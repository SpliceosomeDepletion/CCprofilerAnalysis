#' ---
#' title: "Proteoform annotation"
#' author: "Isabell Bludau"
#' date: "October 21th, 2018"
#' output:
#'   html_document:
#'     toc: true
#'   html_notebook:
#'     toc: true
#'   pdf_document:
#'     toc: true
#' ---

#' # Overview
#' In this script we ...
#'
#' # Step-by-step workflow
#'
#' ## Load CCprofiler package and set working directory:
library(devtools)
library(nFactors)
#install_github("CCprofiler/CCprofiler", ref = "DA_module")
#library(CCprofiler)
if (length(grep("nas21.ethz.ch",getwd()))>0) {
  setwd("~/mysonas/CCprofiler")
  load_all()
  setwd("~/mysonas/PRPF8/analysis/output")
  knitr::opts_knit$set(root.dir = '~/mysonas/PRPF8/analysis/output')
} else {
  setwd("/Volumes/ibludau-1/CCprofiler")
  load_all()
  setwd("/Volumes/ibludau-1/PRPF8/analysis/output")
  knitr::opts_knit$set(root.dir = '/Volumes/ibludau-1/PRPF8/analysis/output')
}

#' ## Load data
traces_list <- readRDS("pepTracesImpSPCmultipep.rda")

#' ## Combine traces across conditions
#' Combine traces across conditions by appending fractions to generate a
#' large traces object with twice as many fractions. This is done in order to
#' allow consistent proteoform detection across conditions.
traces_combined <- combineTracesMutiCond(traces_list)

#' ## Remove outlier peptides
#' Outlier peptides are defined as peptides that do not have any high correlating
#' sibling peptide. Outliers are removed iteratively untill no more peptides
#' fall below the cutoff.
traces_maxCorr <- filterByMaxCorr(traces_combined,
  cutoff = 0.8, plot = T, PDF=T, name="maxCorrHist")

#' ## Remove proteins that only contain a single peptide after outlier filtering
traces_maxCorr_multi <-  filterSinglePeptideHits(traces_maxCorr) # actually this doesn't make a difference anymore
saveRDS(traces_maxCorr_multi,"traces_maxCorr_multi.rda")

#' ## Compute the minimum correlation between any two sibling peptides
traces_maxCorr_multi_minCorr <- calculateMinCorr(traces_maxCorr_multi,
  plot = TRUE, PDF=TRUE)

#' ## Estimate proteoform p-values for each gene/protein
#' The null hypothesis is that each gene/protein only consists of one
#' proteoform, meaning that all its sibling peptides have a high correlation
#' along the fractionation dimension. A gamma distribution is fitted to the
#' high-correlating end of the minimum correlation value distribution and
#' p-values are estimated accordingly.
traces_maxCorr_multi_minCorr_pval <- estimateProteoformPval(traces_maxCorr_multi_minCorr,
  plot = TRUE, PDF=TRUE)

#' ## Determine proteoforms by peptide clustering
#' Peptides of each gene/protein are clustered based on the correlation of all sibling peptides.
traces_maxCorr_multi_minCorr_pval_clustered <- clusterPeptides(traces_maxCorr_multi_minCorr_pval,
  clusterN = NULL, clusterH = 0.6, nFactorAnalysis = FALSE,
  plot = TRUE, PDF=TRUE)

#' ## Evaluate clustering results
#' ##### Differentiate genes/proteins with one vs. multiple proteoforms based on minimum peptide correlation criterion.
single_proteoform <- traces_maxCorr_multi_minCorr_pval_clustered$trace_annotation[proteoform_pval_adj > 0.05]
multi_proteoforms <- traces_maxCorr_multi_minCorr_pval_clustered$trace_annotation[proteoform_pval_adj <= 0.05]

single_proteoform[,n_proteoforms:=length(unique(proteoform_id)),by="protein_id"]
multi_proteoforms[,n_proteoforms:=length(unique(proteoform_id)),by="protein_id"]

#' ##### Number of clusters for proteins with 1 vs multiple proteoforms:
#' Contingency table of clusters in genes/proteins with 1 proteoform:
table(unique(subset(single_proteoform,select=c("protein_id","n_proteoforms")))$n_proteoforms)
#' Contingency table of clusters in genes/proteins with >1 proteoform:
table(unique(subset(multi_proteoforms,select=c("protein_id","n_proteoforms")))$n_proteoforms)

#' ## Protein with a proteoform p-value <0.05 should be set to 1 proteoform regardless of clustering
traces_maxCorr_multi_minCorr_pval_clustered_fixed <- copy(traces_maxCorr_multi_minCorr_pval_clustered)
traces_maxCorr_multi_minCorr_pval_clustered_fixed$trace_annotation[,proteoform_id := ifelse(proteoform_pval_adj > 0.05, protein_id, proteoform_id)]
traces_maxCorr_multi_minCorr_pval_clustered_fixed$trace_annotation[,n_proteoforms := length(unique(proteoform_id)), by="protein_id"]

#' ## Plot clusters for some example proteins
test_proteins <- c("P22102","O00468","Q9UBF2",
                   "O15067","O75822","Q9Y266",
                   "P42167","P61978","O95801",
                   "P55060","O14578","O75717","P35221")
for (p in test_proteins){
  plotPeptideCluster(traces_maxCorr_multi_minCorr_pval_clustered_fixed,p)
}


#' ## Transfer proteoform annotation to original traces list
#' The proteoform detection was performed on the combined traces object. To enable
#' condition dependent analysis of he proteoforms the proteoform information is
#' tranferred to the original tarces list object.
pepTraces_proteoformMulti <- annotateProteoformsAcrossConditions(traces_list,
                                                                 traces_maxCorr_multi_minCorr_pval_clustered_fixed)

#' ## Remove peptides with no proteoform annotation
#peptides_annotated_minus <- pepTraces_proteoformMulti$minus$trace_annotation[which(!is.na(pepTraces_proteoformMulti$minus$trace_annotation$proteoform_id))]$id
#peptides_annotated_plus <- pepTraces_proteoformMulti$plus$trace_annotation[which(!is.na(pepTraces_proteoformMulti$plus$trace_annotation$proteoform_id))]$id
#pepTraces_proteoformMulti$minus <- subset(pepTraces_proteoformMulti$minus, trace_subset_ids = peptides_annotated_minus)
#pepTraces_proteoformMulti$plus <- subset(pepTraces_proteoformMulti$plus, trace_subset_ids = peptides_annotated_plus)

#' ## Evaluate analysis results
#' ##### Differentiate genes/proteins with one vs. multiple proteoforms based on minimum peptide correlation criterion.
single_proteoform <- traces_maxCorr_multi_minCorr_pval_clustered_fixed$trace_annotation[proteoform_pval_adj > 0.05]
multi_proteoforms <- traces_maxCorr_multi_minCorr_pval_clustered_fixed$trace_annotation[proteoform_pval_adj <= 0.05]
cat(paste0("Genes/proteins with 1 proteoform: ",length(unique(single_proteoform$protein_id))))
cat(paste0("Genes/proteins with >1 proteoform: ",length(unique(multi_proteoforms$protein_id))))

#' ##### Number of clusters for proteins with 1 vs multiple proteoforms:
#' Contingency table of clusters in genes/proteins with 1 proteoform:
table(unique(subset(single_proteoform,select=c("protein_id","n_proteoforms")))$n_proteoforms)
#' Contingency table of clusters in genes/proteins with >1 proteoform:
table(unique(subset(multi_proteoforms,select=c("protein_id","n_proteoforms")))$n_proteoforms)

#' ## Save objects
saveRDS(traces_maxCorr_multi_minCorr_pval_clustered_fixed,"pepTraces_proteoformMulti_combined.rds")
saveRDS(pepTraces_proteoformMulti,"pepTraces_proteoformMulti.rds")



####
####


#' ## Evaluate agreement of proteoforms with known exons
traces_exon_pval <- evaluateExonLocation(traces_maxCorr_multi_minCorr_pval_clustered_fixed)

#' ## Proteoform stats
proteoform_stats <- copy(traces_exon_pval$trace_annotation)
#' Number of proteins
length(unique(proteoform_stats$protein_id))
#' Number of proteins with > 1 proteoform
length(unique(proteoform_stats[proteoform_pval_adj < 0.05]$protein_id))
#' Number of proteoforms per protein
multi_proteoforms <- subset(proteoform_stats, proteoform_pval_adj < 0.05)
pf_stat <- unique(subset(multi_proteoforms,select=c("protein_id","n_proteoforms","proteoform_pval_adj")))
pdf("n_proteoforms.pdf", height=3, width=3)
ggplot(pf_stat,aes(x=n_proteoforms)) +
  stat_bin(binwidth=1) +
  theme_classic() +
  stat_bin(binwidth=1, geom="text", aes(label=..count..), vjust=-.5,size=2.5)
dev.off()


#' ##### Number of genes/proteins for which clusters are significantly better explained by genommic location compared to random clustering results.
#' Number of proteins with with >1 proteoforms that can be explained by exons:
length(unique(multi_proteoforms[exon_pval_adj <= 0.05]$protein_id))
#' Fraction of proteins with with >1 proteoforms that can be explained by exons:
length(unique(multi_proteoforms[exon_pval_adj <= 0.05]$protein_id))/length(unique(multi_proteoforms$protein_id))

#' ##### Compare number of exons to number of proteoforms
multi_proteoforms[,exon_proteoform_diff:=nExons-n_proteoforms]
diffTable <- unique(subset(multi_proteoforms,select = c("protein_id","exon_proteoform_diff")))
pdf("exon_minus_proteoform_count.pdf")
p <- ggplot(diffTable,aes(x=exon_proteoform_diff)) +
  geom_histogram(bins=50) +
  theme_classic()
print(p)
dev.off()

# nProteoforms minus
getProteoformStats <-function(x){
  stat <- subset(x$trace_annotation,select=c("protein_id","proteoform_id"))
  stat[,nProteoforms := length(unique(proteoform_id)),by=protein_id]
  stat <- unique(stat[,proteoform_id:=NULL])
  return(stat)
}

proteoformStats <- lapply(pepTraces_proteoformMulti,getProteoformStats)
table(proteoformStats$minus$nProteoforms)
table(proteoformStats$plus$nProteoforms)





test <- copy(subset(traces_exon_pval$trace_annotation,n_proteoforms > 1))
test[,minPepPerProteoform := min(table(proteoform_id)),by=protein_id]
test[,meanPepPerProteoform := as.numeric(mean(table(proteoform_id))),by=protein_id]
test <- unique(subset(test,select=c("protein_id","n_peptides","n_proteoforms","exon_pval","nExonMistakes","nExons","exon_pval_adj","minPepPerProteoform","meanPepPerProteoform")))
test[,pepToExon := n_peptides/nExons]
test[,exonToProteoform := nExons/n_proteoforms]
test[,pepToProteoform := n_peptides/n_proteoforms]
test[,pepToExon := ifelse(pepToExon>quantile(test$pepToExon,na.rm=T,probs=c(0.95)),5,pepToExon)]
test[,minPepPerProteoform := ifelse(minPepPerProteoform>quantile(test$minPepPerProteoform,na.rm=T,probs=c(0.95)),10,minPepPerProteoform)]
test[,minCutoff := ifelse(minPepPerProteoform>=3,">=3","<3")]
#test[,n_peptides := ifelse(n_peptides>quantile(test$n_peptides,na.rm=T,probs=c(0.9)),30,n_peptides)]

ggplot(test,aes(x=exon_pval_adj,y=n_peptides,colour=minCutoff)) + geom_point(alpha=0.7)
