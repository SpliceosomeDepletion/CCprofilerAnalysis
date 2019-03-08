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
traces_list <- readRDS("pepTraces_maxCorr_SPC.rda")
# traces_list <- readRDS("traces_coordinates_list")

#' ## Combine traces across conditions
#' Combine traces across conditions by appending fractions to generate a
#' large traces object with twice as many fractions. This is done in order to
#' allow consistent proteoform detection across conditions.
traces_combined <- combineTracesMutiCond(traces_list)

#' ## Compute the minimum correlation between any two sibling peptides
traces_combined_minCorr <- calculateMinCorr(traces_combined,
  plot = TRUE, PDF=TRUE)

#' ## Estimate proteoform p-values for each gene/protein
#' The null hypothesis is that each gene/protein only consists of one
#' proteoform, meaning that all its sibling peptides have a high correlation
#' along the fractionation dimension. A gamma distribution is fitted to the
#' high-correlating end of the minimum correlation value distribution and
#' p-values are estimated accordingly.
traces_combined_minCorr_pval <- estimateProteoformPval(traces_combined_minCorr,
  plot = TRUE, PDF=TRUE)

#' ## Determine proteoforms by peptide clustering
#' Peptides of each gene/protein are clustered based on the correlation of all sibling peptides.
traces_combined_minCorr_pval_clustered <- clusterPeptides(traces_combined_minCorr_pval,
  clusterN = NULL, clusterH = 0.6, nFactorAnalysis = FALSE,
  plot = TRUE, PDF=TRUE)

#' ## Evaluate clustering results
#' ##### Differentiate genes/proteins with one vs. multiple proteoforms based on minimum peptide correlation criterion.
single_proteoform <- traces_combined_minCorr_pval_clustered$trace_annotation[proteoform_pval_adj > 0.05]
multi_proteoforms <- traces_combined_minCorr_pval_clustered$trace_annotation[proteoform_pval_adj <= 0.05]

single_proteoform[,n_proteoforms:=length(unique(proteoform_id)),by="protein_id"]
multi_proteoforms[,n_proteoforms:=length(unique(proteoform_id)),by="protein_id"]

cat(paste0("Genes/proteins with 1 proteoform: ",length(unique(single_proteoform$protein_id))))
cat(paste0("Genes/proteins with >1 proteoform: ",length(unique(multi_proteoforms$protein_id))))
cat(paste0("# proteoforms: ",length(unique(multi_proteoforms$proteoform_id))))

#' ##### Number of clusters for proteins with 1 vs multiple proteoforms:
#' Contingency table of clusters in genes/proteins with 1 proteoform:
table(unique(subset(single_proteoform,select=c("protein_id","n_proteoforms")))$n_proteoforms)
#' Contingency table of clusters in genes/proteins with >1 proteoform:
table(unique(subset(multi_proteoforms,select=c("protein_id","n_proteoforms")))$n_proteoforms)

#' ## Protein with a proteoform p-value <0.05 should be set to 1 proteoform regardless of clustering
traces_combined_minCorr_pval_clustered_fixed <- copy(traces_combined_minCorr_pval_clustered)
traces_combined_minCorr_pval_clustered_fixed$trace_annotation[,proteoform_id := ifelse(proteoform_pval_adj > 0.05, protein_id, proteoform_id)]
traces_combined_minCorr_pval_clustered_fixed$trace_annotation[,n_proteoforms := length(unique(proteoform_id)), by="protein_id"]

#' ## Plot clusters for some example proteins
test_proteins <- c("P22102","O00468","Q9UBF2",
                   "O15067","O75822","Q9Y266",
                   "P42167","P61978","O95801",
                   "P55060","O14578","O75717","P35221")
for (p in test_proteins){
  plotPeptideCluster(traces_combined_minCorr_pval_clustered_fixed,p, PDF=T)
}

#' ## Transfer proteoform annotation to original traces list
#' The proteoform detection was performed on the combined traces object. To enable
#' condition dependent analysis of he proteoforms the proteoform information is
#' tranferred to the original tarces list object.
traces_list_pepClusters <- annotateProteoformsAcrossConditions(traces_list,
                                                                 traces_combined_minCorr_pval_clustered_fixed)

# These should be zero
nrow(subset(traces_list_pepClusters[[1]]$trace_annotation, (n_proteoforms==1) & (proteoform_pval_adj < 0.05)))
nrow(subset(traces_list_pepClusters[[2]]$trace_annotation, (n_proteoforms==1) & (proteoform_pval_adj < 0.05)))


#' ## Save objects
saveRDS(traces_list_pepClusters,"traces_list_pepClusters.rds")
