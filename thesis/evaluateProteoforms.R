#' Load data
reassignedTraces <- readRDS("traces_sum_reassignedProteoforms.rds")
traces_genomic_coordinates <- readRDS("pepTracesSum_filtered_coordinates.rds")
genomic_coordinates <- traces_genomic_coordinates$genomic_coord
reassignedTraces$genomic_coord <- genomic_coordinates[reassignedTraces$traces$id]

reassignedTraces$trace_annotation[, n_peptides := .N, by="protein_id"]

#' ## Evaluate agreement of proteoforms with known exons
traces_exon_pval <- evaluateExonLocation(reassignedTraces)
saveRDS(traces_exon_pval,"traces_exon_pval.rds")

n_significant_proteoforms <- length(unique(traces_exon_pval$trace_annotation[exon_pval <= 0.1]$protein_id))
n_proteoforms <- length(unique(traces_exon_pval$trace_annotation[! is.na(exon_pval)]$protein_id))
n_significant_proteoforms/n_proteoforms

#' ## Evaluate genomic proximity within proteoforms
traces_location_pval <- evaluateProteoformLocation(traces_exon_pval)
saveRDS(traces_location_pval,"traces_location_pval.rds")

n_significant_proteoforms <- length(unique(traces_location_pval$trace_annotation[genomLocation_pval_adj <= 0.1]$proteoform_id))
n_proteoforms <- length(unique(traces_location_pval$trace_annotation[! is.na(genomLocation_pval_adj)]$proteoform_id))
n_significant_proteoforms/n_proteoforms

#' plot exons and clusters
library(ggpubr)
library(RColorBrewer)
plotPeptideCluster(traces_exon_pval,"O75822",PDF=T, closeGaps=T)
plotPeptideCluster(traces_exon_pval,"P27708",PDF=T, closeGaps=T)
plotPeptideCluster(traces_exon_pval,"Q13263",PDF=T, closeGaps=T)
plotPeptideCluster(traces_exon_pval,"Q15056",PDF=T, closeGaps=T)
plotPeptideCluster(traces_exon_pval,"Q9Y6K5",PDF=T, closeGaps=T)
plotPeptideCluster(traces_exon_pval,"Q9NZI8",PDF=T, closeGaps=T)
plotPeptideCluster(traces_exon_pval,"P42167",PDF=T, closeGaps=T)
plotPeptideCluster(traces_exon_pval,"P14618",PDF=T, closeGaps=T)
plotPeptideCluster(traces_exon_pval,"O43390",PDF=T, closeGaps=T)
plotPeptideCluster(traces_exon_pval,"P61978",PDF=T, closeGaps=T)


reassignedFeatures <- readRDS("reassignedProteoformFeatures.rds")
design_matrix <- readRDS("design_matrix.rda")
calibrationFunctions <- readRDS("calibration.rds")
reassignedTraces_list <- readRDS("traces_list_reassignedProteoforms.rds")

id <- "Q9Y6K5"
pdf("proteoform_Q9Y6K5.pdf", height=5, width=8)
sub <- subset(reassignedTraces_list, trace_subset_ids = id, trace_subset_type = "protein_id")
plotPeptideCluster(traces_exon_pval,id, closeGaps=T)
plotFeatures(feature_table = reassignedFeatures,
             traces = reassignedTraces_list,
             calibration=calibrationFunctions,
             feature_id = id,
             design_matrix=design_matrix,
             annotation_label="proteoform_id",
             colour_by="proteoform_id",
             peak_area = T,
             legend = T,
             onlyBest = F,
             monomer_MW=F)
plotFeatures(feature_table = reassignedFeatures,
             traces = reassignedTraces_list,
             calibration=calibrationFunctions,
             feature_id = id,
             design_matrix=design_matrix,
             annotation_label="isoform_id",
             colour_by="isoform_id",
             peak_area = T,
             legend = T,
             onlyBest = F,
             monomer_MW=F)
dev.off()

proteoform_genes <- c("P42167","P14618","P61978")
pdf("proteoform_LAP2AB_PRM_hnrnpr.pdf", height=5, width=8)
for (id in proteoform_genes) {
  sub <- subset(reassignedTraces_list, trace_subset_ids = id, trace_subset_type = "protein_id")
  plotPeptideCluster(traces_exon_pval,id, closeGaps=T)
  plotFeatures(feature_table = reassignedFeatures,
               traces = reassignedTraces_list,
               calibration=calibrationFunctions,
               feature_id = id,
               design_matrix=design_matrix,
               annotation_label="proteoform_id",
               colour_by="proteoform_id",
               peak_area = T,
               legend = T,
               onlyBest = F,
               monomer_MW=F)
  plotFeatures(feature_table = reassignedFeatures,
               traces = reassignedTraces_list,
               calibration=calibrationFunctions,
               feature_id = id,
               design_matrix=design_matrix,
               annotation_label="isoform_id",
               colour_by="isoform_id",
               peak_area = T,
               legend = T,
               onlyBest = F,
               monomer_MW=F)
}
dev.off()


pdf("genomLocation_pval_histogram.pdf")
  hist(traces_location_pval$trace_annotation$genomLocation_pval, breaks = 50)
  hist(traces_location_pval$trace_annotation$genomLocation_pval_adj, breaks = 50)
dev.off()

#' Perform sequence enrichment analysis 
allDetectedProteins <- unique(traces_genomic_coordinates$trace_annotation$protein_id)
write.table(allDetectedProteins,"allDetectedProteins.txt",quote = F, col.names = F, row.names = F, sep="\t")
proteoformProteins <- unique(reassignedTraces$trace_annotation[grep("_",proteoform_id)]$protein_id)
write.table(proteoformProteins,"proteoformProteins.txt",quote = F, col.names = F, row.names = F, sep="\t")

david_proteoformProteins_vs_allDetectedProteins_UpSeqFeature <- fread("david_proteoformProteins_vs_allDetectedProteins_UpSeqFeature.txt")

plotEnrichment <- function(david_enrichment, enrichment_cutoff = 1.5, Bonferroni_cutoff = 0.05, name="david_enrichment", PDF=T) {
  names(david_enrichment) <- gsub(" ", "_", names(david_enrichment))
  david_enrichment_sig <- subset(david_enrichment, Fold_Enrichment >= enrichment_cutoff & Bonferroni <= Bonferroni_cutoff)
  david_enrichment_sig[, neglog10Bonferroni := -log10(Bonferroni)]
  cat <- david_enrichment_sig$Term
  david_enrichment_sig$Term <- factor(david_enrichment_sig$Term, levels=cat)
  david_enrichment_sig[,name:=gsub("^GO:.*~","",Term)]
  david_enrichment_sig <- david_enrichment_sig[order(Fold_Enrichment)]
  david_enrichment_sig$name <- factor(david_enrichment_sig$name, levels = david_enrichment_sig$name)
  if(nrow(david_enrichment_sig) > 10) {
    david_enrichment_sig <- david_enrichment_sig[1:10]
  }
  
  q <- ggplot(data=david_enrichment_sig,aes(x=name,y=Fold_Enrichment, fill=neglog10Bonferroni)) +
    geom_bar(stat = "identity") +
    theme_classic() +
    theme(axis.text.y  = element_text(angle = 20, hjust = 1)) +
    labs(fill='-log10Bonferroni',x='Category',y='Fold Enrichment') +
    coord_flip() +
    geom_text(aes(label=round(Fold_Enrichment,digits=1)), hjust=1.2,color="white") 
  if (PDF){
    pdf(paste0(name,"_enrichment_cutoff_",enrichment_cutoff,"_Bonferroni_cutoff_",Bonferroni_cutoff,".pdf"),height=4,width=6)
  }
  plot(q)
  if (PDF){
    dev.off()
  }
}

plotEnrichment(david_proteoformProteins_vs_allDetectedProteins_UpSeqFeature,
               name="david_proteoformProteins_vs_allDetectedProteins_UpSeqFeature", 
               enrichment_cutoff = 1.25, Bonferroni_cutoff = 0.05, PDF=T)


###########################
###########################
###########################
###########################
###########################






#allProts <- unique(reassignedTraces$trace_annotation$protein_id)
#reassignedTraces <- subset(reassignedTraces, allProts[1:10], trace_subset_type="protein_id")

#' ## Evaluate agreement of proteoforms with known exons
#traces_exon_pval <- evaluateExonLocation(reassignedTraces)
#saveRDS(traces_exon_pval,"traces_exon_pval.rds")

#n_significant_proteoforms <- length(unique(traces_exon_pval$trace_annotation[exon_pval <= 0.1]$protein_id))
#n_proteoforms <- length(unique(traces_exon_pval$trace_annotation[! is.na(exon_pval)]$protein_id))
#n_significant_proteoforms/n_proteoforms


#' ## Proteoform stats
#proteoform_stats <- copy(traces_exon_pval$trace_annotation)
#' Number of proteins
#length(unique(proteoform_stats$protein_id))
#' Number of proteins with > 1 proteoform
#length(unique(proteoform_stats[proteoform_pval_adj < 0.05]$protein_id))
#' Number of proteoforms per protein
#multi_proteoforms <- subset(proteoform_stats, proteoform_pval_adj < 0.05)
#pf_stat <- unique(subset(multi_proteoforms,select=c("protein_id","n_proteoforms","proteoform_pval_adj")))
#pdf("n_proteoforms.pdf", height=3, width=3)
#ggplot(pf_stat,aes(x=n_proteoforms)) +
#  stat_bin(binwidth=1) +
#  theme_classic() +
#  stat_bin(binwidth=1, geom="text", aes(label=..count..), vjust=-.5,size=2.5)
#dev.off()


#' ##### Number of genes/proteins for which clusters are significantly better explained by genommic location compared to random clustering results.
#' Number of proteins with with >1 proteoforms that can be explained by exons:
#length(unique(multi_proteoforms[exon_pval_adj <= 0.05]$protein_id))
#' Fraction of proteins with with >1 proteoforms that can be explained by exons:
#length(unique(multi_proteoforms[exon_pval_adj <= 0.05]$protein_id))/length(unique(multi_proteoforms$protein_id))

#' ##### Compare number of exons to number of proteoforms
#multi_proteoforms[,exon_proteoform_diff:=nExons-n_proteoforms]
#diffTable <- unique(subset(multi_proteoforms,select = c("protein_id","exon_proteoform_diff")))
#pdf("exon_minus_proteoform_count.pdf")
#p <- ggplot(diffTable,aes(x=exon_proteoform_diff)) +
#  geom_histogram(bins=50) +
#  theme_classic()
#print(p)
#dev.off()

# nProteoforms minus
#getProteoformStats <-function(x){
#  stat <- subset(x$trace_annotation,select=c("protein_id","proteoform_id"))
#  stat[,nProteoforms := length(unique(proteoform_id)),by=protein_id]
#  stat <- unique(stat[,proteoform_id:=NULL])
#  return(stat)
#}

#proteoformStats <- lapply(pepTraces_proteoformMulti,getProteoformStats)
#table(proteoformStats$minus$nProteoforms)
#table(proteoformStats$plus$nProteoforms)





#test <- copy(subset(traces_exon_pval$trace_annotation,n_proteoforms > 1))
#test[,minPepPerProteoform := min(table(proteoform_id)),by=protein_id]
#test[,meanPepPerProteoform := as.numeric(mean(table(proteoform_id))),by=protein_id]
#test <- unique(subset(test,select=c("protein_id","n_peptides","n_proteoforms","exon_pval","nExonMistakes","nExons","exon_pval_adj","minPepPerProteoform","meanPepPerProteoform")))
#test[,pepToExon := n_peptides/nExons]
#test[,exonToProteoform := nExons/n_proteoforms]
#test[,pepToProteoform := n_peptides/n_proteoforms]
#test[,pepToExon := ifelse(pepToExon>quantile(test$pepToExon,na.rm=T,probs=c(0.95)),5,pepToExon)]
#test[,minPepPerProteoform := ifelse(minPepPerProteoform>quantile(test$minPepPerProteoform,na.rm=T,probs=c(0.95)),10,minPepPerProteoform)]
#test[,minCutoff := ifelse(minPepPerProteoform>=3,">=3","<3")]
##test[,n_peptides := ifelse(n_peptides>quantile(test$n_peptides,na.rm=T,probs=c(0.9)),30,n_peptides)]

#ggplot(test,aes(x=exon_pval_adj,y=n_peptides,colour=minCutoff)) + geom_point(alpha=0.7)
