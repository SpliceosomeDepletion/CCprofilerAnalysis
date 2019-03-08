#' ---
#' title: "Proteoform feature finding"
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
peptide_traces_list <- readRDS("pepTraces_proteoformMulti.rds")
design_matrix <- readRDS("design_matrix.rda")
pepTracesSum <- readRDS("pepTracesSum.rda")
calibrationFunctions <- readRDS("calibration.rds")
filteredDataProteoformResolved <- readRDS("filteredDataProteoformResolved.rds")

proteoforms <- unique(pepTracesSum$trace_annotation[grep("_",proteoform_id)]$proteoform_id)

proteoforms_withFeature <- unique(filteredDataProteoformResolved[grep("_",proteoform_ids)]$proteoform_ids)

library(splitstackshape)

#install_github("js229/Vennerable"); library(Vennerable);
library("Vennerable")

#install.packages('mgsub')
library('mgsub')

#traces <- copy(pepTracesSum)
#features <- copy(filteredDataProteoformResolved)


reassigned <- refineProteoformsByDetectedFeatures(traces = pepTracesSum, 
                                                  features = filteredDataProteoformResolved)

reassignedTraces <- reassigned[[1]]
reassignedFeatures <- reassigned[[2]]

length(unique(reassignedTraces$trace_annotation$protein_id))
length(unique(reassignedTraces$trace_annotation[n_proteoforms>1]$protein_id))
length(unique(reassignedTraces$trace_annotation[n_proteoforms>1]$protein_id))

#old_traces <- copy(peptide_traces_list)

#peptide_traces_list$minus$trace_annotation[,proteoform_id:=NULL]
#peptide_traces_list$minus$trace_annotation[,n_proteoforms:=NULL]
#peptide_traces_list$minus$trace_annotation[,cluster:=NULL]
#peptide_traces_list$plus$trace_annotation[,proteoform_id:=NULL]
#peptide_traces_list$plus$trace_annotation[,n_proteoforms:=NULL]
#peptide_traces_list$plus$trace_annotation[,cluster:=NULL]

reassignedTraces_multi <- annotateProteoformsAcrossConditions(peptide_traces_list,
                                                              reassignedTraces)


saveRDS(reassignedTraces_multi,"reassignedTraces_multi.rds")
saveRDS(reassignedFeatures,"reassignedFeatures.rds")

#' ## Remove peptides with no proteoform annotation
#removeUnassignedPeptides <- function(traces){
#  proteoforms_assigned <- unique(traces$trace_annotation[!is.na(proteoform_id)]$id)
#  traces_new <- subset(traces, trace_subset_ids = proteoforms_assigned)
#  return(traces_new)
#}

#reassignedTraces_multi <- lapply(reassignedTraces_multi, removeUnassignedPeptides)


#' ## Plot clusters for some example proteins
test_proteins <- c("P22102","O00468","Q9UBF2",
                   "O15067","O75822","Q9Y266",
                   "P42167","P61978","O95801",
                   "P55060","O14578","O75717","P35221")
for (p in test_proteins){
  plotPeptideCluster(reassignedTraces_multi$minus,p)
}

pdf(paste0("reassignedTraces.pdf"),width=8,height=4)
for (id in test_proteins){
  test <- subset(reassignedTraces_multi, trace_subset_ids = id, trace_subset_type = "protein_id")
  plot(traces = test,
       colour_by="proteoform_id",
       legend = T)
}
dev.off()


#test_proteins <- unique(reassignedTraces_multi$minus$trace_annotation[n_proteoforms_beforeReassignment>n_proteoforms]$protein_id)
test_proteins <- unique(reassignedTraces_multi$minus$trace_annotation[n_proteoforms>1]$protein_id)[1:100]
testTraces <- subset(reassignedTraces_multi$minus, trace_subset_ids = test_proteins, trace_subset_type = "protein_id")

testTraces <- copy(reassignedTraces_multi$minus)

#traces_exon_pval <- evaluateExonLocation(testTraces, adj.method = "fdr", optional_filter = F)
traces_exon_pval_filter <- evaluateExonLocation(testTraces, adj.method = "fdr", optional_filter = T)
saveRDS(traces_exon_pval_filter,"traces_exon_pval_filter.rds")

max05 <- subset(traces_exon_pval_filter$trace_annotation, (min_possible_pval<=0.05))
no_minZero_max05 <- subset(traces_exon_pval_filter$trace_annotation, (min_possible_pval!=0) & (min_possible_pval<=0.05))
no_minZero <- subset(traces_exon_pval_filter$trace_annotation, min_possible_pval!=0)

hist(max05, breaks = 30)
hist(no_minZero_max05, breaks = 30)
hist(no_minZero, breaks = 30)

#all tested proteins
length(unique(traces_exon_pval$trace_annotation$protein_id))

# proteins with only one proteoform
length(unique(traces_exon_pval$trace_annotation[is.na(exon_pval_adj)]$protein_id))
# proteins with significant exon location p-value
length(unique(traces_exon_pval$trace_annotation[exon_pval_adj <= 0.05]$protein_id))
# proteins reaching minimum possible exon location p-value
length(unique(traces_exon_pval$trace_annotation[(exon_pval_adj > 0.05) & (exon_pval<=min_possible_pval)]$protein_id))
# proteins not explainable by exon location!
length(unique(traces_exon_pval$trace_annotation[(exon_pval_adj > 0.05) & (exon_pval>min_possible_pval)]$protein_id))




pdf(paste0("reassignedFeatures_conditionSum.pdf"),width=8,height=4)
for (id in unique(reassignedFeatures[n_proteoform_ids>1]$protein_id)[1:40]){
  plotFeatures(feature_table = reassignedFeatures,
               traces = reassignedTraces,
               calibration=calibrationFunctions,
               feature_id = id,
               annotation_label="proteoform_id",
               colour_by="proteoform_id",
               peak_area = T,
               legend = T,
               onlyBest = F)
}
dev.off()

pdf(paste0("reassignedFeatures_ORIGINAL_conditionSum.pdf"),width=8,height=4)
for (id in unique(reassignedFeatures[n_proteoform_ids>1]$protein_id)[1:40]){
  plotFeatures(feature_table = filteredDataProteoformResolved,
               traces = pepTracesSum,
               calibration=calibrationFunctions,
               feature_id = id,
               annotation_label="proteoform_id",
               colour_by="proteoform_id",
               peak_area = T,
               legend = T,
               onlyBest = F)
}
dev.off()
