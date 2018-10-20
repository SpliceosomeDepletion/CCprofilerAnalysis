library(devtools)
setwd("~/mysonas/CCprofiler")
load_all()
setwd("~/mysonas/PRPF8/analysis/output")
#install_github("CCprofiler/CCprofiler", ref = "proteoformAnnotation")
#library(CCprofiler)

traces_list <- readRDS("pepTracesImpSPCmultipep.rda")
design_matrix <- readRDS("design_matrix.rda")

# in case only one iteratioon of outlier removal should be performed
# traces_maxCorr <- filterByMaxCorr(traces_list,
#   cutoff = 0.85, plot = T, PDF=T, name="maxCorrHist")

traces_maxCorrItr <- iterativeMaxCorrFilter(traces_list,
  cutoff = 0.85, plot = T, PDF=T, name="maxCorrHist")

traces_maxCorr_multi <-  filterSinglePeptideHits(traces_maxCorrItr)

traces_maxCorr_multi_minCorr <- calculateMinCorr(traces_maxCorr_multi,
  plot = TRUE, PDF=TRUE)

traces_maxCorr_multi_minCorr_pval <- estimateProteoformPval(traces_maxCorr_multi_minCorr,
  plot = TRUE, PDF=TRUE)

traces_maxCorr_multi_minCorr_pval_clustered_multiCond <- clustPepMultiCond(traces_maxCorr_multi_minCorr_pval,
  clusterN = NULL, clusterH = 0.75, plot = TRUE, PDF=TRUE)

saveRDS(traces_maxCorr_multi_minCorr_pval_clustered_multiCond,"pepTraces_proteoform.rds")

# in case clustering should be performed for each condition separately
# traces_maxCorr_multi_minCorr_pval_clustered <- clusterPeptides(traces_maxCorr_multi_minCorr_pval,
#   plot = TRUE, PDF=TRUE)

# How many of the proteins that were significant in the clustering,
# actually cluster in multiple proteoforms?
sig_minus <- traces_maxCorr_multi_minCorr_pval_clustered_multiCond$minus$trace_annotation[proteoform_pval_adj <= 0.05]
sig_minus[, n_proteoforms := length(unique(proteoform_id)), by=c("protein_id")]
length(unique(sig_minus$protein_id))
table(unique(subset(sig_minus,select=c("protein_id","n_proteoforms")))$n_proteoforms)

non_sig_minus <- traces_maxCorr_multi_minCorr_pval_clustered_multiCond$minus$trace_annotation[proteoform_pval_adj > 0.05]
non_sig_minus[, n_proteoforms := length(unique(proteoform_id)), by=c("protein_id")]
length(unique(non_sig_minus$protein_id))
table(unique(subset(non_sig_minus,select=c("protein_id","n_proteoforms")))$n_proteoforms)

sig_plus <- traces_maxCorr_multi_minCorr_pval_clustered_multiCond$plus$trace_annotation[proteoform_pval_adj <= 0.05]
sig_plus[, n_proteoforms := length(unique(proteoform_id)), by=c("protein_id")]
length(unique(sig_plus$protein_id))
table(unique(subset(sig_plus,select=c("protein_id","n_proteoforms")))$n_proteoforms)

non_sig_plus <- traces_maxCorr_multi_minCorr_pval_clustered_multiCond$plus$trace_annotation[proteoform_pval_adj > 0.05]
non_sig_plus[, n_proteoforms := length(unique(proteoform_id)), by=c("protein_id")]
length(unique(non_sig_plus$protein_id))
table(unique(subset(non_sig_plus,select=c("protein_id","n_proteoforms")))$n_proteoforms)
