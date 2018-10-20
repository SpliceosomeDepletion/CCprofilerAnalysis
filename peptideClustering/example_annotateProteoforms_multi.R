library(devtools)
setwd("~/mysonas/CCprofiler")
load_all()
setwd("~/mysonas/PRPF8/analysis/output")
#install_github("CCprofiler/CCprofiler", ref = "proteoformAnnotation")
#library(CCprofiler)

traces_list <- readRDS("pepTracesImpSPCmultipep.rda")
design_matrix <- readRDS("design_matrix.rda")


traces_combined <- combineTracesMutiCond(traces_list)

traces_maxCorrItr <- iterativeMaxCorrFilter(traces_combined,
  cutoff = 0.8, plot = T, PDF=T, name="maxCorrHist")

traces_maxCorr_multi <-  filterSinglePeptideHits(traces_maxCorrItr)

traces_maxCorr_multi_minCorr <- calculateMinCorr(traces_maxCorr_multi,
  plot = TRUE, PDF=TRUE)

traces_maxCorr_multi_minCorr_pval <- estimateProteoformPval(traces_maxCorr_multi_minCorr,
  plot = TRUE, PDF=TRUE)

library(nFactors)

traces_maxCorr_multi_minCorr_pval_clustered <- clusterPeptides(traces_maxCorr_multi_minCorr_pval,
  clusterN = NULL, clusterH = NULL, nFactorAnalysis = TRUE,
  plot = TRUE, PDF=TRUE)

traces_list_annotated <- annotateProteoformsAcrossConditions(traces_list,
  traces_maxCorr_multi_minCorr_pval_clustered)

saveRDS(traces_list_annotated,"pepTraces_proteoformMulti.rds")
saveRDS(traces_maxCorr_multi_minCorr_pval_clustered,"pepTraces_proteoformMulti_combined.rds")
