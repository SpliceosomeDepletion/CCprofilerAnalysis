#' ## Load data
pepTracesList_filtered <- readRDS("pepTracesList_filtered.rda")
design_matrix <- readRDS("design_matrix.rda")
calibrationFunctions <- readRDS("calibration.rds")

#' ## Quantify on protein level
protein_traces_list <- proteinQuantification(pepTracesList_filtered,
                                             quantLevel="protein_id",
                                             topN = 1000,
                                             keep_less = TRUE,
                                             rm_decoys = TRUE,
                                             use_sibPepCorr = FALSE,
                                             use_repPepCorr = FALSE,
                                             full_intersect_only = FALSE,
                                             verbose = FALSE)

#' ## Update fraction annotation for protein traces
protein_traces_list <- updateTraces(protein_traces_list)

#' ## Summarize protein traces across conditions
protein_traces_sum <- integrateTraceIntensities(protein_traces_list,
                                                design_matrix = NULL,
                                                integrate_within = NULL,
                                                aggr_fun = "sum")

#' ## Save protein traces
saveRDS(protein_traces_list,"protein_traces_list.rds")
saveRDS(protein_traces_sum,"protein_traces_sum.rds")

#corumHypotheses <- readRDS("../../../html/SECpaper/output_data/corum/complexTargetsPlusDecoys.rds")
#saveRDS(corumHypotheses,"corumHypotheses.rds")
corumHypotheses <- readRDS("corumHypotheses.rds")

#' ## Complex feature finding
complexFeatures <-  findComplexFeatures(traces = protein_traces_sum,
                                        complex_hypothesis = corumHypotheses,
                                        corr_cutoff = 0.9,
                                        window_size = 7,
                                        parallelized = F,
                                        n_cores = 1,
                                        collapse_method = "apex_network",
                                        perturb_cutoff = "5%",
                                        rt_height = 2,
                                        smoothing_length = 7)

saveRDS(complexFeatures,"corumComplexFeatures.rda")

complexFeatures <- filterFeatures(complexFeatures,
                                  complex_ids = NULL,
                                  protein_ids = NULL,
                                  min_feature_completeness = NULL,
                                  min_hypothesis_completeness = NULL,
                                  min_subunits = NULL,
                                  min_peak_corr = NULL,
                                  min_monomer_distance_factor = 2)

hypothesis <- "corum"

plotSummarizedMScoverage(hypotheses = corumHypotheses, protein_traces_sum, PDF = T)

complexFeaturesBest <- getBestFeatures(complexFeatures)

complexFeaturesBestFiltered <- scoreFeatures(complexFeaturesBest, FDR=0.05, PDF=T, name="qvalueStats_complexFeatures")

#featuresScored <- calculateCoelutionScore(complexFeaturesBest)
#qvalueFeaturesScored <- calculateQvalue(featuresScored, name="qvalueStats", PDF=F)

scoredDataAll <- appendSecondaryComplexFeatures(scoredPrimaryFeatures = complexFeaturesBestFiltered, 
                                                allFeatures = complexFeatures, peakCorr_cutoff = 0.8)
saveRDS(scoredDataAll,"scoredDataAll.rds")

plotSummarizedComplexes(scoredDataAll, corumHypotheses, protein_traces_sum, PDF=T, name="complex_completeness_pie")

summarizeFeatures(scoredDataAll,
                  plot=TRUE,
                  PDF=T,
                  name="feature_summary_scoredDataAll")

summarizeFeatures(complexFeaturesBestFiltered,
                  plot=TRUE,
                  PDF=T,
                  name="feature_summary_complexFeaturesBestFiltered")

interest <- c("COP9|proteasome|CCT|Septin|SMN|Prefoldin|snRNP|splice|Splice|DISC|CASP|EIF|RNA")
targets <- unique(scoredDataAll[grep(interest,complex_name)]$complex_id)

pdf("complexFeatures_merged.pdf",width=5,height=4)
for(id in targets){
  plotFeatures(
    feature_table=scoredDataAll,
    traces=protein_traces_sum,
    feature_id = id,
    calibration=calibrationFunctions,
    peak_area=TRUE,
    onlyBest = FALSE,
    monomer_MW = TRUE,
    PDF=F,
    name=id
  )
}
dev.off()

pdf("complexFeatures_separate.pdf",width=5,height=4)
for(id in targets){
  plotFeatures(
    feature_table=scoredDataAll,
    traces=protein_traces_list,
    design_matrix = design_matrix,
    feature_id = id,
    calibration=calibrationFunctions,
    peak_area=TRUE,
    onlyBest = FALSE,
    monomer_MW = TRUE,
    PDF=F,
    name=id
  )
}
dev.off()

