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
protTracesSum <- readRDS("protein_traces_sum.rds")
design_matrix <- readRDS("design_matrix.rda")
corumHypotheses <- exampleComplexHypotheses

#' ## Complex feature finding
complexFeatures <-  findComplexFeatures(traces = protTracesSum,
                                        complex_hypothesis = corumHypotheses,
                                        corr_cutoff = 0.9,
                                        window_size = 8,
                                        parallelized = F,
                                        n_cores = 1,
                                        collapse_method = "apex_network",
                                        perturb_cutoff = "5%",
                                        rt_height = 3,
                                        smoothing_length = 7)

saveRDS(complexFeatures,"corumComplexFeatures.rda")

hypothesis <- "corum"
calibrationFunctions <- readRDS("calibration.rds")

plotSummarizedMScoverage(hypotheses = corumHypotheses, protTracesSum, PDF = F)

complexFeaturesBest <- getBestFeatures(complexFeatures)

complexFeaturesBestFiltered <- scoreFeatures(complexFeaturesBest, FDR=0.1, PDF=F, name="qvalueStats_complexFeatures")

scoredDataAll <- appendSecondaryComplexFeatures(scoredPrimaryFeatures = complexFeaturesBestFiltered, allFeatures = complexFeatures, peakCorr_cutoff = 0.5)

plotSummarizedComplexes(scoredDataAll, corumHypotheses, protTracesSum, PDF=F, name="complex_completeness_pie")

summarizeFeatures(scoredDataAll,
                  plot=TRUE,
                  PDF=F,
                  name="feature_summary_scoredDataAll")

summarizeFeatures(complexFeaturesBestFiltered,
                  plot=TRUE,
                  PDF=F,
                  name="feature_summary_complexFeaturesBestFiltered")

targets <- unique(scoredDataAll$complex_id)
#pdf("complexFeatures_merged.pdf",width=8,height=4)
for(id in targets){
  plotFeatures(
    feature_table=scoredDataAll,
    traces=protTracesSum,
    feature_id = id,
    calibration=calibrationFunctions,
    peak_area=TRUE,
    onlyBest = FALSE,
    monomer_MW = TRUE,
    PDF=F,
    name=id
  )
}
#dev.off()

#' ## Extract intensity values for all detected features
#' Load raw peptide-level traces
pepTraces <- readRDS("pepTraces_proteoformMulti.rds")
#' ## Extract feature values
complex_featureVals <- extractFeatureVals(traces = pepTraces,
                                             features = scoredDataAll,
                                             design_matrix = design_matrix,
                                             extract = "subunits_detected",
                                             imputeZero = T,
                                             verbose = F,
                                             perturb_cutoff = "5%")

#' ## Fill feature values 
complex_featureValsFilled <- fillFeatureVals(featureVals = complex_featureVals,
                                                design_matrix = design_matrix)

#' ## Perform peptide-level differential expression testing for all features 
complex_DiffExprPep <- testDifferentialExpression(featureVals = complex_featureValsFilled,
                                                     compare_between = "Condition",
                                                     level = "peptide",
                                                     measuredOnly = FALSE)

saveRDS(complex_DiffExprPep, "complex_DiffExprPep.rda")

#' ## Aggregate differential expression results to the complex level
complex_DiffExprProteoform <- aggregatePeptideTestsToProteoform(complex_DiffExprPep)
saveRDS(complex_DiffExprProteoform, "complex_DiffExprProteoform.rda")

complex_DiffExprProtein <- aggregatePeptideTests(complex_DiffExprPep)
saveRDS(complex_DiffExprProtein, "complex_DiffExprProtein.rda")


complex_DiffExprComplex <- aggregateProteinTests(complex_DiffExprProtein)

#' Make volcano plots
#' General volcano plots
plotVolcano(complex_DiffExprPep, PDF = F, name = "complex_DiffExprPep")
plotVolcano(complex_DiffExprProteoform, PDF = F, name = "complex_DiffExprProteoform")
plotVolcano(complex_DiffExprProtein, PDF = F, name = "complex_DiffExprProtein")
plotVolcano(complex_DiffExprComplex, PDF = F, name = "complex_DiffExprComplex")


#' Volcanoplts highlighting different complexs
plotVolcano(proteoform_DiffExprProteoform, highlight=c("Q6P2Q9_1"), PDF = T, name = "prot_DiffExprProt_PRPF8_Q6P2Q9_1")
plotVolcano(proteoform_DiffExprProteoform, highlight=c("Q6P2Q9_2"), PDF = T, name = "prot_DiffExprProt_PRPF8_Q6P2Q9_2")
plotVolcano(proteoform_DiffExprProteoform, highlight=c("Q6P2Q9_3"), PDF = T, name = "prot_DiffExprProt_PRPF8_Q6P2Q9_3")
