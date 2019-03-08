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
peptide_traces_list <- readRDS("traces_list_reassignedProteoforms.rds")
design_matrix <- readRDS("design_matrix.rda")
calibrationFunctions <- readRDS("calibration.rds")

#' ## Quantify on protein level
protein_traces_list <- proteinQuantification(peptide_traces_list,
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
                                        rt_height = 1,
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


scoredDataAll <- appendSecondaryComplexFeatures(scoredPrimaryFeatures = complexFeaturesBestFiltered, allFeatures = complexFeatures, peakCorr_cutoff = 0.7)

plotSummarizedComplexes(scoredDataAll, corumHypotheses, protein_traces_sum, PDF=T, name="complex_completeness_pie")

summarizeFeatures(scoredDataAll,
                  plot=TRUE,
                  PDF=T,
                  name="feature_summary_scoredDataAll")

summarizeFeatures(complexFeaturesBestFiltered,
                  plot=TRUE,
                  PDF=T,
                  name="feature_summary_complexFeaturesBestFiltered")

targets <- unique(scoredDataAll$complex_id)[1:20]

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

interest <- c("COP9|proteasome|CCT|Septin|SMN|Prefoldin|snRNP|splice|Splice|DISC|CASP|EIF|RNA")
targets <- unique(scoredDataAll[grep(interest,complex_name)]$complex_id)

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

saveRDS(complex_featureVals, "complex_featureVals.rda")

#' ## Fill feature values 
complex_featureValsFilled <- fillFeatureVals(featureVals = complex_featureVals,
                                                design_matrix = design_matrix)

saveRDS(complex_featureValsFilled, "complex_featureValsFilled.rda")

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
saveRDS(complex_DiffExprComplex, "complex_DiffExprComplex.rda")

#' Make volcano plots
#' General volcano plots
plotVolcano(complex_DiffExprPep, PDF = T, name = "complex_DiffExprPep")
plotVolcano(complex_DiffExprProteoform, PDF = T, name = "complex_DiffExprProteoform")
plotVolcano(complex_DiffExprProtein, PDF = T, name = "complex_DiffExprProtein")
plotVolcano(complex_DiffExprComplex, PDF = T, name = "complex_DiffExprComplex")

plotVolcano(complex_DiffExprComplex, highlight=c("1142"), PDF = T, name = "complex_DiffExprComplex_1142_SMN")
plotVolcano(complex_DiffExprComplex, highlight=c("659;685"), PDF = T, name = "complex_DiffExprComplex_659;685_MeCP1")

protTraces <- readRDS("protein_traces_list.rds")
#complex_DiffExprComplex[complex_id=="659;685"]
target <- "1142"
pdf("complexFeatures_1142_SMN.pdf",width=4,height=4)
  plotFeatures(
    feature_table=scoredDataAll,
    traces=protTraces,
    design_matrix = design_matrix,
    feature_id = target,
    calibration=calibrationFunctions,
    annotation_label="Entry_name",
    peak_area=TRUE,
    onlyBest = FALSE,
    monomer_MW = TRUE,
    PDF=F,
    name=id
  )
dev.off()

complex_stats <- complex_DiffExprComplex[complex_id=="1142"]
globalLog2FC <- median(complex_stats$global_sumLog2FC_imp)
setnames(complex_stats,c("qint1","qint2"),c("minus","plus"))
complex_stats <- subset(complex_stats, select=c("complex_id","apex","sumLog2FC","pBHadj"))
complex_stats.m <- melt(complex_stats,id.vars=c("complex_id","apex","pBHadj"))
setnames(complex_stats.m,"value","log2FC")
setnames(complex_stats.m,"apex","fraction")

pdf(paste0("complex_stats",".pdf"), height=2, width=3)
ggplot(complex_stats.m,aes(x=fraction,y=log2FC)) + 
  geom_bar(stat="identity",fill="steelblue") +
  theme_classic() +
  scale_x_continuous(breaks=seq(0, max(protTraces$minus$fraction_annotation$id), 20), 
                     limits = c(0, max(protTraces$minus$fraction_annotation$id))) +
  geom_hline(yintercept = globalLog2FC, linetype="dashed", colour="#E69F00") +
  geom_hline(yintercept = 0, linetype="solid")
dev.off()




target <- "587"
pdf("complexFeatures_587_NuRD.pdf",width=5,height=4)
plotFeatures(
  feature_table=scoredDataAll,
  traces=protTraces,
  design_matrix = design_matrix,
  feature_id = target,
  calibration=calibrationFunctions,
  annotation_label="Entry_name",
  peak_area=TRUE,
  onlyBest = FALSE,
  monomer_MW = TRUE,
  PDF=F,
  name=id
)
dev.off()


target <- "888-1"
pdf("complexFeatures_888_1_MTA2.pdf",width=5,height=4)
plotFeatures(
  feature_table=scoredDataAll,
  traces=protTraces,
  design_matrix = design_matrix,
  feature_id = target,
  calibration=calibrationFunctions,
  annotation_label="Entry_name",
  peak_area=TRUE,
  onlyBest = FALSE,
  monomer_MW = TRUE,
  PDF=F,
  name=id
)
dev.off()
