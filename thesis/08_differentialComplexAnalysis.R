#' Load data
scoredDataAll <- readRDS("scoredDataAll.rds")
pepTracesList <- readRDS("traces_list_reassignedProteoforms.rds")
protTraces <- readRDS("protein_traces_list.rds")
design_matrix <- readRDS("design_matrix.rda")
calibration <- readRDS("calibration.rds")

#' ## Extract feature values
complex_featureVals <- extractFeatureVals(traces = pepTracesList,
                                          features = scoredDataAll,
                                          design_matrix = design_matrix,
                                          extract = "subunits_detected",
                                          imputeZero = T,
                                          verbose = F,
                                          perturb_cutoff = "5%")

saveRDS(complex_featureVals, "complex_featureVals.rda")

#' ## Fill feature values 
complex_featureValsFilled <- fillFeatureVals(featureVals = complex_featureVals,
                                             tracesList = pepTracesList,
                                             design_matrix = design_matrix)

saveRDS(complex_featureValsFilled, "complex_featureValsFilled.rda")

#' ## Perform peptide-level differential expression testing for all features 
complex_DiffExprPep <- testDifferentialExpression(featureVals = complex_featureValsFilled,
                                                  compare_between = "Condition",
                                                  level = "peptide",
                                                  measuredOnly = FALSE)

#' ## Aggregate differential expression results to the complex level
complex_DiffExprProteoform <- aggregatePeptideTestsToProteoform(complex_DiffExprPep)

complex_DiffExprProtein <- aggregatePeptideTests(complex_DiffExprPep)

complex_DiffExprComplex <- aggregateProteinTests(complex_DiffExprProtein)

#' Change sign of fold-change to accomodate the control as the reference (reverse to how test was calculated)
complex_DiffExprPep[,log2FC:=-log2FC]
complex_DiffExprPep[,global_log2FC_imp:=-global_log2FC_imp]
saveRDS(complex_DiffExprPep, "complex_DiffExprPep.rda")

complex_DiffExprProteoform[,medianLog2FC:=-medianLog2FC]
complex_DiffExprProteoform[,global_medianLog2FC:=-global_medianLog2FC]
complex_DiffExprProteoform[,global_medianLog2FC_imp:=-global_medianLog2FC_imp]
saveRDS(complex_DiffExprProteoform, "complex_DiffExprProteoform.rda")

complex_DiffExprProtein[,medianLog2FC:=-medianLog2FC]
complex_DiffExprProtein[,global_medianLog2FC:=-global_medianLog2FC]
complex_DiffExprProtein[,global_medianLog2FC_imp:=-global_medianLog2FC_imp]
saveRDS(complex_DiffExprProtein, "complex_DiffExprProtein.rda")

complex_DiffExprComplex[,medianLog2FC:=-medianLog2FC]
complex_DiffExprComplex[,global_medianLog2FC:=-global_medianLog2FC]
complex_DiffExprComplex[,global_medianLog2FC_imp:=-global_medianLog2FC_imp]
saveRDS(complex_DiffExprComplex, "complex_DiffExprComplex.rda")

#' Make volcano plots
#' General volcano plots
plotVolcano(complex_DiffExprPep, PDF = T, name = "complex_DiffExprPep")
plotVolcano(complex_DiffExprProteoform, PDF = T, name = "complex_DiffExprProteoform")
plotVolcano(complex_DiffExprProtein, PDF = T, name = "complex_DiffExprProtein")
plotVolcano(complex_DiffExprComplex, PDF = T, name = "complex_DiffExprComplex")

#' Volcanoplts highlighting different complexes
library(ggrepel)
plotVolcano(complex_DiffExprComplex, highlight=c("1142"), PDF = T, name = "complex_DiffExprComplex_1142_SMN")
plotVolcano(complex_DiffExprComplex, highlight=c("659;685"), PDF = T, name = "complex_DiffExprComplex_659;685_MeCP1")

######################################################
######################################################
######################################################

target <- "1142"
pdf("complexFeatures_1142_SMN.pdf",width=6,height=4)
plotFeatures(
  feature_table=scoredDataAll,
  traces=protTraces,
  design_matrix = design_matrix,
  feature_id = target,
  calibration=calibration,
  annotation_label="Entry_name",
  peak_area=TRUE,
  onlyBest = FALSE,
  monomer_MW = TRUE,
  PDF=F,
  name=id
)
dev.off()

complex_stats <- complex_DiffExprComplex[complex_id=="1142"]
globalLog2FC <- median(complex_stats$global_medianLog2FC_imp)
#setnames(complex_stats,c("qint1","qint2"),c("control","depleted"))
complex_stats <- subset(complex_stats, select=c("complex_id","apex","medianLog2FC","pBHadj"))
complex_stats.m <- melt(complex_stats,id.vars=c("complex_id","apex","pBHadj"))
setnames(complex_stats.m,"value","log2FC")
setnames(complex_stats.m,"apex","fraction")

pdf(paste0("complex_stats",".pdf"), height=2, width=3)
ggplot(complex_stats.m,aes(x=fraction,y=log2FC)) + 
  geom_bar(stat="identity",fill="steelblue") +
  theme_classic() +
  scale_x_continuous(breaks=seq(0, max(protTraces$control$fraction_annotation$id), 20), 
                     limits = c(0, max(protTraces$control$fraction_annotation$id))) +
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
  calibration=calibration,
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
  calibration=calibration,
  annotation_label="Entry_name",
  peak_area=TRUE,
  onlyBest = FALSE,
  monomer_MW = TRUE,
  PDF=F,
  name=id
)
dev.off()
