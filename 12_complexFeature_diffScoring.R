#'# Finding differentially expressed proteins between plus and minus

#' ## Environment setup
library(devtools)
#setwd("~/mysonas/CCprofiler")
#load_all()
setwd("~/mysonas/PRPF8/analysis/output")
install_github("CCprofiler/CCprofiler", ref = "DA_module")
library(CCprofiler)

#' ## Load data
proteinFeatures <- readRDS("proteinFeatures_005_FDR_integrated.rds")
#complexFeatures <- readRDS("complexFeaturesBestFiltered_005_FDR.rds")
complexFeatures <- fread("corum_scoredDataAll_005_FDR.txt")
pepTraces <- readRDS("pepTracesRaw.rda")
protTraces <- readRDS("protTraces_assembly.rds")
design_matrix <- readRDS("design_matrix.rda")

complexFeaturesUnique <- getUniqueFeatureGroups(complexFeatures,
                                                rt_height = 0,
                                                distance_cutoff = 1.25)

complexFeaturesCollapsed <- callapseByUniqueFeatureGroups(complexFeaturesUnique,
                                                          rm_decoys = TRUE)

saveRDS(complexFeaturesCollapsed, "complexFeaturesCollapsed.rda")

testIDs <- unique(complexFeaturesCollapsed$complex_id)#[1:20]
testComplexFeatures <- complexFeaturesCollapsed[complex_id%in%testIDs]
#testIDs <- unique(complexFeatures$complex_id)
#testComplexFeatures <- complexFeatures[complex_id%in%testIDs]


complex_featureVals <- extractFeatureVals(traces = pepTraces,
                                  features = testComplexFeatures,
                                  design_matrix = design_matrix,
                                  extract = "subunits",
                                  imputeZero = T,
                                  verbose = F,
                                  perturb_cutoff = "5%")

#complex_featureVals <- extractFeatureVals(traces = pepTraces,
#                                          features = testComplexFeatures,
#                                          design_matrix = design_matrix,
#                                          extract = "subunits_detected",
#                                          imputeZero = T,
#                                          verbose = F,
#                                          perturb_cutoff = "5%")

saveRDS(complex_featureVals, "complex_featureVals.rda")

complex_featureValsFilled <- fillFeatureVals(featureVals = complex_featureVals,
                                     design_matrix = design_matrix)

saveRDS(complex_featureValsFilled, "complex_featureValsFilled.rda")

complex_DiffExprPep <- testDifferentialExpression(featureVals = complex_featureValsFilled,
                          compare_between = "Condition",
                          level = "peptide",
                          measuredOnly = FALSE)

saveRDS(complex_DiffExprPep, "complex_DiffExprPep.rda")

#complex_DiffExprPep <- readRDS("complex_DiffExprPep.rda")

complex_DiffExprProt <- aggregatePeptideTests(complex_DiffExprPep)
saveRDS(complex_DiffExprProt, "complex_DiffExprProt.rda")

complex_DiffExprComplex <- aggregateProteinTests(complex_DiffExprProt)
saveRDS(complex_DiffExprComplex, "complex_DiffExprComplex.rda")

# volcano plots
plotVolcano(complex_DiffExprPep, PDF = T, name = "complex_DiffExprPep")
plotVolcano(complex_DiffExprProt, PDF = T, name = "complex_DiffExprProt")
plotVolcano(complex_DiffExprComplex, PDF = T, name = "complex_DiffExprComplex")

plotVolcano(complex_DiffExprComplex, highlight = "1745;1068;1142;1143;835;3118;2757;834;1743;833", PDF = T, name = "complex_DiffExprComplex_SMN")


# differentially behaving
diffPep <- subset(complex_DiffExprPep, (abs(log2FC) > 1) & (pBHadj < 0.01))
diffProt <- subset(complex_DiffExprProt, (abs(sumLog2FC) > 1) & (pBHadj < 0.01))
diffComplex <- subset(complex_DiffExprComplex, (abs(sumLog2FC) > 1) & (pBHadj < 0.01))

diffPep_shift <- subset(diffPep, abs(local_vs_global_log2FC) > 1)
diffProt_shift <- subset(diffProt, abs(local_vs_global_log2FC) > 1)
diffComplex_shift <- subset(diffComplex, abs(local_vs_global_log2FC) > 1)

diffPep_noShift <- subset(diffPep, abs(local_vs_global_log2FC) < 1)
diffProt_noShift <- subset(diffProt, abs(local_vs_global_log2FC) < 1)
diffComplex_noShift <- subset(diffComplex, abs(local_vs_global_log2FC) < 1)

diffComplex <- diffComplex[order(-abs(local_vs_global_log2FC_imp))]
diffComplex_shift <- diffComplex_shift[order(-abs(local_vs_global_log2FC_imp))]
diffComplex_noShift <- diffComplex_noShift[order(abs(local_vs_global_log2FC_imp))]


pdf("diffComplex_SMN.pdf",width=6,height=4)
for (id in "1745;1068;1142;1143;835;3118;2757;834;1743;833"){ #error for 24th ID????
  sub <- subset(complexFeaturesCollapsed, (complex_id==id))
  plotFeatures.tracesList(sub, protTraces, id, design_matrix = design_matrix, annotation_label="Gene_names",
                        onlyBest = F, peak_area = T, legend=T)
}
dev.off()


pdf("diffComplex_SMN_minus.pdf")
for (id in "1745;1068;1142;1143;835;3118;2757;834;1743;833"){ #error for 24th ID????
  sub <- subset(complexFeaturesCollapsed, (complex_id==id))
  protTraces.sub <- subset(protTraces$minus,trace_subset_ids = unlist(strsplit(sub$subunits,split=";")))
  plotFeatures(sub, protTraces.sub, id, annotation_label="Entry_name",
                        onlyBest = F, peak_area = T, legend=T)
}
dev.off()

pdf("diffComplex_SMN_plus.pdf")
for (id in "1745;1068;1142;1143;835;3118;2757;834;1743;833"){ #error for 24th ID????
  sub <- subset(complexFeaturesCollapsed, (complex_id==id))
  protTraces.sub <- subset(protTraces$plus,trace_subset_ids = unlist(strsplit(sub$subunits,split=";")))
  plotFeatures(sub, protTraces.sub, id, annotation_label="Entry_name",
                        onlyBest = F, peak_area = T, legend=T)
}
dev.off()


pdf("diffComplex_sortByShift_imp.pdf",width=6,height=4)
for (id in unique(diffComplex$complex_id)){ #error for 24th ID????
  sub <- subset(complexFeaturesCollapsed, (complex_id==id) & (apex %in% diffComplex[complex_id==id]$apex))
  plotFeatures.tracesList(sub, protTraces, id, design_matrix = design_matrix, annotation_label="Entry_name",
                        onlyBest = F, peak_area = T, legend=T)
}
dev.off()

pdf("diffComplexShift.pdf")
for (id in unique(diffComplex_shift$complex_id)[-24]){ #error for 24th ID????
  sub <- subset(complexFeatures, (complex_id==id) & (apex %in% diffComplex_shift[complex_id==id]$apex))
  plotFeatures.tracesList(sub, protTraces, id, design_matrix = design_matrix, annotation_label="Entry_name",
                        onlyBest = F, peak_area = T, legend=T)
}
dev.off()

pdf("diffComplexNoShift.pdf")
for (id in unique(diffComplex_noShift$complex_id)){
  sub <- subset(complexFeatures, (complex_id==id) & (apex %in% diffComplex_noShift[complex_id==id]$apex))
  plotFeatures.tracesList(sub, protTraces, id, design_matrix = design_matrix, annotation_label="Entry_name",
                          onlyBest = F, peak_area = T, legend=T)
}
dev.off()
