#' # Protein feature finding
#'
#' ## Environment setup
#'
library(devtools)
#setwd("~/mysonas/CCprofiler")
#load_all()
setwd("~/mysonas/PRPF8/analysis/output")
install_github("CCprofiler/CCprofiler", ref = "DA_module")
library(CCprofiler)

#' ## Load data
proteinFeatures <- readRDS("proteinFeatures_005_FDR_integrated.rds")
pepTraces <- readRDS("pepTracesRaw.rda")
design_matrix <- readRDS("design_matrix.rda")

prot_featureVals <- extractFeatureVals(traces = pepTraces,
                                  features = proteinFeatures,
                                  design_matrix = design_matrix,
                                  extract = "subunits_detected",
                                  imputeZero = T,
                                  verbose = F,
                                  perturb_cutoff = "5%")

saveRDS(prot_featureVals, "prot_featureVals.rda")

prot_featureValsFilled <- fillFeatureVals(featureVals = prot_featureVals,
                                     design_matrix = design_matrix)

saveRDS(prot_featureValsFilled, "prot_featureValsFilled.rda")

prot_DiffExprPep <- testDifferentialExpression(featureVals = prot_featureValsFilled,
                          compare_between = "Condition",
                          level = "peptide",
                          measuredOnly = FALSE)

saveRDS(prot_DiffExprPep, "prot_DiffExprPep.rda")

prot_DiffExprProt <- aggregatePeptideTests(prot_DiffExprPep)
saveRDS(prot_DiffExprProt, "prot_DiffExprProt.rda")

# volcano plots
plotVolcano(prot_DiffExprPep, PDF = T, name = "prot_DiffExprPep")
plotVolcano(prot_DiffExprProt, PDF = T, name = "prot_DiffExprProt")

# differentially behaving
diffPep <- subset(prot_DiffExprPep, (abs(log2FC) > 1) & (pBHadj < 0.01))
diffProt <- subset(prot_DiffExprProt, (abs(sumLog2FC) > 1) & (pBHadj < 0.01))

diffPep_shift <- subset(diffPep, abs(local_vs_global_log2FC_imp) > 1)
diffProt_shift <- subset(diffProt, abs(local_vs_global_log2FC_imp) > 1)

diffPep_noShift <- subset(diffPep, abs(local_vs_global_log2FC_imp) < 1)
diffProt_noShift <- subset(diffProt, abs(local_vs_global_log2FC_imp) < 1)

diffProt <- diffProt[order(-abs(local_vs_global_log2FC_imp))]
diffProt_shift <- diffProt_shift[order(-abs(local_vs_global_log2FC_imp))]
diffProt_noShift <- diffProt_noShift[order(abs(local_vs_global_log2FC_imp))]

pdf("diffProt_sortByShift_imp.pdf")
for (id in unique(diffProt$feature_id)[1:100]){ #error for 1st ID????
  sub <- subset(proteinFeatures, (protein_id==id) & (apex %in% diffProt[feature_id==id]$apex))
  plotFeatures.tracesList(sub, pepTraces, id, design_matrix = design_matrix,
                        onlyBest = F, peak_area = T, legend=F)
}
dev.off()

pdf("diffProtShift.pdf")
for (id in unique(diffProt_shift$feature_id)[1:100]){ #error for 1st ID????
  sub <- subset(proteinFeatures, (protein_id==id) & (apex %in% diffProt_shift[feature_id==id]$apex))
  plotFeatures.tracesList(sub, pepTraces, id, design_matrix = design_matrix,
                        onlyBest = F, peak_area = T, legend=F)
}
dev.off()

pdf("diffProtNoShift.pdf")
for (id in unique(diffProt_noShift$feature_id)[1:100]){
  sub <- subset(proteinFeatures, (protein_id==id) & (apex %in% diffProt_noShift[feature_id==id]$apex))
  plotFeatures.tracesList(sub, pepTraces, id, design_matrix = design_matrix,
                          onlyBest = F, peak_area = T, legend=F)
}
dev.off()
