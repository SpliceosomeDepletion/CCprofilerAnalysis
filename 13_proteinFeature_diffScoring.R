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


plotVolcano(prot_DiffExprProt, highlight=c("Q9HB71"), PDF = T, name = "prot_DiffExprProt_CACYBP")


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

n_proteins_with_feature <- length(unique(proteinFeatures$protein_id))
n_features <- nrow(proteinFeatures)
n_features_collapsedApex <- nrow(prot_DiffExprProt)
n_proteins_diff <- length(unique(diffProt$feature_id))
n_features_diff <- nrow(diffProt)

prot_DiffExprProt[, differential := ifelse((abs(sumLog2FC) > 1) & (pBHadj < 0.01), TRUE, FALSE)]
prot_DiffExprProt[, shift := ifelse(abs(local_vs_global_log2FC_imp) > 1, TRUE, FALSE)]
prot_DiffExprProt[is.na(local_vs_global_log2FC_imp)]$shift <- FALSE
prot_DiffExprProt[, sign := sign(sumLog2FC)]

diff <- prot_DiffExprProt[differential==T]
diffPositive <- prot_DiffExprProt[differential & (sign==1)]
diffNegative <- prot_DiffExprProt[differential & (sign==-1)]
shift <- prot_DiffExprProt[differential & shift]
noShift <- prot_DiffExprProt[differential & !(shift)]
noShiftNegative <- prot_DiffExprProt[differential & !(shift) & (sign==-1)]
shiftNegative <- prot_DiffExprProt[differential & (shift) & (sign==-1)]
noShiftPositive <- prot_DiffExprProt[differential & !(shift) & (sign==1)]
shiftPositive <- prot_DiffExprProt[differential & (shift) & (sign==1)]

proteins_all <- unique(prot_DiffExprProt$feature_id)
proteins_diff <- unique(diff$feature_id)
proteins_diff_up <- unique(diffPositive$feature_id)
proteins_diff_down <- unique(diffNegative$feature_id)
proteins_shift <- unique(shift$feature_id)
proteins_noShift <- unique(noShift$feature_id)
proteins_noShiftNegative <- unique(noShiftNegative$feature_id)
proteins_shiftNegative <- unique(shiftNegative$feature_id)
proteins_noShiftPositive <- unique(noShiftPositive$feature_id)
proteins_shiftPositive <- unique(shiftPositive$feature_id)

write.table(proteins_all,"proteins_all.txt",sep="\t",quote=F,row.names = F,col.names=F)
write.table(proteins_diff,"proteins_diff.txt",sep="\t",quote=F,row.names = F,col.names=F)
write.table(proteins_diff_up,"proteins_diff_up.txt",sep="\t",quote=F,row.names = F,col.names=F)
write.table(proteins_diff_down,"proteins_diff_down.txt",sep="\t",quote=F,row.names = F,col.names=F)
write.table(proteins_shift,"proteins_shift.txt",sep="\t",quote=F,row.names = F,col.names=F)
write.table(proteins_noShift,"proteins_noShift.txt",sep="\t",quote=F,row.names = F,col.names=F)
write.table(proteins_noShiftNegative,"proteins_noShiftNegative.txt",sep="\t",quote=F,row.names = F,col.names=F)
write.table(proteins_shiftNegative,"proteins_shiftNegative.txt",sep="\t",quote=F,row.names = F,col.names=F)
write.table(proteins_noShiftPositive,"proteins_noShiftPositive.txt",sep="\t",quote=F,row.names = F,col.names=F)
write.table(proteins_shiftPositive,"proteins_shiftPositive.txt",sep="\t",quote=F,row.names = F,col.names=F)


proteins_global_up <- unique(prot_DiffExprProt[global_sumLog2FC_imp > 1]$feature_id)
proteins_global_down <- unique(prot_DiffExprProt[global_sumLog2FC_imp < -1]$feature_id)

pdf("diffProt_PRPF8.pdf",width=6,height=4)
for (id in "Q6P2Q9"){ #error for 1st ID????
  sub <- subset(proteinFeatures, protein_id==id)
  plotFeatures.tracesList(sub, pepTraces, id, design_matrix = design_matrix,
                        onlyBest = F, peak_area = T, legend=F)
}
dev.off()

pdf("diffProt_PKM.pdf",width=6,height=4)
for (id in "P14618"){ #error for 1st ID????
  sub <- subset(proteinFeatures, protein_id==id)
  plotFeatures.tracesList(sub, pepTraces, id, design_matrix = design_matrix,
                        onlyBest = F, peak_area = T, legend=F)
}
dev.off()


pdf("diffProt_sortByShift_imp.pdf",width=6,height=4)
for (id in unique(diffProt$feature_id)[1:30]){ #error for 1st ID????
  sub <- subset(proteinFeatures, protein_id==id)
  plotFeatures.tracesList(sub, pepTraces, id, design_matrix = design_matrix,
                        onlyBest = F, peak_area = T, legend=F)
}
dev.off()

pdf("diffProtShift.pdf",width=6,height=4)
for (id in unique(diffProt_shift$feature_id)[1:100]){ #error for 1st ID????
  sub <- subset(proteinFeatures, (protein_id==id) & (apex %in% diffProt_shift[feature_id==id]$apex))
  plotFeatures.tracesList(sub, pepTraces, id, design_matrix = design_matrix,
                        onlyBest = F, peak_area = T, legend=F)
}
dev.off()

pdf("diffProtNoShift.pdf",width=6,height=4)
for (id in unique(diffProt_noShift$feature_id)[1:100]){
  sub <- subset(proteinFeatures, (protein_id==id) & (apex %in% diffProt_noShift[feature_id==id]$apex))
  plotFeatures.tracesList(sub, pepTraces, id, design_matrix = design_matrix,
                          onlyBest = F, peak_area = T, legend=F)
}
dev.off()


example_prot <- c("P40763","Q9UJX2","Q12834",
  "Q14674","Q8WVK7",
  "Q9Y266","O43143","P07195")

exampleDiff <- subset(prot_DiffExprProt, feature_id %in% example_prot)

exampleDiff <- merge(exampleDiff,exampleTraceAnnotation, by.x="feature_id", by.y="Entry")

pdf("exampleDiff.pdf",width=6,height=4)
for (id in example_prot){
  sub <- subset(proteinFeatures, (protein_id==id))
  plotFeatures.tracesList(sub, pepTraces, id, design_matrix = design_matrix,
                          onlyBest = F, peak_area = T, legend=F)
}
dev.off()


"P40763""STAT3"

"Q9UJX2""APC8"

"Q12834""CDC20"

"Q14674""Separin"

"Q8IZT6" "ASPM"

 "Q8WVK7" "SKA2"

"Q9Y266""NUDC"

"O43143" "hPrp43"

"P07195" "LDHB"

"P42167" "LAP2"
