#' ## Load data
traces_list_reassignedProteoforms <- readRDS("traces_list_reassignedProteoforms.rds")
reassignedProteoformFeatures <- readRDS("reassignedProteoformFeatures.rds")
design_matrix <- readRDS("design_matrix.rda")

#' ## Extract feature values
proteoform_featureVals <- extractFeatureVals(traces = traces_list_reassignedProteoforms,
                                             features = reassignedProteoformFeatures,
                                             design_matrix = design_matrix,
                                             extract = "subunits_detected",
                                             imputeZero = T,
                                             verbose = F)
saveRDS(proteoform_featureVals, "proteoform_featureVals.rda")

#' ## Fill feature values
proteoform_featureValsFilled <- fillFeatureVals(featureVals = proteoform_featureVals,
                                                tracesList = traces_list_reassignedProteoforms,
                                                design_matrix = design_matrix)
saveRDS(proteoform_featureValsFilled, "proteoform_featureValsFilled.rda")

#' ## Perform peptide-level differential expression testing for all features
proteoform_DiffExprPep <- testDifferentialExpression(featureVals = proteoform_featureValsFilled,
                                                     compare_between = "Condition",
                                                     level = "peptide",
                                                     measuredOnly = FALSE)

#' ## Aggregate differential expression results to the proteoform level
proteoform_DiffExprProteoform <- aggregatePeptideTestsToProteoform(proteoform_DiffExprPep)

#' ## Aggregate differential expression results to the protein level
proteoform_DiffExprProtein <- aggregatePeptideTests(proteoform_DiffExprPep)

#' Change sign of fold-change to accomodate the control as the reference (reverse to how test was calculated)
proteoform_DiffExprPep[,log2FC:=-log2FC]
proteoform_DiffExprPep[,global_log2FC_imp:=-global_log2FC_imp]
saveRDS(proteoform_DiffExprPep, "proteoform_DiffExprPep.rda")

proteoform_DiffExprProteoform[,medianLog2FC:=-medianLog2FC]
proteoform_DiffExprProteoform[,global_medianLog2FC:=-global_medianLog2FC]
proteoform_DiffExprProteoform[,global_medianLog2FC_imp:=-global_medianLog2FC_imp]
saveRDS(proteoform_DiffExprProteoform, "proteoform_DiffExprProteoform.rda")

proteoform_DiffExprProtein[,medianLog2FC:=-medianLog2FC]
proteoform_DiffExprProtein[,global_medianLog2FC:=-global_medianLog2FC]
proteoform_DiffExprProtein[,global_medianLog2FC_imp:=-global_medianLog2FC_imp]
saveRDS(proteoform_DiffExprProtein, "proteoform_DiffExprProtein.rda")

#' Make volcano plots
#' General volcano plots
plotVolcano(proteoform_DiffExprPep, PDF = T, name = "proteoform_DiffExprPep")
plotVolcano(proteoform_DiffExprProteoform, PDF = T, name = "proteoform_DiffExprProteoform")
plotVolcano(proteoform_DiffExprProtein, PDF = T, name = "proteoform_DiffExprProtein")

#' Volcanoplts highlighting different proteoforms
library(ggrepel)

plotVolcano(proteoform_DiffExprProtein, highlight=c("Q6P2Q9"), PDF = T, name = "prot_DiffExprProt_PRPF8_Q6P2Q9")

plotVolcano(proteoform_DiffExprProteoform, highlight=c("Q9HB71"), PDF = T, name = "prot_DiffExprProt_CACYBP_Q9HB71")

diff <- proteoform_DiffExprProtein[(abs(medianLog2FC)>=1) & (pBHadj <= 0.05)]
nrow(proteoform_DiffExprProtein[(abs(medianLog2FC)>=1) & (pBHadj <= 0.05)])
length(unique(proteoform_DiffExprProtein[(abs(medianLog2FC)>=1) & (pBHadj <= 0.05)]$feature_id))

###########################
###########################
###########################
###########################
###########################

#proteoform_stats <- proteoform_DiffExprProtein[feature_id=="Q9HB71"]
#globalLog2FC <- median(proteoform_stats$global_sumLog2FC_imp)
#setnames(proteoform_stats,c("qint1","qint2"),c("minus","plus"))
#proteoform_stats <- subset(proteoform_stats, select=c("feature_id","apex","sumLog2FC","pBHadj"))
#proteoform_stats.m <- melt(proteoform_stats,id.vars=c("feature_id","apex","pBHadj"))
#setnames(proteoform_stats.m,"value","log2FC")
#setnames(proteoform_stats.m,"apex","fraction")

#pdf(paste0("proteoform_stats",".pdf"), height=2, width=3)
#ggplot(proteoform_stats.m,aes(x=fraction,y=log2FC)) +
#  geom_bar(stat="identity",fill="steelblue") +
#  theme_classic() +
#  scale_x_continuous(breaks=seq(0, max(pepTraces$minus$fraction_annotation$id), 20),
#                     limits = c(0, max(pepTraces$minus$fraction_annotation$id))) +
#  geom_hline(yintercept = globalLog2FC, linetype="dashed", colour="#E69F00") +
#  geom_hline(yintercept = 0, linetype="solid")
#dev.off()


###########################
###########################
###########################
###########################
###########################

#proteoform_DiffExprProteoform[grep("Q6P2Q9",feature_id)]

#proteoform_DiffExprPep[feature_id=="P55036_1" & apex==15]
#proteoform_DiffExprPep[feature_id=="P55036_1" & apex==42]


#for (test_protein in c("P55036","A5YKK6","Q13310","P22102","P35221","Q9BQS8")){
#  protTest <- subset(pepTraces, trace_subset_ids = test_protein, trace_subset_type = "protein_id")
#  plot(protTest, design_matrix = design_matrix, log = FALSE, legend = FALSE, PDF = TRUE,
#       name = paste0("PeptideTraces (",paste(test_protein,collapse="_"),")"), plot = TRUE, highlight = NULL, highlight_col = NULL)
#}

#proteoform_DiffExprPep[feature_id=="A5YKK6_1"]
