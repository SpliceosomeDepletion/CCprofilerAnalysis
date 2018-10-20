library(devtools)
setwd("~/mysonas/CCprofiler")
load_all()
setwd("~/mysonas/PRPF8/analysis/output")
#install_github("CCprofiler/CCprofiler", ref = "proteoformLocationMapping")
#library(CCprofiler)

############################
############################
############################

traces <- readRDS("pepTraces_proteoformMulti.rds")
design_matrix <- readRDS("design_matrix.rda")

############################
############################
############################

# proteoform feature finding
pepTracesSum <- integrateTraceIntensities(traces,
                                       design_matrix = NULL,
                                       integrate_within = NULL,
                                       aggr_fun = "sum")

testProts <- unique(traces$minus$trace_annotation$protein_id)[1:10]
subsetTest <- subset(pepTracesSum,trace_subset_ids=testProts,trace_subset_type = "protein_id")

proteinFeatures  <- findProteinFeatures(traces=subsetTest,
                                            corr_cutoff=0.9,
                                            window_size=8,
                                            parallelized=F,
                                            n_cores=1,
                                            collapse_method="apex_only",
                                            perturb_cutoff= "5%",
                                            rt_height=1,
                                            smoothing_length=7,
                                            useRandomDecoyModel=TRUE)

proteoformFeatures  <- findProteinFeatures(traces=subsetTest,
                                            corr_cutoff=0.9,
                                            window_size=8,
                                            parallelized=F,
                                            n_cores=1,
                                            collapse_method="apex_only",
                                            perturb_cutoff= "5%",
                                            rt_height=1,
                                            smoothing_length=7,
                                            useRandomDecoyModel=TRUE,
                                            quantLevel = "proteoform_id")


filteredDataProtein <- scoreFeatures(proteinFeatures, FDR=0.05, PDF=T, name="qvalueStats_proteinFeatures_integrated")
filteredDataProteoform <- scoreFeatures(proteoformFeatures, FDR=0.05, PDF=T, name="qvalueStats_proteoformFeatures_integrated")


summarizeFeatures(filteredDataProtein,plot=TRUE,PDF=TRUE,name="filteredDataProtein_summary_integrated")
summarizeFeatures(filteredDataProteoform,plot=TRUE,PDF=TRUE,name="filteredDataProtein_summary_integrated")

pdf("features_O75822.pdf")
plotFeatures(feature_table = filteredDataProtein,
             traces = pepTracesSum,
             feature_id = "O75822",
             peak_area = TRUE,
             legend = FALSE,
             onlyBest = FALSE)
dev.off()

pdf("features_O75822_1.pdf")
plotFeatures(feature_table = filteredDataProteoform,
             traces = pepTracesSum,
             feature_id = "O75822_1",
             peak_area = TRUE,
             legend = FALSE,
             onlyBest = FALSE)
dev.off()

pdf("features_O75822_2.pdf")
plotFeatures(feature_table = filteredDataProteoform,
             traces = pepTracesSum,
             feature_id = "O75822_2",
             peak_area = TRUE,
             legend = FALSE,
             onlyBest = FALSE)
dev.off()

#Q9Y2Z0

pdf("features_Q9Y2Z0.pdf")
plotFeatures(feature_table = filteredDataProtein,
             traces = pepTracesSum,
             feature_id = "Q9Y2Z0",
             peak_area = TRUE,
             legend = FALSE,
             onlyBest = FALSE)
dev.off()

pdf("features_Q9Y2Z0_1.pdf")
plotFeatures(feature_table = filteredDataProteoform,
             traces = pepTracesSum,
             feature_id = "Q9Y2Z0_1",
             peak_area = TRUE,
             legend = FALSE,
             onlyBest = FALSE)
dev.off()

pdf("features_Q9Y2Z0_2.pdf")
plotFeatures(feature_table = filteredDataProteoform,
             traces = pepTracesSum,
             feature_id = "Q9Y2Z0_2",
             peak_area = TRUE,
             legend = FALSE,
             onlyBest = FALSE)
dev.off()


pepTraces <- readRDS("pepTracesRaw.rda")

prot_featureVals <- extractFeatureVals(traces = pepTraces,
                                  features = filteredDataProtein,
                                  design_matrix = design_matrix,
                                  extract = "subunits_detected",
                                  imputeZero = T,
                                  verbose = F,
                                  perturb_cutoff = "5%")

proteoform_featureVals <- extractFeatureVals(traces = pepTraces,
                                  features = filteredDataProteoform,
                                  design_matrix = design_matrix,
                                  extract = "subunits_detected",
                                  imputeZero = T,
                                  verbose = F,
                                  perturb_cutoff = "5%")

prot_featureValsFilled <- fillFeatureVals(featureVals = prot_featureVals,
                                     design_matrix = design_matrix)

proteoform_featureValsFilled <- fillFeatureVals(featureVals = proteoform_featureVals,
                                    design_matrix = design_matrix)

prot_DiffExprPep <- testDifferentialExpression(featureVals = prot_featureValsFilled,
                          compare_between = "Condition",
                          level = "peptide",
                          measuredOnly = FALSE)

proteoform_DiffExprPep <- testDifferentialExpression(featureVals = proteoform_featureValsFilled,
                          compare_between = "Condition",
                          level = "peptide",
                          measuredOnly = FALSE)

proteoform_DiffExprPep[feature_id=="O75822_1" & apex==30]
proteoform_DiffExprPep[feature_id=="O75822_1" & apex==50]

for (test_protein in c("O75822","P68400")){
  protTest <- subset(pepTraces, trace_subset_ids = test_protein, trace_subset_type = "protein_id")
  plot(protTest, design_matrix = design_matrix, log = FALSE, legend = FALSE, PDF = TRUE,
   name = paste0("PeptideTraces (",paste(test_protein,collapse="_"),")"), plot = TRUE, highlight = NULL, highlight_col = NULL)
}

############################
############################
############################

proteoform_available_minus <- unique(traces$minus$trace_annotation$proteoform_id)
proteoform_available_minus <- proteoform_available_minus[!is.na(proteoform_available_minus)]
proteoform_available_plus <- unique(traces$plus$trace_annotation$proteoform_id)
proteoform_available_plus <- proteoform_available_plus[!is.na(proteoform_available_plus)]
proteoform_available <- unique(proteoform_available_minus,proteoform_available_plus)

traces <- subset(traces,trace_subset_ids=proteoform_available,trace_subset_type = "proteoform_id")

proteoformTraces <- proteinQuantification(traces,quantLevel="proteoform_id",
  topN = 1000,
  keep_less = TRUE,
  rm_decoys = TRUE,
  use_sibPepCorr = FALSE,
  use_repPepCorr = FALSE,
  full_intersect_only = FALSE)

proteoformTraces <- updateTraces(proteoformTraces)

proteoformTracesSum <- integrateTraceIntensities(proteoformTraces,
                                       design_matrix = NULL,
                                       integrate_within = NULL,
                                       aggr_fun = "sum")

saveRDS(proteoformTraces,"proteoformTraces.rds")
saveRDS(proteoformTracesSum,"proteoformTracesSum.rds")


## plot some example proteins
test_proteins <- unique(proteoformTraces$minus$trace_annotation[proteoform_pval <= 0.05]$protein_id)[1:30]
#test_proteins <- c("O14578","Q9UBF2","O00468","P22102","Q9Y266","P42167","P61978") ## "P42166",
for (test_protein in test_proteins){
  protTest <- subset(proteoformTraces, trace_subset_ids = test_protein, trace_subset_type = "protein_id")
  plot(protTest, design_matrix = design_matrix, log = FALSE, legend = TRUE, PDF = TRUE,
   name = paste0("ProteoformTraces (",paste(test_protein,collapse="_"),")"), plot = TRUE, highlight = NULL, highlight_col = NULL)
}
#awesome example: Q9Y6E0 >> changes association upon spliceosome inhibition
# another example: Q9UPN9 >. but no change

# evaluate protein mass distribution
proteoformTraces_assembly <- annotateMassDistribution(proteoformTraces)
saveRDS(proteoformTraces_assembly,"proteoformTraces_assembly.rds")

# evaluate differential protein mass distribution
proteoformDiffAssemblyState <- getMassAssemblyChange(proteoformTraces_assembly, design_matrix, quantLevel="proteoform_id")
saveRDS(proteoformDiffAssemblyState,"proteoformDiffAssemblyState.rds")

change_cutoff=0.2
plotMassAssemblyChange(proteoformDiffAssemblyState, change_cutoff=change_cutoff, PDF=T,name="proteoformDiffAssemblyState")

assembly_info <- copy(proteoformDiffAssemblyState)
assembly_info[,protein_id := gsub("\\_.*","",proteoform_id)]
assembly_info[,min_change_perProt := min(change), by=protein_id]
assembly_info[,max_change_perProt := max(change), by=protein_id]
assembly_info[,altProteoforms := ifelse((min_change_perProt < -0.1) & (max_change_perProt > 0.1),1,0)]


interesting_proteins <- unique(assembly_info[altProteoforms==1]$protein_id)
for (test_protein in interesting_proteins){
  protTest <- subset(proteoformTraces, trace_subset_ids = test_protein, trace_subset_type = "protein_id")
  plot(protTest, design_matrix = design_matrix, log = FALSE, legend = TRUE, PDF = TRUE,
   name = paste0("ProteoformTraces (",paste(test_protein,collapse="_"),")"), plot = TRUE, highlight = NULL, highlight_col = NULL)
}
