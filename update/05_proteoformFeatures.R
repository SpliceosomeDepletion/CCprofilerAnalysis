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
peptide_traces_list <- readRDS("pepTraces_proteoformMulti.rds")
design_matrix <- readRDS("design_matrix.rda")

#' ## Remove peptides with no proteoform annotation
peptides_annotated_minus <- peptide_traces_list$minus$trace_annotation[which(!is.na(peptide_traces_list$minus$trace_annotation$proteoform_id))]$id
peptides_annotated_plus <- peptide_traces_list$plus$trace_annotation[which(!is.na(peptide_traces_list$plus$trace_annotation$proteoform_id))]$id
peptide_traces_list_sub <- copy(peptide_traces_list)
peptide_traces_list_sub$minus <- subset(peptide_traces_list_sub$minus, trace_subset_ids = peptides_annotated_minus)
peptide_traces_list_sub$plus <- subset(peptide_traces_list_sub$plus, trace_subset_ids = peptides_annotated_plus)

#' ## Integrate traces across conditions
pepTracesSum <- integrateTraceIntensities(peptide_traces_list_sub,
                                          design_matrix = NULL,
                                          integrate_within = NULL,
                                          aggr_fun = "sum")
# Subset traces for testing purposes
testProts <- c(unique(pepTracesSum$trace_annotation$protein_id)[1:10],"P22102","Q6P2Q9","P14618")
subsetTest <- subset(pepTracesSum,trace_subset_ids=testProts,trace_subset_type = "protein_id")

#' ## Proteoform-specific feature finding
proteoformFeatures  <- findProteinFeatures(traces=subsetTest,
                                           corr_cutoff=0.9,
                                           window_size=8,
                                           parallelized=F,
                                           n_cores=1,
                                           collapse_method="apex_only",
                                           perturb_cutoff= "5%",
                                           rt_height=3,
                                           smoothing_length=7,
                                           useRandomDecoyModel=TRUE,
                                           quantLevel = "proteoform_id")

#' ## Proteoform feature scoring and filtering
filteredDataProteoform <- scoreFeatures(proteoformFeatures, FDR=0.1, PDF=T, name="qvalueStats_proteoformFeatures_integrated")

#' ## Summary of detected proteoform features 
summarizeFeatures(filteredDataProteoform,plot=TRUE,PDF=TRUE,name="filteredDataProtein_summary_integrated")

#' ## Save data
saveRDS(proteoformFeatures, "proteoformFeatures.rda")
saveRDS(filteredDataProteoform, "filteredDataProteoform.rda")

#' ## Plot some example features
pdf("features_P22102_2.pdf")
plotFeatures(feature_table = filteredDataProteoform,
             traces = pepTracesSum,
             feature_id = "P22102_2",
             peak_area = TRUE,
             legend = FALSE,
             onlyBest = FALSE)
dev.off()

pdf("features_P22102_1.pdf")
plotFeatures(feature_table = filteredDataProteoform,
             traces = pepTracesSum,
             feature_id = "P22102_1",
             peak_area = TRUE,
             legend = FALSE,
             onlyBest = FALSE)
dev.off()

#' ## Extract intensity values for all detected features
#' Load raw peptide-level traces
pepTraces <- readRDS("pepTraces_proteoformMulti.rds")
#' ## Extract feature values
proteoform_featureVals <- extractFeatureVals(traces = pepTraces,
                                             features = filteredDataProteoform,
                                             design_matrix = design_matrix,
                                             extract = "subunits_detected",
                                             imputeZero = T,
                                             verbose = F,
                                             perturb_cutoff = "5%")
#' ## Fill feature values 
proteoform_featureValsFilled <- fillFeatureVals(featureVals = proteoform_featureVals,
                                                design_matrix = design_matrix)

#' ## Perform peptide-level differential expression testing for all features 
proteoform_DiffExprPep <- testDifferentialExpression(featureVals = proteoform_featureValsFilled,
                                                     compare_between = "Condition",
                                                     level = "peptide",
                                                     measuredOnly = FALSE)

saveRDS(proteoform_DiffExprPep, "proteoform_DiffExprPep.rda")

#' ## Aggregate differential expression results to the proteoform level
proteoform_DiffExprProteoform <- aggregatePeptideTests(proteoform_DiffExprPep)

saveRDS(proteoform_DiffExprProteoform, "proteoform_DiffExprProteoform.rda")

#' Make volcano plots
#' General volcano plots
plotVolcano(proteoform_DiffExprPep, PDF = T, name = "proteoform_DiffExprPep")
plotVolcano(proteoform_DiffExprProteoform, PDF = T, name = "proteoform_DiffExprProteoform")
#' Volcanoplts highlighting different proteoforms
plotVolcano(proteoform_DiffExprProteoform, highlight=c("Q6P2Q9_1"), PDF = T, name = "prot_DiffExprProt_PRPF8_Q6P2Q9_1")
plotVolcano(proteoform_DiffExprProteoform, highlight=c("Q6P2Q9_2"), PDF = T, name = "prot_DiffExprProt_PRPF8_Q6P2Q9_2")
plotVolcano(proteoform_DiffExprProteoform, highlight=c("Q6P2Q9_3"), PDF = T, name = "prot_DiffExprProt_PRPF8_Q6P2Q9_3")

###########################
###########################
###########################
###########################
###########################

#' ## Combine proteoform features within one protein
#' @TODO 
#' @NEW merge features with max rt diff of +/-1 
#' @NEW aggregte tests results across merged features
#' @NEW Differentiate common and proteoform specific features
#' @NEW Find examples where different proteoforms change in opposing directions
#' 
proteoform_DiffExprProteoform[,protein_id:=gsub("\\_.+","",feature_id)]


getUniqueFeatures <- function(featue_set,rt_height=3,method="complete"){
  if (nrow(featue_set) > 1) {
    apex_dist <- dist(featue_set$apex)
    apex_clust <- hclust(apex_dist,method=method)
    apex_groups <- cutree(apex_clust, h=rt_height+0.01)
    # add 0.01 to rt_height to cut tree right above cutoff (2 features with exact dist rt_height are merged)
    featue_set[,uniqueApex:=apex_groups]
    return(featue_set)
  } else {
    featue_set[,uniqueApex:=1]
    return(featue_set)
  }
}

getUniqueFeatures_list <- function(prot,x){
  y <- subset(x,protein_id==prot)
  getUniqueFeatures(y)}

prots <- unique(proteoform_DiffExprProteoform$protein_id)

res <- lapply(prots, getUniqueFeatures_list, x=proteoform_DiffExprProteoform)
res.all <- do.call(rbind,res)

#' @TODO aggregate properly by function in CCpofiler
res.all[,sumLog2FC_uniqueApex_mean := mean(sumLog2FC),by=c("protein_id","uniqueApex")]
res.all[,pBHadj_uniqueApex_mean := mean(pBHadj),by=c("protein_id","uniqueApex")]
res.all[,local_vs_global_log2FC_imp_uniqueApex_mean := mean(local_vs_global_log2FC_imp),by=c("protein_id","uniqueApex")]
res.all[,apex_mean := round(mean(apex)),by=c("protein_id","uniqueApex")]
res.all[,proteoforms := paste(sort(feature_id),collapse=";"), by=c("protein_id","uniqueApex")]

res_sub <- subset(res.all,select=c("protein_id","proteoforms","apex_mean",
                                   "sumLog2FC_uniqueApex_mean",
                                   "pBHadj_uniqueApex_mean",
                                   "local_vs_global_log2FC_imp_uniqueApex_mean"))
res_sub <- unique(res_sub)

setnames(res_sub,c("sumLog2FC_uniqueApex_mean","pBHadj_uniqueApex_mean","protein_id"),c("sumLog2FC","pBHadj","feature_id"))

plotVolcano(res_sub, PDF = F, name = "proteoform_DiffExprProteoform")
plotVolcano(res_sub, PDF = F, name = "proteoform_DiffExprProteoform",highlight="Q6P2Q9")
plotVolcano(res_sub, PDF = F, name = "proteoform_DiffExprProteoform",highlight="P22102")

###########################
###########################
###########################
###########################
###########################

#' @TODO 
#' @NEW Trick diff exp. function to calculate global FC and p-vals

# on peptide, proteoform and protein level!
# Create filteredDataProteoform for all peptides (names(filteredDataProteoform))
#1) Create list of all peptide-proteoform assignments across conditions
#2) subunits annotated, detected ... = all peptides
#test[,"subunits_annotated":=...]
pepTraces <- readRDS("pepTraces_proteoformMulti.rds")
testProts <- c(unique(pepTraces$minus$trace_annotation$protein_id)[1:50],"P22102","Q6P2Q9","P14618")
pepTraces <- subset(pepTraces,trace_subset_ids=testProts,trace_subset_type = "protein_id")

globalDiffExp <- testGlobalDifferentialExpression(pepTraces,design_matrix=design_matrix)


total_fractions.dt <- lapply(pepTraces, "[[", 4)
total_fractions <- max(unique(unlist(lapply(total_fractions.dt, function(x){x$id}))))

ann <- lapply(pepTraces, "[[", 3)
ann_all <- unique(do.call(rbind,lapply(ann, function(x){subset(x,select=c("id","protein_id"))})))
ann_all[,protein_name := protein_id]
ann_all[,subunits_annotated := paste(id, collapse=';'), by="protein_id"]
mock_features <- unique(ann_all[,id:=NULL])
mock_features[,n_subunits_annotated:=.N, by="protein_id"]
mock_features[,subunits_detected:=subunits_annotated]
mock_features[,n_subunits_detected:=n_subunits_annotated]
mock_features[,completeness:=1]
mock_features[,left_pp:=1]
mock_features[,right_pp:=total_fractions]
mock_features[,apex:=1]
mock_features[,apex_mw:=1]
mock_features[,area:=1]
mock_features[,peak_corr:=1]
mock_features[,monomer_sec:=1]
mock_features[,monomer_mw:=1]
mock_features[,in_complex:=FALSE]
mock_features[,coelution_score:=1]
mock_features[,decoy:=1]
mock_features[,pvalue:=1]
mock_features[,qvalue:=1]

mock_features_featureVals <- extractFeatureVals(traces = pepTraces,
                                             features = mock_features,
                                             design_matrix = design_matrix,
                                             extract = "subunits_detected",
                                             imputeZero = T,
                                             verbose = F,
                                             perturb_cutoff = "5%")
#' ## Fill feature values 
mock_features_featureValsFilled <- fillFeatureVals(featureVals = mock_features_featureVals,
                                                design_matrix = design_matrix)

mock_features_DiffExprPep <- testDifferentialExpression(featureVals = mock_features_featureValsFilled,
                                                     compare_between = "Condition",
                                                     level = "peptide",
                                                     measuredOnly = FALSE)

mock_features_DiffExprProtein <- aggregatePeptideTests(mock_features_DiffExprPep)

mock_features_DiffExprProteoform <- aggregatePeptideTestsToProteoform(mock_features_DiffExprPep)

# Important to use global values here and not local!!! Local uses only fractions where there are values in both conditions!


ids <- mock_features_DiffExprProtein[pBHadj<=0.05 & abs(global_sumLog2FC_imp) > 1]$feature_id

for (test_protein in ids){
  protTest <- subset(pepTraces, trace_subset_ids = test_protein, trace_subset_type = "protein_id")
  plot(protTest, design_matrix = design_matrix, log = FALSE, legend = FALSE, PDF = F,
       name = paste0("PeptideTraces (",paste(test_protein,collapse="_"),")"), plot = TRUE, highlight = NULL, highlight_col = NULL)
}


###########################
###########################
###########################
###########################
###########################

#' @TODO 
#' @NEW traces plotting in general
#' @NEW feature plotting pep, prot, complex & proteoform level
#' @NEW plot traces with exon or proteoform colouring
#' @NEW plot traces with exon annotation >> Max trial version

###########################
###########################
###########################
###########################
###########################

proteoform_DiffExprProteoform[grep("Q6P2Q9",feature_id)]



proteoform_DiffExprPep[feature_id=="P55036_1" & apex==15]
proteoform_DiffExprPep[feature_id=="P55036_1" & apex==42]


for (test_protein in c("P55036","A5YKK6","Q13310","P22102","P35221","Q9BQS8")){
  protTest <- subset(pepTraces, trace_subset_ids = test_protein, trace_subset_type = "protein_id")
  plot(protTest, design_matrix = design_matrix, log = FALSE, legend = FALSE, PDF = TRUE,
       name = paste0("PeptideTraces (",paste(test_protein,collapse="_"),")"), plot = TRUE, highlight = NULL, highlight_col = NULL)
}

proteoform_DiffExprPep[feature_id=="A5YKK6_1"]


