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

#install_github("CCprofiler/CCprofiler", ref = "proteoformLocationMapping")
#library(CCprofiler)
#setwd("~/mysonas/PRPF8/analysis/output")

#setwd("~/Desktop/CCprofiler/CCprofiler")
#load_all()
#setwd("~/Desktop/PRPF8data/")

#' ## Load data
peptide_traces_list <- readRDS("pepTraces_proteoformMulti.rds")
design_matrix <- readRDS("design_matrix.rda")

#' ## Remove peptides with no proteoform annotation
#peptides_annotated_minus <- peptide_traces_list$minus$trace_annotation[which(!is.na(peptide_traces_list$minus$trace_annotation$proteoform_id))]$id
#peptides_annotated_plus <- peptide_traces_list$plus$trace_annotation[which(!is.na(peptide_traces_list$plus$trace_annotation$proteoform_id))]$id
#peptide_traces_list_sub <- copy(peptide_traces_list)
#peptide_traces_list_sub$minus <- subset(peptide_traces_list_sub$minus, trace_subset_ids = peptides_annotated_minus)
#peptide_traces_list_sub$plus <- subset(peptide_traces_list_sub$plus, trace_subset_ids = peptides_annotated_plus)

#' ## Integrate traces across conditions
pepTracesSum <- integrateTraceIntensities(peptide_traces_list,
                                          design_matrix = NULL,
                                          integrate_within = NULL,
                                          aggr_fun = "sum")
saveRDS(pepTracesSum,"pepTracesSum.rda")

#' Subset traces for testing purposes
testProts <- c(unique(pepTracesSum$trace_annotation$protein_id)[1:50],"P22102","Q6P2Q9","P14618","P22102","O00468","Q9UBF2",
               "O15067","O75822","Q9Y266",
               "P42167","P61978","O95801",
               "P55060","O14578","O75717","P35221","Q9HB71","Q9H6S3")
subsetTest <- subset(pepTracesSum,trace_subset_ids=testProts,trace_subset_type = "protein_id")
saveRDS(subsetTest,"subsetTest.rda")

#' ## Proteoform-specific feature finding
proteinFeatures  <- findProteinFeatures(traces=subsetTest,
                                        corr_cutoff=0.9,
                                        window_size=7,
                                        parallelized=F,
                                        n_cores=1,
                                        collapse_method="apex_only",
                                        perturb_cutoff= "5%",
                                        rt_height=1,
                                        smoothing_length=7,
                                        useRandomDecoyModel=TRUE,
                                        quantLevel = "protein_id")

saveRDS(proteinFeatures, "proteinFeatures.rda")

resolveProteoformSpecificFeatures <- function(features,traces){
  idx <- seq_len(nrow(features))
  res <- lapply(idx, function(i){
    x <- features[i,]
    peptides_detected <- unlist(strsplit(x$subunits_detected, split = ";"))
    traces_sub <- subset(traces$trace_annotation,id %in% peptides_detected
    )
    proteoforms_detected <- sort(unique(traces_sub$proteoform_id))
    proteoform_subunits_annotated <- unique(subset(traces$trace_annotation, proteoform_id %in% proteoforms_detected)$id)
    completeness <- length(peptides_detected)/length(proteoform_subunits_annotated)
    report <- copy(x)
    report$subunits_annotated <- paste(proteoform_subunits_annotated,collapse = ";")
    report$n_subunits_annotated <- length(proteoform_subunits_annotated)
    report$subunits_detected <- paste(peptides_detected,collapse = ";")
    report$n_subunits_detected <- length(peptides_detected)
    report$completeness <- completeness
    report[,proteoform_ids := paste(proteoforms_detected, collapse = ";")]
    report[,n_proteoform_ids := length(proteoforms_detected)]
    return(report)
  })
  combiRes <- do.call(rbind,res)
  return(combiRes)
}

proteoformFeaturesResolved <- resolveProteoformSpecificFeatures(proteinFeatures,subsetTest)
saveRDS(proteoformFeaturesResolved,"proteoformFeaturesResolved.rda")

#' ## Proteoform feature scoring and filtering
filteredDataProteoformResolved <- scoreFeatures(proteoformFeaturesResolved, FDR=0.05, PDF=T, name="qvalueStats_proteoformFeaturesResolved")
saveRDS(filteredDataProteoformResolved, "filteredDataProteoformResolved.rda")

#' ## Reassignment of proteoforms
proteoforms_perProtein <- lapply(unique(filteredDataProteoformResolved$protein_id), function(x) {
  d <- subset(filteredDataProteoformResolved,protein_id==x)
  discrete_proteoforms <- unique(d$proteoform_ids)
  n_discrete_proteoforms <- length(unique(d$proteoform_ids))
  unique_peptide_groups <- sort(unique(unlist(strsplit(d$proteoform_ids, ';'))))
  n_unique_peptide_groups <- length(sort(unique(unlist(strsplit(d$proteoform_ids, ';')))))
  return(list(discrete_proteoforms=discrete_proteoforms,
  n_discrete_proteoforms=n_discrete_proteoforms,
  unique_peptide_groups=unique_peptide_groups,
  n_unique_peptide_groups=n_unique_peptide_groups))
  })
names(proteoforms_perProtein)<- unique(filteredDataProteoformResolved$protein_id)

unlist(lapply(names(proteoforms_perProtein),function(x){proteoforms_perProtein[[x]]$n_unique_peptide_groups}))
unlist(lapply(names(proteoforms_perProtein),function(x){proteoforms_perProtein[[x]]$n_discrete_proteoforms}))

#' ## Summary of detected proteoform features
summarizeFeatures(filteredDataProteoformResolved,plot=TRUE,PDF=TRUE,name="filteredDataProtein_summary_filteredDataProteoformResolved")

#' ## Plot some example features
calibrationFunctions <- readRDS("calibration.rds")

#filteredDataProteoform[,protein_id:=gsub("\\_.+","",proteoform_id)]
#filteredDataProteoform[,feature_id:=protein_id]

#filteredDataProteoform_P22102 <- subset(filteredDataProteoform,protein_id=="P22102")
pdf("features_proteoformResolved.pdf",width=4,height=4)
for (id in unique(filteredDataProteoformResolved$protein_id)){
  plotFeatures(feature_table = filteredDataProteoformResolved,
               traces = peptide_traces_list,
               calibration=calibrationFunctions,
               design_matrix = design_matrix,
               feature_id = id,
               annotation_label="proteoform_id",
               colour_by="proteoform_id",
               peak_area = TRUE,
               legend = FALSE,
               onlyBest = FALSE)
}
dev.off()


plotFeatures(feature_table = filteredDataProteoformResolved,
             traces = peptide_traces_list,
             calibration=calibrationFunctions,
             design_matrix = design_matrix,
             feature_id = "O00468",
             annotation_label="proteoform_id",
             colour_by="proteoform_id",
             peak_area = TRUE,
             legend = FALSE,
             onlyBest = FALSE)
filteredDataProteoformResolved[protein_id=="O00468"][order(apex)]

plotFeatures(feature_table = filteredDataProteoformResolved,
             traces = peptide_traces_list,
             calibration=calibrationFunctions,
             design_matrix = design_matrix,
             feature_id = "P22102",
             annotation_label="proteoform_id",
             colour_by="proteoform_id",
             peak_area = TRUE,
             legend = FALSE,
             onlyBest = FALSE)
filteredDataProteoformResolved[protein_id=="P22102"][order(apex)]

plotFeatures(feature_table = filteredDataProteoformResolved,
             traces = peptide_traces_list,
             calibration=calibrationFunctions,
             design_matrix = design_matrix,
             feature_id = "O75717",
             annotation_label="proteoform_id",
             colour_by="proteoform_id",
             peak_area = TRUE,
             legend = FALSE,
             onlyBest = FALSE)
filteredDataProteoformResolved[protein_id=="O75717"][order(apex)]


plotFeatures(feature_table = filteredDataProteoformResolved,
             traces = peptide_traces_list,
             calibration=calibrationFunctions,
             design_matrix = design_matrix,
             feature_id = "P42167",
             annotation_label="proteoform_id",
             colour_by="proteoform_id",
             peak_area = TRUE,
             legend = FALSE,
             onlyBest = FALSE)
filteredDataProteoformResolved[protein_id=="P42167"][order(apex)]



plotFeatures(feature_table = filteredDataProteoformResolved,
             traces = peptide_traces_list,
             calibration=calibrationFunctions,
             design_matrix = design_matrix,
             feature_id = "P61978",
             annotation_label="proteoform_id",
             colour_by="proteoform_id",
             peak_area = TRUE,
             legend = FALSE,
             onlyBest = FALSE)
filteredDataProteoformResolved[protein_id=="P61978"][order(apex)]


plotFeatures(feature_table = filteredDataProteoformResolved,
             traces = peptide_traces_list,
             calibration=calibrationFunctions,
             design_matrix = design_matrix,
             feature_id = "Q9H6S3",
             annotation_label="proteoform_id",
             colour_by="proteoform_id",
             peak_area = TRUE,
             legend = FALSE,
             onlyBest = FALSE)
filteredDataProteoformResolved[protein_id=="Q9H6S3"][order(apex)]


#' ## Extract intensity values for all detected features
#' Load raw peptide-level traces
pepTraces <- readRDS("pepTraces_proteoformMulti.rds")
#' ## Extract feature values
proteoform_featureVals <- extractFeatureVals(traces = pepTraces,
                                             features = filteredDataProteoformCollapsed,
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
#' @TODO change for proteoform traces
proteoform_DiffExprProtein <- aggregatePeptideTests(proteoform_DiffExprPep)
saveRDS(proteoform_DiffExprProtein, "proteoform_DiffExprProtein.rda")

proteoform_DiffExprProteoform <- aggregatePeptideTestsToProteoform(proteoform_DiffExprPep)
saveRDS(proteoform_DiffExprProteoform, "proteoform_DiffExprProteoform.rda")

#' Make volcano plots
#' General volcano plots
plotVolcano(proteoform_DiffExprPep, PDF = F, name = "proteoform_DiffExprPep")
plotVolcano(proteoform_DiffExprProteoform, PDF = F, name = "proteoform_DiffExprProteoform")
plotVolcano(proteoform_DiffExprProtein, PDF = F, name = "proteoform_DiffExprProtein")

#' Volcanoplts highlighting different proteoforms
library(ggrepel)

plotVolcano(proteoform_DiffExprProtein, highlight=c("Q6P2Q9"), PDF = F, name = "prot_DiffExprProt_PRPF8_Q6P2Q9")

plotVolcano(proteoform_DiffExprProteoform, highlight=c("Q9HB71"), PDF = T, name = "prot_DiffExprProt_CACYBP_Q9HB71")

diff <- proteoform_DiffExprProtein[(abs(sumLog2FC)>=1) & (pBHadj <= 0.05)]
nrow(proteoform_DiffExprProtein[(abs(sumLog2FC)>=1) & (pBHadj <= 0.05)])
length(unique(proteoform_DiffExprProtein[(abs(sumLog2FC)>=1) & (pBHadj <= 0.05)]$feature_id))

diff[abs(local_vs_global_log2FC_imp) > 1]

#plotVolcano(proteoform_DiffExprProteoform, highlight=c("Q6P2Q9_1"), PDF = T, name = "prot_DiffExprProt_PRPF8_Q6P2Q9_1")
#plotVolcano(proteoform_DiffExprProteoform, highlight=c("Q6P2Q9_2"), PDF = T, name = "prot_DiffExprProt_PRPF8_Q6P2Q9_2")
#plotVolcano(proteoform_DiffExprProteoform, highlight=c("Q6P2Q9_3"), PDF = T, name = "prot_DiffExprProt_PRPF8_Q6P2Q9_3")

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

proteoform_stats <- proteoform_DiffExprProteoform[feature_id=="Q9HB71"]
globalLog2FC <- median(proteoform_stats$global_sumLog2FC_imp)
setnames(proteoform_stats,c("qint1","qint2"),c("minus","plus"))
proteoform_stats <- subset(proteoform_stats, select=c("feature_id","apex","sumLog2FC","pBHadj"))
proteoform_stats.m <- melt(proteoform_stats,id.vars=c("feature_id","apex","pBHadj"))
setnames(proteoform_stats.m,"value","log2FC")
setnames(proteoform_stats.m,"apex","fraction")

pdf(paste0("proteoform_stats",".pdf"), height=2, width=3)
ggplot(proteoform_stats.m,aes(x=fraction,y=log2FC)) +
  geom_bar(stat="identity",fill="steelblue") +
  theme_classic() +
  scale_x_continuous(breaks=seq(0, max(pepTraces$minus$fraction_annotation$id), 20),
                     limits = c(0, max(pepTraces$minus$fraction_annotation$id))) +
  geom_hline(yintercept = globalLog2FC, linetype="dashed", colour="#E69F00") +
  geom_hline(yintercept = 0, linetype="solid")
dev.off()


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
