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
pepTracesSum <- readRDS("pepTracesSum.rda")
proteinFeatures <- readRDS("proteinFeatures.rda")
calibrationFunctions <- readRDS("calibration.rds")

#' Source function
source("../CCprofilerAnalysis/benchmark/resolveProteoformSpecificFeatures.R")

proteoformFeaturesResolved <- resolveProteoformSpecificFeatures(
    features=proteinFeatures,
    traces=pepTracesSum,
    minProteoformIntensityRatio=0.1,
    perturb_cutoff="5%")

saveRDS(proteoformFeaturesResolved,"proteoformFeaturesResolved.rds")

filteredDataProteoformResolved <- scoreFeatures(
  proteoformFeaturesResolved,
  FDR=0.1, PDF=T,
  name=paste0("qvalueStats_proteoformFeaturesResolved"))

saveRDS(filteredDataProteoformResolved,"filteredDataProteoformResolved.rds")

pdf(paste0("filteredDataProteoformResolved_conditionSum.pdf"),width=4,height=4)
  for (id in unique(filteredDataProteoformResolved$protein_id)[1:50]){
    plotFeatures(feature_table = filteredDataProteoformResolved,
                 traces = pepTracesSum,
                 calibration=calibrationFunctions,
                 feature_id = id,
                 annotation_label="proteoform_id",
                 colour_by="proteoform_id",
                 peak_area = T,
                 legend = F,
                 onlyBest = F)
  }
  dev.off()

pdf(paste0("filteredDataProteoformResolved.pdf"),width=4,height=4)
  for (id in unique(filteredDataProteoformResolved$protein_id)[1:50]){
    plotFeatures(feature_table = filteredDataProteoformResolved,
                 traces = peptide_traces_list,
                 design_matrix=design_matrix,
                 calibration=calibrationFunctions,
                 feature_id = id,
                 annotation_label="proteoform_id",
                 colour_by="proteoform_id",
                 peak_area = T,
                 legend = F,
                 onlyBest = F)
  }
  dev.off()


filteredDataProteoformResolved_single <- filteredDataProteoformResolved[n_proteoform_ids==1]
pdf(paste0("filteredDataProteoformResolved_single.pdf"),width=4,height=4)
  for (id in unique(filteredDataProteoformResolved_single$protein_id)[1:50]){
    plotFeatures(feature_table = filteredDataProteoformResolved_single,
                 traces = peptide_traces_list,
                 design_matrix=design_matrix,
                 calibration=calibrationFunctions,
                 feature_id = id,
                 annotation_label="proteoform_id",
                 colour_by="proteoform_id",
                 peak_area = T,
                 legend = F,
                 onlyBest = F)
  }
  dev.off()


filteredDataProteoformResolved_multi <- filteredDataProteoformResolved[n_proteoform_ids>1]
pdf(paste0("filteredDataProteoformResolved_multi.pdf"),width=4,height=4)
  for (id in unique(filteredDataProteoformResolved_multi$protein_id)[1:50]){
    plotFeatures(feature_table = filteredDataProteoformResolved_multi,
                 traces = peptide_traces_list,
                 design_matrix=design_matrix,
                 calibration=calibrationFunctions,
                 feature_id = id,
                 annotation_label="proteoform_id",
                 colour_by="proteoform_id",
                 peak_area = T,
                 legend = F,
                 onlyBest = F)
  }
  dev.off()
