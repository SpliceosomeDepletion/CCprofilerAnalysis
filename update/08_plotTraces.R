#' ---
#' title: Plot traces
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

#setwd("~/Desktop/CCprofiler/CCprofiler")
#load_all()
#setwd("~/Desktop/PRPF8data/")


#' ## Load data
design_matrix <- readRDS("design_matrix.rda")
pepTraces <- readRDS("pepTraces_proteoformMulti.rds")
#testProts <- c(unique(pepTraces$minus$trace_annotation$protein_id)[1:200],"P22102","Q6P2Q9","P14618")
#pepTraces <- subset(pepTraces,trace_subset_ids=testProts,trace_subset_type = "protein_id")

testTraces <- subset(pepTraces,trace_subset_ids="Q6P2Q9",trace_subset_type = "protein_id")
plot(testTraces,legend = T,design_matrix = design_matrix)
plot(testTraces,legend = T, colour_by="proteoform_id",design_matrix = design_matrix)
plot(testTraces,legend = T, colour_by="protein_id",design_matrix = design_matrix)
plot(testTraces,legend = T, colour_by="isoform_id",design_matrix = design_matrix)

testTraces <- subset(pepTraces,trace_subset_ids="O75717",trace_subset_type = "protein_id")
plot(testTraces,legend = T,design_matrix = design_matrix)
plot(testTraces,legend = T,design_matrix = design_matrix, collapse_conditions = T)
plot(testTraces,legend = T,design_matrix = design_matrix,highlight="IWEDLDDDDPK")
plot(testTraces,legend = T, colour_by="proteoform_id",design_matrix = design_matrix,PDF=T,name="O75717_pepTraces_colouredByProteoform")
plot(testTraces,legend = T, colour_by="proteoform_id",design_matrix = design_matrix, collapse_conditions = T)
plot(testTraces,legend = T, colour_by="protein_id",design_matrix = design_matrix,PDF=T,name="O75717_pepTraces_colouredByProtein")
plot(testTraces,legend = F, colour_by="isoform_id",design_matrix = design_matrix)


####

proteoformTraces <- readRDS("proteoform_traces_list.rds")
testTraces <- subset(proteoformTraces,trace_subset_ids="O75717",trace_subset_type = "protein_id")
plot(testTraces,legend = F,design_matrix = design_matrix)
plot(testTraces,legend = T, colour_by="proteoform_id",design_matrix = design_matrix,highlight="O75717_1")
plot(testTraces,legend = T, colour_by="protein_id",design_matrix = design_matrix, collapse_conditions = T)

