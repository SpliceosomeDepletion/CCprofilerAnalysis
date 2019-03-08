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
traces_list_pepClusters <- readRDS("traces_list_pepClusters.rds")
design_matrix <- readRDS("design_matrix.rda")

#' ## Integrate traces across conditions
traces_sum_pepClusters <- integrateTraceIntensities(traces_list_pepClusters,
                                          design_matrix = NULL,
                                          integrate_within = NULL,
                                          aggr_fun = "sum")

saveRDS(traces_sum_pepClusters,"traces_sum_pepClusters.rds")

#' Subset traces for testing purposes
#testProts <- c(unique(traces_sum_pepClusters$trace_annotation$protein_id)[1:50],"P22102","Q6P2Q9","P14618","P22102","O00468","Q9UBF2",
#               "O15067","O75822","Q9Y266",
#               "P42167","P61978","O95801",
#               "P55060","O14578","O75717","P35221","Q9HB71","Q9H6S3")
#subsetTest <- subset(traces_sum_pepClusters,trace_subset_ids=testProts,trace_subset_type = "protein_id")
#saveRDS(subsetTest,"subsetTest.rda")

#' ## Proteoform-specific feature finding
proteinFeatures  <- findProteinFeatures(traces=traces_sum_pepClusters,
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