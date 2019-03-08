#' ---
#' title: "Proteoform & Protein Quantification"
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

#' ## Remove peptides with no proteoform annotation
#peptides_annotated_minus <- peptide_traces_list$minus$trace_annotation[which(!is.na(peptide_traces_list$minus$trace_annotation$proteoform_id))]$id
#peptides_annotated_plus <- peptide_traces_list$plus$trace_annotation[which(!is.na(peptide_traces_list$plus$trace_annotation$proteoform_id))]$id
#peptide_traces_list_sub <- copy(peptide_traces_list)
#peptide_traces_list_sub$minus <- subset(peptide_traces_list_sub$minus, trace_subset_ids = peptides_annotated_minus)
#peptide_traces_list_sub$plus <- subset(peptide_traces_list_sub$plus, trace_subset_ids = peptides_annotated_plus)

#' ## Quantify on protein level
protein_traces_list <- proteinQuantification(peptide_traces_list,
                                             quantLevel="protein_id",
                                             topN = 1000,
                                             keep_less = TRUE,
                                             rm_decoys = TRUE,
                                             use_sibPepCorr = FALSE,
                                             use_repPepCorr = FALSE,
                                             full_intersect_only = FALSE,
                                             verbose = FALSE)

#' ## Update fraction annotation for protein traces
protein_traces_list <- updateTraces(protein_traces_list)

#' ## Summarize protein traces across conditions
protein_traces_sum <- integrateTraceIntensities(protein_traces_list,
                                                design_matrix = NULL,
                                                integrate_within = NULL,
                                                aggr_fun = "sum")

#' ## Save protein traces
saveRDS(protein_traces_list,"protein_traces_list.rds")
saveRDS(protein_traces_sum,"protein_traces_sum.rds")

#' ## Quantify on proteoform level
proteoform_traces_list <- proteinQuantification(peptide_traces_list,
                                                quantLevel="proteoform_id",
                                                topN = 1000,
                                                keep_less = TRUE,
                                                rm_decoys = TRUE,
                                                use_sibPepCorr = FALSE,
                                                use_repPepCorr = FALSE,
                                                full_intersect_only = FALSE,
                                                verbose = FALSE)

#' ## Update fraction annotation for proteoform traces
proteoform_traces_list <- updateTraces(proteoform_traces_list)

#' ## Summarize proteoform traces across conditions
proteoform_traces_sum <- integrateTraceIntensities(proteoform_traces_list,
                                                   design_matrix = NULL,
                                                   integrate_within = NULL,
                                                   aggr_fun = "sum")

#' ## Save proteoform traces
saveRDS(proteoform_traces_list,"proteoform_traces_list.rds")
saveRDS(proteoform_traces_sum,"proteoform_traces_sum.rds")

