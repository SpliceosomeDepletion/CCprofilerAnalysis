#' ---
#' title: "Proteoform annotation"
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
library(nFactors)
#install_github("CCprofiler/CCprofiler", ref = "DA_module")
#library(CCprofiler)
if (length(grep("nas21.ethz.ch",getwd()))>0) {
  setwd("~/mysonas/CCprofiler")
  load_all()
  setwd("~/mysonas/PRPF8/analysis/output/Depleted")
  knitr::opts_knit$set(root.dir = '~/mysonas/PRPF8/analysis/output/Depleted')
} else {
  setwd("/Volumes/ibludau-1/CCprofiler")
  load_all()
  setwd("/Volumes/ibludau-1/PRPF8/analysis/output/Depleted")
  knitr::opts_knit$set(root.dir = '/Volumes/ibludau-1/PRPF8/analysis/output/Depleted')
}

#' ## Load data
pep_traces <- readRDS("../pepTracesRaw.rda")
calibration <-  readRDS("../calibration.rds")

pep_traces <- pep_traces$depleted

pepTracesMV <- findMissingValues(pep_traces,
                                 bound_left = 1,
                                 bound_right = 1,
                                 consider_borders = TRUE)

#' Impute NA values by fitting a spline:
pepTracesImp <- imputeMissingVals(pepTracesMV, method = "spline")

#' Plot imputation summary:
plotImputationSummary(pepTracesMV, pepTracesImp, PDF = T, plot_traces = T, max_n_traces = 2)

#' ## Filter by consecutive IDs
pepTracesConsIds <- filterConsecutiveIdStretches(pepTracesImp,
                                                 min_stretch_length = 2,
                                                 remove_empty = T)

pepTracesConsIds_multiPep <- filterSinglePeptideHits(pepTracesConsIds)

pepTraces_maxCorr <- filterByMaxCorr(pepTracesConsIds_multiPep,
                                     cutoff = 0.9, plot = T, PDF=T, name="maxCorrHist")

#' ## Calculate Sibling Peptide Correlation
pepTraces_maxCorr_SPC <- calculateSibPepCorr(pepTraces_maxCorr,
                                             plot = T, PDF = T,
                                             name = "SibPepCorr_densityplot_afterMaxCorrFilter")


# traces_list <- readRDS("traces_coordinates_list")

#' ## Combine traces across conditions
#' Combine traces across conditions by appending fractions to generate a
#' large traces object with twice as many fractions. This is done in order to
#' allow consistent proteoform detection across conditions.
traces_combined <- copy(pepTraces_maxCorr_SPC)

#' ## Compute the minimum correlation between any two sibling peptides
traces_combined_minCorr <- calculateMinCorr(traces_combined,
  plot = TRUE, PDF=TRUE)

#' ## Estimate proteoform p-values for each gene/protein
#' The null hypothesis is that each gene/protein only consists of one
#' proteoform, meaning that all its sibling peptides have a high correlation
#' along the fractionation dimension. A gamma distribution is fitted to the
#' high-correlating end of the minimum correlation value distribution and
#' p-values are estimated accordingly.
traces_combined_minCorr_pval <- estimateProteoformPval(traces_combined_minCorr,
  plot = TRUE, PDF=TRUE)

#' ## Determine proteoforms by peptide clustering
#' Peptides of each gene/protein are clustered based on the correlation of all sibling peptides.
traces_combined_minCorr_pval_clustered <- clusterPeptides(traces_combined_minCorr_pval,
  clusterN = NULL, clusterH = 0.6, nFactorAnalysis = FALSE,
  plot = TRUE, PDF=TRUE)

#' ## Evaluate clustering results
#' ##### Differentiate genes/proteins with one vs. multiple proteoforms based on minimum peptide correlation criterion.
single_proteoform <- traces_combined_minCorr_pval_clustered$trace_annotation[proteoform_pval_adj > 0.05]
multi_proteoforms <- traces_combined_minCorr_pval_clustered$trace_annotation[proteoform_pval_adj <= 0.05]

single_proteoform[,n_proteoforms:=length(unique(proteoform_id)),by="protein_id"]
multi_proteoforms[,n_proteoforms:=length(unique(proteoform_id)),by="protein_id"]

cat(paste0("Genes/proteins with 1 proteoform: ",length(unique(single_proteoform$protein_id))))
cat(paste0("Genes/proteins with >1 proteoform: ",length(unique(multi_proteoforms$protein_id))))
cat(paste0("# proteoforms: ",length(unique(multi_proteoforms$proteoform_id))))

#' ##### Number of clusters for proteins with 1 vs multiple proteoforms:
#' Contingency table of clusters in genes/proteins with 1 proteoform:
table(unique(subset(single_proteoform,select=c("protein_id","n_proteoforms")))$n_proteoforms)
#' Contingency table of clusters in genes/proteins with >1 proteoform:
table(unique(subset(multi_proteoforms,select=c("protein_id","n_proteoforms")))$n_proteoforms)

#' ## Protein with a proteoform p-value <0.05 should be set to 1 proteoform regardless of clustering
traces_combined_minCorr_pval_clustered_fixed <- copy(traces_combined_minCorr_pval_clustered)
traces_combined_minCorr_pval_clustered_fixed$trace_annotation[,proteoform_id := ifelse(proteoform_pval_adj > 0.05, protein_id, proteoform_id)]
traces_combined_minCorr_pval_clustered_fixed$trace_annotation[,n_proteoforms := length(unique(proteoform_id)), by="protein_id"]

#' ## Save objects
saveRDS(traces_combined_minCorr_pval_clustered_fixed,"traces_combined_minCorr_pval_clustered_fixed.rds")

test_proteins <- c(
  "P42574", # caspase 3
  "P55212", # caspase 6
  "P55210", # caspase 7
  "Q14790", # caspase 8
  "P55211", # caspase 9
  "P07858", # cathepsin B
  "P07711", # cathepsin L1
  "O60911", # cathepsin L2
  "P07339", # cathepsin D
  "P10144", # granzyme B
  "P0DJD7", # pepsin A4
  "P0DJD9", # pepsin A5
  "P0DJD8", #pepsin A3
  "P07477" # trypsin 1
)

pdf(paste0("PeptideTraces_proteases.pdf"))
for (test_protein in test_proteins){
  pepTest <- subset(traces_combined_minCorr_pval_clustered_fixed, trace_subset_ids = test_protein, trace_subset_type = "protein_id")
  test_protein_name <- pepTest$trace_annotation$Entry_name[1]
  if(nrow(pepTest$trace_annotation) > 0 ){
    #pdf(paste0("PeptideTraces_proteoformAnn_",test_protein,"_",test_protein_name,".pdf"))
    plot(pepTest, log = FALSE, PDF = F,
         name = paste0("PeptideTraces_proteoformAnn_",test_protein,"_",test_protein_name), 
         plot = TRUE, highlight = NULL, highlight_col = NULL,
         colour_by = "proteoform_id")
    #dev.off()
  }
}
dev.off()

##

test_proteins <- c(traces_combined_minCorr_pval_clustered_fixed$trace_annotation[n_proteoforms>1]$protein_id)[1:100]

pdf(paste0("PeptideTraces_proteoforms.pdf"))
for (test_protein in test_proteins){
  pepTest <- subset(traces_combined_minCorr_pval_clustered_fixed, trace_subset_ids = test_protein, trace_subset_type = "protein_id")
  test_protein_name <- pepTest$trace_annotation$Entry_name[1]
  if(nrow(pepTest$trace_annotation) > 0 ){
    #pdf(paste0("PeptideTraces_proteoformAnn_",test_protein,"_",test_protein_name,".pdf"))
    plot(pepTest, log = FALSE, PDF = F,
         name = paste0("PeptideTraces_proteoformAnn_",test_protein,"_",test_protein_name), 
         plot = TRUE, highlight = NULL, highlight_col = NULL,
         colour_by = "proteoform_id")
    #dev.off()
  }
}
dev.off()

