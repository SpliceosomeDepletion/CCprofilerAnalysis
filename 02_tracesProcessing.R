library(devtools)
setwd("~/mysonas/CCprofiler")
load_all()
setwd("~/mysonas/PRPF8/analysis/output")
#install_github("CCprofiler/CCprofiler", ref = "DA_module")
#library(CCprofiler)

traces_list <- readRDS("pepTracesRaw.rda")
design_matrix <- readRDS("design_matrix.rda")

#traces_list <- subset(traces_list,trace_subset_ids="P61978",trace_subset_type="protein_id")


traces_integrated <- integrateTraceIntensities(traces_list,
                                          design_matrix = NULL,
                                          integrate_within = NULL,
                                          aggr_fun = "sum")

map_table <- readRDS("../../data/proteogenomics/HeLa_protgen_combined_1FPKM_header_mapping_v2.rda")

#' ## Annotate the Leading isoform for each peptide
traces_integrated_iso <- annotateLeadingIsoform(traces_integrated,
                                 isoform_col = "isoform_id",
                                 output_col = "LeadingIsoform")

#' ## Find the relative position of each peptide in the protein sequence
#' With the leadingIsoform we can use the mapping table that contains the sequences of the
#' fasta file used for the proteomics search to determine the relative position of each peptide
#' within the protein.
traces_integrated_iso_pos <- annotateRelativePepPos(traces_integrated_iso,
  map_table, multimatch = "first", verbose = T)

#' Select one peptide as representative of peptides all
#' starting at the same starting position
traces_integrated_start <- summarizeAlternativePeptideSequences(
  traces_integrated_iso_pos, topN=1,position="PeptidePositionStart")

traces_integrated_start_end <- summarizeAlternativePeptideSequences(
  traces_integrated_start, topN=1,position="PeptidePositionEnd")

#setkey(traces_integrated_start$trace_annotation,"PeptidePositionStart")
#setkey(traces_integrated_start_end$trace_annotation,"PeptidePositionStart")

# subset traces list to selected peptides
traces_list <- subset(traces_list,
  unique(traces_integrated_start_end$traces$id))

# @TODO genomic position info not transferred at the moment

# traces alignment test
alignTraces(traces_list, plot = T, PDF = T)

# total intensity test
plotGlobalIntensities(traces_list, plot = T, PDF = T)

# Convert 0's in missing value locations to NA
pepTracesMV <- findMissingValues(traces_list,
                                 bound_left = 2,
                                 bound_right = 2,
                                 consider_borders = TRUE)

# Impute NA values
pepTracesImp <- imputeMissingVals(pepTracesMV, method = "spline")

# Plot summary
plotImputationSummary(pepTracesMV, pepTracesImp, PDF = T, plot_traces = T, max_n_traces = 2)

# Median normalization >> does not have large impact on SibPepCore summary statistic
# pepTracesMedianNorm <- smootheTraces(pepTracesImp, method = "median", plot = T, PDF = T, smoothe_on_targets = T, span = 0.1)

# Filter by consecutive IDs
pepTracesConsIds <- filterConsecutiveIdStretches(pepTracesImp,
                                             min_stretch_length = 2,
                                             remove_empty = T)

# Calculate Sibling Peptide Correlation
pepTracesImpSPC <- calculateSibPepCorr(pepTracesConsIds,
                                   plot = T, PDF = T)

# Remove proteins with single peptide
minus_multipep <- unique(subset(pepTracesImpSPC$minus$trace_annotation,! is.na(SibPepCorr))$protein_id)
plus_multipep <- unique(subset(pepTracesImpSPC$plus$trace_annotation,! is.na(SibPepCorr))$protein_id)
multipep <- unique(c(minus_multipep,plus_multipep))
pepTracesImpSPCmultipep <- subset(pepTracesImpSPC,trace_subset_ids=multipep,trace_subset_type="protein_id")

# Filter by Sibling Peptide Correlation
# pepTracesImpSPC <- filterBySibPepCorr(pepTracesImpSPC, absolute_spcCutoff = 0 ,fdr_cutoff = NULL, plot = T, PDF = T)

# Add mock eplicate peptide correlation
# @TODO: Update functions not to depend on this!
pepTracesImpSPCmultipep$minus$trace_annotation$RepPepCorr=1
pepTracesImpSPCmultipep$plus$trace_annotation$RepPepCorr=1

pepTracesImpSPCmultipep <- updateTraces(pepTracesImpSPCmultipep)

saveRDS(pepTracesImpSPCmultipep, "pepTracesImpSPCmultipep.rda")

# Trace merging for combined protein feature finding
pepTracesImpSPCmultipepSum <- integrateTraceIntensities(pepTracesImpSPCmultipep,
                                          design_matrix = NULL,
                                          integrate_within = NULL,
                                          aggr_fun = "sum")
saveRDS(pepTracesImpSPCmultipepSum,"pepTracesImpSPCmultipepSum.rds")
