library(devtools)
#setwd("~/mysonas/CCprofiler")
#load_all()
setwd("~/mysonas/PRPF8/analysis/output")
install_github("CCprofiler/CCprofiler", ref = "DA_module")
library(CCprofiler)

traces_list <- readRDS("pepTracesRaw.rda")
design_matrix <- readRDS("design_matrix.rda")

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
# pepTracesConsIds <- filterConsecutiveIdStretches(pepTracesMedianNorm,
#                                             min_stretch_length = 2,
#                                             remove_empty = T)

# Calculate Sibling Peptide Correlation
pepTracesImpSPC <- calculateSibPepCorr(pepTracesImp,
                                   plot = T, PDF = T)

# Remove proteins with single peptide
minus_multipep <- unique(subset(pepTracesImpSPC$minus$trace_annotation,! is.na(SibPepCorr))$protein_id)
plus_multipep <- unique(subset(pepTracesImpSPC$plus$trace_annotation,! is.na(SibPepCorr))$protein_id)
multipep <- unique(c(minus_multipep,plus_multipep))
pepTracesImpSPCmultipep <- subset(pepTracesImpSPC,trace_subset_ids=multipep,trace_subset_type="protein_id")

# Filter by Sibling Peptide Correlation
# pepTracesImpSPC <- filterBySibPepCorr(pepTracesImpSPC, absolute_spcCutoff = 0 ,fdr_cutoff = NULL, plot = T, PDF = T)

# Add mock eplicate peptide correlation
# Update functions not to depend on this!
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
