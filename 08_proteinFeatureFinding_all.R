#' # Protein feature finding
#'
#' ## Environment setup
#'
library(devtools)
#setwd("~/mysonas/CCprofiler")
#load_all()
setwd("~/mysonas/PRPF8/analysis/output")
install_github("CCprofiler/CCprofiler", ref = "DA_module")
library(CCprofiler)

#' ## Load data
pepTracesSum <- readRDS("pepTracesImpSPCmultipepSum.rds")
design_matrix <- readRDS("design_matrix.rda")

#' #### Feature finding
#' The protein centric analysis focusses on differences in co-elution events along the SEC separation dimension. We use the feature finding functionality to define the borders of these elution events. This part will be run on euler.
#+ cache = T
proteinFeatures  <- findProteinFeatures(traces=pepTracesSum,
                                             corr_cutoff=0.9,
                                             window_size=8,
                                             parallelized=TRUE,
                                             n_cores=10,
                                             collapse_method="apex_only",
                                             perturb_cutoff= "5%",
                                             rt_height=1,
                                             smoothing_length=7,
                                             useRandomDecoyModel=TRUE)
getwd()
saveRDS(proteinFeatures, "integrated_proteinFeatures.rda")

#' Protein features are filtered with an FDR cutoff based on the detection of decoy proteins.
#+ cache = T, message =F
filteredData <- scoreFeatures(proteinFeatures, FDR=0.05, PDF=T, name="qvalueStats_proteinFeatures_integrated")
summarizeFeatures(filteredData,plot=TRUE,PDF=TRUE,name="proteinFeature_summary_integrated")
# plotGlobalProteomeAssemblyState(filteredData, PDF=TRUE, name="globalProteomeAssemblyState_integrated")

#' Saving the feature data
saveRDS(filteredData, paste0("proteinFeatures_005_FDR_integrated.rds"))
write.table(filteredData,paste0("proteinFeatures_005_FDR_integrated.txt"),sep="\t",quote=F,row.names = F)
summarizeFeatures(filteredData,plot=TRUE,PDF=TRUE,name="proteinFeature_summary_integrated")
