#' # Protein feature finding
#'
#' ## Environment setup
#'
library(devtools)
#setwd("~/mysonas/CCprofiler")
#load_all()
setwd("/nfs/nas21.ethz.ch/nas/fs2102/biol_ibt_usr_s1/ibludau/PRPF8/analysis/output")
install_github("CCprofiler/CCprofiler", ref = "DA_module")
library(CCprofiler)

#' ## Load data
pepTracesSum <- readRDS("pepTracesImpSPCmultipepSum.rds")
pepTraces <- readRDS("pepTracesImpSPCmultipep.rda")
design_matrix <- readRDS("design_matrix.rda")

#' #### Feature finding
#' The protein centric analysis focusses on differences in co-elution events along the SEC separation dimension. We use the feature finding functionality to define the borders of these elution events. This part will be run on euler.
#' We perform feature finding both on the traces from the individual conditions as well as on the integrated traces. The individual features are used mainly for isoform scoring, while the summed
#' features are used for differential analysis between conditions.
#+ cache = T
#proteinFeatures  <- findProteinFeatures(traces=pepTracesSum,
#                                             corr_cutoff=0.9,
#                                             window_size=8,
#                                             parallelized=TRUE,
#                                             n_cores=10,
#                                             collapse_method="apex_only",
#                                             perturb_cutoff= "5%",
#                                             rt_height=1,
#                                             smoothing_length=7,
#                                             useRandomDecoyModel=TRUE)
#
#print("Protein feature finding on summed traces done.")
#getwd()
#saveRDS(proteinFeatures, "integrated_proteinFeatures.rda")

proteinFeaturesMinus  <- findProteinFeatures(traces=pepTraces$minus,
                                             corr_cutoff=0.9,
                                             window_size=8,
                                             parallelized=TRUE,
                                             n_cores=10,
                                             collapse_method="apex_only",
                                             perturb_cutoff= "5%",
                                             rt_height=1,
                                             smoothing_length=7,
                                             useRandomDecoyModel=TRUE)

#proteinFeaturesPlus  <- findProteinFeatures(traces=pepTraces$plus,
#                                             corr_cutoff=0.9,
#                                             window_size=8,
#                                             parallelized=TRUE,
#                                             n_cores=10,
#                                             collapse_method="apex_only",
#                                             perturb_cutoff= "5%",
#                                             rt_height=1,
#                                             smoothing_length=7,
#                                             useRandomDecoyModel=TRUE)

#proteinFeaturesInd  <- lapply(pepTraces, findProteinFeatures,
#                                             corr_cutoff=0.9,
#                                             window_size=8,
#                                             parallelized=TRUE,
#                                             n_cores=5,
#                                             collapse_method="apex_only",
#                                             perturb_cutoff= "5%",
#                                             rt_height=1,
#                                             smoothing_length=7,
#                                             useRandomDecoyModel=TRUE)
print("Protein feature finding on individual traces done.")
getwd()

saveRDS(proteinFeaturesMinus, "proteinFeaturesMinus.rda")

#' Protein features are filtered with an FDR cutoff based on the detection of decoy proteins.
#+ cache = T, message =F
filteredDataMinus <- scoreFeatures(proteinFeaturesMinus, FDR=0.05, PDF=T, name="qvalueStats_proteinFeaturesMinus")
summarizeFeatures(filteredDataMinus,plot=TRUE,PDF=TRUE,name="proteinFeature_summaryMinus")
# plotGlobalProteomeAssemblyState(filteredData, PDF=TRUE, name="globalProteomeAssemblyState_integrated")


pdf("proteinFeature_summary_individual_poster.pdf",height=4.5,width=3)
summarizeFeatures(filteredDataMinus,plot=TRUE,PDF=FALSE,name="proteinFeature_summaryMinus")
dev.off()
