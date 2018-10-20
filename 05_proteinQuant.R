# protein quentification
library(devtools)
setwd("~/mysonas/CCprofiler")
load_all()
setwd("~/mysonas/PRPF8/analysis/output")
#install_github("CCprofiler/CCprofiler", ref = "DA_module")
#library(CCprofiler)

pepTraces <- readRDS("pepTracesImpSPCmultipep.rda")
design_matrix <- readRDS("design_matrix.rda")

#protTraces <- proteinQuantification(pepTraces,
#                       topN = 2,
#                       keep_less = FALSE,
#                       rm_decoys = TRUE,
#                       use_sibPepCorr = FALSE,
#                       use_repPepCorr = FALSE,
#                       full_intersect_only = FALSE)

protTraces <- proteinQuantification(pepTraces,
                      topN = 1000,
                      keep_less = TRUE,
                      rm_decoys = TRUE,
                      use_sibPepCorr = FALSE,
                      use_repPepCorr = FALSE,
                      full_intersect_only = FALSE)

protTraces <- updateTraces(protTraces)

protTracesSum <- integrateTraceIntensities(protTraces,
                                       design_matrix = NULL,
                                       integrate_within = NULL,
                                       aggr_fun = "sum")

saveRDS(protTraces,"protTraces.rds")
saveRDS(protTracesSum,"protTracesSum.rds")
