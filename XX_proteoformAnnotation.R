# proteoform annotation
library(devtools)
setwd("~/mysonas/CCprofiler")
load_all()
setwd("~/mysonas/PRPF8/analysis/output")
#install_github("CCprofiler/CCprofiler", ref = "DA_module")
#library(CCprofiler)

pepTraces <- readRDS("pepTracesImpSPCmultipep.rda")
design_matrix <- readRDS("design_matrix.rda")
