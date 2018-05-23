# protein quentification
library(devtools)
#setwd("~/mysonas/CCprofiler")
#load_all()
setwd("~/mysonas/PRPF8/analysis/output")
install_github("CCprofiler/CCprofiler", ref = "DA_module")
library(CCprofiler)

protTracesSum <- readRDS("protTracesSum.rds")
corumHypotheses <- readRDS("~/mysonas/html/SECpaper/output_data/corum/complexTargetsPlusDecoys.rds")

#testHyp <- subset(corumHypotheses,complex_id %in% unique(corumHypotheses$complex_id)[1:20])
complexFeatures <-  findComplexFeatures(traces = protTracesSum,
                                            complex_hypothesis = corumHypotheses,
                                            corr_cutoff = 0.9,
                                            window_size = 8,
                                            parallelized = T,
                                            n_cores = 10,
                                            collapse_method = "apex_network",
                                            perturb_cutoff = "5%",
                                            rt_height = 1,
                                            smoothing_length = 7)

saveRDS(complexFeatures,"corumComplexFeatures.rda")
