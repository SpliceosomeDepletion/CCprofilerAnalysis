#' # Load CCprofiler package and set working directory
library(devtools)
options(warn=-1)
if (length(grep("nas22.ethz.ch",getwd()))>0) {
  setwd("~/myspectrumscale/CCprofiler")
  suppressPackageStartupMessages(load_all())
  setwd("~/myspectrumscale/PRPF8/analysis/output")
  knitr::opts_knit$set(root.dir = '~/myspectrumscale/PRPF8/analysis/output')
} else if (length(grep("/Volumes/ibludau-1/",getwd()))>0) {
  setwd("/Volumes/ibludau-1/CCprofiler")
  load_all()
  suppressPackageStartupMessages(setwd("/Volumes/ibludau-1/PRPF8/analysis/output"))
  knitr::opts_knit$set(root.dir = '/Volumes/ibludau-1/PRPF8/analysis/output')
} else {
  setwd("~/Desktop/CCprofiler_analysis/CCprofiler")
  load_all()
  setwd("~/Desktop/CCprofiler_analysis/PRPF8")
}

