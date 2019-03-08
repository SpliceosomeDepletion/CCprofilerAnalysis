#' # Benchmark proteoform detection
#' ## Negative set: proteins with only one proteoform
#' proteins with min sib pep corr > 0.95
#' ## Positive set: proteins with only >= 2 proteoforms
#' Mixture of peptides of protein pairs in negative set
#' maybe exclude protein pairs with min STRING pathway length <= 2
library(devtools)

setwd("~/Desktop/CCprofiler/CCprofiler")
load_all()
setwd("~/Desktop/PRPF8data/benchmark")

pepTracesSum <- readRDS("pepTracesSum.rda")

proteins_minPepCorr_095 <- unique(subset(pepTracesSum$trace_annotation, mincorr >= 0.95)$protein_id)

negative_traces <- subset(pepTracesSum, proteins_minPepCorr_095,trace_subset_type = "protein_id")

odd_idx <- seq(1,length(proteins_minPepCorr_095),2)
even_idx <- seq(2,length(proteins_minPepCorr_095),2)

combinations <- unlist(lapply(seq_len(floor(length(proteins_minPepCorr_095)/2)),function(i){
  paste(c(proteins_minPepCorr_095[odd_idx[i]],proteins_minPepCorr_095[even_idx[i]]), collapse = ";")
}))

res <- lapply(combinations,function(x){
  proteins <- unlist(strsplit(x, split = ";"))
  traces_ann <- subset(negative_traces$trace_annotation, protein_id %in% proteins)
  traces_ann[,protein_id:=x]
  traces_ann[,protein_name:=x]
})
resDT <- do.call(rbind,res)
resDT <- resDT[negative_traces$traces$id, on="id"]

positive_traces <- copy(negative_traces)
positive_traces$trace_annotation <- resDT
positive_traces$trace_annotation[,id:=paste0(id,"_mixed")]
positive_traces$traces[,id:=paste0(id,"_mixed")]



