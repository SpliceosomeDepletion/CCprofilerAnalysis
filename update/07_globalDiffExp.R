#' ---
#' title: "Global differential expression analysis"
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
#install_github("CCprofiler/CCprofiler", ref = "DA_module")
#library(CCprofiler)
if (length(grep("nas21.ethz.ch",getwd()))>0) {
  setwd("~/mysonas/CCprofiler")
  load_all()
  setwd("~/mysonas/PRPF8/analysis/output")
  knitr::opts_knit$set(root.dir = '~/mysonas/PRPF8/analysis/output')
} else {
  setwd("/Volumes/ibludau-1/CCprofiler")
  load_all()
  setwd("/Volumes/ibludau-1/PRPF8/analysis/output")
  knitr::opts_knit$set(root.dir = '/Volumes/ibludau-1/PRPF8/analysis/output')
}

#' ## Load data
design_matrix <- readRDS("design_matrix.rda")
pepTraces <- readRDS("pepTraces_proteoformMulti.rds")
#testProts <- c(unique(pepTraces$minus$trace_annotation$protein_id)[1:200],"P22102","Q6P2Q9","P14618")
#pepTraces <- subset(pepTraces,trace_subset_ids=testProts,trace_subset_type = "protein_id")

#' Perform differential expression analysis
globalDiffExp <- testGlobalDifferentialExpression(pepTraces,design_matrix=design_matrix)
saveRDS(globalDiffExp, "globalDiffExp.rda")

plotVolcano(globalDiffExp$diffPeptides, PDF = T, name = "globalDiffPeptides")
plotVolcano(globalDiffExp$diffProteins, PDF = T, name = "globalDiffProteins")
plotVolcano(globalDiffExp$diffProteoforms, PDF = T, name = "globalDiffProteoforms")

plotVolcano(globalDiffExp$diffProteins, highlight=c("Q6P2Q9"), PDF = T, name = "prot_DiffExprProt_PRPF8")
plotVolcano(globalDiffExp$diffProteoforms, highlight=c("Q6P2Q9"), PDF = T, name = "prot_DiffExprProteoform_PRPF8")
plotVolcano(globalDiffExp$diffProteins, highlight=c("O75717"), PDF = T, name = "prot_DiffExprProt_O75717")
plotVolcano(globalDiffExp$diffProteoforms, highlight=c("O75717"), PDF = T, name = "prot_DiffExprProteoform_O75717")


proteoform_stats_O75717 <- globalDiffExp$diffProteoforms[protein_id=="O75717"]
setnames(proteoform_stats_O75717,c("global_int1_imp","global_int2_imp"),c("minus","plus"))
proteoform_stats_O75717.m <- melt(proteoform_stats_O75717,id.vars=c("proteoform_id","protein_id","sumLog2FC","pBHadj","qVal","pVal"))
pdf(paste0("proteoform_stats_O75717",".pdf"), height=4, width=4)
ggplot(proteoform_stats_O75717.m,aes(x=variable,y=value,fill=proteoform_id,group=variable)) + 
  geom_bar(stat="identity") +
  facet_grid(.~proteoform_id) +
  theme_classic() +
  scale_fill_manual(values=c("#E69F00", "#56B4E9"))
dev.off()

proteins_up <- globalDiffExp$diffProteins[(sumLog2FC <= -1) & (pBHadj < 0.01)]
proteins_down <- globalDiffExp$diffProteins[(sumLog2FC >= 1) & (pBHadj < 0.01)]
write.table(proteins_up,"globalProteinsUp.tsv",sep="\t",quote=F,row.names=F)
write.table(proteins_down,"globalProteinsDown.tsv",sep="\t",quote=F,row.names=F)



