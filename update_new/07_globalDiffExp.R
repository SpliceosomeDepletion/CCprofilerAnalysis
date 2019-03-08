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
setwd("~/Desktop/CCprofiler/CCprofiler")
load_all()
setwd("~/Desktop/PRPF8data/")

#' ## Load data
design_matrix <- readRDS("design_matrix.rda")
pepTraces <- readRDS("pepTraces_proteoformMulti.rds")
#testProts <- c(unique(pepTraces$minus$trace_annotation$protein_id)[500:550],"P22102","Q6P2Q9","P14618")
#pepTraces <- subset(pepTraces,trace_subset_ids=testProts,trace_subset_type = "protein_id")

#' ## Perform differential expression analysis
globalDiffExp <- testGlobalDifferentialExpression(pepTraces,design_matrix=design_matrix)
saveRDS(globalDiffExp, "globalDiffExp.rda")

globalDiffExp <- readRDS("globalDiffExp.rda")


#' ## Plot global volcano plots
plotVolcano(globalDiffExp$diffPeptides, PDF = T, name = "globalDiffPeptides")
plotVolcano(globalDiffExp$diffProteins, PDF = T, name = "globalDiffProteins")
plotVolcano(globalDiffExp$diffProteoforms, PDF = T, name = "globalDiffProteoforms")
# PRPF8 volcano
plotVolcano(globalDiffExp$diffProteins, highlight=c("Q6P2Q9"), PDF = T, name = "prot_DiffExprProt_PRPF8")
plotVolcano(globalDiffExp$diffProteoforms, highlight=c("Q6P2Q9"), PDF = T, name = "prot_DiffExprProteoform_PRPF8")
# Labmeeting volcano
plotVolcano(globalDiffExp$diffProteins, highlight=c("O75717"), PDF = T, name = "prot_DiffExprProt_O75717")
plotVolcano(globalDiffExp$diffProteoforms, highlight=c("O75717"), PDF = T, name = "prot_DiffExprProteoform_O75717")
# Labmeeting figure LAP2 P42166 P42167
plotVolcano(globalDiffExp$diffProteins, highlight=c("P42167"), PDF = T, name = "prot_DiffExprProt_LAP2")
plotVolcano(globalDiffExp$diffProteoforms, highlight=c("P42167"), PDF = T, name = "prot_DiffExprProteoform_LAP2")
# Labmeeting figure HNRNPA1 P09651
plotVolcano(globalDiffExp$diffProteins, highlight=c("P09651"), PDF = T, name = "prot_DiffExprProt_HNRNPA1")
plotVolcano(globalDiffExp$diffProteoforms, highlight=c("P09651"), PDF = T, name = "prot_DiffExprProteoform_HNRNPA1")
# Labmeeting figure HNRNPL P14866
plotVolcano(globalDiffExp$diffProteins, highlight=c("P14866"), PDF = T, name = "prot_DiffExprProt_HNRNPL")
plotVolcano(globalDiffExp$diffProteoforms, highlight=c("P14866"), PDF = T, name = "prot_DiffExprProteoform_HNRNPL")

#' ## Plot intensity barchart
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

proteoform_stats_P14866 <- globalDiffExp$diffProteoforms[protein_id=="P14866"]
setnames(proteoform_stats_P14866,c("global_int1_imp","global_int2_imp"),c("minus","plus"))
proteoform_stats_P14866.m <- melt(proteoform_stats_P14866,id.vars=c("proteoform_id","protein_id","sumLog2FC","pBHadj","qVal","pVal"))
pdf(paste0("proteoform_stats_P14866",".pdf"), height=4, width=4)
ggplot(proteoform_stats_P14866.m,aes(x=variable,y=value,fill=proteoform_id,group=variable)) + 
  geom_bar(stat="identity") +
  facet_grid(.~proteoform_id) +
  theme_classic() +
  scale_fill_manual(values=c("#E69F00", "#56B4E9", "grey"))
dev.off()


#' ## Save differentially abundant proteins for DAVID analysis
all_proteins <- unique(globalDiffExp$diffProteins$protein_id)
proteins_changed <- globalDiffExp$diffProteins[(abs(sumLog2FC) >= 1) & (pBHadj < 0.01)]$protein_id
proteins_unchanged <- all_proteins[which(!all_proteins %in% proteins_changed)]
proteins_up <- unique(globalDiffExp$diffProteins[(sumLog2FC <= -1) & (pBHadj < 0.01)]$protein_id)
proteins_down <- unique(globalDiffExp$diffProteins[(sumLog2FC >= 1) & (pBHadj < 0.01)]$protein_id)
write.table(proteins_up,"globalProteinsUp.tsv",sep="\t",row.names = F,col.names = F,quote=F)
write.table(proteins_down,"globalProteinsDown.tsv",sep="\t",row.names = F,col.names = F,quote=F)
write.table(all_proteins,"all_proteins.txt",sep="\t",row.names = F,col.names = F,quote=F)
write.table(proteins_changed,"proteins_changed.txt",sep="\t",row.names = F,col.names = F,quote=F)
write.table(proteins_unchanged,"proteins_unchanged.txt",sep="\t",row.names = F,col.names = F,quote=F)

all_proteoforms <- unique(globalDiffExp$diffProteoforms$proteoform_id)
proteoforms_changed <- globalDiffExp$diffProteoforms[(abs(sumLog2FC) >= 1) & (pBHadj < 0.01)]$proteoform_id
proteoforms_unchanged <- all_proteoforms[which(!all_proteoforms %in% proteoforms_changed)]
proteoforms_up <- unique(globalDiffExp$diffProteoforms[(sumLog2FC <= -1) & (pBHadj < 0.01)]$proteoform_id)
proteoforms_down <- unique(globalDiffExp$diffProteoforms[(sumLog2FC >= 1) & (pBHadj < 0.01)]$proteoform_id)
length(proteoforms_up)
length(proteoforms_down)
prot_proteoforms_up <- unique(globalDiffExp$diffProteoforms[(sumLog2FC <= -1) & (pBHadj < 0.01)]$protein_id)
prot_proteoforms_down <- unique(globalDiffExp$diffProteoforms[(sumLog2FC >= 1) & (pBHadj < 0.01)]$protein_id)
prot_all_proteoforms <- unique(globalDiffExp$diffProteoforms$protein_id)
length(intersect(prot_proteoforms_up,prot_proteoforms_down))

globalDiffExp$diffProteoforms[protein_id %in% intersect(prot_proteoforms_up,prot_proteoforms_down)]

interesting_prot <- intersect(prot_proteoforms_up,prot_proteoforms_down)
for(p in interesting_prot) {
  testTraces <- subset(pepTraces,trace_subset_ids=p,trace_subset_type = "protein_id")
  plot(testTraces,legend = T, colour_by="proteoform_id",design_matrix = design_matrix,name=paste0(p,"_pepTraces"),PDF=T)
}

#' ## Are proteins with multiple proteoforms more affected by PRPF8 depletion 
proteoformStats <- copy(globalDiffExp$diffProteoforms)
proteoformStats[,nProteoforms:=length(unique(proteoform_id)),by=protein_id]
proteins_one_form <- unique(proteoformStats[nProteoforms==1]$protein_id)
proteins_multiple_forms <- unique(proteoformStats[nProteoforms>1]$protein_id)

diff_exp_multi_proteoforms <- intersect(proteins_multiple_forms,proteins_changed)

m <- length(proteins_multiple_forms)       # Genes IN GO term == genes multiple proteoforms
n <- length(proteins_one_form)       # Genes NOT IN GO term == genes single proteoform
k <- length(proteins_changed)       # Gene hits, that is, differentially expressed == all genes that are diff exp.
x <- c(length(diff_exp_multi_proteoforms):length(proteins_changed))  # Genes both IN GO term and differentially expressed 'hits' == genes with multiple proteoforms that are diff exp

probabilities <- dhyper(x, m, n, k, log = FALSE)
pvalue <- sum(probabilities)
pvalue


#' ## Enrichment analysis 
#' First perform enrichment analysis on DAVID online. 
#' Plot fold-change enrichment plot.
david <- fread("DAVIDenrichment_unchanged_proteins.txt")
david_sig <- subset(david,Benjamini <= 0.05 & `Fold Enrichment` > 1)
david_sig_up <- david_sig[Category=="UP_KEYWORDS"]
pdf("david_unchanged_enrichment.pdf",width=4.5,height=4)
p <- ggplot(david_sig_up, aes(x=Term, y=`Fold Enrichment`)) + 
  geom_bar(stat="identity") +
  theme_classic() + 
  coord_flip() #+ 
  #theme(legend.position="bottom") +
  #scale_fill_gradient2(low = "red", mid = "white", high = "blue", midpoint = 0.05)
print(p)
dev.off()

david <- fread("DAVIDenrichment_changed_proteins.txt")
david_sig <- subset(david,Benjamini <= 0.05 & `Fold Enrichment` > 1)
david_sig_up <- david_sig[Category=="UP_KEYWORDS"]
pdf("david_changed_enrichment.pdf",width=4.5,height=4)
p <- ggplot(david_sig_up, aes(x=Term, y=`Fold Enrichment`)) + 
  geom_bar(stat="identity") +
  theme_classic() + 
  coord_flip() #+ 
  #theme(legend.position="bottom") +
  #scale_fill_gradient2(low = "red", mid = "white", high = "blue", midpoint = 0.05)
print(p)
dev.off()

#' Extract terms from Yansheng's study for enrichment plotting
david_y <- subset(david, Term %in% c("GO:0003743~translation initiation factor activity","mRNA splicing","Ubl conjugation","Cell cycle","RNA-binding","Transcription","Ribosome biogenesis"))
david_y[,Benjamini:=ifelse(Benjamini > 0.15, 0.15, Benjamini)]
david_y$Term <- factor(david_y$Term,
                          levels=c("GO:0003743~translation initiation factor activity",
                                   "mRNA splicing",
                                   "Ubl conjugation","Cell cycle","RNA-binding",
                                   "Transcription","Ribosome biogenesis"),
                          ordered = T)
pdf("y_enrichment.pdf",width=6,height=4)
p <- ggplot(david_y, aes(x=Term, y=`Fold Enrichment`,fill=round(Benjamini,2))) + 
  geom_bar(stat="identity") +
  theme_classic() + 
  scale_x_discrete(limits = rev(levels(david_y$Term))) +
  coord_flip() + 
  theme(legend.position="bottom") +
  scale_fill_gradient2(low = "red", mid = "white", high = "blue", midpoint = 0.05)
print(p)
dev.off()

#' ## Expression comparison to Liu et al.
protein_diff_exp <- globalDiffExp$diffProteins 
protein_diff_exp[,minus:=global_int1_imp]
protein_diff_exp[,plus:=global_int2_imp]

protein_diff_exp_long <- subset(protein_diff_exp,select=c("protein_id","minus","plus"))
protein_diff_exp_long <- melt(protein_diff_exp_long,id.vars = "protein_id")
setnames(protein_diff_exp_long,c("variable","value"),c("condition","protIntensity"))
protein_diff_exp_long[,replicate:=0]

protein_diff_exp_long[protein_id=="Q6P2Q9"]

data_yansheng <- fread("protQuant_yansheng.tsv")
data_yansheng <- subset(data_yansheng,select=names(protein_diff_exp_long))
data_yansheng[,condition:=ifelse(condition=="Control","minus","plus")]

data_all <- rbind(data_yansheng,protein_diff_exp_long)
data_all[,replicate := paste0("rep",replicate)]
data_all$condition <- factor(data_all$condition,levels = c("minus","plus"))

data_tech_rep <- subset(data_all,replicate %in% c("rep1","rep4"))
data_biol_rep <- subset(data_all,replicate %in% c("rep1","rep2"))
data_sec <- subset(data_all,replicate %in% c("rep1","rep0"))

library(plyr)

pdf("tech_rep.pdf",width=4,height=3)
data_tech_rep.c <- dcast(data_tech_rep,protein_id+condition ~ replicate,value.var="protIntensity")
data_tech_rep.c <- na.omit(data_tech_rep.c, cols=seq_along(data_tech_rep.c), invert=FALSE)
p <- ggplot(data_tech_rep.c,aes(x=log(rep1),y=log(rep4))) + geom_point() +
  facet_grid(. ~ condition) +
  theme_classic()
cors <- ddply(data_tech_rep.c, .(condition), summarise, cor = round(cor(rep1, rep4), 2))
p <- p + geom_text(data=cors, aes(label=paste("r=", cor, sep="")), x=min(log(data_tech_rep.c$rep1))+1, y=max(log(data_tech_rep.c$rep4)))
print(p)
dev.off()

pdf("biol_rep.pdf",width=4,height=3)
data_biol_rep.c <- dcast(data_biol_rep,protein_id+condition ~ replicate,value.var="protIntensity")
data_biol_rep.c <- na.omit(data_biol_rep.c, cols=seq_along(data_biol_rep.c), invert=FALSE)
p <- ggplot(data_biol_rep.c,aes(x=log(rep1),y=log(rep2))) + geom_point() +
  facet_grid(. ~ condition) +
  theme_classic()
cors <- ddply(data_biol_rep.c, .(condition), summarise, cor = round(cor(rep1, rep2), 2))
p <- p + geom_text(data=cors, aes(label=paste("r=", cor, sep="")), x=min(log(data_biol_rep.c$rep1))+1, y=max(log(data_biol_rep.c$rep2)))
print(p)
dev.off()

pdf("sec_rep.pdf",width=4,height=3)
data_sec.c <- dcast(data_sec,protein_id+condition ~ replicate,value.var="protIntensity")
data_sec.c <- na.omit(data_sec.c, cols=seq_along(data_sec.c), invert=FALSE)
p <- ggplot(data_sec.c,aes(x=log(rep1),y=log(rep0))) + geom_point() +
  facet_grid(. ~ condition) +
  theme_classic() 
cors <- ddply(data_sec.c, .(condition), summarise, cor = round(cor(rep1, rep0), 2))
p <- p + geom_text(data=cors, aes(label=paste("r=", cor, sep="")), x=min(log(data_sec.c$rep1))+1, y=max(log(data_sec.c$rep0)))
print(p)
dev.off()

