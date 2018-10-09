library(devtools)
#setwd("~/mysonas/CCprofiler")
#load_all()
setwd("~/mysonas/PRPF8/analysis/output")
install_github("CCprofiler/CCprofiler", ref = "DA_module")
library(CCprofiler)

david_moreAssembled <- fread("DAVID_moreAssembled.txt")
names(david_moreAssembled) <- gsub(" ", "_", names(david_moreAssembled))
david_moreAssembled_sig <- subset(david_moreAssembled, Category=="UP_KEYWORDS" & Fold_Enrichment >= 1 & Bonferroni <= 0.05)
david_moreAssembled_sig[, neglog10Bonferroni := -log10(Bonferroni)]
david_moreAssembled_sig <- david_moreAssembled_sig[order(Fold_Enrichment)]
cat <- david_moreAssembled_sig$Term
david_moreAssembled_sig$Term <- factor(david_moreAssembled_sig$Term, levels=cat)

q <- ggplot(data=david_moreAssembled_sig,aes(x=Term,y=Fold_Enrichment, fill=neglog10Bonferroni)) +
  geom_bar(stat = "identity") +
  theme_classic() +
  theme(axis.text.y  = element_text(angle = 30, hjust = 1), legend.position=c(0.8,0.2)) +
  labs(fill='-log10BF',x='UniprotKB keywords',y='Fold Enrichment') +
  coord_flip() +
  geom_text(aes(label=round(Fold_Enrichment,digits=1)), hjust=1.2,color="white")

pdf("david_moreAssembled_5percent.pdf",height=6,width=6)
plot(q)
dev.off()

#############

david_lessAssembled <- fread("DAVID_lessAssembled.txt")
names(david_lessAssembled) <- gsub(" ", "_", names(david_lessAssembled))
david_lessAssembled_sig <- subset(david_lessAssembled, Category=="UP_KEYWORDS" &Fold_Enrichment >= 1 & Bonferroni <= 0.05)
david_lessAssembled_sig[, neglog10Bonferroni := -log10(Bonferroni)]
david_lessAssembled_sig <- david_lessAssembled_sig[order(Fold_Enrichment)]
cat <- david_lessAssembled_sig$Term
david_lessAssembled_sig$Term <- factor(david_lessAssembled_sig$Term, levels=cat)

library(ggplot2)

q <- ggplot(data=david_lessAssembled_sig,aes(x=Term,y=Fold_Enrichment, fill=neglog10Bonferroni)) +
  geom_bar(stat = "identity") +
  theme_classic() +
  theme(axis.text.y  = element_text(angle = 50, hjust = 1), legend.position=c(0.8,0.2)) +
  labs(fill='-log10BF',x='UniprotKB keywords',y='Fold Enrichment') +
  coord_flip() +
  geom_text(aes(label=round(Fold_Enrichment,digits=1)), hjust=1.2,color="white")

pdf("david_lessAssembled_5percent.pdf",height=6,width=6)
plot(q)
dev.off()


#############

david_noChange <- fread("DAVID_noAssemblyChange.txt")
names(david_noChange) <- gsub(" ", "_", names(david_noChange))
david_noChange_sig <- subset(david_noChange, Category=="UP_KEYWORDS" & Fold_Enrichment >= 1 & Bonferroni <= 0.05)
david_noChange_sig[, neglog10Bonferroni := -log10(Bonferroni)]
david_noChange_sig <- david_noChange_sig[order(Fold_Enrichment)]
cat <- david_noChange_sig$Term
david_noChange_sig$Term <- factor(david_noChange_sig$Term, levels=cat)

library(ggplot2)

q <- ggplot(data=david_noChange_sig,aes(x=Term,y=Fold_Enrichment, fill=neglog10Bonferroni)) +
  geom_bar(stat = "identity") +
  theme_classic() +
  theme(axis.text.y  = element_text(angle = 50, hjust = 1), legend.position=c(0.8,0.2)) +
  labs(fill='-log10BF',x='UniprotKB keywords',y='Fold Enrichment') +
  coord_flip() +
  geom_text(aes(label=round(Fold_Enrichment,digits=1)), hjust=1.2,color="white")

pdf("david_noChange_5percent.pdf",height=6,width=6)
plot(q)
dev.off()

#############

david_lessAssembled <- fread("DAVID_lessAssembled.txt")
names(david_lessAssembled) <- gsub(" ", "_", names(david_lessAssembled))
david_lessAssembled_sig <- subset(david_lessAssembled, Category=="UP_KEYWORDS" & Fold_Enrichment >= 1)
david_lessAssembled_sig[, neglog10Bonferroni := -log10(Bonferroni)]
david_lessAssembled_sig <- david_lessAssembled_sig[order(Fold_Enrichment)]
cat <- david_lessAssembled_sig$Term
david_lessAssembled_sig$Term <- factor(david_lessAssembled_sig$Term, levels=cat)

library(ggplot2)

q <- ggplot(data=david_lessAssembled_sig,aes(x=Term,y=Fold_Enrichment, fill=neglog10Bonferroni)) +
  geom_bar(stat = "identity") +
  theme_classic() +
  theme(axis.text.y  = element_text(angle = 50, hjust = 1), legend.position=c(0.8,0.2)) +
  labs(fill='-log10BF',x='UniprotKB keywords',y='Fold Enrichment') +
  coord_flip() +
  geom_text(aes(label=round(Fold_Enrichment,digits=1)), hjust=1.2,color="white")

pdf("david_lessAssembled_NOpercent.pdf",height=6,width=6)
plot(q)
dev.off()

########

david_moreAssembled <- fread("DAVID_moreAssembled.txt")
names(david_moreAssembled) <- gsub(" ", "_", names(david_moreAssembled))
david_moreAssembled_sig <- subset(david_moreAssembled, Category %in% c("GOTERM_BP_DIRECT","GOTERM_MF_DIRECT") & Fold_Enrichment >= 1 & Bonferroni <= 0.05)
david_moreAssembled_sig[, neglog10Bonferroni := -log10(Bonferroni)]
david_moreAssembled_sig <- david_moreAssembled_sig[order(Fold_Enrichment)]
cat <- david_moreAssembled_sig$Term
david_moreAssembled_sig$Term <- factor(david_moreAssembled_sig$Term, levels=cat)

q <- ggplot(data=david_moreAssembled_sig,aes(x=Term,y=Fold_Enrichment, fill=neglog10Bonferroni)) +
  geom_bar(stat = "identity") +
  theme_classic() +
  theme(axis.text.y  = element_text(angle = 30, hjust = 1), legend.position=c(0.8,0.2)) +
  labs(fill='-log10BF',x='GO Biological Process / Molecular Function',y='Fold Enrichment') +
  coord_flip() +
  geom_text(aes(label=round(Fold_Enrichment,digits=1)), hjust=1.2,color="white")

pdf("david_moreAssembled_5percent_GOTERM_BP_MF_DIRECT.pdf",height=6,width=6)
plot(q)
dev.off()
