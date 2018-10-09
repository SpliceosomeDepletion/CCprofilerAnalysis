library(devtools)
#setwd("~/mysonas/CCprofiler")
#load_all()
setwd("~/mysonas/PRPF8/analysis/output")
install_github("CCprofiler/CCprofiler", ref = "DA_module")
library(CCprofiler)

david_diffExpressed <- fread("enrichment_diffProteins.txt")
names(david_diffExpressed) <- gsub(" ", "_", names(david_diffExpressed))
david_diffExpressed_sig <- subset(david_diffExpressed, Category=="GOTERM_BP_DIRECT" & Fold_Enrichment >= 1 & Bonferroni <= 0.05)
david_diffExpressed_sig[, neglog10Bonferroni := -log10(Bonferroni)]
david_diffExpressed_sig <- david_diffExpressed_sig[order(Fold_Enrichment)]
cat <- david_diffExpressed_sig$Term
david_diffExpressed_sig$Term <- factor(david_diffExpressed_sig$Term, levels=cat)
david_diffExpressed_sig[,name:=gsub("^GO:.*~","",Term)]


q <- ggplot(data=david_diffExpressed_sig,aes(x=name,y=Fold_Enrichment, fill=neglog10Bonferroni)) +
  geom_bar(stat = "identity") +
  theme_classic() +
  theme(axis.text.y  = element_text(angle = 30, hjust = 1)) +
  labs(fill='-log10BF',x='GO Biological Process',y='Fold Enrichment') +
  coord_flip() +
  geom_text(aes(label=round(Fold_Enrichment,digits=1)), hjust=1.2,color="white")

pdf("david_diffExpressed_5percent.pdf",height=6,width=6)
plot(q)
dev.off()


david_associationShift <- fread("enrichment_shiftProteins.txt")
names(david_associationShift) <- gsub(" ", "_", names(david_associationShift))
david_associationShift_sig <- subset(david_associationShift, Category=="GOTERM_BP_DIRECT" & Fold_Enrichment >= 1 & Bonferroni <= 0.05)
david_associationShift_sig[, neglog10Bonferroni := -log10(Bonferroni)]
david_associationShift_sig <- david_associationShift_sig[order(Fold_Enrichment)]
cat <- david_associationShift_sig$Term
david_associationShift_sig$Term <- factor(david_associationShift_sig$Term, levels=cat)
david_associationShift_sig[,name:=gsub("^GO:.*~","",Term)]


q <- ggplot(data=david_associationShift_sig,aes(x=name,y=Fold_Enrichment, fill=neglog10Bonferroni)) +
  geom_bar(stat = "identity") +
  theme_classic() +
  theme(axis.text.y  = element_text(angle = 30, hjust = 1)) +
  labs(fill='-log10BF',x='GO Biological Process',y='Fold Enrichment') +
  coord_flip() +
  geom_text(aes(label=round(Fold_Enrichment,digits=1)), hjust=1.2,color="white")

pdf("david_associationShift_5percent.pdf",height=6,width=6)
plot(q)
dev.off()


david_associationShift <- fread("enrichment_shiftProteins.txt")
names(david_associationShift) <- gsub(" ", "_", names(david_associationShift))
david_associationShift_sig <- subset(david_associationShift, Category=="UP_KEYWORDS" & Fold_Enrichment >= 1 & Bonferroni <= 0.05)
david_associationShift_sig[, neglog10Bonferroni := -log10(Bonferroni)]
david_associationShift_sig <- david_associationShift_sig[order(Fold_Enrichment)]
cat <- david_associationShift_sig$Term
david_associationShift_sig$Term <- factor(david_associationShift_sig$Term, levels=cat)
#david_associationShift_sig[,name:=gsub("^GO:.*~","",Term)]


q <- ggplot(data=david_associationShift_sig,aes(x=Term,y=Fold_Enrichment, fill=neglog10Bonferroni)) +
  geom_bar(stat = "identity") +
  theme_classic() +
  theme(axis.text.y  = element_text(angle = 30, hjust = 1)) +
  labs(fill='-log10BF',x='UP_KEYWORDS',y='Fold Enrichment') +
  coord_flip() +
  geom_text(aes(label=round(Fold_Enrichment,digits=1)), hjust=1.2,color="white")

pdf("david_associationShift_5percent_UP_KEYWORDS.pdf",height=6,width=6)
plot(q)
dev.off()
