library(devtools)
#setwd("~/mysonas/CCprofiler")
#load_all()
setwd("~/mysonas/PRPF8/analysis/output")
install_github("CCprofiler/CCprofiler", ref = "DA_module")
library(CCprofiler)

traces_list <- readRDS("pepTracesImpSPCmultipep.rda")
design_matrix <- readRDS("design_matrix.rda")

# count ids
countIDs <- summarize_proteogenomics(traces_list)
getIDs <- getProteogenomicsIds(traces_list)

#isoform <- annotateLeadingIsoform(traces_list)
#isoform_geneLocation <- annotateGenomicCoordinates(isoform)

# plot condition overlaps
library("Vennerable")

venn_genes<-Venn(list("genes_minus"=getIDs$minus$genes,"genes_plus"=getIDs$plus$genes))
print(venn_genes)
pdf("venn_genes_minus_plus.pdf")
  plot(venn_genes, type="circles", doWeights=TRUE)
dev.off()

venn_isotypic<-Venn(list("isotypic_minus"=getIDs$minus$isotypic,"isotypic_plus"=getIDs$plus$isotypic))
print(venn_isotypic)
pdf("venn_isotypic_minus_plus.pdf")
  plot(venn_isotypic, type="circles", doWeights=TRUE)
dev.off()

venn_isoforms<-Venn(list("isoforms_minus"=getIDs$minus$isoforms,"isoforms_plus"=getIDs$plus$isoforms))
print(venn_isoforms)
pdf("venn_isoforms_minus_plus.pdf")
  plot(venn_isoforms, type="circles", doWeights=TRUE)
dev.off()

venn_peptides<-Venn(list("peptides_minus"=getIDs$minus$peptides,"peptides_plus"=getIDs$plus$peptides))
print(venn_peptides)
pdf("venn_peptides_minus_plus.pdf")
  plot(venn_peptides, type="circles", doWeights=TRUE)
dev.off()

genes_uniqueMinus <- length(which(! getIDs$minus$genes %in% getIDs$plus$genes))
genes_uniquePlus <- length(which(! getIDs$plus$genes %in% getIDs$minus$genes))
genes_shared <- length(intersect(getIDs$plus$genes, getIDs$minus$genes))

isoforms_uniqueMinus <- length(which(! getIDs$minus$isoforms %in% getIDs$plus$isoforms))
isoforms_uniquePlus <- length(which(! getIDs$plus$isoforms %in% getIDs$minus$isoforms))
isoforms_shared <- length(intersect(getIDs$plus$isoforms, getIDs$minus$isoforms))

isotypic_uniqueMinus <- length(which(! getIDs$minus$isotypic %in% getIDs$plus$isotypic))
isotypic_uniquePlus <- length(which(! getIDs$plus$isotypic %in% getIDs$minus$isotypic))
isotypic_shared <- length(intersect(getIDs$plus$isotypic, getIDs$minus$isotypic))

peptides_uniqueMinus <- length(which(! getIDs$minus$peptides %in% getIDs$plus$peptides))
peptides_uniquePlus <- length(which(! getIDs$plus$peptides %in% getIDs$minus$peptides))
peptides_shared <- length(intersect(getIDs$plus$peptides, getIDs$minus$peptides))

dt <- data.table(genes=c(genes_shared,genes_uniqueMinus,genes_uniquePlus),
                isoforms=c(isoforms_shared,isoforms_uniqueMinus,isoforms_uniquePlus),
                isotypic=c(isotypic_shared,isotypic_uniqueMinus,isotypic_uniquePlus),
                peptides=c(peptides_shared,peptides_uniqueMinus,peptides_uniquePlus),
                evidence = c("shared","control only","PRPF8 depleted only"))

dt.m <- melt(dt,id.vars="evidence")

pdf("plus_minus_counts.pdf",height=3,width=5.5)
  g <- ggplot(dt.m,aes(x=variable,y=value,fill=evidence)) +
        geom_bar(stat="identity") +
        facet_wrap(~ variable, scales="free",nrow=1) +
        theme_classic() +
        theme(legend.position="bottom") + theme(legend.title = element_blank()) +
        labs(y = "count") +
        theme(axis.title.x=element_blank())
  g
dev.off()

# RNA vs. SWATH comparison
rnaIsoforms <- fread("~/mysonas/PRPF8/data/proteogenomics/rnaIsoforms.txt",sep="\n",header=F)$V1
rnaGenes <- fread("~/mysonas/PRPF8/data/proteogenomics/rnaGenes.txt",sep="\n",header=F)$V1

proteinGenes <- gsub("1/","",unique(c(getIDs$minus$genes,getIDs$plus$genes)))
proteinIsoforms <- unique(c(getIDs$minus$isoforms,getIDs$plus$isoforms))

genes_onlyRNAseq <- length(rnaGenes)
genes_Proteomics <- length(proteinGenes)
isoforms_onlyRNAseq <- length(rnaIsoforms)
isoforms_Proteomics <- length(proteinIsoforms)
genesDB <- data.table(
  count=c(genes_Proteomics,genes_onlyRNAseq),
  evidence=c("SWATH/DIA","RNAseq"))
isoformsDB <- data.table(
  count=c(isoforms_Proteomics,isoforms_onlyRNAseq),
  evidence=c("SWATH/DIA","RNAseq"))

pdf("RnaProt_genes.pdf",height=3,width=2)
  g <- ggplot(genesDB,aes(evidence,count)) +
        geom_bar(stat="identity", fill="steelblue") +
        geom_text(aes(label=count), vjust=1.6, color="white", size=3.5)+
        theme_classic()
  g
dev.off()
pdf("RnaProt_isoforms.pdf",height=3,width=2)
  g <- ggplot(isoformsDB,aes(evidence,count)) +
        geom_bar(stat="identity", fill="steelblue") +
        geom_text(aes(label=count), vjust=1.6, color="white", size=3.5)+
        theme_classic()
  g
dev.off()

venn_RnaProt_isoforms<-Venn(list("isoforms_rna"=rnaIsoforms,"isoforms_protein"=proteinIsoforms))
print(venn_RnaProt_isoforms)
pdf("venn_venn_RnaProt_isoforms.pdf")
  plot(venn_RnaProt_isoforms, type="circles", doWeights=TRUE)
dev.off()

venn_RnaProt_genes<-Venn(list("genes_rna"=rnaGenes,"genes_protein"=proteinGenes))
print(venn_RnaProt_genes)
pdf("venn_venn_RnaProt_genes.pdf")
  plot(venn_RnaProt_genes, type="circles", doWeights=TRUE)
dev.off()
