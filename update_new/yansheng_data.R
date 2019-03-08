library(devtools)
setwd("~/Desktop/CCprofiler/CCprofiler")
load_all()
setwd("~/Desktop/PRPF8data/")

data_old <- fread("feature_alignment_OpenSWATHresult.tsv")
data_old[,file := gsub("/scratch/.*.tmpdir/","",filename)]
data_old[,file := gsub(".mzXML.gz","",file)]
# remove requant
data_old <- subset(data_old,m_score < 2)

ann <- fread("annotation.csv",col.names = c("file","expl","condition","replicate"))

data_annotated <- merge(data_old,ann,by="file",all.x=T) 

data_annotated <- subset(data_annotated,select=c("Sequence","FullPeptideName","ProteinName","Intensity","condition","replicate"))
# summarize charge states
data_annotated[,Intensity:=sum(Intensity),by=c("Sequence","FullPeptideName","ProteinName","condition","replicate")]
data_annotated <- unique(data_annotated)

data_annotated[,prot := gsub("[^/]*/(.*)", "\\1",ProteinName)]
split <- strsplit(data_annotated$prot, split = "/")
newprot <- lapply(split, function(x){gsub("_.*","",x)})
newprot <- lapply(newprot,unique)
newprot <- lapply(newprot,function(x){paste(c(length(x),x), collapse = "/")})
data_annotated[,gene:=newprot]

data_genotypic <- data_annotated[grep("^1/",gene,invert=F)]
data_genotypic[,gene:=gsub("1/","",gene)]

#' ## Map genes to UniProt identifiers
#' Import biomaRt:
library(biomaRt)
ensembl <- useMart("ENSEMBL_MART_ENSEMBL", dataset="hsapiens_gene_ensembl", host = "jul2016.archive.ensembl.org")
ensembl <- useDataset("hsapiens_gene_ensembl",mart=ensembl)
#' Get genes from traces list and convert to ensembl format:
genes <- unique(data_genotypic$gene)
#genes <- gsub("1/","",genes)
#' Create biomaRt mapping table:
bm <- getBM(attributes = c("ensembl_gene_id", "uniprot_swissprot"),
            filters = "ensembl_gene_id",
            values = genes,
            mart = ensembl)
#' Format mapping table for compatability with traces list:
bm <- as.data.table(bm)
bm <- subset(bm, uniprot_swissprot!="")
#bm[,gene_id := paste0("1/",ensembl_gene_id)]
#bm[,ensembl_gene_id := NULL]
setnames(bm,"ensembl_gene_id","gene")
bm[,protein_id := uniprot_swissprot]
#' Keep only one uniprot id per gene:
bm <- bm[!duplicated(bm$gene),]

data_genotypic <- merge(data_genotypic,bm,by="gene",all.x=T,all.y=F)
data_genotypic <- data_genotypic[!is.na(protein_id)]

plot_data <- subset(data_genotypic,protein_id=="Q6P2Q9")
ggplot(plot_data,aes(x=Sequence,y=Intensity,fill=condition,group=condition)) + 
  geom_bar(stat="identity",position=position_dodge()) +
  facet_grid(replicate ~ .)


quant_subset <- subset(data_genotypic,select=c("Intensity","Sequence","protein_id"))
quant_subset[,medianInt:=median(Intensity),by=c("Sequence","protein_id")]
quant_subset[,count := .N,by=c("Sequence","protein_id")]
quant_subset <- unique(quant_subset[,Intensity:=NULL])
quant_subset[,rankCount:=rank(-count,ties.method = "min"),by=c("protein_id")]
quant_subset <- subset(quant_subset,rankCount %in% c(1,2))
quant_subset[,rankInt:=rank(-medianInt,ties.method = "first"),by="protein_id"]
quant_subset <- subset(quant_subset,rankInt %in% c(1,2))
quant_peptides <- unique(quant_subset$Sequence)

data_prot <- subset(data_genotypic,Sequence %in% quant_peptides)
data_prot[,protIntensity := sum(Intensity),by=c("protein_id","condition","replicate")]
data_prot[,quant_peptides := paste0(Sequence,collapse=";"),by=c("protein_id","condition","replicate")]
data_prot[,n_quant_peptides := length(unique(Sequence)),by=c("protein_id","condition","replicate")]
data_prot <- subset(data_prot,select=c("protein_id","condition","replicate","protIntensity","quant_peptides","n_quant_peptides"))
data_prot <- unique(data_prot)

plot_data <- subset(data_prot,protein_id=="Q6P2Q9")
ggplot(plot_data,aes(x=replicate,y=protIntensity,fill=condition)) + 
  geom_bar(stat="identity",position=position_dodge()) +
  theme_classic()

write.table(data_prot,"protQuant_yansheng.tsv",sep="\t",quote=F,col.names=T,row.names = F)

