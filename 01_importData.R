library(devtools)
#setwd("~/mysonas/CCprofiler")
#load_all()
setwd("~/mysonas/PRPF8/analysis/output")
install_github("CCprofiler/CCprofiler", ref = "DA_module")
library(CCprofiler)

#data <- fread("../../data/DIAsearch/output/aligned.csv")
data <- fread("../../data/DIAsearch/output/aligned_filtered.csv")
data[, filename := gsub("/cluster/home/ibludau/mysonas/html/openBIS/.*\\/|\\.mzXML\\.gz", "", filename)]
setkey(data, "filename")

# remove decoys
data <- subset(data, decoy==0)
data <- data[grep("DECOY",data$ProteinName,invert=T)]
# remove reverse from data
# these were in there at some oint because of missing filtering step in spectrast
# idx_rev <- grep("reverse",data$ProteinName)
# if(length(idx_rev) > 0) { data <- data[-idx_rev] }

# remone non genotypic IDs
# all non genotypic IDs start with 2/ or higher >> genotypic ones do not have prefix
data <- data[grep("^ENSG",data$ProteinName,invert=F)]

# maybe at some point add annotation from RNAseq
isoform_annotation <- fread("../../data/DDAsearch/mapping_genotypic.tsv",sep="\t",header=F)
isoform_annotation <- isoform_annotation[grep("^1/",V3)]
isoform_annotation[,gene:=gsub("1/","",V3)]
isoform_subset <- subset(isoform_annotation, select=c("V1","V2","gene"))

# generate naked sequence for mapping
isoform_ex <- copy(isoform_subset)
isoform_ex$V1 <- gsub("\\(.*?\\)", "", isoform_ex$V1)
isoform_ex$V1 <- gsub("\\[.*?\\]", "", isoform_ex$V1)
isoform_ex$V1 <- gsub("\\.", "", isoform_ex$V1)
isoform_ex$V1 <- gsub("\\/.*", "", isoform_ex$V1)

names(isoform_ex) <- c("Sequence","annotation","ProteinName")
isoform_ex_u <- unique(isoform_ex,by=c("Sequence","ProteinName"))

#dups <- isoform_ex_u$Sequence[duplicated(isoform_ex_u$Sequence)]

data.m <- merge(data,isoform_ex_u,all.x=T,all.y=F,by=c("Sequence","ProteinName"))
data.m[,ProteinName:=annotation]
data.m[,annotation:=NULL]

ann <- fread("../../fraction_annotation.csv")

# for some reason important @TODO change this
setkey(data.m, "filename")

traces_list <- importMultipleCondiionsFromOpenSWATH(data=data.m,annotation_table=ann,
  rm_requantified=TRUE,
  rm_decoys = FALSE,
  rm_nonProteotypic = FALSE,
  MS1Quant=FALSE,
  proteogenomicsWF=TRUE,
  verbose=TRUE)

rm(data,data.m)
gc()

summary(traces_list)

# add gene_id
#traces_list$minus$trace_annotation[,gene_id := protein_id]
#traces_list$minus$trace_annotation[,protein_id := NULL]
#traces_list$plus$trace_annotation[,gene_id := protein_id]
#traces_list$plus$trace_annotation[,protein_id := NULL]

# map genes to uniprot
library(biomaRt)
#ensembl=useMart("ensembl")
#ensembl = useDataset("hsapiens_gene_ensembl",mart=ensembl) # run 2 times then it works
ensembl <- useMart("ENSEMBL_MART_ENSEMBL", dataset="hsapiens_gene_ensembl", host = "jul2016.archive.ensembl.org")
ensembl = useDataset("hsapiens_gene_ensembl",mart=ensembl)

genes <- c(traces_list$minus$trace_annotation$gene_id,traces_list$plus$trace_annotation$gene_id)
genes <- unique(genes)
genes <- gsub("1/","",genes)

bm <- getBM(attributes = c("ensembl_gene_id", "uniprot_swissprot"),
            filters = "ensembl_gene_id",
            values = genes,
            mart = ensembl)
bm <- as.data.table(bm)
bm <- subset(bm, uniprot_swissprot!="")
bm[,gene_id := paste0("1/",ensembl_gene_id)]
bm[,ensembl_gene_id := NULL]
bm[,protein_id := uniprot_swissprot]
## keep only one uniprot id per gene
bm <- bm[!duplicated(bm$gene_id),]

traces_list <- annotateTraces(traces_list,
                                bm,
                                traces_id_column = "gene_id",
                                trace_annotation_id_column = "gene_id")

# how many genes do not have uniprot accession available
length(unique(traces_list$minus$trace_annotation$gene_id)) - length(unique(traces_list$minus$trace_annotation$protein_id))

# subset traces to genes that have uniprot accession available
uniprot_available <- unique(c(traces_list$minus$trace_annotation$protein_id,traces_list$plus$trace_annotation$protein_id))
uniprot_available <- uniprot_available[!is.na(uniprot_available)]
traces_list <- subset(traces_list,trace_subset_ids=uniprot_available,trace_subset_type="protein_id")

## Molecular weight calibration
calibrationTable <- fread("../../calibrationTable_bigColumn.txt")
calibration = calibrateMW(calibrationTable,
                          PDF=TRUE,
                          plot=TRUE)
traces_list <- annotateMolecularWeight(traces_list, calibration)

## trace annotation
traces_list <- annotateTraces(traces_list, exampleTraceAnnotation, traces_id_column = "uniprot_swissprot")

## update traces with additional metrics
traces_list <- updateTraces(traces_list)

summary(traces_list)

# design matrix
design_matrix <- data.table(Sample_name = c("minus","plus"),
                            Condition = c("minus","plus"),
                            Replicate = c(1,1))

saveRDS(calibration,"calibration.rds")
saveRDS(design_matrix, "design_matrix.rda")
saveRDS(traces_list, "pepTracesRaw.rda")
