library(devtools)
#setwd("~/mysonas/CCprofiler")
#load_all()
setwd("~/mysonas/PRPF8/analysis/output")
install_github("CCprofiler/CCprofiler", ref = "DA_module")
library(CCprofiler)

#data <- fread("../../data/DIAsearch/output/aligned.csv")
data <- fread("../../data/DIAsearch/output_old/aligned_filtered.csv")
data[, filename := gsub("/cluster/home/ibludau/mysonas/html/openBIS/.*\\/|\\.mzXML\\.gz", "", filename)]
setkey(data, "filename")
setnames(data, "FullUniModPeptideName", "FullPeptideName")
data <- subset(data, peak_group_rank==1)
data <- subset(data, decoy==0)

ann <- fread("../../fraction_annotation.csv")

traces_list <- importMultipleCondiionsFromOpenSWATH(data=data,annotation_table=ann,
  rm_requantified=TRUE,
  rm_decoys = FALSE,
  rm_nonProteotypic = FALSE,
  MS1Quant=FALSE,
  proteogenomicsWF=TRUE,
  verbose=TRUE)

rm(data)
gc()

summary(traces_list)

# remone non genotypic IDs
genes <- c(traces_list$minus$trace_annotation$gene_id,traces_list$plus$trace_annotation$gene_id)
genes <- unique(genes)
genes_typ_idx <- grep("1/",genes)
genes_typ <- genes[genes_typ_idx]
# non-genotypic entries
length(genes) - length(genes_typ_idx)
# subset traces_list
traces_list <- subset(traces_list,trace_subset_ids=genes_typ,trace_subset_type="gene_id")


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
