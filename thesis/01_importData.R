#' ## Read input data and annotate with RNAseq information
#' Read TRIC output file:
data <- fread("../../data/DIAsearch/output/aligned_filtered.csv")
data[, filename := gsub(
  "/cluster/home/ibludau/mysonas/html/openBIS/.*\\/|\\.mzXML\\.gz",
  "", filename)]
setkey(data, "filename")

#' Remove decoys from the TRIC output:
data <- subset(data, decoy==0)
data <- data[grep("DECOY",data$ProteinName,invert=T)]

#' Remove non-genotypic IDs:
# All non genotypic IDs start with 2/ or higher.
# Genotypic IDs do not have a prefix and start with ENSG.
# This also removes all cRAP and iRT proteins.
data <- data[grep("^ENSG",data$ProteinName,invert=F)]

#' Read RNAseq information that was previously removed
#' during library generation:
isoform_annotation <- fread("../../data/DDAsearch/mapping_genotypic.tsv",
                            sep="\t",header=F)
isoform_annotation <- isoform_annotation[grep("^1/",V3)]
isoform_annotation[,gene:=gsub("1/","",V3)]
isoform_subset <- subset(isoform_annotation, select=c("V1","V2","gene"))

#' Generate naked peptide sequence for mapping:
isoform_ex <- copy(isoform_subset)
isoform_ex$V1 <- gsub("\\(.*?\\)", "", isoform_ex$V1)
isoform_ex$V1 <- gsub("\\[.*?\\]", "", isoform_ex$V1)
isoform_ex$V1 <- gsub("\\.", "", isoform_ex$V1)
isoform_ex$V1 <- gsub("\\/.*", "", isoform_ex$V1)

#' Subset annotation table for merging:
names(isoform_ex) <- c("Sequence","annotation","ProteinName")
isoform_ex_u <- unique(isoform_ex,by=c("Sequence","ProteinName"))

#' Merge TRIC output with RNAseq annotation info:
data.m <- merge(data,isoform_ex_u,all.x=T,all.y=F,
                by=c("Sequence","ProteinName"))
data.m[,ProteinName:=annotation]
data.m[,annotation:=NULL]
setkey(data.m, "filename")

saveRDS(data.m,"inputData.rds")

#' ## Convert data to CCprofiler tracesList object
#' Read SEC fraction annotation table:
ann <- fread("../../fraction_annotation.csv")
ann[, condition := gsub("minus", "control", condition)]
ann[, condition := gsub("plus", "depleted", condition)]
ann[, sample := gsub("minus", "control", sample)]
ann[, sample := gsub("plus", "depleted", sample)]

#' Save non-SEC samples
non_sec_samples <- subset(data.m, ! filename %in% ann$filename)
saveRDS(non_sec_samples,"non_sec_samples.rds")

#' Craete tracesList object:
traces_list <- importMultipleCondiionsFromOpenSWATH(data=data.m,
  annotation_table=ann,
  rm_requantified=TRUE,
  rm_decoys = FALSE,
  rm_nonProteotypic = FALSE,
  MS1Quant=FALSE,
  proteogenomicsWF=TRUE,
  verbose=F)
#' Clean memory:
rm(data,data.m)
gc()
#' Inspect traces list:
summary(traces_list)
#' Remove prefix from gene_id
traces_list$control$trace_annotation$gene_id <- gsub("1/","",traces_list$control$trace_annotation$gene_id)
traces_list$depleted$trace_annotation$gene_id <- gsub("1/","",traces_list$depleted$trace_annotation$gene_id)

#' ## Map genes to UniProt identifiers
#' Import biomaRt:
library(biomaRt)
ensembl <- useMart("ENSEMBL_MART_ENSEMBL", dataset="hsapiens_gene_ensembl", host = "www.ensembl.org")
ensembl <- useDataset("hsapiens_gene_ensembl",mart=ensembl)
#' Get genes from traces list and convert to ensembl format:
genes <- c(traces_list$control$trace_annotation$gene_id,traces_list$depleted$trace_annotation$gene_id)
genes <- unique(genes)
#genes <- gsub("1/","",genes)
#' Create biomaRt mapping table:
bm <- getBM(attributes = c("ensembl_gene_id", "uniprotswissprot"),
            filters = "ensembl_gene_id",
            values = genes,
            mart = ensembl)
#' Format mapping table for compatability with traces list:
bm <- as.data.table(bm)
length(which(!(unique(bm[uniprotswissprot==""]$ensembl_gene_id) %in% unique(bm[uniprotswissprot!=""]$ensembl_gene_id))))
bm <- subset(bm, uniprotswissprot!="")
#bm[,gene_id := paste0("1/",ensembl_gene_id)]
#bm[,ensembl_gene_id := NULL]
setnames(bm,"ensembl_gene_id","gene_id")
bm[,protein_id := uniprotswissprot]
#' Keep only one uniprot id per gene:
length(unique(bm[duplicated(bm$gene_id)]$gene_id))
bm <- bm[!duplicated(bm$gene_id),]
#' Annotate traces with biomaRt info:
traces_list <- annotateTraces(traces_list,
                                bm,
                                traces_id_column = "gene_id",
                                trace_annotation_id_column = "gene_id")

#' How many genes do not have uniprot accession available?
length(unique(traces_list$control$trace_annotation$gene_id)) - length(unique(traces_list$control$trace_annotation$protein_id))

#' Subset traces to genes that have uniprot accession available:
uniprot_available <- unique(c(traces_list$control$trace_annotation$protein_id,traces_list$depleted$trace_annotation$protein_id))
uniprot_available <- uniprot_available[!is.na(uniprot_available)]
traces_list <- subset(traces_list,trace_subset_ids=uniprot_available,trace_subset_type="protein_id")

#' ## Molecular weight calibration
#' Read calibration table:
# calibrationTable <- fread("../../calibrationTable_bigColumn.txt")
calibrationTable <- fread("../../../html/HeLa_Kyoto_CCL2/ccl2kyoto_sec_mw_spectronaut.csv")
setnames(calibrationTable,c("run_id","sec_id"),c("filename","fraction_number"))
calibrationTable[, sample := paste(c(condition_id, replicate_id), collapse = "_"), by=c("filename","fraction_number")]
calibrationTable <- unique(subset(calibrationTable, select=c("sec_mw","fraction_number")))
setnames(calibrationTable,c("sec_mw","fraction_number"),c("std_weights_kDa","std_elu_fractions"))
calibrationTable[, std_weights_kDa := std_weights_kDa/1000]

#' Perform molecular weight calibration:
calibration = calibrateMW(calibrationTable,
                          PDF=TRUE,
                          plot=TRUE)
#' Annotate fractions in tares list with according molecular weights:
traces_list <- annotateMolecularWeight(traces_list, calibration)

#' ## Annotate traces
uniprotAnnotation <- fread("../../uniprot_swissprot_annotation.tab")
#' Annotate traces with information from UniProt:
traces_list <- annotateTraces(traces_list,
                              uniprotAnnotation,
                              traces_id_column = "uniprotswissprot")
#' Update traces with additional metrics for each fraction:
traces_list <- updateTraces(traces_list)
#' Inspect traces list:
summary(traces_list)

#' ## Create design matrix
design_matrix <- data.table(Sample_name = c("control","depleted"),
                            Condition = c("control","depleted"),
                            Replicate = c(1,1))

#' ## Save objects
saveRDS(calibration,"calibration.rds")
saveRDS(design_matrix, "design_matrix.rda")
saveRDS(traces_list, "pepTracesRaw.rda")

#' ## QC plots
test_proteins <- c("Q6P2Q9")
pepTest <- subset(traces_list, trace_subset_ids = test_proteins, trace_subset_type = "protein_id")
plot(pepTest, design_matrix = design_matrix, log = FALSE, legend = T, PDF = TRUE,
     name = paste0("RAW_PeptideTraces_","PRPF8"), isoformAnnotation = F, plot = TRUE, highlight = NULL, highlight_col = NULL, colour_by = "Entry_name")

test_proteins <- c("Q15029")
pepTest <- subset(traces_list, trace_subset_ids = test_proteins, trace_subset_type = "protein_id")
plot(pepTest, design_matrix = design_matrix, log = FALSE, legend = T, PDF = TRUE,
     name = paste0("RAW_PeptideTraces_","EFTUD2"), isoformAnnotation = F, plot = TRUE, highlight = NULL, highlight_col = NULL, colour_by = "Entry_name")

test_proteins <- c("Q96DI7")
pepTest <- subset(traces_list, trace_subset_ids = test_proteins, trace_subset_type = "protein_id")
plot(pepTest, design_matrix = design_matrix, log = FALSE, legend = T, PDF = TRUE,
     name = paste0("RAW_PeptideTraces_","SNRNP40"), isoformAnnotation = F, plot = TRUE, highlight = NULL, highlight_col = NULL, colour_by = "Entry_name")

test_proteins <- c("O75643")
pepTest <- subset(traces_list, trace_subset_ids = test_proteins, trace_subset_type = "protein_id")
plot(pepTest, design_matrix = design_matrix, log = FALSE, legend = T, PDF = TRUE,
     name = paste0("RAW_PeptideTraces_","SNRNP200"), isoformAnnotation = F, plot = TRUE, highlight = NULL, highlight_col = NULL, colour_by = "Entry_name")


