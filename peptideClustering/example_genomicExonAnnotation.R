library(devtools)
setwd("~/mysonas/CCprofiler")
load_all()
setwd("~/mysonas/PRPF8/analysis/output")
#install_github("CCprofiler/CCprofiler", ref = "proteoformLocationMapping")
#library(CCprofiler)

traces_annotated <- readRDS("pepTraces_proteoformMulti_combined_genomAnnotation_pvals.rds")
map_table <- readRDS("../../data/proteogenomics/HeLa_protgen_combined_1FPKM_header_mapping_v2.rda")

map_table$IsoformId <- gsub("\\|.*", "", gsub(">.*?\\|","", map_table$header))
isomap <- unique(map_table[,.(LeadingIsoform=IsoformId, LeadingEnsemblProteinSp=protein)])
isomap <- isomap[!duplicated(LeadingIsoform)] # Take the first protein id if the seq is the same

traces_annotated$trace_annotation <- merge(traces_annotated$trace_annotation, isomap,by="LeadingIsoform")
traces_annotated$trace_annotation[, LeadingEnsemblProtein := gsub("_.*", "", LeadingEnsemblProteinSp)]
traces_annotated$trace_annotation <- traces_annotated$trace_annotation[order(id)]

#' Now the genomic position is retrieved with the ensembldb package.
library(EnsDb.Hsapiens.v86)
ensdb <- EnsDb.Hsapiens.v86

traces_coordinates <- annotateGenomicCoordinates(traces_annotated, db=ensdb, proteinIdCol="LeadingEnsemblProtein", verbose=T)

saveRDS(traces_coordinates,"pepTraces_proteoformMulti_combined_genomAnnotation_pvals_GenomicCoordinates.rds")
