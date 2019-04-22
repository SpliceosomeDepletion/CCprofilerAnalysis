#' ## Load data
pepTracesList_filtered <- readRDS("pepTracesList_filtered.rda")
map_table <- readRDS("../../PRPF8_protgen_combined_1FPKM_header_mapping_v2.rda")

# just for testing
traces_withAnnotation <- readRDS("../output_wrongAnnotation/pepTraces_proteoformMulti_combined_genomAnnotation_pvals_GenomicCoordinates.rds")
pepsWithAnnotation <- unique(names(traces_withAnnotation$genomic_coord))

allPeps <- unique(unlist(lapply(pepTracesList_filtered, function(p){p$traces$id})))

pepsNoAnn <- allPeps[! allPeps %in% pepsWithAnnotation]

pepTracesList_filtered <- subset(pepTracesList_filtered, pepsNoAnn)

map_table$IsoformId <- gsub("\\|.*", "", gsub(">.*?\\|","", map_table$header))
isomap <- unique(map_table[,.(LeadingIsoform=IsoformId, LeadingEnsemblProteinSp=protein)])
isomap <- isomap[!duplicated(LeadingIsoform)] # Take the first protein id if the seq is the same

pepTracesList_filtered_sum <- integrateTraceIntensities(pepTracesList_filtered,
                                                        design_matrix = NULL,
                                                        integrate_within = NULL,
                                                        aggr_fun = "sum")

pepTracesList_filtered_sum$trace_annotation <- merge(pepTracesList_filtered_sum$trace_annotation, isomap,by="LeadingIsoform")
pepTracesList_filtered_sum$trace_annotation[, LeadingEnsemblProtein := gsub("_.*", "", LeadingEnsemblProteinSp)]
pepTracesList_filtered_sum$trace_annotation <- pepTracesList_filtered_sum$trace_annotation[order(id)]


#' ##### Load genomic information from the ensembldb package
library(EnsDb.Hsapiens.v86)
ensdb <- EnsDb.Hsapiens.v86

traces_coordinates <- annotateGenomicCoordinates(pepTracesList_filtered_sum, db=ensdb, proteinIdCol="LeadingEnsemblProtein", verbose=T)


#################
#################

#' ## Load data
pepTracesList_filtered <- readRDS("pepTracesList_filtered.rda")
map_table <- readRDS("../../PRPF8_protgen_combined_1FPKM_header_mapping_v2.rda")

# just for testing
traces_withAnnotation <- readRDS("../output_wrongAnnotation/pepTraces_proteoformMulti_combined_genomAnnotation_pvals_GenomicCoordinates.rds")
protsWithAnnotation <- unique(traces_withAnnotation$trace_annotation$protein_id)

allProts <- unique(unlist(lapply(pepTracesList_filtered, function(p){p$trace_annotation$protein_id})))

protsNoAnn <- allProts[! allProts %in% protsWithAnnotation]

pepTracesList_filtered <- subset(pepTracesList_filtered, protsNoAnn, trace_subset_type="protein_id")

map_table$IsoformId <- gsub("\\|.*", "", gsub(">.*?\\|","", map_table$header))
isomap <- unique(map_table[,.(LeadingIsoform=IsoformId, LeadingEnsemblProteinSp=protein)])
isomap <- isomap[!duplicated(LeadingIsoform)] # Take the first protein id if the seq is the same

pepTracesList_filtered_sum <- integrateTraceIntensities(pepTracesList_filtered,
                                                        design_matrix = NULL,
                                                        integrate_within = NULL,
                                                        aggr_fun = "sum")

pepTracesList_filtered_sum$trace_annotation <- merge(pepTracesList_filtered_sum$trace_annotation, isomap,by="LeadingIsoform")
pepTracesList_filtered_sum$trace_annotation[, LeadingEnsemblProtein := gsub("_.*", "", LeadingEnsemblProteinSp)]
pepTracesList_filtered_sum$trace_annotation <- pepTracesList_filtered_sum$trace_annotation[order(id)]


#' ##### Load genomic information from the ensembldb package
library(EnsDb.Hsapiens.v86)
ensdb <- EnsDb.Hsapiens.v86

traces_coordinates <- annotateGenomicCoordinates(pepTracesList_filtered_sum, db=ensdb, proteinIdCol="LeadingEnsemblProtein", verbose=T)
