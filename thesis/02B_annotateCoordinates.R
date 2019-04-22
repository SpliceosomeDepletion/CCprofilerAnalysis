#' ## Load data
pepTracesList_filtered <- readRDS("pepTracesList_filtered.rda")
map_table <- readRDS("../../PRPF8_protgen_combined_1FPKM_header_mapping_v2.rda")

# just for testing
#traces_withAnnotation <- readRDS("../output_wrongAnnotation/pepTraces_proteoformMulti_combined_genomAnnotation_pvals_GenomicCoordinates.rds")
#pepsWithAnnotation <- unique(names(traces_withAnnotation$genomic_coord))
#pepTracesList_filtered <- subset(pepTracesList_filtered, pepsWithAnnotation)

#' ## Anntate traces with genomic coordinates
#' ##### Combine traces with RNAseq information
#' Format RNAseq mapping table:
map_table$IsoformId <- gsub("\\|.*", "", gsub(">.*?\\|","", map_table$header))
isomap <- unique(map_table[,.(LeadingIsoform=IsoformId, LeadingEnsemblProteinSp=protein)])
isomap <- isomap[!duplicated(LeadingIsoform)] # Take the first protein id if the seq is the same

pepTracesList_filtered_sum <- integrateTraceIntensities(pepTracesList_filtered,
                                                        design_matrix = NULL,
                                                        integrate_within = NULL,
                                                        aggr_fun = "sum")

#' Update trace annotation:
pepTracesList_filtered_sum$trace_annotation <- merge(pepTracesList_filtered_sum$trace_annotation, isomap,by="LeadingIsoform")
pepTracesList_filtered_sum$trace_annotation[, LeadingEnsemblProtein := gsub("_.*", "", LeadingEnsemblProteinSp)]
pepTracesList_filtered_sum$trace_annotation <- pepTracesList_filtered_sum$trace_annotation[order(id)]

#' ##### Load genomic information from the ensembldb package
library(EnsDb.Hsapiens.v86)
ensdb <- EnsDb.Hsapiens.v86

#pepTracesList_filtered_sum <- subset(pepTracesList_filtered_sum,
#                 trace_subset_ids = c("P22102","O00468","Q9UBF2",
#                                        "O15067","O75822","Q9Y266",
#                                        "P42167","P61978","O95801",
#                                        "P55060","O14578","O75717","P35221"),
#                 trace_subset_type = "protein_id")

#pepTracesList_filtered_sum <- subset(pepTracesList_filtered_sum,
#                                    trace_subset_ids = pepTracesList_filtered[[1]]$trace_annotation[id != Sequence]$protein_id[1000:1020],
#                                    trace_subset_type = "protein_id")

#' ## Annotate the genomic coordinates
#' Annotate summed traces
traces_coordinates <- annotateGenomicCoordinates(pepTracesList_filtered_sum, db=ensdb, proteinIdCol="LeadingEnsemblProtein", verbose=T)
#' Transfer annotation to traces list
traces_coordinates_list <- copy(pepTracesList_filtered)
traces_coordinates_list$control$genomic_coord <- traces_coordinates$genomic_coord[traces_coordinates_list$control$trace_annotation$id]
traces_coordinates_list$depleted$genomic_coord <- traces_coordinates$genomic_coord[traces_coordinates_list$depleted$trace_annotation$id]

#' ## Save objects
saveRDS(traces_coordinates_list, "pepTracesList_filtered_coordinates.rda")
saveRDS(traces_coordinates,"pepTracesSum_filtered_coordinates.rds")
