#' ## Load data
traces_list_pepClusters <- readRDS("traces_list_pepClusters.rds")
traces_sum_pepClusters <- readRDS("traces_sum_pepClusters.rds")
proteinFeatures <- readRDS("proteinFeatures.rda")
design_matrix <- readRDS("design_matrix.rda")
calibrationFunctions <- readRDS("calibration.rds")

#' ## Resolve protein features by clusters
proteoformFeaturesResolved <- resolveProteoformSpecificFeatures(
    features=proteinFeatures,
    traces=traces_sum_pepClusters,
    minProteoformIntensityRatio=0.1,
    perturb_cutoff="5%")

saveRDS(proteoformFeaturesResolved,"proteoformFeaturesResolved.rds")

#' ## Score protein features
filteredDataProteoformResolved <- scoreFeatures(
  proteoformFeaturesResolved,
  FDR=0.1, PDF=T,
  name=paste0("qvalueStats_proteoformFeaturesResolved"))

saveRDS(filteredDataProteoformResolved,"filteredDataProteoformResolved.rds")

#' ## Reassign proteoform annotation in tarces and features
library("splitstackshape")
library("Vennerable")
library('mgsub')
reassigned <- refineProteoformsByDetectedFeatures(traces = traces_sum_pepClusters, 
                                                  features = filteredDataProteoformResolved)

reassignedTraces <- reassigned[[1]]
reassignedFeatures <- reassigned[[2]]

reassignedTraces_list <- annotateProteoformsAcrossConditions(traces_list_pepClusters,
                                                              reassignedTraces)

summarizeFeatures(reassignedFeatures, PDF=T, name="proteoformReassiged_feature_summary")

saveRDS(reassignedTraces,"traces_sum_reassignedProteoforms.rds")
saveRDS(reassignedTraces_list,"traces_list_reassignedProteoforms.rds")
saveRDS(reassignedFeatures,"reassignedProteoformFeatures.rds")

#reassignedTraces <- readRDS("traces_sum_reassignedProteoforms.rds")
#reassignedTraces_list <- readRDS("traces_list_reassignedProteoforms.rds")
#reassignedFeatures <- readRDS("reassignedProteoformFeatures.rds")

genes_noProteoforms <- unique(reassignedFeatures[grep("_",proteoform_ids, invert = T)]$protein_id)
genes_withProteoforms <- unique(reassignedFeatures[grep("_",proteoform_ids)]$protein_id)
proteofroms <- unique(reassignedFeatures[grep("_",proteoform_ids)]$proteoform_ids)
proteoform_clusters <- unique(unlist(strsplit(proteofroms, split = ";")))

proteoform_DT <- subset(reassignedFeatures, select = c("protein_id","n_unique_proteoforms"))
proteoform_DT <- unique(proteoform_DT)

proteoform_count <- as.data.table(table(proteoform_DT$n_unique_proteoforms))
names(proteoform_count) <- c("n_unique_proteoforms","proteins")
proteoform_count[,proteoform:=ifelse(n_unique_proteoforms==1,"single proteoform", "multiple proteoforms")]
proteoform_count$proteoform <- factor(proteoform_count$proteoform, levels = c("single proteoform", "multiple proteoforms"))

pdf("proteoform_count.pdf", width=6, height=4)
p <- ggplot(proteoform_count , aes(x=n_unique_proteoforms, y=proteins, group=proteoform)) +
  facet_wrap(.~proteoform, scales = "free") + 
  geom_bar(stat="identity") +
  geom_text(aes(label=round(proteins, digits = 2)), vjust=-0.1, color="black",
            position = position_dodge(0.9), size=3.5) +
  theme_classic() +
  xlab(label = "N proteoforms")
print(p)
dev.off()


#########

#id = "P42167" # LAP2B
#id = "Q15149" # PLEC
#id = "P14618" # PKM
#id = "O94925" # GLS >> looks promising but was not found
#id = "P20020" # ATP2B1 >> works
#id = "Q96HC4" # PDLIM5 >> works
#id = "O43390" # hnrnpr >> works??
#id = "Q9UJU6" # dbnl >> works
#id = "P67936" # tpm4 >> looks promising but was not found
#id = "Q9Y6E0" # >> brilliant example >> 2 terminal peptides are corresponding to isoform 2

proteoform_genes <- c("P42167","P14618","O43390")
pdf("proteoform_LAP2B_PKM_hnrnpr.pdf", height=5, width=8)
for (id in proteoform_genes) {
  sub <- subset(reassignedTraces_list, trace_subset_ids = id, trace_subset_type = "protein_id")
  plotPeptideCluster(reassignedTraces,id)
  plotFeatures(feature_table = reassignedFeatures,
               traces = reassignedTraces_list,
               calibration=calibrationFunctions,
               feature_id = id,
               design_matrix=design_matrix,
               annotation_label="proteoform_id",
               colour_by="proteoform_id",
               peak_area = T,
               legend = T,
               onlyBest = F,
               monomer_MW=T)
  plot(traces = sub,
       design_matrix=design_matrix,
       colour_by="proteoform_id",
       legend = T)
  plot(traces = sub,
       design_matrix=design_matrix,
       colour_by="isoform_id",
       legend = T)
  plot(traces = sub,
       design_matrix=design_matrix,
       colour_by="ensembl_protein_id",
       legend = T)
}
dev.off()


proteoform_genes <- unique(reassignedTraces$trace_annotation[n_proteoforms>1]$protein_id)[1:100]
#proteoform_genes <- c("P42167","P14618","O43390")

pdf("proteoform_traces.pdf", height=5, width=8)
for (id in proteoform_genes) {
  sub <- subset(reassignedTraces_list, trace_subset_ids = id, trace_subset_type = "protein_id")
  plotPeptideCluster(reassignedTraces,id)
  plotFeatures(feature_table = reassignedFeatures,
               traces = reassignedTraces_list,
               calibration=calibrationFunctions,
               feature_id = id,
               design_matrix=design_matrix,
               annotation_label="proteoform_id",
               colour_by="proteoform_id",
               peak_area = T,
               legend = T,
               onlyBest = F,
               monomer_MW=T)
}
dev.off()

pdf("proteoform_isoform_traces.pdf")
for (id in proteoform_genes) {
  sub <- subset(reassignedTraces_list, trace_subset_ids = id, trace_subset_type = "protein_id")
  plotPeptideCluster(reassignedTraces,id)
  plot(traces = sub,
       design_matrix=design_matrix,
       colour_by="proteoform_id",
       legend = T)
  plot(traces = sub,
       design_matrix=design_matrix,
       colour_by="isoform_id",
       legend = T)
  plot(traces = sub,
       design_matrix=design_matrix,
       colour_by="ensembl_protein_id",
       legend = T)
}
dev.off()


get_id_idx <- function(y){
  y = unlist(strsplit(y,"/"))
  y = y[-1]
  y = gsub("([[:alnum:]]*[-])(.*)","\\2",y)
  paste(sort(y), collapse = ";")}

pdf("proteoform_isoform_traces_201_250_isoIdx.pdf")
for (id in proteoform_genes[201:250]) {
  sub <- subset(reassignedTraces_list, trace_subset_ids = id, trace_subset_type = "protein_id")
  sub <- lapply(sub, function(x) {
    x$trace_annotation <- x$trace_annotation[,isoform_id_idx := lapply(isoform_id,get_id_idx), by = "id"]
    return(x)
  })
  class(sub) <- "tracesList"
  plotPeptideCluster(reassignedTraces,id)
  plot(traces = sub,
       design_matrix=design_matrix,
       colour_by="proteoform_id",
       legend = T)
  plot(traces = sub,
       design_matrix=design_matrix,
       colour_by="isoform_id",
       legend = T)
  plot(traces = sub,
       design_matrix=design_matrix,
       colour_by="isoform_id_idx",
       legend = T)
}
dev.off()


  