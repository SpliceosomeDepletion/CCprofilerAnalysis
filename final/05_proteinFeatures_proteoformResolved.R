#' ---
#' title: "Proteoform feature finding"
#' author: "Isabell Bludau"
#' date: "October 21th, 2018"
#' output:
#'   html_document:
#'     toc: true
#'   html_notebook:
#'     toc: true
#'   pdf_document:
#'     toc: true
#' ---

#' # Overview
#' In this script we ...
#'
#' # Step-by-step workflow
#'
#' ## Load CCprofiler package and set working directory:
library(devtools)
#install_github("CCprofiler/CCprofiler", ref = "DA_module")
#library(CCprofiler)
if (length(grep("nas21.ethz.ch",getwd()))>0) {
  setwd("~/mysonas/CCprofiler")
  load_all()
  setwd("~/mysonas/PRPF8/analysis/output")
  knitr::opts_knit$set(root.dir = '~/mysonas/PRPF8/analysis/output')
} else {
  setwd("/Volumes/ibludau-1/CCprofiler")
  load_all()
  setwd("/Volumes/ibludau-1/PRPF8/analysis/output")
  knitr::opts_knit$set(root.dir = '/Volumes/ibludau-1/PRPF8/analysis/output')
}

#' ## Load data
traces_list_pepClusters <- readRDS("traces_list_pepClusters.rds")
traces_sum_pepClusters <- readRDS("traces_sum_pepClusters.rds")
proteinFeatures <- readRDS("proteinFeatures.rda")
design_matrix <- readRDS("design_matrix.rda")

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

saveRDS(reassignedTraces,"traces_sum_reassignedProteoforms.rds")
saveRDS(reassignedTraces_list,"traces_list_reassignedProteoforms.rds")
saveRDS(reassignedFeatures,"reassignedProteoformFeatures.rds")

##### 

length(unique(reassignedTraces$trace_annotation$protein_id))
length(unique(reassignedTraces$trace_annotation[n_proteoforms>1]$protein_id))

proteoform_DT <- subset(reassignedTraces$trace_annotation, select = c("protein_id","n_proteoforms"))
proteoform_DT <- unique(proteoform_DT)

proteoform_count <- as.data.table(table(proteoform_DT$n_proteoforms))
names(proteoform_count) <- c("n_proteoforms","genes")
proteoform_count[,proteoform:=ifelse(n_proteoforms==1,"single functional proteoform", "multiple functional proteoforms")]
proteoform_count$proteoform <- factor(proteoform_count$proteoform, levels = c("single functional proteoform", "multiple functional proteoforms"))

pdf("proteoform_count.pdf", width=6, height=4)
p <- ggplot(proteoform_count , aes(x=n_proteoforms, y=genes, group=proteoform)) +
  facet_wrap(.~proteoform, scales = "free") + 
  geom_bar(stat="identity") +
  geom_text(aes(label=round(genes, digits = 2)), vjust=-0.1, color="black",
            position = position_dodge(0.9), size=3.5) +
  theme_classic() 
print(p)
dev.off()

calibrationFunctions <- readRDS("calibration.rds")
design_matrix <- readRDS("design_matrix.rda")

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


proteoform_genes <- unique(reassignedTraces$trace_annotation[n_proteoforms>1]$protein_id)

pdf("proteoform_traces.pdf")
for (id in proteoform_genes[1:50]) {
  sub <- subset(reassignedTraces_list, trace_subset_ids = id, trace_subset_type = "protein_id")
  plotPeptideCluster(reassignedTraces,id)
  plot(traces = sub,
       design_matrix=design_matrix,
       colour_by="proteoform_id",
       legend = T)
  plotFeatures(feature_table = reassignedFeatures,
               traces = reassignedTraces,
               calibration=calibrationFunctions,
               feature_id = id,
               annotation_label="proteoform_id",
               colour_by="proteoform_id",
               peak_area = T,
               legend = T,
               onlyBest = F)
  plotFeatures(feature_table = reassignedFeatures,
               traces = reassignedTraces_list,
               calibration=calibrationFunctions,
               feature_id = id,
               design_matrix=design_matrix,
               annotation_label="proteoform_id",
               colour_by="proteoform_id",
               peak_area = T,
               legend = T,
               onlyBest = F)
}
dev.off()

proteoform_genes <- unique(reassignedTraces$trace_annotation[n_proteoforms>1]$protein_id)

pdf("proteoform_isoform_traces.pdf")
for (id in proteoform_genes[1:50]) {
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


  