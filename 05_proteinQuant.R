# protein quentification
library(devtools)
#setwd("~/mysonas/CCprofiler")
#load_all()
setwd("~/mysonas/PRPF8/analysis/output")
install_github("CCprofiler/CCprofiler", ref = "DA_module")
library(CCprofiler)

pepTraces <- readRDS("pepTracesImpSPCmultipep.rda")
design_matrix <- readRDS("design_matrix.rda")

protTraces <- proteinQuantification(pepTraces,
                       topN = 2,
                       keep_less = FALSE,
                       rm_decoys = TRUE,
                       use_sibPepCorr = FALSE,
                       use_repPepCorr = FALSE,
                       full_intersect_only = FALSE)

protTraces <- updateTraces(protTraces)

protTracesSum <- integrateTraceIntensities(protTraces,
                                       design_matrix = NULL,
                                       integrate_within = NULL,
                                       aggr_fun = "sum")

saveRDS(protTraces,"protTraces.rds")
saveRDS(protTracesSum,"protTracesSum.rds")

# evaluate protein mass distribution
protTraces_assembly <- annotateMassDistribution(protTraces)
saveRDS(protTraces_assembly,"protTraces_assembly.rds")


diffAssemblyState <- getMassAssemblyChange(protTraces_assembly, design_matrix)
saveRDS(diffAssemblyState,"diffAssemblyState.rds")

fold_change_cutoff=1

plotMassAssemblyChange(diffAssemblyState, FC_cutoff=fold_change_cutoff, PDF=T)

assembly_info <- copy(diffAssemblyState)

all_proteins <- subset(assembly_info, !is.na(change))$protein_id
no_change_proteins <- subset(assembly_info, abs(change) <= fold_change_cutoff)$protein_id
more_assembled_proteins <- subset(assembly_info, change < -fold_change_cutoff)$protein_id
less_assembled_proteins <- subset(assembly_info, change > fold_change_cutoff)$protein_id

write.table(all_proteins,paste0("all_proteins.txt"),sep="\t",quote=F,row.names = F,col.names=F)
write.table(no_change_proteins,paste0("no_change_proteins.txt"),sep="\t",quote=F,row.names = F,col.names=F)
write.table(more_assembled_proteins,paste0("more_assembled_proteins.txt"),sep="\t",quote=F,row.names = F,col.names=F)
write.table(less_assembled_proteins,paste0("less_assembled_proteins.txt"),sep="\t",quote=F,row.names = F,col.names=F)
