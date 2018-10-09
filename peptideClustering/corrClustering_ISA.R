library(devtools)
library(matrixStats)
library(fitdistrplus)
library(NbClust)
library('clusteval')
#setwd("~/mysonas/CCprofiler")
#load_all()
setwd("~/mysonas/PRPF8/analysis/output")
install_github("CCprofiler/CCprofiler", ref = "DA_module")
library(CCprofiler)

traces_list <- readRDS("pepTracesImpSPCmultipep.rda")
design_matrix <- readRDS("design_matrix.rda")

traces <- traces_list$minus
#traces <- traces_list$plus
traces_n3 <- filterConsecutiveIdStretches(traces,
                                            min_stretch_length = 3)

traces_n3_maxCorr <- filterByMaxCorr(traces_n3, cutoff = 0.85, plot = T, PDF=T, name="maxCorrHist")

traces_n3_maxCorr <- calculateSibPepCorr(traces_n3_maxCorr)

# Remove proteins with single peptide
multipep <- unique(subset(traces_n3_maxCorr$trace_annotation,! is.na(SibPepCorr))$protein_id)
traces_filtCorr <- subset(traces_n3_maxCorr,trace_subset_ids=multipep,trace_subset_type="protein_id")

genePepList <- getGenePepList(traces_filtCorr)

MinCorrMatrices <- calculateMinCorr(genePepList, plot = TRUE, PDF=TRUE)

FittedDistr <- fitDistr(NULL, MinCorrMatrices, plot = TRUE, PDF=TRUE)

SplicePval <- calculateSplicePval(NULL, MinCorrMatrices, FittedDistr, plot = TRUE, PDF=TRUE)

SplicePval[adj_pval <= 0.05]
# here we choose 0.05 as the cutoff for pval, providing the idxs/indexs of selected protein
# for further analysis in proteogenomicsVisualization.R

SigGenes <- SplicePval[adj_pval <= 0.05]$protein_id

PeptideClust <- clusterPeptides(genePepList[SigGenes], plot = TRUE, PDF=TRUE)

nClust <- lapply(PeptideClust, function(cluster){
  n <- max(cluster$cluster)
  min_pep_count <- min(table(cluster$cluster))
  data.table(n=n,min_pep_count=min_pep_count)})
nClustDT <- do.call(rbind, nClust)
nClustDT[,protein_id := names(nClust)]

sigClust <- subset(nClustDT,(n!=1) & (min_pep_count!=1))

SigGenes_clust <- intersect(sigClust$protein_id,SigGenes)

for (id in SigGenes_clust[1:25]) {
  test_proteins <- id

  pepTest <- subset(traces_filtCorr, trace_subset_ids = test_proteins, trace_subset_type = "protein_id")
  highlight_peps = subset(PeptideClust[[test_proteins]],cluster==2)$id
  pdf(paste0("PeptideTraces_",test_proteins,".pdf"),width=10,height=5)
  plot(pepTest,log = FALSE, legend = T, PDF = FALSE,
   name = paste0("PeptideTraces_",test_proteins), plot = TRUE, highlight = highlight_peps, highlight_col = NULL)
  dev.off()
}


clustering_minus <- PeptideClust
clustering_plus <- PeptideClust

prot_minus <- names(clustering_minus)
prot_plus <- names(clustering_plus)
prot_both <- intersect(prot_minus,prot_plus)

clust_res <- lapply(prot_both, function(gene){
  plus <- clustering_plus[[gene]]
  minus <- clustering_minus[[gene]]
  names(plus) <- c("id","plus")
  names(minus) <- c("id","minus")
  common <- merge(plus,minus,by="id",all=F)
  jaccard <- cluster_similarity(common$plus, common$minus, similarity = "jaccard", method = "independence")
  list(clusters=common,jaccard=jaccard)
})
names(clust_res) <- prot_both

clust_res_filtered <- clust_res[lapply(clust_res, function(l){l$jaccard}) == 1]

proteins_identical_clusters <- names(clust_res_filtered)


pdf(paste0("PeptideTraces_significantArossConditions.pdf"),width=10,height=5)
for (id in proteins_identical_clusters) {
  test_proteins <- id
  pepTest <- subset(traces_list, trace_subset_ids = test_proteins, trace_subset_type = "protein_id")
  isoforms <- clust_res[[id]]$clusters
  traces_minus <- pepTest$minus$traces
  traces_minus[,condition:="Control"]
  traces_plus <- pepTest$plus$traces
  traces_plus[,condition:="PRPF8 depleted"]
  PKM_traces_all <- rbind(traces_minus,traces_plus)
  PKM_traces_all <- merge(PKM_traces_all,isoforms,by="id")
  PKM_traces_all.m <- melt(PKM_traces_all,id.vars=c("id","condition","plus","minus"))
  PKM_traces_all.m$variable <- as.numeric(PKM_traces_all.m$variable)
  PKM_traces_all.m[,line:=paste0(plus,id)]
  g <- ggplot(PKM_traces_all.m,aes(x=variable, y=value, color=factor(plus), group=line)) +
    geom_line() +
    facet_wrap(~ condition,nrow = 2) +
    theme_classic() + scale_colour_manual(values=c("grey","#0085ff","#FFA500"))
  print(g)
}
dev.off()




sig_both_cond <- intersect(SigGenes_minus,SigGenes_plus)

for (id in sig_both_cond) {
  test_proteins <- id
  pepTest <- subset(traces_list, trace_subset_ids = test_proteins, trace_subset_type = "protein_id")
  highlight_peps = subset(PeptideClust[[test_proteins]],cluster==2)$id
  pdf(paste0("PeptideTraces_",test_proteins,".pdf"),width=10,height=5)
  plot(pepTest,log = FALSE, legend = T, PDF = FALSE,
   name = paste0("PeptideTraces_",test_proteins), plot = TRUE, highlight = highlight_peps, highlight_col = NULL)
  dev.off()
}

## Interesting
#Q15691 >> https://www.uniprot.org/uniprot/Q15691 >> https://www.uniprot.org/citations/25217626
#Q676U5 >> multiple alternative splicing variants annotated
# O60524 >> ubiquitinylation >> downregulated upon depletion
# Q9NQG6

##changes upon splicing depletion
#Q9UKG1 >> terminal amino acids in separate cluster
# Q96CW5 >> only one isoform changes assembly state
# Q13263 >> only one isoform downregulated
# O75717
# Q9Y4E8 >> ubiquin conjigation related
