library(devtools)
#library(matrixStats)
#library(fitdistrplus)
#library(NbClust)
#library('clusteval')
#setwd("~/mysonas/CCprofiler")
#load_all()
setwd("~/mysonas/PRPF8/analysis/output")
install_github("CCprofiler/CCprofiler", ref = "proteoformAnnotation")
library(CCprofiler)

traces_list <- readRDS("pepTracesImpSPCmultipep.rda")
design_matrix <- readRDS("design_matrix.rda")

#traces_maxCorr <- filterByMaxCorr(traces_list, cutoff = 0.85, plot = T, PDF=T, name="maxCorrHist")

traces_maxCorrItr <- iterativeMaxCorrFilter(traces_list,
  cutoff = 0.85, plot = T, PDF=T, name="maxCorrHist")

traces_maxCorr_multi <-  filterSinglePeptideHits(traces_maxCorrItr)

#test_prot <- unique(traces_maxCorr_multi$minus$trace_annotation$protein_id)[1:500]
#test <- subset(traces_maxCorr_multi, trace_subset_ids=test_prot,trace_subset_type="protein_id")

traces_maxCorr_multi_minCorr <- calculateMinCorr(traces_maxCorr_multi,
  plot = TRUE, PDF=TRUE)

traces_maxCorr_multi_minCorr_pval <- estimateProteoformPval(traces_maxCorr_multi_minCorr,
  plot = TRUE, PDF=TRUE)

traces_maxCorr_multi_minCorr_pval_clustered_multiCond <- clustPepMultiCond(traces_maxCorr_multi_minCorr_pval, clusterN = NULL, clusterH = 0.75, plot = TRUE, PDF=TRUE)

# in case clustering should be performed for each condition separately
#traces_maxCorr_multi_minCorr_pval_clustered <- clusterPeptides(traces_maxCorr_multi_minCorr_pval,
#  plot = TRUE, PDF=TRUE)

# check examples
PKM <- c("P14618")
PKM_traces <- subset(traces_maxCorr_multi_minCorr_pval_clustered_multiCond,
  trace_subset_ids=PKM,trace_subset_type="protein_id")
plot(PKM_traces, design_matrix = design_matrix, PDF = TRUE,
 name = paste0("PeptideTraces_","PKM"), plot = TRUE)


GART <- c("P22102")
GART_traces <- subset(traces_maxCorr_multi_minCorr_pval_clustered_multiCond,
  trace_subset_ids=GART,trace_subset_type="protein_id")
plot(GART_traces, design_matrix = design_matrix, PDF = TRUE,
 name = paste0("PeptideTraces_","GART"), plot = TRUE)

 # How many of the proteins that were significant in the clustering,
 # actually cluster in multiple proteoforms?
 sig_minus <- traces_maxCorr_multi_minCorr_pval_clustered_multiCond$minus$trace_annotation[proteoform_pval_adj <= 0.05]
 sig_minus[, n_proteoforms := length(unique(proteoform_id)), by=c("protein_id")]
 table(unique(subset(sig_minus,select=c("protein_id","n_proteoforms")))$n_proteoforms)
 sig_minus_prot <- unique(sig_minus$protein_id)

 non_sig_minus <- traces_maxCorr_multi_minCorr_pval_clustered_multiCond$minus$trace_annotation[proteoform_pval_adj > 0.05]
 non_sig_minus[, n_proteoforms := length(unique(proteoform_id)), by=c("protein_id")]
 table(unique(subset(non_sig_minus,select=c("protein_id","n_proteoforms")))$n_proteoforms)


 sig_plus <- traces_maxCorr_multi_minCorr_pval_clustered_multiCond$plus$trace_annotation[proteoform_pval_adj <= 0.05]
 sig_plus[, n_proteoforms := length(unique(proteoform_id)), by=c("protein_id")]
 table(unique(subset(sig_plus,select=c("protein_id","n_proteoforms")))$n_proteoforms)
 sig_plus_prot <- unique(sig_plus$protein_id)

 non_sig_plus <- traces_maxCorr_multi_minCorr_pval_clustered_multiCond$plus$trace_annotation[proteoform_pval_adj > 0.05]
 non_sig_plus[, n_proteoforms := length(unique(proteoform_id)), by=c("protein_id")]
 table(unique(subset(non_sig_plus,select=c("protein_id","n_proteoforms")))$n_proteoforms)

####



traces_maxCorr_multi_minCorr_pval_clustered_multiCond$minus$trace_annotation[, multi_avail := ifelse(is.na(cluster_multiCond),0,1)]
traces_maxCorr_multi_minCorr_pval_clustered_multiCond$minus$trace_annotation[, cluster_multiCond := ifelse(is.na(cluster_multiCond),0,cluster_multiCond)]
traces_maxCorr_multi_minCorr_pval_clustered_multiCond$minus$trace_annotation[,jaccard :=
  cluster_similarity(cluster, cluster_multiCond,similarity = "jaccard", method = "independence"), by=c("protein_id","multi_avail")]
pdf("jaccard_minus.pdf")
  hist(traces_maxCorr_multi_minCorr_pval_clustered_multiCond$minus$trace_annotation$jaccard, breaks = 50)
dev.off()

traces_maxCorr_multi_minCorr_pval_clustered_multiCond$plus$trace_annotation[, multi_avail := ifelse(is.na(cluster_multiCond),0,1)]
traces_maxCorr_multi_minCorr_pval_clustered_multiCond$plus$trace_annotation[, cluster_multiCond := ifelse(is.na(cluster_multiCond),0,cluster_multiCond)]
traces_maxCorr_multi_minCorr_pval_clustered_multiCond$plus$trace_annotation[,jaccard :=
  cluster_similarity(cluster, cluster_multiCond,similarity = "jaccard", method = "independence"), by=c("protein_id","multi_avail")]
pdf("jaccard_plus.pdf")
  hist(traces_maxCorr_multi_minCorr_pval_clustered_multiCond$plus$trace_annotation$jaccard, breaks = 50)
dev.off()



getPepDistMatrix <- function(traces) {
  UseMethod("getPepDistMatrix", traces)
}

#' @describeIn getDistanceMatrix calculate pairwise correlation between peptides
#' @export
getPepDistMatrix.traces <- function(traces) {
  genePeptideList <- getGenePepList(traces)
  idx <- seq_along(genePeptideList)
  distMatrix <- lapply(idx, function(i){
    gene <- genePeptideList[[i]]
    gene_name <- names(genePeptideList)[[i]]
    genecorr <- cor(gene)
    return(genecorr)
  })
  names(distMatrix) <- names(genePeptideList)
  return(distMatrix)
}

#' @describeIn getDistanceMatrix calculate pairwise correlation between peptides
#' @export
getPepDistMatrix.tracesList <- function(tracesList) {
  #.tracesListTest(tracesList)
  res <- lapply(tracesList, getDistanceMatrix.traces)
  return(res)
}


getPepDistMultiCond <- function(tracesList) {
  #.tracesListTest(tracesList)
  dist_list <- getPepDistMatrix(tracesList)
  proteins.l <- lapply(dist_list, names)
  proteins <- intersect(proteins.l[[1]],proteins.l[[2]])
  idx <- seq_along(proteins)
  res <- lapply(idx, function(i){
    #print(i)
    gene <- proteins[[i]]
    x <- dist_list[[1]][[gene]]
    y <- dist_list[[2]][[gene]]
    x.dt <- as.data.table(x,keep.rownames = "id")
    x.m <- melt(x.dt,id.vars="id")
    y.dt <- as.data.table(y,keep.rownames = "id")
    y.m <- melt(y.dt,id.vars="id")
    z <- merge(x.m,y.m,by=c("id","variable"),all=F)
    if(nrow(z) == 0) {
      return(NA)
    }
    z.u <- z[,.(meanCorr = rowMeans(.SD)), by = c("id","variable")]
    a <- dcast(z.u,id ~ variable,value.var="meanCorr")
    a.m <- as.matrix(a, rownames="id")
    a.d <- as.dist(a.m)
    return(a.d)
  })
  names(res) <- proteins
  # remove NAs
  no_NA <- which(!is.na(res))
  res <- res[no_NA]
  singlePep <- which(sapply(res, nrow) > 1)
  res <- res[singlePep]
  return(res)
}

clusterPeptidesMultiCond <- function(tracesList, method = c("complete","single"), clusterH = 0.5,
                    clusterN = NULL, index = "silhouette",
                    plot = FALSE, PDF=FALSE, name="hclustMultiCond", ...) {
  method <- match.arg(method)
  if (PDF) {
    pdf(paste0(name,".pdf"))
  }
  combi_dist <- getPepDistMultiCond(tracesList)
  proteins <- names(combi_dist)
  idx <- seq_along(proteins)
  clustMatricesMethod <- lapply(idx, function(i){
    genecorr <- combi_dist[[i]]
    gene_name <- names(combi_dist)[[i]]
    if (is.null(clusterH)) {
      if (is.null(clusterN)) {
        nbClust_res <- NbClust(data=genecorr, diss=as.dist(1-genecorr),
        distance = NULL, min.nc=2,
        max.nc=min(4,nrow(genecorr)-1)[1],method = method, index = index)
        clusterN <- nbClust_res$Best.nc[1]
      }
    }
    cl <- hclust(as.dist(1-genecorr), method)
    if (plot == TRUE) {
      #plot(cl, labels = FALSE)
      plot(as.dendrogram(cl), ylim = c(0,1))
    }
    groups <- cutree(cl, k = clusterN, h = clusterH)
    groupsDT <- data.table(id=names(groups),cluster_multiCond=groups,protein_id=gene_name,
                          proteoform_id_multiCond = paste0(gene_name,"_",groups))
    return(groupsDT)
  })
  if (PDF) {
    dev.off()
  }
  #return(clustMatricesMethod)
  clustMatrices <- do.call(rbind, clustMatricesMethod)
  tracesList[[1]]$trace_annotation <- merge(tracesList[[1]]$trace_annotation,
    clustMatrices,by=c("id","protein_id"),sort=F, all.x=T, all.y=F)
  tracesList[[2]]$trace_annotation <- merge(tracesList[[2]]$trace_annotation,
    clustMatrices,by=c("id","protein_id"),sort=F, all.x=T, all.y=F)
  return(tracesList)
}



combineTracesMutiCond <- function(tracesList){
  .tracesListTest(tracesList)
  cond <- names(tracesList)
  idx <- seq_along(cond)
  traces <- lapply(idx, function(i){
    t <- tracesList[[i]]
    trac <- t$traces
    trac.m <- melt(trac,id.vars="id")
    trac.m[, cond := cond[[i]]]
    return(trac.m)
    })
  traces_combi <- rbindlist(traces)
  traces_all <- dcast(traces_combi,id ~ variable+cond,value.var="value",fill=0)
  names(traces_all) <- c("id",seq(1,ncol(traces_all)-1,1))
  setcolorder(traces_all, c(seq(1,ncol(traces_all)-1,1),"id"))

  trace_type_all <- tracesList[[1]]$trace_type

  trace_annotation <- lapply(tracesList, function(t){
    res <- subset(t$trace_annotation,select=c("id","protein_id"))
    return(res)
    })
  trace_annotation_combi <- rbindlist(trace_annotation)
  trace_annotation_combi <- unique(trace_annotation_combi)
  trace_annotation_combi[order(match(id, traces_all$id))]

  fraction_annotation_all <- data.table(id=seq(1,ncol(traces_all)-1,1))

  combi_tracesList <- list(
    traces = traces_all,
    trace_type = trace_type_all,
    trace_annotation = trace_annotation_combi,
    fraction_annotation = fraction_annotation_all
  )
  class(combi_tracesList) <- "tracesList"
  .tracesListTest(combi_tracesList)
  return(combi_tracesList)
}


getPepCorrMultiCond <- function(tracesList) {
  #.tracesListTest(tracesList)
  genePeptideList <- getGenePepList(traces)
  idx <- seq_along(genePeptideList)
  distMatrix <- lapply(idx, function(i){
    gene <- genePeptideList[[i]]
    gene_name <- names(genePeptideList)[[i]]
    genecorr <- cor(gene)
    return(genecorr)
  })
  names(distMatrix) <- names(genePeptideList)
  return(distMatrix)
}

clustPepMultiCond <- function(tracesList, tracesListRaw, method = c("complete","single"), clusterH = 0.5,
                    clusterN = NULL, index = "silhouette",
                    plot = FALSE, PDF=FALSE, name="hclustMultiCond", ...) {
  method <- match.arg(method)
  if (PDF) {
    pdf(paste0(name,".pdf"))
  }
  combi_dist <- getPepDistMultiCond(tracesList)
  proteins <- names(combi_dist)
  idx <- seq_along(proteins)
  clustMatricesMethod <- lapply(idx, function(i){
    genecorr <- combi_dist[[i]]
    gene_name <- names(combi_dist)[[i]]
    if (is.null(clusterH)) {
      if (is.null(clusterN)) {
        nbClust_res <- NbClust(data=genecorr, diss=as.dist(1-genecorr),
        distance = NULL, min.nc=2,
        max.nc=min(4,nrow(genecorr)-1)[1],method = method, index = index)
        clusterN <- nbClust_res$Best.nc[1]
      }
    }
    cl <- hclust(as.dist(1-genecorr), method)
    if (plot == TRUE) {
      #plot(cl, labels = FALSE)
      plot(as.dendrogram(cl), ylim = c(0,1))
    }
    groups <- cutree(cl, k = clusterN, h = clusterH)
    groupsDT <- data.table(id=names(groups),cluster_multiCond=groups,protein_id=gene_name,
                          proteoform_id_multiCond = paste0(gene_name,"_",groups))
    return(groupsDT)
  })
  if (PDF) {
    dev.off()
  }
  #return(clustMatricesMethod)
  clustMatrices <- do.call(rbind, clustMatricesMethod)
  tracesList[[1]]$trace_annotation <- merge(tracesList[[1]]$trace_annotation,
    clustMatrices,by=c("id","protein_id"),sort=F, all.x=T, all.y=F)
  tracesList[[2]]$trace_annotation <- merge(tracesList[[2]]$trace_annotation,
    clustMatrices,by=c("id","protein_id"),sort=F, all.x=T, all.y=F)
  return(tracesList)
}









####


# check examples:
PKM <- c("P14618")
PKM_traces <- subset(traces_maxCorr_multi_minCorr_pval_clustered,
  trace_subset_ids=PKM,trace_subset_type="protein_id")
plot(PKM_traces, design_matrix = design_matrix, PDF = TRUE,
 name = paste0("PeptideTraces_","PKM"), plot = TRUE)


GART <- c("P22102")
GART_traces <- subset(traces_maxCorr_multi_minCorr_pval_clustered,
  trace_subset_ids=GART,trace_subset_type="protein_id")
plot(GART_traces, design_matrix = design_matrix, PDF = TRUE,
 name = paste0("PeptideTraces_","GART"), plot = TRUE)


##

minus <- subset(traces_maxCorr_multi_minCorr_pval_clustered$minus$trace_annotation,
  select=c("id","protein_id","proteoform_pval_adj", "proteoform_id", "cluster"))

plus <- subset(traces_maxCorr_multi_minCorr_pval_clustered$plus$trace_annotation,
  select=c("id","protein_id","proteoform_pval_adj", "proteoform_id", "cluster"))

combi <- merge(minus,plus,by=c("id","protein_id"),suffixes=c("_minus","_plus"),all=F)

combi[,jaccard := cluster_similarity(cluster_minus, cluster_plus,
  similarity = "jaccard", method = "independence"), by="protein_id"]

jaccard_unique <- unique(subset(combi, select=c("protein_id","jaccard")))

pdf("jaccard.pdf")
  hist(jaccard_unique$jaccard, breaks = 50)
dev.off()

sig_jaccard <- subset(combi,(proteoform_pval_adj_minus <= 0.05) & (proteoform_pval_adj_plus <= 0.05))
sig_jaccard[, n_proteoforms_minus := length(unique(proteoform_id_minus)), by=c("protein_id")]
sig_jaccard[, n_proteoforms_plus := length(unique(proteoform_id_plus)), by=c("protein_id")]

sig_jaccard_unique <- unique(subset(sig_jaccard, select=c("protein_id","jaccard","n_proteoforms_minus","n_proteoforms_plus")))

pdf("sig_jaccard.pdf")
  hist(sig_jaccard_unique$jaccard, breaks = 50)
  hist(sig_jaccard_unique$n_proteoforms_minus-sig_jaccard_unique$n_proteoforms_plus, breaks = 50)
dev.off()

sig_jaccard_unique[jaccard==1]
pdf("jaccard_proteins.pdf")
for (p in sig_jaccard_unique$protein_id){
  test_traces <- subset(traces_maxCorr_multi_minCorr_pval_clustered,
    trace_subset_ids=p,trace_subset_type="protein_id")
  plot(test_traces, design_matrix = design_matrix, PDF = F,
   name = paste0("PeptideTraces_",p), plot = TRUE)
}
dev.off()
# Q9NQG6 interesting case

####
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



####
plot.tracesList <- function(traces,design_matrix=NULL,PDF=T,plot=T,name="plot") {

  tracesList <- lapply(names(traces), function(tr){
    res <- toLongFormat(traces[[tr]]$traces)
    res$Condition <- design_matrix[Sample_name == tr, Condition]
    res$Replicate <- design_matrix[Sample_name == tr, Replicate]
    res
  })
  traces_long <- do.call("rbind", tracesList)

  isoform_annotation <- lapply(names(traces), function(tr){
    t <- subset(traces[[tr]]$trace_annotation,select=c("id","proteoform_id"))
    t[,Condition := tr]})
    isoform_annotation <- unique(do.call("rbind", isoform_annotation))
    #isoform_annotation$proteoform_id <- gsub("ENSG[0-9]+-","",isoform_annotation$proteoform_id)
  traces_long <- merge(traces_long,isoform_annotation, by=c("id","Condition"),all.x=T)

  traces_long[,isoform_id:=paste0(id,"_",proteoform_id)]
  traces_long[,line:=proteoform_id]

  ## Create a common fraction annotation
  traces_frac <- unique(do.call("rbind", lapply(traces, "[[", "fraction_annotation")))
  traces_frac <- unique(subset(traces_frac, select = names(traces_frac) %in% c("id","molecular_weight")))
  traces_long <- merge(traces_long,traces_frac,by.x="fraction",by.y="id")

  p <- ggplot(traces_long) +
    xlab('fraction') +
    ylab('intensity') +
    theme_bw() +
    theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
    theme(plot.margin = unit(c(1,.5,.5,.5),"cm")) +
    ggtitle(name) +
    #scale_color_manual(values=colorMap) +
    theme(plot.title = element_text(vjust=19,size=10))

  p <- p + facet_grid(Condition ~ Replicate) +
    #geom_line(aes_string(x='fraction', y='intensity', color='proteoform_id', group='line'))
    geom_line(aes_string(x='fraction', y='intensity', color='line', group='isoform_id'))


  if(PDF){
    pdf(paste0(name,".pdf"))
  }
  if(plot){
    plot(p)
  }else{
    return(p)
  }
  if(PDF){
    dev.off()
  }
}
