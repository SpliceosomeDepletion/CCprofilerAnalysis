library(nFactors)

method = "complete"
index = "silhouette"

traces <- copy(traces_maxCorr_multi_minCorr_pval)

genePeptideList <- getGenePepList(traces)
idx <- seq_along(genePeptideList)


idx <- idx[4:50]
factorNumbers <- lapply(idx, function(i){
  print(i)
  gene <- genePeptideList[[i]]
  gene_name <- names(genePeptideList)[[i]]
  genecorr <- cor(gene)
  ev <- eigen(genecorr) # get eigenvalues
  ap <- parallel(subject=nrow(gene),var=ncol(gene),
    rep=100,cent=.05)
  if(length(ev$values) > 2) {
    nS <- nScree(x=ev$values, aparallel=ap$eigen$qevpea)
    pdf(paste0(gene_name,"nFactors.pdf"))
      plotnScree(nS) 
    dev.off()
    fit <- factanal(gene, nS$Components$nparallel, rotation="varimax")
    load <- fit$loadings[,1:ncol(fit$loadings)]
    #print(fit, digits=2, cutoff=.3, sort=TRUE)
    #nbClust_res <- NbClust(data=genecorr, diss=as.dist(1-genecorr),
    #distance = NULL, min.nc=2,
    #max.nc=min(4,nrow(genecorr)-1)[1],method = method, index = index)
    #clusterN <- nbClust_res$Best.nc[1]
  }
  if(length(ev$values) > 2) {
    #res <- data.table(nS=nS$Components$nparallel,clusterN=clusterN)
    res <- load
  } else {
    #res <- data.table(nS=1,clusterN=1)
    res <- 1
  }
  return(res)
})

res <- do.call(rbind,factorNumbers)


  cl <- hclust(as.dist(1-genecorr), method)
  if (plot == TRUE) {
    #plot(cl, labels = FALSE)
    plot(as.dendrogram(cl), ylim = c(0,1))
  }
  groups <- cutree(cl, k = clusterN, h = clusterH)
  groupsDT <- data.table(id=names(groups),cluster=groups,protein_id=gene_name,
                        proteoform_id = paste0(gene_name,"_",groups))
  return(groupsDT)
