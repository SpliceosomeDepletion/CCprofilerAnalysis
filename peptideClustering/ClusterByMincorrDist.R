calculateMinCorr <- function(traces, corrMatrices = NULL) {
  if (is.null(corrMatrices)) {
    corrMatrices <- getCorrelationMatrices(traces)
  }
  minCorrMatrices <- lapply(genePepList, function(gene){
    genecorr <- cor(gene)
    mincorr <- min(genecorr)
    minpeps <- rownames(which(genecorr == min(genecorr), arr.ind=T))
    ints <- apply(gene, 2, max)
    minint <- min(ints)
    return(data.table(mincorr = mincorr, minint = minint,
                      minpep1 = minpeps[1], minpep2 = minpeps[2]))
  })
  return(minCorrMatrices)
}


fitDistr <- function(prot_names = NULL, mincorrtable, distr = "gamma", split_value = 0.2, plot = FALSE) {
  if (is.null(prot_names)) {
    mincorrs <- mincorrtable$mincorr
  } else {
    prot_idxs <- which(prot_names %in% mincorrtable$protein_id)
    mincorrs <- mincorrtable$mincorr[prot_idxs]
  }
  GammaFitted <- fitdist(1 - mincorrs[mincorrs > split_value], distr)
  if (plot == TRUE) {
    plot(GammaFitted)
  }
  return(GammaFitted)
}


calculateSplicePval <- function(prot_names= NULL, mincorrtable, distr_fitted, 
                                adj.method = "fdr", plot = FALSE) {
  if (is.null(prot_names)) {
    mincorrs <- mincorrtable$mincorr
  } else {
    prot_idxs <- which(prot_names %in% mincorrtable$protein_id)
    mincorrs <- mincorrtable$mincorr[prot_idxs]
  }
  
# here we use pgamma directly in the function for our gamma distr example, if we want to
# use other distr, we need a general way to represent all ditribution, e.g.,
# "p + distr abbr."
  pval <- 1-(pgamma(1 - mincorrs, shape = distr_fitted$estimate[1][[1]],
                    rate = distr_fitted$estimate[2][[1]]))
  pval_adj <- p.adjust(pval, adj.method)
  if (plot == TRUE) {
    hist(pval, breaks = 50)
    hist(pval_adj, breaks = 100)
  }
  return(pval_adj)
}

# in cluster function, here we return cut the tree by number of clustering, we could also
# do cutoff. We return cl here, can also return groups by adding an option in the function.
clusterPeptides <- function(traces, corrMatrices = NULL, method = "single", clusterN = 2, ...) {
  if (is.null(corrMatrices)) {
    corrMatrices <- getCorrelationMatrices(traces)
  }
  clustMatricesMethod <- lapply(corrMatrices, function(gene){
    cl <- hclust(as.dist(1-gene), method)
    groups <- cutree(cl, k = clusterN)
    return(cl)
  })
  return(clustMatricesMethod)
}
