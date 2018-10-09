## Hierarchical clustering based detection of splicing peptides
library(CCprofiler)
library(matrixStats)
library(fitdistrplus)
source("C:/Users/Yujia Cai/Desktop/analysis-20180523T183534Z-001/analysis/SummerProject/filterByMaxCorr.R")
source("C:/Users/Yujia Cai/Desktop/analysis-20180523T183534Z-001/analysis/SummerProject/ClusterByMincorrDist.R")
setwd("C:/Users/Yujia Cai/Desktop/analysis-20180523T183534Z-001/analysis/SummerProject")
traces <- readRDS('../output/pepTracesRaw.rda')
design_matrix <- readRDS('../output/design_matrix.rda')

## Focus on one condition first
traces <- traces$plus

traces_filt <- filterConsecutiveIdStretches(traces, 
                                            min_stretch_length = 3)
traces_filtCorr <- filterByMaxCorr(traces_filt, cutoff = 0.85)

genes <- unique(traces_filtCorr$trace_annotation$protein_id)
intMat <- getIntensityMatrix(traces_filtCorr)
filtCorrMatrices <- getCorrelationMatrices(traces_filtCorr)
genePepList <- lapply(genes, FUN = function(gene){
  peps <- traces_filtCorr$trace_annotation[protein_id == gene, id]
  res <- intMat[peps,]
  return(t(res))
})
names(genePepList) <- genes

MincorrMatrices <- calculateMinCorr(traces_filtCorr, filtCorrMatrices)

mincorrtable = do.call(rbind, MincorrMatrices)
mincorrtable$protein_id <- genes

FittedDistr <- fitDistr(NULL, mincorrtable)

SplicePval <- calculateSplicePval(NULL, mincorrtable, FittedDistr)

# here we choose 0.05 as the cutoff for pval, providing the idxs/indexs of selected protein
# for further analysis in proteogenomicsVisualization.R
idxs <- which(TRUE == (SplicePval <= 0.05))

PeptideClust <- clusterPeptides(traces_filtCorr, filtCorrMatrices)
