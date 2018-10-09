#' ---	
#' title: 'Proteogenomics visualization with CCprofiler'	
#' author: "Max Frank"	
#' date: '2017-11-21'	
#' output:	
#'   html_document:	
#'     toc: true	
#'   html_notebook:	
#'     toc: true	
#'   pdf_document:	
#'     toc: true	
#' ---	
#'
#' # Introduction
#' Proteogenomics is a powerful approach for the detection of novel genes/isoforms by combining evidence
#' on the protein level (usually from Mass spectrometry) with knowledge from genomics (usually from
#' RNA sequencing or genome sequencing). While it is possible to detect the presence or absence of particular
#' genomic variants (e.g. protein isoforms, mutations, InDels) on the genomic level, it is usually difficult
#' to predict the phenotypic effect of these variants from large scale sequencing experiments. Protein
#' correllation profiling allows the detection of protein interaction networks on a large scale which are
#' often immediately relevant for phenotypic effects. Here we try to integrate the power of both approaches,
#' by visulizing individual peptides in the context of known genomic features.
#'
#' # Step by step guide to visualize the plotted peptides in their genomic context
#'
#' ## Data import and Setup

library(data.table)
library(devtools)
# install_github("CCprofiler/CCprofiler", ref="protgen_vis")
library(CCprofiler)
## Read data
project_root <- "C:/Users/Yujia Cai/Desktop/analysis-20180523T183534Z-001/analysis/SummerProject"
map_table <- readRDS("HeLa_protgen_combined_1FPKM_header_mapping_v2.rda")
## map_table <- readRDS(paste0(project_root, "data/proteogenomics/HeLa_protgen_combined_1FPKM_header_mapping_v2.rda"))
#' traces <- readRDS(paste0(project_root, "/../output/pepTracesRaw.rda"))
#' 
#' #' ## Annotate the Leading isoform for each peptide
#' traces <- annotateLeadingIsoform(traces,
#'                                  isoform_col = "isoform_id",
#'                                  output_col = "LeadingIsoform")
#' 
#' #' ## Find the relative position of each peptide in the protein sequence
#' #' With the leadingIsoform we can use the mapping table that contains the sequences of the
#' #' fasta file used for the proteomics search to determine the relative position of each peptide
#' #' within the protein.
#' traces <- annotateRelativePepPos(traces, map_table, multimatch = "first", verbose = T)
#' 
#' #' We do see quite some instan0ces where multiple sequences are found. This is due
#' #' to the fact, that sometimes a single isoform id in the mapping table corresponds
#' #' to multiple protein ids. However in that case they will have the same sequence, so the relative
#' #' position will be the same in every case.
#' 
#' #' Map the genomic position of every peptide
#' 
#' #' For finding the genomic position of all peptides we need to know the
#' #' ENSEMBL protein id of the isoform they originate from. Therefore
#' #' we map the leading isoform to its corresponding protein id.
#' 
#' map_table$IsoformId <- gsub("\\|.*", "", gsub(">.*?\\|","", map_table$header))
#' isomap <- unique(map_table[,.(LeadingIsoform=IsoformId, LeadingEnsemblProteinSp=protein)])
#' isomap <- isomap[!duplicated(LeadingIsoform)] # Take the first protein id if the seq is the same
#' traces$minus$trace_annotation <- merge(traces$minus$trace_annotation, isomap)
#' traces$minus$trace_annotation[, LeadingEnsemblProtein := gsub("_.*", "", LeadingEnsemblProteinSp)]
#' traces$minus$trace_annotation <- traces$minus$trace_annotation[order(id)]
#' 
#' traces$plus$trace_annotation <- merge(traces$plus$trace_annotation, isomap)
#' traces$plus$trace_annotation[, LeadingEnsemblProtein := gsub("_.*", "", LeadingEnsemblProteinSp)]
#' traces$plus$trace_annotation <- traces$plus$trace_annotation[order(id)]

#' Now the genomic position is retrieved with the ensembldb package.
library(EnsDb.Hsapiens.v86)

ensdb <- EnsDb.Hsapiens.v86
# Fetched_data <- annotateGenomicCoordinates(traces$minus, db=ensdb, proteinIdCol="LeadingEnsemblProtein", verbose=T)

#' Note this function is currently very inefficient and takes several hours to run
#' It is recommended to run this on a full traces object over night

#' Here we use a pre-annotated traces object
traces_minus <- readRDS("tracesMinusGenomeAnn.rda")
head(traces_minus$genomic_coord)

traces_plus <- readRDS("tracesPlusGenomeAnn.rda")
head(traces_plus$genomic_coord)

#' ## Plotting peptides in their genomic context
#' Here we use the genomic annotations that are present in the traces object to show the relative
#' position of the peptide sequences

library(ggbio)
plotPeptidesInGenome(traces_plus, gene_id="ENSG00000166855", ensdb=ensdb)
plotPeptidesInGenome2(traces_minus, gene_id="ENSG00000166855", highlight = "LLEGTIVNVPEK", ensdb=ensdb, name="awesomeplot")

#'Also objects that do not have the annotation can be plotted. In that case the position
#' is fetched as needed

plotPeptidesInGenome(traces$plus, gene_id="ENSG00000166855", ensdb=ensdb)

#' If we don't supply an EnsDb object the genomic context is not displayed
plotPeptidesInGenome(traces_plus, gene_id="ENSG00000166855", ensdb=NULL)

#' We can also highlight peptides which is especially powerful to visualize splice isoforms
plotPeptidesInGenome(traces_plus, gene_id="ENSG00000166855", ensdb=ensdb, highlight="LLQDANYNVEK")

#' ## A few examples

geneId <- unique(traces_plus$trace_annotation[protein_id == "P42167", gene_id])
hl <- c("NRPPLPAGTNSK", "SELVANNVTLPAGEQR", "PEFLEDPSVLT", "LKSELVANNVTLPAGEQR", "GPPDFSSDEEREPTPVLGSGAAAAGR", "DVYVQLYLQHLTAR")

geneId <- unique(traces_plus$trace_annotation[protein_id == "Q6P2E9", gene_id])
hl <- c("SSHSTWPVDVSQIK", "QPEGTPLNHFR", "FLITGADQNR", "ALQDVQIR", "VISVSTSER")
hl <- c("EAFQSVVLPAFEK", "EPVLAQLR", "GEVSVALK", "HSQEELLQR", "LGTQEYLQQLESHMK", "SLEPMAGQLSNSVATK")

geneId <- unique(traces_plus$trace_annotation[protein_id == "Q15459", gene_id])
hl <- c("AQEPSAAIPK", "IRQNEINNPK", "PAGPVQAVPPPPPVPTEPK", "QNEINNPK")
hl <- c("GLVPEDDTK", "GLVPEDDTKEK", "LQYEGIFIK", "NKGPVSIK", "QKLQYEGIFIK", "TEDSLMPEEEFLRR", "VQVPNMQDK")

plotPeptidesInGenome(traces_plus, gene_id=geneId, ensdb=ensdb, highlight=hl)


source("C:/Users/Yujia Cai/Desktop/analysis-20180523T183534Z-001/analysis/SummerProject/proteogenomics_plotting_custom.R")
environment(plotPeptidesInGenome2) <- environment(plotPeptidesInGenome) 
environment(plot_tracesCustom) <- environment(CCprofiler:::plot.traces)
pdf("ProteogenomicsVisualization_minus0.85.pdf")
prot_ids <- names(clustMatricesSingle)
filtedprot_ids <- prot_ids[idxs]
maxheights <- sapply(clustMatricesPlSingle, function(x) max(x$height))

for(prot_id in filtedprot_ids[order(maxheights, decreasing = F)]){
  peps <- clustMatricesSingle[[prot_id]]
  print(prot_id)
  geneId <- unique(traces_minus$trace_annotation[protein_id == prot_id, gene_id])
  for(group in unique(peps)){
    # if(nrow(proteinFeaturesfiltMerged[protein_id == prot_id]) > 0) {
      p <- plotPeptidesInGenome2(traces_minus, gene_id=geneId, ensdb=ensdb, 
                                 highlight = names(peps[peps == group]),
                                 legend = F, name = prot_id)
    # }else{
    #   p <- plot(subset(traces_filt, prot_id, "protein_id"), highlight = names(peps[peps == group]))
    # }
    plot(p)
  }
}


dev.off()

pdf("ProteogenomicsVisualization_plus0.85_singletest1.pdf")
prot_ids <- names(clustMatricesSingle)
filtedprot_ids <- prot_ids[idxs]
orderedprot_ids <- names(clustMatricesSingle)[order(corrtable$mincorr, decreasing = T)]
subcorrtable <- corrtable[,prot_ids[idxs] %in% corrtable$protein_id]
Newprot_ids <- intersect(orderedprot_ids, filtedprot_ids)

for(prot_id in Newprot_ids){
  peps <- clustMatricesSingle[[prot_id]]
  geneId <- unique(traces_plus$trace_annotation[protein_id == prot_id, gene_id])
  for(group in unique(peps)){
    # if(nrow(proteinFeaturesfiltMerged[protein_id == prot_id]) > 0) {
    p <- plotPeptidesInGenome(traces_plus, geneId=geneId, ensdb=ensdb, 
                                 highlight = names(peps[peps == group]),
                                 legend = F, name = prot_id)
    # }else{
      # p <- plot(subset(traces_filt, prot_id, "protein_id"), highlight = names(peps[peps == group])
    # }
    print(p)
  }
}


dev.off()
