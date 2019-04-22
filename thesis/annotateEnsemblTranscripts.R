traces_exon_pval <- readRDS("traces_exon_pval.rds")

#' Import biomaRt:
library(biomaRt)
ensembl <- useMart("ENSEMBL_MART_ENSEMBL", dataset="hsapiens_gene_ensembl", host = "www.ensembl.org")
ensembl <- useDataset("hsapiens_gene_ensembl",mart=ensembl)

traces_location_pval <- readRDS("traces_location_pval.rds")
subset(traces_location_pval,"P61978", trace_subset_type = "protein_id")[]

x = subset(traces_exon_pval,"P61978", trace_subset_type = "protein_id")$genomic_coord

exon_peptide_mapping <- lapply(x, function(y){y$exon_id})
exon_peptide_mapping <- melt(exon_peptide_mapping)
exon_peptide_mapping <- as.data.table(exon_peptide_mapping)
setnames(exon_peptide_mapping,c("value","L1"),c("ensembl_exon_id","peptide"))

exons <- unique(unlist(lapply(x, function(y){y$exon_id})))

bm <- getBM(attributes = c("ensembl_exon_id", "ensembl_transcript_id"),
            filters = "ensembl_exon_id",
            values = exons,
            mart = ensembl)

exon_peptide_mapping <- merge(bm,exon_peptide_mapping,all.x=T,all.y=T,by="ensembl_exon_id")
exon_peptide_mapping <- as.data.table(exon_peptide_mapping)


peptides_ENST00000360384 <- unique(exon_peptide_mapping[ensembl_transcript_id=="ENST00000360384"]$peptide)
peptides_ENST00000457156 <- unique(exon_peptide_mapping[ensembl_transcript_id=="ENST00000457156"]$peptide)

shared_peptides <- intersect(peptides_ENST00000360384,peptides_ENST00000457156)
unique_peptides_ENST00000360384 <- peptides_ENST00000360384[!peptides_ENST00000360384 %in% shared_peptides]
unique_peptides_ENST00000457156 <- peptides_ENST00000457156[!peptides_ENST00000457156 %in% shared_peptides]

proteoform_DiffExprPep <- readRDS("proteoform_DiffExprPep.rda")

global_P61978_diff <- unique(subset(proteoform_DiffExprPep[id %in% exon_peptide_mapping$peptide], select=c(id,global_log2FC,global_log2FC_imp,global_pBHadj)))

reassignedTraces_list <- readRDS("traces_list_reassignedProteoforms.rds")
reassignedTraces_list_sub <- subset(reassignedTraces_list, trace_subset_ids = "P61978", trace_subset_type = "protein_id")

reassignedTraces_list_sub <- lapply(reassignedTraces_list_sub, function(t){
  t$trace_annotation[,transcript:=ifelse(id %in% shared_peptides, "ENST00000360384/ENST00000457156","ENST00000360384"), by="id"]
  t$trace_annotation[,transcript:=ifelse(id %in% c(shared_peptides,unique_peptides_ENST00000360384), transcript, "other"), by="id"]
  return(t)
  })
class(reassignedTraces_list_sub) <- "tracesList"

reassignedFeatures <- readRDS("reassignedProteoformFeatures.rds")
reassignedFeatures_sub <- subset(reassignedFeatures, protein_id=="P61978")

design_matrix <- readRDS("design_matrix.rda")
calibrationFunctions <- readRDS("calibration.rds")

plotFeatures(feature_table = reassignedFeatures_sub,
             traces = reassignedTraces_list_sub,
             calibration=calibrationFunctions,
             feature_id = "P61978",
             design_matrix=design_matrix,
             annotation_label="transcript",
             colour_by="transcript",
             peak_area = F,
             legend = T,
             onlyBest = F,
             monomer_MW=F)


traces <- traces_exon_pval
protein="P61978"
PDF=FALSE
closeGaps=T

plotPeptideCluster(traces_exon_pval,"O75822",PDF=T, closeGaps=F)
plotPeptideCluster(traces_exon_pval,"P27708",PDF=T, closeGaps=F)
plotPeptideCluster(traces_exon_pval,"Q13263",PDF=T, closeGaps=F)
plotPeptideCluster(traces_exon_pval,"Q15056",PDF=T, closeGaps=F)
plotPeptideCluster(traces_exon_pval,"Q9Y6K5",PDF=T, closeGaps=F)
plotPeptideCluster(traces_exon_pval,"Q9NZI8",PDF=T, closeGaps=F)
plotPeptideCluster(traces_exon_pval,"P42167",PDF=T, closeGaps=F) # TMPO / LAP2Beta
plotPeptideCluster(traces_exon_pval,"P61978",PDF=T, closeGaps=F) # HNRNPK
plotPeptideCluster(traces_exon_pval,"Q9H6S3",PDF=T, closeGaps=F) # EPS8L2

pdf("P14618_sequence_cluster.pdf",width=10,height=10)
plotPeptideCluster(traces_exon_pval,"P14618",PDF=F, closeGaps=T)
dev.off()

reassignedFeatures <- readRDS("reassignedProteoformFeatures.rds")
design_matrix <- readRDS("design_matrix.rda")
calibrationFunctions <- readRDS("calibration.rds")
reassignedTraces_list <- readRDS("traces_list_reassignedProteoforms.rds")

proteoform_genes <- c("P42167","P61978","Q9H6S3","P14618","Q9Y6E0","Q8IWV8")
pdf("proteoform_LAP2AB_HNRNPK_EPS8L2_PKM_etc.pdf", height=5, width=8)
for (id in proteoform_genes) {
  sub <- subset(reassignedTraces_list, trace_subset_ids = id, trace_subset_type = "protein_id")
  plotPeptideCluster(traces_exon_pval,id, closeGaps=F)
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
               monomer_MW=F)
  plotFeatures(feature_table = reassignedFeatures,
               traces = reassignedTraces_list,
               calibration=calibrationFunctions,
               feature_id = id,
               design_matrix=design_matrix,
               annotation_label="isoform_id",
               colour_by="isoform_id",
               peak_area = T,
               legend = T,
               onlyBest = F,
               monomer_MW=F)
}
dev.off()

proteoform_genes <- c("Q9Y6E0")
pdf("proteoform_proteolytic.pdf", height=3, width=8)
for (id in proteoform_genes) {
  sub <- subset(reassignedTraces_list, trace_subset_ids = id, trace_subset_type = "protein_id")
  plotPeptideCluster(traces_exon_pval,id, closeGaps=F)
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
               monomer_MW=F)
  plotFeatures(feature_table = reassignedFeatures,
               traces = reassignedTraces_list,
               calibration=calibrationFunctions,
               feature_id = id,
               design_matrix=design_matrix,
               annotation_label="isoform_id",
               colour_by="isoform_id",
               peak_area = T,
               legend = T,
               onlyBest = F,
               monomer_MW=F)
}
dev.off()

 plotPeptideCluster <- function(traces,protein, PDF=FALSE, closeGaps=FALSE){
  traces <- subset(traces, protein, trace_subset_type="protein_id")
  if ("genomic_coord" %in% names(traces)) {
    traces$trace_annotation[,exon_id:=lapply(id,getGenomicCoord, traces=traces), by="id"]
    traces$trace_annotation[,min_exon_id_start:=min(PeptidePositionStart), by="exon_id"]
    getExonLevels <- unique(subset(traces$trace_annotation, select=c("exon_id","min_exon_id_start")))
    setorder(getExonLevels, min_exon_id_start)
    traces$trace_annotation$exon_id <- factor(traces$trace_annotation$exon_id, levels=getExonLevels$exon_id)
    getPalette = colorRampPalette(brewer.pal(9, "Spectral"))
    
    exons <- unique(unlist(lapply(traces$genomic_coord, function(y){y$exon_id})))
    ensembl <- useMart("ENSEMBL_MART_ENSEMBL", dataset="hsapiens_gene_ensembl", host = "www.ensembl.org")
    ensembl <- useDataset("hsapiens_gene_ensembl",mart=ensembl)
    bm <- getBM(attributes = c("ensembl_exon_id", "ensembl_transcript_id"),
                filters = "ensembl_exon_id",
                values = exons,
                mart = ensembl)
    uniqueTranscripts <- unique(bm$ensembl_transcript_id)
    for (t in uniqueTranscripts){
      exons <- subset(bm, ensembl_transcript_id==t)$ensembl_exon_id
      traces$trace_annotation[, eval(t) := ifelse(exon_id %in% exons,cluster,NA)]
    }
  }
  dt <- subset(traces$trace_annotation,protein_id==protein)
  setkeyv(dt, c("protein_id","PeptidePositionStart"))
  dt[,PeptidePositionStartRank := seq_len(.N), by="protein_id"]
  cbPalette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7","#756bb1","#1c9099")
  dt$cluster <- as.factor(dt$cluster)
  if (PDF){
    if ("genomic_coord" %in% names(traces)) {
      pdf(paste0(protein,"_sequence_cluster.pdf"),width=10,height=5)
    } else {
      pdf(paste0(protein,"_sequence_cluster.pdf"),width=10,height=3)
    }
  }
  if (closeGaps) {
    q <- ggplot(dt,aes(x=PeptidePositionStartRank,
                       y=1,
                       fill=cluster)) +
      geom_bar(stat="identity") + theme_classic() +
      theme(axis.text = element_blank(),
            axis.ticks = element_blank(),
            axis.title= element_blank(),
            axis.line = element_blank()) +
      theme(legend.position="bottom") +
      scale_fill_manual(values=cbPalette) 
    if ("genomic_coord" %in% names(traces)) {
      e <- ggplot(dt,aes(x=PeptidePositionStartRank,
                         y=1,
                         fill=exon_id)) +
        geom_bar(stat="identity") + theme_classic() +
        theme(axis.text = element_blank(),
              axis.ticks = element_blank(),
              axis.title= element_blank(),
              axis.line = element_blank()) +
        theme(legend.position="top") +
        scale_fill_manual(values = getPalette(length(unique(dt$exon_id))))
      f <- ggarrange(e, q, 
                     labels = c("", ""),
                     ncol = 1, nrow = 2)
      print(annotate_figure(f, fig.lab = paste0(protein," : ",unique(dt$Gene_names))))
      
      transcriptAnn <- subset(dt, select=c("id","cluster","PeptidePositionStartRank",uniqueTranscripts))
      transcriptAnn.m <- melt(transcriptAnn, id.vars = c("id","cluster","PeptidePositionStartRank"), variable.name = "transcript")
      transcriptAnn.m$value <- factor(transcriptAnn.m$value)
      #tfPalette <- c("white", "#0072B2")
      t <-  ggplot(transcriptAnn.m,aes(x=PeptidePositionStartRank,
                    y=1,
                    fill=value)) +
        geom_bar(stat="identity") + theme_classic() +
        theme(axis.text = element_blank(),
              axis.ticks = element_blank(),
              axis.title= element_blank(),
              axis.line = element_blank()) +
        theme(legend.position="bottom") +
        scale_fill_manual(values=cbPalette, na.value="white") +
        facet_wrap(. ~ transcript, ncol = 1)
      print(t)
    } else {
      print(q + ggtitle(paste0(protein," : ",unique(dt$Gene_names)))) 
    }
  } else {
    q <- ggplot(dt,aes(x=PeptidePositionStart,
                       y=1,
                       fill=cluster)) +
      geom_bar(stat="identity") + theme_classic() +
      theme(axis.text = element_blank(),
            axis.ticks = element_blank(),
            axis.title= element_blank(),
            axis.line = element_blank()) +
      theme(legend.position="bottom") +
      scale_fill_manual(values=cbPalette) 
    if ("genomic_coord" %in% names(traces)) {
      e <- ggplot(dt,aes(x=PeptidePositionStart,
                         y=1,
                         fill=exon_id)) +
        geom_bar(stat="identity") + theme_classic() +
        theme(axis.text = element_blank(),
              axis.ticks = element_blank(),
              axis.title= element_blank(),
              axis.line = element_blank()) +
        theme(legend.position="top") +
        scale_fill_manual(values = getPalette(length(unique(dt$exon_id))))
      f <- ggarrange(e, q, 
                     labels = c("", ""),
                     ncol = 1, nrow = 2)
      print(annotate_figure(f, fig.lab = paste0(protein," : ",unique(dt$Gene_names))))
      
      transcriptAnn <- subset(dt, select=c("id","cluster","PeptidePositionStart",uniqueTranscripts))
      transcriptAnn.m <- melt(transcriptAnn, id.vars = c("id","cluster","PeptidePositionStart"), variable.name = "transcript")
      transcriptAnn.m$value <- factor(transcriptAnn.m$value)
      #tfPalette <- c("white", "#0072B2")
      t <-  ggplot(transcriptAnn.m,aes(x=PeptidePositionStart,
                                       y=1,
                                       fill=value)) +
        geom_bar(stat="identity") + theme_classic() +
        theme(axis.text = element_blank(),
              axis.ticks = element_blank(),
              axis.title= element_blank(),
              axis.line = element_blank()) +
        theme(legend.position="bottom") +
        scale_fill_manual(values=cbPalette, na.value="white") +
        facet_wrap(. ~ transcript, ncol = 1)
      print(t)
    } else {
      print(q + ggtitle(paste0(protein," : ",unique(dt$Gene_names)))) 
    }
  }
  if (PDF){
    dev.off()
  }
}
