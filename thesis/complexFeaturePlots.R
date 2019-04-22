complex_DiffExprComplex <- readRDS("complex_DiffExprComplex.rda")

complex_DiffExprProteoform <- readRDS("complex_DiffExprProteoform.rda")
complex_DiffExprProteoform[, has_proteoform := ifelse(any(grep("_",proteoform_id)),TRUE,FALSE), by = c("complex_id")]

n_proteoform_table <- unique(subset(complex_DiffExprProteoform, select=c("complex_id","apex","has_proteoform")))
complex_DiffExprComplex <- merge(complex_DiffExprComplex, n_proteoform_table, by=c("complex_id","apex"))

noProteoforms_DiffExprComplex <- subset(complex_DiffExprComplex, has_proteoform==FALSE)
withProteoforms_DiffExprComplex <- subset(complex_DiffExprComplex, has_proteoform==TRUE)

getCount <- function(data){
  n_proteins <- length(unique(data$complex_id))
  n_differential_proteins <- length(unique(data[pBHadj<0.05][abs(medianLog2FC)>1]$complex_id))
  up <- unique(data[pBHadj < 0.05][medianLog2FC > 1]$complex_id)
  down <- unique(data[pBHadj < 0.05][medianLog2FC < -1]$complex_id)
  n_both <- length(intersect(up,down))
  n_up <- length(up[!up %in% down])
  n_down <- length(down[!down %in% up])
  n_unchanged <- n_proteins-n_differential_proteins
  data.table(
    level=c("complexes","complexes","complexes","complexes"),
    quant=c("unchanged","up & down","up","down"),
    count=c(n_unchanged, n_both, n_up, n_down)
  )
}

complexStats_noProteoforms <- getCount(noProteoforms_DiffExprComplex)
complexStats_noProteoforms[,proteoform:=1]
complexStats_withProteoforms <- getCount(withProteoforms_DiffExprComplex)
complexStats_withProteoforms[,proteoform:=2]

proteinFeatureDiffStats <- do.call("rbind", list(complexStats_noProteoforms, 
                                                 complexStats_withProteoforms))

proteinFeatureDiffStats$quant <- factor(proteinFeatureDiffStats$quant, levels = c("unchanged","up & down","up","down"))
proteinFeatureDiffStats$proteoform <- factor(proteinFeatureDiffStats$proteoform, levels = c(2,1))

pdf("complexDiffStats.pdf", width=5, height=3)
ggplot(proteinFeatureDiffStats, aes(x=quant, y=count, fill=proteoform)) + 
  geom_col() + 
  theme_classic() +
  theme(axis.title.x=element_blank()) +
  scale_fill_manual(values=c("#3288BD","#999999"), labels = c(expression(phantom(x)>= "2"),"1"), name="N proteoforms") #+
  #theme(legend.position = "bottom")
dev.off()

all_proteins_noProteoforms <- sum(proteinFeatureDiffStats[proteoform == 1][level=="proteins"]$count)
diff_proteins_noProteoforms <- sum(proteinFeatureDiffStats[proteoform == 1][level=="proteins"][quant != "unchanged"]$count)

all_proteins_noProteoforms
diff_proteins_noProteoforms
diff_proteins_noProteoforms/all_proteins_noProteoforms

all_proteins_withProteoforms <- sum(proteinFeatureDiffStats[proteoform == 2][level=="proteins"]$count)
diff_proteins_withProteoforms <- sum(proteinFeatureDiffStats[proteoform == 2][level=="proteins"][quant != "unchanged"]$count)

all_proteins_withProteoforms
diff_proteins_withProteoforms
diff_proteins_withProteoforms/all_proteins_withProteoforms

###########################
###########################
###########################

#up
length(unique(complex_DiffExprComplex[pBHadj < 0.05][medianLog2FC < -1]$complex_id))
#down
length(unique(complex_DiffExprComplex[pBHadj < 0.05][medianLog2FC > 1]$complex_id))

###########################
###########################
###########################

complex_DiffExprComplex[, n_features := .N, by = c("complex_id")]

complex_DiffExprComplex_countDT <- unique(subset(complex_DiffExprComplex, select=c("complex_id","n_features")))

pdf("complex_nFeatures.pdf", width=4, height=4)
q <- ggplot(complex_DiffExprComplex_countDT,aes(x=n_features)) +
  stat_bin(binwidth=1) +
  stat_bin(binwidth=1, geom="text", aes(label=..count..), vjust=-0.5) +
  labs(x="N co-elution features",y="N complexes") +
  theme_classic() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
print(q)
dev.off()

###########################
###########################
###########################

complex_DiffExprProteoform[, is_proteoform := ifelse(any(grep("_",proteoform_id)),TRUE,FALSE), by = c("feature_id")]
proteoform_includingComplexes <- complex_DiffExprProteoform[has_proteoform==T][is_proteoform==T]

proteoform_includingComplexes[,n_proteoforms_total := length(unique(proteoform_id)), by=c("feature_id")]
proteoform_includingComplexes[,n_proteoforms_inFeature := length(unique(proteoform_id)), by=c("feature_id","complex_id","apex")]

proteoform_includingComplexes[n_proteoforms_total!=n_proteoforms_inFeature]

proteoform_includingComplexes[,max_proteoform_FC:=max(medianLog2FC), by=c("feature_id","complex_id","apex")]
proteoform_includingComplexes[,min_proteoform_FC:=min(medianLog2FC), by=c("feature_id","complex_id","apex")]

proteoform_includingComplexes[,diff_proteoform_FC:=max_proteoform_FC-min_proteoform_FC]

proteoform_includingComplexes_countDT <- unique(subset(proteoform_includingComplexes, select=c("feature_id","complex_id","apex","diff_proteoform_FC")))

interestingIDs <- unique(proteoform_includingComplexes_countDT[diff_proteoform_FC>2]$feature_id)

reassignedTraces_list <- readRDS("traces_list_reassignedProteoforms.rds")
design_matrix <- readRDS("design_matrix.rda")
calibrationFunctions <- readRDS("calibration.rds")
reassignedFeatures <- readRDS("reassignedProteoformFeatures.rds")
protTraces <- readRDS("protein_traces_list.rds")
scoredDataAll <- readRDS("scoredDataAll.rds")

pdf("proteoformSpecificComplexes_interestingIDs.pdf", height=5, width=8)
for (id in interestingIDs) {
  complex <- proteoform_includingComplexes_countDT[feature_id==id]$complex_id[1]
  sub <- subset(reassignedTraces_list, trace_subset_ids = id, trace_subset_type = "protein_id")
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
  plotFeatures(
    feature_table=scoredDataAll,
    traces=protTraces,
    design_matrix = design_matrix,
    feature_id = complex,
    calibration=calibrationFunctions,
    annotation_label="Entry_name",
    peak_area=TRUE,
    onlyBest = FALSE,
    monomer_MW = TRUE,
    PDF=F,
    name=id
  )
  
}
dev.off()

id = "Q99613"
complex <- proteoform_includingComplexes_countDT[feature_id==id]$complex_id[1]
sub <- subset(reassignedTraces_list, trace_subset_ids = id, trace_subset_type = "protein_id")

plotPeptideCluster(reassignedTraces,"Q99613",PDF = F)

pdf("proteoformSpecificComplexes_Q99613_peptides.pdf", height=5, width=8)
plotFeatures(feature_table = reassignedFeatures,
             traces = reassignedTraces_list,
             calibration=calibrationFunctions,
             feature_id = id,
             design_matrix=design_matrix,
             annotation_label="proteoform_id",
             colour_by="proteoform_id",
             peak_area = F,
             apex = F,
             legend = T,
             onlyBest = F,
             monomer_MW=T)
dev.off()

pdf("proteoformSpecificComplexes_Q99613_complex.pdf", height=7, width=8)
plotFeatures(
  feature_table=scoredDataAll,
  traces=protTraces,
  design_matrix = design_matrix,
  feature_id = complex,
  calibration=calibrationFunctions,
  annotation_label="Entry_name",
  peak_area=TRUE,
  onlyBest = FALSE,
  monomer_MW = TRUE,
  PDF=F,
  name=id
)
dev.off()

library(ggrepel)
plotVolcano <- function(testResults, 
         highlight=NULL, 
         FC_cutoff=2, 
         pBHadj_cutoff=0.01,
         name="volcanoPlot", 
         PDF=FALSE,
         level=c("feature","global")) {
  level <- match.arg(level)
  if (level=="feature"){
    if (PDF) {
      pdf(paste0(name,".pdf"), height=4, width=4)
    }
    if ("medianLog2FC" %in% names(testResults)) {
      p <- ggplot(testResults, aes(x=medianLog2FC,y=-log10(pBHadj)))
    } else {
      p <- ggplot(testResults, aes(x=log2FC,y=-log10(pBHadj)))
    }
    p <- p +
      geom_point(size=1) +
      theme_classic() +
      geom_hline(yintercept=-log10(pBHadj_cutoff), colour="red", linetype="dashed") +
      geom_vline(xintercept=-log2(FC_cutoff), colour="red", linetype="dashed") +
      geom_vline(xintercept=log2(FC_cutoff), colour="red", linetype="dashed")
    if (! is.null(highlight)){
      if ("feature_id" %in% names(testResults)) {
        sub <- subset(testResults,feature_id %in% highlight)
        col <- "feature_id"
      } else if ("complex_id" %in% names(testResults)) {
        sub <- subset(testResults,complex_id %in% highlight)
        col <- "complex_id"
      } else if (highlight %in% testResults$protein_id) {
        sub <- subset(testResults,protein_id %in% highlight)
        col <- "protein_id"
      } else if (highlight %in% testResults$proteoform_id) {
        sub <- subset(testResults,proteoform_id %in% highlight)
        col <- "proteoform_id"
      } else {
        stop("The testResults do not have the proper format. Input should be the result from testDifferentialExpression.")
      }
      if ("medianLog2FC" %in% names(testResults)) {
        p <- p + geom_point(data=sub, aes(x=medianLog2FC,y=-log10(pBHadj)), colour="red", fill="red", size=3, shape=23) #+
          #geom_text_repel(data=sub, aes(label=get(col)), size=4, vjust=0, hjust=-0.1, colour="red")
      } else {
        p <- p + geom_point(data=sub, aes(x=log2FC,y=-log10(pBHadj)), colour="red", fill="red", size=3, shape=23)#+
         # geom_text_repel(data=sub, aes(label=get(col)), size=4, vjust=0, hjust=-0.1, colour="red")
      }
    }
  } else if (level=="global"){
    if (PDF) {
      pdf(paste0(name,"_",level,".pdf"), height=4, width=4)
    }
    if ("medianLog2FC" %in% names(testResults)) {
      p <- ggplot(testResults, aes(x=global_medianLog2FC,y=-log10(global_pBHadj)))
    } else {
      p <- ggplot(testResults, aes(x=global_log2FC,y=-log10(global_pBHadj)))
    }
    p <- p +
      geom_point(size=1) +
      theme_classic() +
      geom_hline(yintercept=-log10(pBHadj_cutoff), colour="red", linetype="dashed") +
      geom_vline(xintercept=-log2(FC_cutoff), colour="red", linetype="dashed") +
      geom_vline(xintercept=log2(FC_cutoff), colour="red", linetype="dashed")
    if (! is.null(highlight)){
      if ("feature_id" %in% names(testResults)) {
        sub <- subset(testResults,feature_id %in% highlight)
        col <- "feature_id"
      } else if ("complex_id" %in% names(testResults)) {
        sub <- subset(testResults,complex_id %in% highlight)
        col <- "complex_id"
      } else if (highlight %in% testResults$protein_id) {
        sub <- subset(testResults,protein_id %in% highlight)
        col <- "protein_id"
      } else if (highlight %in% testResults$proteoform_id) {
        sub <- subset(testResults,proteoform_id %in% highlight)
        col <- "proteoform_id"
      } else {
        stop("The testResults do not have the proper format. Input should be the result from testDifferentialExpression.")
      }
      if ("global_medianLog2FC" %in% names(testResults)) {
        p <- p + geom_point(data=sub, aes(x=global_medianLog2FC,y=-log10(global_pBHadj)), colour="red", fill="red", size=3, shape=23) #+
         # geom_text_repel(data=sub, aes(label=get(col)), size=4, vjust=0, hjust=-0.1, colour="red")
      } else {
        p <- p + geom_point(data=sub, aes(x=global_log2FC,y=-log10(global_pBHadj)), colour="red", fill="red", size=3, shape=23)#+
          #geom_text_repel(data=sub, aes(label=get(col)), size=4, vjust=0, hjust=-0.1, colour="red")
      }
    }
  } else {
    stop("Level can only be feature or global.")
  }
  print(p)
  if (PDF) {
    dev.off()
  }
}
plotVolcano(complex_DiffExprProteoform, highlight=c("Q99613"), PDF = T, name = "complex_DiffExprComplex_Q99613")


