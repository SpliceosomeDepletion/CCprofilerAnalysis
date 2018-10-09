library(devtools)
#setwd("~/mysonas/CCprofiler")
#load_all()
setwd("~/mysonas/PRPF8/analysis/output")
install_github("CCprofiler/CCprofiler", ref = "DA_module")
library(CCprofiler)

protTraces <- readRDS("protTracesSum.rds")
conditionTraces <- readRDS("protTraces.rds")
complexHypotheses <- readRDS("~/mysonas/html/SECpaper/output_data/corum/complexTargetsPlusDecoys.rds")
complexFeatures <- readRDS("corumComplexFeatures.rda")
hypothesis <- "corum"
calibrationFunctions <- readRDS("calibration.rds")

plotSummarizedMScoverage(hypotheses = complexHypotheses, protTraces, PDF = TRUE)

# filter complex features to only keep features at higher MW than 2x monomer MW
complexFeatures <- filterFeatures(complexFeatures,
                                  complex_ids = NULL,
                                  protein_ids = NULL,
                                  min_feature_completeness = NULL,
                                  min_hypothesis_completeness = NULL,
                                  min_subunits = NULL,
                                  min_peak_corr = NULL,
                                  min_monomer_distance_factor = 2
                                  )

# subset to best feature per hypothesis only
# (secondary features do not need to follow similar stringency)
complexFeaturesBest <- getBestFeatures(complexFeatures)
saveRDS(complexFeaturesBest,"complexFeaturesBest.rds")

complexFeaturesBestFiltered <- scoreFeatures(complexFeaturesBest, FDR=0.05, PDF=T, name="qvalueStats_complexFeatures")
saveRDS(complexFeaturesBestFiltered,paste0("complexFeaturesBestFiltered_005_FDR.rds"))
write.table(complexFeaturesBestFiltered,paste0("complexFeaturesBestFiltered_005_FDR.txt"),sep="\t",quote=F,row.names = F)

# select peak-correlation cutoff for secondary features
# all features are filtered for best feature FDR and secondary feature peak corr cutoff
scoredDataAll <- appendSecondaryComplexFeatures(scoredPrimaryFeatures = complexFeaturesBestFiltered, allFeatures = complexFeatures, peakCorr_cutoff = 0.5)
saveRDS(scoredDataAll, "scoredDataAll_005_FDR.rds")
write.table(scoredDataAll, paste0(hypothesis,"_scoredDataAll_005_FDR.txt"),sep="\t",quote=F,row.names=F)

plotSummarizedComplexes(scoredDataAll, complexHypotheses, protTraces, PDF=TRUE, name="complex_completeness_pie")

plotSummarizedComplexes_emptyPie <- function(complexFeatures,hypotheses,protTraces,PDF=FALSE,name="complex_completeness_pie"){
  targetFeatures <- copy(complexFeatures)
  targetHypotheses <- copy(hypotheses)
  targetFeatures <- targetFeatures[grep("DECOY",targetFeatures$complex_id,invert=TRUE)]
  targetHypotheses <- targetHypotheses[grep("DECOY",targetHypotheses$complex_id,invert=TRUE)]
  proteins_in_targetHypotheses <- unique(targetHypotheses$protein_id)
  proteins_in_traces <- unique(protTraces$traces$id)
  targetHypotheses[,annotated:=1]
  targetHypotheses[,detected:=ifelse(protein_id %in% proteins_in_traces, 1, 0)]
  #targetHypotheses[,protein_collapsed := paste(protein_id,collapse=";"),by=complex_id]
  targetHypotheses[,annotated_collapsed := sum(annotated),by=complex_id]
  targetHypotheses[,detected_collapsed := sum(detected),by=complex_id]
  targetHypotheses[,ms_completeness := detected_collapsed/annotated_collapsed,by=complex_id]
  unique_targetHypotheses <- unique(targetHypotheses,by="complex_id")
  #unique_targetHypotheses_50 <- subset(unique_targetHypotheses,ms_completeness >= 0.5)

  targetFeatures <- getBestFeatures(targetFeatures)
  targetFeatures_min50 <- subset(targetFeatures,(completeness>=0.5) & (completeness<1))
  targetFeatures_100 <- subset(targetFeatures,completeness==1)
  targetFeatures_lower50 <- subset(targetFeatures,completeness<0.5)

  #complexFeatures_noDecoys <- complexFeatures[grep("DECOY",complexFeatures$complex_id,invert=TRUE)]
  complexCompletenessSummary <- data.table(name=c("no co-elution","co-elution feature"),
                                        count=c(
                                          sum(!(unique_targetHypotheses$complex_id %in% targetFeatures$complex_id)),
                                          sum(unique_targetHypotheses$complex_id %in% targetFeatures$complex_id)
                                          )
                                        )


  cbPalette <- c("white", "#56B4E9")
  if(PDF){pdf(gsub("$|\\.pdf$", ".pdf", name))}
    print(pie(x=complexCompletenessSummary$count,labels=paste0(complexCompletenessSummary$name,"\n",complexCompletenessSummary$count),col=cbPalette[1:nrow(complexCompletenessSummary)]))
    print(pie(x=complexCompletenessSummary$count,labels="",col=cbPalette[1:nrow(complexCompletenessSummary)]))
  if(PDF){dev.off()}
}

plotSummarizedComplexes_emptyPie(scoredDataAll, complexHypotheses, protTraces, PDF=TRUE, name="complex_completeness_emptyPie")

summarizeFeatures(scoredDataAll,
                  plot=TRUE,
                  PDF=TRUE,
                  name="feature_summary_scoredDataAll")
summarizeFeatures(complexFeaturesBestFiltered,
                  plot=TRUE,
                  PDF=TRUE,
                  name="feature_summary_complexFeaturesBestFiltered")


targets <- unique(scoredDataAll$complex_id)[1:40]
pdf("complexFeatures_merged.pdf",width=8,height=4)
for(id in targets){
  plotFeatures(
    feature_table=scoredDataAll,
    traces=protTraces,
    feature_id = id,
    calibration=calibrationFunctions,
    peak_area=TRUE,
    onlyBest = FALSE,
    monomer_MW = TRUE,
    PDF=F,
    name=id
    )
}
dev.off()

targets <- unique(scoredDataAll$complex_id)
pdf("complexFeatures_minus.pdf",width=8,height=4)
for(id in targets){
  plotFeatures(
    feature_table=scoredDataAll,
    traces=conditionTraces$minus,
    feature_id = id,
    calibration=calibrationFunctions,
    annotation_label = "Entry_name",
    peak_area=TRUE,
    onlyBest = FALSE,
    monomer_MW = TRUE,
    PDF=F,
    name=paste0("minus_",id)
    )
}
dev.off()

pdf("complexFeatures_plus.pdf",width=8,height=4)
for(id in targets){
  plotFeatures(
    feature_table=scoredDataAll,
    traces=conditionTraces$plus,
    feature_id = id,
    calibration=calibrationFunctions,
    annotation_label = "Entry_name",
    peak_area=TRUE,
    onlyBest = FALSE,
    monomer_MW = TRUE,
    PDF=F,
    name=paste0("plus_",id)
    )
}
dev.off()
