library(devtools)
setwd("~/Desktop/CCprofiler/CCprofiler")
load_all()
setwd("~/Desktop/PRPF8data/")

inputData <- readRDS("non_sec_samples.rds")
inputData <- subset(inputData,select=c("Sequence","FullPeptideName","filename","Intensity"))
inputData[,inputData:=sum(Intensity),by=c("Sequence","FullPeptideName","filename")]
inputData <- unique(inputData)
inputData <- subset(inputData, filename %in% c("heuselm_L171221_002","heuselm_L171221_003"))
inputData[,condition:=ifelse(filename=="heuselm_L171221_002","minus","plus")]

pepTraces <- readRDS("pepTraces_proteoformMulti.rds")

pepTracesSum <- integrateTraceIntensities(pepTraces,
                                          design_matrix = NULL,
                                          integrate_within = NULL,
                                          aggr_fun = "sum")


inputData_annotated <- merge(inputData,pepTracesSum$trace_annotation,all.x=F,all.y=F,by="Sequence")

inputData_proteinQuant <- copy(inputData_annotated)
inputData_proteinQuant[,nCond := .N, by=c("protein_id","id","condition")]
#inputData_proteinQuant <- subset(inputData_proteinQuant,nCond==2)
inputData_proteinQuant[,protIntensity:=sum(Intensity),by=c("protein_id","condition")]
inputData_proteinQuant <- unique(subset(inputData_proteinQuant,select=c("protein_id","condition","protIntensity")))
inputData_proteinQuant[,replicate:="input"]

inputData_proteinQuant[protein_id=="Q6P2Q9"]


#' ## Expression comparison to full sec
globalDiffExp <- readRDS("globalDiffExp.rda")

protein_diff_exp <- globalDiffExp$diffProteins 
protein_diff_exp[,minus:=global_int1_imp]
protein_diff_exp[,plus:=global_int2_imp]

protein_diff_exp_long <- subset(protein_diff_exp,select=c("protein_id","minus","plus"))
protein_diff_exp_long <- melt(protein_diff_exp_long,id.vars = "protein_id")
setnames(protein_diff_exp_long,c("variable","value"),c("condition","protIntensity"))
protein_diff_exp_long[,replicate:="SEC"]

protein_diff_exp_long[protein_id=="Q6P2Q9"]


# alternative sec quant
pepTraces_quantMinus <- pepTraces$minus$traces
pepTraces_quantMinus[,Intensity:=rowSums(.SD,na.rm=T),by="id",.SDcols = 1:(ncol(pepTraces_quantMinus)-1)]
pepTraces_quantMinus <- subset(pepTraces_quantMinus, select=c("id","Intensity"))
pepTraces_quantMinus_annotated <- merge(pepTraces$minus$trace_annotation,pepTraces_quantMinus,by="id")
pepTraces_quantMinus_annotated <- subset(pepTraces_quantMinus_annotated,id %in% pepTraces$plus$traces$id)
pepTraces_quantMinus_annotated[,protIntensity:=sum(Intensity),by=c("protein_id")]
pepTraces_quantMinus_annotated <- unique(subset(pepTraces_quantMinus_annotated,select=c("protein_id","protIntensity")))
pepTraces_quantMinus_annotated[,condition:="minus"]
pepTraces_quantMinus_annotated[,replicate:="sumSEC"]

pepTraces_quantPlus <- pepTraces$plus$traces
pepTraces_quantPlus[,Intensity:=rowSums(.SD,na.rm=T),by="id",.SDcols = 1:(ncol(pepTraces_quantPlus)-1)]
pepTraces_quantPlus <- subset(pepTraces_quantPlus, select=c("id","Intensity"))
pepTraces_quantPlus_annotated <- merge(pepTraces$plus$trace_annotation,pepTraces_quantPlus,by="id")
pepTraces_quantPlus_annotated <- subset(pepTraces_quantPlus_annotated,id %in% pepTraces$minus$traces$id)
pepTraces_quantPlus_annotated[,protIntensity:=sum(Intensity),by=c("protein_id")]
pepTraces_quantPlus_annotated <- unique(subset(pepTraces_quantPlus_annotated,select=c("protein_id","protIntensity")))
pepTraces_quantPlus_annotated[,condition:="plus"]
pepTraces_quantPlus_annotated[,replicate:="sumSEC"]

sumSECquant <- rbind(pepTraces_quantMinus_annotated,pepTraces_quantPlus_annotated)

#' Merge data
data_quant <- rbind(inputData_proteinQuant,protein_diff_exp_long) 

data_sec.c <- dcast(data_quant,protein_id+condition ~ replicate,value.var="protIntensity")
data_sec.c <- na.omit(data_sec.c, cols=seq_along(data_sec.c), invert=FALSE)

library(plyr)

pdf("quant_inputVsSEQ.pdf",width=4,height=3)
p <- ggplot(data_sec.c,aes(x=log(SEC),y=log(input))) + geom_point(alpha=.3) +
  facet_grid(. ~ condition) +
  theme_classic() 
cors <- ddply(data_sec.c, .(condition), summarise, cor2 = round(cor(SEC, input)^2, 4))
p <- p + geom_text(data=cors, aes(label=paste("r2=", cor2, sep="")), x=min(log(data_sec.c[which(SEC > 0)]$SEC))+2, y=max(log(data_sec.c$input)))

print(p)
dev.off()


data_quant <- rbind(sumSECquant,inputData_proteinQuant) 

data_sec.c <- dcast(data_quant,protein_id+condition ~ replicate,value.var="protIntensity")
data_sec.c <- na.omit(data_sec.c, cols=seq_along(data_sec.c), invert=FALSE)

library(plyr)

pdf("quant_inputVsSumSEQ.pdf",width=4,height=3)
p <- ggplot(data_sec.c,aes(x=log(sumSEC),y=log(input))) + geom_point(alpha=.3) +
  facet_grid(. ~ condition) +
  theme_classic() 
cors <- ddply(data_sec.c, .(condition), summarise, cor2 = round(cor(sumSEC, input)^2, 4))
p <- p + geom_text(data=cors, aes(label=paste("r2=", cor2, sep="")), x=min(log(data_sec.c[which(sumSEC > 0)]$sumSEC))+2, y=max(log(data_sec.c$input)))

print(p)
dev.off()

# diff
data_sec.d <- dcast(data_quant,protein_id+replicate ~ condition,value.var="protIntensity")
data_sec.d[,log2FC := log2(minus/plus)]
data_sec.d <- subset(data_sec.d,select=c("protein_id","replicate","log2FC"))
data_sec.d <- dcast(data_sec.d,protein_id ~ replicate ,value.var="log2FC")

p <- ggplot(data_sec.d,aes(x=input,y=SEC)) + geom_point(alpha=.3) +
  theme_classic() +
  geom_vline(xintercept = 1,colour="red") +
  geom_vline(xintercept = -1,colour="red") +
  geom_hline(yintercept = 1,colour="red") +
  geom_hline(yintercept = -1,colour="red") +
  ylim(c(-5,5)) + xlim(c(-5,5))

p
