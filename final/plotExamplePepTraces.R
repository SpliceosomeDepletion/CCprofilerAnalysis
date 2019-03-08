#' ---
#' title: "Import data to CCprofiler"
#' author: "Isabell Bludau"
#' date: "October 20th, 2018"
#' output:
#'   html_document:
#'     toc: true
#'   html_notebook:
#'     toc: true
#'   pdf_document:
#'     toc: true
#' ---

#' # Overview
#' In this script we import the TRIC output data into CCprofiler and create
#' a tracesList object. The traces are annotated with information from
#' biomaRt and UniProt. The SEC is calibrated according to molecular weight
#' standarts and a design matrix for the dataset is created.
#'
#' # Step-by-step workflow
#'
#' ## Load CCprofiler package and set working directory:
library(devtools)
#install_github("CCprofiler/CCprofiler", ref = "DA_module")
#library(CCprofiler)
if (length(grep("nas21.ethz.ch",getwd()))>0) {
  setwd("~/mysonas/CCprofiler")
  load_all()
  setwd("~/mysonas/PRPF8/analysis/output")
  knitr::opts_knit$set(root.dir = '~/mysonas/PRPF8/analysis/output')
} else {
  setwd("/Volumes/ibludau-1/CCprofiler")
  load_all()
  setwd("/Volumes/ibludau-1/PRPF8/analysis/output")
  knitr::opts_knit$set(root.dir = '/Volumes/ibludau-1/PRPF8/analysis/output')
}

reassignedTraces_list <- readRDS("traces_list_reassignedProteoforms.rds")
reassignedTraces<- readRDS("traces_sum_reassignedProteoforms.rds")
design_matrix <- readRDS("design_matrix.rda")
calibration <- readRDS("calibration.rds")

proteins_with_proteoforms <- unique(reassignedTraces$trace_annotation[n_proteoforms>1]$protein_id)
all_proteins <- unique(reassignedTraces$trace_annotation$protein_id)

write.table(proteins_with_proteoforms,"proteins_with_proteoforms.txt",sep="/t",quote=F, col.names = F, row.names = F)
write.table(all_proteins,"all_proteins.txt",sep="/t",quote=F, col.names = F, row.names = F)

proteins_with_proteoforms_005 <- unique(reassignedTraces$trace_annotation[(n_proteoforms>1) & (proteoform_pval_adj<=0.05)]$protein_id)
write.table(proteins_with_proteoforms_005,"proteins_with_proteoforms_005.txt",sep="/t",quote=F, col.names = F, row.names = F)

rand_idx <- sample(length(all_proteins),length(proteins_with_proteoforms_005),replace = F)

proteins_with_proteoforms_005_random <- unique(reassignedTraces$trace_annotation[rand_idx]$protein_id)
write.table(proteins_with_proteoforms_005_random,"proteins_with_proteoforms_005_random.txt",sep="/t",quote=F, col.names = F, row.names = F)

########### 

david_enrichment <- fread("proteoforms_vs_allProteins_davidEnrichment.txt")
names(david_enrichment) <- gsub(" ", "_", names(david_enrichment))

david_enrichment_sig <- subset(david_enrichment, Category=="UP_SEQ_FEATURE" & Fold_Enrichment >= 1 & Bonferroni <= 0.01)
david_enrichment_sig[, neglog10Bonferroni := -log10(Bonferroni)]
david_enrichment_sig <- david_enrichment_sig[order(Fold_Enrichment)]
cat <- david_enrichment_sig$Term
david_enrichment_sig$Term <- factor(david_enrichment_sig$Term, levels=cat)
david_enrichment_sig[,name:=Term]

q <- ggplot(data=david_enrichment_sig,aes(x=name,y=Fold_Enrichment, fill=neglog10Bonferroni)) +
  geom_bar(stat = "identity") +
  theme_classic() +
  theme(axis.text.y  = element_text(angle = 30, hjust = 1)) +
  labs(fill='-log10BF',x='UP SEQ FEATURE',y='Fold Enrichment') +
  coord_flip() +
  geom_text(aes(label=round(Fold_Enrichment,digits=1)), hjust=1.2,color="white")

pdf("proteoforms_vs_allProteins_davidEnrichment_1percent.pdf",height=6,width=6)
plot(q)
dev.off()

david_enrichment <- fread("proteoforms_005_vs_allProteins_davidEnrichment.txt")
names(david_enrichment) <- gsub(" ", "_", names(david_enrichment))

david_enrichment_sig <- subset(david_enrichment, Category=="UP_SEQ_FEATURE" & Fold_Enrichment >= 1 & Bonferroni <= 0.01)
david_enrichment_sig[, neglog10Bonferroni := -log10(Bonferroni)]
david_enrichment_sig <- david_enrichment_sig[order(Fold_Enrichment)]
cat <- david_enrichment_sig$Term
david_enrichment_sig$Term <- factor(david_enrichment_sig$Term, levels=cat)
david_enrichment_sig[,name:=Term]

q <- ggplot(data=david_enrichment_sig,aes(x=name,y=Fold_Enrichment, fill=neglog10Bonferroni)) +
  geom_bar(stat = "identity") +
  theme_classic() +
  theme(axis.text.y  = element_text(angle = 30, hjust = 1)) +
  labs(fill='-log10BF',x='UP SEQ FEATURE',y='Fold Enrichment') +
  coord_flip() +
  geom_text(aes(label=round(Fold_Enrichment,digits=1)), hjust=1.2,color="white")

pdf("proteoforms_005_vs_allProteins_davidEnrichment_1percent.pdf",height=6,width=6)
plot(q)
dev.off()

david_enrichment <- fread("proteoforms_005_random_vs_allProteins_davidEnrichment.txt")
names(david_enrichment) <- gsub(" ", "_", names(david_enrichment))

david_enrichment_sig <- subset(david_enrichment, Category=="UP_SEQ_FEATURE" & Fold_Enrichment >= 1 & Bonferroni <= 0.01)
david_enrichment_sig[, neglog10Bonferroni := -log10(Bonferroni)]
david_enrichment_sig <- david_enrichment_sig[order(Fold_Enrichment)]
cat <- david_enrichment_sig$Term
david_enrichment_sig$Term <- factor(david_enrichment_sig$Term, levels=cat)
david_enrichment_sig[,name:=Term]

nrow(david_enrichment_sig)
# negative control check passed :-)

########### 

david_enrichment <- fread("proteoforms_vs_allProteins_davidEnrichment.txt")
names(david_enrichment) <- gsub(" ", "_", names(david_enrichment))

david_enrichment_sig <- subset(david_enrichment, Category=="GOTERM_BP_DIRECT" & Fold_Enrichment >= 1 & Bonferroni <= 0.05)
david_enrichment_sig[, neglog10Bonferroni := -log10(Bonferroni)]
david_enrichment_sig <- david_enrichment_sig[order(Fold_Enrichment)]
cat <- david_enrichment_sig$Term
david_enrichment_sig$Term <- factor(david_enrichment_sig$Term, levels=cat)
david_enrichment_sig[,name:=gsub("^GO:.*~","",Term)]

q <- ggplot(data=david_enrichment_sig,aes(x=name,y=Fold_Enrichment, fill=neglog10Bonferroni)) +
  geom_bar(stat = "identity") +
  theme_classic() +
  theme(axis.text.y  = element_text(angle = 30, hjust = 1)) +
  labs(fill='-log10BF',x='GOTERM BP DIRECT',y='Fold Enrichment') +
  coord_flip() +
  geom_text(aes(label=round(Fold_Enrichment,digits=1)), hjust=1.2,color="white")

q

pdf("proteoforms_vs_allProteins_davidEnrichment_5percent_GOTERM_BP_DIRECT.pdf",height=6,width=6)
plot(q)
dev.off()

########### 

test_proteins <- c(#"Q6P2Q9", # PRP8
                   #"Q96DI7", # SNR40
                   #"Q13435", # SF3B2
                   #"O95905", # ECD
                   #"Q9UMS4", # PRP19
                   #"Q9BRX9", # WDR83
                   #"Q96NB3", # ZN830
                   "Q9UJX2", # CDC23 / APC8
                   "Q9UJX5", # APC4 / ANAPC4
                   "P30260", # CDC27
                   "Q9UJX6", # ANAPC2 >> check
                   "Q9UJX5", # ANAPC4
                   "Q13042", # APC6 / CDC16 / ANAPC6
                   "P55072", # VPC
                   "P55036", # PSMD4
                   "P28070", # PSMB4
                   "Q9UNM6", # PSMD13
                   "P51665", # PSMD7
                   "P60900", # PSMA6
                   "O43502", # RAD51C
                   "P38398", # BRCA1 
                   "Q8IY92", # SLX4
                   "Q9BXW9", # FANC >> check
                   "P61978", # HNRNPK >> check
                   "Q96BD6", # SSB1
                   "P14618", # KPYM
                   "P04075", # ALDOA
                   "Q92643", # GPI8
                   "P07195", # LDHB
                   "P60174", # TPI1
                   "P40926", # MDH2
                   "P06733", # ENO1
                   #"P40763", # STAT3
                   #"Q12834", # CDC20 
                   #"Q14674", # ESPL1
                   #"Q8IZT6", # "ASPM"
                   #"Q8WVK7", # SKA2
                   #"Q9Y266", # NUDC
                   #"O43143", # DHX15
                   #"P07195", # LDHB
                   #"P42167"  #LAP2B
                   "O75717" # WD repeat and HMG-box DNA-binding protein 1
                   )

get_id_idx <- function(y){
  y = unlist(strsplit(y,"/"))
  y = y[-1]
  y = gsub("([[:alnum:]]*[-])(.*)","\\2",y)
  paste(sort(y), collapse = ";")}

pdf(paste0("PeptideTraces_proteoformAnn.pdf"))
for (test_protein in test_proteins){
  pepTest <- subset(reassignedTraces_list, trace_subset_ids = test_protein, trace_subset_type = "protein_id")
  test_protein_name <- pepTest[[1]]$trace_annotation$Entry_name[1]
  pepTest <- lapply(pepTest, function(x) {
    x$trace_annotation <- x$trace_annotation[,isoform_id_idx := lapply(isoform_id,get_id_idx), by = "id"]
    return(x)
  })
  class(pepTest) <- "tracesList"
  if(nrow(pepTest$control$trace_annotation) > 0 ){
    #pdf(paste0("PeptideTraces_proteoformAnn_",test_protein,"_",test_protein_name,".pdf"))
    plot(pepTest, design_matrix=design_matrix, log = FALSE, PDF = F,
         name = paste0("PeptideTraces_proteoformAnn_",test_protein,"_",test_protein_name), 
         plot = TRUE, highlight = NULL, highlight_col = NULL,
         colour_by = "proteoform_id")
    plot(pepTest, design_matrix=design_matrix, log = FALSE, PDF = F,
         name = paste0("PeptideTraces_",test_protein,"_",test_protein_name), 
         plot = TRUE, highlight = NULL, highlight_col = NULL, legend=T,
         colour_by = "isoform_id_idx")
    #dev.off()
  }
}
dev.off()

test_proteins <- c(
  "Q9UJX6", # ANAPC2 >> check
  "Q9BXW9", # FANC >> check
  "P61978", # HNRNPK >> check
  "O75717" # WD repeat and HMG-box DNA-binding protein 1
)


subset(pepTest$control$trace_annotation, select=c("Sequence","proteoform_id"))

pdf(paste0("PeptideTraces_proteoformAnn_seqInfo.pdf"))
for (test_protein in test_proteins){
  plotPeptideCluster(reassignedTraces,test_protein)
}
dev.off()


###########

test_proteins <- c(
  "P42574", # caspase 3
  "P55212", # caspase 6
  "P55210", # caspase 7
  "Q14790", # caspase 8
  "P55211", # caspase 9
  "P07858", # cathepsin B
  "P07711", # cathepsin L1
  "O60911", # cathepsin L2
  "P07339", # cathepsin D
  "P10144", # granzyme B
  "P0DJD7", # pepsin A4
  "P0DJD9", # pepsin A5
  "P0DJD8", #pepsin A3
  "P07477" # trypsin 1
)

pdf(paste0("PeptideTraces_proteases.pdf"))
for (test_protein in test_proteins){
  pepTest <- subset(reassignedTraces_list, trace_subset_ids = test_protein, trace_subset_type = "protein_id")
  test_protein_name <- pepTest[[1]]$trace_annotation$Entry_name[1]
  if(nrow(pepTest$control$trace_annotation) > 0 ){
    #pdf(paste0("PeptideTraces_proteoformAnn_",test_protein,"_",test_protein_name,".pdf"))
    plot(pepTest, design_matrix=design_matrix, log = FALSE, PDF = F,
         name = paste0("PeptideTraces_proteoformAnn_",test_protein,"_",test_protein_name), 
         plot = TRUE, highlight = NULL, highlight_col = NULL,
         colour_by = "proteoform_id")
    #dev.off()
  }
}
dev.off()


test_proteins <- c(
  "Q14790", # caspase 8
  "Q13158", # FADD
  "O15519", # CFLAR
  "Q15121", #PEA15
  "P19838"
)

pdf(paste0("PeptideTraces_caspase8.pdf"))
for (test_protein in test_proteins){
  pepTest <- subset(reassignedTraces_list, trace_subset_ids = test_protein, trace_subset_type = "protein_id")
  test_protein_name <- pepTest[[1]]$trace_annotation$Entry_name[1]
  if(nrow(pepTest$control$trace_annotation) > 0 ){
    #pdf(paste0("PeptideTraces_proteoformAnn_",test_protein,"_",test_protein_name,".pdf"))
    plot(pepTest, design_matrix=design_matrix, log = FALSE, PDF = F,
         name = paste0("PeptideTraces_proteoformAnn_",test_protein,"_",test_protein_name), 
         plot = TRUE, highlight = NULL, highlight_col = NULL,
         colour_by = "proteoform_id")
    #dev.off()
  }
}
dev.off()


test_proteins <- unique(reassignedTraces$trace_annotation[n_proteoforms>1]$protein_id)[1:300]
pdf(paste0("PeptideTraces_clusters.pdf"))
for (test_protein in test_proteins){
  pepTest <- subset(reassignedTraces_list, trace_subset_ids = test_protein, trace_subset_type = "protein_id")
  plotPeptideCluster(reassignedTraces,test_protein)
}
dev.off()

#>> gap peptide cluster
test_proteins <- c(
"O43399", # : TPD52L2
"P16615", # : ATP2A2 ATP2B !!!!!!!
"Q8TD16", # : BICD2 KIAA0699
"O75116", # : ROCK2 KIAA0619
"P49959", # : MRE11 HNGS1 MRE11A >> no cluster gap, but also no cleavage gap
"Q16891", # : IMMT HMP MIC60 MINOS2 PIG4 PIG52 !!!!
"Q14683", # : SMC1A DXS423E KIAA0178 SB1.8 SMC1 SMC1L1
"Q01082", # : SPTBN1 SPTB2
"Q9H6T3", # : RPAP3
"P15311", # : EZR VIL2
"Q15811", # : ITSN1 ITSN SH3D1A
"P42704", # : LRPPRC LRP130
"Q6P1J9", # : CDC73 C1orf28 HRPT2
"O75717", # : WDHD1 AND1
"Q9NR30", # : DDX21
"P52701", # : MSH6 GTBP
"Q8ND24" # : RNF214
)

pdf(paste0("PeptideTraces_clusterGap.pdf"))
for (test_protein in test_proteins){
  pepTest <- subset(reassignedTraces_list, trace_subset_ids = test_protein, trace_subset_type = "protein_id")
  test_protein_name <- pepTest[[1]]$trace_annotation$Entry_name[1]
  if(nrow(pepTest$control$trace_annotation) > 0 ){
    #pdf(paste0("PeptideTraces_proteoformAnn_",test_protein,"_",test_protein_name,".pdf"))
    plot(pepTest, design_matrix=design_matrix, log = FALSE, PDF = F,
         name = paste0("PeptideTraces_proteoformAnn_",test_protein,"_",test_protein_name), 
         plot = TRUE, highlight = NULL, highlight_col = NULL,
         colour_by = "proteoform_id")
    #dev.off()
  }
}
dev.off()



