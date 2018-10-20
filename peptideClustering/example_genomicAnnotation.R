library(devtools)
setwd("~/mysonas/CCprofiler")
load_all()
setwd("~/mysonas/PRPF8/analysis/output")
#install_github("CCprofiler/CCprofiler", ref = "proteoformLocationMapping")
#library(CCprofiler)

############################
############################
############################

traces <- readRDS("pepTraces_proteoformMulti_combined.rds")
design_matrix <- readRDS("design_matrix.rda")

map_table <- readRDS("../../data/proteogenomics/HeLa_protgen_combined_1FPKM_header_mapping_v2.rda")

#' ## Annotate the Leading isoform for each peptide
traces <- annotateLeadingIsoform(traces,
                                 isoform_col = "isoform_id",
                                 output_col = "LeadingIsoform")

#' ## Find the relative position of each peptide in the protein sequence
#' With the leadingIsoform we can use the mapping table that contains the sequences of the
#' fasta file used for the proteomics search to determine the relative position of each peptide
#' within the protein.
traces <- annotateRelativePepPos(traces, map_table, multimatch = "first", verbose = T)
#' We do see quite some instances where multiple sequences are found. This is due
#' to the fact, that sometimes a single isoform id in the mapping table corresponds
#' to multiple protein ids. However in that case they will have the same sequence, so the relative
#' position will be the same in every case.

saveRDS(traces,"pepTraces_proteoformMulti_combined_genomAnnotation.rds")

############################
############################
############################

traces <- readRDS("pepTraces_proteoformMulti_combined_genomAnnotation.rds")

traces_annotated <- evaluateProteoformLocation(traces)

saveRDS(traces_annotated,"pepTraces_proteoformMulti_combined_genomAnnotation_pvals.rds")
#traces_annotated <- readRDS("pepTraces_proteoformMulti_combined_genomAnnotation_pvals.rds")

test_proteins <- c("P22102","O00468","Q9UBF2","O15067","O75822","Q9Y266","P42167","P61978","O95801","P55060") # "P42166",
for (p in test_proteins){
  plotPeptideCluster(traces_annotated,p)
}

lap2 <- traces_annotated$trace_annotation[protein_id=="P42167"]
setkey(lap2,PeptidePositionStart)

hnrpk <- traces_annotated$trace_annotation[protein_id=="P61978"]
setkey(hnrpk,PeptidePositionStart)


############################
############################
############################

# evaluate analysis results #
single_proteoform <- traces_annotated$trace_annotation[proteoform_pval_adj > 0.05]
multi_proteoforms <- traces_annotated$trace_annotation[proteoform_pval_adj <= 0.05]

# number of proteins with >1 proteoform
proteins_single_proteoform <- unique(single_proteoform$protein_id)
length(proteins_single_proteoform)
proteins_multi_proteoforms <- unique(multi_proteoforms$protein_id)
length(proteins_multi_proteoforms)

# number of clusters
table(unique(subset(single_proteoform,select=c("protein_id","n_proteoforms")))$n_proteoforms)
table(unique(subset(multi_proteoforms,select=c("protein_id","n_proteoforms")))$n_proteoforms)

# proteins where all clusters can be explained by genomic location
multi_proteoforms[,max_genomLocation_pval := lapply(.SD,max), by="protein_id", .SDcols="genomLocation_pval_adj"]
multi_proteoforms_allExplainedByLocation <- subset(multi_proteoforms,max_genomLocation_pval <= 0.05)
proteins_allExplainedByLocation <- unique(multi_proteoforms_allExplainedByLocation$protein_id)
length(proteins_allExplainedByLocation)

test_proteins <- proteins_allExplainedByLocation[1:30]
for (p in test_proteins){
  plotPeptideCluster(traces_annotated,p)
}

multi_proteoforms[,min_genomLocation_pval := lapply(.SD,min), by="protein_id", .SDcols="genomLocation_pval_adj"]
multi_proteoforms_oneExplainedByLocation <- subset(multi_proteoforms,min_genomLocation_pval <= 0.05)
proteins_oneExplainedByLocation <- unique(multi_proteoforms_oneExplainedByLocation$protein_id)
length(proteins_oneExplainedByLocation)

test_proteins <- proteins_oneExplainedByLocation[which(!proteins_oneExplainedByLocation %in% proteins_allExplainedByLocation)][1:30]
for (p in test_proteins){
  plotPeptideCluster(traces_annotated,p)
}

# awesome example:https://www.uniprot.org/uniprot/O14578
# example for skipped exon:https://www.uniprot.org/uniprot/O75717

############################
############################
############################

traces_list <- readRDS("pepTracesImpSPCmultipep.rda")

## CSF
test_proteins <- c("P35613")
pepTest <- subset(traces_list, trace_subset_ids = test_proteins, trace_subset_type = "protein_id")
plot(pepTest, design_matrix = design_matrix, log = FALSE, legend = F, PDF = TRUE,
 name = paste0("PeptideTraces_","CSF_P35613"), isoformAnnotation = F, plot = TRUE, highlight = NULL, highlight_col = NULL)

 ## NFKB1A
 test_proteins <- c("P19838")
 pepTest <- subset(traces_list, trace_subset_ids = test_proteins, trace_subset_type = "protein_id")
 plot(pepTest, design_matrix = design_matrix, log = FALSE, legend = F, PDF = TRUE,
  name = paste0("PeptideTraces_","NFKB1A_P19838"), isoformAnnotation = F, plot = TRUE, highlight = NULL, highlight_col = NULL)

## MDM2 = Q00987
test_proteins <- c("Q00987")
## nope

## p21
test_proteins <- c("P38936")
pepTest <- subset(traces_list, trace_subset_ids = test_proteins, trace_subset_type = "protein_id")
plot(pepTest, design_matrix = design_matrix, log = FALSE, legend = F, PDF = TRUE,
 name = paste0("PeptideTraces_","p21_P38936"), isoformAnnotation = F, plot = TRUE, highlight = NULL, highlight_col = NULL)


## "P61978"= HNRPK
test_proteins <- c("P61978")
pepTest <- subset(traces_list, trace_subset_ids = test_proteins, trace_subset_type = "protein_id")
plot(pepTest, design_matrix = design_matrix, log = FALSE, legend = F, PDF = TRUE,
 name = paste0("PeptideTraces_","HNRPK_P61978"), isoformAnnotation = F, plot = TRUE, highlight = NULL, highlight_col = NULL)

## P42166 = LAP2alpha
test_proteins <- c("P42166")
traces_list <- readRDS("pepTracesImpSPCmultipep.rda")
pepTest <- subset(traces_list, trace_subset_ids = test_proteins, trace_subset_type = "protein_id")
plot(pepTest, design_matrix = design_matrix, log = FALSE, legend = F, PDF = TRUE,
 name = paste0("PeptideTraces_","LAP2alpha_P42166"), isoformAnnotation = F, plot = TRUE, highlight = NULL, highlight_col = NULL)

# P42167 = LAP2beta
test_proteins <- c("P42167")
traces_list <- readRDS("pepTracesImpSPCmultipep.rda")
pepTest <- subset(traces_list, trace_subset_ids = test_proteins, trace_subset_type = "protein_id")
plot(pepTest, design_matrix = design_matrix, log = FALSE, legend = F, PDF = TRUE,
 name = paste0("PeptideTraces_","LAP2beta_P42167"), isoformAnnotation = F, plot = TRUE, highlight = NULL, highlight_col = NULL)

## "Q9Y266" = NUDC
test_proteins <- c("Q9Y266")
traces_list <- readRDS("pepTracesImpSPCmultipep.rda")
pepTest <- subset(traces_list, trace_subset_ids = test_proteins, trace_subset_type = "protein_id")
plot(pepTest, design_matrix = design_matrix, log = FALSE, legend = F, PDF = TRUE,
 name = paste0("PeptideTraces_","NUDC_Q9Y266"), isoformAnnotation = F, plot = TRUE, highlight = NULL, highlight_col = NULL)

## Awesome example: "O00468" Agarin with many known isoforms
test_proteins <- c("O00468")
traces_list <- readRDS("pepTracesImpSPCmultipep.rda")
pepTest <- subset(traces_list, trace_subset_ids = test_proteins, trace_subset_type = "protein_id")
plot(pepTest, design_matrix = design_matrix, log = FALSE, legend = F, PDF = TRUE,
 name = paste0("PeptideTraces_","Agarin"), isoformAnnotation = F, plot = TRUE, highlight = NULL, highlight_col = NULL)

## Awesome example: "Q9UBF2" 2 isoforms annotated in uniprot
test_proteins <- c("Q9UBF2")
#traces_list <- readRDS("pepTracesImpSPCmultipep.rda")
pepTest <- subset(traces_list, trace_subset_ids = test_proteins, trace_subset_type = "protein_id")
plot(pepTest, design_matrix = design_matrix, log = FALSE, legend = F, PDF = TRUE,
 name = paste0("PeptideTraces_","Q9UBF2"), isoformAnnotation = F, plot = TRUE, highlight = NULL, highlight_col = NULL)

 test_proteins <- c("P22102")
 #traces_list <- readRDS("pepTracesImpSPCmultipep.rda")
 pepTest <- subset(traces_list, trace_subset_ids = test_proteins, trace_subset_type = "protein_id")
 plot(pepTest, design_matrix = design_matrix, log = FALSE, legend = F, PDF = TRUE,
  name = paste0("PeptideTraces_","GART"), isoformAnnotation = F, plot = TRUE, highlight = NULL, highlight_col = NULL)

test_proteins <- c("Q6P2E9")
#traces_list <- readRDS("pepTracesImpSPCmultipep.rda")
pepTest <- subset(traces_list, trace_subset_ids = test_proteins, trace_subset_type = "protein_id")
plot(pepTest, design_matrix = design_matrix, log = FALSE, legend = F, PDF = TRUE,
 name = paste0("PeptideTraces_","EDC4"), isoformAnnotation = F, plot = TRUE, highlight = NULL, highlight_col = NULL)

 test_proteins <- c("O75717")
 #traces_list <- readRDS("pepTracesImpSPCmultipep.rda")
 pepTest <- subset(traces_list, trace_subset_ids = test_proteins, trace_subset_type = "protein_id")
 plot(pepTest, design_matrix = design_matrix, log = FALSE, legend = F, PDF = TRUE,
  name = paste0("PeptideTraces_","O75717"), isoformAnnotation = F, plot = TRUE, highlight = NULL, highlight_col = NULL)

  test_proteins <- c("O00468")
  #traces_list <- readRDS("pepTracesImpSPCmultipep.rda")
  pepTest <- subset(traces_list, trace_subset_ids = test_proteins, trace_subset_type = "protein_id")
  plot(pepTest, design_matrix = design_matrix, log = FALSE, legend = F, PDF = TRUE,
   name = paste0("PeptideTraces_","O00468"), isoformAnnotation = F, plot = TRUE, highlight = NULL, highlight_col = NULL)


 ############################
 ############################
 ############################



getProbability <- function(n,N,c){
  choose(n,c)*choose(N-n,c)/choose(N,n)
}

testIfPerfectGenomicExplanation <- function(v,n_pep){
  if((max(v)-min(v)+1) == n_pep[1]) {
    perfectExp=1
  } else {
    perfectExp=0
  }
  perfectExp
}

x[, perfectProb := getProbability(n_pepPerProteoform,n_peptides,0)]
x[, perfectExp := lapply(.SD,testIfPerfectGenomicExplanation,n_peptides), by=c("proteoform_id"), .SDcols = "PeptidePositionStartRank"]

multi <- x[n_proteoforms!=1]
multi <- x[proteoform_pval_adj <= 0.05]
multi <- multi[!is.na(proteoform_id)]
#multi <- multi[n_pepPerProteoform > 3]
#multi <- multi[n_peptides > 1.5*n_pepPerProteoform]
multi <- subset(multi,select=c("protein_id","proteoform_id",
  "n_proteoforms","n_peptides","n_pepPerProteoform",
  "n_notOne","n_notOneRandom_1","n_notOneRandom_2","n_notOneRandom_3","n_notOneRandom_4","n_notOneRandom_5",
  "f_notOne","f_notOneRandom_1","f_notOneRandom_2","f_notOneRandom_3","f_notOneRandom_4","f_notOneRandom_5",
  "normsdStartRank","normsdRandomRank_1","normsdRandomRank_2","normsdRandomRank_3","normsdRandomRank_4","normsdRandomRank_5",
  "perfectProb","perfectExp"))
setkey(multi, "proteoform_id")
multi <- unique(multi)

multi[, .(n_notOneRandom_mean = rowMeans(.SD)), by = "proteoform_id", .SDcols=c("n_notOneRandom_1","n_notOneRandom_2","n_notOneRandom_3","n_notOneRandom_4","n_notOneRandom_5"), drop=F]

# no good histogram
#multi[,pval_genomicEvidence := empPvals(stat=normsdStartRank,stat0=normsdRandomRank,pool=TRUE)]
#pdf("genomePosition_pval_normsdRandomRank.pdf")
#  p <- ggplot(multi, aes(x=pval_genomicEvidence)) +
#  geom_histogram(bins = 50) +
#  theme_classic()
#  plot(p)
#dev.off()

multi.n <- melt(multi, id.vars=c("protein_id","proteoform_id","n_proteoforms"),measure.vars=c("n_notOne","n_notOneRandom_1","n_notOneRandom_2","n_notOneRandom_3","n_notOneRandom_4","n_notOneRandom_5"))
pdf("genomePosition_stats_number.pdf")
  p <- ggplot(multi.n, aes(x=value, colour=variable)) +
  geom_histogram(position="identity",fill="white", alpha=0.1,bins = 50) +
  theme_classic()
  plot(p)
dev.off()

multi.f <- melt(multi, id.vars=c("protein_id","proteoform_id","n_proteoforms"),measure.vars=c("f_notOne","f_notOneRandom_1","f_notOneRandom_2","f_notOneRandom_3","f_notOneRandom_4","f_notOneRandom_5"))
pdf("genomePosition_stats_fraction.pdf")
  p <- ggplot(multi.f, aes(x=value, colour=variable)) +
  geom_histogram(position="identity",fill="white", alpha=0.1,bins = 50)+
  theme_classic()
  plot(p)
dev.off()

multi.sd <- melt(multi, id.vars=c("protein_id","proteoform_id","n_proteoforms"),measure.vars=c("normsdStartRank","normsdRandomRank_1","normsdRandomRank_2","normsdRandomRank_3","normsdRandomRank_4","normsdRandomRank_5"))
pdf("genomePosition_stats_normsd.pdf")
  p <- ggplot(multi.sd, aes(x=value, colour=variable)) +
  geom_histogram(position="identity",fill="white", alpha=0.1,bins = 80)+
  theme_classic()
  plot(p)
  q <- ggplot(multi.sd, aes(x=value, colour=variable)) +
  geom_histogram(position="identity",fill="white", alpha=0.1,bins = 80) +
  xlim(c(0.9,3.5))+
  theme_classic()
  plot(q)
dev.off()


#' Map the genomic position of every peptide

#' For finding the genomic position of all peptides we need to know the
#' ENSEMBL protein id of the isoform they originate from. Therefore
#' we map the leading isoform to its corresponding protein id.


traces_annotated <- readRDS("pepTraces_proteoformMulti_combined_genomAnnotation_pvals.rds")
map_table <- readRDS("../../data/proteogenomics/HeLa_protgen_combined_1FPKM_header_mapping_v2.rda")

map_table$IsoformId <- gsub("\\|.*", "", gsub(">.*?\\|","", map_table$header))
isomap <- unique(map_table[,.(LeadingIsoform=IsoformId, LeadingEnsemblProteinSp=protein)])
isomap <- isomap[!duplicated(LeadingIsoform)] # Take the first protein id if the seq is the same

traces_annotated$trace_annotation <- merge(traces_annotated$trace_annotation, isomap,by="LeadingIsoform")
traces_annotated$trace_annotation[, LeadingEnsemblProtein := gsub("_.*", "", LeadingEnsemblProteinSp)]
traces_annotated$trace_annotation <- traces_annotated$trace_annotation[order(id)]

#' Now the genomic position is retrieved with the ensembldb package.
library(EnsDb.Hsapiens.v86)
ensdb <- EnsDb.Hsapiens.v86

test_proteins <- unique(c(unique(traces_annotated$trace_annotation$protein_id)[1:250],
  c("P22102","O00468","Q08211","Q15149",
  "Q9UBF2","O15067","O75822","Q9Y266", "Q00341",
  "P42167","P61978","Q9UDY2","Q9P0K7","Q04637")))

test <- subset(traces_annotated, trace_subset_ids=test_proteins,
  trace_subset_type="protein_id")

traces <- annotateGenomicCoordinates(test, db=ensdb, proteinIdCol="LeadingEnsemblProtein", verbose=T)

saveRDS(traces,"testTracesWithGenomicCoordinates.rds")
# traces <- readRDS("testTracesWithGenomicCoordinates.rds")

# if multiple exons are annotated, the first one is taken
#exonRes <- evaluateExonSupport(traces)
#exonStats <- plotRealVsRandomSwaps(exonRes)

#traces$trace_annotation <- merge(traces$trace_annotation,exonStats,by="protein_id")

############################
############################
############################

traces <- readRDS("pepTraces_proteoformMulti_combined_genomAnnotation_pvals_GenomicCoordinates.rds")
traces_exon_pval <- evaluateExonLocation(traces)
saveRDS(traces_exon_pval,"traces_exon_pval.rds")

############################
############################
############################

proteoformTraces <- readRDS("pepTraces_proteoformMulti.rds")
design_matrix <- readRDS("design_matrix.rda")

proteoformTraces <- proteinQuantification(proteoformTraces,quantLevel="proteoform_id",
  topN = 1000,
  keep_less = TRUE,
  rm_decoys = TRUE,
  use_sibPepCorr = FALSE,
  use_repPepCorr = FALSE,
  full_intersect_only = FALSE)

test_proteins <- c("P22102","O00468","Q9UBF2","O15067","O75822","Q9Y266","P42167","P61978","O95801","P55060") ## "P42166",
for (test_protein in test_proteins){
  protTest <- copy(proteoformTraces)
  protTest$plus$genomic_coord <- NULL
  protTest$minus$genomic_coord <- NULL
  protTest <- subset(protTest, trace_subset_ids = test_protein, trace_subset_type = "protein_id")
  plot(protTest, design_matrix = design_matrix, log = FALSE, legend = TRUE, PDF = TRUE,
   name = paste0("ProteoformTraces (",paste(test_protein,collapse="_"),")"), plot = TRUE, highlight = NULL, highlight_col = NULL)
}

pepTraces <- readRDS("pepTraces_proteoformMulti.rds")
test_proteins <- c("P22102","O00468","Q9UBF2","O15067","O75822","Q9Y266","P42167","P61978","O95801","P55060") ## "P42166",
for (test_protein in test_proteins){
  protTest <- copy(pepTraces)
  protTest$plus$genomic_coord <- NULL
  protTest$minus$genomic_coord <- NULL
  protTest <- subset(protTest, trace_subset_ids = test_protein, trace_subset_type = "protein_id")
  plot(protTest, design_matrix = design_matrix, log = FALSE, legend = FALSE, PDF = TRUE,
   name = paste0("ProteoformTraces (",paste(test_protein,collapse="_"),")"), plot = TRUE, highlight = NULL, highlight_col = NULL)
}
