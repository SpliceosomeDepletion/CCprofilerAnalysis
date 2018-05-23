library(devtools)
#setwd("~/mysonas/CCprofiler")
#load_all()
setwd("~/mysonas/PRPF8/analysis/output")
install_github("CCprofiler/CCprofiler", ref = "DA_module")
library(CCprofiler)

traces_list <- readRDS("pepTracesImpSPCmultipep.rda")
design_matrix <- readRDS("design_matrix.rda")
# pepTracesConsIdsSPC <- readRDS("pepTracesConsIdsSPC.rda")

test_proteins <- c("Q6P2Q9")
pepTest <- subset(traces_list, trace_subset_ids = test_proteins, trace_subset_type = "protein_id")
plot(pepTest, design_matrix = design_matrix, log = FALSE, legend = T, PDF = TRUE,
 name = paste0("PeptideTraces_","PRPF8"), isoformAnnotation = TRUE, plot = TRUE, highlight = NULL, highlight_col = NULL)

test_proteins <- c("Q96DI7")
pepTest <- subset(traces_list, trace_subset_ids = test_proteins, trace_subset_type = "protein_id")
plot(pepTest, design_matrix = design_matrix, log = FALSE, legend = T, PDF = TRUE,
 name = paste0("PeptideTraces_","SNRNP40"), isoformAnnotation = TRUE, plot = TRUE, highlight = NULL, highlight_col = NULL)

test_proteins <- c("P14618")
pepTest <- subset(traces_list, trace_subset_ids = test_proteins, trace_subset_type = "protein_id")
plot(pepTest, design_matrix = design_matrix, log = FALSE, legend = F, PDF = TRUE,
 name = paste0("PeptideTraces_","PKMall"), isoformAnnotation = TRUE, plot = TRUE, highlight = NULL, highlight_col = NULL)

PKM_traces <- subset(traces_list,trace_subset_ids="P14618",trace_subset_type="protein_id")
traces_minus <- PKM_traces$minus$traces
traces_minus[,isoform := PKM_traces$minus$trace_annotation$isoform_id]
traces_minus[,isoform_prot := PKM_traces$minus$trace_annotation$ensembl_protein_id]
traces_minus[,condition:="Control"]
traces_plus <- PKM_traces$plus$traces
traces_plus[,isoform := PKM_traces$plus$trace_annotation$isoform_id]
traces_plus[,isoform_prot := PKM_traces$plus$trace_annotation$ensembl_protein_id]
traces_plus[,condition:="PRPF8 depleted"]
PKM_traces_all <- rbind(traces_minus,traces_plus)
PKM_traces_all[,PKM1:=""]
PKM_traces_all[,PKM2:=""]
PKM_traces_all$PKM1[grep("ENSP00000320171",PKM_traces_all$isoform_prot)] <- "M1"
PKM_traces_all$PKM2[grep("ENSP00000334983",PKM_traces_all$isoform_prot)] <- "M2"
PKM_traces_all[,PKM:=paste0(PKM1,PKM2)]
PKM_traces_all[,PKM1:=NULL]
PKM_traces_all[,PKM2:=NULL]
PKM_traces_all[,isoform_prot:=NULL]
PKM_traces_all[,isoform:=NULL]
PKM_traces_all.m <- melt(PKM_traces_all,id.vars=c("id","condition","PKM"))
PKM_traces_all.m$PKM[PKM_traces_all.m$PKM=="M1M2"] <- "M_shared"
PKM_traces_all.m$variable <- as.numeric(PKM_traces_all.m$variable)
PKM_traces_all.m$condition <- as.factor(PKM_traces_all.m$condition)

PKM_traces_all.m[,line:=paste0(PKM,id)]

pdf("PKM_traces_all.pdf",height=5,width=5)
g <- ggplot(PKM_traces_all.m,aes(x=variable, y=value, color=factor(PKM), group=line)) +
  geom_line() +
  facet_wrap(~ condition,nrow = 2) +
  theme_classic() + scale_colour_manual(values=c("grey","#0085ff","#FFA500"))
g
dev.off()

PKM_traces_M1M2only.m <- subset(PKM_traces_all.m,PKM!="M_shared")
pdf("PKM_traces_M1M2only.pdf",height=5,width=5)
g <- ggplot(PKM_traces_M1M2only.m,aes(x=variable, y=value, color=factor(PKM), group=id)) +
  geom_line() +
  facet_wrap(~ condition,nrow = 2) +
  theme_classic() + scale_colour_manual(values=c("#0085ff", "#FFA500"))
g
dev.off()

PKM_traces_M1M2only.m[,logInt:=log(value)]
pdf("PKM_traces_M1M2only_log.pdf",height=5,width=5)
g <- ggplot(PKM_traces_M1M2only.m,aes(x=variable, y=logInt, color=factor(PKM), group=id)) +
  geom_line() +
  facet_wrap(~ condition,nrow = 2) +
  theme_classic() + scale_colour_manual(values=c("#0085ff", "#FFA500"))
g
dev.off()


PKM_traces_M1M2only.m[,sum:=sum(value),by=c("condition","PKM","id")]

uniqueM1M2only <- unique(PKM_traces_M1M2only.m,by=c("condition","PKM","id"))
uniqueM1M2only[,logSum := log2(sum)]
uniqueM1M2only[,log2FC := diff(logSum), by=c("PKM","id")]

pdf("PKM_FC.pdf",width=2,height=3)
  ggplot(subset(uniqueM1M2only,condition=="PRPF8 depleted"),aes(x=factor(PKM),y=log2FC,colour=PKM)) + geom_point() + scale_colour_manual(values=c("#0085ff", "#FFA500"))
dev.off()

#ggplot(uniqueM1M2only,aes(x=factor(condition),y=logSum,colour=PKM, group=id)) + geom_point() + geom_line() + facet_wrap(~ PKM) + scale_colour_manual(values=c("grey","#0085ff","#FFA500"))

header <- readRDS("~/mysonas/PRPF8/data/proteogenomics/HeLa_protgen_combined_1FPKM_header_mapping_v2.rda")

PKM_header <- subset(header,protein %in% c("ENSP00000320171", "ENSP00000334983"))
PKM_header[,PKM1:=""]
PKM_header[,PKM2:=""]
PKM_header$PKM1[grep("ENSP00000320171",PKM_header$header)] <- "M1"
PKM_header$PKM2[grep("ENSP00000334983",PKM_header$header)] <- "M2"
PKM_header[,PKM:=paste0(PKM1,PKM2)]
PKM_header[,PKM1:=NULL]
PKM_header[,PKM2:=NULL]
PKM_header <- subset(PKM_header, PKM != "")

PKM1_seq <- unique(subset(PKM_traces_M1M2only.m,PKM=="M1")$id)
PKM2_seq <- unique(subset(PKM_traces_M1M2only.m,PKM=="M2")$id)
PKM_all_seq <- unique(subset(PKM_traces_all.m,PKM=="M_shared")$id)

PKM1_full_sequence <- subset(PKM_header,PKM=="M1")$sequence
PKM1_full_sequence_detected <- PKM1_full_sequence
PKM2_full_sequence <- subset(PKM_header,PKM=="M2")$sequence
PKM2_full_sequence_detected <- PKM2_full_sequence

PKM1_match <- unlist(lapply(PKM1_seq,function(x){grepRaw(x, PKM1_full_sequence_detected)}))
PKM1_len <- unlist(lapply(PKM1_seq,nchar))
idx_v <- seq(1,length(PKM1_match))
PKM1_detected <- unique(unlist(lapply(idx_v,function(x){PKM1_match[x]+seq(0,PKM1_len[x]-1)})))
PKM1_seq_dt <- data.table(AS=seq(1,nchar(PKM1_full_sequence_detected)),evidence="no")
PKM1_seq_dt$evidence[PKM1_detected] <- "PKM1"

PKM2_match <- unlist(lapply(PKM2_seq,function(x){grepRaw(x, PKM2_full_sequence_detected)}))
PKM2_len <- unlist(lapply(PKM2_seq,nchar))
idx_v <- seq(1,length(PKM2_match))
PKM2_detected <- unique(unlist(lapply(idx_v,function(x){PKM2_match[x]+seq(0,PKM2_len[x]-1)})))
PKM2_seq_dt <- data.table(AS=seq(1,nchar(PKM2_full_sequence_detected)),evidence="no")
PKM2_seq_dt$evidence[PKM2_detected] <- "PKM2"

all_PKM1_match <- unlist(lapply(PKM_all_seq,function(x){grepRaw(x, PKM1_full_sequence_detected)}))
all_PKM1_len <- unlist(lapply(PKM_all_seq,nchar))
idx_v <- seq(1,length(all_PKM1_match))
all_PKM1_detected <- unique(unlist(lapply(idx_v,function(x){all_PKM1_match[x]+seq(0,all_PKM1_len[x]-1)})))
all_PKM1_detected <- all_PKM1_detected[! all_PKM1_detected %in% PKM1_detected]
PKM1_seq_dt$evidence[all_PKM1_detected] <- "PKM"

all_PKM2_match <- unlist(lapply(PKM_all_seq,function(x){grepRaw(x, PKM2_full_sequence_detected)}))
all_PKM2_len <- unlist(lapply(PKM_all_seq,nchar))
idx_v <- seq(1,length(all_PKM2_match))
all_PKM2_detected <- unique(unlist(lapply(idx_v,function(x){all_PKM2_match[x]+seq(0,all_PKM2_len[x]-1)})))
all_PKM2_detected <- all_PKM2_detected[! all_PKM2_detected %in% PKM2_detected]
PKM2_seq_dt$evidence[all_PKM2_detected] <- "PKM"


pdf("PKM1_sequence_coverage.pdf",width=5,height=1)
ggplot(PKM1_seq_dt,aes(x=AS,y=1,fill=evidence)) + geom_bar(stat="identity",alpha=1) + theme_classic() +
  theme(axis.text = element_blank(),axis.ticks = element_blank(),axis.title= element_blank(), axis.line = element_blank()) +
  theme(legend.position="bottom") + scale_fill_manual(values=c("grey80","grey30", "#0085ff"))
dev.off()

pdf("PKM2_sequence_coverage.pdf",width=5,height=1)
ggplot(PKM2_seq_dt,aes(x=AS,y=1,fill=evidence)) + geom_bar(stat="identity",alpha=1) + theme_classic() +
  theme(axis.text = element_blank(),axis.ticks = element_blank(),axis.title= element_blank(), axis.line = element_blank()) +
  theme(legend.position="bottom") + scale_fill_manual(values=c("grey80","grey30", "#FFA500"))
dev.off()

##

test_proteins <- c("P40763")
pepTest <- subset(traces_list, trace_subset_ids = test_proteins, trace_subset_type = "protein_id")
plot(pepTest, design_matrix = design_matrix, log = FALSE, legend = T, PDF = TRUE,
     name = paste0("PeptideTraces_","STAT3"), isoformAnnotation = TRUE, plot = TRUE, highlight = NULL, highlight_col = NULL)

 # mitotic Ashok
 test_proteins <- c("Q9UJX2")
 pepTest <- subset(traces_list, trace_subset_ids = test_proteins, trace_subset_type = "protein_id")
 plot(pepTest, design_matrix = design_matrix, log = FALSE, legend = T, PDF = TRUE,
      name = paste0("PeptideTraces_","APC8"), isoformAnnotation = TRUE, plot = TRUE, highlight = NULL, highlight_col = NULL)

 test_proteins <- c("Q12834")
 pepTest <- subset(traces_list, trace_subset_ids = test_proteins, trace_subset_type = "protein_id")
 plot(pepTest, design_matrix = design_matrix, log = TRUE, legend = T, PDF = TRUE,
      name = paste0("PeptideTraces_","CDC20"), isoformAnnotation = TRUE, plot = TRUE, highlight = NULL, highlight_col = NULL)

 test_proteins <- c("Q14674")
 pepTest <- subset(traces_list, trace_subset_ids = test_proteins, trace_subset_type = "protein_id")
 plot(pepTest, design_matrix = design_matrix, log = FALSE, legend = T, PDF = TRUE,
      name = paste0("PeptideTraces_","Separin"), isoformAnnotation = TRUE, plot = TRUE, highlight = NULL, highlight_col = NULL)

 test_proteins <- c("Q8IZT6") # "ASPM"
 pepTest <- subset(traces_list, trace_subset_ids = test_proteins, trace_subset_type = "protein_id")

 test_proteins <- c("Q8WVK7")
 pepTest <- subset(traces_list, trace_subset_ids = test_proteins, trace_subset_type = "protein_id")
 plot(pepTest, design_matrix = design_matrix, log = FALSE, legend = T, PDF = TRUE,
      name = paste0("PeptideTraces_","SKA2"), isoformAnnotation = TRUE, plot = TRUE, highlight = NULL, highlight_col = NULL)

 test_proteins <- c("P49450")
 pepTest <- subset(traces_list, trace_subset_ids = test_proteins, trace_subset_type = "protein_id")

 test_proteins <- c("Q9Y266")
 pepTest <- subset(traces_list, trace_subset_ids = test_proteins, trace_subset_type = "protein_id")
 plot(pepTest, design_matrix = design_matrix, log = FALSE, legend = T, PDF = TRUE,
      name = paste0("PeptideTraces_","NUDC"), isoformAnnotation = TRUE, plot = TRUE, highlight = NULL, highlight_col = NULL)

test_proteins <- c("O43143")
pepTest <- subset(traces_list, trace_subset_ids = test_proteins, trace_subset_type = "protein_id")
plot(pepTest, design_matrix = design_matrix, log = FALSE, legend = T, PDF = TRUE,
     name = paste0("PeptideTraces_","hPrp43"), isoformAnnotation = TRUE, plot = TRUE, highlight = NULL, highlight_col = NULL)

test_proteins <- c("P07195")
pepTest <- subset(traces_list, trace_subset_ids = test_proteins, trace_subset_type = "protein_id")
plot(pepTest, design_matrix = design_matrix, log = FALSE, legend = T, PDF = TRUE,
     name = paste0("PeptideTraces_","LDHB"), isoformAnnotation = TRUE, plot = TRUE, highlight = NULL, highlight_col = NULL)

test_proteins <- c("P42167")
pepTest <- subset(traces_list, trace_subset_ids = test_proteins, trace_subset_type = "protein_id")
plot(pepTest, design_matrix = design_matrix, log = FALSE, legend = T, PDF = TRUE,
 name = paste0("PeptideTraces_","LAP2"), isoformAnnotation = TRUE, plot = TRUE, highlight = NULL, highlight_col = NULL)
