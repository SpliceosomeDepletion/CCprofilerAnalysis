# protein quentification
library(devtools)
#setwd("~/mysonas/CCprofiler")
#load_all()
setwd("~/mysonas/PRPF8/analysis/output")
install_github("CCprofiler/CCprofiler", ref = "DA_module")
library(CCprofiler)

protTraces <- readRDS("protTraces.rds")
protTracesSum <- readRDS("protTracesSum.rds")
design_matrix <- readRDS("design_matrix.rda")

complexHypotheses <- readRDS("~/mysonas/html/SECpaper/output_data/corum/complexTargetsPlusDecoys.rds")

test_proteins <- c("Q6P2Q9","O75643","Q15029")
protTest <- subset(protTraces, trace_subset_ids = test_proteins, trace_subset_type = "id")
plot(protTest, design_matrix = design_matrix, log = FALSE, legend = TRUE, PDF = TRUE,
 name = paste0("ProteinTraces_PRPF8_BRR2_EFTUD2 (",paste(test_proteins,collapse="_"),")"), plot = TRUE, highlight = NULL, highlight_col = NULL)

 test_proteins <- c("P35998","P43686","P49720")
 protTest <- subset(protTraces, trace_subset_ids = test_proteins, trace_subset_type = "id")
 plot(protTest, design_matrix = design_matrix, log = FALSE, legend = FALSE, PDF = TRUE,
  name = paste0("ProteinTraces_",paste(test_proteins,collapse="_")), plot = TRUE, highlight = NULL, highlight_col = NULL)

test_proteins <- c("P42167","P42166")
protTest <- subset(protTraces, trace_subset_ids = test_proteins, trace_subset_type = "id")
plot(protTest, design_matrix = design_matrix, log = FALSE, legend = TRUE, PDF = TRUE,
  name = paste0("ProteinTraces_",paste(test_proteins,collapse="_")), plot = TRUE, highlight = NULL, highlight_col = NULL)

test_proteins <- c("Q6P2Q9","Q96DI7")
protTest <- subset(protTraces, trace_subset_ids = test_proteins, trace_subset_type = "id")
plot(protTest, design_matrix = design_matrix, log = F, legend = T, PDF = TRUE,
  name = paste0("ProteinTraces_PRPF8_SNRNP40"), plot = TRUE, highlight = NULL, highlight_col = NULL)

spliceosome_ids <- complexHypotheses$protein_id[grep("spliceosome",complexHypotheses$complex_name)]
protTest <- subset(protTraces, trace_subset_ids = spliceosome_ids, trace_subset_type = "id")
plot(protTest, design_matrix = design_matrix, log = F, legend = F, PDF = TRUE,
  name = paste0("ProteinTraces_Spliceosome"), plot = TRUE, highlight = NULL, highlight_col = NULL)

PRPF8complex_ids <- fread("../../PRPF8complex.csv",header=F)$V1
protTest <- subset(protTraces, trace_subset_ids = PRPF8complex_ids, trace_subset_type = "id")
plot(protTest, design_matrix = design_matrix, log = T, legend = T, PDF = TRUE,
  name = paste0("ProteinTraces_PRPF8complex_ids"), plot = TRUE, highlight = NULL, highlight_col = NULL)

APC_ids <- complexHypotheses$protein_id[grep("APC",complexHypotheses$complex_name)]
protTest <- subset(protTraces, trace_subset_ids = APC_ids, trace_subset_type = "id")
plot(protTest, design_matrix = design_matrix, log = F, legend = T, PDF = TRUE,
  name = paste0("ProteinTraces_APC"), plot = TRUE, highlight = NULL, highlight_col = NULL)

PKM_ids <- "P14618"
protTest <- subset(protTraces, trace_subset_ids = PKM_ids, trace_subset_type = "id")
plot(protTest, design_matrix = design_matrix, log = F, legend = T, PDF = TRUE,
  name = paste0("ProteinTraces_PKM"), plot = TRUE, highlight = NULL, highlight_col = NULL)

plot(protTest$minus, log = F, legend = T, PDF = TRUE,
  name = paste0("ProteinTraces_PKM"), plot = TRUE, highlight = NULL, highlight_col = NULL)



SF3_ids <- c("Q9Y266",complexHypotheses$protein_id[grep("SF3b complex",complexHypotheses$complex_name)])
protTest <- subset(protTraces, trace_subset_ids = SF3_ids, trace_subset_type = "id")
plot(protTest, design_matrix = design_matrix, log = F, legend = T, PDF = TRUE,
  name = paste0("ProteinTraces_SF3"), plot = TRUE, highlight = NULL, highlight_col = NULL)

  plotSMN <- function(traces,
                              design_matrix = NULL,
                              collapse_conditions = FALSE,
                              log=FALSE,
                              legend = TRUE,
                              PDF=FALSE,
                              name="Traces",
                              annotation_label="id",
                              plot = TRUE,
                              highlight=NULL,
                              highlight_col=NULL) {
    if(!is.null(design_matrix)){
      if(!all(design_matrix$Sample_name %in% names(traces))){
        stop("Invalid design matrix")
      }
    }else{
      design_matrix <- data.table(Sample_name = names(traces),
                                  Condition = "",
                                  Replicate = 1:length(traces))
    }
    tracesList <- lapply(names(traces), function(tr){
      res <- toLongFormat(traces[[tr]]$traces)
      res$Condition <- design_matrix[Sample_name == tr, Condition]
      res$Replicate <- design_matrix[Sample_name == tr, Replicate]
      if (annotation_label != "id") {
        if (annotation_label %in% names(traces[[tr]]$trace_annotation)){
          id_mapping <- subset(traces[[tr]]$trace_annotation, select=c(get(annotation_label),id))
          res_merged <- merge(res,id_mapping)
          res_merged[,id:=get(annotation_label)]
          res <- res_merged
        }
      }
      res
    })

    traces_long <- do.call("rbind", tracesList)
    traces_frac <- unique(do.call("rbind", lapply(traces, "[[", "fraction_annotation")))
    traces_frac <- unique(subset(traces_frac, select = names(traces_frac) %in% c("id","molecular_weight")))
    traces_long <- merge(traces_long,traces_frac,by.x="fraction",by.y="id")
    if(!is.null(highlight)){
      traces_long$outlier <- gsub("\\(.*?\\)","",traces_long$id) %in% gsub("\\(.*?\\)","",highlight)
      if(!any(traces_long$outlier)) highlight <- NULL
    }

    p <- ggplot(traces_long) +
      xlab('fraction') +
      ylab('intensity') +
      theme_bw() +
      theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
      theme(plot.margin = unit(c(1,.5,.5,.5),"cm")) +
      ggtitle(name) +
      theme(plot.title = element_text(vjust=19,size=10))
    if(collapse_conditions){
      p <- p + facet_grid(~ Replicate) +
        geom_line(aes_string(x='fraction', y='intensity', color='id', lty = 'Condition'))
    }else{
      p <- p + facet_grid(Condition ~ Replicate) +
        geom_line(aes_string(x='fraction', y='intensity', color='id'))

    }
    if(!is.null(highlight)){
      legend_peps <- unique(traces_long[outlier == TRUE, id])
      if(is.null(highlight_col)){
        if(collapse_conditions){
          p <- p +
            geom_line(data = traces_long[outlier == TRUE], aes_string(x='fraction', y='intensity', color='id', lty = 'Condition'), lwd=2) +
            scale_color_discrete(breaks = legend_peps)
        }else{
          p <- p +
            geom_line(data = traces_long[outlier == TRUE], aes_string(x='fraction', y='intensity', color='id'), lwd=2) +
            scale_color_discrete(breaks = legend_peps)
        }

      }else{
        legend_map <- unique(ggplot_build(p)$data[[1]]$colour)
        names(legend_map) <- unique(p$data$id)
        legend_map[legend_peps] <- highlight_col
        legend_vals <- rep(highlight_col, ceiling(length(legend_peps)/ length(highlight_col)))[1:length(legend_peps)]
        if(collapse_conditions){
          p <- p +
            geom_line(data = traces_long[outlier == TRUE],
                      aes(x=fraction, y=intensity, lty = Condition, group = interaction(Condition, id), color = id),
                       lwd=2) +
            # scale_color_discrete(guide = F)
            scale_color_manual(values = legend_map, limits = legend_peps)
          # guides(lty = FALSE)
          # scale_color_manual(limits = legend_peps, values = rep(highlight_col, length(legend_peps))) +
          # geom_line(aes_string(x='fraction', y='intensity', color='id'))
        }else{
          p <- p +
            geom_line(data = traces_long[outlier == TRUE], aes_string(x='fraction', y='intensity', color = 'id'),
                      lwd=2) +
            # scale_color_discrete(guide = F)
            scale_color_manual(values = legend_map, limits = legend_peps)
          # guides(lty = FALSE)
          # scale_color_manual(limits = legend_peps, values = rep(highlight_col, length(legend_peps))) +
          # geom_line(aes_string(x='fraction', y='intensity', color='id'))
        }

      }
    }


    if (log) {
      p <- p + scale_y_log10('log(intensity)')
    }
    if (!legend) {
      p <- p + theme(legend.position="none")
    }
    if(PDF){
      pdf(paste0(name,".pdf"))
    }
    if(plot){
      plot(p)
    }else{
      return(p)
    }
    if(PDF){
      dev.off()
    }
  }

SMN_ids <- subset(complexHypotheses,complex_id=="1142")$protein_id
protTest <- subset(protTraces, trace_subset_ids = SMN_ids, trace_subset_type = "id")
pdf("ProteinTraces_SMN.pdf",width=8,height=4)
plotSMN(protTest, design_matrix = design_matrix, log = F, legend = T, PDF = FALSE,
  name = paste0("ProteinTraces_SMN"), plot = TRUE, annotation_label="Gene_names", highlight = NULL, highlight_col = NULL)
dev.off()

SF3_ids <- c("Q9Y266","P26368","Q9Y3B4","O75533","P43034")
protTest <- subset(protTraces, trace_subset_ids = SF3_ids, trace_subset_type = "id")
plot(protTest, design_matrix = design_matrix, log = F, legend = T, PDF = TRUE,
  name = paste0("ProteinTraces_NUDC_U2AF2_SF3B1_SF3B6_LIS1"), plot = TRUE, highlight = NULL, highlight_col = NULL)

SF3_ids <- c("Q9Y266","P26368","Q9Y3B4","O75533","P43034","P08238")
protTest <- subset(protTraces, trace_subset_ids = SF3_ids, trace_subset_type = "id")
plot(protTest, design_matrix = design_matrix, log = T, legend = T, PDF = TRUE,
  name = paste0("ProteinTraces_NUDC_U2AF2_SF3B1_SF3B6_LIS1_HSP90"), plot = TRUE, highlight = NULL, highlight_col = NULL)

EIF3_Ribosomal40S_ids <- c(complexHypotheses$protein_id[grep("40S ribosomal subunit, cytoplasmic",complexHypotheses$complex_name)],complexHypotheses$protein_id[grep("1097",complexHypotheses$complex_id)])
Ribosomal40S_ids <- complexHypotheses$protein_id[grep("40S ribosomal subunit, cytoplasmic",complexHypotheses$complex_name)]
EIF3_ids <- complexHypotheses$protein_id[grep("1097",complexHypotheses$complex_id)]

EIF3_Ribosomal40S_traces <- subset(protTraces, trace_subset_ids = EIF3_Ribosomal40S_ids, trace_subset_type = "id")

traces_minus <- EIF3_Ribosomal40S_traces$minus$traces
traces_minus[,type := ifelse(EIF3_Ribosomal40S_traces$minus$trace_annotation$id %in% EIF3_ids, "EIF3", "Ribosomal40S")]
traces_minus[,condition:="Control"]
traces_plus <- EIF3_Ribosomal40S_traces$plus$traces
traces_plus[,type := ifelse(EIF3_Ribosomal40S_traces$plus$trace_annotation$id %in% EIF3_ids, "EIF3", "Ribosomal40S")]
traces_plus[,condition:="PRPF8 depleted"]
EIF3_Ribosomal40S_traces_all <- rbind(traces_minus,traces_plus)
EIF3_Ribosomal40S_traces_all.m <- melt(EIF3_Ribosomal40S_traces_all,id.vars=c("id","type","condition"))
EIF3_Ribosomal40S_traces_all.m$variable <- as.numeric(EIF3_Ribosomal40S_traces_all.m$variable)
EIF3_Ribosomal40S_traces_all.m$condition <- as.factor(EIF3_Ribosomal40S_traces_all.m$condition)
pdf("EIF3_Ribosomal40S_traces_all.pdf",height=10,width=10)
  g <- ggplot(EIF3_Ribosomal40S_traces_all.m,aes(x=variable, y=value, color=type, group=id)) +
        geom_line() +
        facet_wrap(~ condition, nrow=2) +
        theme_classic()
  g
dev.off()


# no real interactors
NUDC_string_interactors <- paste0(c("NUDC","CENPE","PLK1","ZW10","CLIP1","NUDC","NUP85","MAD1L1","ZWILCH","NDEL1","PAFAH1B1","NDE1"),"_HUMAN")
protTest <- subset(protTraces, trace_subset_ids = NUDC_string_interactors, trace_subset_type = "Entry_name")
plot(protTest, design_matrix = design_matrix, log = F, legend = T, PDF = TRUE,
  name = paste0("ProteinTraces_NUDC_string_interactors"), plot = TRUE, highlight = NULL, highlight_col = NULL)
