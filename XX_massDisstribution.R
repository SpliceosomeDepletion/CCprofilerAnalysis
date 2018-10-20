# protein quentification
library(devtools)
setwd("~/mysonas/CCprofiler")
load_all()
setwd("~/mysonas/PRPF8/analysis/output")
#install_github("CCprofiler/CCprofiler", ref = "DA_module")
#library(CCprofiler)

design_matrix <- readRDS("design_matrix.rda")
protTraces <- readRDS("protTraces.rds")


# evaluate protein mass distribution
protTraces_assembly <- annotateMassDistribution(protTraces)
saveRDS(protTraces_assembly,"protTraces_assembly.rds")


diffAssemblyState <- getMassAssemblyChange(protTraces_assembly, design_matrix)
saveRDS(diffAssemblyState,"diffAssemblyState.rds")

fold_change_cutoff=1

plotMassAssemblyChange(diffAssemblyState, FC_cutoff=fold_change_cutoff, PDF=T)


plotMassAssemblyChange_empty <- function(assamblyTest, FC_cutoff=2, name="massAssemblyChange", PDF=FALSE){
  if (PDF) {
    pdf(paste0(name,".pdf"))
  }
  no_change <- nrow(subset(assamblyTest, abs(change) <= FC_cutoff))
  more_assembled <- nrow(subset(assamblyTest, change < -FC_cutoff))
  less_assembled <- nrow(subset(assamblyTest, change > FC_cutoff))
  pie(c(no_change,more_assembled,less_assembled),
      labels="",
      main = paste0(unique(assamblyTest$testOrder)))

  h <- ggplot(assamblyTest, aes(x=change)) + geom_histogram(binwidth = 1) + theme_classic()
  print(h)
  if (PDF) {
    dev.off()
  }
}

plotMassAssemblyChange_empty(diffAssemblyState, FC_cutoff=fold_change_cutoff, PDF=T, name="massAssemblyChange_empty")

assembly_info <- copy(diffAssemblyState)

all_proteins <- subset(assembly_info, !is.na(change))$protein_id
no_change_proteins <- subset(assembly_info, abs(change) <= fold_change_cutoff)$protein_id
more_assembled_proteins <- subset(assembly_info, change < -fold_change_cutoff)$protein_id
less_assembled_proteins <- subset(assembly_info, change > fold_change_cutoff)$protein_id

write.table(all_proteins,paste0("all_proteins.txt"),sep="\t",quote=F,row.names = F,col.names=F)
write.table(no_change_proteins,paste0("no_change_proteins.txt"),sep="\t",quote=F,row.names = F,col.names=F)
write.table(more_assembled_proteins,paste0("more_assembled_proteins.txt"),sep="\t",quote=F,row.names = F,col.names=F)
write.table(less_assembled_proteins,paste0("less_assembled_proteins.txt"),sep="\t",quote=F,row.names = F,col.names=F)


plot.traces <- function(traces,
                        log=FALSE,
                        legend = TRUE,
                        PDF=FALSE,
                        name="Traces",
                        plot = TRUE,
                        highlight=NULL,
                        highlight_col=NULL,
                        type=c("traces","heatmap"),
                        heatmap_id=FALSE) {

  .tracesTest(traces)
  type <- match.arg(type)
  traces.long <- toLongFormat(traces$traces)
  traces.long <- merge(traces.long,traces$fraction_annotation,by.x="fraction",by.y="id")
  if(!is.null(highlight)){
    traces.long$outlier <- gsub("\\(.*?\\)","",traces.long$id) %in% gsub("\\(.*?\\)","",highlight)
    if(!any(traces.long$outlier)) highlight <- NULL
  }

  if (type=="traces") {
    p <- ggplot(traces.long) +
      geom_line(aes_string(x='fraction', y='intensity', color='id')) +
      xlab('fraction') +
      ylab('intensity') +
      theme_bw() +
      theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
      theme(plot.margin = unit(c(1,.5,.5,.5),"cm")) +
      ggtitle(name)#+
      #theme(plot.title = element_text(vjust=19,size=10))
    if (log) {
      p <- p + scale_y_log10('log(intensity)')
    }
    if(!is.null(highlight)){
      legend_peps <- unique(traces.long[outlier == TRUE, id])
      if(is.null(highlight_col)){
        p <- p +
          geom_line(data = traces.long[outlier == TRUE], aes_string(x='fraction', y='intensity', color='id'), lwd=2) +
          scale_color_discrete(breaks = legend_peps)
      }else{
        legend_map <- unique(ggplot_build(p)$data[[1]]$colour)
        names(legend_map) <- unique(p$data$id)
        legend_map[legend_peps] <- highlight_col
        legend_vals <- rep(highlight_col, ceiling(length(legend_peps)/ length(highlight_col)))[1:length(legend_peps)]
        p <- p +
          geom_line(data = traces.long[outlier == TRUE], aes_string(x='fraction', y='intensity', lty = 'id'), color = highlight_col, lwd=2) +
          # scale_color_discrete(guide = F)
          scale_color_manual(values = legend_map, limits = legend_peps)
        # guides(lty = FALSE)
        # scale_color_manual(limits = legend_peps, values = rep(highlight_col, length(legend_peps))) +
        # geom_line(aes_string(x='fraction', y='intensity', color='id'))
      }
    }
  } else if (type == "heatmap") {
    traces.long[, intensity_zscore := (intensity - mean(intensity))/sd(intensity), by=c("id")]
    traces.int <- getIntensityMatrix(traces)
    if (traces$trace_type == "protein") {
      traces.long <- merge(traces.long,traces$trace_annotation,by="id",all.x=T,all.y=F)
      clustM <- hclust(as.dist((1-cor(t(traces.int), use = "pairwise.complete.obs"))/2), method = "ward.D")
      orderN <- row.names(traces.int)[clustM$order]
      traces.long$id = factor(traces.long$id, levels = orderN)
    } else {
      traces.long <- merge(traces.long,traces$trace_annotation,by="id",all.x=T,all.y=F)
      protTraces.t <- proteinQuantification(traces,
                             topN = 2,
                             keep_less = FALSE,
                             rm_decoys = TRUE,
                             use_sibPepCorr = FALSE,
                             use_repPepCorr = FALSE,
                             full_intersect_only = FALSE)
      protTraces.int <- getIntensityMatrix(protTraces.t)
      clustM <- hclust(as.dist((1-cor(t(protTraces.int)))/2), method = "ward.D")
      orderN <- row.names(protTraces.int)[clustM$order]
      traces.long$protein_id = factor(traces.long$protein_id, levels = orderN)
      setorderv(traces.long,cols="protein_id")
    }
    p <-  ggplot(traces.long) +
      geom_tile(aes_string(x='fraction', y='id', fill='intensity_zscore')) +
      xlab('fraction') +
      ylab(traces$trace_type) +
      scale_fill_gradient(name = 'intensity z-score') +
      theme_bw() +
      theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
      theme(plot.margin = unit(c(1,.5,.5,.5),"cm")) #+
      #ggtitle(name)
      if (!(heatmap_id)) {
        p <- p + theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())
      }
  } else {
    stop("Type not recognized. Type should be either traces or heatmap.")
  }
  if (!legend) {
    p <- p + theme(legend.position="none")
  }

  if ("molecular_weight" %in% names(traces$fraction_annotation)) {
    p2 <- p
    p <- p + scale_x_continuous(name="fraction",
                                breaks=seq(min(traces$fraction_annotation$id),
                                           max(traces$fraction_annotation$id),10),
                                labels=seq(min(traces$fraction_annotation$id),
                                           max(traces$fraction_annotation$id),10))
    p2 <- p2 + scale_x_continuous(name="molecular weight (kDa)",
                                  breaks=seq(min(traces$fraction_annotation$id),max(traces$fraction_annotation$id),10),
                                  labels=round(traces$fraction_annotation$molecular_weight,digits=0)[seq(1,length(traces$fraction_annotation$id),10)]
    )
    ## extract gtable
    g1 <- ggplot_gtable(ggplot_build(p))
    g2 <- ggplot_gtable(ggplot_build(p2))
    ## overlap the panel of the 2nd plot on that of the 1st plot
    pp <- c(subset(g1$layout, name=="panel", se=t:r))

    g <- gtable_add_grob(g1, g2$grobs[[which(g2$layout$name=="panel")]], pp$t, pp$l, pp$b, pp$l)
    ## steal axis from second plot and modify
    ia <- which(g2$layout$name == "axis-b")
    ga <- g2$grobs[[ia]]
    ax <- ga$children[[2]]
    ## switch position of ticks and labels
    ax$heights <- rev(ax$heights)
    ax$grobs <- rev(ax$grobs)
    ## modify existing row to be tall enough for axis
    g$heights[[2]] <- g$heights[g2$layout[ia,]$t]
    ## add new axis
    g <- gtable_add_grob(g, ax, 2, 4, 2, 4)
    ## add new row for upper axis label
    g <- gtable_add_rows(g, g2$heights[1], 1)
    ## steal axis label from second plot
    ia2 <- which(g2$layout$name == "xlab-b")
    ga2 <- g2$grobs[[ia2]]
    g <- gtable_add_grob(g,ga2, 3, 4, 2, 4)

    if(PDF){
      pdf(paste0(name,".pdf"))
    }
    if(plot){
      grid.draw(g)
    }else{
      return(g)
    }

    if(PDF){
      dev.off()
    }
  }else{
    if(PDF){
      pdf(paste0(name,".pdf"))
    }
    if(plot){
      plot(p)
    }else{
      return(ggplot_gtable(ggplot_build(p)))
    }

    if(PDF){
      dev.off()
    }
  }
}


library(grid)
library(gtable)

plot(protTraces$minus,name="protTracesMinus",PDF=TRUE,legend=F)
plot(protTraces$minus,name="protHeatmapMinus",PDF=TRUE,legend=F, type="heatmap")

plot(pepTraces$minus,name="pepTracesMinus",PDF=TRUE,legend=F)
plot(pepTraces$minus,name="pepHeatmapMinus",PDF=TRUE,legend=F, type="heatmap")

plot(protTraces$plus,name="protTracesPlus",PDF=TRUE,legend=F)
plot(protTraces$plus,name="protHeatmapPlus",PDF=TRUE,legend=F, type="heatmap")

plot(pepTraces$plus,name="pepTracesPlus",PDF=TRUE,legend=F)
plot(pepTraces$plus,name="pepHeatmapPlus",PDF=TRUE,legend=F, type="heatmap")
