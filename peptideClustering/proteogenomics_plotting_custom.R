plotPeptidesInGenome2 <- function (traces, gene_id, ensdb = NULL, colorMap = NULL, highlight = NULL, 
          plot = T, ...) 
{
  if (!is.null(ensdb)) {
    edb <- validateEnsDB(ensdb)
    ## TODO explore truncate.gaps to make the exons bigger
    ## see https://bioconductor.org/packages/release/bioc/vignettes/biovizBase/inst/doc/intro.pdf
    gene_subset <- GeneIdFilter(gene_id)
    p_gm <- try(ggbio::autoplot(edb, which = gene_subset) + theme_bw(), silent = T)
    if(class(p_gm) != "try-error"){
      tr <- tracks(Gene_model = p_gm, heights = 4)
    }
  }
  tr_sub <- subset.traces(traces, gene_id, "gene_id")
  if (is.null(tr_sub[["genomic_coord"]])) {
    message("No genomic Coordinates found. Trying to fetch them...")
    tr_sub <- annotateGenomicCoordinates(tr_sub, db = edb) 
  }
  pepRanges <- tr_sub$genomic_coord
  gr <- unlist(pepRanges)
  mcols(gr)$id <- names(gr)
  pepRanges <- relist(gr, pepRanges)
  if (!is.null(highlight)) {
    mcols(pepRanges)$outlier <- gsub("\\(.*?\\)", "", names(pepRanges)) %in% 
      gsub("\\(.*?\\)", "", highlight)
    mcols(pepRanges)$outlier <- sapply(mcols(pepRanges)$outlier, 
                                       ifelse, "a", "b")
    if (!any(mcols(pepRanges)$outlier == "a")) 
      highlight <- NULL
  }
  else {
    mcols(pepRanges)$outlier <- "b"
  }
  if (!is.null(colorMap)) {
    if (!all(unique(names(pepRanges)) %in% names(colorMap))) {
      stop("Invalid colorMap specified. Not all traces to be plotted are contained in the colorMap")
    }
  }
  else {
    colorMap <- createGGplotColMap(names(pepRanges))
  } 
  
  p_pep <- ggplot() + geom_alignment(pepRanges[mcols(pepRanges)$outlier == "b",], aes(fill = id, 
                                                                                      alpha = outlier),color = NA, type = "model", group.selfish = F) + 
    theme_bw() +
    theme(legend.position = "none") + 
    scale_fill_manual(values = colorMap) + 
    scale_alpha_manual(values = c(a = 1, b = 0.3))
  
  if (!is.null(ensdb) & class(p_gm) != "try-error") {
    h <- tr@heights
    tr <- tr + tracks(Peptides = p_pep)
    tr@heights <- c(h,1)
  }else {
    tr <- tracks(Peptides = p_pep, heights = 1)
  }
  
  if(!is.null(highlight)){
    p_pepO <- ggplot() + geom_alignment(pepRanges[mcols(pepRanges)$outlier == "a",], aes(fill = id),
                                        color = NA, type = "model", group.selfish = F) + 
      theme_bw() +
      theme(legend.position = "none") +
      scale_fill_manual(values = colorMap)
    h <- tr@heights
    tr <- tr + tracks(OutlPeptides = p_pepO, heights = 1)
    tr@heights <- c(h,1)
  }
    
  p_tr <- plot_tracesCustom(tr_sub, plot = F, highlight = highlight, colorMap = colorMap, ...)
  tr <- as(tr, "grob")
  if (plot) {
    grid.arrange(p_tr, tr, ncol = 1)
  }
  else {
    return(grid.arrange(p_tr, tr, ncol = 1))
  }
}

## TODO Change the name of this and the function call in plotPeptidesInGenome
## when putting into package
plot_tracesCustom2 <- function (traces, log = FALSE, legend = TRUE, PDF = FALSE, name = "Traces", 
          plot = TRUE, highlight = NULL, highlight_col = NULL) 
{
  .tracesTest(traces)
  traces.long <- toLongFormat(traces$traces)
  traces.long <- merge(traces.long, traces$fraction_annotation, 
                       by.x = "fraction", by.y = "id")
  if (!is.null(highlight)) {
    traces.long$outlier <- gsub("\\(.*?\\)", "", traces.long$id) %in% 
      gsub("\\(.*?\\)", "", highlight)
    if (!any(traces.long$outlier)) 
      highlight <- NULL
  }
  p <- ggplot(traces.long) + geom_line(aes_string(x = "fraction", 
                                                  y = "intensity", color = "id")) + xlab("fraction") + 
    ylab("intensity") + theme_bw() + theme(panel.grid.major = element_blank(), 
                                           panel.grid.minor = element_blank()) + theme(plot.margin = unit(c(1, 
                                                                                                            0.5, 0.5, 0.5), "cm")) + ggtitle(name)
  if (log) {
    p <- p + scale_y_log10("log(intensity)")
  }
  if (!legend) {
    p <- p + theme(legend.position = "none")
  }
  if (!is.null(highlight)) {
    legend_peps <- unique(traces.long[outlier == TRUE, id])
    if (is.null(highlight_col)) {
      p <- p + geom_line(data = traces.long[outlier == 
                                              TRUE], aes_string(x = "fraction", y = "intensity", 
                                                                color = "id"), lwd = 2) + scale_color_discrete(breaks = legend_peps)
    }
    else {
      legend_map <- unique(ggplot_build(p)$data[[1]]$colour)
      names(legend_map) <- unique(p$data$id)
      legend_map[legend_peps] <- highlight_col
      legend_vals <- rep(highlight_col, ceiling(length(legend_peps)/length(highlight_col)))[1:length(legend_peps)]
      p <- p + geom_line(data = traces.long[outlier == 
                                              TRUE], aes_string(x = "fraction", y = "intensity", 
                                                                lty = "id"), color = highlight_col, lwd = 2) + 
        scale_color_manual(values = legend_map, limits = legend_peps)
    }
  }
  
  if ("molecular_weight" %in% names(traces$fraction_annotation)) {
    mwtransform <- getMWcalibration(traces)
    mw <- traces$fraction_annotation$molecular_weight
    # frbreaks <-ggplot_build(p)$layout$panel_ranges[[1]]$x.major_source
    # breaks <- mw[frbreaks]
    breaks <- mw[seq(1,length(mw), length.out = 8)]
    p <- p + scale_x_continuous(sec.axis = sec_axis(trans = mwtransform,
                                                    breaks = breaks))
  }
  if (PDF) {
    pdf(paste0(name, ".pdf"))
  }
  if (plot) {
    plot(p)
  }
  else {
    return(ggplot_gtable(ggplot_build(p)))
  }
  if (PDF) {
    dev.off()
  }
}

getMWcalibration <- function(traces){
  fr_ann <- traces$fraction_annotation
  tr <- lm(log(fr_ann$molecular_weight) ~ fr_ann$id)
  intercept <- as.numeric(tr$coefficients[1])
  slope <- as.numeric(tr$coefficients[2])
  FractionToMW <- as.formula(paste0("~exp(",slope, "* . + ", intercept, ")"))
  return(FractionToMW)
}

getMWbreaks <- function(traces){
  fractions <- traces$trace_annotation$id
  tr <- getMWcalibration(traces)
  M
}


plot_tracesCustom <- function(traces, log = FALSE, legend = TRUE, PDF = FALSE, name = "Traces", 
          plot = TRUE, highlight = NULL, highlight_col = NULL, colorMap = NULL) 
{
  .tracesTest(traces)
  traces.long <- toLongFormat(traces$traces)
  traces.long <- merge(traces.long, traces$fraction_annotation, 
                       by.x = "fraction", by.y = "id")
  if (!is.null(highlight)) {
    traces.long$outlier <- gsub("\\(.*?\\)", "", traces.long$id) %in% 
      gsub("\\(.*?\\)", "", highlight)
    if (!any(traces.long$outlier)) 
      highlight <- NULL
  }
  if (!is.null(colorMap)) {
    if (!all(unique(traces.long$id) %in% names(colorMap))) {
      stop("Invalid colorMap specified. Not all traces to be plotted are contained in the colorMap")
    }
  }
  else {
    colorMap <- createGGplotColMap(unique(traces.long$id))
  }
  p <- ggplot(traces.long) + geom_line(aes_string(x = "fraction", 
                                                  y = "intensity", color = "id")) + xlab("fraction") + 
    ylab("intensity") + theme_bw() + theme(panel.grid.major = element_blank(), 
                                           panel.grid.minor = element_blank()) + theme(plot.margin = unit(c(1, 
                                                                                                            0.5, 0.5, 0.5), "cm")) + ggtitle(name) + scale_color_manual(values = colorMap)
  if (log) {
    p <- p + scale_y_log10("log(intensity)")
  }
  if (!legend) {
    p <- p + theme(legend.position = "none")
  }
  if (!is.null(highlight)) {
    legend_peps <- unique(traces.long[outlier == TRUE, id])
    if (is.null(highlight_col)) {
      p <- p + geom_line(data = traces.long[outlier == 
                                              TRUE], aes_string(x = "fraction", y = "intensity", 
                                                                color = "id"), lwd = 2) + scale_color_manual(values = colorMap, 
                                                                                                             breaks = legend_peps)
    }
    else {
      p <- p + geom_line(data = traces.long[outlier == 
                                              TRUE], aes_string(x = "fraction", y = "intensity", 
                                                                lty = "id"), color = highlight_col, lwd = 2)
    }
  }
  if ("molecular_weight" %in% names(traces$fraction_annotation)) {
    mwtransform <- getMWcalibration(traces)
    mw <- traces$fraction_annotation$molecular_weight
    # frbreaks <-ggplot_build(p)$layout$panel_ranges[[1]]$x.major_source
    # breaks <- mw[frbreaks]
    breaks <- mw[seq(1,length(mw), length.out = 8)]
    p <- p + scale_x_continuous(sec.axis = sec_axis(trans = mwtransform,
                                                    breaks = breaks))
  }
  if (PDF) {
    pdf(paste0(name, ".pdf"))
  }
  if (plot) {
    plot(p)
  }
  else {
    return(ggplot_gtable(ggplot_build(p)))
  }
  if (PDF) {
    dev.off()
  }
}
