generateMixtureTraces <- function(traces, trueInteractions=NULL, n_proteins = 2, n_mixtureTraces = 500, peptide_ratio = 0.5, min_peptide_count = 2, seed = 123, n_iterations=10) {
  unique_proteins <- unique(traces$trace_annotation$protein_id)
  max_mixtureTraces <- choose(length(unique_proteins), n_proteins)
  if (max_mixtureTraces < n_mixtureTraces) {
    message(paste0("Can only generate maximally ", length(max_mixtureTraces), " mixture traces."))
  }
  n_mixtureTraces <- min(n_mixtureTraces, ncol(max_mixtureTraces))
  
  set.seed(seed)
  if (is.null(trueInteractions)) {
    protein_combinations_noInteraction_max <- replicate(n_mixtureTraces, sample(unique_proteins, n_proteins, replace = FALSE))
  } else {
    #trueInteractions[,idx:=.I]
    #trueInteractions[,id:= paste(sort(c(a,b)),collapse = ";"), by=idx]
    protein_combinations_noInteraction <- matrix()
    it <- 0
    while((ncol(protein_combinations_noInteraction) < n_mixtureTraces) & (it < n_iterations)) {
      it <- it+1
      protein_combinations <- replicate(n_mixtureTraces*2, sample(unique_proteins, n_proteins, replace = FALSE))
      protein_combinations_id <- lapply(c(1:ncol(protein_combinations)),
                                        function(j){
                                          binary_comb <- combn(protein_combinations[,j], m = 2)
                                          lapply(c(1:ncol(binary_comb)),function(x){paste(sort(binary_comb[,x]),collapse = ";")})
                                        })
      noInteraction_idxs_TF <- lapply(protein_combinations_id,function(x){any(unlist(x) %in% trueInteractions$id)})
      
      noInteraction_idxs <- which(noInteraction_idxs_TF=="FALSE")
      protein_combinations_noInteraction <- protein_combinations[,noInteraction_idxs]
    }
    if(ncol(protein_combinations_noInteraction) < n_mixtureTraces) {
      message(paste0("Maximum number of iterations (",it,") reached. Only ",
                     ncol(protein_combinations_noInteraction)," mixtur proteins could be generated."))
      protein_combinations_noInteraction_max <- protein_combinations_noInteraction
    } else {
      protein_combinations_noInteraction_max <- protein_combinations_noInteraction[,1:n_mixtureTraces]
    }
  }
  
  set.seed(seed)
  traces_subsets <- lapply(c(1:ncol(protein_combinations_noInteraction_max)), function(n){
    proteins <- protein_combinations_noInteraction_max[,n]
    proteins <- sort(proteins)
    peptides <- lapply(proteins, function(p){
      traces_sub <- subset(traces, p, trace_subset_type = "protein_id")
      all_peptides <- unique(traces_sub$trace_annotation$id)
      max_peptide_count <- ceiling(length(all_peptides)*peptide_ratio)
      n_peptides_to_select <- max(min_peptide_count,max_peptide_count)
      peptide_subset <- sample(all_peptides, n_peptides_to_select, replace = FALSE)
      return(peptide_subset)
    })
    traces_subset <- subset(traces, unlist(peptides))
    prot_names <- paste(proteins, collapse = "_")
    traces_subset$trace_annotation[,id := paste(paste(prot_names,id, sep = "_"), protein_id, sep="_")]
    traces_subset$traces[,id := traces_subset$trace_annotation$id]
    traces_subset$trace_annotation[,protein_id := prot_names]
    traces_subset$trace_annotation[,n_mixed_proteins := n_proteins]
    traces_subset$trace_annotation[,is_mixed:=TRUE]
    if ("genomic_coord" %in% names(traces_subset)) {
      traces_subset$genomic_coord <- NULL
    }
    return(traces_subset)
  })
  mixed_traces <- list(traces = do.call(rbind,lapply(traces_subsets, `[[`, "traces")),
                       trace_type = traces$trace_type,
                       trace_annotation = do.call(rbind,lapply(traces_subsets, `[[`, "trace_annotation")),
                       fraction_annotation = traces$fraction_annotation)
  class(mixed_traces) <- "traces"
  .tracesTest(mixed_traces)
  return(mixed_traces)
}

evaluateProteoformSensitivityOfMixtureTraces <- function(traces, pvalue_cutoff = 0.05, by=0.01, plot=TRUE, PDF=FALSE, name="TPR_vs_proteoform_pval_adj") {
  pval_range <- seq(0,1,by)  
  stats <- lapply(pval_range, function(p){
    sub <- subset(traces$trace_annotation, proteoform_pval_adj <= p)
    n_positives <- length(unique(subset(traces$trace_annotation,is_mixed == TRUE)$protein_id))
    n_true_positives <- length(unique(subset(sub,is_mixed==TRUE)$protein_id))
    n_negatives <- length(unique(subset(traces$trace_annotation,is_mixed == FALSE)$protein_id))
    n_false_positives <- length(unique(subset(sub,is_mixed==FALSE)$protein_id))
    TPR <- n_true_positives/n_positives
    FPR <- n_false_positives/n_negatives
    data.table(n_positives=n_positives,n_true_positives=n_true_positives,
               n_negatives=n_negatives,n_false_positives=n_false_positives,
               TPR=TPR,FPR=FPR,proteoform_pval_adj=p)
  })
  statsDT <- do.call(rbind,stats)
  if (plot==TRUE) {
    if (PDF==TRUE) {
      pdf(paste0(name,".pdf"),width = 3,height = 3)
    }
    report <- statsDT[proteoform_pval_adj==pvalue_cutoff]
    p <- ggplot(statsDT,aes(x=proteoform_pval_adj,y=TPR)) + 
      geom_line() + 
      theme_classic() +
      geom_vline(xintercept = pvalue_cutoff, colour = "red", linetype="dashed") +
      annotate("text", x=pvalue_cutoff + 0.1, y=0, label= pvalue_cutoff, colour="red") + 
      ggtitle(paste0("N positives = ", report$n_positives, "\nN true positives = ",report$n_true_positives, "\nTPR = ",report$TPR))
    plot(p)
    if (PDF==TRUE) {
    dev.off()
    }
  }
  return(statsDT)
}


evaluateProteoformClusteringOfMixtureTraces <- function(traces){
  all_proteins <- unique(traces$trace_annotation$protein_id)
  mistakes <- lapply(all_proteins, function(x){
    sub <- subset(traces$trace_annotation,protein_id==x)
    unique_proteoforms <- unique(sub$proteoform_id)
    n_proteoforms <- length(unique_proteoforms)
    unique_proteins <- unique(sub$Entry_name)
    n_proteins <- length(unique_proteins)
    n_mistakes <- sum(unlist(lapply(unique_proteoforms,function(p){
      subSub <- subset(sub,proteoform_id==p)
      nrow(subSub)-max(table(subSub$Entry_name))
    })))
    return(data.table(protein=x,n_proteins=n_proteins,n_proteoforms=n_proteoforms,n_mistakes=n_mistakes))
  })
  mistakes <- do.call(rbind,mistakes)
  return(mistakes)
}
                                                         



