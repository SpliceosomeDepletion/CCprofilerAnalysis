complexHypotheses <- readRDS("~/mysonas/html/SECpaper/output_data/corum/complexTargetsPlusDecoys.rds")
APC_ids <- complexHypotheses$protein_id[grep("APC",complexHypotheses$complex_name)]
protTest <- subset(protTracesTopAll, trace_subset_ids = APC_ids, trace_subset_type = "id")
plot(protTest, design_matrix = design_matrix, log = F, legend = T, PDF = TRUE,
  name = paste0("ProteinTracesTopAll_APC"), plot = TRUE, highlight = NULL, highlight_col = NULL)

  SMN_ids <- subset(complexHypotheses,complex_id=="1142")$protein_id
  protTest <- subset(protTracesTopAll, trace_subset_ids = SMN_ids, trace_subset_type = "id")
  plot(protTest, design_matrix = design_matrix, log = F, legend = T, PDF = TRUE,
    name = paste0("ProteinTracesTopAll_SMN"), plot = TRUE, highlight = NULL, highlight_col = NULL)

test_proteins <- c("Q6P2Q9","O75643","Q15029")
    protTest <- subset(protTracesTopAll, trace_subset_ids = test_proteins, trace_subset_type = "id")
    plot(protTest, design_matrix = design_matrix, log = F, legend = T, PDF = TRUE,
      name = paste0("ProteinTracesTopAll_PRPF8_BRR2_EFTUD2 (",paste(test_proteins,collapse="_"),")"), plot = TRUE, highlight = NULL, highlight_col = NULL)

#######

prot_tracesList <- protTracesTopAll

intMats <- lapply(prot_tracesList, getIntensityMatrix)
allIds <- unique(do.call("c", lapply(intMats, rownames)))

corrPairs <- combn(names(prot_tracesList), m = 2)

cross_corr <- apply(corrPairs, 2, function(pair){
    # print(pair)
    tr1 <- intMats[[pair[1]]]
    tr2 <- intMats[[pair[2]]]
    ids1 <- rownames(tr1)
    ids2 <- rownames(tr2)
    common_ids <- intersect(ids1, ids2)

    t1 <- tr1[common_ids,]
    t2 <-  tr2[common_ids,]
    # sapply(1:nrow(t1), function(i) cor(t1[i,], t2[i,]))
    cA <- t1 - rowMeans(t1)
    cB <- t2 - rowMeans(t2)
    sA <- sqrt(rowMeans(cA^2))
    sB <- sqrt(rowMeans(cB^2))

    rowMeans(cA * cB) / (sA * sB)
  })

repCorrDTprotTopAll <- data.table(id=row.names(cross_corr),corr=cross_corr[,1])


prot_tracesList <- protTraces

intMats <- lapply(prot_tracesList, getIntensityMatrix)
allIds <- unique(do.call("c", lapply(intMats, rownames)))

corrPairs <- combn(names(prot_tracesList), m = 2)

cross_corr <- apply(corrPairs, 2, function(pair){
    # print(pair)
    tr1 <- intMats[[pair[1]]]
    tr2 <- intMats[[pair[2]]]
    ids1 <- rownames(tr1)
    ids2 <- rownames(tr2)
    common_ids <- intersect(ids1, ids2)

    t1 <- tr1[common_ids,]
    t2 <-  tr2[common_ids,]
    # sapply(1:nrow(t1), function(i) cor(t1[i,], t2[i,]))
    cA <- t1 - rowMeans(t1)
    cB <- t2 - rowMeans(t2)
    sA <- sqrt(rowMeans(cA^2))
    sB <- sqrt(rowMeans(cB^2))

    rowMeans(cA * cB) / (sA * sB)
  })

repCorrDTprot <- data.table(id=row.names(cross_corr),corr=cross_corr[,1])

repCorrDTprot[,type:="protein"]
repCorrDTprotTopAll[,type:="proteinTopAll"]

corrDT <- rbind(repCorrDTprot, repCorrDTprotTopAll)


pdf("repCorrBox_annotation.pdf", height=3,width=3)
  p <- ggplot(corrDT, aes(x=type,y=corr)) +
  geom_violin() +
  geom_boxplot(width=0.1,outlier.size = 0.05) +
  theme_classic() +
  theme(legend.position="bottom") +
  theme(axis.title.x = element_blank(),legend.title =  element_blank()) +
  theme(axis.text.x = element_text(angle = 25, hjust = 1)) +
  scale_colour_manual(values=cbPalette)
  print(p)
dev.off()

pdf("repCorrDist.pdf", height=3,width=3)
  p <- ggplot(corrDT, aes(x=corr, colour=type)) +
  geom_density() +
  theme_classic() +
  theme(legend.position="bottom") +
  theme(axis.title.x = element_blank(),legend.title =  element_blank()) +
  theme(axis.text.x = element_text(angle = 25, hjust = 1))
  print(p)
dev.off()
