proteoform_DiffExprProtein <- readRDS("proteoform_DiffExprProtein.rda")

proteoform_DiffExprProtein[, n_features := .N, by = c("feature_id")]
proteoform_DiffExprProtein[, n_features := ifelse(n_features>=7, 7, n_features)]

proteoform_DiffExprProtein_countDT <- unique(subset(proteoform_DiffExprProtein, select=c("feature_id","n_features")))

pdf("protein_nFeatures.pdf", width=4, height=4)
q <- ggplot(proteoform_DiffExprProtein_countDT,aes(x=n_features)) +
  stat_bin(binwidth=1) +
  stat_bin(binwidth=1, geom="text", aes(label=..count..), vjust=-0.5) +
  labs(x="N co-elution features",y="N proteins") +
  theme_classic() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  scale_x_discrete(limits=seq(1,7,1), 
                   labels=c("1","2","3","4","5","6",expression(phantom(x) >= "7")))
print(q)
dev.off()

###########################
###########################
###########################

proteoform_DiffExprProteoform <- readRDS("proteoform_DiffExprProteoform.rda")
proteoform_DiffExprProteoform[, n_proteoforms_perProtein := length(unique(proteoform_id)), by = c("feature_id")]
proteoform_DiffExprProteoform[, n_proteoforms_perProtein := ifelse(n_proteoforms_perProtein>=7, 7, n_proteoforms_perProtein)]

proteoform_DiffExprProteoform_countDT <- unique(subset(proteoform_DiffExprProteoform, select = c("feature_id","n_proteoforms_perProtein")))
proteoform_DiffExprProteoform_countDT <- merge(proteoform_DiffExprProteoform_countDT,proteoform_DiffExprProtein_countDT, by="feature_id")

max_proteoforms <- max(proteoform_DiffExprProteoform_countDT$n_proteoforms_perProtein)

percentProteoforms <- lapply(unique(proteoform_DiffExprProteoform_countDT$n_features), function(x){
  sub <- subset(proteoform_DiffExprProteoform_countDT, n_features==x)
  total_features <- nrow(sub)
  p <- c()
  for (i in seq(1,max_proteoforms,1)) {
    proteoforms_perProtein <- nrow(subset(sub, n_proteoforms_perProtein == i))
    p <- c(p,proteoforms_perProtein)
  }
  count <- p
  p <- round(p/total_features*100, digits = 2)
  data.table(n_features = x, n_proteoforms=seq(1,max_proteoforms,1), percent_features=p, featur_count = count)
})
percentProteoforms <- do.call(rbind,percentProteoforms)

all_multiFeatures <- sum(percentProteoforms[n_features > 1]$featur_count)
explained_multiFeatures <- sum(percentProteoforms[n_features > 1][n_proteoforms > 1]$featur_count)
notExplained_multiFeatures <-sum(percentProteoforms[n_features > 1][n_proteoforms == 1]$featur_count)

# Percent of proteins with multiple features for which also >1 proteoform was detected
explained_multiFeatures/all_multiFeatures

percentProteoforms$n_proteoforms <- factor(percentProteoforms$n_proteoforms)

pdf("protein_nFeatures_percentProteoforms.pdf", width=4, height=4)
q <- ggplot(percentProteoforms,aes(x=n_features, y=percent_features, fill=n_proteoforms)) +
  geom_col() +
  labs(x="N co-elution features",y="percent of features explained by N proteoforms") +
  theme_classic() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  scale_x_discrete(limits=seq(1,7,1), 
                   labels=c("1","2","3","4","5","6",expression(phantom(x) >= "7"))) +
  scale_fill_brewer(palette="Spectral", labels=c("1","2","3","4","5","6",expression(phantom(x) >= "7")), name="N proteoforms")
  
print(q)
dev.off()

###########################
###########################
###########################

all_proteins <- unique(proteoform_DiffExprProteoform_countDT$feature_id)
proteins_multipleFeatures <- unique(proteoform_DiffExprProteoform_countDT[n_features>=2]$feature_id)
proteins_multipleProteoforms <- unique(proteoform_DiffExprProteoform_countDT[n_features>=2][n_proteoforms_perProtein>=2]$feature_id)
differential_proteins <- unique(proteoform_DiffExprProtein[(abs(medianLog2FC)>=1) & (pBHadj <= 0.05)]$feature_id)

write.table(all_proteins,"all_proteins.txt",sep="\t",quote=F, row.names = F, col.names = F)
write.table(proteins_multipleFeatures,"proteins_multipleFeatures.txt",sep="\t",quote=F, row.names = F, col.names = F)
write.table(proteins_multipleProteoforms,"proteins_multipleProteoforms.txt",sep="\t",quote=F, row.names = F, col.names = F)
write.table(differential_proteins,"differential_proteins.txt",sep="\t",quote=F, row.names = F, col.names = F)

david_multipleFeatures_vs_multipleProteoforms_GO_MF_BP <- fread("david_multipleFeatures_vs_multipleProteoforms_GO_MF_BP.txt")
david_allProteins_vs_multipleFeatures_GO_MF_BP  <- fread("david_allProteins_vs_multipleFeatures_GO_MF_BP.txt")
david_allProteins_vs_differentialProteins_GO_MF_BP <- fread("david_allProteins_vs_differentialProteins_GO_MF_BP.txt")
david_allProteins_vs_differentialProteins_KEGG <- fread("david_allProteins_vs_differentialProteins_KEGG.txt")

plotEnrichment <- function(david_enrichment, enrichment_cutoff = 1.5, Bonferroni_cutoff = 0.05, name="david_enrichment", PDF=T) {
  names(david_enrichment) <- gsub(" ", "_", names(david_enrichment))
  david_enrichment_sig <- subset(david_enrichment, Fold_Enrichment >= enrichment_cutoff & Bonferroni <= Bonferroni_cutoff)
  david_enrichment_sig[, neglog10Bonferroni := -log10(Bonferroni)]
  cat <- david_enrichment_sig$Term
  david_enrichment_sig$Term <- factor(david_enrichment_sig$Term, levels=cat)
  david_enrichment_sig[,name:=gsub("^GO:.*~","",Term)]
  david_enrichment_sig <- david_enrichment_sig[order(Fold_Enrichment)]
  david_enrichment_sig$name <- factor(david_enrichment_sig$name, levels = david_enrichment_sig$name)
  if(nrow(david_enrichment_sig) > 10) {
    david_enrichment_sig <- david_enrichment_sig[1:10]
  }
  
  q <- ggplot(data=david_enrichment_sig,aes(x=name,y=Fold_Enrichment, fill=neglog10Bonferroni)) +
    geom_bar(stat = "identity") +
    theme_classic() +
    theme(axis.text.y  = element_text(angle = 20, hjust = 1)) +
    labs(fill='-log10Bonferroni',x='Category',y='Fold Enrichment') +
    coord_flip() +
    geom_text(aes(label=round(Fold_Enrichment,digits=1)), hjust=1.2,color="white") 
  if (PDF){
    pdf(paste0(name,"_enrichment_cutoff_",enrichment_cutoff,"_Bonferroni_cutoff_",Bonferroni_cutoff,".pdf"),height=4,width=6)
  }
  plot(q)
  if (PDF){
    dev.off()
  }
}

plotEnrichmentKEGG <- function(david_enrichment, enrichment_cutoff = 1.5, Bonferroni_cutoff = 0.05, name="david_enrichment", PDF=T) {
  names(david_enrichment) <- gsub(" ", "_", names(david_enrichment))
  david_enrichment_sig <- subset(david_enrichment, Fold_Enrichment >= enrichment_cutoff & Bonferroni <= Bonferroni_cutoff)
  david_enrichment_sig[, neglog10Bonferroni := -log10(Bonferroni)]
  cat <- david_enrichment_sig$Term
  david_enrichment_sig$Term <- factor(david_enrichment_sig$Term, levels=cat)
  david_enrichment_sig[,name:=gsub("^hsa.*:","",Term)]
  david_enrichment_sig <- david_enrichment_sig[order(Fold_Enrichment)]
  david_enrichment_sig$name <- factor(david_enrichment_sig$name, levels = david_enrichment_sig$name)
  if(nrow(david_enrichment_sig) > 10) {
    david_enrichment_sig <- david_enrichment_sig[1:10]
  }
  
  q <- ggplot(data=david_enrichment_sig,aes(x=name,y=Fold_Enrichment, fill=neglog10Bonferroni)) +
    geom_bar(stat = "identity") +
    theme_classic() +
    theme(axis.text.y  = element_text(angle = 20, hjust = 1)) +
    labs(fill='-log10Bonferroni',x='KEGG pathway',y='Fold Enrichment') +
    coord_flip() +
    geom_text(aes(label=round(Fold_Enrichment,digits=1)), hjust=1.2,color="white") 
  if (PDF){
    pdf(paste0(name,"_enrichment_cutoff_",enrichment_cutoff,"_Bonferroni_cutoff_",Bonferroni_cutoff,".pdf"),height=3,width=5)
  }
  plot(q)
  if (PDF){
    dev.off()
  }
}
 

plotEnrichment(david_multipleFeatures_vs_multipleProteoforms_GO_MF_BP,
                   name="david_multipleFeatures_vs_multipleProteoforms_GO_MF_BP", 
                   enrichment_cutoff = 1.3, Bonferroni_cutoff = 0.1, PDF=T)

plotEnrichmentKEGG(david_allProteins_vs_differentialProteins_KEGG,
               name="david_allProteins_vs_differentialProteins_KEGG", 
               enrichment_cutoff = 1.5, Bonferroni_cutoff = 0.05, PDF=T)

pdf("david_allProteins_vs_differentialProteins_GO_MF_BP.pdf",height=4,width=7)
plotEnrichment(david_allProteins_vs_differentialProteins_GO_MF_BP,
               name="david_allProteins_vs_differentialProteins_GO_MF_BP", 
               enrichment_cutoff = 1.3, Bonferroni_cutoff = 0.1, PDF=F)
dev.off()

pdf("david_allProteins_vs_multipleFeatures_GO_MF_BP.pdf",height=4,width=7)
plotEnrichment(david_allProteins_vs_multipleFeatures_GO_MF_BP,
                   name="david_allProteins_vs_multipleFeatures_GO_MF_BP", 
                   enrichment_cutoff = 1.3, Bonferroni_cutoff = 0.1, PDF=F)
dev.off()

###########################
###########################
###########################

proteoform_DiffExprProtein 

proteoform_DiffExprProteoform

n_proteoform_table <- unique(subset(proteoform_DiffExprProteoform, select=c("feature_id","apex","n_proteoforms_perProtein")))
proteoform_DiffExprProtein <- merge(proteoform_DiffExprProtein, n_proteoform_table, by=c("feature_id","apex"))

noProteoforms_DiffExprProtein <- subset(proteoform_DiffExprProtein, n_proteoforms_perProtein==1)
withProteoforms_DiffExprProtein <- subset(proteoform_DiffExprProtein, n_proteoforms_perProtein>=2)

getCount <- function(data){
  n_proteins <- length(unique(data$feature_id))
  n_differential_proteins <- length(unique(data[pBHadj<0.05][abs(medianLog2FC)>1]$feature_id))
  up <- unique(data[pBHadj < 0.05][medianLog2FC > 1]$feature_id)
  down <- unique(data[pBHadj < 0.05][medianLog2FC < -1]$feature_id)
  n_both <- length(intersect(up,down))
  n_up <- length(up[!up %in% down])
  n_down <- length(down[!down %in% up])
  n_unchanged <- n_proteins-n_differential_proteins
  data.table(
    level=c("proteins","proteins","proteins","proteins"),
    quant=c("unchanged","up & down","up","down"),
    count=c(n_unchanged, n_both, n_up, n_down)
  )
}

getCountFeatures <- function(data){
  data_m <- copy(data)
  data_m[, feature_id := paste0(c(feature_id,apex), collapse = "_"), by=c("feature_id","apex")]
  n_proteins <- length(unique(data_m$feature_id))
  n_differential_proteins <- length(unique(data_m[pBHadj<0.05][abs(medianLog2FC)>1]$feature_id))
  up <- unique(data_m[pBHadj < 0.05][medianLog2FC > 1]$feature_id)
  down <- unique(data_m[pBHadj < 0.05][medianLog2FC < -1]$feature_id)
  n_both <- length(intersect(up,down))
  n_up <- length(up[!up %in% down])
  n_down <- length(down[!down %in% up])
  n_unchanged <- n_proteins-n_differential_proteins
  data.table(
    level=c("features","features","features"),
    quant=c("unchanged","up","down"),
    count=c(n_unchanged, n_up, n_down)
  )
}

proteinStats_noProteoforms <- getCount(noProteoforms_DiffExprProtein)
proteinStats_noProteoforms[,proteoform := "1"]
proteinStats_withProteoforms <- getCount(withProteoforms_DiffExprProtein)
proteinStats_withProteoforms[,proteoform := "2"]
proteinFeatureStats_noProteoforms <- getCountFeatures(noProteoforms_DiffExprProtein)
proteinFeatureStats_noProteoforms[,proteoform := "1"]
proteinFeatureStats_withProteoforms <- getCountFeatures(withProteoforms_DiffExprProtein)
proteinFeatureStats_withProteoforms[,proteoform := "2"]

proteinFeatureDiffStats <- do.call("rbind", list(proteinStats_noProteoforms, 
                                                 proteinStats_withProteoforms, 
                                                 proteinFeatureStats_noProteoforms,
                                                 proteinFeatureStats_withProteoforms))

proteinFeatureDiffStats$level <- factor(proteinFeatureDiffStats$level, levels = c("proteins","features"))
proteinFeatureDiffStats$quant <- factor(proteinFeatureDiffStats$quant, levels = c("unchanged","up & down","up","down"))

pdf("proteinFeatureDiffStats.pdf", width=7, height=3.5)
ggplot(proteinFeatureDiffStats, aes(x=quant, y=count, fill=proteoform, group=level)) + 
  geom_col() + 
  facet_wrap(. ~ level, scales = "free") +
  theme_classic() +
  theme(axis.title.x=element_blank()) +
  scale_fill_manual(values=c("#999999", "#3288BD"), labels = c("1",expression(phantom(x) >= "2")), name="N proteoforms") +
  theme(legend.position = "bottom")
dev.off()

all_proteins_noProteoforms <- sum(proteinFeatureDiffStats[proteoform == 1][level=="proteins"]$count)
diff_proteins_noProteoforms <- sum(proteinFeatureDiffStats[proteoform == 1][level=="proteins"][quant != "unchanged"]$count)

all_proteins_noProteoforms
diff_proteins_noProteoforms
diff_proteins_noProteoforms/all_proteins_noProteoforms

all_proteins_withProteoforms <- sum(proteinFeatureDiffStats[proteoform == 2][level=="proteins"]$count)
diff_proteins_withProteoforms <- sum(proteinFeatureDiffStats[proteoform == 2][level=="proteins"][quant != "unchanged"]$count)

all_proteins_withProteoforms
diff_proteins_withProteoforms
diff_proteins_withProteoforms/all_proteins_withProteoforms

###########################
###########################
###########################

example_ids <- c("Q9Y6K5","P61978","Q9Y6E0")
example_subset <- subset(proteoform_DiffExprProtein, feature_id %in% example_ids)

###########################
###########################
###########################

#' ## Are proteins with multiple proteoforms more affected by PRPF8 depletion 
m <- all_proteins_withProteoforms      # Genes IN GO term == genes multiple proteoforms
n <- all_proteins_noProteoforms       # Genes NOT IN GO term == genes single proteoform
k <- diff_proteins_withProteoforms + diff_proteins_noProteoforms    # Gene hits, that is, differentially expressed == all genes that are diff exp.
x <- diff_proteins_withProteoforms  # Genes both IN GO term and differentially expressed 'hits' == genes with multiple proteoforms that are diff exp

probabilities <- dhyper(x, m, n, k, log = FALSE)
pvalue <- sum(probabilities)
pvalue

# Fisher
notDiff_noProteoform <- all_proteins_noProteoforms - diff_proteins_noProteoforms
diff_noProteoform <- diff_proteins_noProteoforms
notDiff_withProteoform <- all_proteins_withProteoforms-diff_proteins_withProteoforms
diff_withProteoform <- diff_proteins_withProteoforms

cMatrix <- matrix(c(notDiff_noProteoform,
         diff_noProteoform,
         notDiff_withProteoform,
         diff_withProteoform),
       nrow=2,
       dimnames=list(proteoform=c("noProteoform","Proteoform"),
                     differential=c("notDiff","Diff")))

fisher.test(cMatrix, alternative = "greater")


