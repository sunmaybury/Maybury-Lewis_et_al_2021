# ==================================================================
# Violin plots
# Divide genes by chromatin (dynamic-AA/AQ, and the rest as stable)
# Look at differences in gene expression changes
# Martynoga et al RNA-seq data
# v2 9929 DBsites with -1/+1kb from TSS proximal gene assignments
# October 2019
# ==================================================================

setwd("~/Dropbox/Desktop/ATAC-seq-PEAKS/AvsQ_analysis_2019_v2/Expression_distribution")

library(dplyr)
library(ggplot2)

# Load chromatin data
chrom <- read.table("AvsQ.DBsites.proximal+distal.genelist.txt", sep="\t", stringsAsFactors=F)    # 4074
chrom <- chrom[, c(4,7,5)]
colnames(chrom) <- c("DBsite", "gene", "Fold")

# Load RNA-seq data
expres <- read.csv("~/Dropbox/Desktop/RNA-seq-data/Martynoga_RNA-seq_complete.csv", stringsAsFactors = F, header = T)
# since log2FC is (quiescent/activated), recalculate and reorder columns
expres$log2FC_AQ <- log2(expres$FPKM_aNSCs / expres$FPKM_qNSCs)
expres <- expres[,c(2,7,6)]
colnames(expres) <- c("gene", "log2FC", "padj")

# Merge chromatin and RNA-seq data
dynamic_genes <- merge(expres, chrom, by="gene")
# Filter by significant DE
dynamic_genes_simple <- dynamic_genes[,c(1,2,3)]
stable_genes <- setdiff(expres, dynamic_genes_simple)
stable_genes <- stable_genes[!duplicated(stable_genes$gene), ]

# Skipped dividng stable genes into open and closed chromatin for now

# Divide up dynamic genes

# Dynamic AA and AQ
dynamic_AA <- dynamic_genes[dynamic_genes$Fold > 0, ]
dynamic_AQ <- dynamic_genes[dynamic_genes$Fold < 0, ]

dynamic_AA <- dynamic_AA[,c(1,2,3)]
dynamic_AQ <- dynamic_AQ[,c(1,2,3)]

dynamic_AA <- dynamic_AA[!duplicated(dynamic_AA$gene), ]
dynamic_AQ <- dynamic_AQ[!duplicated(dynamic_AQ$gene), ]

# ========================================================================

# Plotting

# Add groups for plot
stable_genes$group <- "1:stable"
dynamic_AA$group <- "2:dynamic-AA"
dynamic_AQ$group <- "3:dynamic-AQ"

# Make plots
df_all <- rbind(stable_genes, dynamic_AA, dynamic_AQ)

p1 <- ggplot(df_all, aes(x=group, y=log2FC)) +
  geom_violin(trim=F) +
  geom_boxplot(outlier.shape=NA, width=0.1) +
  theme_bw()
p1

# ========================================================================

# Statistical testing

wilcox.test(stable_genes$log2FC, dynamic_AA$log2FC, paired=F)
# Wilcoxon rank sum test with continuity correction
# 
# data:  stable_genes$log2FC and dynamic_AA$log2FC
# W = 1128342, p-value = 1.777e-06
# alternative hypothesis: true location shift is not equal to 0

result2 <- wilcox.test(stable_genes$log2FC, dynamic_AA$log2FC, paired=F)
result2$p.value
# [1] 1.777482e-06

wilcox.test(stable_genes$log2FC, dynamic_AQ$log2FC, paired=F)
# Wilcoxon rank sum test with continuity correction
# 
# data:  stable_genes$log2FC and dynamic_AQ$log2FC
# W = 442811, p-value = 5.945e-09
# alternative hypothesis: true location shift is not equal to 0

result2 <- wilcox.test(stable_genes$log2FC, dynamic_AQ$log2FC, paired=F)
result2$p.value
# [1] 5.945195e-09
