# =========================================================================
# Scatterplot of chromatin accessibility (ATAC-seq) 
# and associated gene expression (Leeman et al. RNA-seq)
# v2 9929 DBsites AvsQ DiffBind, -1/+1kb TSS gene assignment
# October 2019
# =========================================================================

setwd("~/Dropbox/Brown_Webb/Desktop/ATAC-seq-PEAKS/AvsQ_analysis_2019_v2/Expression_correlation")

library(dplyr)
library(ggplot2)

# load Leeman RNA-seq data
expres <- read.csv("~/Dropbox/Brown_Webb/Desktop/RNA-seq-data_pub/DEseq_young_aNSC-qNSC_Leeman_et_al.csv", sep=",", stringsAsFactors=F)
expres <- expres[,c(1,3,7)]
colnames(expres) <- c("gene", "log2FC", "padj")

# separate differentially expressed genes
sig_expres <- expres[expres$padj < 0.05, ]
nonsig_expres <- setdiff(expres, sig_expres)

# load chromatin data
chrom <- read.table("AvsQ.DBsites.proximal+distal.genelist.txt", sep="\t", stringsAsFactors=F)

# clean up chromatin data to include DBsite, gene, and log2FC(A/Q)
chrom <- chrom[,c(4,7,5)]
colnames(chrom) <- c("DBsite", "gene", "Fold")

# =========================================================================

# Merge chromatin and expression data
df_sig <- merge(sig_expres, chrom, by="gene")

# Make scatterplot
p <- ggplot(data=df_sig, aes(x=Fold, y=log2FC))
p +
  geom_rect(aes(xmin=-Inf,xmax=Inf,ymin=0,ymax=Inf),alpha=.05,fill="plum2") +
  geom_rect(aes(xmin=-Inf,xmax=Inf,ymin=-Inf,ymax=0),alpha=.05,fill="skyblue1") +
  geom_point() +
  xlab("log2FC (Chromatin accessibility)") +
  ylab("log2FC (Expression)") +
  theme_minimal() +
  theme(panel.grid.major=element_line(color="black", size=0.5)) +
  geom_hline(aes(yintercept=0)) +
  geom_vline(aes(xintercept=0))

res1 <- cor.test(df_sig$Fold, df_sig$log2FC, method="spearman")
res1
# Spearman's rank correlation rho
# 
# data:  df_sig$Fold and df_sig$log2FC
# S = 229500533, p-value = 0.0007376
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
# rho 
# 0.09930482 

# =========================================================================

# calculate number of points in each quadrant
sig_opened <- df_sig[df_sig$Fold > 0, ]
sig_closed <- df_sig[df_sig$Fold < 0, ]

opened_upreg <- sig_opened[sig_opened$log2FC > 0, ]     # 480 points
opened_downreg <- sig_opened[sig_opened$log2FC < 0, ]   # 441 points

closed_upreg <- sig_closed[sig_closed$log2FC > 0, ]     # 72 points
closed_downreg <- sig_closed[sig_closed$log2FC < 0, ]   # 159 points
