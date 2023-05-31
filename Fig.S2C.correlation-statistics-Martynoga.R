# =========================================================================
# Scatterplot with chromatin accessibility(ATAC-seq)
# and associated gene expression (Martynoga et al. RNA-seq)
# v2 9929 DBsites normalized AvsQ
# October 2019
# =========================================================================


setwd("~/Dropbox/Desktop/ATAC-seq-PEAKS/AvsQ_analysis_2019_v2/Expression_correlation")

library(dplyr)
library(ggplot2)

# Load RNA-seq data
martynoga <- read.table("~/Dropbox/Desktop/RNA-seq-data/Martynoga_RNA-seq_complete.txt", sep = "\t", stringsAsFactors = FALSE, header = TRUE)

# Remove Ensembl ID
martynoga[,1] <- NULL

# Re-calculate log2FC (Activated/Quiescent)
colnames(martynoga) <- c("gene", "FPKM_act", "FPKM_quies", "log2FC_Q/A", "padj")
martynoga$log2FC <- log2(martynoga$FPKM_act / martynoga$FPKM_quies)

# Clean up 
martynoga[,4] <- NULL
martynoga <- martynoga[,c(1,5,4)]

# Load chromatin data
chrom <- read.table("AvsQ.DBsites.proximal+distal.genelist.txt", sep="\t", stringsAsFactors=F)

# clean up chromatin data to include DBsite, gene, and log2FC(A/Q)
chrom <- chrom[,c(4,7,5)]
colnames(chrom) <- c("DBsite", "gene", "Fold")

# =========================================================================

# Merge chromatin and expression data
df_sig <- merge(martynoga, chrom, by="gene")

# Make scatter plot
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
# S = 203801260, p-value = 2.587e-10
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
# rho 
# 0.1854045 

# =========================================================================

# calculate number of points in each quadrant

sig_opened <- df_sig[df_sig$Fold > 0, ]
sig_closed <- df_sig[df_sig$Fold < 0, ]

opened_upreg <- sig_opened[sig_opened$log2FC > 0, ]     # 507 points
opened_downreg <- sig_opened[sig_opened$log2FC < 0, ]   # 414 points

closed_upreg <- sig_closed[sig_closed$log2FC > 0, ]     # 53 points
closed_downreg <- sig_closed[sig_closed$log2FC < 0, ]   # 171 points
