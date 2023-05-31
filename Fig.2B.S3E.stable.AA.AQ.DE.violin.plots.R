# ==================================================================
# Violin plots
# Divide genes by chromatin (dynamic-AA/AQ, and the rest as stable)
# Look at differences in gene expression changes
# Leeman et al
# v2 9929 DBsites with -1/+1kb from TSS proximal gene assignments
# October 2019
# ==================================================================

setwd("~/Dropbox/Brown_Webb/Desktop/ATAC-seq-PEAKS/AvsQ_analysis_2019_v2/Expression_distribution_v2")

library(dplyr)
library(ggplot2)

# Load chromatin data
chrom <- read.table("AvsQ.DBsites.proximal+distal.genelist.txt", sep="\t", stringsAsFactors=F)    # 4074
chrom <- chrom[, c(4,7,5)]
colnames(chrom) <- c("DBsite", "gene", "Fold")

# Load Leeman RNA-seq data and filter p<0.05 DE genes
expres <- read.csv("~/Dropbox/Desktop/RNA-seq-data/DEseq_young_aNSC-qNSC_Leeman_et_al.csv", sep=",", stringsAsFactors=F)
expres <- expres[,c(1,3,7)]
expres <- expres[expres$padj < 0.05, ] 
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

# Write up files of differentially expressed genes with AA and AQ chromatin
# AA_upreg <- dynamic_AA[dynamic_AA$log2FC > 0, ][1]
# AA_downreg <- dynamic_AA[dynamic_AA$log2FC < 0, ][1]
# AA_DE <- dynamic_AA[1]
# AQ_DE <- dynamic_AQ[1]

# write.table(AA_DE, file="AA.sig.DE.genes.txt", sep="\t", quote=F, row.names=F, col.names=F)
# write.table(AA_upreg, file="AA.sig.upreg.genes.txt", sep="\t", quote=F, row.names=F, col.names=F)
# write.table(AA_downreg, file="AA.sig.downreg.genes.txt", sep="\t", quote=F, row.names=F, col.names=F)
# write.table(AQ_DE, file="AQ.sig.DE.genes.txt", sep="\t", quote=F, row.names=F, col.names=F)

# ========================================================================

# Plotting #1

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

# Statistical testing

wilcox.test(stable_genes$log2FC, dynamic_AA$log2FC, paired=F)
# Wilcoxon rank sum test with continuity correction
# 
# data:  stable_genes$log2FC and dynamic_AA$log2FC
# W = 1143548, p-value = 0.003401
# alternative hypothesis: true location shift is not equal to 0

result2 <- wilcox.test(stable_genes$log2FC, dynamic_AA$log2FC, paired=F)
result2$p.value
# [1] 0.003401073

wilcox.test(stable_genes$log2FC, dynamic_AQ$log2FC, paired=F)
# Wilcoxon rank sum test with continuity correction
# 
# data:  stable_genes$log2FC and dynamic_AQ$log2FC
# W = 410079, p-value = 8.228e-05
# alternative hypothesis: true location shift is not equal to 0

result2 <- wilcox.test(stable_genes$log2FC, dynamic_AQ$log2FC, paired=F)
result2$p.value
# [1] 8.227501e-05

# ========================================================================

# January 2019
# Adding AR (React overlaps with AA and AQ to violin plot)
# Plotting #2

# # Read data
# react <- read.table("react.AA.AQ.olap.genes.txt", sep="\t", stringsAsFactors=F)
# react[1] <- NULL
# names(react) <- "gene"
# 
# react <- unique(react)

# # Append Leeman DE data
# dynamic_AR <- inner_join(react, expres)

# # Add group
# dynamic_AR$group <- "4:dynamic-AR"
# 
# df_all2 <- rbind(df_all, dynamic_AR)

# p2 <- ggplot(df_all2, aes(x=group, y=log2FC)) +
#   geom_violin(trim=F) +
#   geom_boxplot(outlier.shape=NA, width=0.1) +
#   theme_bw()
# p2

# Statistical testing

# wilcox.test(stable_genes$log2FC, dynamic_AR$log2FC, paired=F)
# Wilcoxon rank sum test with continuity correction
# 
# data:  stable_genes$log2FC and dynamic_AR$log2FC
# W = 1047661, p-value = 0.9553
# alternative hypothesis: true location shift is not equal to 0

# wilcox.test(dynamic_AA$log2FC, dynamic_AR$log2FC, paired=F)
# Wilcoxon rank sum test with continuity correction
# 
# data:  dynamic_AA$log2FC and dynamic_AR$log2FC
# W = 234516, p-value = 0.03871
# alternative hypothesis: true location shift is not equal to 0

# wilcox.test(dynamic_AQ$log2FC, dynamic_AR$log2FC, paired=F)
# Wilcoxon rank sum test with continuity correction
# 
# data:  dynamic_AQ$log2FC and dynamic_AR$log2FC
# W = 53006, p-value = 0.0005619
# alternative hypothesis: true location shift is not equal to 0

# ======== Separate AR to AR:AA and AR:AQ

# Read data
reactAA <- read.table("react.AA.olap.genes.txt", sep="\t", stringsAsFactors=F)
names(reactAA) <- "gene"
reactAA <- unique(reactAA)

reactAQ <- read.table("react.AQ.olap.genes.txt", sep="\t", stringsAsFactors=F)
names(reactAQ) <- "gene"
reactAQ <- unique(reactAQ)

# Append Leeman DE data
AR_AA <- inner_join(reactAA, expres)
AR_AQ <- inner_join(reactAQ, expres)

# Add group
AR_AA$group <- "4:AR:AA"
AR_AQ$group <- "5:AR:AQ"

df_all2 <- rbind(df_all, AR_AA, AR_AQ)

# ====================================================================================

# Plotting #3

# unique <- read.table("React.unique.GREAT.genes.txt", sep="\t", stringsAsFactors=F)
# names(unique) <- "gene"
# 
# unique <- unique(unique)
# 
# reactOnly <- inner_join(unique, expres)
# 
# test <- inner_join(reactOnly, df_all)
# 
# reactOnly$group <- "6:AR:unique"

p2 <- ggplot(df_all2, aes(x=group, y=log2FC)) +
  geom_violin(trim=F, width=1.15) +
  geom_boxplot(outlier.shape=NA, width=0.1) +
  theme_bw()
p2


# Statistical testing

wilcox.test(stable_genes$log2FC, AR_AA$log2FC, paired=F)
# Wilcoxon rank sum test with continuity correction
# 
# data:  stable_genes$log2FC and AR_AA$log2FC
# W = 1101832, p-value = 0.00261
# # alternative hypothesis: true location shift is not equal to 0


wilcox.test(stable_genes$log2FC, AR_AQ$log2FC, paired=F)
# Wilcoxon rank sum test with continuity correction
# 
# data:  stable_genes$log2FC and AR_AQ$log2FC
# W = 238120, p-value = 2.197e-06
# alternative hypothesis: true location shift is not equal to 0




