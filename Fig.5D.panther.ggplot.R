
# ==========================================================
# Plot PANTHER GO-slim and complete biological processes
# quiescent age-DBsites genes
# for Fig. 7 Aging Cell submission
# December 2020
# ==========================================================

setwd("~/Dropbox/Desktop/ATAC-seq-PEAKS-AGE/Panther/Benjamini-Hochberg/")

library(ggplot2)

# ===================================
# GO-slim
# ===================================

# biol. processes in quiescent age-DBsite genes
df <- read.table("quiescent.age-DBsites.panther.GO-slim-biol.BH.corrected.txt", sep="\t", header=T)

# filter by padj < 0.05
df <- df[df$`bh.pvalue`<0.05, ]
### no significantly enriched GO-slim biological processes

# # filter by padj < 0.05
# df <- df[df$`bh.pvalue`<0.05, ]
# 
# # calculate -log10(padj)
# df$`-Log10(P-value)` <- -log10(df$bh.pvalue)
# 
# # change the order of factor levels for biological process
# df$biological_process <- factor(df$biological_process, levels= df$biological_process[order(df$`-Log10(P-value)`)])

# ===================================
# GO complete biological processes
# ===================================

# PANTHER complete biological processes in quiescent age-DBsite genes
df <- read.table("quiescent.age-DBsites.panther.GO-bio.BH.corrected.txt", sep="\t", header=T)

# filter by padj < 0.05
df <- df[df$`bh.pvalue`<0.05, ]

# calculate -log10(padj)
df$`-Log10(P-value)` <- -log10(df$bh.pvalue)

# change the order of factor levels for biological process
df$biological_process <- factor(df$biological_process, levels= df$biological_process[order(df$`-Log10(P-value)`)])

# take top biol. processes by padj
dat <- df[1:16, ]

# plot
p1 <- ggplot(dat, aes(`-Log10(P-value)`, biological_process)) +
    geom_point(aes(colour = gene_number, size = gene_number)) +
    scale_colour_gradient(low="navyblue", high="red") +
    scale_size_continuous(range = c(2,6)) +
    guides(color = guide_legend(), size = guide_legend()) +
    theme_bw() +
    theme(text = element_text(size=6),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank()) +
    xlim(5,8) +
    coord_fixed()
p1
