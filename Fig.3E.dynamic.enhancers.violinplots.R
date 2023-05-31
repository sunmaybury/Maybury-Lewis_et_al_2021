# ============================================================
# Violin plots
# Martynoga et al RNA-seq data
# Presence of active enhancers in distal dynamic chromatin
# AA and AQ separately done for aNSCs and qNSCs marks
# v2 9929 DBsites AvsQ normalized 
# October 2019
# ============================================================


library(dplyr)
library(ggplot2)

setwd("~/Dropbox/Brown_Webb/Desktop/ATAC-seq-PEAKS/AvsQ_analysis_2019_v2/Expression_distribution")

"%!in%" <- function(x,y) ! x %in% y

# Load RNA-seq data
expres <- read.csv("~/Dropbox/Desktop/RNA-seq-data/Martynoga_RNA-seq_complete.csv", stringsAsFactors = F, header = T)
# since log2FC is (quiescent/activated), recalculate and reorder columns
expres$log2FC_AQ <- log2(expres$FPKM_aNSCs / expres$FPKM_qNSCs)
expres <- expres[,c(2,7,6)]
colnames(expres) <- c("gene", "log2FC", "padj")

# separate differentially expressed genes
sig <- expres[expres$padj < 0.05, ]
nonsig <- setdiff(expres, sig)

# =============================================================

# Load distal dynamic sites with active enhancer marks (AA/AQ)
AAenhancers <- read.table("AA.distal.active.enhancers.genelist.txt", sep="\t", stringsAsFactors=F)
AQenhancers <- read.table("AQ.distal.active.enhancers.genelist.txt", sep="\t", stringsAsFactors=F)

# Clean up chromatin data to include DBsites and genes
AAenhancers <- AAenhancers[,c(4,7)]
AQenhancers <- AQenhancers[,c(4,7)]
colnames(AAenhancers) <- c("site", "gene")
colnames(AQenhancers) <- c("site", "gene")

# Add group columns
AAenhancers$group <- "1:AA"
AQenhancers$group <- "2:AQ"

# Merge expression and chromatin data
sig_AAenhancers <- merge(sig, AAenhancers, by="gene")
sig_AQenhancers <- merge(sig, AQenhancers, by="gene")
sig_AAenhancers[4] <- NULL
sig_AQenhancers[4] <- NULL

# Remove duplicate genes and expression
sig_AAenhancers <- sig_AAenhancers[!duplicated(sig_AAenhancers$gene), ]
sig_AQenhancers <- sig_AQenhancers[!duplicated(sig_AQenhancers$gene), ]

sig_dynamic_enhancers <- rbind(sig_AAenhancers, sig_AQenhancers)

# =============================================================

# All sigDE genes - proximal/distal dynamic genes = stable sigDE genes

# Load all genes with dynamic chromatin
all_dynamic <- read.table("AvsQ.DBsites.proximal+distal.genelist.txt", sep="\t", stringsAsFactors=F)
all_dynamic <- all_dynamic[, c(4,7)]
colnames(all_dynamic) <- c("site", "gene")

# Subtract dynamic chromatin genes from sig DE genes
sig_stable <- sig[sig$gene %!in% all_dynamic$gene, ]
sig_stable <- sig_stable[!duplicated(sig_stable$gene), ]   # 3423
sig_stable$group <- "0:stable"


# Violin plot
df_sig <- rbind(sig_stable, sig_AAenhancers, sig_AQenhancers)

p1 <- ggplot(df_sig, aes(x=group, y=log2FC)) +
  geom_violin(trim=F) +
  geom_boxplot(outlier.shape=NA, width=0.1)+
  theme_bw() 
p1

wilcox.test(sig_stable$log2FC, sig_AAenhancers$log2FC, paired=F)
# Wilcoxon rank sum test with continuity correction
# 
# data:  sig_stable$log2FC and sig_AAenhancers$log2FC
# W = 305577, p-value = 0.001545
# alternative hypothesis: true location shift is not equal to 0

wilcox.test(sig_stable$log2FC, sig_AQenhancers$log2FC, paired=F)
# Wilcoxon rank sum test with continuity correction
# 
# data:  sig_stable$log2FC and sig_AQenhancers$log2FC
# W = 197764, p-value = 7.637e-07
# alternative hypothesis: true location shift is not equal to 0



# Bar graph of mean log2FC

# calculate means (boxplot shows median)
means <- aggregate(log2FC~group, df_sig, mean)

p2 <- ggplot(data=means, aes(x=group, y=log2FC, fill=group)) +
  geom_bar(width=.8, stat="identity") +
  coord_cartesian(ylim=c(-1.5,1.5)) + 
  geom_hline(yintercept=0) + 
  theme_bw()

p2

t.test(sig_stable$log2FC, sig_AAenhancers$log2FC)
t.test(sig_stable$log2FC, sig_AQenhancers$log2FC)

summary(sig_stable$log2FC)
summary(sig_AAenhancers$log2FC)
summary(sig_AQenhancers$log2FC)
