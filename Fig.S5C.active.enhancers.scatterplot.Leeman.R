# ======================================================
# Scatterplot
# Leeman et al RNA-seq FDR<0.05 DE genes
# Presence of active enhancers in dynamic chromatin
# AA and AQ separately done for aNSCs and qNSCs marks
# v2 9929 DBsites AvsQ normalized
# October 2019
# ======================================================

library(dplyr)
library(ggplot2)

setwd("~/Dropbox/Desktop/ATAC-seq-PEAKS/AvsQ_analysis_2019_v2/Expression_distribution_v2")

# Load chromatin data
AA <- read.table("AA.distal.active.enhancers.genelist.txt", sep="\t", stringsAsFactors=F)
AA <- AA[,c(4,5,7)]
colnames(AA) <- c("site", "fold", "gene")

AQ <- read.table("AQ.distal.active.enhancers.genelist.txt", sep="\t", stringsAsFactors=F)
AQ <- AQ[,c(4,5,7)]
colnames(AQ) <- c("site", "fold", "gene")

# =======================================================================

# Load RNA-seq data
expres <- read.csv("~/Dropbox/Desktop/RNA-seq-data/DEseq_young_aNSC-qNSC_Leeman_et_al.csv", sep=",", stringsAsFactors=F)
expres <- expres[,c(1,3,7)]
colnames(expres) <- c("gene", "log2FC", "padj")

# separate differentially expressed genes
sig <- expres[expres$padj < 0.05, ]
nonsig <- setdiff(expres, sig)

# =======================================================================

# Merge chromatin and RNA-seq ata
AA_sig <- inner_join(sig, AA, by="gene")    #249 opened enhancers associated with sig DE genes
AQ_sig <- inner_join(sig, AQ, by="gene")    #77 closed enhancers associated with sig DE genes

# Combine opened and closed enhancers associated with DE genes for scatterplot
df_sig <- rbind(AA_sig, AQ_sig)

# Plot
p <- ggplot(data=df_sig, aes(x=fold, y=log2FC))
p +
  geom_rect(aes(xmin=-Inf,xmax=Inf,ymin=0,ymax=Inf),alpha=.05,fill="plum2") +
  geom_rect(aes(xmin=-Inf,xmax=Inf,ymin=-Inf,ymax=0),alpha=.05,fill="skyblue1") +
  geom_point() +
  xlab("log2FC (Chromatin accessibility)") +
  ylab("log2FC (Expression)") +
  xlim(-4,4) +
  ylim(-25,25) +
  theme_minimal() +
  theme(panel.grid.major=element_line(color="black", size=0.5)) +
  geom_hline(aes(yintercept=0)) +
  geom_vline(aes(xintercept=0))

# =======================================================================

# Count points in each quadrant
AA_upreg <- AA_sig[AA_sig$log2FC>0, ]    # 110
AA_downreg <- AA_sig[AA_sig$log2FC<0, ]  # 139

AQ_upreg <- AQ_sig[AQ_sig$log2FC>0, ]     # 21
AQ_downreg <- AQ_sig[AQ_sig$log2FC<0, ]   # 56
