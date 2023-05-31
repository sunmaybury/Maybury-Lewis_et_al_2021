# Filter iRegulon target genes for interested predictors
# AA significant upreg genes
# Building GRN for upreg genes in opening chromatin with NSC activation
# Combine all transcriptional regulators
# September 2020

# Target genes of predicted regulators from iRegulon results panel, Transcription Factors, saved as txt files

setwd("~/Dropbox/Desktop/ATAC-seq-PEAKS/AvsQ_analysis_2019_v2/iRegulon")

library(tidyr)
library(dplyr)

# Load ComplexHeatmap package
library(ComplexHeatmap)
library(circlize)
library(tidyverse)

# ============================================================================
# Load normalized expression values and corresponding differential expression
# ============================================================================

# Load VST values
vst <- read.csv("~/Dropbox/Desktop/RNA-Seq-data/Leeman_normExpressionValues.csv", stringsAsFactors = F)
vst <- vst[,c(1,11,12,5,6,7)]
colnames(vst) <- c("gene", "qNSCrep1", "qNSCrep3", "aNSCrep1", "aNSCrep2", "aNSCrep3")

# Load differential expression
expr <- read.csv("~/Dropbox/Desktop/RNA-Seq-data/DEseq_young_aNSC-qNSC_Leeman_et_al.csv", stringsAsFactors = F)[,c(1,3,7)]
colnames(expr) <- c("gene", "log2FC", "padj")

# Filter by padj
expr <- expr[expr$padj < 0.05, ]

# Get VST values of sig. diff genes 
de <- merge(vst, expr, by="gene")      # 1 extra because of 2-March gene name change from excel, keep both.

# Convert VST values to matrix with genes as row names
count.mat <- as.matrix(de)[,c(1,2,3,4,5,6)]
rownames(count.mat) <- de[,1]
count.mat <- count.mat[,-1]

# ============================================================================
# Load iRegulon predicted regulator/targets results
# ============================================================================

# E2F1
target1 <- read.table("E2f1.TF.targets.txt", sep="\t", stringsAsFactors=F)    # 241
names(target1) <- "gene"
target1$regulator <- "E2F1"

# Sox
target2 <- read.table("Sox.TF.targets.txt", sep="\t", stringsAsFactors=F)    # 79
names(target2) <- "gene"
target2$regulator <- "SOX"

# KLF
target3 <- read.table("Klf.TF.targets.txt", sep="\t", stringsAsFactors=F)    # 138
names(target3) <- "gene"
target3$regulator <- "KLF"

# NFYC
target4 <- read.table("Nfyc.TF.targets.txt", sep="\t", stringsAsFactors=F)    # 168
names(target4) <- "gene"
target4$regulator <- "NFYC"

# POU3F1
target5 <- read.table("Pou3f1.TF.targets.txt", sep="\t", stringsAsFactors=F)    # 53
names(target5) <- "gene"
target5$regulator <- "POU3F1"

# ==============================================
# Look for shared targets among regulators
# ==============================================

# E2f_Sox<- as.data.frame(intersect(target1$gene, target2$gene))    # 75
# E2f_Klf <- as.data.frame(intersect(target1$gene, target3$gene))    # 128
# E2f_Nfyc <- as.data.frame(intersect(target1$gene, target4$gene))    # 154
# E2f_Pou3f1 <- as.data.frame(intersect(target1$gene, target5$gene))    # 50

# Since majority of the targets are shared, make a dataframe of all targets with duplicates removed
# Then add on regulator bars to the right of the heatmap

# ===================================================
# Generate expression matrix for all target genes
# ===================================================

# Combine all target genes and regulators
targetall <- rbind(target1, target2, target3, target4, target5)

# Get a duplicate-removed list of target genes for expression matrix
targetgenes <- as.data.frame(targetall$gene)
targetgenes <- unique(targetgenes)
names(targetgenes) <- "gene"

# Convert expression values from TF target genes into numeric matrix
# Merge list of all targets by de
dfall <- merge(targetgenes, de, by="gene")

# As matrix
mat <- as.matrix(dfall[, grep("NSC", colnames(dfall))])
rownames(mat) <- dfall[,1]

# Scale matrix
mat_scaled <- t(apply(mat, 1, scale))
colnames(mat_scaled) <- c("qNSCrep1", "qNSCrep3", "aNSCrep1", "aNSCrep2", "aNSCrep3")

# ===========================
# Annotations - regulators
# ===========================

# Make data frame of target genes with yes/no for each regulator
targetgenes$E2F1 <- "no"
targetgenes$SOX <- "no"
targetgenes$KLF <- "no"
targetgenes$NFYC <- "no"
targetgenes$POU3F1 <- "no"

# Change to yes if target of each regulator
targetgenes$E2F1[targetgenes$gene %in% target1$gene] <- "yes"
targetgenes$SOX[targetgenes$gene %in% target2$gene] <- "yes"
targetgenes$KLF[targetgenes$gene %in% target3$gene] <- "yes"
targetgenes$NFYC[targetgenes$gene %in% target4$gene] <- "yes"
targetgenes$POU3F1[targetgenes$gene %in% target5$gene] <- "yes"

# Annotation for each gene
ha1 <- targetgenes$E2F1
ha1_col <- c("yes"="hotpink3", "no"="white")

ha2 <- targetgenes$SOX
ha2_col <- c("yes"="hotpink4", "no"="white")

ha3 <- targetgenes$KLF
ha3_col <- c("yes"="indianred", "no"="white")

ha4 <- targetgenes$NFYC
ha4_col <- c("yes"="indianred1", "no"="white")

ha5 <- targetgenes$POU3F1
ha5_col <- c("yes"="lightpink3", "no"="white")


# Heatmap
ht <- Heatmap(mat_scaled, column_order=1:5, column_km=2, name="Target gene expression", clustering_distance_rows="pearson",
               col=colorRamp2(c(-3, 0, 3), c("navyblue", "white", "darkred")),
               show_row_names=F, border=T) +
      Heatmap(ha1, name="E2F1", col=ha1_col, width=unit(0.5, "cm")) +
      Heatmap(ha2, name="SOX", col=ha2_col, width=unit(0.5, "cm")) +
      Heatmap(ha3, name="KLF", col=ha3_col, width=unit(0.5, "cm")) +
      Heatmap(ha4, name="NFYC", col=ha4_col, width=unit(0.5, "cm")) +
      Heatmap(ha5, name="POU3F1", col=ha5_col, width=unit(0.5, "cm"))

pdf("iRegulon.AA.upreg.regulators.targets.heatmap.pdf")
draw(ht)
dev.off()

