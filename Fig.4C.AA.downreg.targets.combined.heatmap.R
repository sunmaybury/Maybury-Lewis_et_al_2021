# Filter iRegulon target genes for interested predictors
# AA significant downreg genes
# Building GRN for downreg genes in opening chromatin with NSC activation
# August 2020

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

# FOXO
target1 <- read.table("Foxo.TF.targets.txt", sep="\t", stringsAsFactors=F)    # 203
names(target1) <- "gene"
target1$regulator <- "FOXO"

# NFI
target2 <- read.table("NFI.TF.targets.txt", sep="\t", stringsAsFactors=F)    # 227
names(target2) <- "gene"
target2$regulator <- "NFI"

# SMAD1
target3 <- read.table("Smad1.TF.targets.txt", sep="\t", stringsAsFactors=F)    # 30
names(target3) <- "gene"
target3$regulator <- "SMAD1"

# ISL1
target4 <- read.table("Isl1.TF.targets.txt", sep="\t", stringsAsFactors=F)    # 84
names(target4) <- "gene"
target4$regulator <- "ISL1"

# ETV
target5 <- read.table("Etv.TF.targets.txt", sep="\t", stringsAsFactors=F)    # 227
names(target5) <- "gene"
target5$regulator <- "ETV"

# RFX
target6 <- read.table("Rfx.TF.targets.txt", sep="\t", stringsAsFactors=F)    # 41
names(target6) <- "gene"
target6$regulator <- "RFX"

# =======================================================
# Generate expression matrix for all target genes
# =======================================================

# Combine all target genes and regulators
targetall <- rbind(target1, target2, target3, target4, target5, target6)

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
targetgenes$FOXO <- "no"
targetgenes$NFI <- "no"
targetgenes$SMAD1 <- "no"
targetgenes$ISL1 <- "no"
targetgenes$ETV <- "no"
targetgenes$RFX <- "no"

# Change to yes if target of each regulator
targetgenes$FOXO[targetgenes$gene %in% target1$gene] <- "yes"
targetgenes$NFI[targetgenes$gene %in% target2$gene] <- "yes"
targetgenes$SMAD1[targetgenes$gene %in% target3$gene] <- "yes"
targetgenes$ISL1[targetgenes$gene %in% target4$gene] <- "yes"
targetgenes$ETV[targetgenes$gene %in% target5$gene] <- "yes"
targetgenes$RFX[targetgenes$gene %in% target6$gene] <- "yes"

# Annotation for each gene
ha1 <- targetgenes$FOXO
ha1_col <- c("yes"="blue4", "no"="white")

ha2 <- targetgenes$NFI
ha2_col <- c("yes"="blue1", "no"="white")

ha3 <- targetgenes$SMAD1
ha3_col <- c("yes"="dodgerblue2", "no"="white")

ha4 <- targetgenes$ISL1
ha4_col <- c("yes"="royalblue4", "no"="white")

ha5 <- targetgenes$ETV
ha5_col <- c("yes"="slateblue3", "no"="white")

ha6 <- targetgenes$RFX
ha6_col <- c("yes"="steelblue4", "no"="white")

# Heatmap
ht <- Heatmap(mat_scaled, column_order=1:5, column_km=2, name="Target gene expression", clustering_distance_rows="pearson",
              col=colorRamp2(c(-3, 0, 3), c("navyblue", "white", "darkred")),
              show_row_names=F, border=T) +
  Heatmap(ha1, name="FOXO", col=ha1_col, width=unit(0.5, "cm")) +
  Heatmap(ha2, name="NFI", col=ha2_col, width=unit(0.5, "cm")) +
  Heatmap(ha3, name="SMAD1", col=ha3_col, width=unit(0.5, "cm")) +
  Heatmap(ha4, name="ISL1", col=ha4_col, width=unit(0.5, "cm")) +
  Heatmap(ha5, name="ETV", col=ha5_col, width=unit(0.5, "cm")) +
  Heatmap(ha6, name="RFX", col=ha6_col, width=unit(0.5, "cm"))

pdf("iRegulon.AA.downreg.regulators.targets.heatmap.pdf")
draw(ht)
dev.off()