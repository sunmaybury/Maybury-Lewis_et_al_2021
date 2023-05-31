# =============================================================================
# AvsQ proximal stable genes that are upregulated by pvalue (DE)
# Leeman et al RNA-seq data
# Took top quartile of upregulated genes
# Background is all truStable proximal genes
# Then each bed file was adjusted to overlaps with -50/+50bp window around TSS
# Based on Tammy Gershon's feedback
# All files were converted to fa and run on ElemeNT
# This script is to determine whether there are enriched elements
# v2 DiffBind (9929 DBsites) AvsQ
# January 2019
# =============================================================================

library(dplyr)

setwd("~/Dropbox/Brown_Webb/Desktop/ATAC-seq-PEAKS/AvsQ_analysis_2019_v2/Core_promoter/expression_v2_50bp")

# Load ElemeNT output from total population (all genes with stable promoters)
population <- read.table("AvsQ.truStable.proximal.TSS100bp.ElemeNT.output.txt",
                         sep="\t", stringsAsFactors=F, header=T)    #11514

# Get numbers of core promoter elements in the population
TATA <- population[!grepl("no", population$TATA), ]    #285
hInr <- population[!grepl("no", population$hInr), ]    #11368
dInr <- population[!grepl("no", population$dInr), ]    #6735
DPE <- population[!grepl("no", population$DPE), ]    #3562
MTE <- population[!grepl("no", population$MTE), ]    #63
hTCT <- population[!grepl("no", population$hTCT), ]    #2967
dTCT <- population[!grepl("no", population$dTCT), ]    #354

TATA <- TATA[1]
write.table(TATA, file="stable.TSS100bp.TATA.txt", sep="\t", quote=F, row.names=F, col.names=F) 
# This is included as background in STRING analysis All genes with stable open promoters containing TATA

# ============================================
# AvsQ truStable upregulated top quartile (1)
# ============================================

sample <- read.table("AvsQ.truStable.proximal.top.quartile.upreg.TSS100bp.ElemeNT.output.txt",
                     sep="\t", stringsAsFactors=F, header=T)  #475

# sample2 <- sample[, c(1,3,10)]
# write.table(sample2, "Avsq.truStable.proximal.top.quartile.upreg.TSS100bp.TATA.TCT.output.txt",
#             sep="\t", quote=F, row.names=F, col.names=T)

# Get numbers of core promoter elements in the sample
TATA <- sample[!grepl("no", sample$TATA), ]    #30
hInr <- sample[!grepl("no", sample$hInr), ]    #474
dInr <- sample[!grepl("no", sample$dInr), ]    #282
DPE <- sample[!grepl("no", sample$DPE), ]    #132
MTE <- sample[!grepl("no", sample$MTE), ]    #2
hTCT <- sample[!grepl("no", sample$hTCT), ]    #162
dTCT <- sample[!grepl("no", sample$dTCT), ]    #42

# Test enrichment
population <- 11514
true_pop <- c(285,11368,3562,63,2967)
sample <- 475
true_sam <- c(30,474,132,2,162)
p.val <- phyper(true_sam - 1, true_pop, population-true_pop, sample, lower.tail=F)
p.adjust(p.val, method="bonferroni")

# TATA, hInr, DPE, MTE, hTCT
# [1] 0.0000114077 0.0754359774 1.0000000000 1.0000000000 0.0001137164

TATA <- TATA[1]
write.table(TATA, file="top.quartile.upreg.stable.TSS100bp.TATA.txt", sep="\t", quote=F, row.names=F, col.names=F)

TCT <- hTCT[1]
write.table(TCT, file="top.quartile.upreg.stable.TSS100bp.hTCT.txt", sep="\t", quote=F, row.names=F, col.names=F)

# ====================================================
# AvsQ truStable upregulated bottom quartile (2+3+4)
# ====================================================

sample <- read.table("AvsQ.truStable.proximal.bottom.quartile.upreg.TSS100bp.ElemeNT.output.txt",
                     sep="\t", stringsAsFactors=F, header=T)    #1431

# Get numbers of core promoter elements in the sample
TATA <- sample[!grepl("no", sample$TATA), ]    #42
hInr <- sample[!grepl("no", sample$hInr), ]    #1424
dInr <- sample[!grepl("no", sample$dInr), ]    #879
DPE <- sample[!grepl("no", sample$DPE), ]    #434
MTE <- sample[!grepl("no", sample$MTE), ]    #12
hTCT <- sample[!grepl("no", sample$hTCT), ]    #375
dTCT <- sample[!grepl("no", sample$dTCT), ]    #42

# Test enrichment
population <- 11514
true_pop <- c(285,11368,3562,63,2967)
sample <- 1431
true_sam <- c(42,1424,434,12,375)
p.val <- phyper(true_sam - 1, true_pop, population-true_pop, sample, lower.tail=F)
p.adjust(p.val, method="bonferroni")

# TATA, hInr, DPE, MTE, hTCT
# [1] 0.676601658 0.007653916 1.000000000 0.426068851 1.000000000

TATA <- TATA[1]
write.table(TATA, file="bottom.quartile.upreg.stable.TSS100bp.TATA.txt", sep="\t", quote=F, row.names=F, col.names=F)

TCT <- hTCT[1]
write.table(TCT, file="bottom.quartile.upreg.stable.TSS100bp.hTCT.txt", sep="\t", quote=F, row.names=F, col.names=F)


# ====================================================
# AvsQ truStable upregulated second quartile (2)
# ====================================================

sample <- read.table("AvsQ.truStable.proximal.second.quartile.upreg.TSS100bp.ElemeNT.output.txt",
                     sep="\t", stringsAsFactors=F, header=T)    #458

# Get numbers of core promoter elements in the sample
TATA <- sample[!grepl("no", sample$TATA), ]    #19
hInr <- sample[!grepl("no", sample$hInr), ]    #457
dInr <- sample[!grepl("no", sample$dInr), ]    #290
DPE <- sample[!grepl("no", sample$DPE), ]    #142
MTE <- sample[!grepl("no", sample$MTE), ]    #4
hTCT <- sample[!grepl("no", sample$hTCT), ]    #130
dTCT <- sample[!grepl("no", sample$dTCT), ]    #16

# Test enrichment
population <- 11514
true_pop <- c(285,11368,3562,63,2967)
sample <- 458
true_sam <- c(19,457,142,4,130)
p.val <- phyper(true_sam - 1, true_pop, population-true_pop, sample, lower.tail=F)
p.adjust(p.val, method="bonferroni")

# TATA, hInr, DPE, MTE, hTCT
# [1] 0.09701444 0.09153105 1.00000000 1.00000000 0.53064800

# ====================================================
# AvsQ truStable upregulated second quartile (3)
# ====================================================

sample <- read.table("AvsQ.truStable.proximal.third.quartile.upreg.TSS100bp.ElemeNT.output.txt",
                     sep="\t", stringsAsFactors=F, header=T)    #489

# Get numbers of core promoter elements in the sample
TATA <- sample[!grepl("no", sample$TATA), ]    #16
hInr <- sample[!grepl("no", sample$hInr), ]    #487
dInr <- sample[!grepl("no", sample$dInr), ]    #293
DPE <- sample[!grepl("no", sample$DPE), ]    #134
MTE <- sample[!grepl("no", sample$MTE), ]    #4
hTCT <- sample[!grepl("no", sample$hTCT), ]    #109
dTCT <- sample[!grepl("no", sample$dTCT), ]    #16

# Test enrichment
population <- 11514
true_pop <- c(285,11368,3562,63,2967)
sample <- 489
true_sam <- c(16,487,134,4,109)
p.val <- phyper(true_sam - 1, true_pop, population-true_pop, sample, lower.tail=F)
p.adjust(p.val, method="bonferroni")

# TATA, hInr, DPE, MTE, hTCT
# [1] 0.7788032 0.2455996 1.0000000 1.0000000 1.0000000


