
# ==============================================
# Making figures for DiffBind young vs. aged
# Quiescent and activated NSCs
# Code tested from markdown
# 2 replicates each, with batch correction
# Fig 5 for Aging Cell submission
# December 2020
# ==============================================

setwd("~/Dropbox/Brown_Webb/Desktop/ATAC-seq-PEAKS-AGE/Diffbind/files")
library(DiffBind)

# ============================
# Quiescent
# ============================

# load data
all <- dba(sampleSheet="aged.ATAC.quiescent.2.csv")

# count reads in bam files, and establish contrasts
all <- dba.count(all, minOverlap=2)
all <- dba.contrast(all, categories=DBA_TISSUE, block=DBA_REPLICATE, minMembers=2)

# differential analysis
all <- dba.analyze(all)

# PCA with batch correction
# pdf("diffbind.aged.quiescent.PCA.pdf")
# dba.plotPCA(all, contrast=1, method=DBA_DESEQ2_BLOCK, label=DBA_REPLICATE)
# dev.off()

# plot Pearson correlation heatmap
values <- dba.plotHeatmap(all)    ## Fig 5A

# =====================================================
# get consensus peaks between quiescent ages
all <- dba(sampleSheet="aged.ATAC.quiescent.2.csv")
Qpeaks <- dba(all, mask=all$masks$Quiescent, minOverlap=4)
Qsites <- dba.peakset(Qpeaks, bRetrieve=T, DataType=DBA_DATA_FRAME)[,c(1:3)]
write.table(Qsites, file="Quies.young.old.consPeaks.bed", sep="\t", quote=F, row.names=F, col.names=F)
# =====================================================


# ============================
# Activated
# ============================

# load data
all <- dba(sampleSheet="aged.ATAC.activated.2.csv")

# count reads in bam files, and establish contrasts
all <- dba.count(all, minOverlap=2)
all <- dba.contrast(all, categories=DBA_TISSUE, block=DBA_REPLICATE, minMembers=2)

# differential analysis
all <- dba.analyze(all)

# PCA with batch correction
# pdf("diffbind.aged.activated.PCA.pdf")
# dba.plotPCA(all, contrast=1, method=DBA_DESEQ2_BLOCK, label=DBA_REPLICATE)
# dev.off()

# plot Pearson correlation heatmap
values <- dba.plotHeatmap(all)     # Fig 5A


