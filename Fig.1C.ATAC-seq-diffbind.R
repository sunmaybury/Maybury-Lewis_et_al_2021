

setwd("~/Dropbox/Brown_Webb/Desktop/ATAC-seq-PEAKS/DiffBind")
library(DiffBind)
library(rtracklayer)
# library(rgl)

# =============================================================
## Comparison across quiescent, activated
## First compare quiescent and activated
## Then compare activated and reactivated
## Reactivated downsampled to 40 million reads
## October 25 2019
## R version 3.3.1
# =============================================================

### Part I - Activated vs Quiescent

# ===================================================================================
# Save data with read counts as DBA object in R workspace

# Load data
# all <- dba(sampleSheet="AvsQ.postnatal.analysis.csv")

# Save peaks as DB object
# savefile <- dba.save(all, 'AvsQPeaks')

# Count reads
# all <- dba.count(all, minOverlap=2)

# Save read count DB object
# savefile <- dba.save(all, 'AvsQReads')

# ===================================================================================

# Load data
all <- dba.load('AvsQReads')

# Get pearson's correlation across conditions
## Fig. 1C
values <- dba.plotHeatmap(all)
          # aNSC_rep3 aNSC_rep2 qNSC_rep1 qNSC_rep2
# aNSC_rep3      1.00      0.94      0.74      0.77
# aNSC_rep2      0.94      1.00      0.74      0.77
# qNSC_rep1      0.74      0.74      1.00      0.95
# qNSC_rep2      0.77      0.77      0.95      1.00

# Establish contrast across conditions
all <- dba.contrast(all, categories=DBA_CONDITION, minMembers=2)

# Perform differential analysis
all <- dba.analyze(all)
# 4 Samples, 46156 sites in matrix:
#   ID    Tissue    Factor Condition  Treatment Replicate Caller Intervals FRiP
# 1 aNSC_rep2 Postnatal Chromatin Activated Full-Media         2 counts     46156 0.07
# 2 aNSC_rep3 Postnatal Chromatin Activated Full-Media         3 counts     46156 0.06
# 3 qNSC_rep1 Postnatal Chromatin Quiescent BMP4-Media         1 counts     46156 0.06
# 4 qNSC_rep2 Postnatal Chromatin Quiescent BMP4-Media         2 counts     46156 0.05

# 1 Contrast:
#   Group1 Members1    Group2 Members2 DB.DESeq2
# 1 Activated        2 Quiescent        2      9929

# Generate plots
dba.plotPCA(all, label=DBA_ID)
dba.plotPCA(all, label=DBA_ID, b3D=T)

# Retrieve differential sites
all.DB.analysis <- dba.report(all)
write.csv(`all.DB.analysis`, "AvsQ.DB.analysis.csv", row.names=FALSE, quote=FALSE)
write.table(`all.DB.analysis`, "AvsQ.DB.analysis.txt", sep="\t", row.names=FALSE, quote=FALSE)


# ===================================================================================

# Get consensus peaks for each condition

# Load data
all <- dba.load('AvsQPeaks')

# Get consensus peaksets
Qpeaks <- dba(all, mask=all$masks$Quiescent, minOverlap=2)
Qsites <- dba.peakset(Qpeaks, bRetrieve=T)
export(Qsites, "Quies.consPeaks.bed")

Apeaks <- dba(all, mask=all$masks$Activated, minOverlap=2)
Asites <- dba.peakset(Apeaks, bRetrieve=T)
export(Asites, "Act.consPeaks.bed")


# Add consensus peak name for each condition

# Quiescent
df <- read.table("Quies.consPeaks.bed", sep="\t")

dfnew <- cbind(df, 1:nrow(df))
colnames(dfnew) <- c("chr", "start", "end", "site")
dfnew$site <- sub("^", "Quies_site", dfnew$site)
write.table(dfnew, file="Quies.consPeaks.bed", sep="\t", quote=F, row.names=F, col.names=F)

# Activated
df <- read.table("Act.consPeaks.bed", sep="\t")

dfnew <- cbind(df, 1:nrow(df))
colnames(dfnew) <- c("chr", "start", "end", "site")
dfnew$site <- sub("^", "Act_site", dfnew$site)
write.table(dfnew, file="Act.consPeaks.bed", sep="\t", quote=F, row.names=F, col.names=F)

# Retrieve overlapping accessible sites (co-accessible, i.e. stable across conditions)
all_shared <- dba.overlap(all, all$masks$Chromatin)
all_shared$inAll # shows GRanges object, 23,130 sites




# ===================================================================================

### Part II - Reactivated vs Activated

# Notes: Reactivated 2 and 3 libraries, downsampled to 40 million reads
# Compare with Activated (not downsampled)

# ===================================================================================
# Save data with read counts as DBA object in R workspace

# Load data
# all <- dba(sampleSheet="RvsA.postnatal.analysis.csv")

# Save peaks as DB object
# savefile <- dba.save(all, 'RvsAPeaks')

# Count reads
# all <- dba.count(all, minOverlap=2)

# Save read count DB object
# savefile <- dba.save(all, 'RvsAReads')

# ===================================================================================

# Load data
all <- dba.load('RvsAReads')

# Get pearson's correlation across conditions
values <- dba.plotHeatmap(all)
          # aNSC_rep3 aNSC_rep2 rNSC_rep3 rNSC_rep2
# aNSC_rep3      1.00      0.93      0.82      0.82
# aNSC_rep2      0.93      1.00      0.82      0.82
# rNSC_rep3      0.82      0.82      1.00      0.96
# rNSC_rep2      0.82      0.82      0.96      1.00

# Establish contrast across conditions
all <- dba.contrast(all, categories=DBA_CONDITION, minMembers=2)

# Perform differential analysis
all <- dba.analyze(all)
# 4 Samples, 55830 sites in matrix:
      #   ID    Tissue    Factor   Condition  Treatment Replicate Caller Intervals FRiP
# 1 rNSC_rep2 Postnatal Chromatin Reactivated Full-Media         4 counts     55830 0.15
# 2 rNSC_rep3 Postnatal Chromatin Reactivated Full-Media         4 counts     55830 0.15
# 3 aNSC_rep2 Postnatal Chromatin   Activated Full-Media         2 counts     55830 0.08
# 4 aNSC_rep3 Postnatal Chromatin   Activated Full-Media         3 counts     55830 0.06
# 
# 1 Contrast:
#   Group1 Members1    Group2 Members2 DB.DESeq2
# 1 Reactivated        2 Activated        2     44776

# Generate plots
dba.plotPCA(all, label=DBA_ID)
dba.plotPCA(all, label=DBA_ID, b3D=T)

# ===================================================================================

# Get consensus peaks for each condition

# Load data
all <- dba.load('RvsAPeaks')

# Get consensus peaksets
Rpeaks <- dba(all, mask=all$masks$Reactivated, minOverlap=2)
Rsites <- dba.peakset(Rpeaks, bRetrieve=T)
export(Rsites, "React.consPeaks.bed")

# Add consensus peak name for each condition

# Reactivated
df <- read.table("React.consPeaks.bed", sep="\t")

dfnew <- cbind(df, 1:nrow(df))
colnames(dfnew) <- c("chr", "start", "end", "site")
dfnew$site <- sub("^", "React_site", dfnew$site)
write.table(dfnew, file="React.consPeaks.bed", sep="\t", quote=F, row.names=F, col.names=F)
