# ==========================================================
# Postnatal ATAC-seq analysis
# Adding only reactivated samples SML2 and SML3
# Based on highest FRiP scores (diffbind.react.reps.R)
# Repeating analysis after downsampling Q, A, R libraries
# ~ 40,000,000 reads and peaks called
# DiffBind
# September 2019
# R version 3.6.1
# DiffBind version 2.12.0
# ==========================================================

setwd("~/Dropbox/Brown_Webb/Desktop/ATAC-seq-PEAKS/DiffBind_downsample")
library(DiffBind)
library(rtracklayer)
# library(rgl)

# ========================================================
# Comparing across quiescent, activated, reactivated
# ========================================================

# Load data
# all <- dba(sampleSheet="best.reps.analysis.csv")
# all
# 6 Samples, 74475 sites in matrix (106677 total):
      #   ID    Tissue    Factor   Condition  Treatment Replicate Caller Intervals
# 1 rNSC_rep2 Postnatal Chromatin Reactivated Full-Media         4 narrow    106586
# 2 rNSC_rep3 Postnatal Chromatin Reactivated Full-Media         4 narrow    107137
# 3 aNSC_rep2 Postnatal Chromatin   Activated Full-Media         2 narrow     51918
# 4 aNSC_rep3 Postnatal Chromatin   Activated Full-Media         3 narrow     41480
# 5 qNSC_rep1 Postnatal Chromatin   Quiescent BMP4-Media         1 narrow     45874
# 6 qNSC_rep2 Postnatal Chromatin   Quiescent BMP4-Media         2 narrow     37113

# =====================================
# Save peaks as DBA object
# savefile <- dba.save(all, 'bestPeaks')
# =====================================

# Count reads
# all <- dba.count(all, minOverlap=2)
# all
# 6 Samples, 74475 sites in matrix:
      #   ID    Tissue    Factor   Condition  Treatment Replicate Caller Intervals FRiP
# 1 rNSC_rep2 Postnatal Chromatin Reactivated Full-Media         4 counts     74475 0.17
# 2 rNSC_rep3 Postnatal Chromatin Reactivated Full-Media         4 counts     74475 0.16
# 3 aNSC_rep2 Postnatal Chromatin   Activated Full-Media         2 counts     74475 0.09
# 4 aNSC_rep3 Postnatal Chromatin   Activated Full-Media         3 counts     74475 0.07
# 5 qNSC_rep1 Postnatal Chromatin   Quiescent BMP4-Media         1 counts     74475 0.07
# 6 qNSC_rep2 Postnatal Chromatin   Quiescent BMP4-Media         2 counts     74475 0.07

# ========================================
# Save reads as DBA object
# savefile <- dba.save(all, 'bestReads')
# ========================================

# Correlation heatmap based on read count data
# plot(all)

# Get pearson's correlation across conditions
# values <- dba.plotHeatmap(all)    ## Fig S3A
# values
			    # qNSC_rep1 qNSC_rep2 aNSC_rep3 aNSC_rep2 rNSC_rep2 rNSC_rep3
# qNSC_rep1      1.00      0.91      0.74      0.80      0.77      0.80
# qNSC_rep2      0.91      1.00      0.74      0.80      0.77      0.80
# aNSC_rep3      0.74      0.74      1.00      0.87      0.84      0.83
# aNSC_rep2      0.80      0.80      0.87      1.00      0.88      0.88
# rNSC_rep2      0.77      0.77      0.84      0.88      1.00      0.95
# rNSC_rep3      0.80      0.80      0.83      0.88      0.95      1.00





# ========================================================================================

# Load data
all <- dba.load('bestReads')

# Establish contrast across conditions
all <- dba.contrast(all, categories=DBA_CONDITION, minMembers=2)
all
# 3 Contrasts:
      #   Group1 Members1    Group2 Members2
# 1 Reactivated        2 Activated        2
# 2 Reactivated        2 Quiescent        2
# 3   Activated        2 Quiescent        2

# Perform differential analysis
all <- dba.analyze(all)

# Generate plots
dba.plotPCA(all, label=DBA_ID)     ## Fig S3B


# dba.plotPCA(all, label=DBA_ID, b3D=T)
# dba.plotPCA(all, label=DBA_ID, components=2,3)

# dba.plotPCA(all, label=DBA_ID, bLog=F, b3D=T)

dba.plotVolcano(all, contrast=3)
dba.plotVolcano(all, contrast=1)
dba.plotVolcano(all, contrast=2)


# ========================================================================================


# Retrieve differentially accessible and stable sites

############
# A vs Q
############

# Differential
all <- dba.load('bestReads')
all <- dba.contrast(all, group1=all$masks$Activated, group2=all$masks$Quiescent, name1="Activated", name2="Quiescent")
all <- dba.analyze(all)

all <- dba.report(all)
# write.csv(all, "AvsQ.DB.analysis.best.csv", row.names=FALSE, quote=FALSE)
# write.table(all, "AvsQ.DB.analysis.best.txt", sep="\t", row.names=FALSE, quote=FALSE)

# Stable
all <- dba.load('bestReads')
all <- dba.contrast(all, group1=all$masks$Activated, group2=all$masks$Quiescent, name1="Activated", name2="Quiescent")
all <- dba.analyze(all)

df <- dba.overlap(all, all$masks$Activated | all$masks$Quiescent)
df$inAll
export(df$inAll, "AvsQ.shared.bed")

############
# R vs A
############

# Differential
all <- dba.load('bestReads')
all <- dba.contrast(all, group1=all$masks$Reactivated, group2=all$masks$Activated, name1="Reactivated", name2="Activated")
all <- dba.analyze(all)

all <- dba.report(all)
# write.csv(all, "RvsA.DB.analysis.best.csv", row.names=FALSE, quote=FALSE)
# write.table(all, "RvsA.DB.analysis.best.txt", sep="\t", row.names=FALSE, quote=FALSE)

# Stable
all <- dba.load('bestReads')
all <- dba.contrast(all, group1=all$masks$Reactivated, group2=all$masks$Activated, name1="Reactivated", name2="Activated")
all <- dba.analyze(all)

df <- dba.overlap(all, all$masks$Reactivated | all$masks$Activated)
df$inAll
export(df$inAll, "RvsA.shared.bed")


############
# R vs Q
############

# Differential
all <- dba.load('bestReads')
all <- dba.contrast(all, group1=all$masks$Reactivated, group2=all$masks$Quiescent, name1="Reactivated", name2="Quiescent")
all <- dba.analyze(all)

all <- dba.report(all)
# write.csv(all, "RvsQ.DB.analysis.best.csv", row.names=FALSE, quote=FALSE)
# write.table(all, "RvsQ.DB.analysis.best.txt", sep="\t", row.names=FALSE)

# Stable
all <- dba.load('bestReads')
all <- dba.contrast(all, group1=all$masks$Reactivated, group2=all$masks$Quiescent, name1="Reactivated", name2="Quiescent")
all <- dba.analyze(all)

df <- dba.overlap(all, all$masks$Reactivated | all$masks$Quiescent)
df$inAll
export(df$inAll, "RvsQ.shared.bed")



# ========================================================================================


# Format output files for further analysis

setwd("~/Dropbox/Desktop/ATAC-seq-PEAKS/DiffBind_downsample/results/best_react_reps")

###########################################################
# in bash
# sort -k1,1 -k2,2n AvsQ.shared.bed > AvsQ.shared.sorted.bed
# cut -f1-3 AvsQ.shared.sorted.bed > AvsQ.nondiff.bed
# sort -k1,1 -k2,2n RvsA.shared.bed > RvsA.shared.sorted.bed
# cut -f1-3 RvsA.shared.sorted.bed > RvsA.nondiff.bed
# sort -k1,1 -k2,2n RvsQ.shared.bed > RvsQ.shared.sorted.bed
# cut -f1-3 RvsQ.shared.sorted.bed > RvsQ.nondiff.bed
###########################################################

AQ <- read.table(file="AvsQ.nondiff.bed", sep="\t", header=FALSE, stringsAsFactors=FALSE)
AQ_stable <- cbind(AQ, 1:nrow(AQ))
colnames(AQ_stable) <- c("chr", "start", "end", "stable_site")
AQ_stable$stable_site <- sub("^", "site", AQ_stable$stable_site)

write.table(AQ_stable, file="AvsQ.stable.analysis.bed", sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)

RA <- read.table(file="RvsA.nondiff.bed", sep="\t", header=FALSE, stringsAsFactors=FALSE)
RA_stable <- cbind(RA, 1:nrow(RA))
colnames(RA_stable) <- c("chr", "start", "end", "stable_site")
RA_stable$stable_site <- sub("^", "site", RA_stable$stable_site)

write.table(RA_stable, file="RvsA.stable.analysis.bed", sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)

RQ <- read.table(file="RvsQ.nondiff.bed", sep="\t", header=FALSE, stringsAsFactors=FALSE)
RQ_stable <- cbind(RQ, 1:nrow(RQ))
colnames(RQ_stable) <- c("chr", "start", "end", "stable_site")
RQ_stable$stable_site <- sub("^", "site", RQ_stable$stable_site)

write.table(RQ_stable, file="RvsQ.stable.analysis.bed", sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)




# ========================================================================================


# Get consensus peaks for each condition

setwd("~/Dropbox/Desktop/ATAC-seq-PEAKS/DiffBind_downsample/results/best_react_reps")


# Load data
peaks <- dba.load('bestPeaks')

# Get consensus peak sets
Qpeaks <- dba(peaks, mask=peaks$masks$Quiescent, minOverlap=2)    #28559
Qsites <- dba.peakset(Qpeaks, bRetrieve=TRUE)
export(Qsites, "Quies.consPeaks.bed")      

Apeaks <- dba(peaks, mask=peaks$masks$Activated, minOverlap=2)    #32729
Asites <- dba.peakset(Apeaks, bRetrieve=TRUE)
export(Asites, "Act.consPeaks.bed")

Rpeaks <- dba(peaks, mask=peaks$masks$Reactivated, minOverlap=2)    #69957
Rsites <- dba.peakset(Rpeaks, bRetrieve=TRUE)
export(Rsites, "React.consPeaks.bed")


# ========================================================================================


# Add consensus peak name for each condition

setwd("~/Dropbox/Desktop/ATAC-seq-PEAKS/DiffBind_downsample/results/best_react_reps")

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

# Reactivated
df <- read.table("React.consPeaks.bed", sep="\t")

dfnew <- cbind(df, 1:nrow(df))
colnames(dfnew) <- c("chr", "start", "end", "site")
dfnew$site <- sub("^", "React_site", dfnew$site)
write.table(dfnew, file="React.consPeaks.bed", sep="\t", quote=F, row.names=F, col.names=F)