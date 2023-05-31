setwd("~/Dropbox/Desktop/ATAC-seq-PEAKS/AvsQ_analysis_2019_v2/Reactivated")

library(tidyverse)
library(dplyr)

# ===============================================================================
# Part 1
# Subtracting React.AA and React.AQ sites from Reactivated consensus peaks
# ===============================================================================

consR <- read.table("React.consPeaks.bed", sep="\t", stringsAsFactors=F)

reactAA <- read.table("React.AA.rmdup.bed", sep="\t", stringsAsFactors=F)
reactAQ <- read.table("React.AQ.rmdup.bed", sep="\t", stringsAsFactors=F)

tmp <- setdiff(consR, reactAQ)
tmp2 <- setdiff(tmp, reactAA)

write.table(tmp2, file="React.consPeaks.noDB.bed", sep="\t", quote=F, row.names=F, col.names=F)

# ===============================================================================
# Part 2
# Subtracting React.stable, React.AA, React.AQ peaks to get React.unique sites
# React.stable, React.AA, React.AQ duplicate rows were removed first
# ===============================================================================

# Load consPeaks without DBsites 
consR.noDB <- read.table("React.consPeaks.noDB.bed", sep="\t", stringsAsFactors=F)
# test <- distinct(`consR.noDB`)

# Load truStable overlaps
react.stable <- read.table("React.truStable.rmdup.bed", sep="\t", stringsAsFactors=F)
# test <- distinct(`react.stable`)

R.unique <- setdiff(`consR.noDB`, `react.stable`)
# test <- distinct(`R.unique`)
write.table(`R.unique`, file="React.unique.bed", sep="\t", quote=F, row.names=F, col.names=F)

# ===============================================================================
# Part 3
# Checking overlaps between React.stable and AA and AQ
# Total number React consPeaks 51939
# ===============================================================================

reactAA <- read.table("React.AA.rmdup.bed", sep="\t", stringsAsFactors=F)
reactAQ <- read.table("React.AQ.rmdup.bed", sep="\t", stringsAsFactors=F)
react.stable <- read.table("React.truStable.rmdup.bed", sep="\t", stringsAsFactors=F)
react.unique <- read.table("React.unique.bed", sep="\t", stringsAsFactors=F)

all <- rbind(reactAA,reactAQ,`react.stable`, `react.unique`)    # 52030
dup <- all[duplicated(all), ]   # 91
all2 <- distinct(all)    # 51939

olap1 <- inner_join(dup, reactAA)    # 87 overlaps
olap2 <- inner_join(dup, reactAQ)    # 4 overlaps

# Subtract duplicates that overlap with reactAA and reactAQ from React.stable and rewrite file
new.stable <- setdiff(react.stable, dup)
write.table(`new.stable`, file="React.truStable.rmdup.bed", sep="\t", quote=F, row.names=F, col.names=F)
