setwd("~/Dropbox/Brown_Webb/Desktop/ATAC-seq-PEAKS/AvsQ_analysis_2019_v2/GREAT_assign")
library(dplyr)

# Genes associated with truStable sites
df <- read.table("AvsQ.truStable.proximal+distal.genes.txt", sep="\t", stringsAsFactors=F)  # 21079 genes
df2 <- distinct(df)    # 12779 genes

# Genes associated with DBsites
dat <- read.table("AvsQ.DBsites.proximal+distal.genes.txt", sep="\t", stringsAsFactors=F)   # 4074 genes
dat2 <- distinct(dat)    # 3252 genes

stableOnly <- setdiff(df2, dat2)
# 10683 genes are exclusively associated with stable open chromatin
