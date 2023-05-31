
library(dplyr)

setwd("~/Dropbox/Brown_Webb/Desktop/ATAC-seq-PEAKS/AvsQ_analysis_2019_v2/Expression_distribution_v2")

stable <- read.table("AvsQ.truStable.proximal+distal.genes.only.txt", sep="\t", stringsAsFactors=F)
dynamic <- read.table("AvsQ.DBsites.proximal+distal.genes.only.txt", sep="\t", stringsAsFactors=F)

# remove duplicates
stable <- distinct(stable)
dynamic <- distinct(dynamic)

stable_only <- setdiff(stable, dynamic)
names(stable_only) <- "gene"

# expression in vivo
expres <- read.csv("~/Dropbox/Brown_Webb/Desktop/RNA-seq-data_pub/DEseq_young_aNSC-qNSC_Leeman_et_al.csv", sep=",", stringsAsFactors=F)
expres <- expres[,c(1,3,7)]
expres <- expres[expres$padj < 0.05, ] 
colnames(expres) <- c("gene", "log2FC", "padj")

upreg <- expres[expres$log2FC>0, ]
downreg <- expres[expres$log2FC<0, ]

# merge expression and stable only
stable_only_exp <- inner_join(stable_only, expres)
stable_only_upreg <- inner_join(stable_only, upreg)
stable_only_downreg <- inner_join(stable_only, downreg)
