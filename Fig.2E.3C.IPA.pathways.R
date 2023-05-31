# ===================================
# Figure for IPA analysis
# v2 DiffBind AvsQ 9929 DBsites
# November 2019
# ===================================


setwd("~/Dropbox/Desktop/ATAC-seq-PEAKS/AvsQ_analysis_2019_v2/IPA")

library(ggplot2)
library(grid)
library(dplyr)

# Load canonical pathways IPA results
# For genes with only stable chromatin 
# and for genes with stable+dynamic chromatin
stable <- read.table("truStable.IPA.txt", sep="\t", stringsAsFactors=F, header=T)
dynamic <- read.table("DBsites.IPA.txt", sep="\t", stringsAsFactors=F, header=T)

stable <- stable[, 1:2]
dynamic <- dynamic[, 1:2]

colnames(stable) <- c("pathway", "-log(pval)")
colnames(dynamic) <- c("pathway", "-log(pval)")

# List by descending order of -log(p-value)
stable <- stable[with(stable, order(-`-log(pval)`)), ]
stable$pathway <- factor(stable$pathway, levels=stable$pathway[order(stable$`-log(pval)`)])

dynamic <- dynamic[with(dynamic, order(-`-log(pval)`)), ]
dynamic$pathway <- factor(dynamic$pathway, levels=dynamic$pathway[order(dynamic$`-log(pval)`)])

# Take top 12 pathways
stable <- stable[-c(10), ]   # delete row 10 since NER and nucleotide excision repair was repeated
stable <- stable[1:12, ]
dynamic <- dynamic[1:12, ]

# Make bar plots
p1 <- ggplot(stable, aes(x=pathway,y=`-log(pval)`,fill=pathway)) +
  geom_bar(stat="identity", width=0.2) +
  coord_flip() +
  theme_bw() +
  theme(legend.position='none') +
  ggtitle("DE genes with stable chromatin")
p1

p2 <- ggplot(dynamic, aes(x=pathway,y=`-log(pval)`,fill=pathway)) +
  geom_bar(stat="identity", width=0.2) +
  coord_flip() +
  theme_bw() +
  theme(legend.position='none') +
  ggtitle("DE genes with dynamic chromatin")
p2