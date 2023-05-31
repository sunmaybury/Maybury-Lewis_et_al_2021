## Correlation: Leeman RNA-seq data and SML RNA-seq data
## VST-normalized read counts for both
## Young NSCs comparison
## May 2021

setwd("~/Dropbox/Brown_Webb/RNA-seq_SML/correlation/")
library(corrplot)
library(PerformanceAnalytics)

## Load data
leeman <- read.csv("Leeman_normExpressionValues.csv")
maybur <- read.csv("VST-normalized.read.counts.csv")

names(leeman)[1] <- "gene"
names(maybur)[1] <- "gene"

## Leeman data has some genes converted to dates, take them out
leeman <- leeman[!grepl("-Mar",leeman$gene), ]
leeman <- leeman[!grepl("-Sep",leeman$gene), ]

## Get intersecting genes
common <- merge(leeman, maybur)

# =====================================================
## All young NSCs
# =====================================================

## Subset all young samples
young <- common[,c(1,11,12,5,6,7,16,17,18,25,26,27)]
rownames(young) <- young$gene
young[1] <- NULL

## Rename samples
colnames(young) <- c("Leeman_q1","Leeman_q2","Leeman_a1","Leeman_a2","Leeman_a3",
                     "Maybury_q1","Maybury_q2","Maybury_q3","Maybury_a1","Maybury_a2","Maybury_a3")

## Get average
young$Leeman_qui <- rowMeans(subset(young[,c("Leeman_q1","Leeman_q2")]))
young$Leeman_act <- rowMeans(subset(young[,c("Leeman_a1","Leeman_a2","Leeman_a3")]))
young$maybury_qui <- rowMeans(subset(young[,c("Maybury_q1","Maybury_q2","Maybury_q3")]))
young$maybury_act <- rowMeans(subset(young[,c("Maybury_a1","Maybury_a2","Maybury_a3")]))

young <- young[,c(12:15)]

# Get correlation
res <- cor(young)

# =====================================================
## Correlation plots!
# =====================================================

## Different types of plots
# 
# ## All young NSCs
# corrplot(res, type="upper", tl.col="black", tl.srt=45)
# corrplot(res, method="color", tl.col="black")

## PDF plots with pearson correlation numbers

pdf("young.NSCs.combined.mean.corrplot.pdf")
corrplot(res, method="color", cl.cex=0.7, order="original", addrect=2, cl.lim=c(0,1),
         addCoef.col="white", number.cex = 1.25,
         tl.col="black", tl.srt=90,
         title="All young NSCs")
dev.off()




