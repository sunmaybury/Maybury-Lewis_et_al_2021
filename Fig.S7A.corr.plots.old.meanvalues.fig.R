## Correlation: Leeman RNA-seq data and SML RNA-seq data
## VST-normalized read counts for both
## Old NSCs comparison
## May 2021

setwd("~/Dropbox/Brown_Webb/RNA-seq_SML/correlation/")
library(corrplot)

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
## All old NSCs
# =====================================================

## Subset all old samples
old <- common[,c(1,8,9,10,2,3,4,19,20,21,28,29,30)]
rownames(old) <- old$gene
old[1] <- NULL

## Rename samples
colnames(old) <- c("Leeman_q1","Leeman_q2","Leeman_q3","Leeman_a1","Leeman_a2","Leeman_a3",
                     "Maybury_q1","Maybury_q2","Maybury_q3","Maybury_a1","Maybury_a2","Maybury_a3")

## Get average
old$Leeman_qui <- rowMeans(subset(old[,c("Leeman_q1","Leeman_q2","Leeman_q3")]))
old$Leeman_act <- rowMeans(subset(old[,c("Leeman_a1","Leeman_a2","Leeman_a3")]))
old$maybury_qui <- rowMeans(subset(old[,c("Maybury_q1","Maybury_q2","Maybury_q3")]))
old$maybury_act <- rowMeans(subset(old[,c("Maybury_a1","Maybury_a2","Maybury_a3")]))

old <- old[,c(13:16)]

res <- cor(old)

# =====================================================
## Correlation plots!
# =====================================================

## Different types of plots
# 
# ## All old NSCs
# corrplot(res, type="upper", tl.col="black", tl.srt=45)
# corrplot(res, method="color", tl.col="black")

## PDF plots with pearson correlation numbers

pdf("old.NSCs.combined.mean.corrplot.pdf")
corrplot(res, method="color", cl.cex=0.7, order="original", addrect=2, cl.lim=c(0,1),
         addCoef.col="white", number.cex = 1.25,
         tl.col="black", tl.srt=90,
         title="All old NSCs")
dev.off()


