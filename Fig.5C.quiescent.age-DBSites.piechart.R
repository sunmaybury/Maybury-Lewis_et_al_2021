
# =============================================
# Quiescent young vs aged DBsites
# Generate pie chart for genome distribution
# FOr Fig 5 Aging Cell submission
# December 2020
# =============================================

setwd("~/Dropbox/Brown_Webb/Desktop/ATAC-seq-PEAKS-AGE/Genomic_regions_chart/")

library(ChIPseeker)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(clusterProfiler)
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene

# convert diffbind output csv to bed
res <- read.csv("quiescent.aged.two.reps.model.batch.csv")[,1:3]
write.table(res, file="quiescent.aged.two.reps.model.batch.bed", sep="\t", quote=F, row.names=F, col.names=F)

# load bed file
quies <- read.table("quiescent.aged.two.reps.model.batch.bed", sep="\t", stringsAsFactors=F)

# make GRanges
quiesRanges <- GRanges(seqnames=quies$V1, ranges=IRanges(quies$V2, quies$V3), strand=NULL)

# prepare TSS regions to calculate profile 
promoter <- getPromoters(TxDb=txdb, upstream=1000, downstream=1000)

# peak annotation and pie chart
tagMatrix <- getTagMatrix(quiesRanges, windows=promoter)
quiesPeakAnno <- annotatePeak(quiesRanges, tssRegion=c(-1000, 1000), TxDb=txdb)

pdf("quies.age.DBsites.piechart.pdf")
plotAnnoPie(quiesPeakAnno)
dev.off()
