
# =============================================================================
# Activated, Quiescent, Reactivated consensus peaks (consPeaks) from DiffBind
# Using these peaks instead of merged replicates peaks
# DBsites (AvsQ) and truStable sites
# Generating pie charts for genome distribution
# November 2019
# =============================================================================

setwd("~/Dropbox/Desktop/ATAC-seq-PEAKS/AvsQ_analysis_2019_v2/Genomic_regions_chart")

library(ChIPseeker)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(clusterProfiler)
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene

# =================================
# consensusPeaks in each cell state
# =================================

# Load bed files
quies <- read.table("Quies.consPeaks.bed", sep="\t", stringsAsFactors=F)
act <- read.table("Act.consPeaks.bed", sep="\t", stringsAsFactors=F)
react <- read.table("React.consPeaks.bed", sep="\t", stringsAsFactors=F)

# Make Granges
quiesRanges <- GRanges(seqnames=quies$V1, ranges=IRanges(quies$V2, quies$V3), strand=NULL)
actRanges <- GRanges(seqnames=act$V1, ranges=IRanges(act$V2, act$V3), strand=NULL)
reactRanges <- GRanges(seqnames=react$V1, ranges=IRanges(react$V2, react$V3), strand=NULL)

# Prepare TSS regions to calculate the profile of ATAC-seq peaks
# Changing promoter (proximal) designation to 1kb upstream, 1kb downstream from first version of analysis
promoter <- getPromoters(TxDb=txdb, upstream=1000, downstream=1000)

# Peak annotation and make pie chart- Quiescent
tagMatrix <- getTagMatrix(quiesRanges, windows=promoter)
quiesPeakAnno <- annotatePeak(quiesRanges, tssRegion=c(-1000, 1000), TxDb=txdb)
plotAnnoPie(quiesPeakAnno)  

# Peak annotation and make pie chart- Activated
tagMatrix <- getTagMatrix(actRanges, windows=promoter)
actPeakAnno <- annotatePeak(actRanges, tssRegion=c(-1000, 1000), TxDb=txdb)
plotAnnoPie(actPeakAnno)

# Peak annotation and make pie chart - Reactivated
tagMatrix <- getTagMatrix(reactRanges, windows=promoter)
reactPeakAnno <- annotatePeak(reactRanges, tssRegion=c(-1000,1000), TxDb=txdb)
plotAnnoPie(reactPeakAnno)


# =================================
# DBsites and truStable sites AvsQ
# =================================

# Load bed files
dbsites <- read.table("AvsQ.DBsites.bed", sep="\t", stringsAsFactors=F)
stable <- read.table("AvsQ.truStable.bed", sep="\t", stringsAsFactors=F)

# Make GRanges
dbRanges <- GRanges(seqnames=dbsites$V1, ranges=IRanges(dbsites$V2, dbsites$V3), strand=NULL)
stableRanges <- GRanges(seqnames=stable$V1, ranges=IRanges(stable$V2, stable$V3), strand=NULL)

# Prepare TSS regions to calculate the profile of ATAC-seq peaks
# Changing promoter (proximal) designation to 1kb upstream, 1kb downstream from first version of analysis
promoter <- getPromoters(TxDb=txdb, upstream=1000, downstream=1000)

# Peak annotation and make pie chart - DBsites
tagMatrix <- getTagMatrix(dbRanges, windows=promoter)
dbPeakAnno <- annotatePeak(dbRanges, tssRegion=c(-1000,1000), TxDb=txdb)
plotAnnoPie(dbPeakAnno)

# Peak annotation and make pie chart - truStable sites
tagMatrix <- getTagMatrix(stableRanges, windows=promoter)
stablePeakAnno <- annotatePeak(stableRanges, tssRegion=c(-1000,1000), TxDb=txdb)
plotAnnoPie(stablePeakAnno)

