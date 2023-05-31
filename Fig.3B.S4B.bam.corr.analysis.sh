#!/bin/bash

# deepTools
# multiBamSummary computes the read coverage genome-wide for two or more BAM files
# output is a compressed numpy array (.npz) which can be used directly to calculate and visualize pairwise correlation values between read coverages using the tool 'plotCorrelation'

# August 14 2018

cd ~/Dropbox/Desktop/ChIP-seq-PEAKS/bam

# compare activated old (AW) and new (SML) H3K4me3 chIP-seq bam files
multiBamSummary bins --bamfiles AW.H3K4me3.sorted.bam act.H3K4me3.sorted.bam -o activated.old.new.results.npz

# correlation plot
plotCorrelation --corData activated.old.new.results.npz --corMethod pearson --skipZeros --whatToPlot scatterplot --plotFile activated.old.new.results.pdf --outFileCorMatrix activated.old.new.results.tsv

## Fig 3B
# compare quiescent and activated new H3K4me3 chIP-seq bam files
multiBamSummary bins --bamfiles quies.H3K4me3.sorted.bam act.H3K4me3.sorted.bam -o quiescent.activated.results.npz

# correlation plot
plotCorrelation --corData quiescent.activated.results.npz --corMethod pearson --skipZeros --whatToPlot scatterplot --plotFile quiescent.activated.results.pdf --outFileCorMatrix quiescent.activated.results.tsv

plotCorrelation --corData quiescent.activated.results.npz --corMethod pearson --skipZeros --log1p --whatToPlot scatterplot --plotFile quiescent.activated.results.log1p.pdf --outFileCorMatrix quiescent.activated.results2.tsv

# compare all 
multiBamSummary bins --bamfiles quies.H3K4me3.sorted.bam act.H3K4me3.sorted.bam AW.H3K4me3.sorted.bam -o all.results.npz

# correlation plot for all
plotCorrelation --corData all.results.npz --corMethod pearson --skipZeros --whatToPlot scatterplot --plotFile all.results.pdf --outFileCorMatrix all.results.tsv



##########################################################


# June 6 2019
## Fig 


# comparison of quies and activated H3K4me3 with other tissues H3K4me3

RenBamFiles="/Users/symaybur/Desktop/bam"


# use merged replicate files from Ren et al ChIP-seq datasets 


cd ~/Dropbox/Desktop/ChIP-seq-PEAKS/bam

# compare quies and activated NSPCs H3K4me3 with heart H3K4me3
multiBamSummary bins --bamfiles quies.H3K4me3.sorted.bam act.H3K4me3.sorted.bam ${RenBamFiles}/heart.H3K4me3.merged.bam ${RenBamFiles}/intestine.H3K4me3.merged.bam ${RenBamFiles}/kidney.H3K4me3.merged.bam ${RenBamFiles}/liver.H3K4me3.merged.bam ${RenBamFiles}/lung.H3K4me3.merged.bam -o NSCs.adult.tissues.H3K4me3.results.npz

# correlation plot for all
plotCorrelation --corData NSCs.adult.tissues.H3K4me3.results.npz --corMethod pearson --skipZeros --log1p --whatToPlot scatterplot --plotFile NSCs.adult.tissues.H3K4me3.results.pdf --outFileCorMatrix NSCs.adult.tissues.H3K4me3.results.tsv



# compare quies and activated NSPCs H3K4me3 with embryonic stem cells and fibroblasts
multiBamSummary bins --bamfiles quies.H3K4me3.sorted.bam act.H3K4me3.sorted.bam ${RenBamFiles}/ESCs.H3K4me3.merged.bam ${RenBamFiles}/MEFs.H3K4me3.merged.bam -o NSCs.embryo.cells.H3K4me3.results.npz

# correlation plot for all
plotCorrelation --corData NSCs.embryo.cells.H3K4me3.results.npz --corMethod pearson --skipZeros --log1p --whatToPlot scatterplot --plotFile NSCs.embryo.cells.H3K4me3.results.pdf --outFileCorMatrix NSCs.embryo.cells.H3K4me3.results.tsv








