#!/bin/bash

# deepTools
# multiBamSummary computes the read coverage genome-wide for two or more bam files
# output is a compressed numpy array (.npz) which can be used directly to calculate and visualize pairwise correlation values between read coverage using 'plotCorrelation'

# December 7 2020

cd ~/Dropbox/Brown_Webb/Desktop/ATAC-seq-PEAKS-AGE/Analysis/plotCorrelation

# ======== Quiescent =========

# merge files
samtools merge quies.young.merged.bam 8_GGACTCCT_L003_sort.bam 20_CGTACTAG_L005_sort.bam
samtools merge quies.old.merged.bam 5_GGACTCCT_L005_sort.bam 14_CGTACTAG_L004_sort.bam

# index merged files
samtools index quies.young.merged.bam quies.young.merged.bai
samtools index quies.old.merged.bam quies.old.merged.bai 

# compare young and old merged
multiBamSummary bins --bamfiles quies.young.merged.bam quies.old.merged.bam -o quies.merged.npz

# correlation plot
plotCorrelation --corData quies.merged.npz --corMethod pearson --skipZeros --log1p --whatToPlot scatterplot --plotFile quies.young.vs.aged.merged.results.pdf --outFileCorMatrix quies.young.vs.aged.merged.results.tsv

# ======== Activated =========

# merge files
samtools merge act.young.merged.bam 19_TAAGGCGA_L005_sort.bam 25_TAAGGCGA_L008_sort.bam
samtools merge act.old.merged.bam 4_TCCTGAGC_L005_sort.bam 13_TAAGGCGA_L004_sort.bam

# index merged files
samtools index act.young.merged.bam act.young.merged.bai
samtools index act.old.merged.bam act.old.merged.bai

# compare young and old merged
multiBamSummary bins --bamfiles act.young.merged.bam act.old.merged.bam -o act.merged.npz

# correlation plot
plotCorrelation --corData act.merged.npz --corMethod pearson --skipZeros --log1p --whatToPlot scatterplot --plotFile act.young.vs.aged.merged.results.pdf --outFileCorMatrix act.young.vs.aged.merged.results.tsv 

# ======================================================
# Move all output to age directory
# ======================================================

OUTDIR="/Users/symaybur/Dropbox/Desktop/ATAC-seq-PEAKS-AGE/Analysis/plotCorrelation"

mv *.npz $OUTDIR
mv *.tsv $OUTDIR
mv *.results.pdf $OUTDIR