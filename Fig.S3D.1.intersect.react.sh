#!/bin/bash

# ================================================
# Script to divide R consensus peaks into:
# AQ, AA, and truStable (between A and Q),
# and R-specific
# v2 DBsites 9929 normalized read counts AvsQ
# Adding in downsampled Reactivated 2 and 3
# React.consPeaks from DiffBind 2019 v2
# ================================================

cd ~/Dropbox/Desktop/ATAC-seq-PEAKS/AvsQ_analysis_2019_v2/Reactivated/

# Intersect reactivated consensus peaks with AQ, AA, stable (A_Q) sites 

bedtools intersect -wa -a React.consPeaks.bed -b AQ.DBsites.bed > React.AQ.bed   

bedtools intersect -wa -a React.consPeaks.bed -b AA.DBsites.bed > React.AA.bed   

bedtools intersect -wa -a React.consPeaks.bed -b AvsQ.truStable.bed > React.truStable.bed 

# ==================================================================================
# Notes

# site number calculations:

# R consensus 51939

# React.AQ.rmdup 1393
# React.AA.rmdup 6466
# React.stable.rmdup 19857
# React.unique 24314
# Total 52030

# Ran check.overlap.4a.R
# Some R consensus peaks had both AA and stable designations
# ==================================================================================


# ==================================================================================
# Notes
# bedtools subtract ends up subtracting more peaks with -A flag

# Subtract React.AA and React.AQ from R consensus peaks first using 
# subtract.DB.react.consPeaks.4b.R (Run part 1)
# then intersect with AvsQ.truStable.analysis.bed to avoid this issue
# ==================================================================================

bedtools intersect -wa -a React.consPeaks.noDB.bed -b AvsQ.truStable.analysis.bed > React.truStable.bed


# ==================================================================================
# Now run subtract.DB.react.consPeaks.4b.R (part 2)
# to get unique Reactivated peaks
# ==================================================================================


# ==================================================================================
# Lastly run intersect.react.rmdup.4c.R
# to remove duplicate rows from React.AQ.bed, React.AA.bed, React.truStable.bed
# ==================================================================================


# Final output:

# React.consPeaks.bed (Reactivated consensus peaks between two replicates)
# React.AQ.rmdup.bed (Reactivated consensus peaks that are AQ)
# React.AA.rmdup.bed (Reactivated consensus peaks that are AA)
# React.truStable.rmdup.bed (Reactivated consensus peaks that are stably open between Quiescent and Activated)
# React.unique.bed

