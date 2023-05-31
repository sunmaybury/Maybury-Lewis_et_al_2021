#!/bin/bash

# ============================================
# Intersecting ATAC-seq and TF ChIP-seq data
# v2 Diffbind AvsQ 9929 DBsites
# October - November 2019
# ============================================

cd ~/Dropbox/Desktop/ATAC-seq-PEAKS/AvsQ_analysis_2019_v2/TF_ATAC-seq/


# Get dynamic Ascl1 binding sites
bedtools intersect -wa -a AvsQ.DBsites.bed -b Ascl1_peaks.bed > DBsites.Ascl1.binding.bed
bedtools intersect -wa -a AA.DBsites.bed -b Ascl1_peaks.bed > AA.Ascl1.binding.bed
bedtools intersect -wa -a AQ.DBsites.bed -b Ascl1_peaks.bed > AQ.Ascl1.binding.bed

# Get dynamic NFI binding sites
bedtools intersect -wa -a AvsQ.DBsites.bed -b NFI.quiNSC_peaks.bed > DBsites.NFI.binding.bed
bedtools intersect -wa -a AA.DBsites.bed -b NFI.quiNSC_peaks.bed > AA.NFI.binding.bed
bedtools intersect -wa -a AQ.DBsites.bed -b NFI.quiNSC_peaks.bed > AQ.NFI.binding.bed

# Get dynamic FOXO3 binding sites
bedtools intersect -wa -a AvsQ.DBsites.bed -b FOXO3_peaks.bed > DBsites.FOXO3.binding.bed
bedtools intersect -wa -a AA.DBsites.bed -b FOXO3_peaks.bed > AA.FOXO3.binding.bed
bedtools intersect -wa -a AQ.DBsites.bed -b FOXO3_peaks.bed > AQ.FOXO3.binding.bed

