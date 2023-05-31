#!/bin/bash

# ==============================================
# Intersects between CTCF and stable chromatin
# without sites that are also p300-binding sites
# ==============================================

cd ~/Dropbox/Brown_Webb/Desktop/ATAC-seq-PEAKS/AvsQ_analysis_2019_v2/CTCF/peaks_bed



bedtools intersect -wa -a CTCF_merged_1E-8_peaks.bed -b AvsQ.truStable.bed > CTCF.stable.intersect.bed
# 4452 sites

bedtools intersect -wa -a CTCF_merged_1E-8_peaks.bed -b AvsQ.DBsites.bed > CTCF.DBsites.intersect.bed
# 512 sites

bedtools intersect -wa -a CTCF_merged_1E-8_peaks.bed -b p300.proNSC_peaks.bed > CTCF.p300.intersect.bed
# 282 sites

bedtools intersect -wa -a CTCF.stable.intersect.bed -b CTCF.p300.intersect.bed > CTCF.stable.p300.intersect.bed 
# 145 sites

# Subtract p300 binding sites from stable chromatin with CTCF binding
bedtools subtract -A -a  CTCF.stable.intersect.bed -b CTCF.stable.p300.intersect.bed > CTCF.stable.insulator.bed


bedtools intersect -wa -a CTCF_merged_1E-8_peaks.bed -b AvsQ.truStable.distal.bed > CTCF.stable.distal.intersect.bed
#2731

bedtools intersect -wa -a AvsQ.truStable.distal.bed -b CTCF_merged_1E-8_peaks.bed > test.bed

# Run GREAT analysis -25kb/+10kb from TSS