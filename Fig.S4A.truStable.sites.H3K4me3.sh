#!/bin/bash

# =====================================
# How many of quies and act H3K4me3 
# are in AvsQ.truStable sites?
# v2 9929 DBsites AvsQ DiffBind
# 19976 AvsQ truStable sites
# November 2019
# =====================================

cd ~/Dropbox/Desktop/ATAC-seq-PEAKS/AvsQ_analysis_2019_v2/H3K4me3


# bedtools intersect -u writes original A entry once if any overlaps are found in B (at least one overlap found in B)

# truStable quies and act H3K4me3 peaks
bedtools intersect -u -a quies_H3K4me3_peaks.bed -b AvsQ.truStable.bed > quies.H3K4me3.truStable.bed
bedtools intersect -u -a act_H3K4me3_peaks.bed -b AvsQ.truStable.bed > act.H3K4me3.truStable.bed


# Quies and act H3K4me3 in AvsQ DBsites
bedtools intersect -u -a quies_H3K4me3_peaks.bed -b AvsQ.DBsites.bed > quies.H3K4me3.DBsites.bed
bedtools intersect -u -a act_H3K4me3_peaks.bed -b AvsQ.DBsites.bed > act.H3K4me3.DBsites.bed



