#!/bin/bash

# ========================================================
# Script to check intersecting peaks between
# DBsites and shared or "stable" sites

# "Stable" sites are from dba.overlap between conditions.
# Since there are intersecting sites between differential
# and "stable" sites, they are really TOTAL overlaps.
# ========================================================


# Activated and Quiescent


cd ~/Dropbox/Desktop/ATAC-seq-PEAKS/DiffBind/results.diffbind.2019/v2_AvsQ

# Intersect "shared" sites with DB sites 
bedtools intersect -wa -wb -a AvsQ.shared.bed -b AvsQ.DBsites.bed -r -f 0.90 > temp.bed

# temp.bed has 3154 instances

# Subtract intersects from pn.stable.bed
# Report b without intersects

bedtools intersect -v -r -f 0.90 -a AvsQ.shared.bed -b AvsQ.DBsites.bed > AvsQ.truStable.bed
rm temp.bed




