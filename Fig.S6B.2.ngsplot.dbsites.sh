#!/bin/bash

# =================================================
# Generate ngsplot in AA and AQ chromatin
# for Ascl1, and NFI ChIP-seq signal
# =================================================

cd ~/Dropbox/Desktop/ATAC-seq-PEAKS/AvsQ_analysis_2019_v2/TF_ATAC-seq

ngs.plot.r -G mm10 -R bed -C config.Ascl1.inp.txt -O DBsites.Ascl1.binding

ngs.plot.r -G mm10 -R bed -C config.NFI.inp.txt -O DBsites.NFI.binding



