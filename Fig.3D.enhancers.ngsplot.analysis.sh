#!/bin/bash

# ==================================================================================
# ngsplots for H3K27ac ChIP from proliferating and quiescent NSCs (Martynoga et al)
# In AA and AQ DBsites
# v2 9929 DBsites AvsQ normalized
# October 2019
# ==================================================================================

# Combine H3K27ac and p300 in configuration files
# Separate quiescent and activated
# Run without clustering

cd ~/Dropbox/Desktop/ATAC-seq-PEAKS/AvsQ_analysis_2019_v2/Enhancers_ngsplot

# ======================
# All AvsQ DBsites
# ======================

# AA
ngs.plot.r -G mm10 -F cortex -R bed -C config.AA.all.inp.txt -O AA.all.enriched.enhancers

# AQ
ngs.plot.r -G mm10 -F cortex -R bed -C config.AQ.all.inp.txt -O AQ.all.enriched.enhancers


# ======================
# Distal AvsQ DBsites
# ======================

# AA distal
ngs.plot.r -G mm10 -F cortex -R bed -C config.AA.DISTAL.inp.txt -O AA.distal.enriched.enhancers

# AQ distal
ngs.plot.r -G mm10 -F cortex -R bed -C config.AQ.DISTAL.inp.txt -O AQ.distal.enriched.enhancers


# ======================
# Proximal AvsQ DBsites
# ======================

# AA proximal
ngs.plot.r -G mm10 -F cortex -R bed -C config.AA.proximal.inp.txt -O AA.proximal.enriched.enhancers

# AQ distal
ngs.plot.r -G mm10 -F cortex -R bed -C config.AQ.proximal.inp.txt -O AQ.proximal.enriched.enhancers

