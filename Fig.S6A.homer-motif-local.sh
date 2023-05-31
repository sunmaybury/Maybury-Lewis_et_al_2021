#!/bin/bash

# motif analysis
# homer

# ==============================
# 10/23/19

cd ~/Dropbox/Brown_Webb/Desktop/ATAC-seq-PEAKS/AvsQ_analysis_2019_v2/Homer_motif_analysis
perl /Users/symaybur/Softwares/homer-4.10/.//configureHomer.pl -install mm10

# AA and AQ DBsites
# Based on DiffBind results using downsampled A, Q, R ATAC-seq results
findMotifsGenome.pl AA.DBsites.bed mm10 AA-DBsites/ -size 200
findMotifsGenome.pl AQ.DBsites.bed mm10 AQ-DBsites/ -size 200

# AvsQ truStable sites
findMotifsGenome.pl AvsQ.truStable.bed mm10 AvsQ-truStable/ -size 200

# ==============================
# 05/27/20

# NFI binding AQ sites
# ASCL1 binding AA sites

cd ~/Dropbox/Brown_Webb/Desktop/ATAC-seq-PEAKS/AvsQ_analysis_2019_v2/Homer_motif_analysis
perl /Users/symaybur/Software/homer/.//configureHomer.pl -install mm10

# NFI binding AQ sites
findMotifsGenome.pl AQ.NFI.target.bed mm10 AQ-NFI/ -size 200

# ASCL1 binding AA sites
findMotifsGenome.pl AA.ASCL1.target.bed mm10 AA-ASCL1/ -size 200

# ==============================
# 06/22/20

# NFI binding AQ sites against AA sites as background
# ASCL1 binding AA sites against AQ sites as background

cd ~/Dropbox/Brown_Webb/Desktop/ATAC-seq-PEAKS/AvsQ_analysis_2019_v2/Homer_motif_analysis

# NFI binding AQ sites against AA sites as background
findMotifsGenome.pl AQ.NFI.target.bed mm10 AQ-NFI-AA/ -size 200 -bg AA.DBsites.bed

# ASCL1 binding AA sites against AQ sites as background
findMotifsGenome.pl AA.ASCL1.target.bed mm10 AA-ASCL1-AQ/ -size 200 -bg AQ.DBsites.bed

# NFI binding AQ sites against AQ sites as background
findMotifsGenome.pl AQ.NFI.target.bed mm10 AQ-NFI-AQ/ -size 200 -bg AQ.DBsites.bed

# ASCL1 binding AA sites against AA sites as background
findMotifsGenome.pl AA.ASCL1.target.bed mm10 AA-ASCL1-AA/ -size 200 -bg AA.DBsites.bed


