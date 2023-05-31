#!/bin/bash

# ngsplot for quiescent young vs old NSCs
# in regions identified as age-DBsites in quiescent cells
# 1389 DiffBind sites

# directory where quiescent young and old merged bam files are
cd ~/Dropbox/Brown_Webb/Desktop/ATAC-seq-PEAKS-AGE/Diffbind/files

# moved quiescent.age-DBsites.bed and input file into this directory
# move these out to ngsplot directory when done

# ==============================================

ngs.plot.r -G mm10 -F cortex -R bed -SC 0,1 -CO blue2 -C config.quies.age.inp.txt -O quies.age-DBsites.ngsplot

# move files to ngsplot directory
mv *.ngsplot.*.pdf ~/Dropbox/Desktop/ATAC-seq-PEAKS-AGE/ngsplot


