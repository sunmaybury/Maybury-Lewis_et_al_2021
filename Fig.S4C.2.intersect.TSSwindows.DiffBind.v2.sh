#!/bin/bash

# Get intersects of TSS windows and AvsQ DBsite and truStable proximal bed files
# 1/6/19
# Readjust window to -/+50bp around TSS 
# Based on Tammy Gershon's feedback

cd ~/Dropbox/Desktop/ATAC-seq-PEAKS/AvsQ_analysis_2019_v2/Core_promoter

# -/+50bp around TSS
bedtools intersect -a AvsQ.truStable.proximal.bed -b TSS.window.100bp.bed > AvsQ.truStable.proximal.TSS100bp.bed



# ============================================================================================

# All upregulated, downregulated, and not DE genes with stable promoter chromatin

# cut -f1-3 AvsQ.truStable.proximal.upregulated.genelist.txt > AvsQ.truStable.proximal.upregulated.bed
# cut -f1-3 AvsQ.truStable.proximal.downregulated.genelist.txt > AvsQ.truStable.proximal.downregulated.bed
# cut -f1-3 AvsQ.truStable.proximal.notDE.genelist.txt > AvsQ.truStable.proximal.notDE.bed

bedtools intersect -a AvsQ.truStable.proximal.upregulated.bed -b TSS.window.100bp.bed > AvsQ.truStable.proximal.upregulated.TSS100bp.bed

bedtools intersect -a AvsQ.truStable.proximal.downregulated.bed -b TSS.window.100bp.bed > AvsQ.truStable.proximal.downregulated.TSS100bp.bed

bedtools intersect -a AvsQ.truStable.proximal.notDE.bed -b TSS.window.100bp.bed > AvsQ.truStable.proximal.notDE.TSS100bp.bed


# ============================================================================================

# Top and bottom quartile of upregulated genes with stable promoter chromatin in vivo (Leeman)

# cut -f1-3 AvsQ.truStable.proximal.top.quartile.upreg.genelist.txt > AvsQ.truStable.proximal.top.quartile.upreg.bed

bedtools intersect -a AvsQ.truStable.proximal.top.quartile.upreg.bed -b TSS.window.100bp.bed > AvsQ.truStable.proximal.top.quartile.upreg.TSS100bp.bed 

# cut -f1-3 AvsQ.truStable.proximal.second.quartile.upreg.genelist.txt > AvsQ.truStable.proximal.second.quartile.upreg.bed
# cut -f1-3 AvsQ.truStable.proximal.third.quartile.upreg.genelist.txt > AvsQ.truStable.proximal.third.quartile.upreg.bed
# cut -f1-3 AvsQ.truStable.proximal.bottom.quartile.upreg.genelist.txt > AvsQ.truStable.proximal.bottom.quartile.upreg.bed

bedtools intersect -a AvsQ.truStable.proximal.second.quartile.upreg.bed -b TSS.window.100bp.bed > AvsQ.truStable.proximal.second.quartile.upreg.TSS100bp.bed

bedtools intersect -a AvsQ.truStable.proximal.third.quartile.upreg.bed -b TSS.window.100bp.bed > AvsQ.truStable.proximal.third.quartile.upreg.TSS100bp.bed

bedtools intersect -a AvsQ.truStable.proximal.bottom.quartile.upreg.bed -b TSS.window.100bp.bed > AvsQ.truStable.proximal.bottom.quartile.upreg.TSS100bp.bed   


# ============================================================================================

# Top and bottom quartile of upregulated genes with stable promoter chromatin in vitro (Martynoga)

# cut -f1-3 AvsQ.truStable.proximal.top.quartile.upreg.invitro.genelist.txt > AvsQ.truStable.proximal.top.quartile.upreg.invitro.bed
# cut -f1-3 AvsQ.truStable.proximal.bottom.quartile.upreg.invitro.genelist.txt > AvsQ.truStable.proximal.bottom.quartile.upreg.invitro.bed

bedtools intersect -a AvsQ.truStable.proximal.top.quartile.upreg.invitro.bed -b TSS.window.100bp.bed > AvsQ.truStable.proximal.top.quartile.upreg.invitro.TSS100bp.bed

bedtools intersect -a AvsQ.truStable.proximal.bottom.quartile.upreg.invitro.bed -b TSS.window.100bp.bed > AvsQ.truStable.proximal.bottom.quartile.upreg.invitro.TSS100bp.bed


# ============================================================================================

# Top and bottom quartile of upregulated genes with dynamic promoter chromatin in vivo (Leeman)

# cut -f1-3 AvsQ.DBsites.proximal.top.quartile.upreg.genelist.txt > AvsQ.DBsites.proximal.top.quartile.upreg.bed
# cut -f1-3 AvsQ.DBsites.proximal.bottom.quartile.upreg.genelist.txt > AvsQ.DBsites.proximal.bottom.quartile.upreg.bed

bedtools intersect -a AvsQ.DBsites.proximal.top.quartile.upreg.bed -b TSS.window.100bp.bed > AvsQ.DBsites.proximal.top.quartile.upreg.TSS100bp.bed 

bedtools intersect -a AvsQ.DBsites.proximal.bottom.quartile.upreg.bed -b TSS.window.100bp.bed > AvsQ.DBsites.proximal.bottom.quartile.upreg.TSS100bp.bed
