#!/bin/bash

# run commands on Oscar
module load bedtools/2.26.0

# ====================================
# background
# ====================================

# bed - AvsQ.proximal.TSS100bp regions
# intersect betweeh AvsQ stable proximal and TSS -/+50bp window

bedtools getfasta -fi /users/symaybur/data/ref_genomes/mm10.fa -bed AvsQ.truStable.proximal.TSS100bp.bed -fo AvsQ.truStable.proximal.TSS100bp.fa

bedtools getfasta -fi /users/symaybur/data/ref_genomes/mm10.fa -bed TSS.window.100bp.bed -fo TSS.window.100bp.fa


# =====================================================================================
# AvsQ stable proximal TSS -/+50bp window of upregulated genes by quartile in vivo
# =====================================================================================

bedtools getfasta -fi /users/symaybur/data/ref_genomes/mm10.fa -bed AvsQ.truStable.proximal.top.quartile.upreg.TSS100bp.bed -fo AvsQ.truStable.proximal.top.quartile.upreg.TSS100bp.fa

bedtools getfasta -fi /users/symaybur/data/ref_genomes/mm10.fa -bed AvsQ.truStable.proximal.second.quartile.upreg.TSS100bp.bed -fo AvsQ.truStable.proximal.second.quartile.upreg.TSS100bp.fa

bedtools getfasta -fi /users/symaybur/data/ref_genomes/mm10.fa -bed AvsQ.truStable.proximal.third.quartile.upreg.TSS100bp.bed -fo AvsQ.truStable.proximal.third.quartile.upreg.TSS100bp.fa

bedtools getfasta -fi /users/symaybur/data/ref_genomes/mm10.fa -bed AvsQ.truStable.proximal.bottom.quartile.upreg.TSS100bp.bed -fo AvsQ.truStable.proximal.bottom.quartile.upreg.TSS100bp.fa

# =====================================================================================
# AvsQ stable proximal TSS -/+50bp window of upregulated genes by quartile in vitro
# =====================================================================================

bedtools getfasta -fi /users/symaybur/data/ref_genomes/mm10.fa -bed AvsQ.truStable.proximal.top.quartile.upreg.invitro.TSS100bp.bed -fo AvsQ.truStable.proximal.top.quartile.upreg.invitro.TSS100bp.fa

bedtools getfasta -fi /users/symaybur/data/ref_genomes/mm10.fa -bed AvsQ.truStable.proximal.bottom.quartile.upreg.invitro.TSS100bp.bed -fo AvsQ.truStable.proximal.bottom.quartile.upreg.invitro.TSS100bp.fa




