# ============================================================
# Calculating numbers of upregulated and donwregulated genes 
# with AA and AQ sites
# Leeman et al RNA-seq-data
# v2 AvsQ DBsites 9929
# October 2019
# ============================================================

setwd("~/Dropbox/Desktop/ATAC-seq-PEAKS/AvsQ_analysis_2019_v2/Expression_distribution_v2")

library(dplyr)

# load Leeman RNA-seq data
expres <- read.csv("~/Dropbox/Desktop/RNA-seq-data/DEseq_young_aNSC-qNSC_Leeman_et_al.csv", sep=",", stringsAsFactors=F)
expres <- expres[,c(1,3,7)]
colnames(expres) <- c("gene", "log2FC", "padj")

# separate differentially expressed genes
sig <- expres[expres$padj < 0.05, ]       # 4,313 genes
nonsig <- setdiff(expres, sig)            # 11,016 genes

# upregulated and donwregulated significant DE genes
sig_upreg <- sig[sig$log2FC > 0, ]                           # 2088 genes total upregulated
sig_downreg <- sig[sig$log2FC < 0, ]                         # 2225 genes total downregulated

# ===========================================================

# load chromatin data
chrom <- read.table("AvsQ.DBsites.proximal+distal.genelist.txt", sep="\t", stringsAsFactors=F)

# clean up chromatin data to include DBsite, gene, and log2FC(A/Q)
chrom <- chrom[,c(4,7,5)]
colnames(chrom) <- c("DBsite", "gene", "Fold")

# ===========================================================

# change in chromatin vs no change in chromatin in sig DE genes
sig_dynamic <- merge(sig, chrom, by="gene")                         # 1152 gene-associations with dynamic chromatin

sig_dynamicSimple <- sig_dynamic[,c(1,2,3)]
sig_stable <- setdiff(sig, sig_dynamicSimple)                       # 3,423 gene-associations with stable chromatin

sig_dynamic_nodup <- sig_dynamic[!duplicated(sig_dynamic$gene), ]   # 890 total DE genes with dynamic chromatin
sig_stable_nodup <- sig_stable[!duplicated(sig_stable$gene), ]      # 3,423 total DE genes with stable chromatin

# ===========================================================

# upregulated divided into AA, AQ, or both
upreg_dynamic <- merge(sig_upreg, chrom, by="gene")                         # 552
upreg_dynamic_nodup <- upreg_dynamic[!duplicated(upreg_dynamic$gene), ]     # 419 genes

upreg_AA <- upreg_dynamic[upreg_dynamic$Fold > 0, ]                         # 480
upreg_AQ <- upreg_dynamic[upreg_dynamic$Fold < 0, ]                         # 72

upreg_AA_nodup <- upreg_AA[!duplicated(upreg_AA$gene), ]                    # 362 genes
upreg_AQ_nodup <- upreg_AQ[!duplicated(upreg_AQ$gene), ]                    # 63 genes

upreg_both <- merge(upreg_AA, upreg_AQ, by="gene")
upreg_both_nodup <- upreg_both[!duplicated(upreg_both$gene), ]              # 6 genes

upreg_stable <- merge(sig_upreg, sig_stable, by="gene")                     # 1669 genes

# 356 AA + 57 AQ + 6 both + 1669 stable = 2088 upregulated genes

# ===========================================================

# downregulated divided into AA, AQ, or both
downreg_dynamic <- merge(sig_downreg, chrom, by="gene")                             # 600
downreg_dynamic_nodup <- downreg_dynamic[!duplicated(downreg_dynamic$gene), ]       # 471 genes

downreg_AA <- downreg_dynamic[downreg_dynamic$Fold > 0, ]                           # 441
downreg_AQ <- downreg_dynamic[downreg_dynamic$Fold < 0, ]                           # 159

downreg_AA_nodup <- downreg_AA[!duplicated(downreg_AA$gene), ]                      # 356 genes
downreg_AQ_nodup <- downreg_AQ[!duplicated(downreg_AQ$gene), ]                      # 143 genes

downreg_both <- merge(downreg_AA, downreg_AQ, by="gene")
downreg_both_nodup <- downreg_both[!duplicated(downreg_both$gene), ]                # 28 genes

downreg_stable <- merge(sig_downreg, sig_stable, by="gene")                         # 1754 genes

# 328 AA + 115 AQ + 28 both + 1754 stable = 2225 downregulated genes


# ===========================================================

### Make stacked bar plot in Prism

# Total counts

# 2088 upregulated + 2225 downregulated = 4313 total DE genes

# 890 total DE genes with dynamic chromatin
# 3423 total DE genes with stable chromatin

## Upregulated
# AA 356
# AQ 57
# Both 6
# Total 419

## Downregulated
# AA 328
# AQ 115
# Both 28
# Total 471

# Plot percentage of upregulated and downregulated genes with AA, AQ, or both chromatin states




