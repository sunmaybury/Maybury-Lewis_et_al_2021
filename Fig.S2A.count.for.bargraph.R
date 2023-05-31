# ============================================================
# Calculating numbers of upregulated and donwregulated genes 
# with AA and AQ sites
# Martynoga RNA-seq data
# v2 AvsQ DBsites 9929
# October 2019
# ============================================================

setwd("~/Dropbox/Desktop/ATAC-seq-PEAKS/AvsQ_analysis_2019_v2/Expression_distribution")

library(dplyr)

# Load RNA-seq data
expres <- read.csv("~/Dropbox/Desktop/RNA-seq-data/Martynoga_RNA-seq_complete.csv", stringsAsFactors = F, header = T)
# since log2FC is (quiescent/activated), recalculate and reorder columns
expres$log2FC_AQ <- log2(expres$FPKM_aNSCs / expres$FPKM_qNSCs)
expres <- expres[,c(2,7,6)]
colnames(expres) <- c("gene", "log2FC", "padj")

# separate differentially expressed genes
sig <- expres[expres$padj < 0.05, ]       # 4,435 genes

# upregulated and donwregulated significant DE genes
sig_upreg <- sig[sig$log2FC > 0, ]                           # 1960 genes total upregulated
sig_downreg <- sig[sig$log2FC < 0, ]                         # 2475 genes total downregulated

# ===========================================================

# load chromatin data
chrom <- read.table("AvsQ.DBsites.proximal+distal.genelist.txt", sep="\t", stringsAsFactors=F)

# clean up chromatin data to include DBsite, gene, and log2FC(A/Q)
chrom <- chrom[,c(4,7,5)]
colnames(chrom) <- c("DBsite", "gene", "Fold")

# ===========================================================

# change in chromatin vs no change in chromatin in sig DE genes
sig_dynamic <- merge(sig, chrom, by="gene")                         # 1145 gene-associations with dynamic chromatin

sig_dynamicSimple <- sig_dynamic[,c(1,2,3)]
sig_stable <- setdiff(sig, sig_dynamicSimple)                       # 3,549 gene-associations with stable chromatin

sig_dynamic_nodup <- sig_dynamic[!duplicated(sig_dynamic$gene), ]   # 886 total DE genes with dynamic chromatin
sig_stable_nodup <- sig_stable[!duplicated(sig_stable$gene), ]      # 3,543 total DE genes with stable chromatin

# ===========================================================

# upregulated divided into AA, AQ, or both
upreg_dynamic <- merge(sig_upreg, chrom, by="gene")                         # 560
upreg_dynamic_nodup <- upreg_dynamic[!duplicated(upreg_dynamic$gene), ]     # 414 genes

upreg_AA <- upreg_dynamic[upreg_dynamic$Fold > 0, ]                         # 507
upreg_AQ <- upreg_dynamic[upreg_dynamic$Fold < 0, ]                         # 53

upreg_AA_nodup <- upreg_AA[!duplicated(upreg_AA$gene), ]                    # 377 genes
upreg_AQ_nodup <- upreg_AQ[!duplicated(upreg_AQ$gene), ]                    # 51 genes

upreg_both <- merge(upreg_AA, upreg_AQ, by="gene")
upreg_both_nodup <- upreg_both[!duplicated(upreg_both$gene), ]              # 14 genes

upreg_stable <- merge(sig_upreg, sig_stable_nodup, by="gene")                     # 1546 genes

# 363 AA + 37 AQ + 14 both + 1546 stable = 1960 upregulated genes

# ===========================================================

# downregulated divided into AA, AQ, or both
downreg_dynamic <- merge(sig_downreg, chrom, by="gene")                             # 585
downreg_dynamic_nodup <- downreg_dynamic[!duplicated(downreg_dynamic$gene), ]       # 472 genes

downreg_AA <- downreg_dynamic[downreg_dynamic$Fold > 0, ]                           # 414
downreg_AQ <- downreg_dynamic[downreg_dynamic$Fold < 0, ]                           # 171

downreg_AA_nodup <- downreg_AA[!duplicated(downreg_AA$gene), ]                      # 341 genes
downreg_AQ_nodup <- downreg_AQ[!duplicated(downreg_AQ$gene), ]                      # 150 genes

downreg_both <- merge(downreg_AA, downreg_AQ, by="gene")
downreg_both_nodup <- downreg_both[!duplicated(downreg_both$gene), ]                # 19 genes

downreg_stable <- merge(sig_downreg, sig_stable_nodup, by="gene")                   # 2003 genes

# 322 AA + 131 AQ + 19 both + 2003 stable =  2475 downregulated genes


# ===========================================================

### Make stacked bar plot in Prism

# Total counts

# 1960 upregulated + 2475 downregulated = 4435 total DE genes

# 886 total DE genes with dynamic chromatin
# 3543 total DE genes with stable chromatin

## Upregulated
# AA 363
# AQ 37
# Both 14
# Total 414

## Downregulated
# AA 322
# AQ 131
# Both 19
# Total 472

# Plot percentage of upregulated and downregulated genes with AA, AQ, or both chromatin states




