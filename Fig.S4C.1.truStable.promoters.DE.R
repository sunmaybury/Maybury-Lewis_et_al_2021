
library(dplyr)

setwd("~/Dropbox/Brown_Webb/Desktop/ATAC-seq-PEAKS/AvsQ_analysis_2019_v2/Core_promoter")

# Load AvsQ proximal stable gene list
df <- read.table("AvsQ.truStable.proximal.genelist.txt", sep="\t", stringsAsFactors=F)
df[4] <- NULL
colnames(df) <- c("chr", "start", "end", "gene")

# Load RNA-seq data (Leeman et al in vivo)
expres <- read.csv("~/Dropbox/Desktop/RNA-seq-data/DEseq_young_aNSC-qNSC_Leeman_et_al.csv", 
                   sep=",", stringsAsFactors=F)
expres <- expres[,c(1,3,7)]
colnames(expres) <- c("gene", "log2FC", "padj")

# Separate differentially expressed genes
sig <- expres[expres$padj < 0.05, ]
upreg <- sig[sig$log2FC>0, ]
downreg <- sig[sig$log2FC<0, ]

notDE <- setdiff(expres, sig)

# ================================================================================

# Separate AvsQ proximal stable gene list by differential expression

df.upreg <- inner_join(df, upreg, by="gene")     # 1683
df.downreg <- inner_join(df, downreg, by="gene") # 1133
df.notDE <- inner_join(df, notDE, by="gene")     # 7497

# Write output files
# write.table(`df.upreg`, file="AvsQ.truStable.proximal.upregulated.genelist.txt", sep="\t", 
#             quote=F, row.names=F, col.names=F)
# write.table(`df.downreg`, file="AvsQ.truStable.proximal.downregulated.genelist.txt", sep="\t",
#             quote=F, row.names=F, col.names=F)
# write.table(`df.notDE`, file="AvsQ.truStable.proximal.notDE.genelist.txt", sep="\t",
#             quote=F, row.names=F, col.names=F)

# ================================================================================

# For expression analysis v2

# Get AvsQ proximal stable genes that are upregulated by ascending padj
df.upreg.reorder <- `df.upreg`[order(`df.upreg`$padj), ] #1683

# Take top quartile of upregulation
df.upreg.top <- `df.upreg.reorder`[1:420, ]
write.table(`df.upreg.top`, file="AvsQ.truStable.proximal.top.quartile.upreg.genelist.txt", 
            sep="\t", quote=F, row.names=F, col.names=F)

# Other quartiles of upregulated
df.upreg.second <- `df.upreg.reorder`[421:841, ]
df.upreg.third <- `df.upreg.reorder`[842:1262, ]
df.upreg.bottom <- `df.upreg.reorder`[421:1683, ]   # all excluding top quartile

write.table(`df.upreg.second`, file="AvsQ.truStable.proximal.second.quartile.upreg.genelist.txt",
            sep="\t", quote=F, row.names=F, col.names=F)
write.table(`df.upreg.third`, file="AvsQ.truStable.proximal.third.quartile.upreg.genelist.txt",
            sep="\t", quote=F, row.names=F, col.names=F)
write.table(`df.upreg.bottom`, file="AvsQ.truStable.proximal.bottom.quartile.upreg.genelist.txt",
            sep="\t", quote=F, row.names=F, col.names=F)





