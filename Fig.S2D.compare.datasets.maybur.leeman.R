setwd("~/Dropbox/Brown_webb/RNA-seq_SML/comparison_datasets/")

## version 1
## Load all RNA-seq datasets
maybur <- read.csv("postnatal.quies.vs.activated.results.csv")    ## 6941 genes
leeman <- read.csv("DEseq_young_aNSC-qNSC_Leeman_et_al.csv")
leeman <- leeman[leeman$padj<0.05, ]

## Intersect upregulated genes
maybur_upreg <- maybur[maybur$log2FoldChange>0, ]    ## 3699 genes
leeman_upreg <- leeman[leeman$log2FoldChange>0, ]    ## 2088 genes

upreg <- merge(maybur_upreg, leeman_upreg, by="X")    ## 1026 genes
upregdf <- upreg[1]
# write.table(upregdf, file="upreg.genes.intersect.maybur.leeman.txt", sep="\t", quote=F, row.names=F)

## Intersect downregulated genes
maybur_downreg <- maybur[maybur$log2FoldChange<0, ]    ## 3242 genes
leeman_downreg <- leeman[leeman$log2FoldChange<0, ]    ## 2225 genes 

downreg <- merge(maybur_downreg, leeman_downreg, by="X")  ## 755 genes
downregdf <- downreg[1]
# write.table(downregdf, file="downreg.genes.intersect.maybur.leeman.txt", sep="\t", quote=F, row.names=F)

# ## Total
# upreg <- rbind(maybur_upreg, leeman_upreg)
# upreg_genes <- as.data.frame(upreg$X)
# upreg_num <- unique(upreg_genes)    # 4761
# 
# downreg <- rbind(maybur_downreg, leeman_downreg)
# downreg_genes <- as.data.frame(downreg$X)
# downreg_num <- unique(downreg_genes)    # 4712
# 
# 4761 + 4712

## Total DE genes in vivo and/or NSPCs 9473

## ===========================================================================================

## version 2
## Load all RNA-seq datasets
maybur <- read.csv("postnatal.quies.vs.activated.v2.results.csv")    ## 10292 genes
leeman <- read.csv("DEseq_young_aNSC-qNSC_Leeman_et_al.csv")
leeman <- leeman[leeman$padj<0.05, ]

## Intersect upregulated genes
maybur_upreg <- maybur[maybur$log2FoldChange>0, ]    ## 5136 genes
leeman_upreg <- leeman[leeman$log2FoldChange>0, ]    ## 2088 genes

upreg <- merge(maybur_upreg, leeman_upreg, by="X")    ## 1306 genes

## Intersect downregulated genes
maybur_downreg <- maybur[maybur$log2FoldChange<0, ]    ## 5156 genes
leeman_downreg <- leeman[leeman$log2FoldChange<0, ]    ## 2225 genes 

downreg <- merge(maybur_downreg, leeman_downreg, by="X")  ## 1100 genes
