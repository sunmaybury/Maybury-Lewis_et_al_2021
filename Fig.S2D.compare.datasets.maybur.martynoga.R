setwd("~/Dropbox/Brown_webb/RNA-seq_SML/comparison_datasets/")

## version 1
## Load NSPCs RNA-seq dataset
maybur <- read.csv("postnatal.quies.vs.activated.results.csv")    ## 6941 genes
names(maybur)[1] <- "gene"

## Intersect upregulated genes
maybur_upreg <- maybur[maybur$log2FoldChange>0, ]    ## 3699 genes

martynoga_upreg <- read.csv("Martynoga_RNA-seq_Down_qNSCs.csv")  ## 1960 genes
martynoga_upreg[1] <- NULL
names(martynoga_upreg)[1] <- "gene"

upreg <- merge(maybur_upreg, martynoga_upreg, by="gene")    ## 1054 genes
upregdf <- upreg[1]
# write.table(upregdf, file="upreg.genes.intersect.maybur.martynoga.txt", sep="\t", quote=F, row.names=F)


## Intersect downregulated genes
maybur_downreg <- maybur[maybur$log2FoldChange<0, ]    ## 3242 genes

martynoga_downreg <- read.csv("Martynoga_RNA-seq_Up_qNSCs.csv")  ## 2475 genes
martynoga_downreg[1] <- NULL
names(martynoga_downreg)[1] <- "gene"

downreg <- merge(maybur_downreg, martynoga_downreg, by="gene")  ## 925 genes
downregdf <- downreg[1]
# write.table(downregdf, file="downreg.genes.intersect.maybur.martynoga.txt", sep="\t", quote=F, row.names=F)


## Total DE genes
## For Fisher's exact test calculations
maybur_upreg <- as.data.frame(maybur_upreg$gene)
names(maybur_upreg) <- "gene"
maybur_downreg <- as.data.frame(maybur_downreg$gene)
names(maybur_downreg) <- "gene"
martynoga_upreg <- as.data.frame(martynoga_upreg$gene)
names(martynoga_upreg) <- "gene"
martynoga_downreg <- as.data.frame(martynoga_downreg$gene)
names(martynoga_downreg) <- "gene"

upreg <- rbind(maybur_upreg, martynoga_upreg)
upreg_num <- unique(upreg)    # 4605

downreg <- rbind(maybur_downreg, martynoga_downreg)
downreg_num <- unique(downreg)    # 4789
# 
# 4605 + 4789

## Total DE genes in vivo and/or NSPCs 9394




