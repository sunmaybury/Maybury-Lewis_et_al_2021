
# ==============================================================================
# Benjamini-Hochberg correction on Panther GO biological process analyses
# GO-Slim and complete
# quiescent.age-DBsites genes
# for fig. 5 for Aging Cell submission
# December 2020
# ==============================================================================

setwd("~/Dropbox/Brown_Webb/Desktop/ATAC-seq-PEAKS-AGE/Panther/Benjamini-Hochberg/")

# ===========================
# GO-slim biol processes
# ===========================

# load results
# biological process terms
bh.correct <- read.table("quiescent.age-DBsites.panther.GO-slim-biol.txt", sep="\t", stringsAsFactors=F, header=T)
bh.correct <- `bh.correct`[,c(1,3,6,7)]
colnames(bh.correct) <- c("biological_process", "gene_number", "fold_enrichment", "pvalue")

pvalues <- c(bh.correct$pvalue)
bh.p.vals <- p.adjust(pvalues, method="BH", n=2137)
bh.correct$bh.pvalue <- bh.p.vals

# order by p.adj
bh.correct <- `bh.correct`[with(`bh.correct`, order(`bh.pvalue`)), ]
write.table(`bh.correct`, file="quiescent.age-DBsites.panther.GO-slim-biol.BH.corrected.txt", sep="\t", quote=F, row.names=F,
            col.names=T)

# =================================
# GO complete biological processes
# =================================

# load results
# 596 biological process terms
bh.correct <- read.table("quiescent.age-DBsites.panther.GO-biol.txt", sep="\t", stringsAsFactors=F, header=T)
bh.correct <- `bh.correct`[,c(1,3,6,7)]
colnames(bh.correct) <- c("biological_process", "gene_number", "fold_enrichment", "pvalue")

pvalues <- c(bh.correct$pvalue)
bh.p.vals <- p.adjust(pvalues, method="BH", n=1890)
bh.correct$bh.pvalue <- bh.p.vals

# order by p.adj
bh.correct <- `bh.correct`[with(`bh.correct`, order(`bh.pvalue`)), ]
write.table(`bh.correct`, file="quiescent.age-DBsites.panther.GO-bio.BH.corrected.txt", sep="\t", quote=F, row.names=F,
            col.names=T)



