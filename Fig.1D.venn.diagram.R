###################################################
# venn diagrams for AQ, AA, stable peak #s
# R version 3.6.1
# venneuler 
###################################################

setwd("~/Dropbox/Desktop/ATAC-seq-PEAKS/DiffBind/results.diffbind.2019/v2_AvsQ")

library(venneuler)

# venn diagram for Activated vs Quiescent
dat <- venneuler(c(AA=6777, AQ=3152, "AA&AQ"=19976))
plot(dat)




