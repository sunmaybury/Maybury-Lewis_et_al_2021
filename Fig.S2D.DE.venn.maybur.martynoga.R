setwd("~/Dropbox/Brown_Webb/RNA-seq_SML/comparison_datasets/")

library(venneuler)

## Upregulated genes
## venn diagram for upregulated genes in common b/w NSPC RNA-seq and NS5 RNA-seq
dat <- venneuler(c(NSPCs=2645, NS5s=906, "NSPCs&NS5s"=1054))
plot(dat)

## Fisher's exact test

## Testing whether overlap of upregulated genes in NSPCs with those in NS5 is significant

# upreg_NSPC = 2645
# upreg_NS5 = 906
# overlap = 1054
# total number of DE genes in NSPCs or NS5 = 9394

upreg.olap <- matrix(c(1054,2645,906,4789), nrow=2, dimnames=list(c("invivo","notNS5"), c("NSPC","notNSPC")))

fisher.test(upreg.olap, alternative="greater")
# Fisher's Exact Test for Count Data
# 
# data:  upreg.olap
# p-value < 2.2e-16
# alternative hypothesis: true odds ratio is greater than 1
# 95 percent confidence interval:
#  1.932967      Inf
# sample estimates:
# odds ratio 
#   2.106201

## Exact pval
result2 <- fisher.test(upreg.olap, alternative="greater")
result2$p.value
# [1] 5.851911e-48

## Downregulated genes
## venn diagram for downregulated genes in common b/w NSPC RNA-seq and NS5 RNA-seq
dat <- venneuler(c(NSPCs=2317, NS5s=1550, "NSPCs&NS5s"=925))
plot(dat)

## Fisher's exact test

## Testing whether overlap of downregulated genes in NSPCs with those in NS5 is significant
downreg.olap <- matrix(c(925,2317,1550,4602), nrow=2, dimnames=list(c("invivo","notNS5"), c("NSPC","notNSPC")))

fisher.test(downreg.olap, alternative="greater")
# Fisher's Exact Test for Count Data
# 
# data:  downreg.olap
# p-value = 0.0002769
# alternative hypothesis: true odds ratio is greater than 1
# 95 percent confidence interval:
#  1.092605      Inf
# sample estimates:
# odds ratio 
#   1.185273 

## Exact pval
result2 <- fisher.test(downreg.olap, alternative="greater")
result2$p.value
# [1] 0.0002768905
