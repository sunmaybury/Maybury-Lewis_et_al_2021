setwd("~/Dropbox/Brown_Webb/RNA-seq_SML/comparison_datasets/")

library(venneuler)

# Upregulated genes
## venn diagram for upregulated genes in common b/w NSPC RNA-seq and in vivo RNA-seq
dat <- venneuler(c(NSPCs=2673, NSCs=1062, "NSPCs&NSCs"=1026))
plot(dat)

## Fisher's exact test

## Testing whether overlap of upregulated genes in NSPCs with those in vivo is significant

# upreg_NSPC = 2673
# upreg_invivo = 1062
# overlap = 1026
# total number of DE genes in NSPCs or in vivo = 9473

upreg.olap <- matrix(c(1026,2673,1062,4712), nrow=2, dimnames=list(c("invivo","notinvivo"), c("NSPC","notNSPC")))

fisher.test(upreg.olap, alternative="greater")
# Fisher's Exact Test for Count Data
# 
# data:  upreg.olap
# p-value < 2.2e-16
# alternative hypothesis: true odds ratio is greater than 1
# 95 percent confidence interval:
#  1.56649     Inf
# sample estimates:
# odds ratio 
#   1.703014 

# exact pval 
result2 <- fisher.test(upreg.olap, alternative="greater")
result2$p.value
# [1] 1.589536e-26

## Downregulated genes
## venn diagram for downregulated genes in common b/w NSPC RNA-seq and in vivo RNA-seq
dat <- venneuler(c(NSPCs=2487, NSCs=1470, "NSPCs&NSCs"=755))
plot(dat)

# downreg_NSPC = 2487
# downreg_invivo = 1470
# overlap = 755
# total number of downreg genes in NSPCs or in vivo = 4712

downreg.olap <- matrix(c(755,2487,1470,4761), nrow=2, dimnames=list(c("invivo","notinvivo"), c("NSPC","notNSPC")))

fisher.test(downreg.olap, alternative="greater")
# Fisher's Exact Test for Count Data
# 
# data:  downreg.olap
# p-value = 0.6387
# alternative hypothesis: true odds ratio is greater than 1
# 95 percent confidence interval:
#  0.9025755       Inf
# sample estimates:
# odds ratio 
#  0.9832237 




