# ===================================
# Make weighted venn diagrams
# CTCF binding open chromatin
# CTCF enriched in stable chromatin?
# Fisher's exact test
# v2 DiffBind AvsQ 9929 DBsites
# November 2019
# ===================================

setwd("~/Dropbox/Brown_Webb/Desktop/ATAC-seq-PEAKS/AvsQ_analysis_2019_v2/CTCF")

library(venneuler)

# =========================================
# CTCF total binding 27797 sites
# AvsQ truStable 19976 sites
# DBsites 9929 sites
# CTCF + truStable 4452 sites
# CTCF + DBsites 512 sites
# =========================================

# Make venn diagram

# CTCF in open chromatin
CTCF_stable <- venneuler(c(stable=15524, CTCF=512, "stable&CTCF"=4452))
plot(CTCF_stable)

# Test: are stable chromatin sites enriched with CTCF binding?
table <- matrix(c(4452,15524,512,9417), nrow=2, dimnames = list(c("CTCF", "notCTCF"), c("stable", "notstable")))
        # stable notstable
# CTCF      4452       512
# notCTCF  15524      9417                        

# Odds ratio
odds_CTCF <- 4452 / 512
odds_notCTCF <- 15524 / 9417
odds_ratio <- odds_CTCF / odds_notCTCF   # 5.27

# Fisher's test
fisher.test(table, alternative="greater")
# Fisher's Exact Test for Count Data
# 
# data:  table
# p-value < 2.2e-16
# alternative hypothesis: true odds ratio is greater than 1
# 95 percent confidence interval:
# 4.866973      Inf
# sample estimates:
# odds ratio 
# 5.274318 

# To get exact p-value
result2 <- fisher.test(table, alternative="greater")
result2$p.value
# [1] 0