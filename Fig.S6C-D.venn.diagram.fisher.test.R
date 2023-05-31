# ===============================================
# Make weighted venn diagrams
# AA sites that are also Ascl1 binding sites
# AQ sites that are also NFI binding sites
# Fisher's exact test
# v2 DiffBind AvsQ 9929 DBsites
# November 2019
# ===============================================

setwd("~/Dropbox/Desktop/ATAC-seq-PEAKS/AvsQ_analysis_2019_v2/TF_ATAC-seq")

library(venneuler)

# Counting Ascl1 binding sites with dynamic chromatin 
AA_Ascl1 <- venneuler(c(AA=750, dynamic_Ascl1=951, "AA&dynamic_Ascl1"=6027))    
plot(AA_Ascl1)

# counting NFI binding sites with dynamic chromatin
AQ_NFI <- venneuler(c(AQ=1941, dynamic_NFI=1840, "AQ&dynamic_NFI"=1211))
plot(AQ_NFI)

AA_NFI <- venneuler(c(AA=4937, dynamic_NFI=1211, "AA&dynamic_NFI"=1840))
plot(AA_NFI)

# ==============================================================

# AA sites = 6777
# AQ sites = 3152
# (dynamic sites = AA+ AQ)
# Ascl1-binding dynamic sites = 6978
# NFI-binding dynamic sites = 3051

# ==============================================================
# Testing whether dynamic sites are enriched with TF binding
# ==============================================================

### Ascl1

AA_Ascl1 <- matrix(c(6027,750,951,2201), nrow=2, dimnames = list(c("Ascl1", "notAscl1"), c("AA", "notAA")))
           # AA notAA
# Ascl1    6027   951
# notAscl1  750  2201

# Odds ratio:
odds_Ascl1 <- 6027 / 951
odds_not_Ascl1 <- 750 / 2201
odds_ratio <- odds_Ascl1 / odds_not_Ascl1     # 18.60

# Fisher's test:
fisher.test(AA_Ascl1, alternative="greater")
# Fisher's Exact Test for Count Data
# 
# data:  AA_Ascl1
# p-value < 2.2e-16
# alternative hypothesis: true odds ratio is greater than 1
# 95 percent confidence interval:
# 16.96492      Inf
# sample estimates:
# odds ratio 
# 18.58983 

# To get the exact p-value:
result2 <- fisher.test(AA_Ascl1, alternative="greater")
result2$p.value
# [1] 0

### NFI

AQ_NFI <- matrix(c(1211,1941,1840,4937), nrow=2, dimnames = list(c("NFI", "notNFI"), c("AQ", "notAQ")))
         # AQ notAQ
# NFI    1211  1840
# notNFI 1941  4937

# Odds ratio:
odds_NFI <- 1211 / 1840
odds_not_NFI <- 1941 / 4937
odds_ratio <- odds_NFI / odds_not_NFI    # 1.67

# Fisher's test:
fisher.test(AQ_NFI, alternative="greater")
# Fisher's Exact Test for Count Data
# 
# data:  AQ_NFI
# p-value < 2.2e-16
# alternative hypothesis: true odds ratio is greater than 1
# 95 percent confidence interval:
# 1.551128      Inf
# sample estimates:
# odds ratio 
# 1.673998 

# To get the exact p-value:
result2 <- fisher.test(AQ_NFI, alternative="greater")
result2$p.value
# [1] 2.120805e-29

AA_NFI <- matrix(c(1840,4937,1211,1941), nrow=2, dimnames = list(c("NFI", "notNFI"), c("AA", "notAA")))
         # AA notAA
# NFI    1840  1211
# notNFI 4937  1941

# Odds ratio:
odds_NFI <- 1840 / 1211
odds_not_NFI <- 4937 / 1941
odds_ratio <- odds_NFI / odds_not_NFI    # 0.59

# Fisher's test:
fisher.test(AA_NFI, alternative="greater")
# Fisher's Exact Test for Count Data
# 
# data:  AA_NFI
# p-value = 1
# alternative hypothesis: true odds ratio is greater than 1
# 95 percent confidence interval:
# 0.5536121       Inf
# sample estimates:
# odds ratio 
# 0.5973724 

# To get the exact p-value:
result2 <- fisher.test(AA_NFI, alternative="greater")
result2$p.value
