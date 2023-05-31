# ========================================================
# Venn diagrams for H3K4me3 in qNSCs and aNSCs
# that overlap with AvsQ.truStable sites
# v2 AvsQ DiffBind 9929 DBsites
# 19976 AvsQ truStable sites
# October-November 2019
# ========================================================

setwd("~/Dropbox/Brown_Webb/Desktop/ATAC-seq-PEAKS/AvsQ_analysis_2019_v2/H3K4me3")


library(venneuler)

# quiescent H3K4me3 and truStable sites
quies <- venneuler(c(H3K4me3=6247, stable=9037, "stable&H3K4me3"=10939))
plot(quies)

# activated H3K4me3 and truStable sites
act <- venneuler(c(H3K4me3=5357, stable=9353,"stable&H3K4me3"=10623))
plot(act)



# ==================================================
# Fisher's exact test
# Are AvsQ truStable sites enriched with H3K4me3?
# ==================================================

# ====================
# Quiescent H3K4me3
# ====================

# Not stable + Not H3K4me3 = DBsites - DB.quies.H3K4me3

quies_test <- matrix(c(10939,9037,6247,7868), nrow=2,
                     dimnames = list(c("H3K4me3", "notH3K4me3"), c("stable", "notStable")))
           # stable notStable
# H3K4me3     10939      6247
# notH3K4me3   9037      7868

# Odds ratio
odds_H3K4me3 <- 10939 / 6247 
odds_notH3K4me3 <- 9037 / 7868
odds_ratio <- odds_H3K4me3/odds_notH3K4me3    #1.52

# Fisher's test
fisher.test(quies_test, alternative="greater")
# Fisher's Exact Test for Count Data
# 
# data:  quies_test
# p-value < 2.2e-16
# alternative hypothesis: true odds ratio is greater than 1
# 95 percent confidence interval:
# 1.469741      Inf
# sample estimates:
# odds ratio 
# 1.524509 

# To get the exact p-value:
result2 <- fisher.test(quies_test, alternative="greater")
result2$p.value
# [1] 1.03646e-81



# ====================
# Activated H3K4me3
# ====================

# Not stable + Not H3K4me3 = DBsites - DB.act.H3K4me3

act_test <- matrix(c(10623,9353,5357,8053), nrow=2,
                   dimnames = list(c("H3K4me3", "notH3K4me3"), c("stable", "notStable")))
         
           # stable notStable
# H3K4me3     10623      5357
# notH3K4me3   9353      8053

# Odds ratio
odds_H3K4me3 <- 10623 / 5357
odds_notH3K4me3 <- 9353 / 8063
odds_ratio <- odds_H3K4me3/odds_notH3K4me3    # 1.71


# Fisher's test
fisher.test(act_test, alternative="greater")
# Fisher's Exact Test for Count Data
# 
# data:  act_test
# p-value < 2.2e-16
# alternative hypothesis: true odds ratio is greater than 1
# 95 percent confidence interval:
# 1.644617      Inf
# sample estimates:
# odds ratio 
# 1.707394 

# To get the exact p-value:
result2 <- fisher.test(act_test, alternative="greater")
result2$p.value
# [1] 3.339426e-125

