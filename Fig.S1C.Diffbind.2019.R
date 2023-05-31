

setwd("~/Dropbox/Brown_Webb/Desktop/Ren_ATAC-seq")
library(DiffBind)
library(rtracklayer)
library(rgl)

# ======================================================================
# ATAC-seq from the Ren lab
# Comparison of aNSCs ATAC-seq with other tissues
# DiffBind
# November 2019
# R version 3.3.1
# ======================================================================

# Load data
all <- dba(sampleSheet="Ren.ATAC.tissues.csv")
all
# 14 Samples, 204886 sites in matrix (328626 total):
          #   ID    Tissue    Factor     Condition  Treatment Replicate Caller Intervals
# 1  intestine_1 intestine Chromatin Proliferating Full-Media         1 narrow    100719
# 2  intestine_2 intestine Chromatin Proliferating Full-Media         2 narrow     85748
# 3     kidney_1    kidney Chromatin Proliferating Full-Media         1 narrow    100922
# 4     kidney_2    kidney Chromatin Proliferating Full-Media         2 narrow     74991
# 5      liver_1     liver Chromatin Proliferating Full-Media         1 narrow    139683
# 6      liver_2     liver Chromatin Proliferating Full-Media         2 narrow    109099
# 7       lung_1      lung Chromatin Proliferating Full-Media         1 narrow    159871
# 8       lung_2      lung Chromatin Proliferating Full-Media         2 narrow    155130
# 9    stomach_1   stomach Chromatin Proliferating Full-Media         1 narrow    118933
# 10   stomach_2   stomach Chromatin Proliferating Full-Media         2 narrow    123387
# 11      aNSC_1      NSCs Chromatin Proliferating Full-Media         1 narrow     52503
# 12      aNSC_2      NSCs Chromatin Proliferating Full-Media         2 narrow     41480
# 13   qNSC_rep1      NSCs Chromatin     Quiescent BMP4-Media         1 narrow     47056
# 14   qNSC_rep2      NSCs Chromatin     Quiescent BMP4-Media         2 narrow     52076


# Count reads
# all <- dba.count(all, minOverlap=2)

# Save read count DB object
# savefile <- dba.save(all, 'allReads')

# Load read count
all <- dba.load('allReads')
all
# 14 Samples, 204886 sites in matrix:
          #   ID    Tissue    Factor     Condition  Treatment Replicate Caller Intervals FRiP
# 1  intestine_1 intestine Chromatin Proliferating Full-Media         1 counts    204886 0.25
# 2  intestine_2 intestine Chromatin Proliferating Full-Media         2 counts    204886 0.22
# 3     kidney_1    kidney Chromatin Proliferating Full-Media         1 counts    204886 0.21
# 4     kidney_2    kidney Chromatin Proliferating Full-Media         2 counts    204886 0.17
# 5      liver_1     liver Chromatin Proliferating Full-Media         1 counts    204886 0.46
# 6      liver_2     liver Chromatin Proliferating Full-Media         2 counts    204886 0.37
# 7       lung_1      lung Chromatin Proliferating Full-Media         1 counts    204886 0.33
# 8       lung_2      lung Chromatin Proliferating Full-Media         2 counts    204886 0.37
# 9    stomach_1   stomach Chromatin Proliferating Full-Media         1 counts    204886 0.20
# 10   stomach_2   stomach Chromatin Proliferating Full-Media         2 counts    204886 0.24
# 11      aNSC_1      NSCs Chromatin Proliferating Full-Media         1 counts    204886 0.12
# 12      aNSC_2      NSCs Chromatin Proliferating Full-Media         2 counts    204886 0.11
# 13   qNSC_rep1      NSCs Chromatin     Quiescent BMP4-Media         1 counts    204886 0.11
# 14   qNSC_rep2      NSCs Chromatin     Quiescent BMP4-Media         2 counts    204886 0.10

# Get pearson's correlation among ATAC-seq
values <- dba.plotHeatmap(all)
values

# =======================
# Without qNSCs
# =======================
            # aNSC_1 aNSC_2 lung_2 lung_1 liver_2 liver_1 intestine_2 intestine_1 kidney_2 kidney_1 stomach_1 stomach_2
# aNSC_1        1.00   0.85   0.56   0.57    0.46    0.45        0.49        0.49     0.60     0.61      0.54      0.61
# aNSC_2        0.85   1.00   0.57   0.58    0.46    0.46        0.50        0.50     0.60     0.61      0.56      0.61
# lung_2        0.56   0.57   1.00   0.95    0.49    0.50        0.51        0.52     0.62     0.65      0.65      0.72
# lung_1        0.57   0.58   0.95   1.00    0.48    0.49        0.51        0.52     0.63     0.66      0.64      0.72
# liver_2       0.46   0.46   0.49   0.48    1.00    0.92        0.56        0.58     0.59     0.58      0.60      0.59
# liver_1       0.45   0.46   0.50   0.49    0.92    1.00        0.57        0.59     0.59     0.58      0.61      0.61
# intestine_2   0.49   0.50   0.51   0.51    0.56    0.57        1.00        0.91     0.66     0.66      0.74      0.72
# intestine_1   0.49   0.50   0.52   0.52    0.58    0.59        0.91        1.00     0.67     0.67      0.77      0.74
# kidney_2      0.60   0.60   0.62   0.63    0.59    0.59        0.66        0.67     1.00     0.87      0.70      0.73
# kidney_1      0.61   0.61   0.65   0.66    0.58    0.58        0.66        0.67     0.87     1.00      0.71      0.74
# stomach_1     0.54   0.56   0.65   0.64    0.60    0.61        0.74        0.77     0.70     0.71      1.00      0.89
# stomach_2     0.61   0.61   0.72   0.72    0.59    0.61        0.72        0.74     0.73     0.74      0.89      1.00

# =======================
# With qNSCs
# =======================

            # lung_2 lung_1 aNSC_1 aNSC_2 qNSC_rep1 qNSC_rep2 liver_2 liver_1 intestine_2 intestine_1 kidney_2 kidney_1 stomach_1 stomach_2
# lung_2        1.00   0.95   0.55   0.55      0.54      0.55    0.50    0.52        0.52        0.53     0.63     0.65      0.65      0.73
# lung_1        0.95   1.00   0.56   0.56      0.55      0.56    0.49    0.50        0.52        0.53     0.64     0.66      0.65      0.73
# aNSC_1        0.55   0.56   1.00   0.85      0.82      0.84    0.44    0.44        0.47        0.47     0.59     0.60      0.53      0.59
# aNSC_2        0.55   0.56   0.85   1.00      0.81      0.83    0.45    0.45        0.49        0.49     0.59     0.60      0.55      0.60
# qNSC_rep1     0.54   0.55   0.82   0.81      1.00      0.88    0.44    0.44        0.48        0.47     0.59     0.60      0.54      0.59
# qNSC_rep2     0.55   0.56   0.84   0.83      0.88      1.00    0.46    0.45        0.49        0.49     0.60     0.61      0.55      0.60
# liver_2       0.50   0.49   0.44   0.45      0.44      0.46    1.00    0.92        0.57        0.59     0.59     0.58      0.60      0.60
# liver_1       0.52   0.50   0.44   0.45      0.44      0.45    0.92    1.00        0.58        0.60     0.59     0.59      0.62      0.62
# intestine_2   0.52   0.52   0.47   0.49      0.48      0.49    0.57    0.58        1.00        0.91     0.66     0.66      0.75      0.72
# intestine_1   0.53   0.53   0.47   0.49      0.47      0.49    0.59    0.60        0.91        1.00     0.67     0.67      0.77      0.75
# kidney_2      0.63   0.64   0.59   0.59      0.59      0.60    0.59    0.59        0.66        0.67     1.00     0.87      0.70      0.73
# kidney_1      0.65   0.66   0.60   0.60      0.60      0.61    0.58    0.59        0.66        0.67     0.87     1.00      0.71      0.75
# stomach_1     0.65   0.65   0.53   0.55      0.54      0.55    0.60    0.62        0.75        0.77     0.70     0.71      1.00      0.89
# stomach_2     0.73   0.73   0.59   0.60      0.59      0.60    0.60    0.62        0.72        0.75     0.73     0.75      0.89      1.00