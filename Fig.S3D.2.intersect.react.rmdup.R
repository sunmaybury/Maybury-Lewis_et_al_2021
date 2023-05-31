

setwd("~/Dropbox/Desktop/ATAC-seq-PEAKS/AvsQ_analysis_2019_v2/Reactivated")

library(tidyverse)

# Remove duplicate rows from Reactivated consensus peaks that are AQ
react.AQ <- read.table("React.AQ.bed", sep="\t")
react.AQ <- as_tibble(`react.AQ`)

react.AQ.rmdup <- unique(`react.AQ`)
write.table(`react.AQ.rmdup`, file="React.AQ.rmdup.bed", sep="\t", quote=F, row.names=F, col.names=F)

# Remove duplicate rows from Reactivated consensus peaks that are AA
react.AA <- read.table("React.AA.bed", sep="\t")
react.AA <- as_tibble(`react.AA`)

react.AA.rmdup <- unique(`react.AA`)
write.table(`react.AA.rmdup`, file="React.AA.rmdup.bed", sep="\t", quote=F, row.names=F, col.names=F)

# Remove duplicate rows from Reactivated consensus peaks that are stable between Activated and Quiescent
react.stable <- read.table("React.truStable.bed", sep="\t")
react.stable <- as_tibble(`react.stable`)

react.stable.rmdup <- unique(`react.stable`)
write.table(`react.stable.rmdup`, file="React.truStable.rmdup.bed", sep="\t", quote=F, row.names=F, col.names=F)



