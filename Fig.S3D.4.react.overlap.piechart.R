# =====================================================
# Make a pichart to show breakdown of overlaps
# between Reactivated consensus peaks and
# AQ, AA, stable, and react-specific sites
# v2 DBsites 9929 normalized read counts b/w AvsQ
# Reactivated downsampled, MACS q 0.0001
# November 2019
# =====================================================


setwd("~/Dropbox/Desktop/ATAC-seq-PEAKS/AvsQ_analysis_2019_v2/Reactivated")

# ==================================================
# Part 1 - all Reactivated consPeaks
# ==================================================

library(ggplot2)

df <- data.frame(
  group=c("AQ", "AA", "Stable", "Unique"),
  value=c(1393, 6466, 19766, 24314)
)

head(df)

bp <- ggplot(df, aes(x="", y=value, fill=group)) +
  geom_bar(width=1, stat="identity")

pie <- bp + coord_polar("y", start=0)
pie

ggsave("React.consPeaks.all.pie.pdf", width=5, height=5)

# ==================================================
# Part 2 - Only overlapping  Reactivated consPeaks
# ==================================================

library(ggplot2)

df <- data.frame(
  group=c("AQ", "AA", "Stable"),
  value=c(1393, 6466, 19766)
)

head(df)

bp <- ggplot(df, aes(x="", y=value, fill=group)) +
  geom_bar(width=1, stat="identity")

pie <- bp + coord_polar("y", start=0)
pie

ggsave("React.consPeaks.onlyAvsQ.overlap.pie.pdf", width=5, height=5)
