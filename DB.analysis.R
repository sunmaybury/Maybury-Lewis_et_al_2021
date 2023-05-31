
# ========================================================================
# Taking differential analysis from DiffBind and generating DBsite lists
# To divide the DBsites into proximal and distal for gene assignment (GREAT)
# and compare with diffentially expressed gene lists
# Includes postnatal quiescent, activated samples
# Not downsampled, without including reactivated for normalizing read counts
# October 2019
# =========================================================================

setwd("~/Dropbox/Brown_Webb/Desktop/ATAC-seq-PEAKS/DiffBind/results.diffbind.2019")

# Create a column that numbers each DBsite and write output to a bed file
# All csv files are diffbind results FDR < 0.05

# Load Activated vs Quiescent
df <- read.csv("v2_AvsQ/AvsQ.DB.analysis.csv")

# Add column with row numbers
new.df <- cbind(df, 1:nrow(df))
new.df <- `new.df`[, c(1,2,3,12,9,11)]
colnames(`new.df`) <- c("chr", "start", "end", "DBsite", "Fold", "FDR")

`new.df`$DBsite <- sub("^", "DBsite", `new.df`$DBsite)

write.table(`new.df`, file="AvsQ.DBsites.bed", sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)

# Separate AA and AQ DBsites
AA <- `new.df`[`new.df`$Fold>0, ]
AQ <- `new.df`[`new.df`$Fold<0, ]

write.table(AA, file="AA.DBsites.bed", sep="\t", quote=F, row.names=F, col.names=F)
write.table(AQ, file="AQ.DBsites.bed", sep="\t", quote=F, row.names=F, col.names=F)



