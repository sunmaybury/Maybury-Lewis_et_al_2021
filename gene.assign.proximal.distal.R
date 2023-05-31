
# =============================================================
# This script combines the description of gene assignment to 
# reactivated analysis bed files using GREAT,
# and formating the results to bed files for further analysis

# v2 redone with AvsQ DBsites (reactivated to be added later)
# 9929 DBsites
# 19976 truStable sites

# October 2019
# =============================================================

setwd("~/Dropbox/Brown_Webb/Desktop/ATAC-seq-PEAKS/AvsQ_analysis_2019_v2/GREAT_assign")

# First ran GREAT on AvsQ.DBsites and AvsQ.truStable sites
# with -1kb/+1kb to get proximal gene assignments
# Then ran GREAT again with -25kb/+10kb to get proximal + distal assigments

library(tidyr)
library(dplyr)


# =================================================================================
# Part 1
# Append chromosome positions to reactivated consensus sites with assigned genes
# Separate rows with multiple gene associations
# Write out files to be edited by bash
# =================================================================================


# =======================
# AvsQ.DB
# =======================

### Proximal + distal

# Results from GREAT, reactivated consensus site and associated gene
db <- read.table("AvsQ.DBsites.proximal+distal.txt", sep="\t", stringsAsFactors=F)
colnames(db) <- c("DBsite", "gene")
chrpos <- read.table("AvsQ.DBsites.bed", sep="\t", stringsAsFactors=F)
colnames(chrpos) <- c("chr", "start", "end", "DBsite", "Fold", "FDR")
df <- inner_join(chrpos, db, by="DBsite")

# Select only rows that are associated with genes
df_genes <- df[!grepl("NONE", df$gene), ]

# Separate multiple genes that are associated with same DBsite into different rows
df_genes_sep <- df_genes %>%
  mutate(gene = strsplit(as.character(gene), ",")) %>%
  unnest(gene)

write.table(df_genes_sep, file="AvsQ.DBsites.proximal+distal.clean.txt", sep="\t", quote=F, row.names=F, col.names=F)

### Proximal 

db <- read.table("AvsQ.DBsites.proximal.txt", sep="\t", stringsAsFactors=F)
colnames(db) <- c("DBsite", "gene")
df <- inner_join(chrpos, db, by="DBsite")

# Select only rows that are associated with genes
df_genes <- df[!grepl("NONE", df$gene), ]

# Separate multiple genes that are associated with same DBsite into different rows
df_genes_sep <- df_genes %>%
  mutate(gene = strsplit(as.character(gene), ",")) %>%
  unnest(gene)

write.table(df_genes_sep, file="AvsQ.DBsites.proximal.clean.txt", sep="\t", quote=F, row.names=F, col.names=F)


# =======================
# AvsQ.truStable
# =======================

### Proximal + distal

db <- read.table("AvsQ.truStable.proximal+distal.txt", sep="\t", stringsAsFactors=F)
colnames(db) <- c("site", "gene")
chrpos <- read.table("AvsQ.truStable.bed", sep="\t", stringsAsFactors=F)
colnames(chrpos) <- c("chr", "start", "end", "site")
df <- inner_join(chrpos, db, by="site")

# Select only rows that are associated with genes
df_genes <- df[!grepl("NONE", df$gene), ]

# Separate multiple genes that are associated with same DBsite into different rows
df_genes_sep <- df_genes %>%
  mutate(gene = strsplit(as.character(gene), ",")) %>%
  unnest(gene)

write.table(df_genes_sep, file="AvsQ.truStable.proximal+distal.clean.txt", sep="\t", quote=F, row.names=F, col.names=F)

### Proximal

db <- read.table("AvsQ.truStable.proximal.txt", sep="\t", stringsAsFactors=F)
colnames(db) <- c("site", "gene")
df <- inner_join(chrpos, db, by="site")

# Select only rows that are associated with genes
df_genes <- df[!grepl("NONE", df$gene), ]

# Separate multiple genes that are associated with same DBsite into different rows
df_genes_sep <- df_genes %>%
  mutate(gene = strsplit(as.character(gene), ",")) %>%
  unnest(gene)

write.table(df_genes_sep, file="AvsQ.truStable.proximal.clean.txt", sep="\t", quote=F, row.names=F, col.names=F)

# =================================================================================
# Next run modify.GREAT.genelist.1a.sh
# to omit distance to TSS from bed files
# =================================================================================


# =================================================================================
# Part 2
# Take lists of genes that with TSS distances omitted
# Append back to the React.AQ, React.AA, React.truStable, React.unique bed files
# =================================================================================


# ===================
# AvsQ.DBsites
# ===================

df <- read.table("AvsQ.DBsites.proximal+distal.clean.txt", sep="\t", stringsAsFactors=F)
genes <- read.table("AvsQ.DBsites.proximal+distal.genes.txt")
df2 <- cbind(df, genes)
df2$V7 <- NULL
write.table(df2, file="AvsQ.DBsites.proximal+distal.genelist.txt", sep="\t", quote=F, row.names=F, col.names=F)

df <- read.table("AvsQ.DBsites.proximal.clean.txt", sep="\t", stringsAsFactors=F)
genes <- read.table("AvsQ.DBsites.proximal.genes.txt")
df2 <- cbind(df, genes)
df2$V7 <- NULL
write.table(df2, file="AvsQ.DBsites.proximal.genelist.txt", sep="\t", quote=F, row.names=F, col.names=F)

# ===================
# AvsQ.truStable
# ===================

df <- read.table("AvsQ.truStable.proximal+distal.clean.txt", sep="\t", stringsAsFactors=F)
genes <- read.table("AvsQ.truStable.proximal+distal.genes.txt")
df2 <- cbind(df, genes)
df2$V5 <- NULL
write.table(df2, file="AvsQ.truStable.proximal+distal.genelist.txt", sep="\t", quote=F, row.names=F, col.names=F)

df <- read.table("AvsQ.truStable.proximal.clean.txt", sep="\t", stringsAsFactors=F)
genes <- read.table("AvsQ.truStable.proximal.genes.txt")
df2 <- cbind(df, genes)
df2$V5 <- NULL
write.table(df2, file="AvsQ.truStable.proximal.genelist.txt", sep="\t", quote=F, row.names=F, col.names=F)


# =================================================================================
# Subtract Proximal from Proximal + Distal
# to get distal only lists
# =================================================================================

# =============
# AvsQ DBsites
# =============

all <- read.table("AvsQ.DBsites.proximal+distal.genelist.txt", sep="\t", stringsAsFactors=F)
prox <- read.table("AvsQ.DBsites.proximal.genelist.txt", sep="\t", stringsAsFactors=F)
distal <- setdiff(all, prox)
write.table(distal, file="AvsQ.DBsites.distal.genelist.txt", sep="\t", quote=F, row.names=F, col.names=F)

# =============
# AvsQ truStable
# =============

all <- read.table("AvsQ.truStable.proximal+distal.genelist.txt", sep="\t", stringsAsFactors=F)
prox <- read.table("AvsQ.truStable.proximal.genelist.txt", sep="\t", stringsAsFactors=F)
distal <- setdiff(all, prox)
write.table(distal, file="AvsQ.truStable.distal.genelist.txt", sep="\t", quote=F, row.names=F, col.names=F)

