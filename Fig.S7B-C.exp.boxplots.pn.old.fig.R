# =============================================================
## Benchmarking primary activated and quiescent NSCs in culture
## RNA-seq relative gene expression (VST-normalized)
## Figure for Aging Cell Manuscript
## May 2021
## SML
# =============================================================

setwd("~/Dropbox/Brown_Webb/RNA-seq_SML/benchmark/")

library(ggplot2)
library(reshape2)
library(tidyverse)
library(ggpubr)
library(RColorBrewer)
library(extrafont)
loadfonts(device="win")

## Load relative gene expression
vstval <- read.csv("VST-normalized.read.counts.csv")
names(vstval)[1] <- "gene"

## Keep postnatal and old data
vstval <- vstval[,c(1,2,3,4,8,9,10,11,12,13,17,18,19)]

## ======================================================================
## Quiescent NSCs markers
## Fgfr3, Glast, Glt-1, Id2, S100b, Sox9 from Leeman et al
## Glast, Glt-1, Id2, S100b, Nfix, Gfap
# qui.genes <- c("Slc1a3","Slc1a2","Id2","S100b","Vcam1","Gfap")
qui.genes <- c("Dner","Slc1a2","Id2","S100b","Vcam1","Gfap")
qui.expres <- subset(vstval, vstval$gene %in% qui.genes)

## Melt data
df <- melt(qui.expres)
## Separate sample ID into age/celltype and replicates
df2 <- df %>%
    separate(variable, into=c("group","rep"), sep="_(?=[^_]+$)")
## Set groups as factors
df2$group <- factor(df2$group,
                       levels = c('postnat_qNSC','old_qNSC',
                                  'postnat_aNSC','old_aNSC') ,ordered = TRUE)
## Scale expression values 0 to 1
range01 <- function(x){(x-min(x))/(max(x)-min(x))}
df2$normval <- range01(df2$value)

## Boxplot
ggplot(df2, aes(x=group, y=normval, fill=factor(gene))) +
    geom_boxplot() +
    labs(x="NSCs", y="Relative expression", family="Arial") +
    facet_wrap(. ~ gene, ncol=2) +
    theme(panel.grid.major=element_blank(), 
          panel.grid.minor=element_blank(),
          panel.background=element_rect(colour="black", size=0.5, fill=NA),
          strip.background=element_rect(colour="black", fill="white", size=0.5, linetype="solid"),
          strip.text.x=element_text(size=8, color="black", face="bold.italic"),
          axis.text.x = element_text(angle=50, hjust=1, size=8))

## ======================================================================

## ======================================================================
## Activated NSCs markers
## Ascl1, Cdk6, Dlx1, Dlx2, Mcm2 from Leeman et al
## Cdk6, Dlx1, Dlx2, Mcm2
# act.genes <- c("Cdk6","Egfr","Dlx2","Mcm2")
act.genes <- c("Egfr","Sdc1","Mki67","Cdk6","Mcm2","Nes")
act.expres <- subset(vstval, vstval$gene %in% act.genes)

## Melt data
df <- melt(act.expres)
## Separate sample ID into age/celltype and replicates
df2 <- df %>%
    separate(variable, into=c("group","rep"), sep="_(?=[^_]+$)")
## Set groups as factors
df2$group <- factor(df2$group,
                    levels = c('postnat_qNSC','old_qNSC',
                               'postnat_aNSC','old_aNSC') ,ordered = TRUE)
## Scale expression values 0 to 1
df2$normval <- range01(df2$value)

## Boxplot
ggplot(df2, aes(x=group, y=normval, fill=factor(gene))) +
    geom_boxplot() +
    labs(x="NSCs", y="Relative expression", family="Arial") +
    facet_wrap(. ~ gene, ncol=2) +
    theme(panel.grid.major=element_blank(), 
          panel.grid.minor=element_blank(),
          panel.background=element_rect(colour="black", size=0.5, fill=NA),
          strip.background=element_rect(colour="black", fill="white", size=0.5, linetype="solid"),
          strip.text.x=element_text(size=8, color="black", face="bold.italic"),
          axis.text.x = element_text(angle=50, hjust=1, size=8))

