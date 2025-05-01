rm(list = ls())
setwd("D:/R.data/DHFR_project/6. OXPHOS meta analysis")

library(data.table)
library(limma)
library(dplyr)
library(ggplot2)

#load metal meta result
dat1 = fread("Meta_result1.tbl")
dat1 = data.frame(dat1, stringsAsFactors = F)

#calculate adjusted P value
dat1$padj = p.adjust(dat1$P.value, method = "fdr")

#choose significant gene
idx1 = which(dat1$padj < 0.05)
dat1 = dat1[idx1, ]

#save file
save(file = "OXPHOS_meta_gene.rdata", dat1)
