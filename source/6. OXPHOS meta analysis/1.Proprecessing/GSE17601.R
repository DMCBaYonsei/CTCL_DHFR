rm(list=ls())
options(stringsAsFactors = F)
setwd("D:/R.data/OXPHOS")

library(data.table)
library(limma)
library(plyr)
library(WGCNA)

#load dataset
dat1 <- fread("D:/R.data/GEO/GSE17601_series_matrix.txt", fill = T)
dat1 <- data.frame(dat1, stringsAsFactors = F)
dat1 <- dat1[-nrow(dat1), ]

sample <- fread("GSE17601_DHFR_Phenotype.txt", fill = T)
sample = t(sample)
sample = data.frame(sample)
sample = rbind(colnames(sample), sample)
sample = unlist(sample$X3)
table(sample)


#select expression data
dat1 <- dat1[-(1:65), ]
colnames(dat1)[1] <- "ID"

#reorganize group
idx.h = which(sample == "HIGH")
idx.l = which(sample == "LOW")

dat1 <- dat1[, c(1, idx.l, idx.h)]
write.table(dat1, "temp.txt", col.names = F, row.names = F, sep = "\t")

#read data again with header
dat1 = fread("temp.txt", header = T, stringsAsFactors = F)
dat1 = data.frame(dat1)

# scale and normalize
dat1[, -1] <- normalizeQuantiles(as.matrix(dat1[, -1]))
boxplot(dat1[sample(1:nrow(dat1), 2000), -1])
rownames(dat1)=dat1$ID_REF
dat1 = dat1[,-1]

#ID to gene symbol
source("D:/R.data/GPL/GPL96.R")

#probe selection
rownames(dat1) = paste0("s", (1:nrow(dat1)))
dat2 = collapseRows(dat1[, -1], rowID = rownames(dat1), rowGroup = dat1[, 1])
datExpr = dat2$datETcollapsed
dx = paste0("High_", c(1:length(idx.h)))
cn = paste0("Low_", c(1:length(idx.l)))
colnames(datExpr) = c(cn, dx)
datExpr = data.frame(datExpr)

#select OXPHOS gene
datExpr$symbol = rownames(datExpr)
gene_list <- fread("Oxphos_gene_list.txt", header = T, fill = T)
idx = intersect(datExpr$symbol,gene_list$Symbol)
dat3 = datExpr[idx,]
dat3 = dat3[,-17]

#export normalization data
write.table(dat3, "GSE17601_oxphos.txt", sep = "\t")
