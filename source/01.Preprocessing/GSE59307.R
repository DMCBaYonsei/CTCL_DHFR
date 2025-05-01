rm(list=ls())
options(stringsAsFactors = F)
setwd("D:/R.data/DHFR_project/data")

library(data.table)
library(limma)
library(plyr)
library(WGCNA)

#load dataset
dat1 <- fread("D:/R.data/GEO/GSE59307_series_matrix.txt", fill = T)
dat1 <- data.frame(dat1, stringsAsFactors = F)
dat1 <- dat1[-nrow(dat1), ]

#select and export sample information
dat2 = data.frame(t(dat1[c(31, 42, 39), -1]))
colnames(dat2) = c("id", "tissue", "dx")
sample_storage = paste0("GSE59307_sample.txt")
write.table(dat2, sample_storage, sep = "\t", row.names = F)

#select compare condition
sample = unlist(dat1[39, ])
table(sample)

#select expression data
dat1 <- dat1[-(1:62), ]
colnames(dat1)[1] <- "ID"

#reorganize group
idx.c = which(sample == "condition: Healthy")
idx.d = which(sample == "condition: Cutaneus Tcell Lymphoma")

dat1 <- dat1[, c(1, idx.c, idx.d)]
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
source("D:/R.data/GPL/GPL570.R")

#probe selection
rownames(dat1) = paste0("s", (1:nrow(dat1)))
dat2 = collapseRows(dat1[, -1], rowID = rownames(dat1), rowGroup = dat1[, 1])
datExpr = dat2$datETcollapsed
cn = paste0("Normal_", c(1:length(idx.c)))
dx = paste0("CTCL_", c(1:length(idx.d)))
colnames(datExpr) = c(cn, dx)

#export normalization data
write.table(datExpr, "GSE59307_norm.txt", sep = "\t")
