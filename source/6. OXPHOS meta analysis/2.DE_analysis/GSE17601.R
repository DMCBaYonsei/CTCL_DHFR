rm(list = ls())
setwd("D:/R.data/OXPHOS")
library(data.table)
library(limma)
library(ggplot2)
library(dplyr)
library(ggrepel)

#load normalization data
data = "GSE17601"
storage = paste0(data, "_OXPHOS.txt")
dat1 = fread(storage)
dat1 = data.frame(dat1, stringsAsFactors = F)
rownames(dat1) = dat1[, 1]
dat1 = dat1[, -1]

#Differential expression analysis with limma
Disease = c(rep(0, 8), rep(1, 8))
design = model.matrix(~ Disease)
fit = lmFit(dat1, design)
fit = eBayes(fit, trend = TRUE)
result = topTable(fit, coef=2, n = nrow(dat1), adjust="fdr", p = 1)
result = result[rownames(dat1), ]

#export result
storage = paste0(data, "_deg.txt")
write.table(result, storage, sep = "\t")


#export data for metal
result$Symbol = rownames(result)
DEG <- data.frame(Symbol = result$Symbol,
                  log2FC = result$logFC,
                  SE = result$logFC/result$t)
storage = paste0(data, "_deg_2.txt")
write.table(DEG, storage, row.names = FALSE, sep = "\t")
