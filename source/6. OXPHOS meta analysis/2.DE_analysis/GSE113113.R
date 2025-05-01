rm(list=ls())
options(stringsAsFactors = F)
setwd("D:/R.data/OXPHOS")

library(data.table)
library(limma)
library(plyr)
library(dplyr)
library(DESeq2)
library(ggplot2)
library(ggrepel)


#load count data
expr <- fread("GSE113113_OXPHOS.txt")
expr <- data.frame(expr, row.names = 1)

Disease = c(rep("high", 13), rep("low", 12))
sample = colnames(expr)
sample = data.frame(sample, Disease)
sample$Disease = as.factor(sample$Disease)
levels(sample$Disease)
sample$Disease = relevel(sample$Disease, ref = "low")

#Differential expression analysis with DESeq2
dds = DESeqDataSetFromMatrix(countData = expr,
                             colData = sample,
                             design = ~ Disease)

dds <- DESeq(dds)

results <- results(dds)


write.csv(as.data.frame(results), file = "DESeq2_results.csv")

result = read.csv("DESeq2_results.csv")

#export result
write.table(result, "GSE113113_deg.txt", sep = "\t", row.names = F)


#export data for metal
DEG <- data.frame(Symbol = result$X,
                  log2FC = result$log2FoldChange,
                  SE = result$lfcSE)
write.table(DEG, "GSE113113_deg_2.txt", row.names = F, sep = "\t")
