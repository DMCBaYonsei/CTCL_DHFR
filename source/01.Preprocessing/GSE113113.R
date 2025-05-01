rm(list=ls())
options(stringsAsFactors = F)
setwd("D:/R.data/DHFR_project/data")

library(data.table)
library(limma)
library(plyr)
library(WGCNA)

#Load dataset
expr <- fread("D:/R.data/GEO/GSE113113_data_final.txt", header = T, fill = T)
expr <- data.frame(expr, stringsAsFactors = F)
sample <- fread("D:/R.data/GEO/GSE113113_series_matrix.txt", fill = T)
sample = data.frame(t(sample[c(30, 39, 41), ]))
colnames(sample) = c("id", "dx", "stage")
idx1 = match(colnames(expr[,-1]), sample[-1,]$id)
sample = sample[-1,]
sample = sample[idx1,]
write.table(sample, "GSE113113_sample.txt", sep = "\t", row.names = F)

#exp data
expr[,-1] = exp(expr[,-1])

#gene selection
expr = data.frame(expr)
rownames(expr) = paste0("s", (1:nrow(expr)))
gene_row_count = table(expr[,1])
gene_row_count = data.frame(gene_row_count)

#genes only one value
idx1 = which(gene_row_count$Freq ==1)
gene_conserved = as.character(gene_row_count$Var1[idx1])

#device one value gene and duplication gene
idx1 = which(expr$REFID %in% gene_conserved)
expr_cons = expr[idx1,]
expr_dup = expr[-idx1,]

#duplicated gene selection
expr = expr_dup
rownames(expr) = paste0("s", (1:nrow(expr)))
dat2 = collapseRows(expr[, -1], rowID = rownames(expr), rowGroup = expr[, 1])
datExpr = dat2$datETcollapsed
expr_cons = data.frame(expr_cons, row.names = 1)

#combine final data expression
expr_final = rbind(expr_cons, datExpr)

#tranform to gene count
expr_final = round(expr_final, 0)

#export count data
write.table(expr_final, "GSE113113_count.txt", sep = "\t")
