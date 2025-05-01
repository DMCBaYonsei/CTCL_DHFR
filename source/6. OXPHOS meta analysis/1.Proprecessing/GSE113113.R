rm(list=ls())
options(stringsAsFactors = F)
setwd("D:/R.data/OXPHOS")

library(data.table)
library(limma)
library(plyr)
library(WGCNA)


#Load dataset
expr <- fread("D:/R.data/GEO/GSE113113_data_final.txt", header = T, fill = T)
expr <- data.frame(expr, stringsAsFactors = F)

sample <- fread("GSE113113_DHFR_Phenotype.txt", fill = T)
sample = t(sample)
sample = data.frame(sample)
sample = rbind(colnames(sample), sample)
sample = unlist(sample$X3)


#reorganize group
idx.h = which(sample == "HIGH")
idx.l = which(sample == "LOW")

expr <- expr[, c(1, idx.h, idx.l)]

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

#
cn = paste0("High_", c(1:length(idx.h)))
dx = paste0("Low_", c(1:length(idx.l)))
colnames(expr_final) = c(cn, dx)

#select OXPHOS gene
expr_final$symbol = rownames(expr_final)
gene_list <- fread("Oxphos_gene_list.txt", header = T, fill = T)
idx = intersect(expr_final$symbol,gene_list$Symbol)
dat3 = expr_final[idx,]
dat3 = dat3[,-26]

#export count data
write.table(dat3, "GSE113113_OXPHOS.txt", sep = "\t")
