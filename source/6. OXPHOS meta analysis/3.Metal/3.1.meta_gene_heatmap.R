rm(list = ls())
setwd("D:/R.data/DHFR_project/6. OXPHOS meta analysis")
library(data.table)
library(limma)
library(dplyr)
library(ggplot2)

#
load("OXPHOS_meta_gene.rdata")

#
dat1$t = dat1$Effect/dat1$StdErr
oxphos = fread("oxphos_genes.txt", header = F)
oxphos = data.frame(oxphos)
rownames(dat1) = dat1$MarkerName
idx = intersect(dat1$MarkerName, oxphos$V1)
dat1 = dat1[idx,]
dat1 = dat1[oxphos$V1,]

#
df = data.frame(Gene = dat1$MarkerName, Meta = dat1$t)

#
candidate = c("GSE17601_deg.txt", "GSE9479_deg.txt")

#
for(i in 1:length(candidate)){
  de_summary = read.table(candidate[i], header = T)  
  de_summary$Symbol = rownames(de_summary)
  idx1 = match(df$Gene, de_summary$Symbol)
  df[, i+2] = de_summary$t[idx1]
  dataset_name = gsub("_deg.txt", "", candidate[i])
  colnames(df)[i+2] = dataset_name
}

dat = read.table("GSE113113_deg.txt", header = T)
idx1 = match(df$Gene, dat$X)
dat = dat[idx1,]
dat$t = dat$log2FoldChange/dat$lfcSE
df$GSE113113 = dat$t

#
df = data.frame(df, row.names = 1)

df = na.omit(df)
#
library(plot.matrix)
num_group = 20
cut_off = seq(-2, 5, 1)
cut_off = c(min(df),cut_off, max(df))

#
pal_neg = colorRampPalette(c("#80C4C8", "#e9f2f2"))
pal_pos = colorRampPalette(c("#f3dde8", "#fb4e53"))

#
idx1 = which(cut_off < 0)
idx1 = length(idx1)

#
col_neg = pal_neg(idx1)
col_pos = pal_pos(length(cut_off) - idx1 - 1)
col = c(col_neg, col_pos)

df = t(df)

# Open WMF graphics device
win.metafile("F.meta_OXPHOS.wmf", width = 20, height = 5)
par(mar = c(9, 12, 4, 5))
plot(df, 
     breaks = cut_off, 
     col = col,
     axis.col = list(side = 1, las = 2, cex.axis = 1.8),
     axis.row = list(side = 2, las = 1, cex.axis = 1.8),
     cex.main = 1.8,
     main = "", 
     xlab = "",
     ylab = "",
     border = T)
dev.off()