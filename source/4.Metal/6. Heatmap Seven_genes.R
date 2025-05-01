rm(list = ls())
setwd("D:/R.data/DHFR_project/data")
library(data.table)
library(limma)
library(dplyr)
library(ggplot2)
library(ggrepel)

#Load meta result
dat1 = fread("D:/R.data/DHFR_project/Data/Meta_result1.tbl")
dat1 = data.frame(dat1, stringsAsFactors = F)
rownames(dat1) = dat1[, 1]
dat1 = dat1[, -c(2,3)]
dat1$padj = p.adjust(dat1$P.value, method = "fdr")
colnames(dat1) = c("Gene", "logFC", "SE", "pvalue", "direction", "Padj")

idx1 = grep("\\?", dat1[,5])
dat1 = dat1[-idx1,]

idx1 = grep("\\+-", dat1[,5])
dat1 = dat1[-idx1,]
idx1 = which(dat1$direction == "-++" | dat1$direction == "--+")
dat1 = dat1[-idx1,]
# visualize and export significant results
dat1$sig = ifelse(dat1$direction == "+++" & dat1$Padj<0.05, "Up regulation",
                  ifelse(dat1$direction == "---" & dat1$Padj<0.05, "Down regulation", "Not Significant"))
rownames(dat1) = dat1$Gene

#select ONCOin panel
gene_list <- fread("seven_genes.txt", header = T, fill = T)
idx = intersect(dat1$Gene,gene_list$Gene_name_2)
dat1 = dat1[idx,]
dat1 = dat1[gene_list$Gene_name_2,]
#
dat1$t = dat1$logFC/dat1$SE

#
df = data.frame(Gene = dat1$Gene, Meta = dat1$t)

#
candidate = c("GSE59307_deg.txt", "GSE9479_deg.txt")

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

#
library(plot.matrix)
num_group = 20
cut_off = seq(-6, 6, 2)
cut_off = c(cut_off, max(df))

#
pal_neg = colorRampPalette(c("#0183b1", "#e9f2f2"))
pal_pos = colorRampPalette(c("#f3dde8", "#fb4e53"))

#
idx1 = which(cut_off < 0)
idx1 = length(idx1)

#
col_neg = pal_neg(idx1)
col_pos = pal_pos(length(cut_off) - idx1 - 1)
col = c(col_neg, col_pos)

df = t(df)

#
win.metafile("F.seven genes.wmf", width = 8, height = 4)
par(mar = c(6, 9, 1, 3))
plot(df, 
     breaks = cut_off, 
     col = col,
     axis.col = list(side = 1, las = 2, cex.axis = 1.7),
     axis.row = list(side = 2, las = 1, cex.axis = 1.7),
     cex.main = 1.7,
     main = "", 
     xlab = "",
     ylab = "",
     border = T)
dev.off()
