rm(list = ls())
setwd("D:/R.data/DHFR_project/data")
library(data.table)
library(limma)
library(dplyr)
library(ggplot2)

#Load the meta-analysis result data
load("lymphoma_meta_gene.rdata")

#Calculate t-statistics for each gene
dat1$t = dat1$Effect/dat1$StdErr

#Sort the data by decreasing absolute t-value
dat1 = dat1[order(-1*dat1$t), ]

#Determine the cutoff t-value corresponding to the top 150-200 genes
cut = abs(dat1$t)
cut = cut[order(-1*cut)]
cut_2 = cut[150]
cut_1 = cut[201]

# Filter genes with absolute t-value above the cutoff
dat1 = dat1[(abs(dat1$t) > cut_1) & (abs(dat1$t) < cut_2), ]
dat1 = dat1[order(-1*dat1$t), ]

# Create a dataframe with gene names and meta-analysis t-values
df = data.frame(Gene = dat1$MarkerName, Meta = dat1$t)

# Define DEG summary files to extract individual dataset t-values
candidate = c("GSE59307_deg.txt", "GSE9479_deg.txt")

# For each DEG dataset, extract t-values corresponding to the top meta genes
for(i in 1:length(candidate)){
  de_summary = read.table(candidate[i], header = T)  
  de_summary$Symbol = rownames(de_summary)
  idx1 = match(df$Gene, de_summary$Symbol)
  df[, i+2] = de_summary$t[idx1]
  dataset_name = gsub("_deg.txt", "", candidate[i])
  colnames(df)[i+2] = dataset_name
}

# Read in a third DEG dataset with slightly different column names
dat = read.table("GSE113113_deg.txt", header = T)
idx1 = match(df$Gene, dat$X)
dat = dat[idx1,]
dat$t = dat$log2FoldChange/dat$lfcSE
df$GSE113113 = dat$t

# Set row names to genes and remove gene name column
df = data.frame(df, row.names = 1)

# Load plot.matrix for heatmap-style visualization
library(plot.matrix)

# Define number of color groups and breakpoints for coloring
num_group = 20
cut_off = seq(-4, 12, 2)
cut_off = c(min(df), cut_off, max(df))

# Create color palettes for negative and positive values
pal_neg = colorRampPalette(c("#0183b1", "#e9f2f2"))
pal_pos = colorRampPalette(c("#f3dde8", "#fb4e53"))

# Determine how many colors are needed on each side of zero
idx1 = which(cut_off < 0)
idx1 = length(idx1)

# Generate the actual color vectors
col_neg = pal_neg(idx1)
col_pos = pal_pos(length(cut_off) - idx1 - 1)
col = c(col_neg, col_pos)

# Transpose the dataframe for plotting
df = t(df)

# Save the heatmap to a PDF file
win.metafile("F.meta_lymphoma_genes_4.wmf", width = 20, height = 5.5)
par(mar = c(11, 10, 5, 5))

# Plot the matrix
plot(df, 
     breaks = cut_off, 
     col = col,
     axis.col = list(side = 1, las = 2, cex.axis = 1.8),
     axis.row = list(side = 2, las = 1, cex.axis = 1.8),
     cex.main = 1.8,
     main = "Top 151 - 200 CTCL_meta_genes", 
     xlab = "",
     ylab = "",
     border = T)
dev.off()
