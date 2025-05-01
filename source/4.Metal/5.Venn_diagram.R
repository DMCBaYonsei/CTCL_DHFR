rm(list = ls())
setwd("D:/R.data/DHFR_project/data")
library(data.table)
library(limma)
library(ggplot2)
library(ggvenn)

#load GSE113113 result
dat1 = fread("GSE113113_deg2.txt")
dat1 = data.frame(dat1, stringsAsFactors = F)
rownames(dat1) = dat1[, 1]
idx1 = which(dat1$sig == "Up regulation")
dat1 = dat1[idx1,]
GSE113113 = dat1[,1]

#load GSE59307 result
dat1 = fread("GSE59307_deg2.txt")
dat1 = data.frame(dat1, stringsAsFactors = F)
rownames(dat1) = dat1[, 1]
idx1 = which(dat1$sig == "Up regulation")
dat1 = dat1[idx1,]
GSE59307 = dat1[,1]

#load GSE9479 result
dat1 = fread("GSE9479_deg2.txt")
dat1 = data.frame(dat1, stringsAsFactors = F)
rownames(dat1) = dat1[, 1]
idx1 = which(dat1$sig == "Up regulation")
dat1 = dat1[idx1,]
GSE9479 = dat1[,1]

#Load meta result
dat1 = fread("meta_sign.txt")
dat1 = data.frame(dat1, stringsAsFactors = F)
rownames(dat1) = dat1[, 1]
idx1 = which(dat1$sig == "Up regulation")
dat1 = dat1[idx1,]
Meta = dat1[,1]


#Combine all upregulated gene lists into a named list
deg_list = list(
  Meta = Meta,
  GSE113113 = GSE113113,
  GSE59307 = GSE59307,
  GSE9479 = GSE9479)

#Plot Venn diagram to visualize overlaps in upregulated genes
plot <- ggvenn(
  deg_list, 
  fill_color = c("#fb4e53", "#f3dde8", "#EAFFD0", "#0183b1"),
  stroke_size = 0.5, 
  set_name_size = 6,
  text_size = 5.5)

# Copy current plot to a WMF file
dev.copy(win.metafile, filename = "Venn_diagram.wmf", width = 9, height = 9)
dev.off()

#Find intersection of upregulated genes across all datasets
common = intersect(GSE113113, GSE59307)
common = intersect(common, GSE9479)
common = intersect(common, Meta)

#Convert the result into a data frame and save to file
common = as.data.frame(common)
write.table(common, file = "common.txt", sep = "\t")
