rm(list=ls())
options(stringsAsFactors = F)
setwd("D:/R.data/DHFR_project/data")

library(data.table)
library(limma)
library(plyr)
library(WGCNA)
library(dplyr)

#load dataset
dat1 <- fread("D:/R.data/GEO/GSE9479_series_matrix.txt", fill = T)
dat1 <- data.frame(dat1, stringsAsFactors = F)
dat1 <- dat1[-nrow(dat1), ]

#select compare condition
stage = unlist(dat1[42, ])
table(stage)
stage = as.data.frame(stage)
stage<-stage %>% 
  mutate(stage_2=recode(stage,
                        "Stage IA"="local", 
                        "Stage IB"="local", 
                        "Stage IIB"="advance", 
                        "Stage III"="remove"))

stage = unlist(stage[,2])
#reorganize group
idx.local = which(stage == "local")
idx.advance = which(stage == "advance")

dat1 <- dat1[, c(1, idx.local, idx.advance)]


#select and export sample information
dat2 = data.frame(t(dat1[c(34, 42), -1]))
colnames(dat2) = c("id", "stage")
sample_storage = paste0("GSE9479_sample.txt")
write.table(dat2, sample_storage, sep = "\t", row.names = F)

#select expression data
dat1 <- dat1[-(1:63), ]
colnames(dat1)[1] <- "ID"

write.table(dat1, "temp.txt", col.names = F, row.names = F, sep = "\t")

#read data again with header
dat1 = fread("temp.txt", header = T, stringsAsFactors = F)
dat1 = data.frame(dat1)

# scale and normalize
dat1[,-1] = log10(dat1[,-1])
dat1[, -1] <- normalizeQuantiles(as.matrix(dat1[, -1]))
boxplot(dat1[sample(1:nrow(dat1), 2000), -1])
rownames(dat1)=dat1$ID_REF
dat1 = dat1[,-1]

#ID to gene symbol
source("D:/R.data/GPL/GPL4685.R")

#probe selection
rownames(dat1) = paste0("s", (1:nrow(dat1)))
dat2 = collapseRows(dat1[, -1], rowID = rownames(dat1), rowGroup = dat1[, 1])
datExpr = dat2$datETcollapsed
cn = paste0("Local_", c(1:length(idx.local)))
dx = paste0("Advance_", c(1:length(idx.advance)))
colnames(datExpr) = c(cn, dx)

#export normalization data
write.table(datExpr, "GSE9479_norm.txt", sep = "\t")
