rm(list = ls())
setwd("D:/R.data/DHFR_project/data/")
library(data.table);library(matrixStats);library(limma);library(WGCNA)
library(ggplot2)
library(dplyr)
library(pheatmap)
library(plyr)
library(WGCNA)

#load count data
expr <- fread("GSE113113_count.txt")
expr <- data.frame(expr, row.names = 1)

#load sample data
sample <- fread("GSE113113_sample.txt", fill = T)

#match two data
idx1 = match(colnames(expr), sample[-1,]$id)
sample = sample[-1,]
sample = sample[idx1,]

#select stage
stage = unlist(sample[, 3])
table(stage)

#choose MF
stage = as.data.frame(stage)
stage<-stage %>% 
  mutate(stage_2=recode(stage,
                        "tumor stage: IA"="MF", 
                        "tumor stage: IB"="MF", 
                        "tumor stage: IIA"="MF", 
                        "tumor stage: IIB"="MF", 
                        "tumor stage: IVA2"="SS", 
                        "tumor stage: Normal"="Normal"))

stage = unlist(stage[,2])
idx.c = which(stage == "Normal")
idx.mf = which(stage == "MF")
idx.ss = which(stage == "SS")
expr = expr[, c(idx.mf)]
sample = sample[c(idx.mf),]

#Classify stage
expr = as.matrix(expr)
sample = sample[,-2]
rownames(sample) = sample$id
sample<-sample %>% 
  mutate(dx=recode(stage,
                   "tumor stage: IA"="local", 
                   "tumor stage: IB"="local", 
                   "tumor stage: IIA"="local", 
                   "tumor stage: IIB"="advance"))
sample = unlist(sample[,3])
idx.local = which(sample == "local")
idx.advance = which(sample == "advance")
expr = expr[c(idx.local, idx.advance), ]

#Label stage
local = paste0("Local_", c(1:length(idx.local)))
advance = paste0("Advance_", c(1:length(idx.advance)))
colnames(expr) = c(local, advance)

#PCA
dat1 = expr
dat1 = data.frame(dat1)
dx = c(rep("Early-stages",29),rep("Stage IIb",6))
dx = factor(dx, levels = c("Early-stages", "Stage IIb"))


#high_variance
sd = rowSds(as.matrix(dat1))
cut_off = quantile(sd, 0.8)
dat1 = dat1[sd > cut_off, ]

#PCA analysis
dat2 = t(dat1)
prcomp = prcomp(dat2)
pc = prcomp$x
pc = pc[,c(1:2)]
pc = data.frame(pc)
pc$group = dx

p1 = ggplot(pc, aes(x=PC1, y = PC2, color = group)) + 
  geom_point(size = 4) + 
  ggtitle("Principal components analysis GSE113113") +
  theme_bw() + 
  theme(legend.title = element_blank(),
        legend.text = element_blank(),
        axis.title.x = element_text(size = 24),
        axis.title.y = element_text(size = 24),
        axis.text = element_text(size = 20),
        plot.title = element_text(hjust = 0.5, size = 26)) +
  scale_color_manual(values = c("Early-stages" = "#0183b1", "Stage IIb" = "#eb6841"))

# Print the plot to screen (important!)
print(p1)

# Copy current plot to WMF
dev.copy(win.metafile, filename = "F.pca_GSE113113.wmf", width = 8, height = 6)
dev.off()

#save PCA result
library(ggpubr)
ggsave("F.pca_GSE113113.pdf", width = 8, height = 6)

#corellation matrix
corMatrix <- cor(dat1,use="c")
pheatmap(corMatrix, fontsize = 8, filename = "F.heatmap_GSE113113.pdf")     