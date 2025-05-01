rm(list = ls())
setwd("D:/R.data/DHFR_project/data/")
library(data.table);library(matrixStats);library(limma);library(WGCNA)
library(ggplot2)
library(dplyr)
library(pheatmap)
library(ggpubr)

#
dat1 = fread("GSE9479_norm.txt")
dat1 = data.frame(dat1)
rownames(dat1) = dat1[,1]
dat1 = dat1[,-1]
dx = c(rep("Early-stages",43),rep("Stage IIb",8))
dx = factor(dx, levels = c("Early-stages","Stage IIb"))

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
  ggtitle("Principal components analysis GSE9479") +
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
dev.copy(win.metafile, filename = "F.pca_GSE9479.wmf", width = 8, height = 6)
dev.off()

#save PCA result
ggsave("F.pca_GSE9479.pdf", width = 8, height = 6)

#corellation matrix
corMatrix <- cor(dat1,use="c")
pheatmap(corMatrix, fontsize = 10, filename = "F.heatmap_GSE9479.pdf")