rm(list = ls())
setwd("D:/R.data/DHFR_project/data/")
library(data.table);library(matrixStats);library(limma);library(WGCNA)
library(ggplot2)
library(dplyr)
library(pheatmap)
library(ggpubr)

#
dat1 = fread("GSE59307_norm.txt")
dat1 = data.frame(dat1)
rownames(dat1) = dat1[,1]
dat1 = dat1[,-1]
dx = c(rep("Normal",8),rep("CTCL",14))
dx = factor(dx, levels = c("Normal","CTCL"))


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
  ggtitle("Principal components analysis GSE59307") +
  theme_bw() + 
  theme(legend.title = element_blank(),
        legend.text = element_blank(),
        axis.title.x = element_text(size = 24),
        axis.title.y = element_text(size = 24),
        axis.text = element_text(size = 20),
        plot.title = element_text(hjust = 0.5, size = 26)) +
  scale_color_manual(values = c("Normal" = "#0183b1", "CTCL" = "#eb6841"))

# Print the plot to screen (important!)
print(p1)

# Copy current plot to WMF
dev.copy(win.metafile, filename = "F.pca_GSE59307.wmf", width = 8, height = 6)
dev.off()

#save PCA result
ggsave("F.pca_GSE59307.pdf", width = 8, height = 6)

#corellation matrix
corMatrix <- cor(dat1,use="c")
pheatmap(corMatrix, fontsize = 10, filename = "F.heatmap_GSE59307.pdf")