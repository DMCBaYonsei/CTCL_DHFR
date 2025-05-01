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

#export meta result
write.table(dat1, "meta_sign.txt", row.names = FALSE, sep = "\t")

#label target genes
list = fread("Target_list.txt", header = F)
dat1$label <- ifelse(dat1$Gene %in% list$V1, dat1$Gene, NA)
dat1$'-log10Padj' = -log10(dat1$Padj)
dat1$log2FoldChange = dat1$logFC

#volcano plot
dat1 %>%
  ggplot(data = ., aes(x = log2FoldChange, y = `-log10Padj`, col = sig)) +  # Use data argument
  geom_point() +
  scale_color_manual(values = c("Up regulation" = "#fb4e53", "Down regulation" = "#0183b1",  "Not Significant" = "grey")) +
  geom_text_repel(aes(label = label), size = 9, box.padding = unit(.8, "lines"), show.legend = FALSE,color = "black") +
  ggtitle("Metal meta-analysis result") +
  theme_bw() +
  theme(legend.title = element_blank(),
        legend.text = element_blank(),
        axis.title.x = element_text(size = 24),
        axis.title.y = element_text(size = 24),
        axis.text = element_text(size = 20),
        plot.title = element_text(hjust = 0.5, size = 26))
#save plot
library(ggpubr)
ggsave("F.DEG_meta.pdf", width = 8.5, height = 8)
