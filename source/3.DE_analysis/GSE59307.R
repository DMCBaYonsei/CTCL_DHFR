rm(list = ls())
setwd("D:/R.data/DHFR_project/data")
library(data.table)
library(limma)
library(ggplot2)
library(dplyr)
library(ggrepel)

#load normalization data
data = "GSE59307"
storage = paste0(data, "_norm.txt")
dat1 = fread(storage)
dat1 = data.frame(dat1, stringsAsFactors = F)
rownames(dat1) = dat1[, 1]
dat1 = dat1[, -1]

#Differential expression analysis with limma
Disease = c(rep(0, 8), rep(1, 14))
design = model.matrix(~ Disease)
fit = lmFit(dat1, design)
fit = eBayes(fit, trend = TRUE)
result = topTable(fit, coef=2, n = nrow(dat1), adjust="fdr", p = 1)
result = result[rownames(dat1), ]

#export result
storage = paste0(data, "_deg.txt")
write.table(result, storage, sep = "\t")

# visualize and export significant results
result$sig = ifelse(result$logFC > 0 & result$adj.P.Val<0.05, "Up regulation",
                    ifelse(result$logFC<0 & result$adj.P.Val<0.05, "Down regulation", "Not Significant"))
result$Symbol = rownames(result)

#label target genes
list = fread("Target_list.txt", header = F)
result$label <- ifelse(result$Symbol %in% list$V1, result$Symbol, NA)

result$'-log10Padj' = -log10(result$adj.P.Val)
result$log2FoldChange = result$logFC

#volcano plot
p1 <- result %>%
  ggplot(data = ., aes(x = log2FoldChange, y = `-log10Padj`, col = sig)) +  # Use data argument
  geom_point() +
  scale_color_manual(values = c("Up regulation" = "#fb4e53", "Down regulation" = "#0183b1", "Not Significant" = "grey")) +
  geom_text_repel(aes(label = label), size = 9, box.padding = unit(.8, "lines"), show.legend = FALSE,color = "black") +
  ggtitle("GSE59307 skin: CTCL vs Normal") +
  theme_bw() +
  theme(legend.title = element_blank(),
        legend.text = element_blank(),
        axis.title.x = element_text(size = 24),
        axis.title.y = element_text(size = 24),
        axis.text = element_text(size = 20),
        plot.title = element_text(hjust = 0.5, size = 26))

# Print the plot to screen (important!)
print(p1)

# Copy current plot to WMF
dev.copy(win.metafile, filename = "F.DEG_GSE59307.wmf", width = 8.5, height = 8)
dev.off()

#save plot
library(ggpubr)
ggsave("F.DEG_GSE59307.pdf", width = 8.5, height = 8)

#export significant result
significant <- filter(result, adj.P.Val<0.05)
storage2 = paste0(data, "_deg2.txt")
write.table(significant, storage2, sep = "\t")

#export data for metal meta-analysis
diffexp <- data.frame(Symbol = result$Symbol,
                      log2FC = result$logFC,
                      SE = result$logFC/result$t,
                      pvalue = result$adj.P.Val,
                      Sample_size = 51)
write.table(diffexp, "D:/R.data/DHFR_project/4.Metal/GSE59307.txt", row.names = FALSE, sep = "\t")