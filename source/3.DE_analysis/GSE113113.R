rm(list=ls())
options(stringsAsFactors = F)
setwd("D:/R.data/DHFR_project/data")

library(data.table)
library(limma)
library(plyr)
library(dplyr)
library(DESeq2)
library(ggplot2)
library(ggrepel)

#load count data
expr <- fread("GSE113113_count.txt")
expr <- data.frame(expr, row.names = 1)

#load sample data
sample <- fread("GSE113113_sample.txt", fill = T)

#select stage
stage = unlist(sample[, 3])
table(stage)

#
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
expr = as.matrix(expr)
sample = sample[,-2]
rownames(sample) = sample$id
sample<-sample %>% 
  mutate(dx=recode(stage,
                   "tumor stage: IA"="local", 
                   "tumor stage: IB"="local", 
                   "tumor stage: IIA"="local", 
                   "tumor stage: IIB"="advance"))
sample = sample[,-2]
sample$dx = as.factor(sample$dx)
levels(sample$dx)
sample$dx = relevel(sample$dx, ref = "local")

#Differential expression analysis with DESeq2
dds = DESeqDataSetFromMatrix(countData = expr,
                             colData = sample,
                             design = ~ dx)

dds <- DESeq(dds)

results <- results(dds)


write.csv(as.data.frame(results), file = "DESeq2_results.csv")

result = read.csv("DESeq2_results.csv")

#export result
write.table(result, "GSE113113_deg.txt", sep = "\t", row.names = F)

# visualize and export significant results
result$sig = ifelse(result$log2FoldChange > 0 & result$padj<0.05, "Up regulation",
                    ifelse(result$log2FoldChange<0 & result$padj<0.05, "Down regulation", "Not Significant"))
result$Symbol = result$X

#label target genes
list = fread("Target_list.txt", header = F)
result$label <- ifelse(result$Symbol %in% list$V1, result$Symbol, NA)
result$'-log10Padj' = -log10(result$padj)

#volcano plot
result %>%
  ggplot(data = ., aes(x = log2FoldChange, y = `-log10Padj`, col = sig)) +  # Use data argument
  geom_point() +
  scale_color_manual(values = c("Up regulation" = "#fb4e53", "Down regulation" = "#0183b1","Not Significant" = "grey")) +
  geom_text_repel(aes(label = label), size = 9, box.padding = unit(.8, "lines"), show.legend = FALSE,color = "black") +
  ggtitle("GSE113113 skin: Stage IIB vs early-stages") +
  theme_bw() +
  theme(legend.title = element_blank(),
        legend.text = element_blank(),
        axis.title.x = element_text(size = 24),
        axis.title.y = element_text(size = 24),
        axis.text = element_text(size = 20),
        plot.title = element_text(hjust = 0.5, size = 26))

#save plot
library(ggpubr)
ggsave("F.DEG_GSE113113.pdf", width = 8.5, height = 8)

#significant filter and export
significant <- filter(result, padj<0.05)
write.table(significant, "GSE113113_deg2.txt", row.names = F, sep = "\t")

#export data for metal meta-analysis
diffexp <- data.frame(Symbol = result$Symbol,
                      log2FC = result$log2FoldChange,
                      SE = result$lfcSE,
                      pvalue = result$padj,
                      Sample_size = 35)
write.table(diffexp, "D:/R.data/DHFR_project/4.Metal/GSE113113.txt", row.names = FALSE, sep = "\t")