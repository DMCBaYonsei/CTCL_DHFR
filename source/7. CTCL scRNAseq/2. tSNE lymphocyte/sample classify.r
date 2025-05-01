# Clear environment
rm(list=ls())
options(stringsAsFactors = FALSE)

# Load libraries
library(Seurat)
library(ggplot2)

# Set working directory
setwd("D:/R.data/DHFR_project/CTCL scRNAseq/data")

# List of dataset file names
file_list <- c("tumor1_CD3_positive.txt", "tumor2_CD3_positive.txt", "tumor3_CD3_positive.txt", 
               "tumor4_CD3_positive.txt", "tumor5_CD3_positive.txt",
               "control1_CD3_positive.txt", "control2_CD3_positive.txt", 
               "control3_CD3_positive.txt", "control4_CD3_positive.txt")

# Create an empty list to store Seurat objects
seurat_list <- list()

# Load and create Seurat objects for each dataset
for (file in file_list) {
  sample_name <- gsub("_CD3_positive.txt", "", file)  # Extract sample name
  
  # Read dataset
  data <- read.table(file, header = TRUE, row.names = 1, sep = "\t")
  
  # Create Seurat object
  seurat_obj <- CreateSeuratObject(counts = data, project = sample_name)
  
  # Assign metadata for group classification (CTCL vs. Healthy)
  seurat_obj$Group <- ifelse(grepl("tumor", sample_name), "CTCL", "Healthy Control")
  
  # Store in list
  seurat_list[[sample_name]] <- seurat_obj
}

# Merge all Seurat objects into one
combined_seurat <- merge(seurat_list[[1]], y = seurat_list[-1], add.cell.ids = names(seurat_list))

# Normalize the data
combined_seurat <- NormalizeData(combined_seurat)

# Identify variable features
combined_seurat <- FindVariableFeatures(combined_seurat)

# Scale the data
combined_seurat <- ScaleData(combined_seurat)

# Run PCA
combined_seurat <- RunPCA(combined_seurat)

# Run t-SNE
combined_seurat <- RunTSNE(combined_seurat, dims = 1:20)

# Save t-SNE plot as PDF
pdf("tsne_classify.pdf", width = 7.8, height = 6)

# Generate and print the t-SNE plot
DimPlot(combined_seurat, reduction = "tsne", group.by = "Group") + 
  ggtitle("t-SNE of lymphocytes")

# Close the PDF device
dev.off()

cat("t-SNE plot saved as tsne_classify.pdf\n")
