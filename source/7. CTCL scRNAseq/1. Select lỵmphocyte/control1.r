# Clear environment
rm(list=ls())
options(stringsAsFactors = FALSE)

# Set working directory
setwd("D:/R.data/DHFR_project/CTCL scRNAseq/data")

# Load libraries
library(Seurat)
library(ggplot2)

# Load the dataset
tumor_data <- read.csv("control1.csv", row.names = 1)  # Assuming gene names are row names

# Create a Seurat object
seurat_obj <- CreateSeuratObject(counts = tumor_data)

# Normalize the data
seurat_obj <- NormalizeData(seurat_obj)

# Identify variable features
seurat_obj <- FindVariableFeatures(seurat_obj)

# Scale the data
seurat_obj <- ScaleData(seurat_obj)

# Get all gene names
gene_names <- rownames(GetAssayData(seurat_obj, assay = "RNA", layer = "counts"))

# Define CD3-related genes for T cell identification
cd3_genes <- c("CD3E", "CD3D", "CD3G")

# Check which CD3 genes exist in the dataset
cd3_genes_present <- intersect(cd3_genes, gene_names)

# Stop if none of the CD3 genes are found
if (length(cd3_genes_present) == 0) {
  stop("No CD3 genes found in dataset! Check available gene names using head(gene_names).")
} else {
  print(paste("Using CD3 genes:", paste(cd3_genes_present, collapse=", ")))
}

# Extract CD3 gene expression and sum across selected genes
cd3_expression <- colSums(GetAssayData(seurat_obj, assay = "RNA", layer = "counts")[cd3_genes_present, , drop = FALSE])

# Select cells where summed CD3 expression is greater than zero
cd3_positive_cells <- names(cd3_expression)[cd3_expression > 0]

# Subset the Seurat object to keep only CD3-positive cells
seurat_cd3_positive <- subset(seurat_obj, cells = cd3_positive_cells)

# Save the filtered cell matrix
write.table(GetAssayData(seurat_cd3_positive, assay = "RNA", layer = "counts"), 
            "control1_CD3_positive.txt", sep = "\t", row.names = TRUE, col.names = NA, quote = FALSE)

# Print success message
cat("Filtered data saved to control1_CD4_positive.txt\n")
