# Clear environment
rm(list=ls())
options(stringsAsFactors = FALSE)

# Set working directory
setwd("D:/R.data/DHFR_project/CTCL scRNAseq/data")

# Load the tumor1_CD3_positive.txt data
tumor1_cd3_positive_data <- read.table("tumor5_CD3_positive.txt", header = TRUE, row.names = 1)

# Check if DHFR is present in the row names (genes) of the dataset
if ("DHFR" %in% rownames(tumor1_cd3_positive_data)) {
  print("DHFR gene is found in the dataset.")
} else {
  stop("DHFR gene not found in the dataset!")
}

# Extract DHFR expression levels for CD3-positive cells (from tumor1 dataset)
dhfr_expression <- tumor1_cd3_positive_data["DHFR", , drop = FALSE]

# Define DHFR-positive cells (expression > 0)
dhfr_positive_cells <- colnames(dhfr_expression)[dhfr_expression > 0]

# Calculate the percentage of cells with DHFR-positive expression
total_cells <- ncol(tumor1_cd3_positive_data)  # Total number of cells
dhfr_positive_percentage <- (length(dhfr_positive_cells) / total_cells) * 100

# Print the percentage of DHFR-positive cells
print(paste("Percentage of DHFR-positive cells:", round(dhfr_positive_percentage, 2), "%"))
