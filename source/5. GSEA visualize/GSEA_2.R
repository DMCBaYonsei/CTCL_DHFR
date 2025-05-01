# Load the ggplot2 library
library(ggplot2)
library(data.table)
library(viridis)  # For a clear, perceptually uniform color palette
library(ggtext)

# Import data frame
data = fread("data.txt", header = TRUE)

data = as.data.frame(data)
data$Enrichment_Score = as.numeric(data$Enrichment_Score)
data$P_value = as.numeric(data$P_value)
data$LogP = log10(data$P_value)


# Base plot
p1 <- ggplot(data, aes(x = GSE_Dataset, y = Gene_set)) +
  geom_point(aes(size = Enrichment_Score, color = LogP)) +
  scale_color_viridis(option = "D", name = "-log10(P-value)") +
  scale_size_continuous(range = c(5, 15)) +
  labs(
    title = "Gene Set Enrichment",
    x = NULL,
    y = NULL,
    size = "Normalized Enrichment Score",
    color = "-log10(P-value)"
  ) +
  theme_minimal() +
  theme(
    text = element_text(size = 16),
    axis.text.x = element_text(angle = 0, hjust = 0.5, size = 16)
  ) +
  # Add label ONLY for "Oxidative phosphorylation"
  geom_text(
    data = subset(data, Gene_set == "Oxidative phosphorylation"),
    aes(x = GSE_Dataset, y = Gene_set, label = Gene_set),
    size = 5,
    color = "red",
    vjust = 0.5
  )

# Print
print(p1)

#Save as Windows Metafile
dev.copy(win.metafile, filename = "DotPlot_GeneSet_Enrichment.wmf", width = 14, height = 7)
dev.off()