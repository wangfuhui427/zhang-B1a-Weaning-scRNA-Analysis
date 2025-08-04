# Load necessary libraries
library(Seurat)
library(ggplot2)
library(dplyr)

# --- Figure 1b: Quality Control Violin Plots ---
# Assumes 'sample_id' is a column in the metadata
VlnPlot(seurat_object_peritoneal, 
        features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
        group.by = "sample_id", 
        pt.size = 0, 
        ncol = 3)

# --- Figure 1c: UMAP of Transcriptional Clusters ---
DimPlot(seurat_object_peritoneal, 
        reduction = "umap", 
        group.by = "seurat_clusters", 
        label = TRUE, 
        label.size = 5) + NoLegend()

# --- Figure 1d: Marker Gene Dot Plot ---
# List of key marker genes for the 12 clusters
marker_genes <- c("Fos", "Junb", "Cd25", "Hspa1a", "Ighd", "Cdkn1a", "Mki67", ...) # Add your top markers
DotPlot(seurat_object_peritoneal, 
        features = marker_genes, 
        group.by = "seurat_clusters") + 
        RotatedAxis()

# --- Figure 1e: Cluster Proportion Bar Plot ---
metadata <- seurat_object_peritoneal@meta.data
metadata %>%
  group_by(timepoint, seurat_clusters) %>%
  summarise(n = n()) %>%
  mutate(freq = n / sum(n)) %>%
  ggplot(aes(x = timepoint, y = freq, fill = seurat_clusters)) +
  geom_bar(stat = "identity", position = "fill") +
  theme_classic() +
  scale_y_continuous(labels = scales::percent) +
  ylab("Proportion")
