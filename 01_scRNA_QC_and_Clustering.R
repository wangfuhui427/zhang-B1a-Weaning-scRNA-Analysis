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
#
library(devtools)
all_markers <- FindAllMarkers(seurat_object,
                              only.pos = TRUE,     
                              min.pct = 0.25,      
                              logfc.threshold = 0.25)
write.csv(all_markers,"10group_markers.csv")

top20_markers <- all_markers %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 20)
valid_genes <- intersect(top20$gene, rownames(seurat_object@assays$RNA@scale.data))

heatmap_plot <- DoHeatmap(seurat_object, 
                          features = valid_genes, 
                          size = 4) +
  theme(axis.text.y = element_text(size = 4))
heatmap_plot <- DoHeatmap(seurat_object, 
                          features = valid_genes, 
                          group.by = "seurat_clusters",assay = "RNA",
                          group.colors = c("0" = "#82B5E9", "1" = "#EAE939", "2" = "#ACD030",  "3"="#3e4f94","4"="#8BD98A",
                                           "5"="#35b0ab","6"="#FF6155", "7" ="#b54764"),
                          size = 2) +
  theme(axis.text.y = element_text(size = 4))

heatmap_plot
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
# --- Figure 1f: Cluster marker umap Featureplot---
library(Seurat)
library(ggplot2)
features <- c("Fos", "Fosb", "Jun", "Junb", "Lifr", "Il2ra", 
              "Hspa1a", "Hspa1b", "Ighd", "Fcer2a", "Cdkn1a", "Mki67")
#FeaturePlot
plots <- FeaturePlot(seurat_object,
                     features = features,
                     cols = c("lightgrey", "#FF6155"),
                     reduction = "umap")
plots <- lapply(plots, function(p) {
  p + theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    axis.line = element_blank(),
    panel.grid = element_blank(),
    panel.border = element_blank()
  )
})
