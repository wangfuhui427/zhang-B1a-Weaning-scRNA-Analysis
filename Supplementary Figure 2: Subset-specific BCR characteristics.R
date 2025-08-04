# --- Supp Fig 2a: VH-VL Pairing Heatmap ---
# This requires custom code to count pairs and ggplot2 with geom_tile
vh_vl_pairs <- bcr_data %>%
  group_by(timepoint, v_gene_heavy, v_gene_light) %>%
  summarise(count = n()) %>%
  # ... further processing to get top 15 and create a matrix ...

# --- Supp Fig 2b: Volcano Plot of DEGs ---
# Requires differential expression results between Ighv11-2 and Ighv12-3 cells
# A package like EnhancedVolcano is excellent for this
library(EnhancedVolcano)
EnhancedVolcano(deg_results,
                lab = rownames(deg_results),
                x = 'avg_log2FC',
                y = 'p_val_adj')

# --- Supp Fig 2c: UMAP of Specific Clone Expansion ---
weaning_clone_barcodes <- bcr_data %>%
  filter(CTaa == "CMRYGNYWYFDVW_CLQHGESPYTF", timepoint == "5wk") %>%
  pull(barcode)
DimPlot(seurat_object_peritoneal, 
        cells.highlight = weaning_clone_barcodes,
        cols.highlight = "blue")

# --- Supp Fig 2d: Rare vs. Dominant Clone Proportions ---
# Requires calculating clone sizes and plotting with ggplot2
clone_sizes <- bcr_data %>%
  group_by(seurat_clusters, CTaa) %>%
  summarise(clone_size = n())
# ... further processing to categorize as rare/dominant and plot proportions ...
