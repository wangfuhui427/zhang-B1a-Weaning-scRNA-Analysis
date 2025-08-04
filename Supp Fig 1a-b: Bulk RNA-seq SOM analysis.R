# --- Supp Fig 1a-b: Bulk RNA-seq SOM analysis ---
# Code for SOM analysis is complex and specific to the package used (e.g., 'kohonen').
# This script would be separate.
# Please refer to the script `S1_Bulk_RNAseq_Analysis.R` in the repository.

# --- Supp Fig 1c-e: Neonatal Splenic Analysis ---
# Assumes `seurat_neonatal` is the Seurat object for P7 splenic B-1a cells
# UMAP plot
DimPlot(seurat_neonatal, reduction = "umap", label = TRUE)

# Feature plot for progenitor genes
FeaturePlot(seurat_neonatal, features = c("Kit", "Il7r"))

# CytoTRACE analysis
# Assumes CytoTRACE has been run and scores are in metadata
FeaturePlot(seurat_neonatal, features = "CytoTRACE_score") +
    scale_color_viridis_c()
