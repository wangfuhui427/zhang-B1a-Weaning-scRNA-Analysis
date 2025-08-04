# --- Figure 2a: Integrated UMAP by Time Point ---
DimPlot(seurat_object_integrated, 
        reduction = "umap", 
        group.by = "timepoint")

# --- Figure 2b: UMAP of Cell Cycle Phase ---
# Assumes CellCycleScoring has been run and "Phase" is in metadata
DimPlot(seurat_object_integrated, 
        reduction = "umap", 
        group.by = "Phase")

# --- Figure 2c: Feature Plots of Key Genes ---
FeaturePlot(seurat_object_integrated, 
            features = c("Mki67", "Isg15", "Zbtb32"), 
            ncol = 3)

# --- Figure 2d: ClusterGVis Heatmap Preparation ---
# This code prepares the data for ClusterGVis
pseudo_bulk_data <- AggregateExpression(seurat_object_integrated, 
                                        group.by = c("timepoint"), 
                                        assays = "RNA")$RNA
# The resulting matrix 'pseudo_bulk_data' is used as input for the ClusterGVis tool.

# --- Figure 2e: Functional Enrichment Plot ---
# Assumes 'go_results' is a dataframe from Enrichr analysis of ClusterGVis modules
ggplot(go_results, aes(x = reorder(Term, -log10(p.value)), y = -log10(p.value))) +
  geom_bar(stat = "identity", fill = "steelblue") +
  coord_flip() +
  theme_minimal() +
  xlab("GO Biological Process")
