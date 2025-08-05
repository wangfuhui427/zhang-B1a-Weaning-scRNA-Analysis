# --- Figure 2a: Integrated UMAP by Time Point ---
DimPlot(seurat_object_integrated, 
        reduction = "umap", 
        group.by = "timepoint")
# --- Figure 2c: Dotplot of markers ---
library(Seurat)
library(ggplot2)
library(dplyr)

# --- Step 1: Define your list of selected marker genes ---
# Using a named list or ordering the vector can help control the gene order on the plot.
# Let's ensure the order is exactly as we want it.
marker_genes_list <- c("Rhob","Fos", "Fosb","Junb",
                       "H2-Q7","Pgap1","Sspn",
                       "Fbxw13", "Lifr","Fcgr2b",
                      "Lmo2", "Cfp","Hsp90b1",
                      "Ighd","Sell","Cd55","Abca1",
                      "Apoe","Plaur",
                      "Egr3", "Marcksl1","Nfkbid",
                      "Ahnak2","Zcwpw1","Ccr9","Pou2af2",
                      "Mki67", "Top2a", "Cenpe","Ccnb2",
  "Bhlhe41","Cdkn1a", "Tnfsf8","Psrc1","Ifi206","Slfn5", "Ifit3" ,"Hspa1a","Hspa1b","Hsph1" # Add more genes as needed
)
# Make sure the gene names are unique and present in your Seurat object.
marker_genes_list <- unique(marker_genes_list) 

# --- Step 2: Generate and Customize the Dot Plot ---
p_dotplot_custom <- DotPlot(seurat_object, 
                            features = marker_genes_list,
                            group.by = "seurat_clusters",
                            scale = TRUE) + # `scale = TRUE` ensures expression is z-scored within each gene
  scale_colour_gradient2(low = "#5381ec", mid = 'lightgrey', high = "#dd6d60" ) + # Custom color scale for expression
  labs(
    title = "Expression of Key Marker Genes Across Clusters",
    x = "Genes",
    y = "Transcriptional Clusters"
  ) +
  theme_bw() + 
  theme(
   
    panel.border = element_rect(colour = "black", fill = NA, size = 1), 
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(), 
    
    
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 10, colour = "black"), 
    axis.text.y = element_text(size = 10, colour = "black"),
    axis.title = element_text(size = 12, face = "bold", colour = "black"),
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold", colour = "black"),
    
    
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 8) 
  )
  
# Display the plot
p_dotplot_custom


# --- Figure : Feature Plots of Key Genes ---
FeaturePlot(seurat_object_integrated, 
            features = c("Mki67", "Isg15", "Zbtb32"), 
            ncol = 3)
# --- Figure 2d :ssGSEA Plots of cell cycle Genesets ---

####
gene_sets <- lapply(gene_sets, function(genes) {
  genes[genes != "" & !is.na(genes)]  
})


if (is.null(rownames(data_matrix))) {
  stop("Row names of the expression matrix must be gene names.")
}

gene_set_scores <- list()
for (gene_set_name in names(gene_sets)) {

  selected_genes <- gene_sets[[gene_set_name]]
  common_genes <- intersect(rownames(data), selected_genes)
  if (length(common_genes) == 0) {
    warning(paste("No matching genes found for", gene_set_name))
    next
  }
  subset_expr <- data[common_genes, , drop = FALSE]
  gene_set_scores[[gene_set_name]] <- colMeans(subset_expr, na.rm = TRUE)
}
gene_set_scores_df <- as.data.frame(gene_set_scores)
gene_set_scores_df$Sample <- metadata $Sample
gene_set_scores_df$cluster <-metadata$seurat_clusters
gene_set_scores_df$IGHV_group <-metadata$group
gene_set_scores_df$S.score <-metadata$S.Score
gene_set_scores_df$G2M.score <-metadata$G2M.Score
gene_set_scores_df$highlight<-metadata$highlight
gene_set_scores_df$CTaa<-metadata$CTaa

ggplot(sub_umap_data, aes(x = UMAP_1, y = UMAP_2, color = geneset)) +
  geom_point(alpha = 1, size = 0.5) +
  facet_wrap(~ Sample, nrow = 2, ncol = 5) + 
  scale_color_gradient(low = "white", high = "red") +
  theme_minimal() +
  #theme(
  #panel.grid = element_blank(),    
  #panel.background = element_blank()
  #) +
  labs(title = "KEGG_apoptosis Expression in UMAP by Sample",
       x = "UMAP 1", y = "UMAP 2", color = "Expression")
# --- Figure 2e: ClusterGVis Heatmap Preparation ---
# This code prepares the data for ClusterGVis
pseudo_bulk_data <- AggregateExpression(seurat_object_integrated, 
                                        group.by = c("timepoint"), 
                                        assays = "RNA")$RNA
# The resulting matrix 'pseudo_bulk_data' is used as input for the ClusterGVis tool.
fpkm_matrix<-fpkm_matrix[complete.cases(fpkm_matrix),]
fpkm_matrix =fpkm_matrix[apply(fpkm_matrix, 1, function(x) sd(x)!=0),] 
fpkm_matrix<- fpkm_matrix[rowMeans(fpkm_matrix) > 0.2, ]
cm <- clusterData(obj  = fpkm_matrix,
                  cluster.method = "mfuzz",
                  cluster.num = 6)
# --- Figure 2f: Functional Enrichment Plot ---
# Assumes 'go_results' is a dataframe from Enrichr analysis of ClusterGVis modules
                               gene_df<-merged_df
gene_df$geneCluster<-gene_df$Cluster
# 数据库
dbs <- c("GO_Biological_Process_2021", 
         "GO_Cellular_Component_2021", 
         "GO_Molecular_Function_2021", 
         "KEGG_2021_Human")
httr::set_config(httr::config(connecttimeout = 30))
# 分组分析
enrich_results <- gene_df %>%
  group_by(geneCluster) %>%
  group_split() %>%
  map_dfr(function(df) {
    genes <- unique(df$Gene)
    if (length(genes) < 5) return(NULL)  # 跳过太小的集合
    enriched <- enrichr(genes, dbs)
    
    # 筛选并合并结果
    bind_rows(lapply(names(enriched), function(dbname) {
      enriched[[dbname]] %>%
        filter(P.value < 0.1) %>%
        mutate(Database = dbname, Cluster = unique(df$geneCluster))
    }))
  })
enrich_results$Overlap<-as.character(enrich_results$Overlap)

enrich_results$Overlap <- paste0('"', enrich_results$Overlap, '"')
write.csv(enrich_results,"4scKO_clustergvis_enrichr.csv")
ggplot(go_results, aes(x = reorder(Term, -log10(p.value)), y = -log10(p.value))) +
  geom_bar(stat = "identity", fill = "steelblue") +
  coord_flip() +
  theme_minimal() +
  xlab("GO Biological Process")
