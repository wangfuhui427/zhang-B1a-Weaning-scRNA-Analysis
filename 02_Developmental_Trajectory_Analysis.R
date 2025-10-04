# ---Integrated UMAP by Time Point ---
sub_obj <- WhichCells(seurat_object, idents = c("WT_1WK","WT_3WK",""WT_5WK",""WT_7WK","WT_10WK"))
sub_obj <- subset(seurat_object, cells = sub_obj)
B1A.list <- SplitObject(sub_obj, split.by = "Sample")
B1A.list <- lapply(X = B1A.list, FUN = function(x) {
  DefaultAssay(x) <- "RNA"
  x@assays$RNA@scale.data <- matrix()
  x@assays$RNA@meta.features <- data.frame(row.names = rownames(x@assays$RNA@data))  # <<<<<< 修复点
  x@assays$RNA@var.features <- character()
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
  
  return(x)
})

features <- SelectIntegrationFeatures(object.list = B1A.list)
cell.anchors <- FindIntegrationAnchors(object.list = B1A.list,anchor.features = features, dims = 1:20)
cell.combined <- IntegrateData(anchorset = cell.anchors, dims = 1:20)

subset_object<-cell.combined
subset_object <- ScaleData(subset_object)
subset_object <- RunPCA(subset_object, npcs = 30)
subset_object <- RunUMAP(subset_object, dims = 1:30)
subset_object <- RunTSNE(subset_object,  dims = 1:30)
subset_object<- FindNeighbors(subset_object,  dims = 1:30)
subset_object <- FindClusters(subset_object, resolution = 0.5)

# Visualization
p1 <- DimPlot(subset_object, reduction = "umap", label=TRUE, group.by = "seurat_clusters",raster = FALSE)
p2 <- DimPlot(subset_object, reduction = "umap", label = TRUE, repel = TRUE,split="Sample",raster = FALSE)
p1
DimPlot(seurat_object_integrated, 
        reduction = "umap", 
        group.by = "Sample")
DimPlot(seurat_object_integrated, 
        reduction = "umap", 
        group.by = "timepoint")
# --- Dotplot of markers ---
library(Seurat)
library(ggplot2)
library(dplyr)

# --- Step 1: Define marker genes ---
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


# ---  Feature Plots of Key Genes ---
FeaturePlot(seurat_object_integrated, 
            features = c("Mki67", "Isg15", "Zbtb32"), 
            ncol = 3)

# ---cytotrace2
devtools::install_github("digitalcytometry/cytotrace2", subdir ="cytotrace2_r")
library(CytoTRACE2)
result_sce <- cytotrace2(seurat_object,
                         is_seurat = TRUE,
                         slot_type = "counts",
                         species = "mouse")
head(result_sce,2)
annotation <- data.frame(phenotype=result_sce@meta.data$seurat_clusters) %>% set_rownames(., colnames(result_sce))
plots <- plotData(cytotrace2_result=result_sce, annotation=annotation, is_seurat=TRUE)
p1 <- plots$CytoTRACE2_UMAP
p2 <- plots$CytoTRACE2_Potency_UMAP
p3 <- plots$CytoTRACE2_Relative_UMAP
p4 <- plots$CytoTRACE2_Boxplot_byPheno
p5<-plots$Phenotype_UMAP
library(patchwork)
(p1+p2+p3+p4) + plot_layout(ncol = 2)
p5
# ---monocle3 pseudotime analysis---
library(monocle3)
data <- GetAssayData(seurat_object, assay = 'RNA', slot = 'counts')
cell_metadata <- seurat_object@meta.data
gene_annotation <- data.frame(gene_short_name = rownames(data))
rownames(gene_annotation) <- rownames(data)
cds <- new_cell_data_set(data,
                         cell_metadata = cell_metadata,
                         gene_metadata = gene_annotation)

cds <- preprocess_cds(cds, num_dim = 50)

cds <- reduce_dimension(cds, preprocess_method = "PCA",reduction_method = "UMAP")

library(ggplot2)
p1 <- plot_cells(cds, reduction_method="UMAP", color_cells_by="seurat_clusters",
                 group_label_size=5) + ggtitle('cds.umap')+
  scale_color_manual(values = new_colors) 
new_colors <- c("0" = "#82B5E9", "1" = "#EAE939", "2" = "#ACD030",  "3"="#3e4f94","4"="#8BD98A",
                "5"="#35b0ab","6"="#FF6155", "7" ="#b54764")
p1

cds.embed <- cds@int_colData$reducedDims$UMAP
int.embed <- Embeddings(seurat_object, reduction = "umap")
int.embed <- int.embed[rownames(cds.embed),]
cds@int_colData$reducedDims$UMAP <- int.embed
p2 <- plot_cells(cds, reduction_method="UMAP", color_cells_by="seurat_clusters",
                 group_label_size=5) + ggtitle('int.umap')+scale_color_manual(values = new_colors)
p1+p2

cds<-cluster_cells(cds)
plot_cells(cds,color_cells_by = "partition")
cds<-learn_graph(cds,use_partition=FALSE)                                                                              
p = plot_cells(cds, color_cells_by = "seurat_clusters", label_groups_by_cluster = TRUE, label_leaves = TRUE, 
               label_branch_points = TRUE)
p
cds <- order_cells(cds)
plot_cells(cds, color_cells_by = "pseudotime", label_cell_groups = FALSE, label_leaves = FALSE,
           label_branch_points = FALSE, graph_label_size = 1.5)
get_earliest_principal_node <- function(cds, time_bin = "2") {
  cell_ids <- which(colData(cds)[, "seurat_clusters"] == time_bin)
  
  closest_vertex <- cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
  root_pr_nodes <- igraph::V(principal_graph(cds)[["UMAP"]])$name[as.numeric(names(which.max(table(closest_vertex[cell_ids,
  ]))))]
  
  root_pr_nodes
}
cds <- order_cells(cds, root_pr_nodes = get_earliest_principal_node(cds))
plot_cells(cds,color_cells_by = "pseudotime",label_cell_groups=FALSE,
           label_leaves=FALSE,label_branch_points=FALSE,graph_label_size=1.5)

library(magrittr)
library(dplyr)
Track_genes <- graph_test(cds, neighbor_graph="principal_graph", cores=4)
write.csv(Track_genes,"1wkWT-trackgenes.csv")
Track_genes_sig <- Track_genes %>% top_n(n=40, morans_I) %>%
  pull(gene_short_name) %>% as.character()

data1<-cds[Track_genes_sig,]
Track_genes_sig_self<-c("Igll1","Vpreb1","Il7r","Lig1","Top2a","Mcm5","Pclaf","Mcm6","Ccna2","Stmn1")
plot_genes_in_pseudotime(cds[Track_genes_sig_self,], color_cells_by="seurat_clusters", cell_size = 0.25,
                         min_expr=0.2, ncol = 2)+scale_color_manual(values = new_colors)

plot_cells(cds, genes=Track_genes_sig, show_trajectory_graph=FALSE,
           label_cell_groups=FALSE,  label_leaves=FALSE)

genelist <- pull(Track_genes, gene_short_name) %>% as.character()
gene_module <- find_gene_modules(cds[genelist,], resolution=1e-2, cores = 10)
cell_group <- tibble::tibble(cell=row.names(colData(cds)), 
                             cell_group=colData(cds)$seurat_clusters)
agg_mat <- aggregate_gene_expression(cds, gene_module, cell_group)
row.names(agg_mat) <- stringr::str_c("Module ", row.names(agg_mat))
pheatmap::pheatmap(agg_mat, scale="column", clustering_method="ward.D2",fontsize = 3)
write.csv(gene_module,"5group_pseudotime_monocle3_module.csv")

genes = c("Il7r", "Mki67", "Hmga2","Ttk","Plk1","Cenpe","Ccnb2")
lineage_cds = cds[rowData(cds)$gene_short_name %in% genes,
                      clusters(cds) %in% c("6","2","4","7","5","3","1","0")]
lineage_cds<-order_cells(lineage_cds,root_pr_nodes = get_earliest_principal_node(cds))
lineage_cds<-order_cells(lineage_cds,root_pr_nodes = "2")
plot_genes_in_pseudotime(lineage_cds,
                         color_cells_by="seurat_clusters",
                         min_expr=0.5)
                         facet_by="clusters")

samples <- unique(colData(lineage_cds)$Sample)

for (gene in genes) {
  for (s in samples) {
    cat("Plotting", gene, "for sample:", s, "\n")
    # 子集化 CellDataSet
    sub_cds <- lineage_cds[, colData(lineage_cds)$Sample == s]
    p <- plot_genes_in_pseudotime(sub_cds[gene, ],
                                  color_cells_by = "Sample",
                                  min_expr = 0.5) +
      ggtitle(paste0(gene, " - ", s))
    print(p)
  }
}
pdf("pseudotime_by_sample.pdf", width = 6, height = 4)
# loop 内 print(p)
dev.off()

ciliated_genes = c("Ighm", "Mki67", "Hmga2","Ttk","Setd2","Dnmt3a")
cds_subset = cds[rowData(cds)$gene_short_name %in% ciliated_genes,]
gene_fits = fit_models(cds_subset, model_formula_str = "~seurat_clusters")

fit_coefs = coefficient_table(gene_fits)

emb_time_terms = fit_coefs %>% filter(term == "seurat_clusters")
# coefficient_table()use Benjamini and Hochberg（BH）
emb_time_terms %>% filter (q_value < 0.05) %>%
  select(gene_short_name, term, q_value, estimate)

plot_genes_violin(cds_subset, group_cells_by="seurat_clusters", ncol=2) +
  theme(axis.text.x=element_text(angle=45, hjust=1))
plot_genes_in_pseudotime(cds_subset,
                         color_cells_by="seurat_clusters",
                         min_expr=0.5)
plot_cells(cds,
           genes=gene_module %>% filter(module %in% c(100,74,171,120,41,137)),
           label_cell_groups=FALSE,
           show_trajectory_graph=FALSE)

# ---geneset analysis---
# --- ssGSEA Plots of cell cycle Genesets ---
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
# ---ClusterGVis Heatmap Preparation ---
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
# ---Functional Enrichment Plot ---
# Assumes 'go_results' is a dataframe from Enrichr analysis of ClusterGVis modules
                               gene_df<-merged_df
gene_df$geneCluster<-gene_df$Cluster

dbs <- c("GO_Biological_Process_2021", 
         "GO_Cellular_Component_2021", 
         "GO_Molecular_Function_2021", 
         "KEGG_2021_Human")
httr::set_config(httr::config(connecttimeout = 30))
enrich_results <- gene_df %>%
  group_by(geneCluster) %>%
  group_split() %>%
  map_dfr(function(df) {
    genes <- unique(df$Gene)
    if (length(genes) < 5) return(NULL)  
    enriched <- enrichr(genes, dbs)
    
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
