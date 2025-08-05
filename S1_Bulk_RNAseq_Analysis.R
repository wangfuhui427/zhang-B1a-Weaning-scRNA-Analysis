# --- Supp Fig 1a-b: Bulk RNA-seq SOM analysis ---
# Code for SOM analysis is complex and specific to the package used (e.g., 'kohonen').
# --- 1. Load Libraries ---
library(ggplot2)
library(clusterProfiler)
library(ggrepel)
install.packages(c("kohonen")) # Installs the 'kohonen' package if not already present
library(kohonen)
library(RColorBrewer)

# --- 2. Data Loading and Preprocessing ---
# Load the data from a CSV file. Assumes the first column contains gene names.
data <- read.csv("/Users/wangfuhui/Desktop/data.csv", header = TRUE, check.names = FALSE)

# Remove duplicate gene names to ensure unique rownames
data <- data[!duplicated(data$gene_name), ] 

# Set gene names as rownames and remove the original gene name column
rownames(data) <- data[, 1] 
data <- data[, -1] 

# Get a summary of the data
summary(data)

# Filter out genes with zero variance (no change across samples)
data <- data[apply(data, 1, function(x) sd(x) != 0), ] 

# Filter out genes with very low expression (row sum < 0.4)
data <- data[rowSums(data) >= 0.4, ]

# --- 3. Data Scaling ---
# Scale the data so that each gene (row) has a mean of 0 and a standard deviation of 1
data_scale <- as.matrix(t(scale(t(data))))
names(data_scale) <- names(data) # Retain original column names
head(data_scale)

# --- 4. SOM Training ---
# Define the size and shape of the SOM grid (e.g., 10x10 rectangular)
som_grid <- somgrid(xdim = 10, ydim = 10, topo = "rectangular")

# Train the SOM model using the scaled data
som_model <- supersom(data_scale, grid = som_grid, keep.data = TRUE)

# --- 5. SOM Quality and Training Visualization ---
# Visualize the training process. The plot shows the mean distance to the closest unit
# for each iteration. The distance should decrease and stabilize, indicating convergence.
plot(som_model, type = "changes")

# Visualize the number of genes mapped to each node (unit) in the SOM grid.
# This helps to see the data distribution on the map.
coolBlueHotRed <- function(n, alpha = 0.7) {
  rainbow(n, end = 4/6, alpha = alpha)[n:1]
}
plot(som_model, type = "counts", main = "Node Counts", palette.name = coolBlueHotRed)

# Visualize the quality of the SOM map (mean distance of objects in a unit to the unit's codebook vector).
# Lower values (bluer colors) indicate better quality/homogeneity within the node.
plot(som_model, type = "quality", main = "Node Quality/Distance", palette.name = coolBlueHotRed)

# Visualize the U-Matrix (Unified Distance Matrix).
# Darker colors indicate larger distances between neighboring nodes, suggesting potential cluster boundaries.
plot(som_model, type = "dist.neighbours", main = "SOM Neighbour Distances", palette.name = grey.colors)

# Visualize the codebook vectors for each node.
# Each line represents a variable (sample) and shows its trend across the SOM grid.
plot(som_model, type = "codes", codeRendering = "lines")

# --- 6. Extracting Gene Information from the SOM ---
# Get the mapping of each gene to its corresponding SOM node (unit).
table(som_model$unit.classif)

# Create a data frame containing each gene and the SOM node it was mapped to.
# The nodes are numbered from bottom-left to top-right.
som_model_code_class <- data.frame(name = rownames(data_scale), code_class = som_model$unit.classif)
head(som_model_code_class)

# --- 7. Clustering the SOM Nodes ---
# Determine the optimal number of clusters using the Elbow Method (WCSS).
mydata <- as.matrix(as.data.frame(som_model$codes))
wss <- (nrow(mydata) - 1) * sum(apply(mydata, 2, var))
for (i in 2:20) wss[i] <- sum(kmeans(mydata, centers = i)$withinss)
par(mar = c(5.1, 4.1, 4.1, 2.1))
plot(1:20, wss, type = "b", xlab = "Number of Clusters",
     ylab = "Within groups sum of squares", main = "Within Cluster Sum of Squares (WCSS)")

# Perform hierarchical clustering on the codebook vectors to group the SOM nodes.
# Here we choose 7 clusters based on the elbow plot or biological knowledge.
som_cluster <- cutree(hclust(dist(mydata)), 7)

# Define a color palette for the clusters
cluster_palette <- function(x, alpha = 0.6) {
  n <- length(unique(x)) * 2
  rainbow(n, start = 2/6, end = 6/6, alpha = alpha)[seq(n, 0, -2)]
}
cluster_palette_init <- cluster_palette(som_cluster)
bgcol <- cluster_palette_init[som_cluster]

# Plot the SOM map with nodes colored by their cluster assignment.
plot(som_model, type = "codes", bgcol = bgcol, main = "Clusters", codeRendering = "lines")
add.cluster.boundaries(som_model, som_cluster)

# Add the final cluster assignment to our gene information dataframe.
som_model_code_class_cluster <- som_model_code_class
som_model_code_class_cluster$cluster <- som_cluster[som_model_code_class$code_class]
head(som_model_code_class_cluster)

# --- 8. Mapping Properties (e.g., Gene Expression) onto the SOM ---
# This section demonstrates how to visualize the expression of a specific sample across the SOM grid.
# This can be done for any property, such as pathway scores, as long as it's a numeric vector.
color_by_var <- names(data_scale)[4] # Example: Select the 4th sample
color_by <- data_scale[, color_by_var]
unit_colors <- aggregate(color_by, by = list(som_model$unit.classif), FUN = mean, simplify = TRUE)
plot(som_model, type = "property", property = unit_colors[, 2], main = color_by_var, palette.name = coolBlueHotRed)

# --- 9. Loop to Generate Property Maps for All Samples ---
# This loop will create and save a property map for each sample in the dataset.
# Define the color ramp palette.
coolBlueHotRed <- colorRampPalette(c("#1e95d3", "white", "#f65f72")) # Blue-White-Red

# Get all sample names from the column names of the scaled data
sample_names <- colnames(data_scale)

# Loop through each sample name
for (i in seq_along(sample_names)) {
  color_by_var <- sample_names[i]
  color_by <- data_scale[, color_by_var]
  
  # Ensure the property vector length matches the number of genes mapped
  color_by <- color_by[1:length(som_model$unit.classif)]
  
  # Calculate the mean property value for each SOM unit
  unit_colors <- aggregate(color_by, by = list(som_model$unit.classif), FUN = mean, simplify = TRUE)
  
  # Determine the global color range for this specific sample to keep the scale consistent
  global_range <- range(color_by, na.rm = TRUE)
  
  # Export the plot to a PNG file
  png(filename = paste0("SOM_", color_by_var, ".png"), width = 800, height = 800)
  plot(
    som_model,
    type = "property",
    property = unit_colors[, 2],
    main = color_by_var,
    palette.name = coolBlueHotRed,
    zlim = global_range # Fix the color scale to the range of the current sample
  )
  dev.off() # Close the PNG device
}

# --- 10. Export Gene Lists ---
# To export the final table with genes, their mapped SOM node, and final cluster assignment.
head(som_model_code_class_cluster)
# write.csv(som_model_code_class_cluster, "SOM_gene_clusters.csv", row.names = FALSE)
# --- Supp Fig 1c-e: Neonatal Splenic Analysis ---
# Assumes `seurat_neonatal` is the Seurat object for P7 splenic B-1a cells
# UMAP plot
DimPlot(seurat_neonatal, reduction = "umap", label = TRUE)


# CytoTRACE analysis
# Assumes CytoTRACE has been run and scores are in metadata
FeaturePlot(seurat_neonatal, features = "CytoTRACE_score") +
    scale_color_viridis_c()
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
saveRDS(result_sce,"1wk_Cytotrace2_result_sce.rds")
# --- Supp Fig 1f: Neonatal Splenic pseudotime Analysis ---
seurat_object<-readRDS("wt_1wk_alone.rds")
seurat_object$seurat_clusters
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
