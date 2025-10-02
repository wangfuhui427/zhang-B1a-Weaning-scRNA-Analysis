# Load necessary libraries
library(Seurat)
library(ggplot2)
library(dplyr)

# --- Quality Control Violin Plots ---
seurat_object[["percent.mt"]] <- PercentageFeatureSet(seurat_object, pattern = "^mt-")
seurat_object[["percent.ribo"]] <- PercentageFeatureSet(seurat_object, pattern = "^Rps|^Rpl")
ery_genes <- c("Hbb-bt", "Hba-a1", "Hba-a2","Hbb-bs","Klf1","Tal1")
mye_genes<-c("Cd68","Ly6c","Ly6g","Csf1r","Cd34")
T_genes<-c("Cd3e","Cd8a","Cd4")
CLP_genes<-c("Il7r","Cd19","Cd79a","Cd79b","Vpreb1","Dntt","Rag1","Rag2")
Idents(seurat_object)<-"Sample"
ery_genes <- ery_genes[ery_genes %in% rownames(seurat_object)]
mye_genes <- mye_genes[mye_genes %in% rownames(seurat_object)]
T_genes <- T_genes[T_genes %in% rownames(seurat_object)]
CLP_genes <- CLP_genes[CLP_genes %in% rownames(seurat_object)]
seurat_object[["ery_expr"]] <- Matrix::colSums(seurat_object@assays$RNA@data[ery_genes, , drop = FALSE])
seurat_object[["mye_expr"]] <- Matrix::colSums(seurat_object@assays$RNA@data[mye_genes, , drop = FALSE])
seurat_object[["T_expr"]] <- Matrix::colSums(seurat_object@assays$RNA@data[T_genes, , drop = FALSE])
seurat_object[["CLP_expr"]] <- Matrix::colSums(seurat_object@assays$RNA@data[CLP_genes, , drop = FALSE])
VlnPlot(seurat_object, features = c("mye_expr", "T_expr","ery_expr"), ncol = 3,raster=FALSE)
VlnPlot(seurat_object_peritoneal, 
        features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
        group.by = "sample_id", 
        pt.size = 0, 
        ncol = 3)

# --- CCA and UMAP of Transcriptional Clusters ---
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
DimPlot(seurat_object_peritoneal, 
        reduction = "umap", 
        group.by = "seurat_clusters", 
        label = TRUE, 
        label.size = 5) + NoLegend()

# --- Marker Gene Dot Plot ---
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
