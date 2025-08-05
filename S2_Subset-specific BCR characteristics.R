# --- Supp Fig 2a: VH-VL Pairing Heatmap ---
library(ggplot2)
library(dplyr)
library(tidyr)
library(cowplot)
library(RColorBrewer)

my_palette <- colorRampPalette(c("#FFDEAD", "darkred"))
cell_totals <- c(WT_3WK = 15206, WT_5WK = 16272, WT_7WK = 13381, WT_10WK = 16020)

bcr_ratio_data <- top20_combos_per_sample %>%
separate(IGHV_IGKV, into = c("IGHV", "IGKV"), sep = "\\+", remove = FALSE) %>%
mutate(Ratio = TotalCount  / cell_totals[Sample])



global_max <- max(bcr_ratio_data$Ratio)
plot_bcr_heatmap <- function(df, sample_name, scale_max) {
df_sub <- df %>%
         filter(Sample == sample_name)
    
       # 
       p_main <- ggplot(df_sub, aes(x = IGKV, y = IGHV, fill = Ratio)) +
           geom_tile(color = "white") +
          scale_fill_gradientn(colors = my_palette(100), limits = c(0, scale_max), name = "Ratio") +
          theme_minimal(base_size = 9) +
          theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
                             axis.text.y = element_text(size = 7),
                             axis.title = element_blank(),
                            panel.grid = element_blank(),
                             legend.position = "none")
       
         # up
         p_top <- df_sub %>%
             group_by(IGKV) %>%
             summarise(total = sum(Ratio)) %>%
             ggplot(aes(x = IGKV, y = total)) +
             geom_bar(stat = "identity", fill = "#98C6EA") +
             theme_minimal(base_size = 7) +
             theme(axis.text.x = element_blank(),
                               axis.title = element_blank(),
                              axis.ticks.x = element_blank(),
                               panel.grid = element_blank())
         
           # right
           p_right <- df_sub %>%
               group_by(IGHV) %>%
               summarise(total = sum(Ratio)) %>%
               ggplot(aes(x = total, y = IGHV)) +
               geom_bar(stat = "identity", fill = "#98C6EA") +
               theme_minimal(base_size = 7) +
               theme(axis.text.y = element_blank(),
                                 axis.title = element_blank(),
                                 axis.ticks.y = element_blank(),
                                 panel.grid = element_blank())
           
             # combine
             heatmap_plot <- ggdraw() +
                 draw_plot(p_top, x = 0.21, y = 0.8, width = 0.59, height = 0.2) +
                 draw_plot(p_main, x = 0.1, y = 0, width = 0.7, height = 0.8) +
                 draw_plot(p_right, x = 0.8, y = 0.15, width = 0.2, height = 0.64) +
                 draw_plot_label(label = sample_name, x = 0.1, y = 1, size = 4)
            
              return(heatmap_plot)
           }
p1 <- plot_bcr_heatmap(bcr_ratio_data, "WT_3WK", global_max)
p2 <- plot_bcr_heatmap(bcr_ratio_data, "WT_5WK", global_max)
p3 <- plot_bcr_heatmap(bcr_ratio_data, "WT_7WK", global_max)
p4 <- plot_bcr_heatmap(bcr_ratio_data, "WT_10WK", global_max)
plot_grid(p1, p2, p3, p4, ncol = 2)

# --- Supp Fig 2b: Volcano Plot of DEGs ---
# Requires differential expression results between Ighv11-2 and Ighv12-3 cells
# A package like EnhancedVolcano is excellent for this
library(EnhancedVolcano)
EnhancedVolcano(deg_results,
                lab = rownames(deg_results),
                x = 'avg_log2FC',
                y = 'p_val_adj')
library(ggplot2)
library(dplyr)

df <- df %>%
  mutate(delta_pct = (pct.1 - pct.2)*100)

df <- df %>%
  mutate(Regulation = case_when(
    avg_log2FC >= 0.25  & p_val < 0.05 ~ "11-2_Up",
    avg_log2FC <= -0.25  & p_val < 0.05 ~ "12-3_Up",
    TRUE ~ "None"
  ))


ggplot(df, aes(x = delta_pct, y = avg_log2FC, color = Regulation)) +
  geom_point(alpha = 0.8, size = 1) +
  scale_color_manual(values = c("11-2_Up" = "#E64B35", "12-3_Up" = "#1e95d3", "None" = "gray")) +
  geom_vline(xintercept = c(-5, 5), linetype = "dashed", color = "black") +
  geom_hline(yintercept = c(-0.25, 0.25), linetype = "dashed", color = "black") +
  theme_bw(base_size = 10) +
  labs(x = "Delta Percent (pct.1 - pct.2)",
       y = "Average log2 Fold Change",
       title = "Differential Expression: log2FC vs Delta pct")
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

