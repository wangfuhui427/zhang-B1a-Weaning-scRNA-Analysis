# Load necessary libraries
library(immunarch)
library(dplyr)
library(ggplot2)
library(Seurat)

# --- Data Preprocessing Step ---
# Assumption: Your BD Rhapsody VDJ output file is named "bcr_annotations.csv"
# and it contains data for all samples, with a column named 'sample_id' to differentiate them.
# If each sample is in a separate file, you will need to read them individually and combine them using `rbind`.

# 1. Read the VDJ data from BD Rhapsody
# Note: Choose the correct function (e.g., read.csv or read.tsv) based on your file's delimiter.
raw_bcr_data <- read.csv("path/to/your/bcr_annotations.csv", stringsAsFactors = FALSE)

# 2. "Translate" BD Rhapsody column names to the standard immunarch format
# This is an example; you must adjust it to match the actual column names in your file.
# Core immunarch columns: Clones, CDR3.nt, CDR3.aa, V.name, D.name, J.name
# Other important columns: Barcode, ContigID, Sample
bcr_data_formatted <- raw_bcr_data %>%
  rename(
    Clones = Clone.ID,           # 'Clones' is required by immunarch for clone counting
    CDR3.aa = CDR3.Amino.Acid,   # CDR3 amino acid sequence
    CDR3.nt = CDR3.Nucleotide,   # CDR3 nucleotide sequence
    V.name = V.Gene,             # V-gene segment
    D.name = D.Gene,             # D-gene segment (if present)
    J.name = J.Gene,             # J-gene segment
    Barcode = Cell.Label,        # Cell barcode
    Sample = sample_id           # Sample identifier
  ) %>%
  # Ensure only productive, paired heavy and light chain sequences are retained
  filter(Is.Productive == TRUE, Chain.Type %in% c("IGH", "IGK", "IGL"))

# For simplicity, we often start by analyzing only the heavy chains (IGH) for V-gene usage, etc.
igh_data <- bcr_data_formatted %>% filter(Chain.Type == "IGH")

# --- Figure 3a: Ighv Gene Usage Pie Charts ---
# 3. Create an immunarch object
# We will load directly from our pre-processed data frame
immdata_list <- list("your_data_name" = igh_data)
immdata <- repLoad(immdata_list)

# 4. Generate V-gene usage statistics
# immunarch typically recognizes the "V.name" column
gene_usage <- geneUsage(immdata$data, .gene = "V", .norm = TRUE)

# 5. Plot the pie chart
# The `.plot = "pie"` option may not be directly supported in newer immunarch versions,
# so we will use ggplot2 for more control.
# First, get the data
gene_usage_df <- gene_usage %>% select(Gene, Freq)
# Plot
ggplot(gene_usage_df, aes(x = "", y = Freq, fill = Gene)) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar("y", start = 0) +
  theme_void() + 
  labs(fill = "IGHV Gene")

# --- Figure 3b: UMAP of Dominant Clones ---
# 6. Get the cell barcodes for a specific V-gene
# Ensure the barcode format matches the one in your Seurat object exactly!
ighv11_2_barcodes <- igh_data %>% 
  filter(V.name == "IGHV11-2") %>% 
  pull(Barcode) %>%
  unique() # Get each cell's barcode only once

# 7. Highlight these cells in the Seurat object
DimPlot(seurat_object_peritoneal, 
        cells.highlight = ighv11_2_barcodes, 
        cols.highlight = "red",
        sizes.highlight = 1) # Make the highlighted points more visible

# --- Figure 3c & 3d: D50 Diversity Boxplots ---
# 8. Calculate diversity by sample
# The `.col = "Sample"` argument tells immunarch to group the data by the 'Sample' column
div_by_time <- repDiversity(immdata$data, .method = "d50", .col = "Sample")

# 9. Plot the boxplot
ggplot(div_by_time, aes(x = Sample, y = d50)) + 
    geom_boxplot() +
    geom_jitter(width = 0.1, alpha = 0.7) + # Add jittered points to show data distribution
    theme_classic() +
    xlab("Timepoint") +
    ylab("D50 Index")

# To calculate diversity by cluster, you first need to merge Seurat cluster information into the igh_data dataframe.
# cluster_info <- seurat_object_peritoneal@meta.data %>% select(Barcode = "cell_id_col", seurat_clusters)
# igh_data_with_clusters <- left_join(igh_data, cluster_info, by = "Barcode")
# Then, load and calculate diversity using a similar approach.

# --- Figure 3e: Specific Clonotype Tracking ---
# 10. Track the dynamics of specific CDR3 amino acid sequences (CTaa)
# Ideally, a unique clonotype is defined by both heavy and light chain CDR3s.
# Here, we assume you can find unique CDR3.aa combinations via Clone.ID or paired barcodes.
# For this simplified example, we continue using the heavy chain CDR3.aa.
clonotype_dynamics <- igh_data %>%
  filter(CDR3.aa %in% c("CMRYSNYWYFDVW", "CMRYGNYWYFDVW")) %>% # Example simplified CDR3 sequences
  group_by(Sample, CDR3.aa) %>%
  summarise(count = n()) %>%
  # Calculate frequency
  left_join(igh_data %>% group_by(Sample) %>% summarise(total_count = n()), by = "Sample") %>%
  mutate(frequency = count / total_count * 100)

# 11. Plot the bar chart
ggplot(clonotype_dynamics, aes(x = Sample, y = frequency, fill = CDR3.aa)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_classic() +
  xlab("Timepoint") +
  ylab("Clonotype Frequency (%)")
