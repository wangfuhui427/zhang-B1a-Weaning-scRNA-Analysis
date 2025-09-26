# zhang-B1a-Weaning-scRNA-Analysis
R scripts for the analysis and figure generation of the manuscript "Weaning-driven Dietary Transition Shapes Transcriptome and BCR Repertoire of Mouse Peritoneal B1a Cells: Insights from Single-cell Analysis"
# Code for: Weaning-driven Dietary Transition Shapes Transcriptome and BCR Repertoire of Mouse Peritoneal B1a Cells: Insights from Single-cell Analysis

This repository contains the R scripts used for the analysis and figure generation in our manuscript submitted to *Scientific Data*.

## Data Availability

All raw and processed data have been deposited at NODE with the project ID: OEP00006590 (https://www.biosino.org/node/). 

## Software and Environment

The analyses were performed in R (v4.3.0) on a macOS system. The following key R packages are required to run the scripts:

- Seurat (v5.3.0)
- immunarch (v0.9.1)
- dplyr (v1.1.0)
- ggplot2 (v3.4.0)
- patchwork (v1.1.2)
- and others as listed in the scripts.

## How to Reproduce the Analysis

1.  Download the processed Seurat objects and BCR data from the GEO repository linked above.
2.  Place the data files into a designated `data/` subfolder.
3.  Run the scripts in the following order to reproduce the main figures:
    - `01_scRNA_QC_and_Clustering.R`: Generates Figure 1.
    - `02_Developmental_Trajectory_Analysis.R`: Generates Figure 2.
    - `03_BCR_Repertoire_Analysis.R`: Generates Figure 3 and Supplementary Figure 2.
    - `S1_Bulk_RNAseq_Analysis.R`: Generates Supplementary Figure 1.
    - `S2_Subset-specific BCR characteristics.R`: Generates Supplementary Figure 2.
Please refer to the comments within each script for detailed explanations of each step.
