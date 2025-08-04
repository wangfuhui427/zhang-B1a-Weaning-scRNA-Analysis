# zhang-B1a-Weaning-scRNA-Analysis
R scripts for the analysis and figure generation of the manuscript "A single-cell transcriptomic and B-cell receptor repertoire dataset of mouse peritoneal B-1a cells across the weaning transition"
# Code for: A single-cell transcriptomic and B-cell receptor repertoire dataset of mouse peritoneal B-1a cells across the weaning transition

This repository contains the R scripts used for the analysis and figure generation in our manuscript submitted to *Scientific Data*.

## Data Availability

The raw and processed data are available at the Gene Expression Omnibus (GEO) under accession number: [请在这里填写您的GEO登录号].

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
    - `S1_Bulk_RNAseq_Analysis.R`: Generates Supplementary Figure 1a-b.

Please refer to the comments within each script for detailed explanations of each step.
