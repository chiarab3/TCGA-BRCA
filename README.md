# RNA-seq Differential Expression Analysis: Breast Cancer Tumor vs. Normal

This repository contains the code and documentation for a differential expression analysis pipeline comparing breast cancer tumor tissue with matched normal tissue using TCGA-BRCA RNA-seq data. The workflow leverages Bioconductor packages to identify genes with altered expression in cancer and explore their functional significance.

## Data Source
The analysis uses RNA-seq data from The Cancer Genome Atlas Breast Cancer (TCGA-BRCA) project, accessed programmatically via the TCGAbiolinks package. Instead of processing raw FASTQ files, it is utilized the pre-computed STAR-Counts provided by the GDC, which enables to focus on downstream analysis while maintaining reproducibility.

## Environment Setup
The analysis is performed in a dedicated conda environment running in WSL (Ubuntu 22.04):

```bash
# Create and activate environment
conda create -n rnaseq python=3.8
conda activate rnaseq

# Configure bioconda channels
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge

# Install core tools
conda install -y r-base
conda install -y bioconductor-deseq2 bioconductor-enhancedvolcano bioconductor-gsea
conda install -y r-ggplot2 r-pheatmap

# Install TCGAbiolinks
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("TCGAbiolinks")
BiocManager::install("SummarizedExperiment")
```

## Analysis Workflow
### 1. Data Acquisition
- Query TCGA-BRCA RNA-seq data
- Select matched tumor-normal pairs
- Download and prepare count matrices

### 2. Exploratory Data Analysis
- Quality control and filtering
- Principal Component Analysis (PCA)
- Sample clustering and visualization

### 3. Differential Expression Analysis
- Normalization using DESeq2
- Statistical modeling of gene expression changes
- Identification of significantly differentially expressed genes

### 4. Functional Analysis
- Pathway enrichment analysis
- Gene set enrichment analysis (GSEA)
- Biological interpretation of results

### 5. Visualization
- Volcano plots of differential expression
- Heatmaps of top differentially expressed genes
- Pathway visualization

## Repository Structure
- `data/`: Contains downloaded and processed data
- `scripts/`: R scripts for each analysis step
- `results/`: Output files and visualizations
- `docs/`: Documentation and notes

## Dependencies
- R >= 4.0
- Bioconductor packages:
  - TCGAbiolinks
  - SummarizedExperiment
  - DESeq2
  - EnhancedVolcano
  - clusterProfiler
- Additional R packages:
  - ggplot2
  - pheatmap
  - RColorBrewer
