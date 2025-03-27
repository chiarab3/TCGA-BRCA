# Load required libraries
library(TCGAbiolinks)
library(SummarizedExperiment)
library(GenomicRanges)
library(DESeq2)
library(pheatmap) 
library(ggplot2)  

# Load the dataset
load("tcga_brca_subset.RData")

# Basic dataset information
print("Dataset overview:")
print(data)
print(dim(data))  # Dimensions: genes x samples
print("Sample metadata:")
print(colData(data)[, c("patient", "sample_type", "shortLetterCode")])

# Inspect gene annotation
print("First few genes:")
print(head(rowData(data)))

# Check raw counts distribution
counts <- assay(data)  # Extract raw counts
mean_counts <- rowMeans(counts)  # Compute mean expression per gene

# Histogram of mean counts (log10 transformed)
pdf("02_EDA.pdf")
hist(log10(mean_counts + 1), main="Distribution of Mean Counts", 
     xlab="log10(mean counts + 1)", col="lightblue", border="black")

# Sample correlation heatmap
sample_cor <- cor(counts)
pheatmap(sample_cor, main="Sample Correlation Heatmap", clustering_method="ward.D2")

# Filter low-expression genes (mean counts < 10)
keep_genes <- mean_counts > 10
data_filtered <- data[keep_genes, ]
print("Number of genes after filtering:")
print(dim(data_filtered))

# Principal Component Analysis (PCA)
dds <- DESeqDataSetFromMatrix(countData = assay(data_filtered), 
                              colData = colData(data_filtered), 
                              design = ~ sample_type)  # Sample type as grouping variable

vsd <- vst(dds, blind=TRUE)  # Variance Stabilizing Transformation

pca_data <- plotPCA(vsd, intgroup="sample_type", returnData=TRUE)
percentVar <- round(100 * attr(pca_data, "percentVar"))

ggplot(pca_data, aes(PC1, PC2, color=sample_type)) +
  geom_point(size=4) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  ggtitle("PCA of TCGA-BRCA Samples") +
  theme_minimal()
dev.off()

# Save the filtered dataset for differential expression analysis
save(data_filtered, file = "tcga_brca_filtered.RData")
