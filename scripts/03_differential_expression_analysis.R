# Load required libraries
library(DESeq2)
library(TCGAbiolinks)
library(SummarizedExperiment)
library(ggplot2)
library(EnhancedVolcano)

# Load filtered data
load("tcga_brca_filtered.RData")

# Prepare DESeq2 dataset
dds <- DESeqDataSetFromMatrix(
  countData = assay(data_filtered),
  colData = colData(data_filtered),
  design = ~ sample_type)

# Normalization and pre-processing
dds <- DESeq(dds)  # Run DESeq2 pipeline
normalized_counts <- counts(dds, normalized=TRUE)  # Extract normalized counts
write.csv(normalized_counts, "normalized_counts.csv")  # Save normalized counts

# Extract Differential Expression Results
res <- results(dds, contrast=c("sample_type", "Primary Tumor", "Solid Tissue Normal"))
# Reduction of variability in log2FoldChange estimates for better stability in the results
res <- lfcShrink(dds, coef=2, type="apeglm")

# Order results by adjusted p-value
res <- res[order(res$padj), ]
# Save results to file
write.csv(as.data.frame(res), "differential_expression_results.csv")

# Summary of DE analysis
print(summary(res))

# Identify significantly differentially expressed genes (DEGs)
# Apply cutoff for significance (padj < 0.05 and log2FoldChange > 1 or < -1)
sig_genes <- subset(res, padj < 0.05 & abs(log2FoldChange) > 1)
print(dim(sig_genes))
head(sig_genes)
write.csv(as.data.frame(sig_genes), "significant_DEGs.csv")

# Volcano plot to visualize DEGs
pdf("03_DEA.pdf")
EnhancedVolcano(res,
  lab = rownames(res),
  x = "log2FoldChange",
  y = "padj",
  title = "Differential Expression Analysis",
  pCutoff = 0.05,
  FCcutoff = 1,
  pointSize = 2.0,
  labSize = 3.0)
dev.off()
