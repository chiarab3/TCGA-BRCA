# Load required libraries
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(DOSE)

# Load DEGs from previous analysis
deg_data <- read.csv("significant_DEGs.csv", row.names = 1)

# Convert gene IDs from Ensembl to Entrez
gene_list_clean <- gsub("\\..*", "", rownames(deg_data))  # Remove Ensembl version numbers
entrez_ids <- mapIds(org.Hs.eg.db, keys = gene_list_clean, column = "ENTREZID",
                     keytype = "ENSEMBL", multiVals = "first")
entrez_ids <- na.omit(entrez_ids)  # Remove NA values
print(length(entrez_ids))

# GO Enrichment Analysis
go_enrich <- enrichGO(gene = entrez_ids,
                      OrgDb = org.Hs.eg.db,
                      keyType = "ENTREZID",
                      ont = "BP",  # "BP" = Biological Process, "MF" = Molecular Function, "CC" = Cellular Component
                      pAdjustMethod = "BH",
                      pvalueCutoff = 0.05,
                      qvalueCutoff = 0.05)

pdf("04_FA.pdf")

# Visualize the top 10 enriched biological processes
barplot(go_enrich, showCategory = 10, title = "Top 10 Enriched GO Terms")
dotplot(go_enrich, showCategory = 10, title = "GO Enrichment Analysis")

# KEGG Pathway Analysis
kegg_enrich <- enrichKEGG(gene = entrez_ids,
                          organism = "hsa",  # Homo sapiens
                          pAdjustMethod = "BH",
                          pvalueCutoff = 0.05)

# Visualizing KEGG pathways
barplot(kegg_enrich, showCategory = 10, title = "Top 10 Enriched KEGG Pathways")
dotplot(kegg_enrich, showCategory = 10, title = "KEGG Pathway Enrichment")

# Create a ranked gene list for GSEA
gene_ranks <- deg_data$log2FoldChange
names(gene_ranks) <- entrez_ids

# Remove missing values
gene_ranks <- na.omit(gene_ranks)

# Ensure unique gene names
gene_ranks <- gene_ranks[!duplicated(names(gene_ranks))]

# Filter out NA values
valid_idx <- !is.na(names(gene_ranks))  # Keep only valid entries
gene_ranks <- gene_ranks[valid_idx]

# Sort genes by log2FoldChange
gene_ranks <- sort(gene_ranks, decreasing = TRUE)

# Perform GSEA
gsea_res <- gseGO(geneList = gene_ranks,
                  OrgDb = org.Hs.eg.db,
                  keyType = "ENTREZID",
                  ont = "BP",
                  pAdjustMethod = "BH",
                  pvalueCutoff = 0.05)

print(gsea_res)
# Visualize GSEA results
gseaplot2(gsea_res, geneSetID = 1, title = gsea_res@result$Description[1])
dotplot(gsea_res, showCategory = 3, title = "Top Enriched Pathways (GSEA)")
dev.off()
