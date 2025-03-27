library(TCGAbiolinks)
   
# Query TCGA for RNA-seq data
query <- GDCquery(
  project = "TCGA-BRCA",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts",
  sample.type = c("Primary Tumor", "Solid Tissue Normal"), 
  barcode = c("TCGA-A7-A0CE", "TCGA-A7-A0CE-11A", "TCGA-A7-A0D9", "TCGA-A7-A0D9-11A", 
              "TCGA-BH-A0B3", "TCGA-BH-A0B3-11A", "TCGA-BH-A0B8", "TCGA-BH-A0B8-11A",
              "TCGA-BH-A0C0", "TCGA-BH-A0C0-11A")  # 5 tumor-normal pairs)
   
# Download data
GDCdownload(query, method = "api", files.per.chunk = 5)
   
# Prepare the data
data <- GDCprepare(query)
   
# Save to file
save(data, file = "tcga_brca_subset.RData")

# Print summary to verify data
print(dim(data))
print(head(colnames(data)))
