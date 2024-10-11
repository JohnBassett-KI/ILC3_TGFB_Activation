#before running:
#module load R_packages/4.1.1
#module load Rstudio/2022.07.1-554
if (!requireNamespace("BiocManager", quietly=TRUE)) {
  install.packages("BiocManager")
}
BiocManager::install("zellkonverter")
library('zellkonverter')
library('Seurat')
#Raw Data
sce <- readH5AD("data/gca_adult_large_intestine/gca_adult_largeInt.h5ad")

gcA_adult_largeInt <- as.Seurat(sce, counts = "X", data = NULL)

save(gcA_adult_largeInt, file = "data/gca_adult_large_intestine/gca_adult_largeInt_Seurat.robj")


#Normalized data
sce <- readH5AD("data/gca_adult_large_intestine/gca_adult_largeInt_norm.h5ad")

gcA_adult_largeInt_norm <- as.Seurat(sce, counts = NULL, data = "X")

save(gcA_adult_largeInt_norm, file = "data/gca_adult_large_intestine/gca_adult_largeInt_Seurat_norm.robj")
