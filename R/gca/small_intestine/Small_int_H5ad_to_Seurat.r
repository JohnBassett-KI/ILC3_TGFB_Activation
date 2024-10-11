#gut cell atlas data (gca) has been processed on a HPC and reduced to only data 
#pertaining to the small intestine. The gca data is originally processed
#as an AnnData object in python and saved in the .h5ad format. This file
#converts the Anndata object to a seurat object and saves it in the .robj format.

BiocManager::install("zellkonverter")
BiocManager::install("SingleCellExperiment")
library("zellkonverter")
library("Seurat")

sce <- zellkonverter::readH5AD("data/gca_adult_small_intestine/adult_smallInt.h5ad")
gca_adult_smallInt <- as.Seurat(sce, counts = "X", data = NULL)
save(gca_adult_smallInt, file = "data/gca_adult_small_intestine/gca_SI_Seurat.robj")


sce <- readH5AD("data/gca_adult_small_intestine/gca_adult_smallInt_norm.h5ad")
gcA_adult_smallInt_norm <- as.Seurat(sce, counts = NULL, data = "X")
save(gcA_adult_smallInt_norm, file = "data/gca_adult_small_intestine/gca_adult_smallInt_Seurat_norm.robj")

