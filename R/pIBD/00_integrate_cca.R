#set up environment
###################
library('Seurat')
library('harmony')
library('dplyr')
###################
pIBD <- readRDS(file =  "data/pIBD/pIBD.RDS")

#spit the data into experiment based layers
pIBD[["RNA"]] <- split(pIBD[["RNA"]], f = pIBD$experiment)


# Check factors for integration
pIBD <- NormalizeData(pIBD)
pIBD <- FindVariableFeatures(pIBD)
pIBD <- ScaleData(pIBD)
pIBD <- RunPCA(pIBD)
pIBD <- FindNeighbors(pIBD, dims = 1:30, reduction = "pca")
pIBD <- FindClusters(pIBD, resolution = 2, cluster.name = "unintegrated_clusters")
pIBD <- RunUMAP(pIBD, dims = 1:30, reduction = "pca", reduction.name = "umap.unintegrated")

#check inflammation based differences: Very Few
DimPlot(pIBD, reduction = "umap.unintegrated", group.by = c("Tissue", "donor"))
#check for donor based diff
DimPlot(pIBD, reduction = "umap.unintegrated", group.by = c("donor", "seurat_clusters"))
#cell cycle phase: Not substantial
DimPlot(pIBD, reduction = "umap.unintegrated", group.by = c("phase", "seurat_clusters"))
#10x experiment based differences: FLAG for Integration
DimPlot(pIBD, reduction = "umap.unintegrated", group.by = c("experiment", "seurat_clusters"))
#check to see if donor and experiment are same: experiment integration accounts for Donor based differences
DimPlot(pIBD, reduction = "umap.unintegrated", group.by = c("experiment", "donor"))

pIBD <- IntegrateLayers(object = pIBD, method = CCAIntegration, orig.reduction = "pca", new.reduction = "cca",
                        verbose = FALSE)

# re-join layers after integration
pIBD[["RNA"]] <- JoinLayers(pIBD[["RNA"]])

pIBD <- FindNeighbors(pIBD, reduction = "cca", dims = 1:30)
pIBD <- FindClusters(pIBD, resolution = 0.7)
pIBD <- RunUMAP(pIBD, dims = 1:30, reduction = "cca")
# Visualization
DimPlot(pIBD, reduction = "umap", group.by = c("donor", "experiment"))
DimPlot(pIBD, reduction = "umap", group.by = c("Tissue", "annot"))

#save integrated object
saveRDS(pIBD, file="data/pIBD/batchCorr_PIBD.rds")
