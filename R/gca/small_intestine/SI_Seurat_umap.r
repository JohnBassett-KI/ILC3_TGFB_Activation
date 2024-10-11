library('Seurat')
library('harmony')
library('magrittr')
load("data/gca_adult_small_intestine/gca_adult_smallIntQC_downsampled.robj")
gcA_samp <- gca_adult_smallIntQC_downsampled
rm(gca_adult_smallIntQC_downsampled)
#gcA_samp <- RenameAssays(object = gcA_samp, originalexp = "RNA")
gcA_samp <- NormalizeData(gcA_samp, normalization.method ="LogNormalize") %>%
  FindVariableFeatures() %>%
  ScaleData(vars.to.regress= 'pct_counts_mt') %>%
  RunPCA(verbose = FALSE) #TODO regressions in scale data
gcA_samp <- RunHarmony(gcA_samp, group.by.vars = c("batch", "Fraction"))
#Visualization and Clustering
# These are now standard steps in the Seurat workflow for visualization and clustering
gcA_samp <- RunUMAP(gcA_samp, dims = 1:30,
                    verbose = FALSE, reduction = "harmony")

gcA_samp <- FindNeighbors(gcA_samp, dims = 1:30,
                          verbose = FALSE, reduction = "harmony")
gcA_samp <- FindClusters(gcA_samp, verbose = FALSE)
#FeaturePlot(gcA_samp, features = c("ADAR")) 
DimPlot(gcA_samp, label = T, group.by = "Integrated_05",
        repel = T, label.size = 5, shuffle = T) + NoLegend()
#Highlight cell types
Idents(gcA_samp) <- gcA_samp$Integrated_05
HighlightCells = WhichCells(gcA_samp, idents = "ILC3")
DimPlot(gcA_samp, label = T, group.by = "Integrated_05",
        repel = T, label.size = 3, shuffle = T, cells.highlight = HighlightCells) + NoLegend()
#remove all ILCs which do not cluster in cluster 24
Idents(gcA_samp)<- gcA_samp$Integrated_05
namesILC3 <- WhichCells(gcA_samp, idents= "ILC3")
Idents(gcA_samp) <- gcA_samp$seurat_clusters
namesCluster24 <- WhichCells(gcA_samp, idents= "24")
removeUNCLUST <- namesILC3[which(!namesILC3 %in% namesCluster24)]

gcA_samp <- subset(gcA_samp, cells = removeUNCLUST, invert = T)

# FeaturePlot(gcA_samp, features = c("CD2", "CD58",  "HLA-B", "KIR3DL2"), cells = HighlightCells)
# FeatureScatter(gcA_samp, feature1 = 'total_counts', feature2 = 'n_genes_by_counts',
#                group.by = 'Integrated_05', pt.size = 0.5)
# FeatureScatter(gcA_samp, feature1 = 'total_counts', feature2 = 'pct_counts_mt',
#                group.by = 'Integrated_05', pt.size = 0.5)
# VlnPlot(gcA_samp, features = 'pct_counts_mt' , pt.size = 0, group.by ='Integrated_05', sort = T ) + NoLegend()

#save(gcA_samp, file = "data/gca_adult_small_intestine/gca_adult_SI_QC_downsamp_final.robj")
writeLines(Cells(gcA_samp), con = "data/gca_adult_small_intestine/cell_names.txt")

