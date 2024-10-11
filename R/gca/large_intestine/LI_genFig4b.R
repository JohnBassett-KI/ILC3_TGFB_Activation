library('Seurat')
library('pheatmap')
#generate figures for figure 4b

load(file = "data/gca_adult_large_intestine/gca_adult_largeInt_Seurat_norm.robj")
DATA <- gcA_adult_largeInt_norm
load(file = "data/gca_adult_large_intestine/gca_adult_LI_QC_downsamp_final.robj")
rm(gcA_adult_largeInt_norm)
#Classify cell types
cellTypes <- unique(gcA_samp$Integrated_05)
#Manually group cells
names(cellTypes) <- c("Lymphoid", "Lymphoid", "Enteric Glia", "Epithelial", "Epithelial",
                    "Lymphoid", "Lymphoid", "Epithelial", "Endothelial", "Lymphoid",
                    "Enteroendocrine","Epithelial", "Lymphoid", "Lymphoid", "Lymphoid",
                    "Enteroendocrine","Endothelial", "Endothelial", "Endothelial", "Endothelial",
                    "Myeloid", "Lymphoid", "Lymphoid", "Myeloid", "Myeloid",
                    "Endothelial", "Endothelial", "Lymphoid", "Mesothelial","Epithelial",
                    "Myeloid","Lymphoid","Lymphoid","Epithelial", "Endothelial",
                    "Lymphoid", "Lymphoid","Stem Cell", "Stromal", "Stromal",
                    "Stromal", "Stromal", "Fibroblast", "Stem Cell", "Lymphoid",
                    "Lymphoid", "Lymphoid", "Lymphoid", "Stromal", "Lymphoid",
                    "Epithelial", "Endothelial", "Myeloid","Myeloid", "Lymphoid",
                    "Stromal","Fibroblast", "Fibroblast")
rm(gcA_samp)
cellTypes <- cellTypes[order(names(cellTypes))]
cellTypedf <- data.frame(Cell = cellTypes, Class = names(cellTypes))


TGFBsig <- c("ITGAV", "ITGB1", "ITGB3", "ITGB5","ITGB8", "MMP9", "MMP14", "SMAD2", "SMAD3", "SMAD7", "TGFB1", "TGFBR1", "TGFBR2", "SDC4", "TGM2", "FN1", "FBN1", "LTBP1", "LTBP3")  
save(TGFBsig, file = "data/TGFBsig.robj")
load(file = "data/TGFBsig.robj")
pheatMat2 <- AverageExpression(DATA, group.by = "Integrated_05",
                               assays = 'originalexp', layer = 'data',
                               features = TGFBsig, return.seurat = T)
HeatMatrix2<- pheatMat2$originalexp$data #average log counts
HeatMatrix2<-HeatMatrix2[, cellTypes]

#Heatmap for all cell types
# All_cell_types <- pheatmap::pheatmap(HeatMatrix3, cluster_rows = T, cluster_cols = T, clustering_method =  "complete",
#                                angle_col = 45, cellwidth = 18, cellheight = 20, fontsize = 12,
#                                #gaps_col = c(9,10,12,19,22,43,44,50,52),
#                                main = "Large Intestine TGFb Signature")
#Heatmap for selected cell types
selectedCells <- cellTypedf[,"Cell"][c(which(cellTypedf[,"Class"] == "Stromal"),
                      which(cellTypedf[,"Class"] == "Epithelial"))]
selectedCells <- as.character(selectedCells)
selectedCells <- c(selectedCells, c("Activated CD4 T", "Activated CD8 T", "CD8 Tmem",
                                    "Th17", "gdT","Treg","NK cell","ILC3"))

selectedCells <- selectedCells[-5] 
selectedCells<- selectedCells[-5]
save(selectedCells, file="data/gca_adult_large_intestine/selectedCellTypes.robj")

HeatMatrix4 <- HeatMatrix2[, selectedCells]
Selected_cell_types <- pheatmap::pheatmap(HeatMatrix4, cluster_rows = F, cluster_cols = F, color = colorRampPalette(c("#5A619D", "white", "#F15F22"))(100),
                                          angle_col = 45, cellwidth = 20, cellheight = 20, fontsize = 12, scale = "row",
                                          gaps_row = seq_along(rownames(HeatMatrix4))[-19],
                                          main = "Large Intestine TGFb Signature")

