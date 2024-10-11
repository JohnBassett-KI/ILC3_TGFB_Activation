library('pheatmap')
library('magrittr')
library('Seurat')
#generate figures for figure 4a

load(file = "data/gca_adult_small_intestine/gca_adult_smallInt_Seurat_norm.robj")
DATA <- gcA_adult_smallInt_norm
rm(gcA_adult_smallInt_norm)
load(file = "data/gca_adult_small_intestine/gca_adult_SI_QC_downsamp_final.robj")
#Classify cell types
cellTypes <- unique(gcA_samp$Integrated_05)
#Manually group cells
names(cellTypes) <- c("Lymphoid", "Lymphoid", "Enteric Glia", "Epithelial", "Epithelial",
                    "Lymphoid", "Endothelial", "Lymphoid", "Epithelial", "Epithelial",
                    "Lymphoid", "Lymphoid", "Lymphoid", "Endothelial", "Endothelial",
                    "Endothelial", "Endothelial", "Myeloid", "Myeloid", "Myeloid",
                    "Endothelial", "Endothelial", "Lymphoid", "Epithelial", "Myeloid",
                    "Lymphoid", "Epithelial", "Stem Cell", "Stromal", "Stromal",
                    "Stromal", "Stem Cell", "Lymphoid", "Lymphoid", "Lymphoid",
                    "Stromal", "Lymphoid", "Epithelial", "Endothelial", "Myeloid",
                    "Myeloid", "Lymphoid", "Fibroblast")
cellTypes <- cellTypes[order(names(cellTypes))]
cellTypedf <- data.frame(Cell = cellTypes, Class = names(cellTypes))


#TGFBsig <- c("ITGAV", "ITGB8", "ITGB1", "ITGB3", "ITGB5", "MMP9", "MMP14", "SMAD2", "SMAD3", "SMAD7", "TGFB1", "TGFBR1", "TGFBR2", "SDC4", "TGM2", "FN1", "FBN1", "LTBP1", "LTBP3")  
#pheatobj3<- DATA[c(TGFBsig),]
load(file= "data/TGFBsig.robj")

pheatMat3 <- AverageExpression(DATA, group.by = "Integrated_05",
                               assays = 'originalexp', layer = 'data',
                               features = TGFBsig, return.seurat = T)
HeatMatrix3<- pheatMat3$originalexp$data #average log counts
HeatMatrix3<-HeatMatrix3[, cellTypes]

#Heatmap for all cell types
# All_cell_types <- pheatmap::pheatmap(HeatMatrix3, cluster_rows = T, cluster_cols = T, clustering_method =  "complete",
#                                angle_col = 45, cellwidth = 18, cellheight = 20, fontsize = 12,
#                                gaps_col = c(8,9,16,17,31,37,39),
#                                main = "Small Intestine TGFb Signature")
#Heatmap for selected cell types
selectedCells <- cellTypedf[,"Cell"][c(which(cellTypedf[,"Class"] == "Stromal"),
                      which(cellTypedf[,"Class"] == "Epithelial"))]
selectedCells <- as.character(selectedCells)
selectedCells <- c(selectedCells, c("Activated CD4 T", "Activated CD8 T", "CD8 Tmem",
                     "gdT","Treg", "ILC3"))
selectedCells <- selectedCells[-4]
#save(selectedCells, file="data/gca_adult_small_intestine/SelectedCells_si.robj")



HeatMatrix4 <- HeatMatrix3[, selectedCells]
Selected_cell_types <- pheatmap::pheatmap(HeatMatrix4, cluster_rows = F, cluster_cols = F, color = colorRampPalette(c("#5957A6", "white", "#F26222"))(100),
                                          angle_col = 45, cellwidth = 20, cellheight = 20, fontsize = 14, scale = "row",
                                          gaps_row = seq_along(rownames(HeatMatrix4))[-19],
                                          main = "Small Intestine TGFb Signature")
#VlnPlot(DATA,features = TGFBsig, stack = T, group.by = "Integrated_05") + NoLegend()
