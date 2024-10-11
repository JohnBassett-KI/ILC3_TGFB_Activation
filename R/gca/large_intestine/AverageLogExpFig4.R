#generate figures for figure 4b


#load(file = "data/gca_adult_large_intestine/gca_adult_LI_QC_downsamp_final.robj")
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
cellTypes <- cellTypes[order(names(cellTypes))]
cellTypedf <- data.frame(Cell = cellTypes, Class = names(cellTypes))


TGFBsig <- c("ITGAV", "ITGB8", "ITGB3", "MMP9", "MMP14", "SMAD2", "SMAD3", "SMAD7", "TGFB1", "TGFBR1", "TGFBR2", "SDC4", "TGM2", "FN1", "FBN1", "LTBP1", "LTBP3")  
Idents(gcA_adult_largeInt_norm) <- gcA_adult_largeInt_norm$Integrated_05
pheatobj3<- subset(gcA_adult_largeInt_norm, features = TGFBsig)

pheatMat3<- AverageExpression(pheatobj3, return.seurat = FALSE)
HeatMatrix3 <- pheatMat3$originalexp
HeatMatrix3<-HeatMatrix3[, cellTypes]

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
selectedCells
selectedCells <- selectedCells[-5] 
selectedCells<- selectedCells[-5]


HeatMatrix4 <- HeatMatrix3[, selectedCells]
Selected_cell_types <- pheatmap::pheatmap(HeatMatrix4, cluster_rows = F, cluster_cols = F, color = colorRampPalette(c("#5A619D", "white", "#F15F22"))(100),
                                          angle_col = 45, cellwidth = 20, cellheight = 20, fontsize = 12, scale = "row",
                                          gaps_row = seq_along(rownames(HeatMatrix4))[-17],
                                          main = "Large Intestine TGFb Signature")
