library("Seurat")
library("pheatmap")
library("ggplot2")
library("scales")
library('dplyr')
#pIBD data

pIBD <- readRDS(file = "data/pIBD/combatCorr_pIBD.rds")

onlyILC3s <- subset(pIBD, subset = annot == "ILC")
  
#Signature genes
#TGFBsig <- c("ITGAV", "ITGB8", "ITGB3", "MMP9", "MMP14", "SMAD2", "SMAD3", "SMAD7", "TGFB1", "TGFBR1", "TGFBR2", "SDC4", "TGM2", "FN1", "FBN1", "LTBP1", "LTBP3")  
load(file = "data/TGFBsig.robj")
Idents(pIBD) <- pIBD$annot

#HeatMap
pheatMat2 <- AverageExpression(pIBD, group.by = "annot",
                               assays = 'comBat', layer = 'data',
                               features = TGFBsig, return.seurat = T)
HeatMatrix2<-pheatMat2@assays$comBat$data
HeatMatrix2 <- HeatMatrix2[,c(1:8,12,10,11,9)]
TC_heat2 <- pheatmap::pheatmap(HeatMatrix2, cluster_rows = F, cluster_cols = F, color = colorRampPalette(c("#5A619D", "white", "#F15F22"))(100),
                               angle_col = 45, cellwidth = 20, cellheight = 20, fontsize = 12, scale = "row",
                               gaps_row = seq_along(rownames(HeatMatrix2)[-17]),
                               main = "pIBD TGFb Signature")





#Violin Plot
# Idents(onlyILC3s) <- onlyILC3s$Tissue
# #TGFBsig1 <- c("SMAD2", "SMAD3", "SMAD7", "TGFB1", "TGFBR1", "TGFBR2")
# TGFBsig2 <- c("SDC4", "TGM2", "FN1", "FBN1", "LTBP1","ITGAV", "ITGB8", "MMP9", "MMP14")
#VlnPlot(pIBD, features = c("MMP9", "SDC4", "TGM2", "ITGAV", "MMP14", "ITGB3"),raster = F, add.noise = F, log = F, flip = F, stack = F, group.by = "annot", assay = "comBat", layer = "data") + NoLegend()
# #VlnPlot(onlyILC3s, features = TGFBsig1,raster = T, add.noise = F, stack = T, log = T) + NoLegend()
# VlnPlot(onlyILC3s, features = TGFBsig2,raster = F, add.noise = F, stack  = F, log = F, assay = "RNA", layer = "counts") + NoLegend()
# VlnPlot(pIBD, features = "MMP9", stack = F, log = F,  assay = "comBat")

#DotPlot

Idents(onlyILC3s) <- onlyILC3s$Tissue
pct_exp <- DotPlot(onlyILC3s, features = TGFBsig, scale = F, idents = c("Noninflamed", "Inflamed"), cols = c("#5957A6","#F26222"), scale.by = "radius", assay = "comBat") 
pct_exp  
#DotPlot(pIBD, features = TGFBsig, scale = T)

#Dot plot requires percent expressing data to come from RNA layer and Average expression to come from Combat integrated data

percent_expr <- DotPlot(onlyILC3s, features = TGFBsig, scale = F, idents = c("Noninflamed", "Inflamed"), cols = c("#5957A6","#F26222"), scale.by = "radius", assay = "RNA")$data %>%
  dplyr::select(features.plot, id, pct.exp)

expr_integrated <- DotPlot(onlyILC3s, features = TGFBsig, scale = F, idents = c("Noninflamed", "Inflamed"), cols = c("#5957A6","#F26222"), scale.by = "radius", assay = "comBat")$data %>%
  dplyr::select(features.plot, id, avg.exp)
# Merge the percent expression and average expression data
dot_data <- merge(percent_expr, expr_integrated, by = c("features.plot", "id"))

# Plot the dot plot
ggplot(dot_data, aes(x = id, y = features.plot)) +
  geom_point(aes(size = pct.exp, color = avg.exp)) +
  scale_color_gradient(low = "#5957A6", high = "#F26222") +  # Adjust color scheme as needed
  scale_size_continuous(range = c(0.5, 5), name = "Percent Expressed") +  # Adjust the size range
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid = element_blank(),
        axis.line = element_line(),
        axis.ticks = element_line(color = "black")) +
  labs(title = "ILC3: Inflamed v. Non-Inflamed Tissue",
       x = "Clusters", y = "Genes", color = "Avg. Expression", size = "Percent Expressed")+
  coord_flip()

