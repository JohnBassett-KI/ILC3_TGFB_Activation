library('Seurat')
load("data/gca_adult_small_intestine/gca_adult_smallInt_Seurat_QC.robj")
#gca_final is the seurat object corresponding to the final qc output
gca_finalQC <- gca_final #rename to gca_finalQC for clarity
rm(gca_final)

# ##########
# #CHECK QC#
# ##########
# 
#Cell frequencies after QC
Cells <- table(gca_finalQC@meta.data$Integrated_05)
order <- sort(Cells, decreasing = T)
names(order)

#Sort cells according to their pct_counts_mt
SortMT <- cbind(as.character(gca_finalQC@meta.data$Integrated_05), gca_finalQC@meta.data$pct_counts_mt)
sortLevels <- levels(factor(SortMT[,1]))
medOUT <- NULL
for(l in sortLevels){
  CellMT <-  SortMT[SortMT[,1] == l,]
  if(!is.null(dim(CellMT))){
    res <- median(as.numeric(CellMT[,2]))
    MTmed <- cbind(CellMT[1,1], res)
  }else{
    res <- as.numeric(CellMT[2])
    MTmed <- cbind(CellMT[1], res)
  }
  if(is.null(medOUT)){
    medOUT <- MTmed
  }else{
    medOUT <- rbind(medOUT,MTmed)
  }
}
medOUT <- as.data.frame(medOUT)
rownames(medOUT) <- medOUT[,1]
medOUT[,2] <- as.numeric(medOUT[,2])
medOUT <- medOUT[order(medOUT[,'res'], decreasing = T),]
f_order <- rownames(medOUT)
empties <- levels(gca_finalQC@meta.data$Integrated_05)[!(levels(gca_finalQC@meta.data$Integrated_05) %in% f_order)]
f_order <- c(f_order,empties)
gca_finalQC@meta.data$Integrated_05 <- factor(gca_finalQC@meta.data$Integrated_05, levels = f_order)


#Visualize Mitochondrial QC
VlnPlot(gca_finalQC,features=c('pct_counts_mt'), pt.size = 0, group.by = 'Integrated_05', adjust = 2) + NoLegend()
VlnPlot(gca_finalQC,features=c('n_genes'), pt.size = 0, group.by = 'Integrated_05', adjust = 2) + NoLegend()
VlnPlot(gca_finalQC,features=c('total_counts'), pt.size = 0, group.by = 'Integrated_05', adjust = 2) + NoLegend()
VlnPlot(gca_finalQC,features=c('pct_counts_mt','n_genes','total_counts'), pt.size = 0, group.by = 'Integrated_05', stack = T, flip = T) + NoLegend()

############
#Downsample#
############
#Downsample Cell types with excessive representation compared with ILC target
gcA_samp = NULL
sample.size = 50
NotAugmented = data.frame()
for(cell in names(Cells)){
  sampled = F
  n_cells <- as.numeric(Cells[cell]) #number of cells in group
  if(n_cells == 0){
    next
  }
  intermediate_Seurat_Object <- gca_finalQC[,which(gca_finalQC[['Integrated_05']] == cell)]
  if(n_cells > sample.size){
    cat('sampling: ', cell, ' n = ', sample.size,  '\n')
    downSample_Seurat_obj <- intermediate_Seurat_Object[, sample(colnames(intermediate_Seurat_Object), size = sample.size, replace = F)]
    sampled = T
  }else{
    NAsize = dim(intermediate_Seurat_Object)[2]
    cat('Small sample size', NAsize, '\n')
    NotAugmented[cell,1]  <-  data.frame(NAsize)
  }
  if(is.null(gcA_samp)){
    message('gcA_samp init')
    if(sampled == T){
      gcA_samp <- downSample_Seurat_obj
    }else{
      gcA_samp <- intermediate_Seurat_Object
    }
  }else if(sampled == T){
    gcA_samp <- merge(gcA_samp, downSample_Seurat_obj)
  }else{
    gcA_samp <- merge(gcA_samp, intermediate_Seurat_Object)
  }
}
#Cell types which have not been downsampled due to low sample size
cat("Cell types not downsampled due to low samp size \n")
NotAugmented
#free up memory
rm(intermediate_Seurat_Object)
rm(downSample_Seurat_obj)
#rm(gca_finalQC)
gc()
#manually check proper merging of downsampled datasets
dim(gcA_samp)
table(gcA_samp@meta.data$Integrated_05)
table(gcA_samp@meta.data$Integrated_05)['ILC3']

#check data processing state
unique(gcA_samp@assays$originalexp@counts[,111])

gcA_samp@meta.data$Integrated_05 <- factor(gcA_samp@meta.data$Integrated_05, levels = f_order)
VlnPlot(gcA_samp,features=c('pct_counts_mt','n_genes','total_counts'),
        pt.size = 0, group.by = 'Integrated_05', stack = T, flip = T) + NoLegend()

#remove batches with cell counts too low to be able to integrate
table(gcA_samp$batch)
badBatches <- names(which(table(gcA_samp$batch) <10))
badCells <- which(gcA_samp$batch %in%badBatches)
badCellNames <- names(gcA_samp$batch[badCells])
retain <-  Cells(gcA_samp)[!Cells(gcA_samp) %in% badCellNames]

gca_adult_smallIntQC_downsampled <- subset(gcA_samp, cells = retain)
table(gca_adult_smallIntQC_downsampled@meta.data$Integrated_05)['ILC3']
table(gca_adult_smallIntQC_downsampled@meta.data$Integrated_05)
#write.csv(Cells(gca_adult_largeIntQC_downsampled), "data/gca_downsampled_cell_names.csv",
#          row.names = F, )

gca_adult_smallIntQC_downsampled <- RenameAssays(object = gca_adult_smallIntQC_downsampled, originalexp = "rna")
#gca_adult_smallIntQC_downsampled <- NormalizeData(gca_adult_smallIntQC_downsampled, normalization.method ="LogNormalize")
save(gca_adult_smallIntQC_downsampled, file='data/gca_adult_small_intestine/gca_adult_smallIntQC_downsampled.robj')

#prepare for cellphonedb

#DEGS####
#########
#optional for DEGs analysis
#find DEGs
# DEGS <- FindAllMarkers(gca_adult_largeIntQC_downsampled,
#                        test.use = 'LR',
#                        verbose = F,
#                        only.pos = T,
#                        random.seed =1,
#                        logfc.threshold = 0.2,
#                        min.pct = 0.1,
#                        return.thresh = 0.05)
# 
# cond <- DEGS[,'p_val_adj'] < 0.05 & DEGS[,"avg_log2FC"]> 0.1
# fDEGs = DEGS[which(cond),]
# 
# fDEGs <- fDEGs[,c('cluster','gene','p_val_adj', 'p_val', 'avg_log2FC', 'pct.1', 'pct.2')]
# write.table(fDEGs, file = 'gca_expanded_DEGs.tsv', sep='\t', quote = F, row.names = F)

# #write Metadata for statistical analysis
# Idents(gca_adult_largeIntQC_downsampled) <- gca_adult_largeIntQC_downsampled$Integrated_05
# meta <- cbind(colnames(gca_adult_largeIntQC_downsampled), as.character(Idents(gca_adult_largeIntQC_downsampled)))
# colnames(meta) <- c("Cell", "cell_type")
# write.table(x = meta, file = 'gca_Meta.tsv',  row.names = F, quote = F, sep = '\t')
# 
# 
# #write counts file to a format readab\le by cellphone db
# #NOTE: due to internal server storage problems h5ad files are not accessible
# library("Matrix")
# writeMM(gca_adult_largeIntQC_downsampled@assays$RNA@data, file ='adult_largeInt_cellphonedb_mtx/matrix.mtx')
# write(x=rownames(gca_adult_largeIntQC_downsampled@assays$RNA@data), file = 'adult_largeInt_cellphonedb_mtx/features.tsv')
# write(x = colnames(gca_adult_largeIntQC_downsampled@assays$RNA@data), file = 'adult_largeInt_cellphonedb_mtx/barcodes.tsv')
