library('Seurat')
#library('ggplot2')
load('data/gca_adult_large_intestine/gca_adult_largeInt_Seurat.robj')


##split into individual 10x runs
gca_list <- SplitObject(gcA_adult_largeInt, split.by='sample_name')

#Count TA cells
sumInds = 0
for(run in 1:length(gca_list)){
  inds <- which(gca_list[[run]]$Integrated_05 == 'TA')
  sumInds <- sumInds+ length(inds)
}
print(sumInds) 

#Iterate through experiments and threshold according to Mt strategy, remove TA cells, remove poor experiments
for(run in 1:length(gca_list)){
  ExperName <- names(gca_list[run])
  cat("Experiment: ", ExperName, "\n")
  #remove all TA cells
  #inds <- which(!gca_list[[run]]$Integrated_05 == 'TA') #index of cells which not TA
  #cat(" TA Cells: ", length(inds), "\n")
  #tempSeurat <- subset(gca_list[[run]], cells = inds)
  
  #Keep TA cells
  tempSeurat <- gca_list[[run]]
  
  #check how many cells are in exp and flag low cell exps
  if(dim(tempSeurat)[2] < 100){
    message('warning ', names(gca_list[run]),
            ' has low cell number: ', dim(tempSeurat)[2] )
    low.cells <- table(tempSeurat$Integrated_05)
    low.cells <- low.cells[low.cells>0]
    print(low.cells)
  }
  
  #Variable Threshold of %MT
  retain = NULL
  cats = names(table(tempSeurat$Integrated_05))
  for(type in cats){
    cell.inds = which(tempSeurat$Integrated_05==type)
    if(length(cell.inds)<2){
      next
    }
    cell.mt <- tempSeurat$pct_counts_mt[cell.inds]
    cell.mt.q <- quantile(cell.mt, probs =.75)
    
    if(type %in% c('LEC3 (ADGRG3+)','LEC5 (CLDN11+)', 'Th17', 'Mesothelium(PRG4+)', 'cDC1')){
      retain <- c(retain, names(which(cell.mt < 30)))
    }else if(cell.mt.q < 14 | (type %in% c('IgA plasma cell',
                                     'Activated CD4 T',
                                     'Adult Glia',
                                     'Stromal 1 (ADAMDEC1+)',
                                     'TA'))){ #set threshold to 15
      retain <- c(retain, names(which(cell.mt < 15)))
    }else if(cell.mt.q < 40){ # set threshold to 30
      retain <- c(retain, names(which(cell.mt < 20)))
    }else if(type %in% c('Macrophages', 'Stem cells', 'gdT', 'cDC2')){
      retain <- c(retain, names(which(cell.mt < 30)))
    }else{ #set threshold to 50
      retain <- c(retain, names(which(cell.mt < 50)))
    }
    
  }
  
  tempSeurat = subset(tempSeurat, cells = retain)
  gca_list[[names(gca_list[run])]] <- tempSeurat
}

#find all cells to retain
retain.cell <- NULL
for(i in 1:length(gca_list)){
  retain.cell = c(retain.cell, Cells(gca_list[[i]]))
}
#Remove cells which failed qc
gca_final <- subset(gcA_adult_largeInt, cells = retain.cell)
names.remove = names(which(table(gca_final$Integrated_05) < 3))
#Remove cells which have an n which is less than 3
gca_final <- subset(gca_final, cells = which(!gca_final$Integrated_05 %in% names.remove))
diff <- cbind(table(gca_final$Integrated_05), table(gcA_adult_largeInt$Integrated_05))
gca_final$Integrated_05 <- factor(gca_final$Integrated_05)

#check celltype representation
colnames(diff) <- c("final", "original")
diff

#Visual finalized data after qc
VlnPlot(gca_final, features = 'pct_counts_mt', group.by = 'Integrated_05', pt.size = 0) + NoLegend()
p2 = FeatureScatter(gca_final,
                    feature2 = 'pct_counts_mt',
                    feature1 = 'total_counts',
                    group.by = 'Integrated_05' ) + NoLegend()
p3 = FeatureScatter(gca_final,
                    feature2 = 'n_genes_by_counts',
                    feature1 = 'total_counts',
                    group.by = 'Integrated_05' ) + NoLegend()

p2|p3

cells.in.orig <- names(table(gcA_adult_largeInt$Integrated_05))
cells.in.final <- names(table(gca_final$Integrated_05))

#check which cells types were removed
cells.types.removed <- cells.in.orig[!(cells.in.orig %in% cells.in.final)]
cells.types.removed
#check cell types of interest that were removed for value in analysis
table(gcA_adult_largeInt$Integrated_05)['Stromal 2 (CH25H+)']
table(gcA_adult_largeInt$Integrated_05)['LEC5 (CLDN11+)']
table(gca_final$Integrated_05)['LEC5 (CLDN11+)']
table(gcA_adult_largeInt$Integrated_05)['LEC3 (ADGRG3+)']
table(gca_final$Integrated_05)['LEC3 (ADGRG3+)']
#no conflict with removing stromal 2 cell group. There were only 2 in the 
#original analysis

#save final data as robj
save(gca_final, file = 'data/gca_adult_large_intestine/gca_adult_LargInt_Seurat_QC.robj')
