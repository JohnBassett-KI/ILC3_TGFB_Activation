library('Seurat')


load('data/gca_adult_small_intestine/gca_SI_Seurat.robj')
##split into individual 10x runs
gca_list <- SplitObject(gca_adult_smallInt, split.by='sample_name')

#Calculate the number of TA cells in the dataset:
#Transit-Amplifying cells are an undifferentiated population in transition 
#between SCs and differentiated cells. In the Large Intestine they are 
#disproportionately enriched in the data. Here we check to see if this is the 
#same for the small intestine. 
sumInds = 0
for(run in 1:length(gca_list)){
  inds <- which(gca_list[[run]]$Integrated_05 == 'TA')
  sumInds <- sumInds+ length(inds)
}
print(sumInds) 
#2698 TA cells found: no issue.

#Manual look at percent.mt
#VlnPlot(gca_adult_smallInt, features = "pct_counts_mt", group.by = "Integrated_05", pt.size = 0)+
#geom_boxplot(width = 0.1, fill = "white", color = "black", outlier.shape = NA) +
#NoLegend()

#load outlier thresholds (This will be regenerated in the following loop)
load(file = "data/gca_adult_small_intestine/medianOutlierThresholds.robj")
#Iterate through experiments and threshold according to Mt strategy, remove TA cells, remove poor experiments
MitoByCell_Outliers = list()
for(run in 1:length(gca_list)){
  experName <- names(gca_list[run])
  cat('QC for Experiment ',experName,'\n')
  #remove all TA cells
  ##inds <- which(!gca_list[[run]]$Integrated_05 == 'TA') #index of cells which not TA
  ##tempSeurat <- subset(gca_list[[run]], cells = inds)

  #Keep all TA cells
  tempSeurat <- gca_list[[run]]
  
  #check how many cells are in experiment and flag low cell experiments
  if(dim(tempSeurat)[2] < 500){
    message('warning ', experName,
            ' has low cell number: ', dim(tempSeurat)[2] )
    low.cells <- table(tempSeurat$Integrated_05)
    low.cells <- low.cells[low.cells>0]
    print(low.cells)
  }else{
    cat('Containing ', dim(tempSeurat)[2],' cells \n')
  }
  
  #Variable Threshold of %MT
  retain = NULL
  ExperOutlierCol <- list()
  cats = names(table(tempSeurat$Integrated_05))
  for(type in cats){
    cell.inds = which(tempSeurat$Integrated_05==type)
    if(length(cell.inds)<2){
      next
    }
    cell.mt <- tempSeurat$pct_counts_mt[cell.inds]
    cell.mt.q <- quantile(cell.mt, probs =c(0.25,.75))
    IQR = cell.mt.q[[2]]-cell.mt.q[[1]]
    OutThresh = 1.5*IQR + cell.mt.q[[2]]
    ExperOutlierCol[type] <- OutThresh
    ############################################################################################
    #MANUAL CHECK
    #######################
    # Generate a VlnPlotx
    #VlnPlot(tempSeurat, features = "pct_counts_mt", group.by = "Integrated_05", pt.size = 0)+ 
    #  geom_boxplot(width = 0.1, fill = "white", color = "black", outlier.shape = NA) + NoLegend()
    ############################################################################################
    cat("-- ", type,
        "\n      Experiment Mitcondrial Threshold: ",
        OutThresh,
        "\n median Outlier threshold for celltype: ", 
        medians[[type]], "\n")
    
    #Set the outlier threshold: 20 unless the median outlier threshold is greater. Then variable thresholding.
    if(medians[[type]]> 50){
      FinalThresh = 50
    }else if(medians[[type]]< 20){
      FinalThresh = 20
    }else{
      FinalThresh = medians[[type]]
    }
    retain <- c(retain, names(which(cell.mt < FinalThresh)))
  }
  MitoByCell_Outliers[experName] <- list(ExperOutlierCol)
  tempSeurat = subset(tempSeurat, cells = retain)
  gca_list[[experName]] <- tempSeurat
}

#Generate a dataframe of all the IQR based outlier thresholds for each experiment and cell type
cellTypes <- unique(gca_adult_smallInt$Integrated_05)
OutlierFrame <- data.frame()
for(name in names(MitoByCell_Outliers)){
  for(cell in cellTypes){
    if(is.null(MitoByCell_Outliers[[name]][[cell]])){
      element = NA
    }else{
      element = MitoByCell_Outliers[[name]][[cell]]
    }
    OutlierFrame[name, cell] <- element
  }
}
medians <- apply(OutlierFrame, 2, mean, na.rm=T)
#save(medians, file = "data/gca_adult_small_intestine/medianOutlierThresholds.robj")

#find all cells to retain
retain.cell <- NULL
for(i in 1:length(gca_list)){
  retain.cell = c(retain.cell, Cells(gca_list[[i]]))
}
#Remove cells which failed qc
gca_final <- subset(gca_adult_smallInt, cells = retain.cell)
names.remove = names(which(table(gca_final$Integrated_05) < 3))
#Remove cells which have an n which is less than 3
gca_final <- subset(gca_final, cells = which(!gca_final$Integrated_05 %in% names.remove))
diff <- cbind(table(gca_final$Integrated_05), table(gca_adult_smallInt$Integrated_05))
gca_final$Integrated_05 <- factor(gca_final$Integrated_05)

#check celltype representation
colnames(diff) <- c("final", "original")
print(diff)

#Visual finalized data after qc
VlnPlot(gca_adult_smallInt, features = 'pct_counts_mt', group.by = 'Integrated_05', pt.size = 0) + NoLegend()
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

cells.in.orig <- names(table(gca_adult_smallInt$Integrated_05))
cells.in.final <- names(table(gca_final$Integrated_05))

#check which cells types were removed
cells.types.removed <- cells.in.orig[!(cells.in.orig %in% cells.in.final)]
cells.types.removed

#Remove Cell types with below 10 cells in the dataset
keep <- (diff > 10)[,1]
keep <- names(which(keep))
Idents(gca_final) <- gca_final$Integrated_05
gca_final <- subset(gca_final, idents = keep, invert = F)
#save final data as robj
save(gca_final, file = 'data/gca_adult_small_intestine/gca_adult_smallInt_Seurat_QC.robj')
