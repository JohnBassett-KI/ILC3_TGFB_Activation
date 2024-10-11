library('sva')
library('Seurat')

pIBD <- readRDS(file =  "data/pIBD/pIBD.RDS")


pIBD <- FindVariableFeatures(pIBD)
#TGFBsig <- c("ITGAV", "ITGB8", "ITGB3", "MMP9", "MMP14", "SMAD2", "SMAD3", "SMAD7", "TGFB1", "TGFBR1", "TGFBR2", "SDC4", "TGM2", "FN1", "FBN1", "LTBP1", "LTBP3")  
load("data/TGFBsig.robj")
var_genes <- c(VariableFeatures(pIBD), TGFBsig)

expr_matrix <- GetAssayData(pIBD, layer = "data")
expr_matSubset <- expr_matrix[rownames(expr_matrix) %in% var_genes, ]
batch_info <- factor(pIBD$experiment)

# Apply ComBat for batch correction
int_pIBD <- ComBat(dat = as.matrix(expr_matSubset), batch = batch_info, par.prior = TRUE, prior.plots = FALSE)


# Store the updated expression matrix back into the Seurat object
batchCorr <- CreateAssayObject(counts = as.matrix(int_pIBD))
pIBD[["comBat"]] <- batchCorr

saveRDS(pIBD, file="data/pIBD/combatCorr_pIBD.rds")
