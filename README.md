# ILC3s modulate immune-epithelial interactions in the gut via TGF-β1 activation

### Authors

Diana Coman¹, John W. Bassett²†, Isabelle Coales¹†, Ainize Peña-Cearra¹†, Emily Read¹, Zuzanna Łukasik³, Robin J. Dart³, Mark A. Travis⁷⁸⁹, Jenny Mjösberg²¹⁰, Luke B. Roberts⁴, Joana F. Neves¹

†These authors contributed equally to this work.

<span style="font-size: 10px;"> 1. Centre for Host-Microbiome Interactions, King’s College London, United Kingdom  
2. Center for Infectious Medicine, Department of Medicine Huddinge, Karolinska Institutet, Karolinska University Hospital Huddinge, Stockholm, Sweden  
3. Centre for Inflammation Biology and Cancer Immunology, School of Immunology & Microbial Sciences, King’s College London at Guy’s Hospital Campus, London, UK  
4. Department of Rheumatology, Faculty of Medicine and Health Sciences, Ghent University and Unit for Molecular Immunology and Inflammation, VIB-UGent Center for Inflammation Research, Ghent, Belgium  
5. Immunosurveillance Laboratory, The Francis Crick Institute, London, UK  
6. Department of Gastroenterology, Guy’s and St Thomas’ Foundation Trust, London, UK  
7. Lydia Becker Institute of Immunology and Inflammation, University of Manchester, United Kingdom  
8. Wellcome Trust Centre for Cell-Matrix Research, University of Manchester, United Kingdom  
9. Division of Immunology, Immunity to Infection and Respiratory Medicine, Faculty of Biology, Medicine and Health, Manchester Academic Health Sciences Centre, University of Manchester, United Kingdom  
10. Clinical Lung- and Allergy Research Unit, Department of Medicine Huddinge, Karolinska Institutet, Stockholm, Sweden & the Department of Respiratory Medicine and Allergy, Karolinska University Hospital Huddinge, Stockholm, Sweden </span>
  

### This repository contains the code for generating Figure 4 for the manuscript.

![](ComanEtAl_Fig4.jpg)


### Data

#### Space-Time Gut Cell Atlas
- **Raw Data**: [Full_obj_log_counts_soupx_v2.h5ad](https://cellgeni.cog.sanger.ac.uk/gutcellatlas/Full_obj_log_counts_soupx_v2.h5ad)
- **Normalized Data**: [Full_obj_raw_counts_nosoupx_v2.h5ad](https://cellgeni.cog.sanger.ac.uk/gutcellatlas/Full_obj_raw_counts_nosoupx_v2.h5ad)
- **Source**: [Nature Article](https://www.nature.com/articles/s41586-021-03852-1)

#### Pediatric IBD
- **Cell Ranger v.2.1.1 Count Matrix**: [GSE169136](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE169136)
- **Annotation and Metadata**: [Kokkinou pIBD GitHub Repository](https://github.com/ramvinay/Kokkinou_pIBD)
- **Source**: [PMC Article](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC10213874/)

### File Overview

This project is organized into two main directories: `HPC` and `R`. Below is a brief description of the contents of each directory.

#### HPC Directory

The `HPC` directory contains Jupyter Notebook files related to processing Gut Cell Atlas data on a high-performance computing (HPC) cluster.

- **Downsample_integrate_checkCounts.ipynb**: A notebook for downsampling and checking counts in the dataset, primarily used for visualization and quality assurance.
- **ImportData_and_subset-Small_intestine.ipynb**: A notebook for importing and subsetting raw data specific to the small intestine.
- **ImportData_and_subset-Large_intestine.ipynb**: A notebook for importing and subsetting raw data specific to the large intestine.
- **ImportData_and_subset_intestine-Norm.ipynb**: A notebook for importing normalized intestine data.

#### R Directory

The `R` directory contains R scripts organized into two subdirectories: `gca` and `pIBD`.

#### gca (Fig 4 A and B)

The `gca` subdirectory contains scripts for analyzing Gut Cell Atlas data, specifically for large and small intestine datasets:

- **large_intestine (Fig 4 B)**:
  - **00_LI_H5ad_to_Seurat.r**: Converts H5ad format data to Seurat objects for large intestine analysis.
  - **01_LI_QC_by_10xrun-cellType.r**: Quality control script for large intestine data, focusing on 10x run and cell type.
  - **02_LI_Downsample_cell_types.r**: Downsampling script for cell types in large intestine data.
  - **03_LI_Seurat_umap.r**: Generates UMAP visualizations for large intestine data using Seurat.
  - **AverageLogExpFig4.R**: Generates figures based on average log expression for large intestine data.

- **small_intestine (Fig 4 A)**:
  - **00_Small_int_H5ad_to_Seurat.r**: Converts H5ad format data to Seurat objects for small intestine analysis.
  - **01_SI_QC_by_10xrun-cellType.r**: Quality control script for small intestine data, focusing on 10x run and cell type.
  - **02_SI_Downsample_cell_types.r**: Downsampling script for cell types in small intestine data.
  - **03_SI_Seurat_umap.r**: Generates UMAP visualizations for small intestine data using Seurat.
  - **SI_genFig4a.R**: Generates figures for small intestine data.

#### pIBD (Fig 4 C)

The `pIBD` subdirectory contains scripts related to the integration and analysis of data for the pediatric inflammatory bowel disease dataset (pIBD):

- **00_integrate_cca.R**: Script for integrating datasets using canonical correlation analysis (CCA).
- **00integrate_comBat.R**: Script for batch effect correction using ComBat. (used for figure)
- **genFigure4e.R**: Generates specific figures related to the IBD analysis.

### Acknowledgments

Operations performed on publicly available data were enabled by resources provided by the National Academic Infrastructure for Supercomputing in Sweden (NAISS), partially funded by the Swedish Research Council through grant agreement no. 2022-06725.
