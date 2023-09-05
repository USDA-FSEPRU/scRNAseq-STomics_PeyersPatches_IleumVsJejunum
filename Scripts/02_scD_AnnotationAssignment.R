### Load required software packages
#library(DropletUtils) 
library(Seurat)
library(SeuratObject)
library(ggplot2)
library(SeuratDisk)
library(dplyr)
library(readxl)
library(writexl)
library(readxl)
#library(scales)
#library(viridis)

## Load data
#Load in h5 Seurat object of our ileum/jejunum data:
seu <- LoadH5Seurat('/home/Jayne.Wiarda/SI_PP_SC_ST/NormalizeIntegrateDimReduc_SConly/AllSamples.h5Seurat')

## Load annotation information:
#To porcine ileum atlas:
CellTypePredictions1 <- read_excel('/home/Jayne.Wiarda/SI_PP_SC_ST/MappingPrediction/WiardaIleumAtlas_CellTypePredictions.xlsx') # read in predictions from mapping to pig ileum atlas data
MappingScores1 <- read_excel('/home/Jayne.Wiarda/SI_PP_SC_ST/MappingPrediction/WiardaIleumAtlas_MappingScores.xlsx')
colnames(CellTypePredictions1) <- paste("IleumAtlas", colnames(CellTypePredictions1), sep = "_")
colnames(MappingScores1) <- paste("IleumAtlas", colnames(MappingScores1), sep = "_")

#To porcine small intestinal epithelial cells:
CellTypePredictions2 <- read_excel('/home/Jayne.Wiarda/SI_PP_SC_ST/MappingPrediction/WiardaSIEpAtlas_CellTypePredictions.xlsx') # read in predictions from mapping to pig ileum atlas data
MappingScores2 <- read_excel('/home/Jayne.Wiarda/SI_PP_SC_ST/MappingPrediction/WiardaSIEpAtlas_MappingScores.xlsx')
colnames(CellTypePredictions2) <- paste("SIEpAtlas", colnames(CellTypePredictions2), sep = "_")
colnames(MappingScores2) <- paste("SIEpAtlas", colnames(MappingScores2), sep = "_")

dif <- setdiff(CellTypePredictions1$IleumAtlas_CellBarcodes, CellTypePredictions2$SIEpAtlas_CellBarcodes) # identify cell barcodes missing from the dataframe because they weren't predicted as epithelial cells
SIEpAtlas_CellBarcodes <- c(CellTypePredictions2$SIEpAtlas_CellBarcodes, dif)
SIEpAtlas_predicted.id <- c(CellTypePredictions2$SIEpAtlas_predicted.id, rep('NA', length(dif)))
SIEpAtlas_prediction.score.Enterocyte   <- c(CellTypePredictions2$SIEpAtlas_prediction.score.Enterocyte, rep('NA', length(dif)))
SIEpAtlas_prediction.score.NEUROD1lo.EE   <- c(CellTypePredictions2$SIEpAtlas_prediction.score.NEUROD1lo.EE, rep('NA', length(dif)))
SIEpAtlas_prediction.score.NEUROD1hi.EE   <- c(CellTypePredictions2$SIEpAtlas_prediction.score.NEUROD1hi.EE, rep('NA', length(dif)))
SIEpAtlas_prediction.score.Goblet   <- c(CellTypePredictions2$SIEpAtlas_prediction.score.Goblet, rep('NA', length(dif)))
SIEpAtlas_prediction.score.BEST4.enterocyte   <- c(CellTypePredictions2$SIEpAtlas_prediction.score.BEST4.enterocyte, rep('NA', length(dif)))
SIEpAtlas_prediction.score.Crypt   <- c(CellTypePredictions2$SIEpAtlas_prediction.score.Crypt, rep('NA', length(dif)))
SIEpAtlas_prediction.score.max   <- c(CellTypePredictions2$SIEpAtlas_prediction.score.max, rep('NA', length(dif)))
CellTypePredictions2 <- data.frame(SIEpAtlas_CellBarcodes, SIEpAtlas_predicted.id, SIEpAtlas_prediction.score.Enterocyte, SIEpAtlas_prediction.score.NEUROD1lo.EE, SIEpAtlas_prediction.score.NEUROD1hi.EE, SIEpAtlas_prediction.score.Goblet, SIEpAtlas_prediction.score.BEST4.enterocyte, SIEpAtlas_prediction.score.Crypt, SIEpAtlas_prediction.score.max)
CellTypePredictions2 <- CellTypePredictions2[ order(match(CellTypePredictions2$SIEpAtlas_CellBarcodes, colnames(seu))), ]

dif <- setdiff(MappingScores1$IleumAtlas_CellBarcodes, MappingScores2$SIEpAtlas_CellBarcodes) # identify cell barcodes missing from the dataframe because they weren't predicted as epithelial cells
SIEpAtlas_CellBarcodes <- c(MappingScores2$SIEpAtlas_CellBarcodes, dif)
SIEpAtlas_MappingScores <- c(MappingScores2$SIEpAtlas_MappingScores, rep('NA', length(dif)))
MappingScores2 <- data.frame(SIEpAtlas_CellBarcodes, SIEpAtlas_MappingScores)
MappingScores2 <- MappingScores2[ order(match(MappingScores2$SIEpAtlas_CellBarcodes, colnames(seu))), ]

## Add metadata:
seu <- AddMetaData(object = seu, 
                   metadata = c(CellTypePredictions1, MappingScores1, CellTypePredictions2, MappingScores2))

## Merge annotations:
#Identify cell barcodes corresponding to different epithelial cell type predictions:
ent <- rownames(seu@meta.data %>% filter(seu$SIEpAtlas_predicted.id == 'Enterocyte'))
cry <- rownames(seu@meta.data %>% filter(seu$SIEpAtlas_predicted.id == 'Crypt'))
nlo <- rownames(seu@meta.data %>% filter(seu$SIEpAtlas_predicted.id == 'NEUROD1lo EE'))
gob <- rownames(seu@meta.data %>% filter(seu$SIEpAtlas_predicted.id == 'Goblet'))
b4e <- rownames(seu@meta.data %>% filter(seu$SIEpAtlas_predicted.id == 'BEST4 enterocyte'))
nhi <- rownames(seu@meta.data %>% filter(seu$SIEpAtlas_predicted.id == 'NEUROD1hi EE'))

#Incorporate epithelial cell annotations into ileum atlas annotations:
bcs <- data.frame(colnames(seu), seu$IleumAtlas_predicted.id)
colnames(bcs) <- c('barcode', 'cellID')
bcs <- bcs %>% mutate(cellID = replace(cellID, barcode %in% ent, 'Enterocytes'))
bcs <- bcs %>% mutate(cellID = replace(cellID, barcode %in% cry, 'Crypt cells'))
bcs <- bcs %>% mutate(cellID = replace(cellID, barcode %in% nlo, 'NEUROD1lo EE cells'))
bcs <- bcs %>% mutate(cellID = replace(cellID, barcode %in% gob, 'Goblet cells'))
bcs <- bcs %>% mutate(cellID = replace(cellID, barcode %in% b4e, 'BEST4 enterocytes'))
bcs <- bcs %>% mutate(cellID = replace(cellID, barcode %in% nhi, 'NEUROD1hi EE cells'))

#Create new meta data slot:
seu$celltype <- bcs$cellID
DimPlot(seu, reduction = 'umap', group.by = 'celltype')

## Further define stromal cells
str <- LoadH5Seurat('/home/Jayne.Wiarda/SI_PP_SC_ST/AnnotatedSeurat/StromalCells.h5seurat')

## Add new metadata:
dif <- setdiff(colnames(seu), colnames(str)) # identify cell barcodes missing from the dataframe because they weren't predicted as stromal cells
Stromal_CellBarcodes <- c(colnames(str), dif)
x <- str$stromal
Stromal_Clusters <- c(as.character(str$stromal), rep('NA', length(dif)))
Stromal <- data.frame(Stromal_CellBarcodes, Stromal_Clusters)
Stromal <- Stromal[ order(match(Stromal$Stromal_CellBarcodes, colnames(seu))), ]
seu <- AddMetaData(object = seu, 
                   metadata = c(Stromal))

## Merge annotations:
#Identify cell barcodes corresponding to different epithelial cell type predictions:
fib <- rownames(seu@meta.data %>% filter(seu$Stromal_Clusters == 'Fibroblasts'))
end <- rownames(seu@meta.data %>% filter(seu$Stromal_Clusters == 'Endothelial cells'))
mus <- rownames(seu@meta.data %>% filter(seu$Stromal_Clusters == 'Muscle cells'))

#Incorporate stromal cell annotations into ileum atlas annotations:
bcs <- data.frame(colnames(seu), seu$celltype)
colnames(bcs) <- c('barcode', 'cellID')
bcs <- bcs %>% mutate(cellID = replace(cellID, barcode %in% fib, 'Fibroblasts'))
bcs <- bcs %>% mutate(cellID = replace(cellID, barcode %in% end, 'Endothelial cells'))
bcs <- bcs %>% mutate(cellID = replace(cellID, barcode %in% mus, 'Muscle cells'))

#Create new meta data slot:
seu$celltype <- bcs$cellID

#Set a logical order for new identities:
Idents(seu) <- seu$celltype
levels(seu) <- c('Activated B cells', 'Cycling B cells', 'Resting B cells', 'Transitioning B cells', 'Antibody-secreting cells',
                 'Cycling CD4 ab T cells', 'Cycling CD8 ab T cells', 'Cycling gd T cells', 'Cycling group 1 ILCs',
                 'Cytotoxic CD8 ab T cells', 'Cytotoxic gd T cells', 'Cytotoxic group 1 ILCs', 
                 'Non-naive CD8 ab T cells', 'Non-naive gd T cells', 'Non-naive group 1 ILCs', 
                 'SELLhi gd T cells', 'CD2neg GD T cells', 
                 'Naive CD4/CD8 ab T cells', 'Non-naive CD4 ab T cells', 'Follicular CD4 ab T cells', 'Group 3 ILCs',
                 'Dendritic cells', 'Macrophages', 'Mast cells',
                 'Crypt cells', 'Enterocytes', 'BEST4 enterocytes', 'Goblet cells', 'NEUROD1lo EE cells', 'NEUROD1hi EE cells',
                 'Endothelial cells', 'Fibroblasts', 'Muscle cells')
seu$celltype <- Idents(seu)

## Query canonical genes:
DefaultAssay(seu) <- 'SCT'
DotPlot(seu,
        features = unique(c('PTPRC',
                            'CD79A', 'CD79B', 'CD19', 'MS4A1', 'PAX5', # B cells
                            'JCHAIN', 'XBP1', 'PRDM1', 'IRF4', # ASCs
                            'CD69', 'CD83', 'SLA-DQB1', 'SLA-DRA', # B activation
                            'GPR183', 'CCR7', 'KLF2', 'SELL', 'FCER2', 'CD40', # resting B
                            'AICDA', 'CD86', 'BCL6', # follicular B
                            'PCLAF', 'BIRC5', 'TOP2A', 'STMN1', # cycling B
                            
                            'CD3E', 'CD3G', 'CD247', # T cell
                            'CD4', 'CD8B', 'TRDC', 'CD2', 'CD8A', # T cell subsets
                            'PCLAF', 'BIRC5', 'TOP2A', 'STMN1', # cycling
                            'CCL5', 'ITGAE', # effector/resident
                            'GZMA-16903', 'GZMB', 'GNLY', # cytotoxic
                            'CTSW', 'XCL1', 'SLA-DRA', 'SLA-DQB1', 'CCR9', 'KLRK1', # activation
                            'FCER1G', 'KLRG1', 'ITGB1', 'ITGB7', 'SELL', # SELLhi gd
                            'ID3', 'RHEX', 'BLK', 'SAMSN1', 'IL26', # CD2- gd
                            'CCR7', 'S1PR1', 'LEF1', 'KLF2', # naive ab
                            'ICOS', 'CTLA4', 'CD40LG', 'IL10', # Non-naive + follicular CD4
                            'CD52', 'IFITM3', 'GPR183', # non-naive CD4
                            'PDCD1', 'CXCR4', 'CD69', # follicular CD4
                            'LTB', 'ID2', 'KIT', 'IL7R', 'IL22', 'KLRB1', 'RORC', 'CXCL8', # group 3 ILC
                            
                            'FLT3', 'SLA-DRA', 'SLA-DQB1', # DC
                            'SIRPA', 'CD68', 'CXCL2', 'C1QA', 'C1QB', 'C1QC', # macrophage
                            'ICAM1', 'CSF2RB', 'MS4A2', 'FCER1A', # mast
                            
                            'RPL5', 'RPS6', 'EEF1B2', 'OLFM4', 'PIGR', 'LYZ', # crypts
                            'FABP2', 'FABP1', 'CLCA4', 'SLC5A1', 'SI', 'ACE2',  # enterocyte
                            'GUCA2A', 'GUCA2B', 'BEST4', 'CFTR', 'NOTCH2', 'OTOP2', # best4 enterocyte
                            'TFF3', 'REG4', 'CLCA1', 'SPINK4', 'MUC2', 'CXCL8', # goblet
                            'PYY', 'GAST', 'SST', 'CCK', 'TTR', 'NTS', # NEUROD1lo EE
                            'NEUROD1', 'CHGA', 'CHGB', 'KRT7', 'SCT', 'PENK', # NEUROD1hi EE
                            
                            'PECAM1', 'CDH5', # endothelial
                            'ECM1', 'COL1A1', 'COL1A2', # fibroblast
                            'TAGLN', 'MYH11', 'ACTG2')), # muscle
                          cols = c('gold', 'red3')) + RotatedAxis() 

## Reclassify annotations
#In dot plot above, most cell markers support annotations, except for cycling group 1 ILCs. We see strong expression of T cell (CD3E, CD3G, CD247) genes and more specifically CD8 ab T cell (CD8A, CD8B) genes. We will go ahead and re-calssify these as cycling CD8 ab T cells.
Idents(seu) <- seu$celltype
seu <- RenameIdents(seu, 'Cycling group 1 ILCs' = 'Cycling CD8 ab T cells')
seu$celltype <- Idents(seu)

Idents(seu) <- seu$celltype
levels(seu) <- c('Activated B cells', 'Cycling B cells', 'Resting B cells', 'Transitioning B cells', 'Antibody-secreting cells',
                 'Cycling CD4 ab T cells', 'Cycling CD8 ab T cells', 'Cycling gd T cells',
                 'Cytotoxic CD8 ab T cells', 'Cytotoxic gd T cells', 'Cytotoxic group 1 ILCs', 
                 'Non-naive CD8 ab T cells', 'Non-naive gd T cells', 'Non-naive group 1 ILCs', 
                 'SELLhi gd T cells', 'CD2neg GD T cells', 
                 'Naive CD4/CD8 ab T cells', 'Non-naive CD4 ab T cells', 'Follicular CD4 ab T cells', 'Group 3 ILCs',
                 'Dendritic cells', 'Macrophages', 'Mast cells',
                 'Crypt cells', 'Enterocytes', 'BEST4 enterocytes', 'Goblet cells', 'NEUROD1lo EE cells', 'NEUROD1hi EE cells',
                 'Endothelial cells', 'Fibroblasts', 'Muscle cells')
seu$celltype <- Idents(seu)

Idents(seu) <- seu$celltype

# Now re-plot:
DefaultAssay(seu) <- 'SCT'
DotPlot(seu,
        features = unique(c('PTPRC',
                            'CD79A', 'CD79B', 'CD19', 'MS4A1', 'PAX5', # B cells
                            'JCHAIN', 'XBP1', 'PRDM1', 'IRF4', # ASCs
                            'CD69', 'CD83', 'SLA-DQB1', 'SLA-DRA', # B activation
                            'GPR183', 'CCR7', 'KLF2', 'SELL', 'FCER2', 'CD40', # resting B
                            'AICDA', 'CD86', 'BCL6', # follicular B
                            'PCLAF', 'BIRC5', 'TOP2A', 'STMN1', # cycling B
                            
                            'CD3E', 'CD3G', 'CD247', # T cell
                            'CD4', 'CD8B', 'TRDC', 'CD2', 'CD8A', # T cell subsets
                            'PCLAF', 'BIRC5', 'TOP2A', 'STMN1', # cycling
                            'CCL5', 'ITGAE', # effector/resident
                            'GZMA-16903', 'GZMB', 'GNLY', # cytotoxic
                            'CTSW', 'XCL1', 'SLA-DRA', 'SLA-DQB1', 'CCR9', 'KLRK1', # activation
                            'FCER1G', 'KLRG1', 'ITGB1', 'ITGB7', 'SELL', # SELLhi gd
                            'ID3', 'RHEX', 'BLK', 'SAMSN1', 'IL26', # CD2- gd
                            'CCR7', 'S1PR1', 'LEF1', 'KLF2', # naive ab
                            'ICOS', 'CTLA4', 'CD40LG', 'IL10', # Non-naive + follicular CD4
                            'CD52', 'IFITM3', 'GPR183', # non-naive CD4
                            'PDCD1', 'CXCR4', 'CD69', # follicular CD4
                            'LTB', 'ID2', 'KIT', 'IL7R', 'IL22', 'KLRB1', 'RORC', 'CXCL8', # group 3 ILC
                            
                            'FLT3', 'SLA-DRA', 'SLA-DQB1', # DC
                            'SIRPA', 'CD68', 'CXCL2', 'C1QA', 'C1QB', 'C1QC', # macrophage
                            'ICAM1', 'CSF2RB', 'MS4A2', 'FCER1A', # mast
                            
                            'RPL5', 'RPS6', 'EEF1B2', 'OLFM4', 'PIGR', 'LYZ', # crypts
                            'FABP2', 'FABP1', 'CLCA4', 'SLC5A1', 'SI', 'ACE2',  # enterocyte
                            'GUCA2A', 'GUCA2B', 'BEST4', 'CFTR', 'NOTCH2', 'OTOP2', # best4 enterocyte
                            'TFF3', 'REG4', 'CLCA1', 'SPINK4', 'MUC2', 'CXCL8', # goblet
                            'PYY', 'GAST', 'SST', 'CCK', 'TTR', 'NTS', # NEUROD1lo EE
                            'NEUROD1', 'CHGA', 'CHGB', 'KRT7', 'SCT', 'PENK', # NEUROD1hi EE
                            
                            'PECAM1', 'CDH5', # endothelial
                            'ECM1', 'COL1A1', 'COL1A2', # fibroblast
                            'TAGLN', 'MYH11', 'ACTG2')), # muscle
        cols = c('gold', 'red3')) + RotatedAxis() 

# More refined gene list:
DefaultAssay(seu) <- 'SCT'
DotPlot(seu, features = c('PTPRC', 'CD79A', 'CD19', 'MS4A1', 
                          'JCHAIN', 'XBP1', 'PRDM1', 'IRF4', 
                          'CD69', 'CD83', 'CD40', 
                          'AICDA', 
                          'PCLAF', 'TOP2A', 
                          'CD3E', 'CD3G', 
                          'CD4', 'CD8B', 'TRDC', 'CD2', 'CD8A', 
                          'CCL5', 'ITGAE', 
                          'GZMB', 
                          'CCR9', 'KLRK1', 
                          'SELL', 
                          'BLK', 
                          'CCR7', 'S1PR1', 'ICOS', 
                          'CTLA4', 'CD40LG', 'PDCD1', 
                          'IL22', 'RORC',
                          'FLT3', 'SLA-DQB1', 
                          'SIRPA', 'CD68',
                          'MS4A2', 'FCER1A', 
                          'OLFM4', 'PIGR', 
                          'FABP2', 'SI', 
                          'GUCA2A', 'BEST4', 
                          'REG4', 'MUC2', 
                          'CHGA', 'CHGB',
                          'PECAM1', 'CDH5', 
                          'ECM1', 'COL1A1', 
                          'MYH11', 'ACTG2'), cols = c('gold', 'red3')) + RotatedAxis()

## Now make some other data plots on UMAP:
# Cell types
DimPlot(seu,
        cols = c('cyan4', 'gold3', 'chartreuse4', 'deeppink4', 'sandybrown',
                 'cornflowerblue', 'navy', 'lightpink', 'salmon', 
                 'deepskyblue2', 'tan4', 'mediumpurple1', 
                 'darkgreen', 'gray50', 'darkmagenta', 'red', 'hotpink', 'khaki', 
                 'orange4', 'limegreen', 'cadetblue3', 'firebrick', 'deepskyblue4', 'darkseagreen', 'burlywood3', 'black',
                 'goldenrod3', 'blue', 'deeppink', 'lightskyblue3', 'mistyrose3', 'purple'),
        shuffle = TRUE) #& NoLegend()

# Cell lineages
DimPlot(seu,
        cols = c(rep('mediumorchid', 5), rep('orange', 15), rep('blue', 3), rep('forestgreen', 6), rep('burlywood4', 3)),
        shuffle = TRUE) #& NoLegend()

# Sample IDs
DimPlot(seu,
        group.by = 'orig.ident',
        cols = c('red', 'red', 'red', 'blue', 'blue', 'blue'),
        shuffle = TRUE)

## Save Seurat objects with updated annotation:
SaveH5Seurat(seu, '/home/Jayne.Wiarda/SI_PP_SC_ST/AnnotatedSeurat/AllSamples_annotated.h5seurat', overwrite = TRUE)

### View session information
sessionInfo()
## R version 4.2.1 (2022-06-23)
## Platform: x86_64-pc-linux-gnu (64-bit)
## Running under: Ubuntu 20.04.4 LTS
## 
## Matrix products: default
## BLAS:   /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3
## LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/liblapack.so.3
## 
## locale:
##  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
##  [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
##  [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
##  [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                 
##  [9] LC_ADDRESS=C               LC_TELEPHONE=C            
## [11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
## 
## attached base packages:
## [1] stats     graphics  grDevices utils     datasets  methods   base     
## 
## other attached packages:
## [1] writexl_1.4.0         readxl_1.4.0          dplyr_1.0.9          
## [4] SeuratDisk_0.0.0.9020 ggplot2_3.3.6         sp_1.5-0             
## [7] SeuratObject_4.1.0    Seurat_4.1.1         
## 
## loaded via a namespace (and not attached):
##   [1] Rtsne_0.16            colorspace_2.0-3      deldir_1.0-6         
##   [4] ellipsis_0.3.2        ggridges_0.5.3        rstudioapi_0.13      
##   [7] spatstat.data_2.2-0   farver_2.1.1          leiden_0.4.2         
##  [10] listenv_0.8.0         bit64_4.0.5           ggrepel_0.9.1        
##  [13] RSpectra_0.16-1       fansi_1.0.3           codetools_0.2-18     
##  [16] splines_4.2.1         knitr_1.39            polyclip_1.10-0      
##  [19] jsonlite_1.8.0        ica_1.0-3             cluster_2.1.3        
##  [22] png_0.1-7             rgeos_0.5-9           uwot_0.1.11          
##  [25] shiny_1.7.2           sctransform_0.3.3     spatstat.sparse_2.1-1
##  [28] compiler_4.2.1        httr_1.4.3            assertthat_0.2.1     
##  [31] Matrix_1.4-1          fastmap_1.1.0         lazyeval_0.2.2       
##  [34] cli_3.3.0             later_1.3.0           htmltools_0.5.3      
##  [37] tools_4.2.1           igraph_1.3.4          gtable_0.3.0         
##  [40] glue_1.6.2            RANN_2.6.1            reshape2_1.4.4       
##  [43] Rcpp_1.0.9            scattermore_0.8       cellranger_1.1.0     
##  [46] vctrs_0.4.1           nlme_3.1-158          progressr_0.10.1     
##  [49] lmtest_0.9-40         spatstat.random_2.2-0 xfun_0.31            
##  [52] stringr_1.4.0         globals_0.15.1        mime_0.12            
##  [55] miniUI_0.1.1.1        lifecycle_1.0.1       irlba_2.3.5          
##  [58] goftest_1.2-3         future_1.27.0         MASS_7.3-58.1        
##  [61] zoo_1.8-10            scales_1.2.0          spatstat.core_2.4-4  
##  [64] promises_1.2.0.1      spatstat.utils_2.3-1  parallel_4.2.1       
##  [67] RColorBrewer_1.1-3    yaml_2.3.5            reticulate_1.25      
##  [70] pbapply_1.5-0         gridExtra_2.3         rpart_4.1.16         
##  [73] stringi_1.7.8         highr_0.9             rlang_1.0.4          
##  [76] pkgconfig_2.0.3       matrixStats_0.62.0    evaluate_0.15        
##  [79] lattice_0.20-45       ROCR_1.0-11           purrr_0.3.4          
##  [82] tensor_1.5            labeling_0.4.2        patchwork_1.1.1      
##  [85] htmlwidgets_1.5.4     bit_4.0.4             cowplot_1.1.1        
##  [88] tidyselect_1.1.2      parallelly_1.32.1     RcppAnnoy_0.0.19     
##  [91] plyr_1.8.7            magrittr_2.0.3        R6_2.5.1             
##  [94] generics_0.1.3        DBI_1.1.3             withr_2.5.0          
##  [97] mgcv_1.8-40           pillar_1.8.0          fitdistrplus_1.1-8   
## [100] survival_3.3-1        abind_1.4-5           tibble_3.1.8         
## [103] future.apply_1.9.0    hdf5r_1.3.5           crayon_1.5.1         
## [106] KernSmooth_2.23-20    utf8_1.2.2            spatstat.geom_2.4-0  
## [109] plotly_4.10.0         rmarkdown_2.14        grid_4.2.1           
## [112] data.table_1.14.2     digest_0.6.29         xtable_1.8-4         
## [115] tidyr_1.2.0           httpuv_1.6.5          munsell_0.5.0        
## [118] viridisLite_0.4.0