library(DropletUtils) 
library(Seurat)
library(SeuratObject)
library(ggplot2)
library(SeuratDisk)
library(dplyr)
library(readxl)
library(writexl)
library(scales)

## Create Seurat object list

#Read in Seurat object with all six samples included:
seu <- readRDS('/home/Jayne.Wiarda/SI_PP_SC_ST/QC_STonly/AllSamples_PostQC.rds')

#See how many spots/genes in Seurat object:
seu

#See how many spots in each sample:
table(seu$SampleID)

#Split Seurat object into list by sample ID:
seu.list <- SplitObject(seu, split.by = "SampleID") # split by sample IDs 
rm(seu)

#Reorganize images appropriately in Seurat list:
seu.list$I2@images$I3 <- NULL
seu.list$I2@images$I4 <- NULL
seu.list$I2@images$J2 <- NULL
seu.list$I2@images$J3 <- NULL
seu.list$I2@images$J4 <- NULL

seu.list$I3@images$I2 <- NULL
seu.list$I3@images$I4 <- NULL
seu.list$I3@images$J2 <- NULL
seu.list$I3@images$J3 <- NULL
seu.list$I3@images$J4 <- NULL

seu.list$I4@images$I2 <- NULL
seu.list$I4@images$I3 <- NULL
seu.list$I4@images$J2 <- NULL
seu.list$I4@images$J3 <- NULL
seu.list$I4@images$J4 <- NULL

seu.list$J2@images$J3 <- NULL
seu.list$J2@images$J4 <- NULL
seu.list$J2@images$I2 <- NULL
seu.list$J2@images$I3 <- NULL
seu.list$J2@images$I4 <- NULL

seu.list$J3@images$J2 <- NULL
seu.list$J3@images$J4 <- NULL
seu.list$J3@images$I2 <- NULL
seu.list$J3@images$I3 <- NULL
seu.list$J3@images$I4 <- NULL

seu.list$J4@images$J2 <- NULL
seu.list$J4@images$J3 <- NULL
seu.list$J4@images$I2 <- NULL
seu.list$J4@images$I3 <- NULL
seu.list$J4@images$I4 <- NULL

#See how many spots/genes in each listed Seurat object:
seu.list

## Normalize data using SCT method
#Performed SCTransform on each sample from list:
for (i in 1:length(seu.list)) { # normalize data using SCTransform method
  seu.list[[i]] <- SCTransform(seu.list[[i]], 
                               return.only.var.genes = FALSE, 
                               verbose = TRUE,
                               assay = 'Spatial') 
}

## Integration and dimensionality reductions on different data combinations
#We don't yet know if we want to integrate together all data versus analyzing jejunum and ileum separately, so we will end up doing integration and dimensionality reduction on three sets of data:
  
#* all samples together (3 jejunum + 3 ileum = 6 samples)
#* only 3 jejunum samples
#* only 3 ileum samples

### Integration & dimensionality reduction for all samples

#### Integrate the data from different samples:
#Since SCTransform took awhile to run, let's store the SCT-corrected seu.list in another object we can tinker with for this combination of data:
seu.list2 <- seu.list

#Find variable features and integrate data:
seu.features <- SelectIntegrationFeatures(seu.list2, # select the genes to use for integration
                                          verbose = TRUE) 
seu.list2 <- PrepSCTIntegration(seu.list2, 
                               anchor.features = seu.features,
                               verbose = TRUE)
seu.anchors <- FindIntegrationAnchors(seu.list2, # identify anchors for integration from top 30 data dimensions
                                      normalization.method = "SCT", 
                                      anchor.features = seu.features, 
                                      dims = 1:30)
seu.integrated <- IntegrateData(seu.anchors, # integrate data
                                normalization.method = "SCT", 
                                dims = 1:30)
rm(seu.list2, seu.features, seu.anchors)

#### Calculate principle components
#Calculate 50 PCs:
seu.integrated <- RunPCA(seu.integrated, # run PCA analysis for 50 dimensions of the data
                         npcs = 50, 
                         verbose = TRUE) 

#Visualize PCs:
#ElbowPlot(seu.integrated,
#          ndims = 50) # look at this plot to find the 'elbow' for significant PCs... use this number of PCs for creating UMAP, tSNE, & spot neighbors & clustering

#Quantitiatively calculate PC cutoff, which we will use to set our data dimensions in most of our subsequent analyses:
pct <- seu.integrated[["pca"]]@stdev / sum(seu.integrated[["pca"]]@stdev) * 100 # find standard deviation for each PC
cumu <- cumsum(pct) # find cumulative percentages for PCs
co1 <- which(cumu > 90 & pct < 5)[1] # find PC representing cumulative percent >90% and less than 5% associated with the single PC
co1 # list PC
co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1 # find last PC where change in percent variation is more than 0.1%
co2 # list PC
pcs <- min(co1, co2) # find the minimum PC from the 2 methods used above
pcs # list PC
plot_df <- data.frame(pct = pct, # put PC values into dataframe for plotting
                      cumu = cumu, 
                      rank = 1:length(pct))
ggplot(plot_df, aes(cumu, pct, label = rank, color = rank > pcs)) + # visualize PCs to use in elbow plot
  geom_text() + 
  geom_vline(xintercept = 90, color = "grey") + 
  geom_hline(yintercept = min(pct[pct > 5]), color = "grey") +
  theme_bw()
PCdims <- 1:pcs # use the minimum PC from the quantitative method above to set the PCs for subsequent steps
rm(pct, cumu, co1, co2, pcs)

#### Perform multidimensional visualization of data:
#Create UMAP dimensions:
seu.integrated <- RunUMAP(seu.integrated, 
                          dims = PCdims, # use our calculated number of PCs
                          reduction = "pca", 
                          dim_embed = 3, # calculate 3 plot dimensions in case we want to try 3D plotting later
                          min.dist = 0.5,
                          spread = 0.2,
                          assay = "SCT") # create UMAP

#Create t-SNE dimensions:
seu.integrated <- RunTSNE(seu.integrated, 
                          dims = PCdims, # use our calculated number of PCs
                          reduction = "pca", 
                          dim_embed = 3) # calculate 3 plot dimensions in case we want to try 3D plotting later


#Visualize UMAP:
DimPlot(seu.integrated,
        reduction = 'umap',
        group.by = 'SampleID')

#Visualize t-SNE:
DimPlot(seu.integrated,
        reduction = 'tsne',
        group.by = 'SampleID')

#### Add normalized/scaled data to RNA assay
#dim(seu.integrated[["Spatial"]]@scale.data) # see that there is no RNA assay scaled data yet
seu.integrated <- NormalizeData(seu.integrated,  # normalize the RNA counts data per spot
                                normalization.method = "LogNormalize", 
                                scale.factor = 10000, 
                                assay = "Spatial")
seu.integrated <- ScaleData(seu.integrated, # scale the RNA counts data relative to other spots
                            assay = "Spatial")
seu.integrated <- ScaleData(seu.integrated, # scale the SCT counts data relative to other spots
                            assay = "SCT")
#dim(seu.integrated[["Spatial"]]@scale.data) # see that all genes are scaled in RNA assay now
#dim(seu.integrated[["SCT"]]@scale.data) # see that all genes are scaled in SCT assay now

#### Save the Seurat object
saveRDS(seu.integrated, '/home/Jayne.Wiarda/SI_PP_SC_ST/NormalizeIntegrateDimReduc_STonly/AllSamples.rds')
rm(seu.integrated)

### Integration & dimensionality reduction for only jejunal samples

#### Integrate the data from different samples:
#Since SCTransform took awhile to run, let's store the SCT-corrected seu.list in another object we can tinker with for this combination of data:
seu.list2 <- seu.list

#Remove unwanted samples from the list:
seu.list2 # see where our jejunal vs. ileal samples are in the list
seu.list2 <- seu.list2[-c(1:3)] # remove unwanted samples
seu.list2 # make sure we have the correct samples left

#Find variable features and integrate data:
seu.features <- SelectIntegrationFeatures(seu.list2, # select the genes to use for integration
                                          verbose = TRUE) 
seu.list2 <- PrepSCTIntegration(seu.list2, 
                                anchor.features = seu.features,
                                verbose = TRUE)
seu.anchors <- FindIntegrationAnchors(seu.list2, # identify anchors for integration from top 30 data dimensions
                                      normalization.method = "SCT", 
                                      anchor.features = seu.features, 
                                      dims = 1:30)
seu.integrated <- IntegrateData(seu.anchors, # integrate data
                                normalization.method = "SCT", 
                                dims = 1:30)
rm(seu.list2, seu.features, seu.anchors)

#### Calculate principle components
# Calculate 50 PCs:
seu.integrated <- RunPCA(seu.integrated, # run PCA analysis for 50 dimensions of the data
                         npcs = 50, 
                         verbose = TRUE) 

#Visualize PCs:
#ElbowPlot(seu.integrated,
#          ndims = 50) # look at this plot to find the 'elbow' for significant PCs... use this number of PCs for creating UMAP, tSNE, & spot neighbors & clustering

#Quantitiatively calculate PC cutoff, which we will use to set our data dimensions in most of our subsequent analyses:
pct <- seu.integrated[["pca"]]@stdev / sum(seu.integrated[["pca"]]@stdev) * 100 # find standard deviation for each PC
cumu <- cumsum(pct) # find cumulative percentages for PCs
co1 <- which(cumu > 90 & pct < 5)[1] # find PC representing cumulative percent >90% and less than 5% associated with the single PC
co1 # list PC
co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1 # find last PC where change in percent variation is more than 0.1%
co2 # list PC
pcs <- min(co1, co2) # find the minimum PC from the 2 methods used above
pcs # list PC
plot_df <- data.frame(pct = pct, # put PC values into dataframe for plotting
                      cumu = cumu, 
                      rank = 1:length(pct))
ggplot(plot_df, aes(cumu, pct, label = rank, color = rank > pcs)) + # visualize PCs to use in elbow plot
  geom_text() + 
  geom_vline(xintercept = 90, color = "grey") + 
  geom_hline(yintercept = min(pct[pct > 5]), color = "grey") +
  theme_bw()
PCdims <- 1:pcs # use the minimum PC from the quantitative method above to set the PCs for subsequent steps
rm(pct, cumu, co1, co2, pcs)

#### Perform multidimensional visualization of data:
#Create UMAP dimensions:
seu.integrated <- RunUMAP(seu.integrated, 
                          dims = PCdims, # use our calculated number of PCs
                          reduction = "pca", 
                          dim_embed = 3, # calculate 3 plot dimensions in case we want to try 3D plotting later
                          min.dist = 0.5,
                          spread = 0.2,
                          assay = "SCT") # create UMAP

#Create t-SNE dimensions:
seu.integrated <- RunTSNE(seu.integrated, 
                          dims = PCdims, # use our calculated number of PCs
                          reduction = "pca", 
                          dim_embed = 3, # calculate 3 plot dimensions in case we want to try 3D plotting later
                          assay = "SCT") # create tSNE plot 

#Visualize UMAP:
DimPlot(seu.integrated,
        reduction = 'umap',
        group.by = 'SampleID')

#Visualize t-SNE:
DimPlot(seu.integrated,
        reduction = 'tsne',
        group.by = 'SampleID')

#### Add normalized/scaled data to RNA assay
#dim(seu.integrated[["Spatial"]]@scale.data) # see that there is no RNA assay scaled data yet
seu.integrated <- NormalizeData(seu.integrated,  # normalize the RNA counts data per spot
                                normalization.method = "LogNormalize", 
                                scale.factor = 10000, 
                                assay = "Spatial")
seu.integrated <- ScaleData(seu.integrated, # scale the RNA counts data relative to other spots
                            assay = "Spatial")
seu.integrated <- ScaleData(seu.integrated, # scale the SCT counts data relative to other spots
                            assay = "SCT")
#dim(seu.integrated[["Spatial"]]@scale.data) # see that all genes are scaled in RNA assay now
#dim(seu.integrated[["SCT"]]@scale.data) # see that all genes are scaled in SCT assay now


#### Save the Seurat object
saveRDS(seu.integrated, '/home/Jayne.Wiarda/SI_PP_SC_ST/NormalizeIntegrateDimReduc_STonly/JejunumOnly.rds')
rm(seu.integrated)

### Integration & dimensionality reduction for only jejunal samples

#### Integrate the data from different samples:

#Since SCTransform took awhile to run, let's store the SCT-corrected seu.list in another object we can tinker with for this combination of data:
seu.list2 <- seu.list

#Remove unwanted samples from the list:
seu.list2 # see where our jejunal vs. ileal samples are in the list
seu.list2 <- seu.list2[-c(4:6)] # remove unwanted samples
seu.list2 # make sure we have the correct samples left

#Find variable features and integrate data:
seu.features <- SelectIntegrationFeatures(seu.list2, # select the genes to use for integration
                                          verbose = TRUE) 
seu.list2 <- PrepSCTIntegration(seu.list2, 
                               anchor.features = seu.features,
                               verbose = TRUE)
seu.anchors <- FindIntegrationAnchors(seu.list2, # identify anchors for integration from top 30 data dimensions
                                      normalization.method = "SCT", 
                                      anchor.features = seu.features, 
                                      dims = 1:30)
seu.integrated <- IntegrateData(seu.anchors, # integrate data
                                normalization.method = "SCT", 
                                dims = 1:30)
rm(seu.list2, seu.features, seu.anchors)

#### Calculate principle components
#Calculate 50 PCs:
seu.integrated <- RunPCA(seu.integrated, # run PCA analysis for 50 dimensions of the data
                         npcs = 50, 
                         verbose = TRUE) 

#Visualize PCs:
#ElbowPlot(seu.integrated,
#          ndims = 50) # look at this plot to find the 'elbow' for significant PCs... use this number of PCs for creating UMAP, tSNE, & spot neighbors & clustering

#Quantitiatively calculate PC cutoff, which we will use to set our data dimensions in most of our subsequent analyses:
pct <- seu.integrated[["pca"]]@stdev / sum(seu.integrated[["pca"]]@stdev) * 100 # find standard deviation for each PC
cumu <- cumsum(pct) # find cumulative percentages for PCs
co1 <- which(cumu > 90 & pct < 5)[1] # find PC representing cumulative percent >90% and less than 5% associated with the single PC
co1 # list PC
co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1 # find last PC where change in percent variation is more than 0.1%
co2 # list PC
pcs <- min(co1, co2) # find the minimum PC from the 2 methods used above
pcs # list PC
plot_df <- data.frame(pct = pct, # put PC values into dataframe for plotting
                      cumu = cumu, 
                      rank = 1:length(pct))
ggplot(plot_df, aes(cumu, pct, label = rank, color = rank > pcs)) + # visualize PCs to use in elbow plot
  geom_text() + 
  geom_vline(xintercept = 90, color = "grey") + 
  geom_hline(yintercept = min(pct[pct > 5]), color = "grey") +
  theme_bw()
PCdims <- 1:pcs # use the minimum PC from the quantitative method above to set the PCs for subsequent steps
rm(pct, cumu, co1, co2, pcs)

#### Perform multidimensional visualization of data:
#Create UMAP dimensions:
seu.integrated <- RunUMAP(seu.integrated, 
                          dims = PCdims, # use our calculated number of PCs
                          reduction = "pca", 
                          dim_embed = 3, # calculate 3 plot dimensions in case we want to try 3D plotting later
                          min.dist = 0.5,
                          spread = 0.2,
                          assay = "SCT") # create UMAP

#Create t-SNE dimensions:
seu.integrated <- RunTSNE(seu.integrated, 
                          dims = PCdims, # use our calculated number of PCs
                          reduction = "pca", 
                          dim_embed = 3, # calculate 3 plot dimensions in case we want to try 3D plotting later
                          assay = "SCT") # create tSNE plot 

#Visualize UMAP:
DimPlot(seu.integrated,
        reduction = 'umap',
        group.by = 'SampleID')

#Visualize t-SNE:
DimPlot(seu.integrated,
        reduction = 'tsne',
        group.by = 'SampleID')

#### Add normalized/scaled data to RNA assay
#dim(seu.integrated[["Spatial"]]@scale.data) # see that there is no RNA assay scaled data yet
seu.integrated <- NormalizeData(seu.integrated,  # normalize the RNA counts data per spot
                                normalization.method = "LogNormalize", 
                                scale.factor = 10000, 
                                assay = "Spatial")
seu.integrated <- ScaleData(seu.integrated, # scale the RNA counts data relative to other spots
                            assay = "Spatial")
seu.integrated <- ScaleData(seu.integrated, # scale the SCT counts data relative to other spots
                            assay = "SCT")
#dim(seu.integrated[["Spatial"]]@scale.data) # see that all genes are scaled in RNA assay now
#dim(seu.integrated[["SCT"]]@scale.data) # see that all genes are scaled in SCT assay now

#### Save the Seurat object
saveRDS(seu.integrated, '/home/Jayne.Wiarda/SI_PP_SC_ST/NormalizeIntegrateDimReduc_STonly/IleumOnly.rds')
rm(seu.integrated)

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
## [1] stats4    stats     graphics  grDevices utils     datasets  methods  
## [8] base     
## 
## other attached packages:
##  [1] scales_1.2.0                writexl_1.4.0              
##  [3] readxl_1.4.0                dplyr_1.0.9                
##  [5] SeuratDisk_0.0.0.9020       ggplot2_3.3.6              
##  [7] sp_1.5-0                    SeuratObject_4.1.0         
##  [9] Seurat_4.1.1                DropletUtils_1.16.0        
## [11] SingleCellExperiment_1.18.0 SummarizedExperiment_1.26.1
## [13] Biobase_2.56.0              GenomicRanges_1.48.0       
## [15] GenomeInfoDb_1.32.2         IRanges_2.30.0             
## [17] S4Vectors_0.34.0            BiocGenerics_0.42.0        
## [19] MatrixGenerics_1.8.1        matrixStats_0.62.0         
## 
## loaded via a namespace (and not attached):
##   [1] plyr_1.8.7                igraph_1.3.4             
##   [3] lazyeval_0.2.2            splines_4.2.1            
##   [5] BiocParallel_1.30.3       listenv_0.8.0            
##   [7] scattermore_0.8           digest_0.6.29            
##   [9] htmltools_0.5.3           fansi_1.0.3              
##  [11] magrittr_2.0.3            tensor_1.5               
##  [13] cluster_2.1.3             ROCR_1.0-11              
##  [15] limma_3.52.2              globals_0.15.1           
##  [17] R.utils_2.12.0            spatstat.sparse_2.1-1    
##  [19] colorspace_2.0-3          ggrepel_0.9.1            
##  [21] xfun_0.31                 crayon_1.5.1             
##  [23] RCurl_1.98-1.7            jsonlite_1.8.0           
##  [25] progressr_0.10.1          spatstat.data_2.2-0      
##  [27] survival_3.3-1            zoo_1.8-10               
##  [29] glue_1.6.2                polyclip_1.10-0          
##  [31] gtable_0.3.0              zlibbioc_1.42.0          
##  [33] XVector_0.36.0            leiden_0.4.2             
##  [35] DelayedArray_0.22.0       Rhdf5lib_1.18.2          
##  [37] future.apply_1.9.0        HDF5Array_1.24.1         
##  [39] abind_1.4-5               DBI_1.1.3                
##  [41] edgeR_3.38.1              spatstat.random_2.2-0    
##  [43] miniUI_0.1.1.1            Rcpp_1.0.9               
##  [45] viridisLite_0.4.0         xtable_1.8-4             
##  [47] reticulate_1.25           spatstat.core_2.4-4      
##  [49] dqrng_0.3.0               bit_4.0.4                
##  [51] htmlwidgets_1.5.4         httr_1.4.3               
##  [53] RColorBrewer_1.1-3        ellipsis_0.3.2           
##  [55] ica_1.0-3                 farver_2.1.1             
##  [57] pkgconfig_2.0.3           R.methodsS3_1.8.2        
##  [59] scuttle_1.6.2             uwot_0.1.11              
##  [61] deldir_1.0-6              locfit_1.5-9.6           
##  [63] utf8_1.2.2                labeling_0.4.2           
##  [65] tidyselect_1.1.2          rlang_1.0.4              
##  [67] reshape2_1.4.4            later_1.3.0              
##  [69] cellranger_1.1.0          munsell_0.5.0            
##  [71] tools_4.2.1               cli_3.3.0                
##  [73] generics_0.1.3            ggridges_0.5.3           
##  [75] evaluate_0.15             stringr_1.4.0            
##  [77] fastmap_1.1.0             goftest_1.2-3            
##  [79] yaml_2.3.5                bit64_4.0.5              
##  [81] knitr_1.39                fitdistrplus_1.1-8       
##  [83] purrr_0.3.4               RANN_2.6.1               
##  [85] nlme_3.1-158              pbapply_1.5-0            
##  [87] future_1.26.1             sparseMatrixStats_1.8.0  
##  [89] mime_0.12                 R.oo_1.25.0              
##  [91] hdf5r_1.3.5               compiler_4.2.1           
##  [93] rstudioapi_0.13           plotly_4.10.0            
##  [95] png_0.1-7                 spatstat.utils_2.3-1     
##  [97] tibble_3.1.7              stringi_1.7.8            
##  [99] highr_0.9                 RSpectra_0.16-1          
## [101] rgeos_0.5-9               lattice_0.20-45          
## [103] Matrix_1.4-1              vctrs_0.4.1              
## [105] pillar_1.8.0              lifecycle_1.0.1          
## [107] rhdf5filters_1.8.0        spatstat.geom_2.4-0      
## [109] lmtest_0.9-40             RcppAnnoy_0.0.19         
## [111] data.table_1.14.2         cowplot_1.1.1            
## [113] bitops_1.0-7              irlba_2.3.5              
## [115] httpuv_1.6.5              patchwork_1.1.1          
## [117] R6_2.5.1                  promises_1.2.0.1         
## [119] KernSmooth_2.23-20        gridExtra_2.3            
## [121] parallelly_1.32.0         codetools_0.2-18         
## [123] MASS_7.3-58               assertthat_0.2.1         
## [125] rhdf5_2.40.0              withr_2.5.0              
## [127] sctransform_0.3.3         GenomeInfoDbData_1.2.8   
## [129] mgcv_1.8-40               parallel_4.2.1           
## [131] rpart_4.1.16              grid_4.2.1               
## [133] beachmat_2.12.0           tidyr_1.2.0              
## [135] rmarkdown_2.14            DelayedMatrixStats_1.18.0
## [137] Rtsne_0.16                shiny_1.7.2
