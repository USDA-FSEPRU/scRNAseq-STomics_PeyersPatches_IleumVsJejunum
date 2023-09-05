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
CellTypePredictions <- read_excel('/home/Jayne.Wiarda/SI_PP_SC_ST/MappingPrediction/WiardaIleumAtlas_CellTypePredictions.xlsx') # read in predictions from mapping to pig ileum atlas data

## Add metadata:
seu <- AddMetaData(object = seu, 
                   metadata = c(CellTypePredictions))

## Further define stromal cells
#Create dataset of only stromal cells:
Idents(seu) <- seu$predicted.id
str <- subset(seu, idents = 'Stromal cells')
DefaultAssay(str) <- 'RNA'
counts <- as.data.frame(str[['RNA']]@counts)
keep <- rowSums(counts) > 0
keep <- rownames(counts[keep,])
str <- DietSeurat(str, 
                  counts = TRUE,
                  data = TRUE,
                  scale.data = FALSE, # remove the scaled data
                  dimreducs = NULL,
                  features = keep, # keep only genes with non-zero counts across all cells
                  assays = 'RNA') # keep only RNA assay and remove SCT and integrated

seu.list <- SplitObject(str, split.by = "orig.ident") # split by sample IDs
seu.list
rm(str)
for (i in 1:length(seu.list)) { # normalize data using SCTransform method
  seu.list[[i]] <- SCTransform(seu.list[[i]], 
                               return.only.var.genes = FALSE, 
                               verbose = TRUE) 
} # use SCT normalization since this was also used on reference data we will soon import

seu.features <- SelectIntegrationFeatures(seu.list, # select the genes to use for integration
                                          verbose = TRUE) 
seu.list <- PrepSCTIntegration(seu.list, 
                               anchor.features = seu.features,
                               verbose = TRUE)
seu.anchors <- FindIntegrationAnchors(seu.list, # identify anchors for integration from top 30 data dimensions
                                      normalization.method = "SCT", 
                                      anchor.features = seu.features, 
                                      k.filter = 11, # reduce to 11 since smallest sample has 12 cells
                                      k.score = 11, # reduce to 11 since smallest sample has 12 cells
                                      dims = 1:11) # reduce to 11 since smallest sample has 12 cells
seu.integrated <- IntegrateData(seu.anchors, # integrate data
                                normalization.method = "SCT", 
                                k.weight = 11, # reduce to 11 since smallest sample has 12 cells
                                dims = 1:11) # reduce to 11 since smallest sample has 12 cells
rm(counts, keep, seu.features, seu.anchors)
seu.integrated <- RunPCA(seu.integrated, # run PCA analysis for 50 dimensions of the data
                         npcs = 50, 
                         verbose = TRUE) 

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
#ggplot(plot_df, aes(cumu, pct, label = rank, color = rank > pcs)) + # visualize PCs to use in elbow plot
#  geom_text() + 
#  geom_vline(xintercept = 90, color = "grey") + 
#  geom_hline(yintercept = min(pct[pct > 5]), color = "grey") +
#  theme_bw()
PCdims <- 1:pcs # use the minimum PC from the quantitative method above to set the PCs for subsequent steps
rm(pct, cumu, co1, co2, pcs)

seu.integrated <- RunTSNE(seu.integrated, 
                          dims = PCdims, # use our calculated number of PCs
                          reduction = "pca", 
                          dim_embed = 3, # calculate 3 plot dimensions in case we want to try 3D plotting later
                          assay = "SCT") # create tSNE plot 
seu.integrated <- RunUMAP(seu.integrated, 
                          dims = PCdims, # use our calculated number of PCs
                          reduction = "pca", 
                          dim_embed = 3, # calculate 3 plot dimensions in case we want to try 3D plotting later
                          min.dist = 0.5,
                          spread = 0.2,
                          assay = "SCT") # create UMAP
#dim(seu.integrated[["RNA"]]@scale.data) # see that there is no RNA assay scaled data yet
seu.integrated <- NormalizeData(seu.integrated,  # normalize the RNA counts data per cell
                                normalization.method = "LogNormalize", 
                                scale.factor = 10000, 
                                assay = "RNA")
seu.integrated <- ScaleData(seu.integrated, # scale the RNA counts data relative to other cells
                            assay = "RNA") # scales all genes instead of just highly variable
seu.integrated <- ScaleData(seu.integrated, # scale the SCT counts data relative to other cells
                            assay = "SCT")
#dim(seu.integrated[["RNA"]]@scale.data) # see that all genes are scaled in RNA assay now
#dim(seu.integrated[["SCT"]]@scale.data) # see that all genes are scaled in SCT assay now

#Define stromal cell clusters:
seu.integrated <- FindNeighbors(seu.integrated, 
                                dims = PCdims, 
                                verbose = FALSE) 
seu.integrated <- FindClusters(seu.integrated, 
                               verbose = FALSE,
                               resolution = 0.8) 
DimPlot(seu.integrated, group.by = 'seurat_clusters', label = TRUE)

##Query expression of canonical genes:
  
#Canonical genes for endothelial cells: PECAM1, CDH5
#Canonical genes for muscle cells: MYH11, ACTG2, TAGLN
#Canonical genes for fibroblasts: COL1A1, COL1A2, ECM1

DefaultAssay(seu.integrated) <- 'SCT'
FeaturePlot(seu.integrated,
            reduction = 'umap',
            features = c('PECAM1', 'CDH5', # endothelial markers
                         'MYH11', 'ACTG2', 'TAGLN', # muscle markers
                         'COL1A1', 'COL1A2', 'ECM1')) # fibroblast markers

#Rename clusters:
Idents(seu.integrated) <- seu.integrated$seurat_clusters
seu.integrated <- RenameIdents(seu.integrated, '0' = 'Fibroblasts')
seu.integrated <- RenameIdents(seu.integrated, '1' = 'Fibroblasts')
seu.integrated <- RenameIdents(seu.integrated, '2' = 'Endothelial cells')
seu.integrated <- RenameIdents(seu.integrated, '3' = 'Fibroblasts')
seu.integrated <- RenameIdents(seu.integrated, '4' = 'Endothelial cells')
seu.integrated <- RenameIdents(seu.integrated, '5' = 'Muscle cells')
seu.integrated <- RenameIdents(seu.integrated, '6' = 'Fibroblasts')
seu.integrated$stromal <- Idents(seu.integrated)
seu.integrated$celltype <- seu.integrated$stromal

#Reset identifity levels:
Idents(seu.integrated) <- seu.integrated$celltype
levels(seu.integrated) <- c('Endothelial cells', 'Fibroblasts', 'Muscle cells')
seu.integrated$celltype <- Idents(seu.integrated)

#Plot new annotations
DimPlot(seu.integrated, group.by = 'celltype', label = TRUE)

#Save data:
SaveH5Seurat(seu.integrated, '/home/Jayne.Wiarda/SI_PP_SC_ST/AnnotatedSeurat/StromalCells.h5seurat', overwrite = TRUE)

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
