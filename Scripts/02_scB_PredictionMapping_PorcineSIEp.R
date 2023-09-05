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

## Load query data

#We will treat our new ileum/jejunum scRNA-seq as the query dataset. 

#Load in h5 Seurat object of our ileum/jejunum data:
seu <- LoadH5Seurat('/home/Jayne.Wiarda/SI_PP_SC_ST/NormalizeIntegrateDimReduc_SConly/AllSamples.h5Seurat')

#Change default assay:
DefaultAssay(seu) <- 'RNA'

#How many cells, genes, gene assays, and dimensional reductions?
seu

## Convert query data gene names 
#This step isn't necessary since we have identical gene nomenclature used across the two separate scRNA-seq datasets.

## Subset to only epithelial cells:

#Take only epithelial cell data:
CellTypePredictions <- read_excel('/home/Jayne.Wiarda/SI_PP_SC_ST/MappingPrediction/WiardaIleumAtlas_CellTypePredictions.xlsx') # read in predictions from mapping to pig ileum atlas data
seu <- AddMetaData(object = seu, 
                  metadata = c(CellTypePredictions)) # add predicted IDs to metadata

Idents(seu) <- seu$predicted.id
seu <- subset(seu, idents = 'Epithelial cells') # subset only cells predicted as epithelial cells
seu

#Re-normalize query data:
seu.list <- SplitObject(seu, split.by = "orig.ident") # split by sample IDs
seu.list
#rm(seu)
for (i in 1:length(seu.list)) { # normalize data using SCTransform method
  seu.list[[i]] <- SCTransform(seu.list[[i]], 
                               return.only.var.genes = FALSE, 
                               verbose = TRUE) 
} # use SCT normalization since this was also used on reference data we will soon import

## Read in reference dataset:
ref <- LoadH5Seurat('/home/Jayne.Wiarda/scRNAseqEpSI/Seurat/GutEpOnlySubset.h5Seurat')
ref
DimPlot(ref, reduction = 'tsne', group.by = 'celltype')

#ref<- UpdateSeuratObject(ref) # skip since already a v4 Seurat object
#ref[['integrated']] <- as(object = ref[['integrated']] , Class = "SCTAssay")

DefaultAssay(ref) <- "integrated"

#Perform cell type predictions:
for(i in 1:length(seu.list)) {
  seu.list[[i]] <- RunPCA(seu.list[[i]], # have to calculate and store PCA to run MappingScore() function later
                          npcs = 23) # reduce to 23 since our smallest sample only has 24 cells
}

CellTypePredictions <- list()
MappingScores <- list()
for(i in 1:length(seu.list)) {
  anchors <- FindTransferAnchors(
    reference = ref,
    query = seu.list[[i]],
    #reduction = "cca", # opted to use cca since the method is recommended for cross-species mapping 
    reduction = 'pcaproject', # not using cca unless cross-species comparison is being performed
    dims = 1:23, # reduce to 23 since this is the number of PCs we calculated
    k.score = 23, # reduce to 23 since our smallest sample only has 24 cells
    normalization.method = "SCT",
    recompute.residuals = FALSE) # have to set to FALSE to use SCT in Seurat v4 
  predictions <- TransferData(anchorset = anchors, 
                              refdata = list(cell_type = ref$celltype), 
                              dims = 1:23,
                              weight.reduction = "pcaproject",
                              k.weight = 23)
  MapScores <- MappingScore(
    anchors = anchors,
    combined.object = anchors@object.list[[1]],
    query.neighbors =  slot(object = seu.list[[i]], name = "neighbors")[["seu.list_ref.nn"]],
    query.weights = Tool(object = seu.list[[i]], slot = "TransferData")$weights.matrix,
    query.embeddings = Embeddings(object = seu.list[[i]]),
    ref.embeddings = Embeddings(object = ref),
    nn.method = "annoy",
    kanchors = 23, # reduce to 23 since our smallest sample only has 24 cells
    ksmooth = 23, # reduce to 23 since our smallest sample only has 24 cells
    ndim = 23 # reduce to 23 since this is the number of PCs we calculated
    # n.trees = n.trees
  )
  CellTypePredictions[[i]] <- predictions
  MappingScores[[i]] <- MapScores
}


CellTypePredictions <- do.call(rbind, CellTypePredictions)
CellTypePredictions <- as.data.frame(CellTypePredictions)
MappingScores <- Reduce(c,MappingScores)
MappingScores <- as.data.frame(MappingScores)

#Save the mapping & prediction results:
CellTypePredictions$CellBarcodes <- rownames(CellTypePredictions)
MappingScores$CellBarcodes <- rownames(MappingScores)
write_xlsx(CellTypePredictions, '/home/Jayne.Wiarda/SI_PP_SC_ST/MappingPrediction/WiardaSIEpAtlas_CellTypePredictions.xlsx')
write_xlsx(MappingScores, '/home/Jayne.Wiarda/SI_PP_SC_ST/MappingPrediction/WiardaSIEpAtlas_MappingScores.xlsx')

#Incorporate cell prediction & mapping scores into original Seurat object of query data:
#seu <- LoadH5Seurat('/home/Jayne.Wiarda/SI_PP_SC_ST/NormalizeIntegrateDimReduc_SConly/AllSamples.h5Seurat')
DefaultAssay(seu) <- 'RNA'
seu <- AddMetaData(object = seu, 
                  metadata = c(CellTypePredictions, MappingScores))

#Visualize prediction & mapping results:
FeaturePlot(seu, reduction = 'tsne', cols = c('gold', 'red'), features = 'MappingScores', min.cutoff = 0.7) # do we see high mapping scores, indicating good representation of query data by reference dataset? (Answer: yes)

VlnPlot(seu, group.by = 'orig.ident', features = 'MappingScores') # any major differences in mapping scores between ileum and jejunum data?
# Since jejunum and ileum samples have similar mapping scores, suggests that the ileum atlas reference dataset annotations are cross-applicable to the jejunal samples used here

predict <- colnames(seu@meta.data %>% select(starts_with("prediction.score."))) # extract list of all prediction score meta data slots
FeaturePlot(seu, # plot all prediction scores
            features = c(predict[28:33]),
            reduction = 'tsne', # change to 'tsne' to overlay onto t-SNE plot
            ncol = 2) & 
  scale_color_gradientn( colours = c('grey90', 'blue3'),  limits = c(0, 1)) & 
  NoAxes() & NoLegend() 

VlnPlot(seu, group.by = 'orig.ident', features = 'prediction.score.max') # any major differences in prediction scores between ileum and jejunum data? 
VlnPlot(seu, group.by = 'predicted.id', features = 'prediction.score.max') # any major differences in prediction scores between ileum and jejunum data? 
# See mostly high prediction scores, giving confidence to cell type identities assigned via reference annotation-based prediction
DimPlot(seu, group.by = 'predicted.id', reduction = 'tsne', label = FALSE) # plot predicted cell IDs (based on highest prediction scores)

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
##   [1] ggbeeswarm_0.6.0      Rtsne_0.16            colorspace_2.0-3     
##   [4] deldir_1.0-6          ellipsis_0.3.2        ggridges_0.5.3       
##   [7] rstudioapi_0.13       spatstat.data_2.2-0   farver_2.1.1         
##  [10] leiden_0.4.2          listenv_0.8.0         bit64_4.0.5          
##  [13] ggrepel_0.9.1         fansi_1.0.3           codetools_0.2-18     
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
##  [70] pbapply_1.5-0         gridExtra_2.3         ggrastr_1.0.1        
##  [73] rpart_4.1.16          stringi_1.7.8         highr_0.9            
##  [76] rlang_1.0.4           pkgconfig_2.0.3       matrixStats_0.62.0   
##  [79] evaluate_0.15         lattice_0.20-45       ROCR_1.0-11          
##  [82] purrr_0.3.4           tensor_1.5            labeling_0.4.2       
##  [85] patchwork_1.1.1       htmlwidgets_1.5.4     bit_4.0.4            
##  [88] cowplot_1.1.1         tidyselect_1.1.2      parallelly_1.32.1    
##  [91] RcppAnnoy_0.0.19      plyr_1.8.7            magrittr_2.0.3       
##  [94] R6_2.5.1              generics_0.1.3        DBI_1.1.3            
##  [97] withr_2.5.0           mgcv_1.8-40           pillar_1.8.0         
## [100] fitdistrplus_1.1-8    survival_3.3-1        abind_1.4-5          
## [103] tibble_3.1.8          future.apply_1.9.0    hdf5r_1.3.5          
## [106] crayon_1.5.1          KernSmooth_2.23-20    utf8_1.2.2           
## [109] spatstat.geom_2.4-0   plotly_4.10.0         rmarkdown_2.14       
## [112] grid_4.2.1            data.table_1.14.2     digest_0.6.29        
## [115] xtable_1.8-4          tidyr_1.2.0           httpuv_1.6.5         
## [118] munsell_0.5.0         beeswarm_0.4.0        viridisLite_0.4.0    
## [121] vipor_0.4.5
