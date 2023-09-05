### Load required software packages

library(Seurat)
library(ggplot2)

## Read in Seurat object
#Read in Seurat object of processed data:
seu <- readRDS('/home/Jayne.Wiarda/SI_PP_SC_ST/NormalizeIntegrateDimReduc_STonly/AllSamples.rds')

#See how many cells/genes in Seurat object:
seu

#See how many cells in each sample:
table(seu$SampleID)

## Define clusters
#Calculate PCs to use:
pct <- seu[["pca"]]@stdev / sum(seu[["pca"]]@stdev) * 100 # find standard deviation for each PC
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
ggsave('/home/Jayne.Wiarda/SI_PP_SC_ST/Figures/st_Clustering_0.jpeg')
PCdims <- 1:pcs # use the minimum PC from the quantitative method above to set the PCs for subsequent steps
rm(pct, cumu, co1, co2, pcs)

#Create clusters:
DefaultAssay(seu) <- 'integrated'
seu <- FindNeighbors(seu, 
                     dims = PCdims, 
                     verbose = FALSE) 
seu <- FindClusters(seu,
                    verbose = FALSE,
                    resolution = 0.08) 
DimPlot(seu, group.by = 'seurat_clusters', label = TRUE)
ggsave('/home/Jayne.Wiarda/SI_PP_SC_ST/Figures/st_Clustering_1.jpeg')
SpatialDimPlot(seu, ncol = 3)
ggsave('/home/Jayne.Wiarda/SI_PP_SC_ST/Figures/st_Clustering_2.jpeg')

#Rename clusters:
Idents(seu) <- seu$seurat_clusters
seu <- RenameIdents(seu, '0' = 'IFZ/PFZ')
seu <- RenameIdents(seu, '1' = 'Villus')
seu <- RenameIdents(seu, '2' = 'Crypt')
seu <- RenameIdents(seu, '3' = 'Follicle')
seu <- RenameIdents(seu, '4' = 'Muscularis')
seu$Region_Clust <- Idents(seu)

Idents(seu) <- seu$Region_Clust
SpatialDimPlot(seu, ncol = 3)
ggsave('/home/Jayne.Wiarda/SI_PP_SC_ST/Figures/st_Clustering_3.jpeg')

#### Save the Seurat object
saveRDS(seu, '/home/Jayne.Wiarda/SI_PP_SC_ST/AnnotatedSeurat/AllSamples_annotated_ST.rds')
rm(seu)

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
## [1] ggplot2_3.3.6      sp_1.5-0           SeuratObject_4.1.0 Seurat_4.1.1      
## 
## loaded via a namespace (and not attached):
##   [1] Rtsne_0.16            colorspace_2.0-3      deldir_1.0-6         
##   [4] ellipsis_0.3.2        ggridges_0.5.3        rstudioapi_0.13      
##   [7] spatstat.data_2.2-0   farver_2.1.1          leiden_0.4.2         
##  [10] listenv_0.8.0         ggrepel_0.9.1         fansi_1.0.3          
##  [13] codetools_0.2-18      splines_4.2.1         knitr_1.39           
##  [16] polyclip_1.10-0       jsonlite_1.8.0        ica_1.0-3            
##  [19] cluster_2.1.3         png_0.1-7             rgeos_0.5-9          
##  [22] uwot_0.1.11           shiny_1.7.2           sctransform_0.3.3    
##  [25] spatstat.sparse_2.1-1 compiler_4.2.1        httr_1.4.3           
##  [28] assertthat_0.2.1      Matrix_1.4-1          fastmap_1.1.0        
##  [31] lazyeval_0.2.2        cli_3.3.0             later_1.3.0          
##  [34] htmltools_0.5.3       tools_4.2.1           igraph_1.3.4         
##  [37] gtable_0.3.0          glue_1.6.2            RANN_2.6.1           
##  [40] reshape2_1.4.4        dplyr_1.0.9           Rcpp_1.0.9           
##  [43] scattermore_0.8       vctrs_0.4.1           nlme_3.1-158         
##  [46] progressr_0.10.1      lmtest_0.9-40         spatstat.random_2.2-0
##  [49] xfun_0.31             stringr_1.4.0         globals_0.15.1       
##  [52] mime_0.12             miniUI_0.1.1.1        lifecycle_1.0.1      
##  [55] irlba_2.3.5           goftest_1.2-3         future_1.26.1        
##  [58] MASS_7.3-58           zoo_1.8-10            scales_1.2.0         
##  [61] spatstat.core_2.4-4   promises_1.2.0.1      spatstat.utils_2.3-1 
##  [64] parallel_4.2.1        RColorBrewer_1.1-3    yaml_2.3.5           
##  [67] reticulate_1.25       pbapply_1.5-0         gridExtra_2.3        
##  [70] rpart_4.1.16          stringi_1.7.8         highr_0.9            
##  [73] rlang_1.0.4           pkgconfig_2.0.3       matrixStats_0.62.0   
##  [76] evaluate_0.15         lattice_0.20-45       ROCR_1.0-11          
##  [79] purrr_0.3.4           tensor_1.5            labeling_0.4.2       
##  [82] patchwork_1.1.1       htmlwidgets_1.5.4     cowplot_1.1.1        
##  [85] tidyselect_1.1.2      parallelly_1.32.0     RcppAnnoy_0.0.19     
##  [88] plyr_1.8.7            magrittr_2.0.3        R6_2.5.1             
##  [91] generics_0.1.3        DBI_1.1.3             withr_2.5.0          
##  [94] mgcv_1.8-40           pillar_1.8.0          fitdistrplus_1.1-8   
##  [97] survival_3.3-1        abind_1.4-5           tibble_3.1.7         
## [100] future.apply_1.9.0    KernSmooth_2.23-20    utf8_1.2.2           
## [103] spatstat.geom_2.4-0   plotly_4.10.0         rmarkdown_2.14       
## [106] grid_4.2.1            data.table_1.14.2     digest_0.6.29        
## [109] xtable_1.8-4          tidyr_1.2.0           httpuv_1.6.5         
## [112] munsell_0.5.0         viridisLite_0.4.0