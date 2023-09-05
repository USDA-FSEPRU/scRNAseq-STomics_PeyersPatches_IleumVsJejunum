library(Seurat)
library(SeuratObject)
library(SeuratDisk)
library(ggplot2)

## Create Seurat object
#Load the spatial data as Seurat objects and merge:
I2 <- Load10X_Spatial(data.dir = '/home/Jayne.Wiarda/SI_PP_SC_ST/SpaceRangerOutputs/I2', 
                      slice = 'I2')
I3 <- Load10X_Spatial(data.dir = '/home/Jayne.Wiarda/SI_PP_SC_ST/SpaceRangerOutputs/I3', 
                      slice = 'I3')
I4 <- Load10X_Spatial(data.dir = '/home/Jayne.Wiarda/SI_PP_SC_ST/SpaceRangerOutputs/I4', 
                      slice = 'I4')
J2 <- Load10X_Spatial(data.dir = '/home/Jayne.Wiarda/SI_PP_SC_ST/SpaceRangerOutputs/J2', 
                      slice = 'J2')
J3 <- Load10X_Spatial(data.dir = '/home/Jayne.Wiarda/SI_PP_SC_ST/SpaceRangerOutputs/J3', 
                      slice = 'J3')
J4 <- Load10X_Spatial(data.dir = '/home/Jayne.Wiarda/SI_PP_SC_ST/SpaceRangerOutputs/J4', 
                      slice = 'J4')
seu <- merge(I2, y = c(I3, I4, J2, J3, J4),
             add.cell.ids = c('I2', 'I3', 'I4', 'J2', 'J3', 'J4'))
rm(I2, I3, I4, J2, J3, J4) # free up space

#See how many cells/genes in Seurat object:
seu

#Also add sample ID in metadata slot:
seu$SampleID <- substr(colnames(seu),start=1,stop=2)
Idents(seu) <- seu$SampleID # set as default identity for later plotting

## Identify any genes with sum zero counts
#Remove non-expressed genes:
keep <- rowSums(seu) > 0 # specify rows with genes that are expressed at least once in any sample
seu <- seu[keep, ] # retain only genes that are expressed at least once
rm(keep) # free up space

#See how many spots/genes in Seurat object:
seu
table(seu$SampleID)

## Plot QC metrics
#Violin plots:
VlnPlot(seu, 
        features = c("nFeature_Spatial", "nCount_Spatial"), # QC metric metadata to plot
        split.by = 'SampleID', # plot each sample separately
        pt.size = 0, # don't include individual cells as points on plots since there are so many
        ncol = 2)

#Number of genes vs. number reads:
meta <- seu@meta.data
ggplot(meta, aes(x=nCount_Spatial, y=nFeature_Spatial, color = SampleID))+
  geom_point() + 
  facet_wrap(~SampleID, nrow =1)+
  theme_get() + 
  xlim(0,70000) +
  ylim(0,8000) +
  ylab("#Genes detected per cell") +
  RotatedAxis()

#Histograms:
ggplot(meta, aes(x=nFeature_Spatial,y=..density..)) + 
  geom_histogram(fill="white",color="black",bins=500) + 
  scale_x_continuous(breaks = seq(0, 8000, 250), lim = c(0, 8000)) + 
  facet_wrap(~SampleID) +
  RotatedAxis() + 
  ggtitle('nFeature_Spatial')

ggplot(meta, aes(x=nCount_Spatial,y=..density..)) + 
  geom_histogram(fill="white",color="black",bins=500) + 
  scale_x_continuous(breaks = seq(0, 70000, 2500), lim = c(0, 70000)) + 
  facet_wrap(~SampleID) +
  RotatedAxis() + 
  ggtitle('nCount_Spatial')

# We see some large differences in distribution of values for the QC metrics across different samples, so next we will overlay the QC values onto each sample's H&E image:
SpatialPlot(seu, 
            features = c("nFeature_Spatial", "nCount_Spatial"))

#We see from this plot that the muscularis layer has lower gene expression; however, we also see that there are areas outside muscularis, such as a chunk of Peyer's patch region in J3, that also have low RNA counts, likely due to the tissue not being properly adhered to slide in that area during the STomics assay. Our goal is to find a way to filter out the areas with poor tissue adherence (technical low RNA count issues) while retaining low-count muscularis regions (biological low RNA count).

## Define a threshold for filtering out low UMI spots
#Based on plots above, we can try some different low intervals for filtering and see how they look on plots.

#Try 1000 UMI threshold:
cells <- WhichCells(seu, expression = nCount_Spatial < 1000)
SpatialPlot(seu, 
            cells.highlight = cells) & # cells with low UMI will be highlighted in red (default color)
  NoLegend()

#Try 3000 UMI threshold:
cells <- WhichCells(seu, expression = nCount_Spatial < 3000)
SpatialPlot(seu, 
            cells.highlight = cells) & # cells with low UMI will be highlighted in red (default color)
  NoLegend()

#Try 5000 UMI threshold:
cells <- WhichCells(seu, expression = nCount_Spatial < 5000)
SpatialPlot(seu, 
            cells.highlight = cells) & # cells with low UMI will be highlighted in red (default color)
  NoLegend()

#Try 7000 UMI threshold:
cells <- WhichCells(seu, expression = nCount_Spatial < 7000)
SpatialPlot(seu, 
            cells.highlight = cells) & # cells with low UMI will be highlighted in red (default color)
  NoLegend()

#Try 9000 UMI threshold:
cells <- WhichCells(seu, expression = nCount_Spatial < 9000)
SpatialPlot(seu, 
            cells.highlight = cells) & # cells with low UMI will be highlighted in red (default color)
  NoLegend()

#Based on comparing plots, we decided to go with a threshold of 5000 UMIs for filtering out low UMI spots:
ggplot(meta, aes(x=nCount_Spatial,y=..density..)) + 
  geom_histogram(fill="white",color="black",bins=500) + 
  scale_x_continuous(breaks = seq(0, 70000, 2500), lim = c(0, 70000)) + 
  facet_wrap(~SampleID) +
  geom_vline(aes(xintercept=5000),color="red",lty="longdash") + # move this cutoff line where you see fit
  RotatedAxis() + 
  ggtitle('nCount_Spatial')
rm(meta)
  
## Filter our spots with low total UMIs
#Identify cells passing each/every QC filter:
keep <- WhichCells(seu, expression = nCount_Spatial > 5000) # cells passing UMI filter
length(keep) # how many cells pass this QC filter?

#Create new Seurat object with only cells passing all QC filters:
seuKeep <- subset(seu, 
                  cells = keep)
rm(seu) # free up space

#How many cells in each sample now?
table(seuKeep$SampleID)

## Look over plots with only cells passing QC
#Violin plots:
VlnPlot(seuKeep, 
        features = c("nFeature_Spatial", "nCount_Spatial"), # QC metric metadata to plot
        split.by = 'SampleID', # plot each sample separately
        pt.size = 0, # don't include individual cells as points on plots since there are so many
        ncol = 2)

#Number of genes vs. number reads:
meta <- seuKeep@meta.data
ggplot(meta, aes(x=nCount_Spatial, y=nFeature_Spatial, color = SampleID))+
  geom_point() + 
  facet_wrap(~SampleID, nrow =1)+
  theme_get() + 
  xlim(0,70000) +
  ylim(0,8000) +
  ylab("#Genes detected per cell") +
  RotatedAxis()

#Histograms:
ggplot(meta, aes(x=nFeature_Spatial,y=..density..)) + 
  geom_histogram(fill="white",color="black",bins=500) + 
  scale_x_continuous(breaks = seq(0, 8000, 250), lim = c(0, 8000)) + 
  facet_wrap(~SampleID) +
  RotatedAxis() + 
  ggtitle('nFeature_Spatial')

ggplot(meta, aes(x=nCount_Spatial,y=..density..)) + 
  geom_histogram(fill="white",color="black",bins=500) + 
  scale_x_continuous(breaks = seq(0, 70000, 2500), lim = c(0, 70000)) + 
  facet_wrap(~SampleID) +
  RotatedAxis() + 
  ggtitle('nCount_Spatial')
rm(meta)

#Overlay the QC values onto each sample's H&E image:
SpatialPlot(seuKeep, 
            features = c("nFeature_Spatial", "nCount_Spatial"))

## Save counts of spots passing QC
saveRDS(seuKeep, 
             '/home/Jayne.Wiarda/SI_PP_SC_ST/QC_STonly/AllSamples_PostQC.rds')

### View session information
sessionInfo()
## R version 4.1.2 (2021-11-01)
## Platform: x86_64-pc-linux-gnu (64-bit)
## Running under: Ubuntu 20.04.3 LTS
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
## [1] ggplot2_3.3.5         SeuratDisk_0.0.0.9019 SeuratObject_4.0.2   
## [4] Seurat_4.0.4         
## 
## loaded via a namespace (and not attached):
##   [1] Rtsne_0.15            colorspace_2.0-2      deldir_0.2-10        
##   [4] ellipsis_0.3.2        ggridges_0.5.3        rstudioapi_0.13      
##   [7] spatstat.data_2.1-0   farver_2.1.0          leiden_0.3.9         
##  [10] listenv_0.8.0         ggrepel_0.9.1         bit64_4.0.5          
##  [13] fansi_0.5.0           codetools_0.2-18      splines_4.1.2        
##  [16] knitr_1.34            polyclip_1.10-0       jsonlite_1.7.2       
##  [19] ica_1.0-2             cluster_2.1.2         png_0.1-7            
##  [22] uwot_0.1.10           shiny_1.6.0           sctransform_0.3.2    
##  [25] spatstat.sparse_2.0-0 compiler_4.1.2        httr_1.4.2           
##  [28] assertthat_0.2.1      Matrix_1.3-4          fastmap_1.1.0        
##  [31] lazyeval_0.2.2        cli_3.0.1             later_1.3.0          
##  [34] htmltools_0.5.2       tools_4.1.2           igraph_1.2.6         
##  [37] gtable_0.3.0          glue_1.4.2            RANN_2.6.1           
##  [40] reshape2_1.4.4        dplyr_1.0.7           Rcpp_1.0.7           
##  [43] scattermore_0.7       vctrs_0.3.8           nlme_3.1-152         
##  [46] lmtest_0.9-38         xfun_0.26             stringr_1.4.0        
##  [49] globals_0.14.0        mime_0.11             miniUI_0.1.1.1       
##  [52] lifecycle_1.0.0       irlba_2.3.3           goftest_1.2-2        
##  [55] future_1.22.1         MASS_7.3-54           zoo_1.8-9            
##  [58] scales_1.1.1          spatstat.core_2.3-0   promises_1.2.0.1     
##  [61] spatstat.utils_2.2-0  parallel_4.1.2        RColorBrewer_1.1-2   
##  [64] yaml_2.2.1            reticulate_1.21       pbapply_1.5-0        
##  [67] gridExtra_2.3         rpart_4.1-15          stringi_1.7.4        
##  [70] highr_0.9             rlang_0.4.11          pkgconfig_2.0.3      
##  [73] matrixStats_0.60.1    evaluate_0.14         lattice_0.20-45      
##  [76] ROCR_1.0-11           purrr_0.3.4           tensor_1.5           
##  [79] labeling_0.4.2        patchwork_1.1.1       htmlwidgets_1.5.4    
##  [82] cowplot_1.1.1         bit_4.0.4             tidyselect_1.1.1     
##  [85] parallelly_1.28.1     RcppAnnoy_0.0.19      plyr_1.8.6           
##  [88] magrittr_2.0.1        R6_2.5.1              generics_0.1.0       
##  [91] DBI_1.1.1             pillar_1.6.2          withr_2.4.2          
##  [94] mgcv_1.8-38           fitdistrplus_1.1-5    survival_3.2-13      
##  [97] abind_1.4-5           tibble_3.1.4          future.apply_1.8.1   
## [100] crayon_1.4.1          hdf5r_1.3.4           KernSmooth_2.23-20   
## [103] utf8_1.2.2            spatstat.geom_2.2-2   plotly_4.9.4.1       
## [106] rmarkdown_2.11        grid_4.1.2            data.table_1.14.0    
## [109] digest_0.6.27         xtable_1.8-4          tidyr_1.1.3          
## [112] httpuv_1.6.3          munsell_0.5.0         viridisLite_0.4.0