library(Seurat)
library(SeuratObject)
library(ggplot2)
library(SeuratDisk)

## Load & process Seurat object

#Load a processed Seurat object of only B cells (from both ileum + jejunum samples combined):
  
seu <- LoadH5Seurat('/home/Jayne.Wiarda/SI_PP_SC_ST/AnnotatedSeurat/B.h5seurat')
DefaultAssay(seu) <- 'SCT'
Idents(seu) <- seu$celltype

#Make a new metadata slot for tissue:
  
seu$tissue <- substr(seu$orig.ident, 1, 1) # I = ileum, J = jejunum

#Add cell cycle information to cells:
  
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
DefaultAssay(seu) <- 'SCT'
seu <- CellCycleScoring(seu, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE) # don't humanize gene names since we did a manual check that all gene symbols lined up or aren't present in our pig annotation
DimPlot(seu, reduction = 'pca')
DimPlot(seu, reduction = 'umap')
FeaturePlot(seu, features = c('G2M.Score', 'S.Score'), cols = c('grey90', 'navy'), reduction = 'pca')
FeaturePlot(seu, features = c('G2M.Score', 'S.Score'), cols = c('grey90', 'navy'), reduction = 'umap')

# Save data:
SaveH5Seurat(seu, '/home/Jayne.Wiarda/SI_PP_SC_ST/AnnotatedSeurat/B.h5seurat', overwrite = TRUE)

sessionInfo()
#R version 4.2.2 Patched (2022-11-10 r83330)
#Platform: x86_64-pc-linux-gnu (64-bit)
#Running under: Ubuntu 22.04.1 LTS

#Matrix products: default
#BLAS:   /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3
#LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/libopenblasp-r0.3.20.so

#locale:
#  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C               LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8     LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8    LC_PAPER=en_US.UTF-8       LC_NAME=C                 
#[9] LC_ADDRESS=C               LC_TELEPHONE=C             LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       

#attached base packages:
#  [1] stats4    stats     graphics  grDevices utils     datasets  methods   base     

#loaded via a namespace (and not attached):
#  [1] scattermore_0.8             R.methodsS3_1.8.2           ragg_1.2.2                  SeuratObject_4.1.1          pkgmaker_0.32.2             tidyr_1.2.1                 ggplot2_3.3.6              
#[8] bit64_4.0.5                 knitr_1.40                  irlba_2.3.5                 DelayedArray_0.22.0         R.utils_2.12.0              data.table_1.14.2           rpart_4.1.16               
#[15] KEGGREST_1.36.3             RCurl_1.98-1.8              doParallel_1.0.17           generics_0.1.3              BiocGenerics_0.42.0         ScaledMatrix_1.4.1          cowplot_1.1.1              
#[22] RSQLite_2.2.17              RANN_2.6.1                  future_1.28.0               bit_4.0.4                   phylobase_0.8.10            spatstat.data_2.2-0         xml2_1.3.3                 
#[29] httpuv_1.6.6                SummarizedExperiment_1.26.1 assertthat_0.2.1            viridis_0.6.2               xfun_0.33                   hms_1.1.2                   evaluate_0.16              
#[36] promises_1.2.0.1            fansi_1.0.3                 progress_1.2.2              igraph_1.3.4                DBI_1.1.3                   htmlwidgets_1.5.4           spatstat.geom_2.4-0        
#[43] purrr_0.3.4                 ellipsis_0.3.2              dplyr_1.0.10                V8_4.2.2                    annotate_1.74.0             gridBase_0.4-7              locfdr_1.1-8               
#[50] deldir_1.0-6                sparseMatrixStats_1.8.0     MatrixGenerics_1.8.1        vctrs_0.4.1                 SingleCellExperiment_1.18.0 Biobase_2.56.0              ROCR_1.0-11                
#[57] abind_1.4-5                 cachem_1.0.6                withr_2.5.0                 ggforce_0.3.4               progressr_0.11.0            sctransform_0.3.4           prettyunits_1.1.1          
#[64] goftest_1.2-3               softImpute_1.4-1            cluster_2.1.4               ape_5.6-2                   lazyeval_0.2.2              crayon_1.5.1                genefilter_1.78.0          
#[71] hdf5r_1.3.5                 labeling_0.4.2              edgeR_3.38.4                pkgconfig_2.0.3             tweenr_2.0.2                GenomeInfoDb_1.32.4         nlme_3.1-159               
#[78] vipor_0.4.5                 rlang_1.0.5                 globals_0.16.1              lifecycle_1.0.2             miniUI_0.1.1.1              registry_0.5-1              rsvd_1.0.5                 
#[85] dichromat_2.0-0.1           cellranger_1.1.0            polyclip_1.10-0             matrixStats_0.62.0          lmtest_0.9-40               rngtools_1.5.2              Matrix_1.5-1               
#[92] Rhdf5lib_1.18.2             zoo_1.8-10                  beeswarm_0.4.0              ggridges_0.5.3              png_0.1-7                   viridisLite_0.4.1           bitops_1.0-7               
#[99] R.oo_1.25.0                 rncl_0.8.6                  KernSmooth_2.23-20          rhdf5filters_1.8.0          Biostrings_2.64.1           blob_1.2.3                  DelayedMatrixStats_1.18.0  
#[106] stringr_1.4.1               zinbwave_1.18.0             parallelly_1.32.1           spatstat.random_2.2-0       S4Vectors_0.34.0            beachmat_2.12.0             scales_1.2.1               
#[113] memoise_2.0.1               magrittr_2.0.3              plyr_1.8.7                  ica_1.0-3                   howmany_0.3-1               zlibbioc_1.42.0             compiler_4.2.2             
#[120] dqrng_0.3.0                 RColorBrewer_1.1-3          fitdistrplus_1.1-8          cli_3.4.0                   ade4_1.7-19                 XVector_0.36.0              listenv_0.8.0              
#[127] patchwork_1.1.2             pbapply_1.5-0               MASS_7.3-58.1               mgcv_1.8-40                 tidyselect_1.1.2            stringi_1.7.8               textshaping_0.3.6          
#[134] yaml_2.3.5                  BiocSingular_1.12.0         locfit_1.5-9.6              ggrepel_0.9.1               grid_4.2.2                  tools_4.2.2                 future.apply_1.9.1         
#[141] parallel_4.2.2              rstudioapi_0.14             uuid_1.1-0                  foreach_1.5.2               RNeXML_2.4.7                gridExtra_2.3               farver_2.1.1               
#[148] Rtsne_0.16                  ggraph_2.0.6                digest_0.6.29               rgeos_0.5-9                 shiny_1.7.2                 Rcpp_1.0.9                  GenomicRanges_1.48.0       
#[155] scuttle_1.6.3               later_1.3.0                 RcppAnnoy_0.0.19            AnnotationDbi_1.58.0        httr_1.4.4                  kernlab_0.9-31              colorspace_2.0-3           
#[162] XML_3.99-0.10               tensor_1.5                  reticulate_1.26             clusterExperiment_2.16.0    IRanges_2.30.1              splines_4.2.2               uwot_0.1.14                
#[169] spatstat.utils_2.3-1        graphlayouts_0.8.1          sp_1.5-0                    mapproj_1.2.8               systemfonts_1.0.4           plotly_4.10.0               xtable_1.8-4               
#[176] jsonlite_1.8.0              tidygraph_1.2.2             R6_2.5.1                    pillar_1.8.1                htmltools_0.5.3             mime_0.12                   NMF_0.24.0                 
#[183] glue_1.6.2                  fastmap_1.1.0               BiocParallel_1.30.3         BiocNeighbors_1.14.0        codetools_0.2-18            maps_3.4.0                  utf8_1.2.2                 
#[190] lattice_0.20-45             spatstat.sparse_2.1-1       tibble_3.1.8                curl_4.3.2                  ggbeeswarm_0.6.0            leiden_0.4.3                gtools_3.9.3               
#[197] survival_3.4-0              limma_3.52.3                rmarkdown_2.16              munsell_0.5.0               rhdf5_2.40.0                GenomeInfoDbData_1.2.8      iterators_1.0.14           
#[204] HDF5Array_1.24.2            reshape2_1.4.4              gtable_0.3.1                spatstat.core_2.4-4   