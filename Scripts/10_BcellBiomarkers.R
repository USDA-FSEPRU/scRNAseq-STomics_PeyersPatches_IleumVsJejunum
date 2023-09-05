library(Seurat)
library(SeuratDisk)

bcell <- LoadH5Seurat('/home/Jayne.Wiarda/SI_PP_SC_ST/PseudotimeTrajectory/B_traj.h5seurat')
DefaultAssay(bcell) -> 'SCT'

seu <- LoadH5Seurat('/home/Jayne.Wiarda/SI_PP_SC_ST/AnnotatedSeurat/AllSamples_annotated.h5seurat')
DefaultAssay(seu) <- 'SCT'

st <- readRDS('/home/Jayne.Wiarda/SI_PP_SC_ST/AnnotatedSeurat/AllSamples_annotated_ST.rds')
DefaultAssay(st) <- 'SCT'
st$tissue <- substr(st$SampleID, 1, 1)   

DotPlot(bcell,
        features = c('PAX5', 'CD79A', 'BCL6', 'MKI67', 'IRF4', 'PRDM1', 'RPL5'), # add last gene so we will scale to 100%
        cols = c('gold', 'red'),
        group.by = 'celltype') + RotatedAxis()

bcell$psuedotime1_bin  <- factor(bcell$psuedotime1_bin, 
                                 levels = rev(c( 'Non-trajectory B', '(0,33.5]', '(33.5,67]', '(67,101]',
                                            '(101,134]', '(134,168]')))
DotPlot(bcell,
        features = c('PAX5', 'CD79A', 'BCL6', 'MKI67', 'IRF4', 'PRDM1', 'RPL5'), # add last gene so we will scale to 100%
        cols = c('gold', 'red'),
        group.by = 'psuedotime1_bin',
        col.min = -1,
        col.max = 1.5) + RotatedAxis()

bcell$psuedotime2_bin  <- factor(bcell$psuedotime2_bin, 
                                 levels = rev(c('Non-trajectory B', '(0,54.5]', '(54.5,109]', '(109,164]',
                                            '(164,218]', '(218,273]')))
DotPlot(bcell,
        features = c('PAX5', 'CD79A', 'BCL6', 'MKI67', 'IRF4', 'PRDM1', 'RPL5'), # add last gene so we will scale to 100%
        cols = c('gold', 'red'),
        group.by = 'psuedotime2_bin',
        col.min = -1,
        col.max = 1.5) + RotatedAxis()

DotPlot(seu,
        features = c('PAX5', 'CD79A', 'BCL6', 'MKI67', 'IRF4', 'PRDM1', 'RPL5'),
        cols = c('gold', 'red'),
        group.by = 'celltype') + RotatedAxis()

st$Region_Clust <- factor(st$Region_Clust, 
                          levels = rev(c('Villus', 'Crypt', 'IFZ/PFZ', 'Follicle', 'Muscularis')))
DotPlot(st,
        features = c('PAX5', 'CD79A', 'BCL6', 'MKI67', 'IRF4', 'PRDM1', 'RPL5'),
        cols = c('gold', 'red'),
        group.by = 'Region_Clust',
        col.min = -1,
        col.max = 1.5) + RotatedAxis()

st$Region_ManAnn <- factor(st$Region_ManAnn, 
                          levels = rev(c('Villus', 'Crypt', 'IFZ/PFZ', 'Follicle', 'Muscularis')))
DotPlot(st,
        features = c('PAX5', 'CD79A', 'BCL6', 'MKI67', 'IRF4', 'PRDM1', 'RPL5'),
        cols = c('gold', 'red'),
        group.by = 'Region_ManAnn',
        col.min = -1,
        col.max = 1.5) + RotatedAxis()

st$comboClus <- paste(st$Region_Clust, st$tissue, sep = '_')
DotPlot(st,
        features = c('PAX5', 'CD79A', 'BCL6', 'MKI67', 'IRF4', 'PRDM1', 'RPL5'),
        cols = c('gold', 'red'),
        group.by = 'comboClus',
        col.min = -1,
        col.max = 1.5) + RotatedAxis()

st$comboMan <- paste(st$Region_ManAnn, st$tissue, sep = '_')
DotPlot(st,
        features = c('PAX5', 'CD79A', 'BCL6', 'MKI67', 'IRF4', 'PRDM1', 'RPL5'),
        cols = c('gold', 'red'),
        group.by = 'comboMan',
        col.min = -1,
        col.max = 1.5) + RotatedAxis()

sessionInfo()
#R version 4.2.2 Patched (2022-11-10 r83330)
#Platform: x86_64-pc-linux-gnu (64-bit)
#Running under: Ubuntu 22.04.2 LTS

#Matrix products: default
#BLAS:   /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3
#LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/libopenblasp-r0.3.20.so

#locale:
#  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C               LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8     LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8    LC_PAPER=en_US.UTF-8       LC_NAME=C                 
#[9] LC_ADDRESS=C               LC_TELEPHONE=C             LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       

#attached base packages:
#  [1] stats4    stats     graphics  grDevices utils     datasets  methods   base     

#other attached packages:
#  [1] SeuratDisk_0.0.0.9020 sp_1.5-0              SeuratObject_4.1.1    Seurat_4.1.1         

#loaded via a namespace (and not attached):
#  [1] rappdirs_0.3.3         scattermore_0.8        coda_0.19-4            pkgmaker_0.32.2        tidyr_1.2.1            knitr_1.40             ggplot2_3.3.6          bit64_4.0.5            irlba_2.3.5           
#[10] data.table_1.14.2      rpart_4.1.16           KEGGREST_1.36.3        RCurl_1.98-1.8         doParallel_1.0.17      generics_0.1.3         BiocGenerics_0.42.0    cowplot_1.1.1          RSQLite_2.2.17        
#[19] RANN_2.6.1             future_1.28.0          bit_4.0.4              tzdb_0.3.0             spatstat.data_2.2-0    xml2_1.3.3             lubridate_1.9.0        httpuv_1.6.6           assertthat_0.2.1      
#[28] gargle_1.2.1           xfun_0.33              hms_1.1.2              evaluate_0.16          promises_1.2.0.1       fansi_1.0.3            progress_1.2.2         dbplyr_2.2.1           igraph_1.3.4          
#[37] DBI_1.1.3              htmlwidgets_1.5.4      spatstat.geom_2.4-0    googledrive_2.0.0      purrr_0.3.4            ellipsis_0.3.2         RSpectra_0.16-1        dplyr_1.0.10           ggpubr_0.4.0          
#[46] backports_1.4.1        gridBase_0.4-7         deldir_1.0-6           vctrs_0.4.1            ggalluvial_0.12.3      Biobase_2.56.0         ROCR_1.0-11            abind_1.4-5            cachem_1.0.6          
#[55] withr_2.5.0            progressr_0.11.0       sctransform_0.3.4      sna_2.7                prettyunits_1.1.1      goftest_1.2-3          svglite_2.1.0          cluster_2.1.4          lazyeval_0.2.2        
#[64] crayon_1.5.1           hdf5r_1.3.5            labeling_0.4.2         pkgconfig_2.0.3        GenomeInfoDb_1.32.4    nlme_3.1-159           rlang_1.0.6            globals_0.16.1         lifecycle_1.0.2       
#[73] miniUI_0.1.1.1         registry_0.5-1         filelock_1.0.2         BiocFileCache_2.4.0    modelr_0.1.10          cellranger_1.1.0       polyclip_1.10-0        matrixStats_0.62.0     lmtest_0.9-40         
#[82] rngtools_1.5.2         Matrix_1.5-1           carData_3.0-5          zoo_1.8-10             reprex_2.0.2           ggridges_0.5.3         GlobalOptions_0.1.2    googlesheets4_1.0.1    png_0.1-7             
#[91] viridisLite_0.4.1      rjson_0.2.21           bitops_1.0-7           KernSmooth_2.23-20     ggnetwork_0.5.10       Biostrings_2.64.1      blob_1.2.3             shape_1.4.6            stringr_1.4.1         
#[100] parallelly_1.32.1      spatstat.random_2.2-0  readr_2.1.2            rstatix_0.7.0          S4Vectors_0.34.0       ggsignif_0.6.3         scales_1.2.1           memoise_2.0.1          magrittr_2.0.3        
#[109] plyr_1.8.7             ica_1.0-3              zlibbioc_1.42.0        compiler_4.2.2         RColorBrewer_1.1-3     clue_0.3-61            fitdistrplus_1.1-8     cli_3.4.0              XVector_0.36.0        
#[118] listenv_0.8.0          patchwork_1.1.2        pbapply_1.5-0          MASS_7.3-58.1          mgcv_1.8-40            tidyselect_1.1.2       stringi_1.7.8          forcats_0.5.2          yaml_2.3.5            
#[127] ggrepel_0.9.1          grid_4.2.2             tools_4.2.2            timechange_0.2.0       future.apply_1.9.1     parallel_4.2.2         circlize_0.4.15        rstudioapi_0.14        foreach_1.5.2         
#[136] gridExtra_2.3          farver_2.1.1           Rtsne_0.16             digest_0.6.29          rgeos_0.5-9            FNN_1.1.3.1            shiny_1.7.2            Rcpp_1.0.9             car_3.1-0             
#[145] broom_1.0.1            later_1.3.0            RcppAnnoy_0.0.19       httr_1.4.4             AnnotationDbi_1.58.0   ComplexHeatmap_2.12.1  colorspace_2.0-3       rvest_1.0.3            XML_3.99-0.10         
#[154] fs_1.5.2               tensor_1.5             reticulate_1.26        IRanges_2.30.1         splines_4.2.2          uwot_0.1.14            spatstat.utils_2.3-1   plotly_4.10.0          systemfonts_1.0.4     
#[163] xtable_1.8-4           jsonlite_1.8.0         R6_2.5.1               pillar_1.8.1           htmltools_0.5.3        mime_0.12              NMF_0.24.0             glue_1.6.2             fastmap_1.1.0         
#[172] BiocParallel_1.30.3    BiocNeighbors_1.14.0   codetools_0.2-18       utf8_1.2.2             lattice_0.20-45        spatstat.sparse_2.1-1  tibble_3.1.8           network_1.17.2         curl_4.3.2            
#[181] leiden_0.4.3           survival_3.4-0         rmarkdown_2.16         statnet.common_4.7.0   munsell_0.5.0          GetoptLong_1.0.5       GenomeInfoDbData_1.2.8 iterators_1.0.14       haven_2.5.1           
#[190] reshape2_1.4.4         gtable_0.3.1           spatstat.core_2.4-4   