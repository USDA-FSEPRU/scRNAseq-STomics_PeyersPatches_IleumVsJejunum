library(Seurat)
library(SeuratDisk)
library(ggplot2)

#  refer to https://kstreet13.github.io/bioc2020trajectories/articles/workshopTrajectories.html#trajectory-inference-1

seu <- LoadH5Seurat('/home/Jayne.Wiarda/SI_PP_SC_ST/PseudotimeTrajectory/B_traj.h5seurat')

## Statistical comparisons of trajectory distributions

#Perform Kolmogorov-Smirnov test to determine if density of pseudotime distribution is different between cells from ileum (I) versus jejunum (J):
  
dat <- data.frame(seu$tissue, seu$pseudotime1, seu$pseudotime2)
jej <- subset(dat, seu.tissue == 'J')
il <- subset(dat, seu.tissue == 'I')

#Trajectory 1:
ks.test(jej$seu.pseudotime1, il$seu.pseudotime1) # Kolmogorov-Smirnov Test assesses if groups have same vs significantly different pseudotime distributions
#Asymptotic two-sample Kolmogorov-Smirnov test
#data:  jej$seu.pseudotime2 and il$seu.pseudotime2
#D = 0.096727, p-value = 1.742e-10
#alternative hypothesis: two-sided

#Trajectory 2:
ks.test(jej$seu.pseudotime2, il$seu.pseudotime2) # Kolmogorov-Smirnov Test assesses if groups have same vs significantly different pseudotime distributions
#Asymptotic two-sample Kolmogorov-Smirnov test
#data:  jej$seu.pseudotime1 and il$seu.pseudotime1
#D = 0.11036, p-value < 2.2e-16
#alternative hypothesis: two-sided

## Plot distribution of cells from different tissue types across trajectories:

ggplot(dat, aes(x=seu.pseudotime1, color = seu.tissue, fill = seu.tissue)) + 
  geom_density(alpha=0.1) +
  scale_color_manual(values=c('red', 'blue')) +
  scale_fill_manual(values=c('red', 'blue')) +
  theme_classic() +
  scale_x_continuous(expand=c(0,0)) +
  scale_y_continuous(expand=c(0,0))
ggsave('/home/Jayne.Wiarda/SI_PP_SC_ST/Figures/sc_Btraj_DifferentialProgression_0.jpeg')

ggplot(dat, aes(x=seu.pseudotime2, color = seu.tissue, fill = seu.tissue)) + 
  geom_density(alpha=0.1) +
  scale_color_manual(values=c('red', 'blue')) +
  scale_fill_manual(values=c('red', 'blue')) +
  theme_classic() +
  scale_x_continuous(expand=c(0,0)) +
  scale_y_continuous(expand=c(0,0))
ggsave('/home/Jayne.Wiarda/SI_PP_SC_ST/Figures/sc_Btraj_DifferentialProgression_1.jpeg')

sessionInfo()
#R version 4.2.2 Patched (2022-11-10 r83330)
#Platform: x86_64-pc-linux-gnu (64-bit)
#Running under: Ubuntu 22.04.1 LTS

#Matrix products: default
#BLAS:   /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3
#LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/libopenblasp-r0.3.20.so

#locale:
#  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C               LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8     LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8    LC_PAPER=en_US.UTF-8      
#[8] LC_NAME=C                  LC_ADDRESS=C               LC_TELEPHONE=C             LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       

#attached base packages:
#  [1] stats4    stats     graphics  grDevices utils     datasets  methods   base     

#other attached packages:
#  [1] ggplot2_3.3.6         SeuratDisk_0.0.0.9020 sp_1.5-0              SeuratObject_4.1.1    Seurat_4.1.1         

#loaded via a namespace (and not attached):
#  [1] scattermore_0.8             coda_0.19-4                 pkgmaker_0.32.2             tidyr_1.2.1                 bit64_4.0.5                 knitr_1.40                  irlba_2.3.5                
#[8] DelayedArray_0.22.0         data.table_1.14.2           rpart_4.1.16                KEGGREST_1.36.3             RCurl_1.98-1.8              doParallel_1.0.17           generics_0.1.3             
#[15] BiocGenerics_0.42.0         ScaledMatrix_1.4.1          cowplot_1.1.1               RSQLite_2.2.17              RANN_2.6.1                  future_1.28.0               bit_4.0.4                  
#[22] phylobase_0.8.10            spatstat.data_2.2-0         xml2_1.3.3                  httpuv_1.6.6                SummarizedExperiment_1.26.1 assertthat_0.2.1            viridis_0.6.2              
#[29] xfun_0.33                   hms_1.1.2                   evaluate_0.16               promises_1.2.0.1            fansi_1.0.3                 progress_1.2.2              igraph_1.3.4               
#[36] DBI_1.1.3                   htmlwidgets_1.5.4           spatstat.geom_2.4-0         purrr_0.3.4                 ellipsis_0.3.2              RSpectra_0.16-1             dplyr_1.0.10               
#[43] ggpubr_0.4.0                backports_1.4.1             annotate_1.74.0             gridBase_0.4-7              locfdr_1.1-8                deldir_1.0-6                MatrixGenerics_1.8.1       
#[50] vctrs_0.4.1                 SingleCellExperiment_1.18.0 ggalluvial_0.12.3           Biobase_2.56.0              ROCR_1.0-11                 abind_1.4-5                 cachem_1.0.6               
#[57] withr_2.5.0                 ggforce_0.3.4               progressr_0.11.0            sctransform_0.3.4           sna_2.7                     prettyunits_1.1.1           goftest_1.2-3              
#[64] softImpute_1.4-1            svglite_2.1.0               cluster_2.1.4               ape_5.6-2                   lazyeval_0.2.2              crayon_1.5.1                genefilter_1.78.0          
#[71] hdf5r_1.3.5                 labeling_0.4.2              edgeR_3.38.4                pkgconfig_2.0.3             tweenr_2.0.2                GenomeInfoDb_1.32.4         nlme_3.1-159               
#[78] vipor_0.4.5                 rlang_1.0.6                 globals_0.16.1              lifecycle_1.0.2             miniUI_0.1.1.1              registry_0.5-1              rsvd_1.0.5                 
#[85] cellranger_1.1.0            polyclip_1.10-0             matrixStats_0.62.0          lmtest_0.9-40               rngtools_1.5.2              Matrix_1.5-1                carData_3.0-5              
#[92] Rhdf5lib_1.18.2             zoo_1.8-10                  beeswarm_0.4.0              ggridges_0.5.3              GlobalOptions_0.1.2         png_0.1-7                   viridisLite_0.4.1          
#[99] rjson_0.2.21                bitops_1.0-7                rhdf5filters_1.8.0          rncl_0.8.6                  KernSmooth_2.23-20          ggnetwork_0.5.10            Biostrings_2.64.1          
#[106] blob_1.2.3                  shape_1.4.6                 stringr_1.4.1               zinbwave_1.18.0             parallelly_1.32.1           spatstat.random_2.2-0       rstatix_0.7.0              
#[113] S4Vectors_0.34.0            ggsignif_0.6.3              beachmat_2.12.0             scales_1.2.1                memoise_2.0.1               magrittr_2.0.3              plyr_1.8.7                 
#[120] ica_1.0-3                   howmany_0.3-1               zlibbioc_1.42.0             compiler_4.2.2              RColorBrewer_1.1-3          clue_0.3-61                 fitdistrplus_1.1-8         
#[127] cli_3.4.0                   ade4_1.7-19                 XVector_0.36.0              listenv_0.8.0               patchwork_1.1.2             pbapply_1.5-0               MASS_7.3-58.1              
#[134] mgcv_1.8-40                 tidyselect_1.1.2            stringi_1.7.8               yaml_2.3.5                  BiocSingular_1.12.0         locfit_1.5-9.6              ggrepel_0.9.1              
#[141] grid_4.2.2                  tools_4.2.2                 future.apply_1.9.1          parallel_4.2.2              circlize_0.4.15             rstudioapi_0.14             uuid_1.1-0                 
#[148] foreach_1.5.2               RNeXML_2.4.7                gridExtra_2.3               farver_2.1.1                Rtsne_0.16                  ggraph_2.0.6                digest_0.6.29              
#[155] rgeos_0.5-9                 FNN_1.1.3.1                 shiny_1.7.2                 Rcpp_1.0.9                  GenomicRanges_1.48.0        car_3.1-0                   broom_1.0.1                
#[162] later_1.3.0                 RcppAnnoy_0.0.19            AnnotationDbi_1.58.0        httr_1.4.4                  ComplexHeatmap_2.12.1       kernlab_0.9-31              colorspace_2.0-3           
#[169] XML_3.99-0.10               tensor_1.5                  reticulate_1.26             clusterExperiment_2.16.0    IRanges_2.30.1              splines_4.2.2               uwot_0.1.14                
#[176] spatstat.utils_2.3-1        graphlayouts_0.8.1          plotly_4.10.0               systemfonts_1.0.4           xtable_1.8-4                jsonlite_1.8.0              tidygraph_1.2.2            
#[183] R6_2.5.1                    pillar_1.8.1                htmltools_0.5.3             mime_0.12                   NMF_0.24.0                  glue_1.6.2                  fastmap_1.1.0              
#[190] BiocParallel_1.30.3         BiocNeighbors_1.14.0        codetools_0.2-18            utf8_1.2.2                  lattice_0.20-45             spatstat.sparse_2.1-1       tibble_3.1.8               
#[197] network_1.17.2              ggbeeswarm_0.6.0            leiden_0.4.3                gtools_3.9.3                survival_3.4-0              limma_3.52.3                rmarkdown_2.16             
#[204] statnet.common_4.7.0        munsell_0.5.0               rhdf5_2.40.0                GetoptLong_1.0.5            GenomeInfoDbData_1.2.8      iterators_1.0.14            HDF5Array_1.24.2           
#[211] reshape2_1.4.4              gtable_0.3.1                spatstat.core_2.4-4   