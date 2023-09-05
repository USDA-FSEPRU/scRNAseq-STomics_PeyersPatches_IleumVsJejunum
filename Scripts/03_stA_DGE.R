### Load required software packages
library(Seurat)
library(writexl)
library(readxl)

## Read in Seurat object

seu <- readRDS('/home/Jayne.Wiarda/SI_PP_SC_ST/AnnotatedSeurat/AllSamples_annotated_ST.rds')
DefaultAssay(seu) <- 'SCT'

# Prep SCT assay for DGE testing:
seu <- PrepSCTFindMarkers(seu)

## Perform DGE testing between regions for all samples together, manual clustering annotation:
Idents(seu) <- seu$Region_ManAnn
de <- FindAllMarkers(seu, 
                     assay = 'SCT',
                     only.pos = FALSE)
de <- subset(de, p_val_adj < 0.05)

#Add Ensembl IDs:
features <- read_excel("/home/Jayne.Wiarda/scRNAseqEpSI/UpdatedGeneNameListForSus97GTF.xlsx") # read in the gene information for the filtered gene list from QC analysis
features <- features[c("ENSID", "FinalAnnot")] # subset only the columns of gene symbols, Ensembl IDs, and the names used for analysis
features$FinalAnnot <-gsub("_", "-", features$FinalAnnot) # replace all underscores with dashes since this occurred when processing data into Seurat
features$FinalAnnot[features$FinalAnnot == "ISG12(A)"] <- "ISG12(A)-ENSSSCG00000035297" # manually fix one mislabeled gene from file (double checked all other genes are fine as-is)
de <- merge(de,
            features, 
            by.x = "gene", 
            by.y = "FinalAnnot") # merge the DE gene lists with the additional gene information
de <- de[order(de$cluster, de$p_val_adj),] # reorder by lowest to highest p-value within each cluster

#Save DEG list:
write_xlsx(de, '/home/Jayne.Wiarda/SI_PP_SC_ST/DGE/SpatialRegions_DEGs_AllSamples_ManualAnnotation.xlsx')

# Save list of background genes used for DGE analyses:
features <- read_excel("/home/Jayne.Wiarda/scRNAseqEpSI/UpdatedGeneNameListForSus97GTF.xlsx") # read in the gene information for the filtered gene list from QC analysis
features <- features[c("ENSID", "FinalAnnot")] # subset only the columns of gene symbols, Ensembl IDs, and the names used for analysis
features$FinalAnnot <-gsub("_", "-", features$FinalAnnot) # replace all underscores with dashes since this occurred when processing data into Seurat
features$FinalAnnot[features$FinalAnnot == "ISG12(A)"] <- "ISG12(A)-ENSSSCG00000035297" # manually fix one mislabeled gene from file (double checked all other genes are fine as-is)
features <- features[features$FinalAnnot %in% rownames(seu),]
write_xlsx(features, '/home/Jayne.Wiarda/SI_PP_SC_ST/DGE/SpatialRegions_DEGs_BackgroundGeneList_AllSamples_ManualAnnotation.xlsx')

# Plot top DE genes:
topgenes <- de %>% group_by(cluster) %>% top_n(20, avg_log2FC) # only plot top 20 genes per cluster, as determined by highest average logFC values
DoHeatmap(subset(seu, downsample = 100), # take only 100 cells per cluster for plotting
          features = as.character(topgenes$gene), 
          assay = "SCT") 

## Perform DGE testing between regions for all samples together, using clustering annotation:
Idents(seu) <- seu$Region_Clust
de <- FindAllMarkers(seu, 
                     assay = 'SCT',
                     only.pos = FALSE)
de <- subset(de, p_val_adj < 0.05)

#Add Ensembl IDs:
features <- read_excel("/home/Jayne.Wiarda/scRNAseqEpSI/UpdatedGeneNameListForSus97GTF.xlsx") # read in the gene information for the filtered gene list from QC analysis
features <- features[c("ENSID", "FinalAnnot")] # subset only the columns of gene symbols, Ensembl IDs, and the names used for analysis
features$FinalAnnot <-gsub("_", "-", features$FinalAnnot) # replace all underscores with dashes since this occurred when processing data into Seurat
features$FinalAnnot[features$FinalAnnot == "ISG12(A)"] <- "ISG12(A)-ENSSSCG00000035297" # manually fix one mislabeled gene from file (double checked all other genes are fine as-is)
de <- merge(de,
            features, 
            by.x = "gene", 
            by.y = "FinalAnnot") # merge the DE gene lists with the additional gene information
de <- de[order(de$cluster, de$p_val_adj),] # reorder by lowest to highest p-value within each cluster

#Save DEG list:
write_xlsx(de, '/home/Jayne.Wiarda/SI_PP_SC_ST/DGE/SpatialRegions_DEGs_AllSamples_ClusteringAnnotation.xlsx')

# Save list of background genes used for DGE analyses:
features <- read_excel("/home/Jayne.Wiarda/scRNAseqEpSI/UpdatedGeneNameListForSus97GTF.xlsx") # read in the gene information for the filtered gene list from QC analysis
features <- features[c("ENSID", "FinalAnnot")] # subset only the columns of gene symbols, Ensembl IDs, and the names used for analysis
features$FinalAnnot <-gsub("_", "-", features$FinalAnnot) # replace all underscores with dashes since this occurred when processing data into Seurat
features$FinalAnnot[features$FinalAnnot == "ISG12(A)"] <- "ISG12(A)-ENSSSCG00000035297" # manually fix one mislabeled gene from file (double checked all other genes are fine as-is)
features <- features[features$FinalAnnot %in% rownames(seu),]
write_xlsx(features, '/home/Jayne.Wiarda/SI_PP_SC_ST/DGE/SpatialRegions_DEGs_BackgroundGeneList_AllSamples_ClusteringAnnotation.xlsx')

# Plot top DE genes:
topgenes <- de %>% group_by(cluster) %>% top_n(20, avg_log2FC) # only plot top 20 genes per cluster, as determined by highest average logFC values
DoHeatmap(subset(seu, downsample = 100), # take only 100 cells per cluster for plotting
          features = as.character(topgenes$gene), 
          assay = "SCT") 

### View session information
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
#  [1] grid      stats4    stats     graphics  grDevices utils     datasets  methods   base     

#other attached packages:
#  [1] readxl_1.4.1       writexl_1.4.0      sp_1.5-0           SeuratObject_4.1.1 Seurat_4.1.1      

#loaded via a namespace (and not attached):
#  [1] utf8_1.2.2             reticulate_1.26        tidyselect_1.1.2       RSQLite_2.2.17         AnnotationDbi_1.58.0   htmlwidgets_1.5.4      Rtsne_0.16             munsell_0.5.0          codetools_0.2-18      
#[10] ica_1.0-3              future_1.28.0          miniUI_0.1.1.1         withr_2.5.0            spatstat.random_2.2-0  colorspace_2.0-3       progressr_0.11.0       Biobase_2.56.0         filelock_1.0.2        
#[19] knitr_1.40             rstudioapi_0.14        ROCR_1.0-11            tensor_1.5             listenv_0.8.0          labeling_0.4.2         GenomeInfoDbData_1.2.8 polyclip_1.10-0        farver_2.1.1          
#[28] bit64_4.0.5            parallelly_1.32.1      vctrs_0.4.1            generics_0.1.3         xfun_0.33              timechange_0.2.0       BiocFileCache_2.4.0    R6_2.5.1               GenomeInfoDb_1.32.4   
#[37] hdf5r_1.3.5            bitops_1.0-7           spatstat.utils_2.3-1   cachem_1.0.6           assertthat_0.2.1       vroom_1.5.7            promises_1.2.0.1       scales_1.2.1           googlesheets4_1.0.1   
#[46] rgeos_0.5-9            gtable_0.3.1           globals_0.16.1         goftest_1.2-3          rlang_1.0.6            splines_4.2.2          lazyeval_0.2.2         gargle_1.2.1           spatstat.geom_2.4-0   
#[55] broom_1.0.1            yaml_2.3.5             reshape2_1.4.4         abind_1.4-5            modelr_0.1.10          backports_1.4.1        httpuv_1.6.6           tools_4.2.2            ggplot2_3.3.6         
#[64] ellipsis_0.3.2         spatstat.core_2.4-4    RColorBrewer_1.1-3     BiocGenerics_0.42.0    ggridges_0.5.3         Rcpp_1.0.9             plyr_1.8.7             progress_1.2.2         zlibbioc_1.42.0       
#[73] purrr_0.3.4            RCurl_1.98-1.8         prettyunits_1.1.1      rpart_4.1.16           deldir_1.0-6           pbapply_1.5-0          cowplot_1.1.1          S4Vectors_0.34.0       zoo_1.8-10            
#[82] haven_2.5.1            ggrepel_0.9.1          cluster_2.1.4          fs_1.5.2               magrittr_2.0.3         data.table_1.14.2      scattermore_0.8        lmtest_0.9-40          reprex_2.0.2          
#[91] RANN_2.6.1             googledrive_2.0.0      fitdistrplus_1.1-8     matrixStats_0.62.0     hms_1.1.2              patchwork_1.1.2        mime_0.12              evaluate_0.16          xtable_1.8-4          
#[100] XML_3.99-0.10          IRanges_2.30.1         gridExtra_2.3          compiler_4.2.2         tibble_3.1.8           KernSmooth_2.23-20     crayon_1.5.1           htmltools_0.5.3        mgcv_1.8-40           
#[109] later_1.3.0            tzdb_0.3.0             tidyr_1.2.1            lubridate_1.9.0        DBI_1.1.3              dbplyr_2.2.1           MASS_7.3-58.1          rappdirs_0.3.3         Matrix_1.5-1          
#[118] readr_2.1.2            cli_3.4.0              parallel_4.2.2         igraph_1.3.4           forcats_0.5.2          pkgconfig_2.0.3        plotly_4.10.0          spatstat.sparse_2.1-1  xml2_1.3.3            
#[127] XVector_0.36.0         rvest_1.0.3            stringr_1.4.1          digest_0.6.29          sctransform_0.3.4      RcppAnnoy_0.0.19       spatstat.data_2.2-0    Biostrings_2.64.1      rmarkdown_2.16        
#[136] cellranger_1.1.0       leiden_0.4.3           uwot_0.1.14            curl_4.3.2             shiny_1.7.2            lifecycle_1.0.2        nlme_3.1-159           jsonlite_1.8.0         viridisLite_0.4.1     
#[145] fansi_1.0.3            pillar_1.8.1           lattice_0.20-45        KEGGREST_1.36.3        fastmap_1.1.0          httr_1.4.4             survival_3.4-0         glue_1.6.2             png_0.1-7             
#[154] bit_4.0.4              stringi_1.7.8          blob_1.2.3             memoise_2.0.1          dplyr_1.0.10           irlba_2.3.5            future.apply_1.9.1   