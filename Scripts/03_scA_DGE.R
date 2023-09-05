### Load required software packages
library(Seurat)
library(writexl)
library(readxl)
library(SeuratDisk)
library(dplyr)
library(ggplot2)
library(scales)

## Read in Seurat object

seu <- LoadH5Seurat('/home/Jayne.Wiarda/SI_PP_SC_ST/AnnotatedSeurat/AllSamples_annotated.h5seurat')
DefaultAssay(seu) <- 'SCT'

# Prep SCT assay for DGE testing:
seu <- PrepSCTFindMarkers(seu)

## Perform DGE testing between all cell types:
Idents(seu) <- seu$celltype
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
write_xlsx(de, '/home/Jayne.Wiarda/SI_PP_SC_ST/DGE/CellTypeDGE_AllSamples.xlsx')

# Save list of background genes used for DGE analyses:
features <- read_excel("/home/Jayne.Wiarda/scRNAseqEpSI/UpdatedGeneNameListForSus97GTF.xlsx") # read in the gene information for the filtered gene list from QC analysis
features <- features[c("ENSID", "FinalAnnot")] # subset only the columns of gene symbols, Ensembl IDs, and the names used for analysis
features$FinalAnnot <-gsub("_", "-", features$FinalAnnot) # replace all underscores with dashes since this occurred when processing data into Seurat
features$FinalAnnot[features$FinalAnnot == "ISG12(A)"] <- "ISG12(A)-ENSSSCG00000035297" # manually fix one mislabeled gene from file (double checked all other genes are fine as-is)
features <- features[features$FinalAnnot %in% rownames(seu),]
write_xlsx(features, '/home/Jayne.Wiarda/SI_PP_SC_ST/DGE/CellTypeDGE_AllSamples_BackgroundGeneList.xlsx')

# Plot top DE genes:
topgenes <- de %>% group_by(cluster) %>% top_n(10, avg_log2FC) # only plot top 10 genes per cluster, as determined by highest average logFC values
DoHeatmap(subset(seu, downsample = 100), # take only 100 cells per cluster for plotting
          features = as.character(topgenes$gene), 
          assay = "SCT") 

## Perform DGE testing between tissues within a cell type:

### Define pairwise comparisons

#Combine tissue and cell identities:
seu$tissue <- substr(seu$orig.ident, 1, 1)   
seu$combo <- paste(seu$celltype, seu$tissue, sep = '_')
table(seu$combo) # identify any anntotation with fewer than 3 cells; we will need to omit this comparison from DGE analysis

#Create pairwise comparisons:
pop1 <- paste(unique(seu$celltype), 'J', sep = '_') # identify all of our cluster IDs
pop2 <- paste(unique(seu$celltype), 'I', sep = '_') # identify all of our cluster IDs
comps <- data.frame(pop1, pop2)
#View(comps) # identify row of cycling gd T cells & BEST4 enterocytes that needs to be removed due to cell numbers being too low in one tissue
comps <- comps[-c(29, 32), ] # remove rows where not enough cells present for comparison

### Perform DGE testing between tissues within a cell type
Idents(seu) <- seu$combo
DefaultAssay(seu) <- 'SCT'
results <- list()
for(i in 1:nrow(comps)) {
  markers <- FindMarkers(seu, 
                         ident.1 = comps[i,1], 
                         ident.2 = comps[i,2],
                         only.pos = FALSE) 
  markers$gene <- rownames(markers)
  markers$pop1 <- paste(comps[i,1])
  markers$pop2 <- paste(comps[i,2])
  markers$comparison <- paste(markers$pop1, markers$pop2, sep = 'v')
  results[[i]] <- markers
} # if any of the comparisons don't turn up DE genes, this function won't work...
pwAll <- do.call(rbind, results)
pwAll <- subset(pwAll, p_val_adj < 0.05)

#Add Ensembl IDs:
features <- read_excel("/home/Jayne.Wiarda/scRNAseqEpSI/UpdatedGeneNameListForSus97GTF.xlsx") # read in the gene information for the filtered gene list from QC analysis
features <- features[c("ENSID", "FinalAnnot")] # subset only the columns of gene symbols, Ensembl IDs, and the names used for analysis
features$FinalAnnot <-gsub("_", "-", features$FinalAnnot) # replace all underscores with dashes since this occurred when processing data into Seurat
features$FinalAnnot[features$FinalAnnot == "ISG12(A)"] <- "ISG12(A)-ENSSSCG00000035297" # manually fix one mislabeled gene from file (double checked all other genes are fine as-is)
pwAll <- merge(pwAll,
               features, 
               by.x = "gene", 
               by.y = "FinalAnnot") # merge the DE gene lists with the additional gene information
pwAll <- pwAll[order(pwAll$comparison, pwAll$p_val_adj),] # reorder by lowest to highest p-value within each cluster

#Save DEG list:
write_xlsx(pwAll, '/home/Jayne.Wiarda/SI_PP_SC_ST/DGE/CellTypeDGE_AllSamples_TissueDEGs.xlsx')

# Show number of DEGs between tissues per cell type:
results <- list()
for(i in 1:nrow(comps)) {
  DE <- sum(pwAll$pop1 == comps[i,1] & pwAll$pop2 == comps[i,2])
  DE$pop1 <- paste(comps[i,1])
  DE$pop2 <- paste(comps[i,2])
  results[[i]] <- DE
}

pwDE <- data.frame(matrix(unlist(results), nrow=length(results), byrow=T))
colnames(pwDE) <- c('DEgenes', 'pop1', 'pop2')
pwDE$DEgenes <- as.numeric(pwDE$DEgenes)
id <- unique(seu$celltype)
id <- id[c(1:28, 30:31)]
pwDE$cellID <- id

ggplot(pwDE, aes(x = 1, y = cellID, fill = DEgenes)) +
  geom_tile(color = 'black')+
  geom_text(aes(label = DEgenes)) +
  scale_fill_gradientn(colours = c('white', 'wheat', 'goldenrod1', 'red'),
                       limits = c(0, 30), oob=squish)+ 
  theme_classic() +
  scale_y_discrete(limits=rev)

# Plot selected genes in B cells:
DefaultAssay(seu) <- 'SCT'
Idents(seu) <- seu$celltype
seu <- subset(seu, idents = c('Activated B cells', 'Cycling B cells', 'Resting B cells', 'Transitioning B cells'))
seu$tissue <- substr(seu$orig.ident, 1, 1)   
seu$combo <- paste(seu$celltype, seu$tissue, sep = '_')
seu$combo <- factor(seu$combo,levels=c('Activated B cells_J', 'Activated B cells_I',
                                       'Cycling B cells_J', 'Cycling B cells_I',
                                       'Resting B cells_J', 'Resting B cells_I', 
                                       'Transitioning B cells_J', 'Transitioning B cells_I'))

Idents(seu) <- seu$combo
VlnPlot(seu, features = c('ENSSSCG00000010077', 'ENSSSCG00000037775', 'ZP3', 
                          'CD74', 'MS4A1', 'SLA-2', 'SLA-DMA', 'SLA-DMB', 'SLA-DQA1', 
                          'SLA-DQB1', 'SLA-DRA', 'SLA-DRB1'), 
        cols = c('blue', 'red', 'blue', 'red', 'blue', 'red', 'blue', 'red'), 
        pt.size = 0.000) & 
  geom_boxplot(width = 0.2, fill = 'white', col = 'black', lwd = 0.5, outlier.shape = NA, coef = 0) &
  stat_summary(fun.y = mean, geom='point', size = 2, colour = "red4") # red dot shows mean; box shows IQR

### View session information
sessionInfo()
#R version 4.2.2 Patched (2022-11-10 r83330)
#Platform: x86_64-pc-linux-gnu (64-bit)
#Running under: Ubuntu 22.04.2 LTS

#Matrix products: default
#BLAS:   /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3
#LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/libopenblasp-r0.3.20.so

#locale:
#  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C               LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8     LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
#[7] LC_PAPER=en_US.UTF-8       LC_NAME=C                  LC_ADDRESS=C               LC_TELEPHONE=C             LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       

#attached base packages:
#  [1] grid      stats4    stats     graphics  grDevices utils     datasets  methods   base     

#other attached packages:
#[1] scales_1.2.1          ggplot2_3.3.6         dplyr_1.0.10          SeuratDisk_0.0.0.9020 readxl_1.4.1          writexl_1.4.0         sp_1.5-0             
#[8] SeuratObject_4.1.1    Seurat_4.1.1         

#loaded via a namespace (and not attached):
#  [1] scattermore_0.8             coda_0.19-4                 pkgmaker_0.32.2             tidyr_1.2.1                 knitr_1.40                 
#[6] bit64_4.0.5                 irlba_2.3.5                 DelayedArray_0.22.0         data.table_1.14.2           rpart_4.1.16               
#[11] KEGGREST_1.36.3             RCurl_1.98-1.8              doParallel_1.0.17           generics_0.1.3              BiocGenerics_0.42.0        
#[16] ScaledMatrix_1.4.1          cowplot_1.1.1               RSQLite_2.2.17              RANN_2.6.1                  future_1.28.0              
#[21] bit_4.0.4                   phylobase_0.8.10            spatstat.data_2.2-0         xml2_1.3.3                  httpuv_1.6.6               
#[26] SummarizedExperiment_1.26.1 assertthat_0.2.1            xfun_0.33                   hms_1.1.2                   evaluate_0.16              
#[31] promises_1.2.0.1            fansi_1.0.3                 progress_1.2.2              igraph_1.3.4                DBI_1.1.3                  
#[36] htmlwidgets_1.5.4           spatstat.geom_2.4-0         purrr_0.3.4                 ellipsis_0.3.2              RSpectra_0.16-1            
#[41] ggpubr_0.4.0                backports_1.4.1             annotate_1.74.0             gridBase_0.4-7              locfdr_1.1-8               
#[46] deldir_1.0-6                MatrixGenerics_1.8.1        vctrs_0.4.1                 SingleCellExperiment_1.18.0 ggalluvial_0.12.3          
#[51] Biobase_2.56.0              ROCR_1.0-11                 abind_1.4-5                 cachem_1.0.6                withr_2.5.0                
#[56] progressr_0.11.0            sctransform_0.3.4           sna_2.7                     prettyunits_1.1.1           goftest_1.2-3              
#[61] softImpute_1.4-1            svglite_2.1.0               cluster_2.1.4               ape_5.6-2                   lazyeval_0.2.2             
#[66] crayon_1.5.1                genefilter_1.78.0           hdf5r_1.3.5                 labeling_0.4.2              edgeR_3.38.4               
#[71] pkgconfig_2.0.3             GenomeInfoDb_1.32.4         vipor_0.4.5                 nlme_3.1-159                rlang_1.0.6                
#[76] globals_0.16.1              lifecycle_1.0.2             miniUI_0.1.1.1              registry_0.5-1              rsvd_1.0.5                 
#[81] ggrastr_1.0.1               cellranger_1.1.0            polyclip_1.10-0             matrixStats_0.62.0          lmtest_0.9-40              
#[86] rngtools_1.5.2              Matrix_1.5-1                carData_3.0-5               Rhdf5lib_1.18.2             zoo_1.8-10                 
#[91] beeswarm_0.4.0              ggridges_0.5.3              GlobalOptions_0.1.2         png_0.1-7                   viridisLite_0.4.1          
#[96] rjson_0.2.21                bitops_1.0-7                rhdf5filters_1.8.0          rncl_0.8.6                  KernSmooth_2.23-20         
#[101] ggnetwork_0.5.10            Biostrings_2.64.1           blob_1.2.3                  shape_1.4.6                 stringr_1.4.1              
#[106] zinbwave_1.18.0             parallelly_1.32.1           spatstat.random_2.2-0       rstatix_0.7.0               S4Vectors_0.34.0           
#[111] ggsignif_0.6.3              beachmat_2.12.0             memoise_2.0.1               magrittr_2.0.3              plyr_1.8.7                 
#[116] ica_1.0-3                   howmany_0.3-1               zlibbioc_1.42.0             compiler_4.2.2              RColorBrewer_1.1-3         
#[121] clue_0.3-61                 fitdistrplus_1.1-8          cli_3.4.0                   ade4_1.7-19                 XVector_0.36.0             
#[126] listenv_0.8.0               patchwork_1.1.2             pbapply_1.5-0               MASS_7.3-58.1               mgcv_1.8-40                
#[131] tidyselect_1.1.2            stringi_1.7.8               yaml_2.3.5                  BiocSingular_1.12.0         locfit_1.5-9.6             
#[136] ggrepel_0.9.1               tools_4.2.2                 future.apply_1.9.1          parallel_4.2.2              circlize_0.4.15            
#[141] rstudioapi_0.14             uuid_1.1-0                  foreach_1.5.2               RNeXML_2.4.7                gridExtra_2.3              
#[146] farver_2.1.1                Rtsne_0.16                  digest_0.6.29               rgeos_0.5-9                 FNN_1.1.3.1                
#[151] shiny_1.7.2                 Rcpp_1.0.9                  GenomicRanges_1.48.0        car_3.1-0                   broom_1.0.1                
#[156] later_1.3.0                 RcppAnnoy_0.0.19            features AnnotationDbi_1.58.0        ComplexHeatmap_2.12.1      
#[161] kernlab_0.9-31              colorspace_2.0-3            XML_3.99-0.10               tensor_1.5                  reticulate_1.26            
#[166] clusterExperiment_2.16.0    IRanges_2.30.1              splines_4.2.2               uwot_0.1.14                 spatstat.utils_2.3-1       
#[171] plotly_4.10.0               systemfonts_1.0.4           xtable_1.8-4                jsonlite_1.8.0              R6_2.5.1                   
#[176] pillar_1.8.1                htmltools_0.5.3             mime_0.12                   NMF_0.24.0                  glue_1.6.2                 
#[181] fastmap_1.1.0               BiocParallel_1.30.3         BiocNeighbors_1.14.0        codetools_0.2-18            utf8_1.2.2                 
#[186] lattice_0.20-45             spatstat.sparse_2.1-1       tibble_3.1.8                network_1.17.2              ggbeeswarm_0.6.0           
#[191] leiden_0.4.3                survival_3.4-0              limma_3.52.3                rmarkdown_2.16              statnet.common_4.7.0       
#[196] munsell_0.5.0               rhdf5_2.40.0                GetoptLong_1.0.5            GenomeInfoDbData_1.2.8      iterators_1.0.14           
#[201] HDF5Array_1.24.2            reshape2_1.4.4              gtable_0.3.1                spatstat.core_2.4-4     