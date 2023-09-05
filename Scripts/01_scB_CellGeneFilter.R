library(DropletUtils) 
library(Seurat)
library(SeuratObject)
library(readxl)
library(ggplot2)
library(scales)

# Specify data directory and check that all files are there:
data_dir <- c(I2 = "/home/Jayne.Wiarda/SI_PP_SC_ST/SoupX_SConly/I2strainedCounts",
              I3 = "/home/Jayne.Wiarda/SI_PP_SC_ST/SoupX_SConly/I3strainedCounts",
              I4 = "/home/Jayne.Wiarda/SI_PP_SC_ST/SoupX_SConly/I4strainedCounts",
              J2 = "/home/Jayne.Wiarda/SI_PP_SC_ST/SoupX_SConly/J2strainedCounts",
              J3 = "/home/Jayne.Wiarda/SI_PP_SC_ST/SoupX_SConly/J3strainedCounts",
              J4 = "/home/Jayne.Wiarda/SI_PP_SC_ST/SoupX_SConly/J4strainedCounts")
lapply(data_dir, dir) # Should show barcodes.tsv.gz, genes.tsv.gz, and matrix.mtx.gz for each sample listed

# Create a Seurat object with all six samples included. This has a built-in step to filter out genes with sum zero expression across cells of all samples, so this also qualifies as our gene filtering step, leaving us with only non-zero count genes.
scRNA_data <- Read10X(data.dir = data_dir) # read the 10X data from all samples into a data matrix
seu = CreateSeuratObject(counts = scRNA_data, 
                         min.cells = 1) # include only genes expressed in at least one cell; THIS SUFFICES AS THE GENE FILTERING STEP FOR OUR DATA
rm(data_dir, scRNA_data) # clear space

# See how many cells/genes in Seurat object:
seu

#See how many cells in each sample:
table(seu$orig.ident)

#We need to calculate the percentage of mitochondrial genes expressed within each cell.

##To start, we need to form a list of mitochondrial genes to look at:
annotKey <- read_excel('/home/Jayne.Wiarda/SI_PP_SC_ST/GeneAnnotationFiles/UpdatedGeneNameListForSus97GTF_06302021_JEW_SKS.xlsx') # read in our key for modified gene annotations
annot <- "/home/Jayne.Wiarda/SI_PP_SC_ST/GeneAnnotationFiles/Sus_scrofa.Sscrofa11.1.97_modified06302021_JEW_SKS.gtf" # specify file path to Sus scrofa 11.1 version 97 annotation file with modified gene names
mitoGenes <-  system2("grep", args = c('^MT', annot, "| grep -o 'ENSSSCG[0-9]*' | uniq"), stdout = TRUE) # extract mitochondrial gene Ensembl IDs from annotation file
mitoGenes <- annotKey[annotKey$ENSID %in% mitoGenes,] # identify rows for mitochondrial genes in the annotation key
mitoGenes <- mitoGenes$FinalList
length(mitoGenes) # make sure the length is 37, the same numebr of mitochondrial genes found in pigs as in humans
rm(annot, annotKey) # clear space

# Now we need to calculate the percentage of mitochondrial read counts in each cell:
counts <- GetAssayData(object = seu, slot = "counts")
mitoCounts <- counts[rownames(counts) %in% mitoGenes,] # identify rows for mitochondrial genes in the annotation key
pctMito <- ((colSums(mitoCounts))/(colSums(counts)))*100
seu <- AddMetaData(seu, pctMito, col.name = "percent_mito")
rm(counts, mitoCounts,pctMito,mitoGenes)

### Plot QC metrics
#Violin plots:
VlnPlot(seu, 
          features = c("nFeature_RNA", "nCount_RNA", "percent_mito"), # QC metric metadata to plot
          split.by = 'orig.ident', # plot each sample separately
          pt.size = 0, # don't include individual cells as points on plots since there are so many
          ncol = 3)

#Number of genes vs. percent mitochondrial reads:
meta <- seu@meta.data
ggplot(meta, aes(x=percent_mito, y=nFeature_RNA, color = orig.ident))+
  geom_point() + 
  facet_wrap(~orig.ident, nrow =1)+
  theme_get() + 
  ylab("#Genes detected per cell")

# Number of genes vs. number reads:
ggplot(meta, aes(x=nCount_RNA, y=nFeature_RNA, color = percent_mito))+
  geom_point() + 
  facet_wrap(~orig.ident, nrow =1)+
  theme_get() + 
  xlim(0,20000) +
  ylim(0,4000) +
  scale_colour_gradient(low = "gold", high = "red", limits=c(0, 20), oob=squish) +
  ylab("#Genes detected per cell")

#Histograms with our desired thresholds shown from looking at all plots:
ggplot(meta, aes(x=percent_mito,y=..density..)) + 
  geom_histogram(fill="white",color="black",bins=500) + 
  scale_x_continuous(breaks = seq(0, 50, 5), lim = c(0, 50)) + 
  facet_wrap(~orig.ident) +
  geom_vline(aes(xintercept=15),color="red",lty="longdash") + # move this cutoff line where you see fit
  RotatedAxis() + 
  ggtitle('percent_mito')

ggplot(meta, aes(x=nFeature_RNA,y=..density..)) + 
  geom_histogram(fill="white",color="black",bins=500) + 
  scale_x_continuous(breaks = seq(0, 4000, 250), lim = c(0, 4000)) + 
  facet_wrap(~orig.ident) +
  geom_vline(aes(xintercept=700),color="red",lty="longdash") + # move this cutoff line where you see fit
  RotatedAxis() + 
  ggtitle('nFeature_RNA')

ggplot(meta, aes(x=nCount_RNA,y=..density..)) + 
  geom_histogram(fill="white",color="black",bins=500) + 
  scale_x_continuous(breaks = seq(0, 20000, 2000), lim = c(0, 20000)) + 
  facet_wrap(~orig.ident) +
  geom_vline(aes(xintercept=1500),color="red",lty="longdash") + # move this cutoff line where you see fit
  RotatedAxis() + 
  ggtitle('nCount_RNA')

rm(meta) 

#Established QC thresholds are:
##Keep cells with <15% mitochondrial reads
##Keep cells with >700 genes detected
##Keep cells with >1500 total reads

##Filter out poor quality cells
#Identify cells passing each/every QC filter:
keepMito <- WhichCells(seu, expression = percent_mito < 15) # cells passing mito filter
length(keepMito) # how many cells pass this QC filter?

keepGenes <- WhichCells(seu, expression = nFeature_RNA > 700) # cells passing gene filter
length(keepGenes) # how many cells pass this QC filter?

keepUMI <- WhichCells(seu, expression = nCount_RNA > 1500) # cells passing UMI filter
length(keepUMI) # how many cells pass this QC filter?

keep <- Reduce(intersect, list(keepMito, keepGenes, keepUMI)) # cells passing all QC filters
length(keep) # how many cells pass this QC filter?

rm(keepMito, keepGenes, keepUMI) # free up space

# Create new Seurat object with only cells passing all QC filters:
seuKeep <- subset(seu, 
                  cells = keep)
rm(seu) # free up space

# How many cells in each sample now?
table(seuKeep$orig.ident)

### Look over plots with only cells passing QC
# Violin plots:
VlnPlot(seuKeep, 
        features = c("nFeature_RNA", "nCount_RNA", "percent_mito"), # QC metric metadata to plot
        split.by = 'orig.ident', # plot each sample separately
        pt.size = 0, # don't include individual cells as points on plots since there are so many
        ncol = 3)

# Number of genes vs. percent mitochondrial reads:
meta <- seuKeep@meta.data
ggplot(meta, aes(x=percent_mito, y=nFeature_RNA, color = orig.ident))+
  geom_point() + 
  facet_wrap(~orig.ident, nrow =1)+
  theme_get() + 
  ylab("#Genes detected per cell")

# Number of genes vs. number reads:
ggplot(meta, aes(x=nCount_RNA, y=nFeature_RNA, color = percent_mito))+
  geom_point() + 
  facet_wrap(~orig.ident, nrow =1)+
  theme_get() + 
  xlim(0,20000) +
  ylim(0,4000) +
  scale_colour_gradient(low = "gold", high = "red", limits=c(0, 20), oob=squish) +
  ylab("#Genes detected per cell")

rm(meta) # free up space

## Save counts of cells passing QC from each sample
Idents(seuKeep) <- seuKeep$orig.ident
sub <- subset(seuKeep, ident = "I2")
write10xCounts(x = sub@assays$RNA@counts, path = "/home/Jayne.Wiarda/SI_PP_SC_ST/QC_SConly/I2onlyFilteredQC", version = "3", overwrite = TRUE) 
sub <- subset(seuKeep, ident = "I3")
write10xCounts(x = sub@assays$RNA@counts, path = "/home/Jayne.Wiarda/SI_PP_SC_ST/QC_SConly/I3onlyFilteredQC", version = "3", overwrite = TRUE) 
sub <- subset(seuKeep, ident = "I4")
write10xCounts(x = sub@assays$RNA@counts, path = "/home/Jayne.Wiarda/SI_PP_SC_ST/QC_SConly/I4onlyFilteredQC", version = "3", overwrite = TRUE) 
sub <- subset(seuKeep, ident = "J2")
write10xCounts(x = sub@assays$RNA@counts, path = "/home/Jayne.Wiarda/SI_PP_SC_ST/QC_SConly/J2onlyFilteredQC", version = "3", overwrite = TRUE) 
sub <- subset(seuKeep, ident = "J3")
write10xCounts(x = sub@assays$RNA@counts, path = "/home/Jayne.Wiarda/SI_PP_SC_ST/QC_SConly/J3onlyFilteredQC", version = "3", overwrite = TRUE) 
sub <- subset(seuKeep, ident = "J4")
write10xCounts(x = sub@assays$RNA@counts, path = "/home/Jayne.Wiarda/SI_PP_SC_ST/QC_SConly/J4onlyFilteredQC", version = "3", overwrite = TRUE) 

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
##  [1] scales_1.2.0                ggplot2_3.3.6              
##  [3] readxl_1.4.0                sp_1.5-0                   
##  [5] SeuratObject_4.1.0          Seurat_4.1.1               
##  [7] DropletUtils_1.16.0         SingleCellExperiment_1.18.0
##  [9] SummarizedExperiment_1.26.1 Biobase_2.56.0             
## [11] GenomicRanges_1.48.0        GenomeInfoDb_1.32.2        
## [13] IRanges_2.30.0              S4Vectors_0.34.0           
## [15] BiocGenerics_0.42.0         MatrixGenerics_1.8.1       
## [17] matrixStats_0.62.0         
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
##  [21] xfun_0.31                 dplyr_1.0.9              
##  [23] crayon_1.5.1              RCurl_1.98-1.7           
##  [25] jsonlite_1.8.0            progressr_0.10.1         
##  [27] spatstat.data_2.2-0       survival_3.3-1           
##  [29] zoo_1.8-10                glue_1.6.2               
##  [31] polyclip_1.10-0           gtable_0.3.0             
##  [33] zlibbioc_1.42.0           XVector_0.36.0           
##  [35] leiden_0.4.2              DelayedArray_0.22.0      
##  [37] Rhdf5lib_1.18.2           future.apply_1.9.0       
##  [39] HDF5Array_1.24.1          abind_1.4-5              
##  [41] DBI_1.1.3                 edgeR_3.38.1             
##  [43] spatstat.random_2.2-0     miniUI_0.1.1.1           
##  [45] Rcpp_1.0.9                viridisLite_0.4.0        
##  [47] xtable_1.8-4              reticulate_1.25          
##  [49] spatstat.core_2.4-4       dqrng_0.3.0              
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
##  [79] yaml_2.3.5                knitr_1.39               
##  [81] fitdistrplus_1.1-8        purrr_0.3.4              
##  [83] RANN_2.6.1                nlme_3.1-158             
##  [85] pbapply_1.5-0             future_1.26.1            
##  [87] sparseMatrixStats_1.8.0   mime_0.12                
##  [89] ggrastr_1.0.1             R.oo_1.25.0              
##  [91] compiler_4.2.1            rstudioapi_0.13          
##  [93] beeswarm_0.4.0            plotly_4.10.0            
##  [95] png_0.1-7                 spatstat.utils_2.3-1     
##  [97] tibble_3.1.7              stringi_1.7.8            
##  [99] highr_0.9                 rgeos_0.5-9              
## [101] lattice_0.20-45           Matrix_1.4-1             
## [103] vctrs_0.4.1               pillar_1.8.0             
## [105] lifecycle_1.0.1           rhdf5filters_1.8.0       
## [107] spatstat.geom_2.4-0       lmtest_0.9-40            
## [109] RcppAnnoy_0.0.19          data.table_1.14.2        
## [111] cowplot_1.1.1             bitops_1.0-7             
## [113] irlba_2.3.5               httpuv_1.6.5             
## [115] patchwork_1.1.1           R6_2.5.1                 
## [117] promises_1.2.0.1          KernSmooth_2.23-20       
## [119] gridExtra_2.3             vipor_0.4.5              
## [121] parallelly_1.32.0         codetools_0.2-18         
## [123] MASS_7.3-58               assertthat_0.2.1         
## [125] rhdf5_2.40.0              withr_2.5.0              
## [127] sctransform_0.3.3         GenomeInfoDbData_1.2.8   
## [129] mgcv_1.8-40               parallel_4.2.1           
## [131] rpart_4.1.16              grid_4.2.1               
## [133] beachmat_2.12.0           tidyr_1.2.0              
## [135] rmarkdown_2.14            DelayedMatrixStats_1.18.0
## [137] Rtsne_0.16                shiny_1.7.2              
## [139] ggbeeswarm_0.6.0