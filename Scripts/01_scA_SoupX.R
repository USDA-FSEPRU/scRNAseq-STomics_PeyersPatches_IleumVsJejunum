library(DropletUtils) 
library(ggplot2) 
library(SoupX) 

## Perform ambient RNA estimation and removal on each individual sample 

### Sample I2
#### Provide file path to Cell Ranger outputs:
sc = load10X('/home/Jayne.Wiarda/SI_PP_SC_ST/CellRangerOutputs/I2') 

#### Check the meta data to make sure we also have cluster information and t-SNE coordinates:
head(sc$metaData, n = 3)

#### Plot data according to Cell Ranger cluster assignments and t-SNE coordinates:
dd = sc$metaData # create an object with all the metadata
mids = aggregate(cbind(tSNE1,tSNE2) ~ clusters,data=dd,FUN=mean) # determine t-SNE coordinates for middle of each cluster
gg = ggplot(dd,aes(tSNE1,tSNE2)) + # make a t-SNE plot
  geom_point(aes(colour=factor(clusters)),size=0.2) +
  geom_label(data=mids,aes(label=clusters)) 
plot(gg) # show plot

#### Check expression patterns for some canonical genes:
dd$CD3E = sc$toc["CD3E", ] # make column of gene expression values for CD3E (T cell gene)
dd$IgLambdaV = sc$toc["ENSSSCG00000038719", ] # make column of gene expression values for gene that codes for Ig lambda V region (B cell gene)
dd$CD79B = sc$toc["CD79B", ] # make column of gene expression values for CD79B (B cell gene)
dd$FABP6 = sc$toc["FABP6", ] # make column of gene expression values for gene that codes for FABP6 (epithelial cell gene)
dd$EPCAM = sc$toc["EPCAM", ] # make column of gene expression values for gene that codes for EPCAM (epithelial cell gene)
dd$GNLY = sc$toc["GNLY", ] # make column of gene expression values for gene that codes for GNLY (cytotoxicty gene)
dd$HBB = sc$toc["HBB", ] # make column of gene expression values for gene that codes for HBB (erythrocyte gene)

a1 <- ggplot(dd, aes(tSNE1,tSNE2)) + geom_point(aes(colour = CD3E > 0)) + ggtitle('CD3E') # which cells express this gene?
a2 <- ggplot(dd, aes(tSNE1,tSNE2)) + geom_point(aes(colour = IgLambdaV > 0)) + ggtitle('IgLambdaV')
a3 <- ggplot(dd, aes(tSNE1,tSNE2)) + geom_point(aes(colour = CD79B > 0)) + ggtitle('CD79B')
a4 <- ggplot(dd, aes(tSNE1,tSNE2)) + geom_point(aes(colour = FABP6 > 0)) + ggtitle('FABP6')
a5 <- ggplot(dd, aes(tSNE1,tSNE2)) + geom_point(aes(colour = EPCAM > 0)) + ggtitle('EPCAM')
a6 <- ggplot(dd, aes(tSNE1,tSNE2)) + geom_point(aes(colour = GNLY > 0)) + ggtitle('GNLY')
a7 <- ggplot(dd, aes(tSNE1,tSNE2)) + geom_point(aes(colour = HBB > 0)) + ggtitle('HBB')
(a1+a2+a3) / (a4+a5+a6) # show 6 plots

a1 <- plotMarkerMap(sc, "CD3E") + ggtitle('CD3E')  # if we assumed all cells were nothing but soup, which cells still show higher than expected expression for the gene (TRUE = expression levels higher than expected if cell was just soup, so likely real expression). This just gives us an idea of soup expression, this is NOT a formal analysis used for removing the soup RNA.
a2 <- plotMarkerMap(sc, "ENSSSCG00000038719") + ggtitle('IgLambdaV')
a3 <- plotMarkerMap(sc, "CD79B") + ggtitle('CD79B')
a4 <- plotMarkerMap(sc, "FABP6") + ggtitle('FABP6')
a5 <- plotMarkerMap(sc, "EPCAM") + ggtitle('EPCAM')
a6 <- plotMarkerMap(sc, "GNLY") + ggtitle('GNLY')
a7 <- plotMarkerMap(sc, "HBB") + ggtitle('HBB')
(a1+a2+a3) / (a4+a5+a6) # show 6 plots

#What we see in these plots is some misplaced gene expression, indicating we have RNA soup to remove!
  
#### Calculate the RNA soup fraction:
sc = autoEstCont(sc) # estimate the fraction of RNAs belonging to soup
out = adjustCounts(sc) # create a corrected count matrix

#### See which genes were most affected by our soup correction:
cntSoggy = rowSums(sc$toc > 0) # list cells with counts greater than 0 before correction for each gene
cntStrained = rowSums(out > 0) # list cells with counts greater than 0 after correction for each gene
tail(sort((cntSoggy - cntStrained)/cntSoggy), n = 10) # list the 10 most affected genes that had expression reduced in total # of cells
tail(sort(rowSums(sc$toc > out)/rowSums(sc$toc > 0)), n = 10) # list the 10 genes that had greatest overall quantities reduced

#### See how soup removal affects the genes we assessed expression patterns for earlier:
a1 <- plotChangeMap(sc, out, "CD3E") + ggtitle('CD3E')
a2 <- plotChangeMap(sc, out, "ENSSSCG00000038719") + ggtitle('IgLambdaV')
a3 <- plotChangeMap(sc, out, "CD79B") + ggtitle('CD79B')
a4 <- plotChangeMap(sc, out, "FABP6") + ggtitle('FABP6')
a5 <- plotChangeMap(sc, out, "EPCAM") + ggtitle('EPCAM')
a6 <- plotChangeMap(sc, out, "GNLY") + ggtitle('GNLY')
a7 <- plotChangeMap(sc, out, "HBB") + ggtitle('HBB')
(a1+a2+a3) / (a4+a5+a6) # show 6 plots

#### Save our strained count matrix to a new location:
write10xCounts("/home/Jayne.Wiarda/SI_PP_SC_ST/SoupX_SConly/I2strainedCounts", out, version = "3", overwrite = TRUE)
rm(dd, gg, mids, out, sc, cntSoggy, cntStrained)

#Now I do the same for the rest of the files, but I won't show most of the output plots to save some space (though I did look at these plots as a double check)

### Sample I3
#### Provide file path to Cell Ranger outputs:
sc = load10X('/home/Jayne.Wiarda/SI_PP_SC_ST/CellRangerOutputs/I3') 

#### Calculate the RNA soup fraction:
sc = autoEstCont(sc) # estimate the fraction of RNAs belonging to soup
out = adjustCounts(sc) # create a corrected count matrix

#### See which genes were most affected by our soup correction:
cntSoggy = rowSums(sc$toc > 0) # list cells with counts greater than 0 before correction for each gene
cntStrained = rowSums(out > 0) # list cells with counts greater than 0 after correction for each gene
tail(sort((cntSoggy - cntStrained)/cntSoggy), n = 10) # list the 10 most affected genes that had expression reduced in total # of cells
tail(sort(rowSums(sc$toc > out)/rowSums(sc$toc > 0)), n = 10) # list the 10 genes that had greatest overall quantities reduced

#### Save our strained count matrix to a new location:
write10xCounts("/home/Jayne.Wiarda/SI_PP_SC_ST/SoupX_SConly/I3strainedCounts", out, version = "3", overwrite = TRUE)
rm(dd, gg, mids, out, sc, cntSoggy, cntStrained)

### Sample I4
#### Provide file path to Cell Ranger outputs:
sc = load10X('/home/Jayne.Wiarda/SI_PP_SC_ST/CellRangerOutputs/I4') 

#### Calculate the RNA soup fraction:
sc = autoEstCont(sc) # estimate the fraction of RNAs belonging to soup
out = adjustCounts(sc) # create a corrected count matrix

#### See which genes were most affected by our soup correction:
cntSoggy = rowSums(sc$toc > 0) # list cells with counts greater than 0 before correction for each gene
cntStrained = rowSums(out > 0) # list cells with counts greater than 0 after correction for each gene
tail(sort((cntSoggy - cntStrained)/cntSoggy), n = 10) # list the 10 most affected genes that had expression reduced in total # of cells
tail(sort(rowSums(sc$toc > out)/rowSums(sc$toc > 0)), n = 10) # list the 10 genes that had greatest overall quantities reduced

#### Save our strained count matrix to a new location:
write10xCounts("/home/Jayne.Wiarda/SI_PP_SC_ST/SoupX_SConly/I4strainedCounts", out, version = "3", overwrite = TRUE)
rm(dd, gg, mids, out, sc, cntSoggy, cntStrained)

### Sample J2
#### Provide file path to Cell Ranger outputs:
sc = load10X('/home/Jayne.Wiarda/SI_PP_SC_ST/CellRangerOutputs/J2') 

#### Calculate the RNA soup fraction:
sc = autoEstCont(sc) # estimate the fraction of RNAs belonging to soup
out = adjustCounts(sc) # create a corrected count matrix

#### See which genes were most affected by our soup correction:
cntSoggy = rowSums(sc$toc > 0) # list cells with counts greater than 0 before correction for each gene
cntStrained = rowSums(out > 0) # list cells with counts greater than 0 after correction for each gene
tail(sort((cntSoggy - cntStrained)/cntSoggy), n = 10) # list the 10 most affected genes that had expression reduced in total # of cells
tail(sort(rowSums(sc$toc > out)/rowSums(sc$toc > 0)), n = 10) # list the 10 genes that had greatest overall quantities reduced

#### Save our strained count matrix to a new location:
write10xCounts("/home/Jayne.Wiarda/SI_PP_SC_ST/SoupX_SConly/J2strainedCounts", out, version = "3", overwrite = TRUE)
rm(dd, gg, mids, out, sc, cntSoggy, cntStrained)

### Sample J3
#### Provide file path to Cell Ranger outputs:
sc = load10X('/home/Jayne.Wiarda/SI_PP_SC_ST/CellRangerOutputs/J3') 

#### Calculate the RNA soup fraction:
sc = autoEstCont(sc) # estimate the fraction of RNAs belonging to soup
out = adjustCounts(sc) # create a corrected count matrix

#### See which genes were most affected by our soup correction:
cntSoggy = rowSums(sc$toc > 0) # list cells with counts greater than 0 before correction for each gene
cntStrained = rowSums(out > 0) # list cells with counts greater than 0 after correction for each gene
tail(sort((cntSoggy - cntStrained)/cntSoggy), n = 10) # list the 10 most affected genes that had expression reduced in total # of cells
tail(sort(rowSums(sc$toc > out)/rowSums(sc$toc > 0)), n = 10) # list the 10 genes that had greatest overall quantities reduced

#### Save our strained count matrix to a new location:
write10xCounts("/home/Jayne.Wiarda/SI_PP_SC_ST/SoupX_SConly/J3strainedCounts", out, version = "3", overwrite = TRUE)
rm(dd, gg, mids, out, sc, cntSoggy, cntStrained)

### Sample J4
#### Provide file path to Cell Ranger outputs:
sc = load10X('/home/Jayne.Wiarda/SI_PP_SC_ST/CellRangerOutputs/J4') 

#### Calculate the RNA soup fraction:
sc = autoEstCont(sc) # estimate the fraction of RNAs belonging to soup
out = adjustCounts(sc) # create a corrected count matrix

#### See which genes were most affected by our soup correction:
cntSoggy = rowSums(sc$toc > 0) # list cells with counts greater than 0 before correction for each gene
cntStrained = rowSums(out > 0) # list cells with counts greater than 0 after correction for each gene
tail(sort((cntSoggy - cntStrained)/cntSoggy), n = 10) # list the 10 most affected genes that had expression reduced in total # of cells
tail(sort(rowSums(sc$toc > out)/rowSums(sc$toc > 0)), n = 10) # list the 10 genes that had greatest overall quantities reduced

#### Save our strained count matrix to a new location:
write10xCounts("/home/Jayne.Wiarda/SI_PP_SC_ST/SoupX_SConly/J4strainedCounts", out, version = "3", overwrite = TRUE)
rm(dd, gg, mids, out, sc, cntSoggy, cntStrained)

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
## [1] stats4    stats     graphics  grDevices utils     datasets  methods  
## [8] base     
## 
## other attached packages:
##  [1] SoupX_1.6.1                 ggplot2_3.3.6              
##  [3] DropletUtils_1.16.0         SingleCellExperiment_1.18.0
##  [5] SummarizedExperiment_1.26.1 Biobase_2.56.0             
##  [7] GenomicRanges_1.48.0        GenomeInfoDb_1.32.2        
##  [9] IRanges_2.30.0              S4Vectors_0.34.0           
## [11] BiocGenerics_0.42.0         MatrixGenerics_1.8.1       
## [13] matrixStats_0.62.0         
## 
## loaded via a namespace (and not attached):
##   [1] plyr_1.8.7                igraph_1.3.4             
##   [3] lazyeval_0.2.2            sp_1.5-0                 
##   [5] splines_4.2.1             BiocParallel_1.30.3      
##   [7] listenv_0.8.0             scattermore_0.8          
##   [9] digest_0.6.29             htmltools_0.5.3          
##  [11] fansi_1.0.3               magrittr_2.0.3           
##  [13] tensor_1.5                cluster_2.1.3            
##  [15] ROCR_1.0-11               limma_3.52.2             
##  [17] globals_0.15.1            R.utils_2.12.0           
##  [19] spatstat.sparse_2.1-1     colorspace_2.0-3         
##  [21] ggrepel_0.9.1             xfun_0.31                
##  [23] dplyr_1.0.9               RCurl_1.98-1.7           
##  [25] jsonlite_1.8.0            progressr_0.10.1         
##  [27] spatstat.data_2.2-0       survival_3.3-1           
##  [29] zoo_1.8-10                glue_1.6.2               
##  [31] polyclip_1.10-0           gtable_0.3.0             
##  [33] zlibbioc_1.42.0           XVector_0.36.0           
##  [35] leiden_0.4.2              DelayedArray_0.22.0      
##  [37] Rhdf5lib_1.18.2           future.apply_1.9.0       
##  [39] HDF5Array_1.24.1          abind_1.4-5              
##  [41] scales_1.2.0              DBI_1.1.3                
##  [43] edgeR_3.38.1              spatstat.random_2.2-0    
##  [45] miniUI_0.1.1.1            Rcpp_1.0.9               
##  [47] viridisLite_0.4.0         xtable_1.8-4             
##  [49] reticulate_1.25           spatstat.core_2.4-4      
##  [51] dqrng_0.3.0               htmlwidgets_1.5.4        
##  [53] httr_1.4.3                RColorBrewer_1.1-3       
##  [55] ellipsis_0.3.2            Seurat_4.1.1             
##  [57] ica_1.0-3                 farver_2.1.1             
##  [59] pkgconfig_2.0.3           R.methodsS3_1.8.2        
##  [61] scuttle_1.6.2             uwot_0.1.11              
##  [63] deldir_1.0-6              locfit_1.5-9.6           
##  [65] utf8_1.2.2                labeling_0.4.2           
##  [67] tidyselect_1.1.2          rlang_1.0.4              
##  [69] reshape2_1.4.4            later_1.3.0              
##  [71] munsell_0.5.0             tools_4.2.1              
##  [73] cli_3.3.0                 generics_0.1.3           
##  [75] ggridges_0.5.3            evaluate_0.15            
##  [77] stringr_1.4.0             fastmap_1.1.0            
##  [79] goftest_1.2-3             yaml_2.3.5               
##  [81] knitr_1.39                fitdistrplus_1.1-8       
##  [83] purrr_0.3.4               RANN_2.6.1               
##  [85] nlme_3.1-158              pbapply_1.5-0            
##  [87] future_1.26.1             sparseMatrixStats_1.8.0  
##  [89] mime_0.12                 R.oo_1.25.0              
##  [91] compiler_4.2.1            rstudioapi_0.13          
##  [93] plotly_4.10.0             png_0.1-7                
##  [95] spatstat.utils_2.3-1      tibble_3.1.7             
##  [97] stringi_1.7.8             highr_0.9                
##  [99] rgeos_0.5-9               lattice_0.20-45          
## [101] Matrix_1.4-1              vctrs_0.4.1              
## [103] pillar_1.8.0              lifecycle_1.0.1          
## [105] rhdf5filters_1.8.0        spatstat.geom_2.4-0      
## [107] lmtest_0.9-40             RcppAnnoy_0.0.19         
## [109] data.table_1.14.2         cowplot_1.1.1            
## [111] bitops_1.0-7              irlba_2.3.5              
## [113] httpuv_1.6.5              patchwork_1.1.1          
## [115] R6_2.5.1                  promises_1.2.0.1         
## [117] KernSmooth_2.23-20        gridExtra_2.3            
## [119] parallelly_1.32.0         codetools_0.2-18         
## [121] MASS_7.3-58               assertthat_0.2.1         
## [123] rhdf5_2.40.0              withr_2.5.0              
## [125] SeuratObject_4.1.0        sctransform_0.3.3        
## [127] GenomeInfoDbData_1.2.8    mgcv_1.8-40              
## [129] parallel_4.2.1            rpart_4.1.16             
## [131] grid_4.2.1                beachmat_2.12.0          
## [133] tidyr_1.2.0               rmarkdown_2.14           
## [135] DelayedMatrixStats_1.18.0 Rtsne_0.16               
## [137] shiny_1.7.2
