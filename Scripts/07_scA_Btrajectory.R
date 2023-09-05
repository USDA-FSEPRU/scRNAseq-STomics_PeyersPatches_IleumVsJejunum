### Load required software packages

library(Seurat)
library(SeuratObject)
library(ggplot2)
library(SeuratDisk)
library(SingleCellExperiment)
library(slingshot)
library(viridis)
library(dplyr)
library(tidyverse)
#library(scales)

## Load & process Seurat object

#Load a processed Seurat object of only B cells (from both ileum + jejunum samples combined):
  
seu <- LoadH5Seurat('/home/Jayne.Wiarda/SI_PP_SC_ST/AnnotatedSeurat/B.h5seurat')
DefaultAssay(seu) <- 'SCT'
Idents(seu) <- seu$celltype
DimPlot(seu, reduction = 'umap', cols = c('cyan4', 'gold3', 'chartreuse4', 'deeppink4', 'sandybrown'))
#DimPlot(seu, reduction = 'pca', cols = c('cyan4', 'gold3', 'chartreuse4', 'deeppink4', 'sandybrown'))

#Find how many PC dimensions to use on dataset:

pct <- seu[["pca"]]@stdev / sum(seu[["pca"]]@stdev) * 100 # find standard deviation for each PC
cumu <- cumsum(pct) # find cumulative percentages for PCs
co1 <- which(cumu > 90 & pct < 5)[1] # find PC representing cumulative percent >90% and less than 5% associated with the single PC
co1 # list PC
co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1 # find last PC where change in percent variation is more than 0.1%
co2 # list PC
pcs <- min(co1, co2) # find the minimum PC from the 2 methods used above
pcs # list PC
PCdims <- 1:pcs # use the minimum PC from the quantitative method above to set the PCs for subsequent steps
rm(co1, co2, cumu, pcs, pct)

#Make a new metadata slot for tissue:
seu$tissue <- substr(seu$orig.ident, 1, 1) # I = ileum, J = jejunum

## Slingshot trajectory analysis

#Convert to a SingleCellExperiment object:

sc <- as.SingleCellExperiment(seu, assay = 'SCT')
# add original dimensionality reductions to sce object
# could opt to re-calculate these for each data subset if desired
reducedDims(sc) <- SimpleList(PCA=Embeddings(seu, 'pca')[,PCdims], 
                              UMAP = Embeddings(seu, 'umap'), 
                              TSNE = Embeddings(seu, 'tsne'))

#Run Slingshot trajectory analysis:

sc <- slingshot(sc, reducedDim = 'PCA', clusterLabels = 'celltype', start.clus = 'Resting B cells')  # leaving clusters unspecified would allow only a single trajectory to be fit to all cells
summary(sc$slingPseudotime_1) # see summary of pseudotime.... may have more than one trajectory to look at (e.g. $slingPseudotime_2) if clusters were specified in the slingshot() command
summary(sc$slingPseudotime_2) 
#summary(sc$slingPseudotime_3) # and so on for as many trajectories as are made

## Plot trajectories

#Plot trajectory 1:

#pt <- slingPseudotime(sc)
#nms <- colnames(pt)
#pal <- mako(100, end = 1)
#colors <- rev(mako(50, alpha = 1)) 
#plotcol <- colors[cut(sc$slingPseudotime_1, breaks=50)]

# plot on pca
#plot(reducedDims(sc)$PCA, col = plotcol, pch=16, asp = 1, cex = 0.2) 

#plot(reducedDims(sc)$PCA, col = plotcol, pch=16, asp = 1, cex = 0.2) 
#lines(SlingshotDataSet(sc), lwd=2, 
#      #type = 'lineages', 
#      col = 'black')

# plot on umap
#plot(reducedDims(sc)$UMAP, col = plotcol, pch=16, asp = 1, cex = 0.2)

#Plot trajectory 2:

#pt <- slingPseudotime(sc)
#nms <- colnames(pt)
#pal <- mako(100, end = 1)
#colors <- rev(mako(50, alpha = 1)) 
#plotcol <- colors[cut(sc$slingPseudotime_2, breaks=50)]

# plot on pca
#plot(reducedDims(sc)$PCA, col = plotcol, pch=16, asp = 1, cex = 0.2) 

#plot(reducedDims(sc)$PCA, col = plotcol, pch=16, asp = 1, cex = 0.2) 
#lines(SlingshotDataSet(sc), lwd=2, 
#      #type = 'lineages', 
#      col = 'black')

# plot on umap
#plot(reducedDims(sc)$UMAP, col = plotcol, pch=16, asp = 1, cex = 0.2)

## Plot with Seurat

seu$pseudotime1 <- sc$slingPseudotime_1
seu$pseudotime2 <- sc$slingPseudotime_2

notcells <- colnames(seu)[is.na(seu$pseudotime1)]
dat <- data.frame(colnames(seu), seu$celltype)
colnames(dat) <- c('barcode', 'celltype')
dat$celltype <- as.character(dat$celltype)
dat$celltype[dat$barcode %in% notcells] <- 'null'
seu$pseudotime1_celltype <- dat$celltype

FeaturePlot(seu, features = c('pseudotime1'), order = TRUE) + 
  scale_colour_viridis(option = 'mako', begin = 0.1, direction = -1, na.value="grey85") 
DimPlot(seu, group.by = 'pseudotime1_celltype',
        cols = c('cyan4', 'gold3', 'grey85', 'chartreuse4', 'deeppink4'))

notcells <- colnames(seu)[is.na(seu$pseudotime2)]
dat <- data.frame(colnames(seu), seu$celltype)
colnames(dat) <- c('barcode', 'celltype')
dat$celltype <- as.character(dat$celltype)
dat$celltype[dat$barcode %in% notcells] <- 'null'
seu$pseudotime2_celltype <- dat$celltype

FeaturePlot(seu, features = c('pseudotime2'), order = TRUE) + 
  scale_colour_viridis(option = 'mako', begin = 0.1, direction = -1, na.value="grey85") 
DimPlot(seu, group.by = 'pseudotime2_celltype',
        cols = c('cyan4', 'sandybrown', 'gold3', 'grey85', 'chartreuse4', 'deeppink4'))
  
dat <- as.data.frame(cbind(slingPseudotime(sc), sc$tissue))
colnames(dat) <- c('Trajectory1', 'Trajectory2', 'Tissue')
dat$Trajectory1 <- as.numeric(dat$Trajectory1)
dat$Trajectory2 <- as.numeric(dat$Trajectory2)

ggplot(dat, aes(x=Trajectory1, color = Tissue, fill = Tissue)) + 
  geom_density(alpha=0.1) +
  scale_color_manual(values=c('magenta4', 'forestgreen')) +
  scale_fill_manual(values=c('magenta4', 'forestgreen')) +
  theme_classic() +
  scale_x_continuous(expand=c(0,0)) +
  scale_y_continuous(expand=c(0,0))

ggplot(dat, aes(x=Trajectory2, color = Tissue, fill = Tissue)) + 
  geom_density(alpha=0.1) +
  scale_color_manual(values=c('magenta4', 'forestgreen')) +
  scale_fill_manual(values=c('magenta4', 'forestgreen')) +
  theme_classic() +
  scale_x_continuous(expand=c(0,0)) +
  scale_y_continuous(expand=c(0,0))

dat <- cbind(dat, sc$orig.ident)
colnames(dat) <- c('Trajectory1', 'Trajectory2', 'Tissue', 'SampleID')

ggplot(dat, aes(x=Trajectory1, color = SampleID)) + 
  geom_density(alpha=0.1) +
  scale_color_manual(values=c('mediumpurple', 'magenta4', 'mediumpurple4', 'palegreen3', 'forestgreen', 'darkolivegreen')) +
  theme_classic() +
  scale_x_continuous(expand=c(0,0)) +
  scale_y_continuous(expand=c(0,0))

ggplot(dat, aes(x=Trajectory2, color = SampleID)) + 
  geom_density(alpha=0.1) +
  scale_color_manual(values=c('mediumpurple', 'magenta4', 'mediumpurple4', 'palegreen3', 'forestgreen', 'darkolivegreen')) +
  theme_classic() +
  scale_x_continuous(expand=c(0,0)) +
  scale_y_continuous(expand=c(0,0))

dat <- cbind(dat, sc$celltype)
colnames(dat) <- c('Trajectory1', 'Trajectory2', 'Tissue', 'SampleID', 'CellType')

ggplot(dat, aes(x=Trajectory1, color = CellType, fill = CellType)) + 
  geom_density(alpha=0.1) +
  scale_color_manual(values=c('cyan4', 'gold3', 'chartreuse4', 'deeppink4', 'sandybrown')) +
  scale_fill_manual(values=c('cyan4', 'gold3', 'chartreuse4', 'deeppink4', 'sandybrown')) +
  theme_classic() +
  scale_x_continuous(expand=c(0,0)) +
  scale_y_continuous(expand=c(0,0))

ggplot(dat, aes(x=Trajectory2, color = CellType, fill = CellType)) + 
  geom_density(alpha=0.1) +
  scale_color_manual(values=c('cyan4', 'gold3', 'chartreuse4', 'deeppink4', 'sandybrown')) +
  scale_fill_manual(values=c('cyan4', 'gold3', 'chartreuse4', 'deeppink4', 'sandybrown')) +
  theme_classic() +
  scale_x_continuous(expand=c(0,0)) +
  scale_y_continuous(expand=c(0,0))

## Break down trajectories into range increments:

dat <- data.frame(seu$pseudotime1, seu$pseudotime2)
max <- max(seu$pseudotime1, na.rm = T)
dat <- dat %>% mutate(psuedotime1_bin = cut(seu.pseudotime1, 
                                            breaks=c(0*max, 
                                                     0.2*max,
                                                     0.4*max,
                                                     0.6*max,
                                                     0.8*max,
                                                     1*max)))
dat$psuedotime1_bin <- as.character(dat$psuedotime1_bin)
dat$psuedotime1_bin <- replace_na(dat$psuedotime1_bin, 'Non-trajectory B')

max <- max(seu$pseudotime2, na.rm = T)
dat <- dat %>% mutate(psuedotime2_bin = cut(seu.pseudotime2, 
                                            breaks=c(0*max, 
                                                     0.2*max,
                                                     0.4*max,
                                                     0.6*max,
                                                     0.8*max,
                                                     1*max)))
dat$psuedotime2_bin <- as.character(dat$psuedotime2_bin)
dat$psuedotime2_bin <- replace_na(dat$psuedotime2_bin, 'Non-trajectory B')

seu <- AddMetaData(seu, metadata = dat)

DimPlot(seu, reduction = 'umap',
        group.by = 'psuedotime1_bin',
        cols = c('red', 'limegreen', 'blue', 'orange', 'gold', 'grey85'))

DimPlot(seu, reduction = 'umap',
        group.by = 'psuedotime2_bin',
        cols = c('red', 'gold', 'limegreen', 'blue', 'orange', 'grey85'))

### Save data:

SaveH5Seurat(seu, '/home/Jayne.Wiarda/SI_PP_SC_ST/PseudotimeTrajectory/B_traj.h5seurat', overwrite = TRUE)
saveRDS(sc, '/home/Jayne.Wiarda/SI_PP_SC_ST/PseudotimeTrajectory/B_traj_sce.rds')

### View session information

sessionInfo()
#R version 4.2.2 Patched (2022-11-10 r83330)
#Platform: x86_64-pc-linux-gnu (64-bit)
#Running under: Ubuntu 22.04.1 LTS

#Matrix products: default
#BLAS:   /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3
#LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/libopenblasp-r0.3.20.so

#locale:
#  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C               LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8     LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8    LC_PAPER=en_US.UTF-8       LC_NAME=C                  LC_ADDRESS=C               LC_TELEPHONE=C             LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       

#attached base packages:
#  [1] stats4    stats     graphics  grDevices utils     datasets  methods   base     

#other attached packages:
#  [1] forcats_0.5.2               stringr_1.4.1               purrr_0.3.4                 readr_2.1.2                 tidyr_1.2.1                 tibble_3.1.8                tidyverse_1.3.2             dplyr_1.0.10                viridis_0.6.2               viridisLite_0.4.1           slingshot_2.4.0             TrajectoryUtils_1.4.0      
#[13] princurve_2.1.6             SingleCellExperiment_1.18.0 SummarizedExperiment_1.26.1 Biobase_2.56.0              GenomicRanges_1.48.0        GenomeInfoDb_1.32.4         IRanges_2.30.1              S4Vectors_0.34.0            BiocGenerics_0.42.0         MatrixGenerics_1.8.1        matrixStats_0.62.0          SeuratDisk_0.0.0.9020      
#[25] ggplot2_3.3.6               sp_1.5-0                    SeuratObject_4.1.1          Seurat_4.1.1               

#loaded via a namespace (and not attached):
#  [1] scattermore_0.8           ragg_1.2.2                pkgmaker_0.32.2           bit64_4.0.5               knitr_1.40                irlba_2.3.5               DelayedArray_0.22.0       data.table_1.14.2         rpart_4.1.16              KEGGREST_1.36.3           RCurl_1.98-1.8            doParallel_1.0.17         generics_0.1.3           
#[14] ScaledMatrix_1.4.1        cowplot_1.1.1             RSQLite_2.2.17            RANN_2.6.1                miloR_1.4.0               future_1.28.0             tzdb_0.3.0                bit_4.0.4                 phylobase_0.8.10          lubridate_1.9.0           spatstat.data_2.2-0       xml2_1.3.3                httpuv_1.6.6             
#[27] assertthat_0.2.1          gargle_1.2.1              xfun_0.33                 hms_1.1.2                 evaluate_0.16             promises_1.2.0.1          fansi_1.0.3               progress_1.2.2            readxl_1.4.1              dbplyr_2.2.1              igraph_1.3.4              DBI_1.1.3                 htmlwidgets_1.5.4        
#[40] spatstat.geom_2.4-0       googledrive_2.0.0         ellipsis_0.3.2            backports_1.4.1           annotate_1.74.0           gridBase_0.4-7            locfdr_1.1-8              deldir_1.0-6              sparseMatrixStats_1.8.0   vctrs_0.4.1               ROCR_1.0-11               abind_1.4-5               cachem_1.0.6             
#[53] withr_2.5.0               ggforce_0.3.4             progressr_0.11.0          sctransform_0.3.4         prettyunits_1.1.1         goftest_1.2-3             softImpute_1.4-1          cluster_2.1.4             ape_5.6-2                 lazyeval_0.2.2            crayon_1.5.1              genefilter_1.78.0         hdf5r_1.3.5              
#[66] edgeR_3.38.4              pkgconfig_2.0.3           labeling_0.4.2            tweenr_2.0.2              nlme_3.1-159              vipor_0.4.5               rlang_1.0.6               globals_0.16.1            lifecycle_1.0.2           miniUI_0.1.1.1            registry_0.5-1            modelr_0.1.10             rsvd_1.0.5               
#[79] cellranger_1.1.0          polyclip_1.10-0           lmtest_0.9-40             rngtools_1.5.2            Matrix_1.5-1              Rhdf5lib_1.18.2           zoo_1.8-10                reprex_2.0.2              beeswarm_0.4.0            googlesheets4_1.0.1       ggridges_0.5.3            png_0.1-7                 bitops_1.0-7             
#[92] rhdf5filters_1.8.0        rncl_0.8.6                KernSmooth_2.23-20        Biostrings_2.64.1         blob_1.2.3                DelayedMatrixStats_1.18.0 zinbwave_1.18.0           parallelly_1.32.1         spatstat.random_2.2-0     beachmat_2.12.0           scales_1.2.1              memoise_2.0.1             magrittr_2.0.3           
#[105] plyr_1.8.7                ica_1.0-3                 howmany_0.3-1             zlibbioc_1.42.0           compiler_4.2.2            RColorBrewer_1.1-3        fitdistrplus_1.1-8        cli_3.4.0                 ade4_1.7-19               XVector_0.36.0            listenv_0.8.0             patchwork_1.1.2           pbapply_1.5-0            
#[118] MASS_7.3-58.1             mgcv_1.8-40               tidyselect_1.1.2          stringi_1.7.8             textshaping_0.3.6         yaml_2.3.5                BiocSingular_1.12.0       locfit_1.5-9.6            ggrepel_0.9.1             grid_4.2.2                timechange_0.2.0          tools_4.2.2               future.apply_1.9.1       
#[131] parallel_4.2.2            rstudioapi_0.14           uuid_1.1-0                foreach_1.5.2             RNeXML_2.4.7              gridExtra_2.3             farver_2.1.1              Rtsne_0.16                ggraph_2.0.6              digest_0.6.29             rgeos_0.5-9               shiny_1.7.2               Rcpp_1.0.9               
#[144] broom_1.0.1               later_1.3.0               RcppAnnoy_0.0.19          httr_1.4.4                AnnotationDbi_1.58.0      kernlab_0.9-31            colorspace_2.0-3          rvest_1.0.3               fs_1.5.2                  XML_3.99-0.10             tensor_1.5                reticulate_1.26           clusterExperiment_2.16.0 
#[157] splines_4.2.2             uwot_0.1.14               spatstat.utils_2.3-1      graphlayouts_0.8.1        plotly_4.10.0             systemfonts_1.0.4         xtable_1.8-4              jsonlite_1.8.0            tidygraph_1.2.2           R6_2.5.1                  pillar_1.8.1              htmltools_0.5.3           mime_0.12                
#[170] NMF_0.24.0                glue_1.6.2                fastmap_1.1.0             BiocParallel_1.30.3       BiocNeighbors_1.14.0      codetools_0.2-18          utf8_1.2.2                lattice_0.20-45           spatstat.sparse_2.1-1     ggbeeswarm_0.6.0          leiden_0.4.3              gtools_3.9.3              survival_3.4-0           
#[183] limma_3.52.3              rmarkdown_2.16            munsell_0.5.0             rhdf5_2.40.0              GenomeInfoDbData_1.2.8    iterators_1.0.14          HDF5Array_1.24.2          haven_2.5.1               reshape2_1.4.4            gtable_0.3.1              spatstat.core_2.4-4   