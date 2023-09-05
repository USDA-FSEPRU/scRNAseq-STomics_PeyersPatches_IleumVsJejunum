library(Seurat)
library(SeuratDisk)
library(miloR)
library(SingleCellExperiment)
library(data.table)
library(ggplot2)
library(readxl)
library(dplyr)


seu <- LoadH5Seurat('/home/Jayne.Wiarda/SI_PP_SC_ST/AnnotatedSeurat/AllSamples_annotated.h5seurat')
DefaultAssay(seu) <- "SCT"
seu$tissue <- substr(seu$orig.ident, 1, 1) # add metadata on tissues; I = ileum, J = jejunum

# Start by looking at cell compositions to get an idea of differences in cellular abundances...
seu$tissue <- substr(seu$orig.ident, 1, 1) 
seu$tissue <- factor(seu$tissue, levels = c('J', 'I'))

# bar plot of cell lineage compositions within each sample
Idents(seu) <- seu$celltype
B <- rep('B', 5)
TILC <- rep('T/ILC', 15)
Myeloid <- rep('Myeloid', 3)
Epithelial <- rep('Epithelial', 6)
Stromal <- rep('Stromal', 3)
CellTypes <- c(B, TILC, Myeloid, Epithelial, Stromal)
length(CellTypes)
seu$cellLin <-seu$celltype
Idents(seu) <- seu$cellLin
names(CellTypes) <- levels(seu) # assign CellTypes to cluster numbers
seu <- RenameIdents(seu, CellTypes) # change dataset identity to cell types in new Seurat object
seu$cellLin <- Idents(seu)

Idents(seu) <- seu$orig.ident
SampleTotalCells <- prop.table(table(seu$cellLin)) # What percent of total cells are from each sample?
SamplePercents <- prop.table(table(Idents(seu),seu$cellLin), 
                             margin = 1) # What percent of cells from each cluster belong to each sample?
SamplePercents <- rbind(SamplePercents, SampleTotalCells) # add row of overall percentages to table
#rowSums(SamplePercents) # make sure all are equal to 1
SamplePercents <- t(SamplePercents) # transpose the table
par(mfrow=c(1, 1), mar=c(5, 5, 4, 8))
barplot(SamplePercents, # create stacked bar plot
        col = c('mediumorchid', 'orange', 'blue', 'forestgreen', 'burlywood4'),
        legend = rownames(SamplePercents),
        xlab = "Cluster #", 
        ylab = "Frequency within cluster", 
        las = 2,
        border = NA,
        space = 0.05,
        legend.text = TRUE, 
        args.legend = list(x = "topright", bty = "n", inset=c(-0.1, 0)))

# barplot of Sample compositions within each cell type
Idents(seu) <- seu$celltype
seu$orig.ident <- factor(seu$orig.ident, levels = c('J2', 'J3', 'J4', 'I2', 'I3', 'I4'))
SampleTotalCells <- prop.table(table(seu$orig.ident)) # What percent of total cells are from each sample?
SamplePercents <- prop.table(table(Idents(seu),seu$orig.ident), 
                             margin = 1) # What percent of cells from each cluster belong to each sample?
SamplePercents <- rbind(SamplePercents, SampleTotalCells) # add row of overall percentages to table
#rowSums(SamplePercents) # make sure all are equal to 1
SamplePercents <- t(SamplePercents) # transpose the table
par(mfrow=c(1, 1), mar=c(15, 5, 4, 8))
barplot(SamplePercents, # create stacked bar plot
        col = c('dodgerblue1', 'dodgerblue3', 'dodgerblue4',
                'indianred1', 'indianred3', 'indianred4'),
        legend = rownames(SamplePercents),
        xlab = "Cluster #", 
        ylab = "Frequency within cluster", 
        las = 2,
        border = NA,
        space = 0.05,
        legend.text = TRUE, 
        args.legend = list(x = "topright", bty = "n", inset=c(-0.15, 0)))

# Next, move on to more formal differential abundance analysis using miloR...
### Create Milo object

#First convert Seurat object to SingleCellExperiment and then to a Milo object:
set.seed(111) # set a seed for reproducibility
milo <- as.SingleCellExperiment(seu, assay = 'SCT')
milo <- Milo(milo)

#Also incorporate the shared nearest neighbors (SNN) graph calculated in Seurat into the Milo object. We will use this to calculate cell neighborhoods.

#First calculate PCs to use for finding neighbors:
  
##Start by again finding number of PCs to use:

pct <- seu[["pca"]]@stdev / sum(seu[["pca"]]@stdev) * 100 # find standard deviation for each PC
cumu <- cumsum(pct) # find cumulative percentages for PCs
co1 <- which(cumu > 90 & pct < 5)[1] # find PC representing cumulative percent >90% and less than 5% associated with the single PC
co1 # list PC
co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1 # find last PC where change in percent variation is more than 0.1%
co2 # list PC
pcs <- min(co1, co2) # find the minimum PC from the 2 methods used above
pcs # list PC
PCdims <- 1:pcs # use the minimum PC from the quantitative method above to set the PCs for subsequent steps
rm(pct, cumu, co1, co2) 

#Next, find and add neighbors:
seu <- FindNeighbors(seu, assay = 'integrated', k.param = 20, dims = PCdims) # calculate nearest neighbors if seu@graphs$integrated_snn does not exist yet
miloR::graph(milo) <- miloR::graph(buildFromAdjacency(seu@graphs$integrated_snn, k=20))


### Create & visualize cell neighborhoods:

#Start by creating cell neighborhoods. Parameters of prop, k, and d may be modified slightly depending on the dataset used. Higher proportions of cells (prop) will take longer to run, but may require up to a value of 0.2 for smaller datasets. We choose to set k and d parameters according to those used to calculate our SNN graph and 'significant' PCs in our previous Seurat analysis.

#Now generate the neighborhoods:

set.seed(111) # set a seed for reproducibility of neighborhood generation
milo <- makeNhoods(milo,
                   prop = 0.2, # sample 20% of cells...probably safe to lower as far as 0.05 for datasets with >30k cells...may consider using proportions up to 0.2 if that helps optimize neighborhood size distribution peak; higher prop will increase run time. We will keep prop at 0.2 since we have some rare cell populations we want to make sure are represented.
                   k = 20, # set to k = 20 because for Seurat FindNeighbors() we used default k.param = 20
                   d = pcs, # set to PCs used to find neighbors in Seurat
                   refined = TRUE, # always use refined unless you use graph-based data batch correction, then consider either-or
                   refinement_scheme="graph") # use graph-based approach so bottlenecking step calcNhoodDistance() isn't required

#Now that we have calculated our cell neighborhoods, let's look at their sizes. Ideally, peak size should fall between 50-100 cells per neighborhood but may be less for extremely small datasets:
plotNhoodSizeHist(milo) # ideally have peak of distribution between 50 and 100...otherwise consider increasing k or prop...peak may be <50 for small datasets

#Now let's move on to look at these cell neighborhoods overlaid onto UMAP coordinates:

milo <- buildNhoodGraph(milo)
plotNhoodGraph(milo, layout = 'UMAP')

### Count cells in each neighborhood

#Now let's do a head count of which cells came from each of our samples within each of our detected cell neighborhoods:

milo <- countCells(milo, meta.data = data.frame(colData(milo)), sample="orig.ident")
head(nhoodCounts(milo))

### Create experimental design

#Create a model of experimental design variables:

milo_design <- data.frame(colData(milo))[,c("orig.ident", "tissue")]
milo_design <- distinct(milo_design)
rownames(milo_design) <- milo_design$orig.ident
milo_design

### Perform DA testing

#Perform DA testing on each neighborhood:

da_results <- testNhoods(milo,
                         design = ~ tissue,
                         design.df = milo_design,
                         fdr.weighting = 'graph-overlap', # make sure to use graph-overlap since we used graph-based approach for makeNhoods()
                         reduced.dim = 'PCA')
head(da_results)

#Make a histogram of p-values found across cell neighborhoods:

ggplot(da_results, aes(PValue)) + geom_histogram(bins=100)

#Make a volcano plot of DA. Each dot is one cell neighborhood:

ggplot(da_results, aes(logFC, -log10(FDR))) +
  geom_point() +
  geom_hline(yintercept = 1) ## Mark significance threshold (10% FDR)

#Overlay logFC scores onto cell neighborhood central coordinates on t-SNE & UMAP plots:

plotNhoodGraphDA(milo, da_results, layout="UMAP",alpha=0.05)

#And we can also look at all cell neighborhoods on a bee swarm plot:

plotDAbeeswarm(da_results, alpha = 0.05)

#The problem right now is that we know which cell neighborhoods are DA and can guess at which cell types these correspond to, but we don't know for certain. Therefore, it may be helpful to annotate our cell neighborhoods and then re-assess DA....

### Annotate cell neighborhoods

#Annotate cell neighborhoods by finding which neighborhoods of most cells belonging to a specific cell ID. 

#Start by calculating the percentage of cells belonging to each annotated cell type within each neighborhood. We will record the annotation that has the largest percentage in each neighborhood:
  
da_results <- annotateNhoods(milo, da_results, coldata_col = "celltype")
head(da_results)

#Create a histogram to look at the largest percentages for a single cell type within each cell neighborhood:
  
ggplot(da_results, aes(celltype_fraction)) + geom_histogram(bins=100)

#Based on this graph, we need to set a cut-off value for what percentage of cells in a neighborhood must share a single ID to be assigned as that cell type annotation. In this case, we will set our cutoff at 0.7.

#Based on this criteria, any cell neighborhood with >70% of cells belonging to a single annotation will be assigned to that identity. Any cell neighborhood with <70% of cells belonging to a single annotation will be catergorized as 'Mixed' cell neighborhoods:
  
da_results$celltype <- ifelse(da_results$celltype_fraction < 0.7, "Mixed", da_results$celltype)
da_results$celltype <- factor(da_results$celltype,
                              levels = rev(c('Mixed', 'Activated B cells', 'Cycling B cells', 'Resting B cells', 'Transitioning B cells', 'Antibody-secreting cells',
                                             'Cycling CD4 ab T cells', 'Cycling CD8 ab T cells', 'Cycling gd T cells', 'Cycling group 1 ILCs',
                                             'Cytotoxic CD8 ab T cells', 'Cytotoxic gd T cells', 'Cytotoxic group 1 ILCs', 
                                             'Non-naive CD8 ab T cells', 'Non-naive gd T cells', 'Non-naive group 1 ILCs', 
                                             'SELLhi gd T cells', 'CD2neg GD T cells', 
                                             'Naive CD4/CD8 ab T cells', 'Non-naive CD4 ab T cells', 'Follicular CD4 ab T cells', 'Group 3 ILCs',
                                             'Dendritic cells', 'Macrophages', 'Mast cells',
                                             'Crypt cells', 'Enterocytes', 'BEST4 enterocytes', 'Goblet cells', 'NEUROD1lo EE cells', 'NEUROD1hi EE cells',
                                             'Endothelial cells', 'Fibroblasts', 'Muscle cells')))

### Plot DA across annotated cell neighborhoods:

#Make a bee swarm plot:
  
plotDAbeeswarm(da_results, group.by = "celltype", alpha = 0.05)

### Further summarization of DA results:

#Let's further summarize only those neighborhoods with an annotated cell type. Start by subsetting only non-mixed neighborhoods, creating a column defining DA significance (FDR < 0.05), and creating a column indicating fold-change towards enrichment in one group vs another:

da_sum <- subset(da_results, celltype_fraction > 0.7)
da_sum <- da_sum %>% mutate(significant = case_when(FDR < 0.05 ~ 'Sig', FDR >= 0.05 ~ 'NonSig'))
da_sum <- da_sum %>% mutate(FC = case_when(logFC > 0 ~ 'Jejunum', logFC < 0 ~ 'Ileum', logFC == 0 ~ 'Neutral'))
da_sum$result <- paste(da_sum$FC, da_sum$significant, sep = '_')
da_sum$result <- replace(da_sum$result, da_sum$result == 'Ileum_NonSig', 'NonSig')
da_sum$result <- replace(da_sum$result, da_sum$result == 'Jejunum_NonSig', 'NonSig')
table(da_sum$result, da_sum$celltype) # see summary of results per cell type

### Save DA results:

write_xlsx(da_results, '/home/Jayne.Wiarda/SI_PP_SC_ST/DifferentialAbundance/DAresults_AllCells.xlsx')

#Save Milo object:

saveRDS(milo, '/home/Jayne.Wiarda/SI_PP_SC_ST/DifferentialAbundance/milo_AllCells.rds')

### View session information

sessionInfo()
#R version 4.2.2 Patched (2022-11-10 r83330)
#Platform: x86_64-pc-linux-gnu (64-bit)
#Running under: Ubuntu 22.04.1 LTS

#Matrix products: default
#BLAS:   /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3
#LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/libopenblasp-r0.3.20.so

#locale:
#  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C               LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8     LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8    LC_PAPER=en_US.UTF-8       LC_NAME=C                  LC_ADDRESS=C              
#[10] LC_TELEPHONE=C             LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       

#attached base packages:
#  [1] stats4    stats     graphics  grDevices utils     datasets  methods   base     

#other attached packages:
#  [1] dplyr_1.0.10                readxl_1.4.1                ggplot2_3.3.6               data.table_1.14.2           SingleCellExperiment_1.18.0 SummarizedExperiment_1.26.1 Biobase_2.56.0              GenomicRanges_1.48.0       
#[9] GenomeInfoDb_1.32.4         IRanges_2.30.1              S4Vectors_0.34.0            BiocGenerics_0.42.0         MatrixGenerics_1.8.1        matrixStats_0.62.0          miloR_1.4.0                 edgeR_3.38.4               
#[17] limma_3.52.3                SeuratDisk_0.0.0.9020       sp_1.5-0                    SeuratObject_4.1.1          Seurat_4.1.1               

#loaded via a namespace (and not attached):
#  [1] scattermore_0.8           R.methodsS3_1.8.2         ragg_1.2.2                pkgmaker_0.32.2           tidyr_1.2.1               bit64_4.0.5               knitr_1.40                irlba_2.3.5               DelayedArray_0.22.0      
#[10] R.utils_2.12.0            rpart_4.1.16              KEGGREST_1.36.3           RCurl_1.98-1.8            doParallel_1.0.17         generics_0.1.3            ScaledMatrix_1.4.1        cowplot_1.1.1             RSQLite_2.2.17           
#[19] RANN_2.6.1                future_1.28.0             bit_4.0.4                 phylobase_0.8.10          spatstat.data_2.2-0       xml2_1.3.3                httpuv_1.6.6              assertthat_0.2.1          viridis_0.6.2            
#[28] xfun_0.33                 hms_1.1.2                 evaluate_0.16             promises_1.2.0.1          fansi_1.0.3               progress_1.2.2            igraph_1.3.4              DBI_1.1.3                 htmlwidgets_1.5.4        
#[37] spatstat.geom_2.4-0       purrr_0.3.4               ellipsis_0.3.2            V8_4.2.2                  annotate_1.74.0           gridBase_0.4-7            locfdr_1.1-8              deldir_1.0-6              sparseMatrixStats_1.8.0  
#[46] vctrs_0.4.1               ROCR_1.0-11               abind_1.4-5               cachem_1.0.6              withr_2.5.0               ggforce_0.3.4             progressr_0.11.0          sctransform_0.3.4         prettyunits_1.1.1        
#[55] goftest_1.2-3             softImpute_1.4-1          cluster_2.1.4             ape_5.6-2                 lazyeval_0.2.2            crayon_1.5.1              genefilter_1.78.0         hdf5r_1.3.5               labeling_0.4.2           
#[64] pkgconfig_2.0.3           tweenr_2.0.2              nlme_3.1-159              vipor_0.4.5               rlang_1.0.5               globals_0.16.1            lifecycle_1.0.2           miniUI_0.1.1.1            registry_0.5-1           
#[73] rsvd_1.0.5                dichromat_2.0-0.1         cellranger_1.1.0          polyclip_1.10-0           lmtest_0.9-40             rngtools_1.5.2            Matrix_1.5-1              Rhdf5lib_1.18.2           zoo_1.8-10               
#[82] beeswarm_0.4.0            ggridges_0.5.3            png_0.1-7                 viridisLite_0.4.1         bitops_1.0-7              R.oo_1.25.0               rncl_0.8.6                KernSmooth_2.23-20        rhdf5filters_1.8.0       
#[91] Biostrings_2.64.1         blob_1.2.3                DelayedMatrixStats_1.18.0 stringr_1.4.1             zinbwave_1.18.0           parallelly_1.32.1         spatstat.random_2.2-0     beachmat_2.12.0           scales_1.2.1             
#[100] memoise_2.0.1             magrittr_2.0.3            plyr_1.8.7                ica_1.0-3                 howmany_0.3-1             zlibbioc_1.42.0           compiler_4.2.2            dqrng_0.3.0               RColorBrewer_1.1-3       
#[109] fitdistrplus_1.1-8        cli_3.4.0                 ade4_1.7-19               XVector_0.36.0            listenv_0.8.0             patchwork_1.1.2           pbapply_1.5-0             MASS_7.3-58.1             mgcv_1.8-40              
#[118] tidyselect_1.1.2          stringi_1.7.8             textshaping_0.3.6         yaml_2.3.5                BiocSingular_1.12.0       locfit_1.5-9.6            ggrepel_0.9.1             grid_4.2.2                tools_4.2.2              
#[127] future.apply_1.9.1        parallel_4.2.2            rstudioapi_0.14           uuid_1.1-0                foreach_1.5.2             RNeXML_2.4.7              gridExtra_2.3             farver_2.1.1              Rtsne_0.16               
#[136] ggraph_2.0.6              digest_0.6.29             rgeos_0.5-9               shiny_1.7.2               Rcpp_1.0.9                scuttle_1.6.3             later_1.3.0               RcppAnnoy_0.0.19          AnnotationDbi_1.58.0     
#[145] httr_1.4.4                kernlab_0.9-31            colorspace_2.0-3          XML_3.99-0.10             tensor_1.5                reticulate_1.26           clusterExperiment_2.16.0  splines_4.2.2             statmod_1.4.37           
#[154] uwot_0.1.14               spatstat.utils_2.3-1      graphlayouts_0.8.1        mapproj_1.2.8             systemfonts_1.0.4         plotly_4.10.0             xtable_1.8-4              jsonlite_1.8.0            tidygraph_1.2.2          
#[163] R6_2.5.1                  pillar_1.8.1              htmltools_0.5.3           mime_0.12                 NMF_0.24.0                glue_1.6.2                fastmap_1.1.0             BiocParallel_1.30.3       BiocNeighbors_1.14.0     
#[172] codetools_0.2-18          maps_3.4.0                utf8_1.2.2                lattice_0.20-45           spatstat.sparse_2.1-1     tibble_3.1.8              curl_4.3.2                ggbeeswarm_0.6.0          leiden_0.4.3             
#[181] gtools_3.9.3              survival_3.4-0            rmarkdown_2.16            munsell_0.5.0             rhdf5_2.40.0              GenomeInfoDbData_1.2.8    iterators_1.0.14          HDF5Array_1.24.2          reshape2_1.4.4           
#[190] gtable_0.3.1              spatstat.core_2.4-4    