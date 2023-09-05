library(slingshot)
library(SingleCellExperiment)
library(tradeSeq)
library(writexl)
library(Seurat)
library(SeuratDisk)
library(ggplot2)
library(dplyr)
library(data.table)
library(clusterExperiment)
library(ComplexHeatmap)
library(circlize)
library(tidyr)

#  refer to https://kstreet13.github.io/bioc2020trajectories/articles/workshopTrajectories.html#trajectory-inference-1

seu <- LoadH5Seurat('/home/Jayne.Wiarda/SI_PP_SC_ST/PseudotimeTrajectory/B_traj.h5seurat')
DefaultAssay(seu) <- 'integrated'
genes <- VariableFeatures(seu) # find most variable genes to test...not testing all genes will severely cut down on run time

sce <- readRDS('/home/Jayne.Wiarda/SI_PP_SC_ST/PseudotimeTrajectory/B_traj_sce.rds')

BPPARAM <- BiocParallel::bpparam()
BPPARAM # lists current options
BPPARAM$workers <- 78 

# Before running DE tests, identify an appropriate k value to use. Refer to link at start of script for info
# look for elbow where decrease begins to stabilize in plots

icMat <- evaluateK(counts = as.matrix(assays(sce)$counts),
                   pseudotime = as.matrix(sce$slingshot@assays@data$pseudotime),
                   cellWeights = as.matrix(sce$slingshot@assays@data$weights),
                   conditions = factor(colData(sce)$tissue),
                   nGenes = 300,
                   k = 3:7)
nknots = 5 # use 5 knots for this data

# Do test without consideration of tissue

set.seed(123)
sce <- fitGAM(counts = as.matrix(assays(sce)$counts),
              sds = sce$slingshot,
              conditions = factor(colData(sce)$tissue),
              nknots = nknots,
              genes = genes,
              parallel = TRUE,
              verbose = TRUE,
              BPPARAM = BPPARAM)
mean(rowData(sce)$tradeSeq$converged)

rowData(sce)$assocRes <- associationTest(sce, lineages = TRUE, l2fc = log2(2))
assocRes <- rowData(sce)$assocRes
assocRes$gene <- rownames(assocRes)
assocRes$padj_lineage1_conditionI <- p.adjust(assocRes$pvalue_lineage1_conditionI, "fdr")
assocRes$padj_lineage1_conditionJ <- p.adjust(assocRes$pvalue_lineage1_conditionJ, "fdr")
assocRes$padj_lineage2_conditionI <- p.adjust(assocRes$pvalue_lineage2_conditionI, "fdr")
assocRes$padj_lineage2_conditionJ <- p.adjust(assocRes$pvalue_lineage2_conditionJ, "fdr")
write_xlsx(assocRes, '/home/Jayne.Wiarda/SI_PP_SC_ST/PseudotimeTrajectory/pseudotimeDE.xlsx')
saveRDS(sce, '/home/Jayne.Wiarda/SI_PP_SC_ST/PseudotimeTrajectory/pseudotime_GAM_model.rds')

# Find DE genes between tissues across pseudotime

condRes <- conditionTest(sce, 
                         global = FALSE,
                         lineage = TRUE,
                         pairwise = TRUE,
                         l2fc = log2(2))
condRes$padj_lineage1 <- p.adjust(condRes$pvalue_lineage1, "fdr")
condRes$padj_lineage2 <- p.adjust(condRes$pvalue_lineage2, "fdr")
condRes$gene <- rownames(condRes)
write_xlsx(condRes, '/home/Jayne.Wiarda/SI_PP_SC_ST/PseudotimeTrajectory/pseudotimeDE_condition.xlsx')

## Plot results in heatmaps

### Trajectory 1
IGenes <-  assocRes$gene[
  which(assocRes$padj_lineage1_conditionI < 0.05)
]
JGenes <-  assocRes$gene[
  which(assocRes$padj_lineage1_conditionJ < 0.05)
]

length(IGenes) # How many DE genes for ileum?
#1125
length(JGenes) # How many DE genes for jejunum?
#1179
length(intersect(IGenes, JGenes)) # How many DE genes are shared?
#958
sameGenes <- intersect(IGenes, JGenes)

dat <- data.frame(seu$pseudotime1, seu$celltype, seu$tissue)
max <- max(seu$pseudotime1, na.rm = T)
dat <- dat %>% mutate(psuedotime1_bin10 = cut(seu.pseudotime1, 
                                              breaks=c(0*max, 
                                                       0.1*max,
                                                       0.2*max,
                                                       0.3*max,
                                                       0.4*max,
                                                       0.5*max,
                                                       0.6*max,
                                                       0.7*max,
                                                       0.8*max,
                                                       0.9*max,
                                                       1*max)))
dat <- na.omit(dat, cols = "pseudotime1")
dat$psuedotime1_bin10 <- factor(dat$pseudotime1_bin10, 
                                levels = c("(0,16.8]",
                                           "(16.8,33.5]",
                                           "(33.5,50.3]",
                                           "(50.3,67]",
                                           "(67,83.8]",
                                           "(83.8,101]",
                                           "(101,117]",
                                           "(117,134]",
                                           "(134,151]",
                                           "(151,168]"))
cells <- rownames(dat)
seu$cellID <- colnames(seu)
Idents(seu) <- seu$cellID
lin1 <- subset(seu, idents = cells)
lin1 <- ScaleData(lin1, assay = 'SCT') # scale data relative to only cells included in the dataset 
identical(rownames(dat), colnames(lin1))
lin1 <- AddMetaData(lin1, dat)
Idents(lin1) <- lin1$psuedotime1_bin10
min(table(lin1$psuedotime1_bin10)) # how many cells in smallest bin? Don't want too small...
av.exp <- AverageExpression(lin1, assays = 'SCT', features = sameGenes, return.seurat = TRUE) # create in-silico bulk RNA-seq dataset for each sample
av.exp <- ScaleData(av.exp, assay = 'SCT') # scale data relative to only cells included in the dataset 

bin <- colnames(av.exp)
heatdata <- as.matrix(av.exp[['SCT']]@scale.data)
selGenes1 <- c('CD74', 'CTSS', 'IFI30', 'MARCHF8', 'PYCARD', 'SLA-DMA', 'SLA-DMB', 'SLA-DQA1', 'SLA-DQB1', 'SLA-DRA')
selGenes2 <- c('ACKR3', 'ANXA1', 'CACNB4', 'CACYBP', 'CCR10', 'CCR7', 'CCR9', 'CD274', 'CD44', 'CD74', 'CREB1', 'CRHBP', 'CSF2RB', 'CXCL10', 'CXCR4',
               'CYLD', 'DUSP1', 'EIF2AK1', 'ENSSSCG00000006418', 'FOXO3', 'GAPDH', 'GBP1', 'GSN', 'HSP90AB1',
               'HSPD1', 'ICAM1', 'IFNGR1', 'IFNGR2', 'IL10RB', 'IL11RA', 'IL17RA', 'IL2RG', 'IL7', 'IRF1', 'IRF5', 'ISG15',
               'JARID2', 'KDM5B', 'MNDA', 'MX1', 'MX2', 'NCL', 'NFKB1', 'NFKBIA', 'NKX3-1', 'NLRC5',
               'OAS2', 'P4HB', 'PARP14', 'PELI3', 'PML', 'PTPRC', 'PYCARD', 'RARG', 'RELB', 'RIPK1', 'RNF125',
               'SIGIRR', 'SLA-DQA1', 'SLC25A5', 'SLC27A1', 'SOCS3', 'SRSF7', 'STAT4', 'TIMP1', 
               'TIMP2', 'TLR4', 'TNFRSF18', 'TNFRSF1B', 'TRAF3', 'TRAF3IP2', 'TRAF5', 'TUBA1B', 'UGCG',
               'VIM', 'XCL1', 'YBX3', 'ZBP1', 'ZFP36', 'ZFP36L1', 'ZFP36L2')
selGenes3 <- c('ATAD5', 'BAZ1A', 'CDT1', 'DBF4', 'E2F8', 'FBXO5', 'MCM3', 'MCM4', 'MCM6', 'MCM7', 'NUCKS1', 'PCLAF', 'RBBP6', 'RFC3', 'RPA3', 'RUVBL1', 'SMC3', 'ZBTB38')
selGenes4 <- c('AURKB', 'BUB1', 'BUB3', 'CCNB1', 'CDT1', 'CENPE', 'FBXO5', 'NDC80', 'NEK6', 'RAD21', 'UBE2C')
selGenes <- unique(c(selGenes1, selGenes2, selGenes3, selGenes4))
av.exp <- AverageExpression(lin1, assays = 'SCT', features = selGenes, return.seurat = TRUE) # create in-silico bulk RNA-seq dataset for each sample
av.exp <- ScaleData(av.exp, assay = 'SCT') # scale data relative to only cells included in the dataset 
bin <- colnames(av.exp)
heatdata <- as.matrix(av.exp[['SCT']]@scale.data)
make_italics <- function(x) {
  as.expression(lapply(x, function(y) bquote(italic(.(y)))))
}
ComplexHeatmap::Heatmap(heatdata, cluster_columns = FALSE, col = c('slateblue', 'black', 'yellow2'),
                        row_labels = make_italics(rownames(heatdata)))


### Trajectory 2
IGenes <-  assocRes$gene[
  which(assocRes$padj_lineage2_conditionI < 0.05)
]
JGenes <-  assocRes$gene[
  which(assocRes$padj_lineage2_conditionJ < 0.05)
]

length(IGenes) # How many DE genes for ileum?
#1002
length(JGenes) # How many DE genes for jejunum?
#1213
length(intersect(IGenes, JGenes)) # How many DE genes are shared?
#917
sameGenes <- intersect(IGenes, JGenes)

dat <- data.frame(seu$pseudotime2, seu$celltype, seu$tissue)
max <- max(seu$pseudotime2, na.rm = T)
dat <- dat %>% mutate(psuedotime2_bin10 = cut(seu.pseudotime2, 
                                              breaks=c(0*max, 
                                                       0.1*max,
                                                       0.2*max,
                                                       0.3*max,
                                                       0.4*max,
                                                       0.5*max,
                                                       0.6*max,
                                                       0.7*max,
                                                       0.8*max,
                                                       0.9*max,
                                                       1*max)))
dat <- na.omit(dat, cols = "pseudotime2")
dat$psuedotime2_bin10 <- factor(dat$pseudotime1_bin10, 
                                levels = c('(0,27.3]', '(27.3,54.5]', 
                                           '(54.5,81.8]', '(81.8,109]',
                                           '(109,136]', '(136,164]',
                                           '(164,191]', '(191,218]',
                                           '(218,245]', '(245,273]'))
cells <- rownames(dat)
seu$cellID <- colnames(seu)
Idents(seu) <- seu$cellID
lin2 <- subset(seu, idents = cells)
lin2 <- ScaleData(lin2, assay = 'SCT') # scale data relative to only cells included in the dataset 
identical(rownames(dat), colnames(lin2))
lin2 <- AddMetaData(lin2, dat)
Idents(lin2) <- lin2$psuedotime2_bin10
min(table(lin2$psuedotime2_bin10)) # how many cells in smallest bin? Don't want too small...
av.exp <- AverageExpression(lin2, assays = 'SCT', features = sameGenes, return.seurat = TRUE) # create in-silico bulk RNA-seq dataset for each sample

bin <- colnames(av.exp)
selGenes1 <- c('CD74', 'CTSS', 'IFI30', 'SLA-DMA', 'SLA-DMB', 'SLA-DOA', 'SLA-DOB', 'SLA-DQA1', 'SLA-DQB1', 'SLA-DRA')
selGenes2 <- c('CD19', 'CD22', 'CD38', 'CD79A', 'CD79B', 'CTLA4', 'ENSSSCG00000033721', 'ENSSSCG00000036224', 'ENSSSCG00000040849', 'FCRL3', 'FOXP1', 'GCSAM', 'MNDA', 'MS4A1', 'PRKCB', 'PTPN6', 'PTPRC')
selGenes3 <- c('C4BPB', 'CCR7', 'CD55', 'CD59', 'CFD', 'CXCL10', 'ENSSSCG00000015294', 'ENSSSCG00000028674', 'ENSSSCG00000033721', 'ENSSSCG00000036224', 
              'ENSSSCG00000037775', 'ENSSSCG00000040849', 'GAPDH', 'GPR183', 'JCHAIN', 'LTA', 'NOTCH2', 'NPY', 'PTPN6', 'PTPRC', 'RGCC', 'TRAF3IP2', 'ZP3')
selGenes4 <- c('AICDA', 'CD22', 'CD40', 'FCRL3', 'HMCES', 'MSPD1', 'IL10', 'MZB1', 'PTPRC', 'TRAF3IP2', 'UNG', 'XBP1')
selGenes5 <- c('CALR', 'CD74', 'CLGN', 'DNAJB11', 'DNAJC3', 'ENSSSCG00000029160', 'ERP44', 'FKBP4', 'HSP90AB1', 'HSP90B1', 
              'HSPA8', 'HSPB1', 'HSPD1', 'HSPE1', 'HYOU1', 'P4HB', 'PDIA4', 'PDIA5',
              'PPIB', 'PPIF', 'PRDX4', 'TXNDC5')

#  refer to https://kstreet13.github.io/bioc2020trajectories/articles/workshopTrajectories.html#trajectory-inference-1

seu <- LoadH5Seurat('/home/Jayne.Wiarda/SI_PP_SC_ST/PseudotimeTrajectory/B_traj.h5seurat')
DefaultAssay(seu) <- 'integrated'
genes <- VariableFeatures(seu) # find most variable genes to test...not testing all genes will severely cut down on run time

sce <- readRDS('/home/Jayne.Wiarda/SI_PP_SC_ST/PseudotimeTrajectory/B_traj_sce.rds')

BPPARAM <- BiocParallel::bpparam()
BPPARAM # lists current options
BPPARAM$workers <- 78 

# Before running DE tests, identify an appropriate k value to use. Refer to link at start of script for info
# look for elbow where decrease begins to stabilize in plots

icMat <- evaluateK(counts = as.matrix(assays(sce)$counts),
                   pseudotime = as.matrix(sce$slingshot@assays@data$pseudotime),
                   cellWeights = as.matrix(sce$slingshot@assays@data$weights),
                   conditions = factor(colData(sce)$tissue),
                   nGenes = 300,
                   k = 3:7)
nknots = 5 # use 5 knots for this data

# Do test without consideration of tissue

set.seed(123)
sce <- fitGAM(counts = as.matrix(assays(sce)$counts),
              sds = sce$slingshot,
              conditions = factor(colData(sce)$tissue),
              nknots = nknots,
              genes = genes,
              parallel = TRUE,
              verbose = TRUE,
              BPPARAM = BPPARAM)
mean(rowData(sce)$tradeSeq$converged)

rowData(sce)$assocRes <- associationTest(sce, lineages = TRUE, l2fc = log2(2))
assocRes <- rowData(sce)$assocRes
assocRes$gene <- rownames(assocRes)
assocRes$padj_lineage1_conditionI <- p.adjust(assocRes$pvalue_lineage1_conditionI, "fdr")
assocRes$padj_lineage1_conditionJ <- p.adjust(assocRes$pvalue_lineage1_conditionJ, "fdr")
assocRes$padj_lineage2_conditionI <- p.adjust(assocRes$pvalue_lineage2_conditionI, "fdr")
assocRes$padj_lineage2_conditionJ <- p.adjust(assocRes$pvalue_lineage2_conditionJ, "fdr")
write_xlsx(assocRes, '/home/Jayne.Wiarda/SI_PP_SC_ST/PseudotimeTrajectory/pseudotimeDE.xlsx')
saveRDS(sce, '/home/Jayne.Wiarda/SI_PP_SC_ST/PseudotimeTrajectory/pseudotime_GAM_model.rds')

# Find DE genes between tissues across pseudotime

condRes <- conditionTest(sce, 
                         global = FALSE,
                         lineage = TRUE,
                         pairwise = TRUE,
                         l2fc = log2(2))
condRes$padj_lineage1 <- p.adjust(condRes$pvalue_lineage1, "fdr")
condRes$padj_lineage2 <- p.adjust(condRes$pvalue_lineage2, "fdr")
condRes$gene <- rownames(condRes)
write_xlsx(condRes, '/home/Jayne.Wiarda/SI_PP_SC_ST/PseudotimeTrajectory/pseudotimeDE_condition.xlsx')

# Plot top DE genes (lowest adjusted p) between tissues for trajectory 1:
g1 <- plotSmoothers(sce, assays(sce)$counts, gene = 'CST3', alpha = 1, border = TRUE) + ggtitle('CST3')
g2 <- plotSmoothers(sce, assays(sce)$counts, gene = 'VIM', alpha = 1, border = TRUE) + ggtitle('VIM')
g3 <- plotSmoothers(sce, assays(sce)$counts, gene = 'DUT', alpha = 1, border = TRUE) + ggtitle('DUT')
g4 <- plotSmoothers(sce, assays(sce)$counts, gene = 'RFC3', alpha = 1, border = TRUE) + ggtitle('RFC3')
g1+g2+g3+g4

# Plot top DE genes (lowest adjusted p) between tissues for trajectory 2:
g1 <- plotSmoothers(sce, assays(sce)$counts, gene = 'ENSSSCG00000039111', alpha = 1, border = TRUE) + ggtitle('ENSSSCG00000039111')
g2 <- plotSmoothers(sce, assays(sce)$counts, gene = 'ZP3', alpha = 1, border = TRUE) + ggtitle('ZP3')
g3 <- plotSmoothers(sce, assays(sce)$counts, gene = 'ZNF70', alpha = 1, border = TRUE) + ggtitle('ZNF70')
g4 <- plotSmoothers(sce, assays(sce)$counts, gene = 'ENSSSCG00000034766', alpha = 1, border = TRUE) + ggtitle('ENSSSCG00000034766')
g1+g2+g3+g4

# GO analysis revealed 7 genes that are part of protein folding that were DE across tissues in trajectory 1. Let's plot:
g1 <- plotSmoothers(sce, assays(sce)$counts, gene = 'HSP90AB1', alpha = 1, border = TRUE) + ggtitle('HSP90AB1')
g2 <- plotSmoothers(sce, assays(sce)$counts, gene = 'HSP90AA1', alpha = 1, border = TRUE) + ggtitle('HSP90AA1')
g3 <- plotSmoothers(sce, assays(sce)$counts, gene = 'HSPA8', alpha = 1, border = TRUE) + ggtitle('HSPA8')
g4 <- plotSmoothers(sce, assays(sce)$counts, gene = 'HSP90AA1', alpha = 1, border = TRUE) + ggtitle('HSP90AA1')
g5 <- plotSmoothers(sce, assays(sce)$counts, gene = 'DNAJA1', alpha = 1, border = TRUE) + ggtitle('DNAJA1')
g6 <- plotSmoothers(sce, assays(sce)$counts, gene = 'FKBP4', alpha = 1, border = TRUE) + ggtitle('FKBP4')
g7 <- plotSmoothers(sce, assays(sce)$counts, gene = 'HSPD1', alpha = 1, border = TRUE) + ggtitle('HSPD1')
g1+g2+g3+g4+g5+g6+g7

# GO analysis revealed 4 genes that are part of humoral immunity that were DE across tissues in trajectory 2. Let's plot:
g1 <- plotSmoothers(sce, assays(sce)$counts, gene = 'ENSSSCG00000010044', alpha = 1, border = TRUE) + ggtitle('ENSSSCG00000010044')
g2 <- plotSmoothers(sce, assays(sce)$counts, gene = 'SLA-DRB1', alpha = 1, border = TRUE) + ggtitle('SLA-DRB1')
g3 <- plotSmoothers(sce, assays(sce)$counts, gene = 'ZP3', alpha = 1, border = TRUE) + ggtitle('ZP3')
g4 <- plotSmoothers(sce, assays(sce)$counts, gene = 'ENSSSCG00000040849', alpha = 1, border = TRUE) + ggtitle('ENSSSCG00000040849')
g1+g2+g3+g4

sessionInfo()
#R version 4.2.2 Patched (2022-11-10 r83330)
#Platform: x86_64-pc-linux-gnu (64-bit)
#Running under: Ubuntu 22.04.2 LTS

#Matrix products: default
#BLAS:   /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3
#LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/libopenblasp-r0.3.20.so

#locale:
#  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C               LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8     LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8    LC_PAPER=en_US.UTF-8       LC_NAME=C                  LC_ADDRESS=C              
#[10] LC_TELEPHONE=C             LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       

#attached base packages:
#  [1] stats4    grid      stats     graphics  grDevices utils     datasets  methods   base     

#other attached packages:
#  [1] tidyr_1.2.1                 circlize_0.4.15             ComplexHeatmap_2.12.1       clusterExperiment_2.16.0    data.table_1.14.2           dplyr_1.0.10                ggplot2_3.3.6               SeuratDisk_0.0.0.9020       sp_1.5-0                   
#[10] SeuratObject_4.1.1          Seurat_4.1.1                writexl_1.4.0               tradeSeq_1.10.0             slingshot_2.4.0             TrajectoryUtils_1.4.0       SingleCellExperiment_1.18.0 SummarizedExperiment_1.26.1 Biobase_2.56.0             
#[19] GenomicRanges_1.48.0        GenomeInfoDb_1.32.4         IRanges_2.30.1              S4Vectors_0.34.0            BiocGenerics_0.42.0         MatrixGenerics_1.8.1        matrixStats_0.62.0          princurve_2.1.6            

#loaded via a namespace (and not attached):
#  [1] scattermore_0.8        coda_0.19-4            pkgmaker_0.32.2        bit64_4.0.5            knitr_1.40             irlba_2.3.5            DelayedArray_0.22.0    rpart_4.1.16           KEGGREST_1.36.3        RCurl_1.98-1.8         doParallel_1.0.17     
#[12] generics_0.1.3         ScaledMatrix_1.4.1     cowplot_1.1.1          RSQLite_2.2.17         RANN_2.6.1             future_1.28.0          bit_4.0.4              phylobase_0.8.10       spatstat.data_2.2-0    xml2_1.3.3             httpuv_1.6.6          
#[23] assertthat_0.2.1       viridis_0.6.2          xfun_0.33              hms_1.1.2              evaluate_0.16          promises_1.2.0.1       fansi_1.0.3            progress_1.2.2         igraph_1.3.4           DBI_1.1.3              htmlwidgets_1.5.4     
#[34] spatstat.geom_2.4-0    purrr_0.3.4            ellipsis_0.3.2         RSpectra_0.16-1        ggpubr_0.4.0           backports_1.4.1        annotate_1.74.0        gridBase_0.4-7         locfdr_1.1-8           deldir_1.0-6           vctrs_0.4.1           
#[45] ggalluvial_0.12.3      Cairo_1.6-0            ROCR_1.0-11            abind_1.4-5            cachem_1.0.6           withr_2.5.0            progressr_0.11.0       sctransform_0.3.4      sna_2.7                prettyunits_1.1.1      goftest_1.2-3         
#[56] softImpute_1.4-1       svglite_2.1.0          cluster_2.1.4          ape_5.6-2              lazyeval_0.2.2         crayon_1.5.1           genefilter_1.78.0      hdf5r_1.3.5            edgeR_3.38.4           pkgconfig_2.0.3        labeling_0.4.2        
#[67] nlme_3.1-159           rlang_1.0.6            globals_0.16.1         lifecycle_1.0.2        miniUI_0.1.1.1         registry_0.5-1         rsvd_1.0.5             cellranger_1.1.0       polyclip_1.10-0        lmtest_0.9-40          rngtools_1.5.2        
#[78] Matrix_1.5-1           carData_3.0-5          Rhdf5lib_1.18.2        zoo_1.8-10             ggridges_0.5.3         GlobalOptions_0.1.2    png_0.1-7              viridisLite_0.4.1      rjson_0.2.21           bitops_1.0-7           rhdf5filters_1.8.0    
#[89] rncl_0.8.6             KernSmooth_2.23-20     ggnetwork_0.5.10       Biostrings_2.64.1      blob_1.2.3             shape_1.4.6            stringr_1.4.1          zinbwave_1.18.0        parallelly_1.32.1      spatstat.random_2.2-0  rstatix_0.7.0         
#[100] ggsignif_0.6.3         beachmat_2.12.0        scales_1.2.1           memoise_2.0.1          magrittr_2.0.3         plyr_1.8.7             ica_1.0-3              howmany_0.3-1          zlibbioc_1.42.0        compiler_4.2.2         RColorBrewer_1.1-3    
#[111] clue_0.3-61            fitdistrplus_1.1-8     cli_3.4.0              ade4_1.7-19            XVector_0.36.0         listenv_0.8.0          patchwork_1.1.2        pbapply_1.5-0          MASS_7.3-58.1          mgcv_1.8-40            tidyselect_1.1.2      
#[122] stringi_1.7.8          yaml_2.3.5             BiocSingular_1.12.0    locfit_1.5-9.6         ggrepel_0.9.1          tools_4.2.2            future.apply_1.9.1     parallel_4.2.2         rstudioapi_0.14        uuid_1.1-0             foreach_1.5.2         
#[133] RNeXML_2.4.7           gridExtra_2.3          farver_2.1.1           Rtsne_0.16             digest_0.6.29          rgeos_0.5-9            FNN_1.1.3.1            shiny_1.7.2            Rcpp_1.0.9             car_3.1-0              broom_1.0.1           
#[144] later_1.3.0            RcppAnnoy_0.0.19       AnnotationDbi_1.58.0   httr_1.4.4             kernlab_0.9-31         colorspace_2.0-3       XML_3.99-0.10          tensor_1.5             reticulate_1.26        splines_4.2.2          uwot_0.1.14           
#[155] spatstat.utils_2.3-1   plotly_4.10.0          systemfonts_1.0.4      xtable_1.8-4           jsonlite_1.8.0         R6_2.5.1               pillar_1.8.1           htmltools_0.5.3        mime_0.12              NMF_0.24.0             glue_1.6.2            
#[166] fastmap_1.1.0          BiocParallel_1.30.3    BiocNeighbors_1.14.0   codetools_0.2-18       utf8_1.2.2             lattice_0.20-45        spatstat.sparse_2.1-1  tibble_3.1.8           network_1.17.2         leiden_0.4.3           survival_3.4-0        
#[177] limma_3.52.3           rmarkdown_2.16         statnet.common_4.7.0   munsell_0.5.0          rhdf5_2.40.0           GetoptLong_1.0.5       GenomeInfoDbData_1.2.8 iterators_1.0.14       HDF5Array_1.24.2       reshape2_1.4.4         gtable_0.3.1          
#[188] spatstat.core_2.4-4   