library(tidyverse)
library(biomaRt)  # used to map GOterms to Ensembl IDs
library(topGO)
library(readxl)
library(ggnewscale)
library(slingshot)
library(SingleCellExperiment)
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

seu <- LoadH5Seurat('/home/Jayne.Wiarda/SI_PP_SC_ST/PseudotimeTrajectory/B_traj.h5seurat')

GO_background_generator <- 
  function(GO_term_universe, seurat_obj, output_path){
    # gene list is enriched genes, 
    # GO_term_universe is already prepared, needs to be filtered
    # seurat obj will be used to filter GO_term_universe
    # browser()
    # need to make sure all genes in seurat object are detected
    counts <- seurat_obj@assays$SCT@counts
    if (!all(rowSums(counts) > 0)){
      print('some genes not detected, filtering genes with 0 counts')
      counts <- counts[rowSums(counts) > 0,]
    }
    
    all_GO <- read_tsv('/home/Jayne.Wiarda/SI_PP_SC_ST/GeneOntology/AllCells_gene_to_GO.tsv')
    new_universe <-  all_GO %>% 
      filter(FinalList %in% rownames(counts)) %>% 
      write_tsv(output_path)
    print(paste('number of genes in all universe',nrow(all_GO)))
    print(paste('number of genes in new universe',nrow(new_universe)))
    
  }

# wrapper function that makes using topGO easier (I think)

# you can also access this function by loading my R package 'funfuns'
# that way you don't need to read in this function every time
# remotes::install_github('jtrachsel/funfuns')
# library(funfuns)


topGO_wrapper <- function(myInterestingGenes, #vector
                          mapping_file,       # two column file
                          ont='BP',
                          algor = 'elim',
                          statistic='Fisher',
                          nodeSize=10, 
                          return_GOdata=F){
  
  require(topGO)
  # browser()
  geneID2GO <- readMappings(mapping_file)
  geneNames <- names(geneID2GO)
  
  # Get the list of genes of interest
  geneList <- factor(as.integer(geneNames %in% myInterestingGenes))
  names(geneList) <- geneNames
  
  #initialize topGOdata object
  GOdata <- new("topGOdata", ontology = ont, allGenes = geneList,
                annot = annFUN.gene2GO, gene2GO = geneID2GO,
                nodeSize=nodeSize)
  
  # contstruct a tibble that maps interesting genes to GO terms
  # this can add a decent amount of time to the function call...
  interesting_genes_in_GOs <- 
    genesInTerm(GOdata) %>% 
    enframe(name='GO.ID', 
            value='genes_in_GO') %>% 
    mutate(involved_genes=map(.x=genes_in_GO, ~ .x[.x %in% myInterestingGenes]), 
           involved_genes=map_chr(.x=involved_genes, ~paste(.x, collapse = '_'))) %>% 
    dplyr::select(GO.ID, involved_genes) 
  
  # Run topGO test
  resultTopGO.elim <- runTest(GOdata, algorithm = algor, statistic = statistic )
  allRes <- GenTable(GOdata, pval = resultTopGO.elim,
                     orderBy = "pval",
                     topNodes = length(GOdata@graph@nodes), #include all nodes
                     numChar=1000)
  
  # clean up results and add in extra info
  allRes <- allRes %>%
    mutate(ont=ifelse(ont=='BP', 'Biological Process',
                      ifelse(ont=='MF', 'Molecular Function', "Cellular Component"))) %>%
    mutate(GO_aspect = ont,
           algorithm = algor,
           statistic = statistic) %>%
    dplyr::select(-ont) %>% 
    left_join(interesting_genes_in_GOs)
  
  if (return_GOdata == TRUE){
    return(list(allRes, GOdata))
  } else {
    return(allRes)
  }
  
  
}


#### Get ensemble gene_IDs for pig genome ###
# This section uses the internet to map these IDs, so sometimes it is very slow
# in fact it is not working right now.  I think the USDA's internet tubes get all 
# blocked up during business hours, what with all the video calls

# select mart and data set
bm <- useMart("ensembl")
bm <- useDataset("sscrofa_gene_ensembl", mart=bm)

# Get ensembl gene ids and GO terms
EG2GO <- getBM(mart=bm, attributes=c('ensembl_gene_id','go_id'))

# examine result
head(EG2GO,15)

# Remove blank entries
EG2GO <- EG2GO[EG2GO$go_id != '',]

### format mapping file for use with topGO wrapper function
EG2GO <- 
  EG2GO %>%
  group_by(ensembl_gene_id) %>% 
  summarise(GO=paste(go_id, sep = ' ', collapse = ',')) %>% 
  transmute(EnsemblID=ensembl_gene_id, 
            GO=GO) 

# scRNA-seq all cells ----

# to map gene 'FinalList' annotation to Ensembl gene_IDs
gene_IDs <- read_excel("/home/Jayne.Wiarda/SI_PP_SC_ST/GeneAnnotationFiles/UpdatedGeneNameListForSus97GTF_06302021_JEW_SKS.xlsx")
gene_IDs$EnsemblID <- gene_IDs$ENSID # make identical to column name in EG2GO

background <- read_excel('/home/Jayne.Wiarda/SI_PP_SC_ST/DGE/CellTypeDGE_AllSamples_BackgroundGeneList.xlsx')
detected_genes <- background$FinalAnnot

# all GO terms detected in this dataset

GO_gene_universe <- 
  gene_IDs %>% 
  filter(FinalList %in% detected_genes) %>% 
  left_join(EG2GO) %>% 
  filter(!is.na(GO))

GO_gene_universe %>%
  dplyr::select(FinalList, GO) %>% 
  write_tsv('/home/Jayne.Wiarda/SI_PP_SC_ST/GeneOntology/AllCells_gene_to_GO.tsv')

# Trajectory B cells only ----
# Modify background generator to take from variable features of integrated assay (since this is what we used for DGE here) rather than SCT:
GO_background_generator <- 
  function(GO_term_universe, seurat_obj, output_path){
    # gene list is enriched genes, 
    # GO_term_universe is already prepared, needs to be filtered
    # seurat obj will be used to filter GO_term_universe
    # browser()
    # need to make sure all genes in seurat object are detected
    counts <- VariableFeatures(seu)
    all_GO <- read_tsv('/home/Jayne.Wiarda/SI_PP_SC_ST/GeneOntology/AllCells_gene_to_GO.tsv')
    new_universe <-  all_GO %>% 
      filter(FinalList %in% genes) %>% 
      write_tsv(output_path)
    print(paste('number of genes in all universe',nrow(all_GO)))
    print(paste('number of genes in new universe',nrow(new_universe)))
  }

seu <- LoadH5Seurat('/home/Jayne.Wiarda/SI_PP_SC_ST/PseudotimeTrajectory/B_traj.h5seurat')
DefaultAssay(seu) <- 'integrated'

GO_background_generator(GO_term_universe = '/home/Jayne.Wiarda/SI_PP_SC_ST/GeneOntology/AllCells_gene_to_GO.tsv',
                        seurat_obj = seu, 
                        output_path = '/home/Jayne.Wiarda/SI_PP_SC_ST/GeneOntology/BCellTrajectory_gene_to_GO.tsv')

## Conserved DGE in ileum and jejunum across pseudotime
BtrajDGE <- read_xlsx('/home/Jayne.Wiarda/SI_PP_SC_ST/PseudotimeTrajectory/pseudotimeDE.xlsx')

IGenesTraj1 <-  BtrajDGE$gene[
  which(BtrajDGE$padj_lineage1_conditionI < 0.05)
]
JGenesTraj1 <-  BtrajDGE$gene[
  which(BtrajDGE$padj_lineage1_conditionJ < 0.05)
]
sameGenesTraj1 <- data.frame(intersect(IGenesTraj1, JGenesTraj1))
colnames(sameGenesTraj1) <- 'gene'
sameGenesTraj1$cluster <- rep('Trajectory1', nrow(sameGenesTraj1))

IGenesTraj2 <-  BtrajDGE$gene[
  which(BtrajDGE$padj_lineage2_conditionI < 0.05)
]
JGenesTraj2 <-  BtrajDGE$gene[
  which(BtrajDGE$padj_lineage2_conditionJ < 0.05)
]
sameGenesTraj2 <- data.frame(intersect(IGenesTraj2, JGenesTraj2))
colnames(sameGenesTraj2) <- 'gene'
sameGenesTraj2$cluster <- rep('Trajectory2', nrow(sameGenesTraj2))

sameGenes <- rbind(sameGenesTraj1, sameGenesTraj2)

Btraj_results <- 
  sameGenes %>% 
  group_by(cluster) %>% 
  nest() %>% 
  mutate(enriched_genes=map(data, ~.x %>% pull(gene))) %>% 
  mutate(GO_results=
           map(enriched_genes, ~topGO_wrapper(myInterestingGenes = .x,
                                              mapping_file = '/home/Jayne.Wiarda/SI_PP_SC_ST/GeneOntology/BCellTrajectory_gene_to_GO.tsv'))) %>% 
  mutate(filt_results=purrr::map(.x = GO_results, .f = ~filter(.x, pval < 0.05))) %>% 
  dplyr::select(cluster, filt_results) %>% 
  unnest(cols=filt_results)

Btraj_results$Fold_enrichment <- Btraj_results$Significant/Btraj_results$Expected # add in stat for fold enrichment
Btraj_results <- subset(Btraj_results, Significant > 1) # remove GO terms with only one enriched genedim

Btraj_results %>% write_tsv('/home/Jayne.Wiarda/SI_PP_SC_ST/GeneOntology/Btraj_ConservedDGE_GOresults.tsv')

## Plot results in heatmaps

### Trajectory 1
assocRes <- BtrajDGE
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

# We selected genes from 4 biological processes to plot:
# 1) antigen processing and presentation of peptide antigen via MHC-II; 2) response to cytokine; 3) DNA replication; 4) mitotic chromatid separation
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

col_fun<-colorRamp2(c(-2, 0, 2), c('slateblue', 'black', 'yellow2'))
ComplexHeatmap::Heatmap(heatdata, cluster_columns = FALSE, col = col_fun,
                        row_labels = make_italics(rownames(heatdata)))

subdata <- heatdata[rownames(heatdata) %in% selGenes1, ] 
h1 <- ComplexHeatmap::Heatmap(subdata, cluster_columns = FALSE, col = col_fun,
                              row_labels = make_italics(rownames(subdata)), name = 'Relative gene expression')

subdata <- heatdata[rownames(heatdata) %in% selGenes2, ] 
h2 <- ComplexHeatmap::Heatmap(subdata, cluster_columns = FALSE, col = col_fun,
                              row_labels = make_italics(rownames(subdata)), name = 'Relative gene expression')

subdata <- heatdata[rownames(heatdata) %in% selGenes3, ] 
h3 <- ComplexHeatmap::Heatmap(subdata, cluster_columns = FALSE, col = col_fun,
                              row_labels = make_italics(rownames(subdata)), name = 'Relative gene expression')

subdata <- heatdata[rownames(heatdata) %in% selGenes4, ] 
h4 <- ComplexHeatmap::Heatmap(subdata, cluster_columns = FALSE, col = col_fun,
                              row_labels = make_italics(rownames(subdata)), name = 'Relative gene expression')

dat <- data.frame(lin1$S.Score, lin1$G2M.Score, lin1$celltype, lin1$psuedotime1_bin10)
colnames(dat) <- c('S.score', 'G2M.score', 'celltype', 'bin')
sScore <- dat %>%
  group_by(bin) %>%
  summarise_at(vars(S.score), list(name = mean))
colnames(sScore) <- c('bin', 'S.score')
G2MScore <- dat %>%
  group_by(bin) %>%
  summarise_at(vars(G2M.score), list(name = mean))
colnames(G2MScore) <- c('bin', 'G2M.score')
scores <- cbind(sScore, G2MScore)
rownames(scores) <- scores$bin
scores <- scores[,c('S.score', 'G2M.score')]
scores <- t(scores)
col_fun<-colorRamp2(c(min(seu$G2M.Score, seu$S.Score),max(seu$G2M.Score, seu$S.Score)), c('grey90', 'darkgreen'))
h5 <- ComplexHeatmap::Heatmap(scores, cluster_columns = FALSE, cluster_rows = FALSE, col = col_fun, name = 'Cell cycle score')

tab <- prop.table(table(dat$celltype, dat$bin), margin = 2)
col_fun<-colorRamp2(c(0,1), c('grey90', 'magenta4'))
h6 <- ComplexHeatmap::Heatmap(tab, cluster_columns = FALSE, cluster_rows = FALSE, col = col_fun, name = 'Cell type frequency')

ht_list <- h6 %v% h5 %v% h1 %v% h2 %v% h3 %v% h4
draw(ht_list) 

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

selGenes1 <- c('CD74', 'CTSS', 'IFI30', 'SLA-DMA', 'SLA-DMB', 'SLA-DOA', 'SLA-DOB', 'SLA-DQA1', 'SLA-DQB1', 'SLA-DRA')
selGenes2 <- c('CD19', 'CD22', 'CD38', 'CD79A', 'CD79B', 'CTLA4', 'ENSSSCG00000033721', 'ENSSSCG00000036224', 'ENSSSCG00000040849', 'FCRL3', 'FOXP1', 'GCSAM', 'MNDA', 'MS4A1', 'PRKCB', 'PTPN6', 'PTPRC')
selGenes3 <- c('C4BPB', 'CCR7', 'CD55', 'CD59', 'CFD', 'CXCL10', 'ENSSSCG00000015294', 'ENSSSCG00000028674', 'ENSSSCG00000033721', 'ENSSSCG00000036224', 
               'ENSSSCG00000037775', 'ENSSSCG00000040849', 'GAPDH', 'GPR183', 'JCHAIN', 'LTA', 'NOTCH2', 'NPY', 'PTPN6', 'PTPRC', 'RGCC', 'TRAF3IP2', 'ZP3')
selGenes4 <- c('AICDA', 'CD22', 'CD40', 'FCRL3', 'HMCES', 'MSPD1', 'IL10', 'MZB1', 'PTPRC', 'TRAF3IP2', 'UNG', 'XBP1')
selGenes5 <- c('CALR', 'CD74', 'CLGN', 'DNAJB11', 'DNAJC3', 'ENSSSCG00000029160', 'ERP44', 'FKBP4', 'HSP90AB1', 'HSP90B1', 
               'HSPA8', 'HSPB1', 'HSPD1', 'HSPE1', 'HYOU1', 'P4HB', 'PDIA4', 'PDIA5',
               'PPIB', 'PPIF', 'PRDX4', 'TXNDC5')
selGenes <- unique(c(selGenes1, selGenes2, selGenes3, selGenes4, selGenes5))
av.exp <- AverageExpression(lin2, assays = 'SCT', features = selGenes, return.seurat = TRUE) # create in-silico bulk RNA-seq dataset for each sample
av.exp <- ScaleData(av.exp, assay = 'SCT') # scale data relative to only cells included in the dataset 
bin <- colnames(av.exp)
heatdata <- as.matrix(av.exp[['SCT']]@scale.data)
make_italics <- function(x) {
  as.expression(lapply(x, function(y) bquote(italic(.(y)))))
}
col_fun<-colorRamp2(c(-2, 0, 2), c('slateblue', 'black', 'yellow2'))
ComplexHeatmap::Heatmap(heatdata, cluster_columns = FALSE, col = col_fun,
                        row_labels = make_italics(rownames(heatdata)))

subdata <- heatdata[rownames(heatdata) %in% selGenes1, ] 
h1 <- ComplexHeatmap::Heatmap(subdata, cluster_columns = FALSE, col = col_fun,
                              row_labels = make_italics(rownames(subdata)), name = 'Relative gene expression')

subdata <- heatdata[rownames(heatdata) %in% selGenes2, ] 
h2 <- ComplexHeatmap::Heatmap(subdata, cluster_columns = FALSE, col = col_fun,
                              row_labels = make_italics(rownames(subdata)), name = 'Relative gene expression')

subdata <- heatdata[rownames(heatdata) %in% selGenes3, ] 
h3 <- ComplexHeatmap::Heatmap(subdata, cluster_columns = FALSE, col = col_fun,
                              row_labels = make_italics(rownames(subdata)), name = 'Relative gene expression')

subdata <- heatdata[rownames(heatdata) %in% selGenes4, ] 
h4 <- ComplexHeatmap::Heatmap(subdata, cluster_columns = FALSE, col = col_fun,
                              row_labels = make_italics(rownames(subdata)), name = 'Relative gene expression')

subdata <- heatdata[rownames(heatdata) %in% selGenes5, ] 
h5 <- ComplexHeatmap::Heatmap(subdata, cluster_columns = FALSE, col = col_fun,
                              row_labels = make_italics(rownames(subdata)), name = 'Relative gene expression')

dat <- data.frame(lin2$S.Score, lin2$G2M.Score, lin2$celltype, lin2$psuedotime2_bin10)
colnames(dat) <- c('S.score', 'G2M.score', 'celltype', 'bin')
sScore <- dat %>%
  group_by(bin) %>%
  summarise_at(vars(S.score), list(name = mean))
colnames(sScore) <- c('bin', 'S.score')
G2MScore <- dat %>%
  group_by(bin) %>%
  summarise_at(vars(G2M.score), list(name = mean))
colnames(G2MScore) <- c('bin', 'G2M.score')
scores <- cbind(sScore, G2MScore)
rownames(scores) <- scores$bin
scores <- scores[,c('S.score', 'G2M.score')]
scores <- t(scores)
col_fun<-colorRamp2(c(min(seu$G2M.Score, seu$S.Score),max(seu$G2M.Score, seu$S.Score)), c('grey90', 'darkgreen'))
h6 <- ComplexHeatmap::Heatmap(scores, cluster_columns = FALSE, cluster_rows = FALSE, col = col_fun, name = 'Cell cycle score')

tab <- prop.table(table(dat$celltype, dat$bin), margin = 2)
col_fun<-colorRamp2(c(0,1), c('grey90', 'magenta4'))
h7 <- ComplexHeatmap::Heatmap(tab, cluster_columns = FALSE, cluster_rows = FALSE, col = col_fun, name = 'Cell type frequency')

ht_list <- h7 %v% h6 %v% h1 %v% h2 %v% h3 %v% h4 %v% h5
draw(ht_list) 



## Tissue-specific DGE
BtrajTissueDGE <- read_xlsx('/home/Jayne.Wiarda/SI_PP_SC_ST/PseudotimeTrajectory/pseudotimeDE_condition.xlsx')

GenesTraj1 <-  BtrajTissueDGE$gene[
  which(BtrajTissueDGE$padj_lineage1 < 0.05)
]
GenesTraj1 <- data.frame(GenesTraj1)
colnames(GenesTraj1) <- 'gene'
GenesTraj1$cluster <- rep('Trajectory1', nrow(GenesTraj1))

GenesTraj2 <-  BtrajTissueDGE$gene[
  which(BtrajTissueDGE$padj_lineage2 < 0.05)
]
GenesTraj2 <- data.frame(GenesTraj2)
colnames(GenesTraj2) <- 'gene'
GenesTraj2$cluster <- rep('Trajectory2', nrow(GenesTraj2))

DEGenes <- rbind(GenesTraj1, GenesTraj2)

Btraj_results <- 
  DEGenes %>% 
  group_by(cluster) %>% 
  nest() %>% 
  mutate(enriched_genes=map(data, ~.x %>% pull(gene))) %>% 
  mutate(GO_results=
           map(enriched_genes, ~topGO_wrapper(myInterestingGenes = .x,
                                              mapping_file = '/home/Jayne.Wiarda/SI_PP_SC_ST/GeneOntology/BCellTrajectory_gene_to_GO.tsv'))) %>% 
  mutate(filt_results=purrr::map(.x = GO_results, .f = ~filter(.x, pval < 0.05))) %>% 
  dplyr::select(cluster, filt_results) %>% 
  unnest(cols=filt_results)

Btraj_results$Fold_enrichment <- Btraj_results$Significant/Btraj_results$Expected # add in stat for fold enrichment

Btraj_results %>% write_tsv('/home/Jayne.Wiarda/SI_PP_SC_ST/GeneOntology/Btraj_TissueSpecificDGE_GOresults.tsv')

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
#  [1] circlize_0.4.15             ComplexHeatmap_2.12.1       clusterExperiment_2.16.0    data.table_1.14.2           SeuratDisk_0.0.0.9020       sp_1.5-0                    SeuratObject_4.1.1          Seurat_4.1.1                writexl_1.4.0              
#[10] slingshot_2.4.0             TrajectoryUtils_1.4.0       SingleCellExperiment_1.18.0 SummarizedExperiment_1.26.1 GenomicRanges_1.48.0        GenomeInfoDb_1.32.4         MatrixGenerics_1.8.1        matrixStats_0.62.0          princurve_2.1.6            
#[19] ggnewscale_0.4.8            readxl_1.4.1                topGO_2.48.0                SparseM_1.81                GO.db_3.15.0                AnnotationDbi_1.58.0        IRanges_2.30.1              S4Vectors_0.34.0            Biobase_2.56.0             
#[28] graph_1.74.0                BiocGenerics_0.42.0         biomaRt_2.52.0              forcats_0.5.2               stringr_1.4.1               dplyr_1.0.10                purrr_0.3.4                 readr_2.1.2                 tidyr_1.2.1                
#[37] tibble_3.1.8                ggplot2_3.3.6               tidyverse_1.3.2            

#loaded via a namespace (and not attached):
#  [1] statnet.common_4.7.0   rsvd_1.0.5             ica_1.0-3              zinbwave_1.18.0        svglite_2.1.0          foreach_1.5.2          lmtest_0.9-40          crayon_1.5.1           spatstat.core_2.4-4    MASS_7.3-58.1          rhdf5filters_1.8.0    
#[12] nlme_3.1-159           backports_1.4.1        reprex_2.0.2           rlang_1.0.6            XVector_0.36.0         ROCR_1.0-11            irlba_2.3.5            limma_3.52.3           phylobase_0.8.10       filelock_1.0.2         BiocParallel_1.30.3   
#[23] rjson_0.2.21           bit64_4.0.5            glue_1.6.2             rngtools_1.5.2         sctransform_0.3.4      parallel_4.2.2         spatstat.sparse_2.1-1  spatstat.geom_2.4-0    haven_2.5.1            tidyselect_1.1.2       fitdistrplus_1.1-8    
#[34] XML_3.99-0.10          zoo_1.8-10             ggpubr_0.4.0           xtable_1.8-4           ggnetwork_0.5.10       magrittr_2.0.3         evaluate_0.16          cli_3.4.0              zlibbioc_1.42.0        rstudioapi_0.14        miniUI_0.1.1.1        
#[45] rpart_4.1.16           locfdr_1.1-8           shiny_1.7.2            BiocSingular_1.12.0    xfun_0.33              clue_0.3-61            cluster_2.1.4          KEGGREST_1.36.3        ggrepel_0.9.1          ape_5.6-2              listenv_0.8.0         
#[56] Biostrings_2.64.1      png_0.1-7              future_1.28.0          withr_2.5.0            bitops_1.0-7           plyr_1.8.7             cellranger_1.1.0       coda_0.19-4            pillar_1.8.1           GlobalOptions_0.1.2    cachem_1.0.6          
#[67] fs_1.5.2               kernlab_0.9-31         hdf5r_1.3.5            GetoptLong_1.0.5       vctrs_0.4.1            ellipsis_0.3.2         generics_0.1.3         NMF_0.24.0             tools_4.2.2            rncl_0.8.6             munsell_0.5.0         
#[78] DelayedArray_0.22.0    fastmap_1.1.0          compiler_4.2.2         abind_1.4-5            httpuv_1.6.6           pkgmaker_0.32.2        plotly_4.10.0          rgeos_0.5-9            GenomeInfoDbData_1.2.8 gridExtra_2.3          edgeR_3.38.4          
#[89] lattice_0.20-45        deldir_1.0-6           utf8_1.2.2             later_1.3.0            BiocFileCache_2.4.0    jsonlite_1.8.0         scales_1.2.1           ScaledMatrix_1.4.1     pbapply_1.5-0          carData_3.0-5          genefilter_1.78.0     
#[100] lazyeval_0.2.2         promises_1.2.0.1       car_3.1-0              doParallel_1.0.17      goftest_1.2-3          spatstat.utils_2.3-1   reticulate_1.26        sna_2.7                rmarkdown_2.16         cowplot_1.1.1          Rtsne_0.16            
#[111] softImpute_1.4-1       uwot_0.1.14            igraph_1.3.4           HDF5Array_1.24.2       survival_3.4-0         yaml_2.3.5             systemfonts_1.0.4      htmltools_0.5.3        memoise_2.0.1          locfit_1.5-9.6         viridisLite_0.4.1     
#[122] digest_0.6.29          assertthat_0.2.1       rappdirs_0.3.3         mime_0.12              registry_0.5-1         RSQLite_2.2.17         future.apply_1.9.1     blob_1.2.3             RNeXML_2.4.7           splines_4.2.2          labeling_0.4.2        
#[133] Rhdf5lib_1.18.2        Cairo_1.6-0            googledrive_2.0.0      RCurl_1.98-1.8         broom_1.0.1            hms_1.1.2              modelr_0.1.10          rhdf5_2.40.0           colorspace_2.0-3       shape_1.4.6            Rcpp_1.0.9            
#[144] RANN_2.6.1             fansi_1.0.3            tzdb_0.3.0             parallelly_1.32.1      R6_2.5.1               ggridges_0.5.3         lifecycle_1.0.2        curl_4.3.2             ggsignif_0.6.3         googlesheets4_1.0.1    leiden_0.4.3          
#[155] Matrix_1.5-1           howmany_0.3-1          RcppAnnoy_0.0.19       RColorBrewer_1.1-3     iterators_1.0.14       htmlwidgets_1.5.4      beachmat_2.12.0        polyclip_1.10-0        network_1.17.2         timechange_0.2.0       rvest_1.0.3           
#[166] mgcv_1.8-40            globals_0.16.1         patchwork_1.1.2        spatstat.random_2.2-0  progressr_0.11.0       codetools_0.2-18       lubridate_1.9.0        FNN_1.1.3.1            prettyunits_1.1.1      dbplyr_2.2.1           gridBase_0.4-7        
#[177] RSpectra_0.16-1        gtable_0.3.1           DBI_1.1.3              ggalluvial_0.12.3      tensor_1.5             httr_1.4.4             KernSmooth_2.23-20     stringi_1.7.8          progress_1.2.2         reshape2_1.4.4         farver_2.1.1          
#[188] uuid_1.1-0             annotate_1.74.0        viridis_0.6.2          xml2_1.3.3             BiocNeighbors_1.14.0   ade4_1.7-19            scattermore_0.8        bit_4.0.4              spatstat.data_2.2-0    pkgconfig_2.0.3        gargle_1.2.1          
#[199] rstatix_0.7.0          knitr_1.40 