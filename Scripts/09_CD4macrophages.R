### Load required software packages

library(Seurat)
library(SeuratObject)
library(ggplot2)
library(SeuratDisk)
library(dplyr)
library(readxl)
library(writexl)
library(readxl)
library(CellChat)
library(tidyverse)
library(biomaRt)  # used to map GOterms to Ensembl IDs
library(topGO)
library(ggrepel)
#library(scales)
#library(viridis)
library(ggnewscale)

## Load data
#Load in h5 Seurat object of our ileum/jejunum data:
seu <- LoadH5Seurat('/home/Jayne.Wiarda/SI_PP_SC_ST/AnnotatedSeurat/AllSamples_annotated.h5seurat')
Idents(seu) <- seu$celltype
DefaultAssay(seu) <- 'SCT'

dat <- data.frame(colnames(seu), seu$celltype, seu[['SCT']]@data['CD4',])
colnames(dat) <- c('barcode', 'celltype', 'CD4')
dat$celltype <- as.character(dat$celltype)

CD4pos <- dat$barcode[dat$celltype == 'Macrophages' & dat$CD4 > 0]
CD4neg <- dat$barcode[dat$celltype == 'Macrophages' & dat$CD4 == 0]

dat$celltype[dat$barcode %in% CD4pos] <- 'CD4pos macrophage'
dat$celltype[dat$barcode %in% CD4neg] <- 'CD4neg macrophage'

seu$celltype2 <- dat$celltype
Idents(seu) <- seu$celltype2

# DGE ----

seu <- PrepSCTFindMarkers(seu)
de <- FindMarkers(seu, 
                  ident.1 = 'CD4pos macrophage', 
                  ident.2 = 'CD4neg macrophage')
de$Gene <- rownames(de)
de <- subset(de, p_val_adj < 0.05)
de$cluster <- ifelse(de$avg_log2FC > 0, 'CD4pos macrophage', 'CD4neg macrophage')
write_xlsx(de, '/home/Jayne.Wiarda/SI_PP_SC_ST/DGE/Macrophages_CD4_DGE.xlsx')

avg <- as.data.frame(log1p(AverageExpression(seu, verbose = FALSE)$SCT))
avg$gene <- rownames(avg)

genes.to.label = c(de$Gene)
lab <- filter(avg, gene %in% c(genes.to.label))

ggplot(avg, aes(`CD4neg macrophage`, `CD4pos macrophage`)) + 
  geom_point() + 
  ggtitle("Macrophages - CD4+ vs. CD4- DEGs") +
  theme_bw() +
  geom_point(data = lab,
             aes(`CD4neg macrophage`, `CD4pos macrophage`),
             color = 'red') +
  theme(axis.text=element_text(size=25),
        axis.title=element_text(size=25)) +
  geom_text_repel(data = lab, aes(label = gene)) # label a subsete of DE genes. #/which genes labeled depends on overlap in graph.

ggplot(avg, aes(`CD4neg macrophage`, `CD4pos macrophage`)) + 
  geom_point() + 
  ggtitle("Macrophages - CD4+ vs. CD4- DEGs") +
  theme_bw() +
  geom_point(data = lab,
             aes(`CD4neg macrophage`, `CD4pos macrophage`),
             color = 'red') +
  theme(axis.text=element_text(size=25),
         axis.title=element_text(size=25))


# GO ----

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

de$cluster <- ifelse(de$avg_log2FC > 0, 'CD4pos macrophage', 'CD4neg macrophage')
results <- 
  de %>% 
  group_by(cluster) %>% 
  nest() %>% 
  mutate(enriched_genes=map(data, ~.x %>% pull(Gene))) %>% 
  mutate(GO_results=
           map(enriched_genes, ~topGO_wrapper(myInterestingGenes = .x,
                                              mapping_file = '/home/Jayne.Wiarda/SI_PP_SC_ST/GeneOntology/AllCells_gene_to_GO.tsv'))) %>% 
  mutate(filt_results=purrr::map(.x = GO_results, .f = ~filter(.x, pval < 0.05))) %>% 
  dplyr::select(cluster, filt_results) %>% 
  unnest(cols=filt_results)

results$Fold_enrichment <- results$Significant/results$Expected # add in stat for fold enrichment
results <- subset(results, Significant > 1) # remove GO terms with only one enriched genedim

results %>% write_tsv('/home/Jayne.Wiarda/SI_PP_SC_ST/GeneOntology/CD4macrophage_GOresults.tsv')

# Create dot plot of selected GO processes for spatial regions
terms <- c('phagocytosis, engulfment',
           'complement activation, classical pathway',
           'positive regulation of catalytic activity',
           'positive regulation of leukocyte mediated cytotoxicity',
           'positive regulation of inflammatory response',
           'cell surface receptor signaling pathway',
           'positive regulation of macrophage activation',
           'acute inflammatory response to antigenic stimulus',
           'positive regulation of translation', 
           'positive regulation of cell motility',
           'positive regulation of locomotion')
results <- results[results$Term %in% terms, ]
results$Term <- factor(results$Term,
                                  levels = rev(c('positive regulation of translation',
                                                 'positive regulation of cell motility',
                                                 'positive regulation of locomotion',
                                                 'positive regulation of macrophage activation',
                                                 'complement activation, classical pathway',
                                                 'phagocytosis, engulfment',
                                                 'acute inflammatory response to antigenic stimulus',
                                                 'positive regulation of inflammatory response',
                                                 'positive regulation of leukocyte mediated cytotoxicity',
                                                 'positive regulation of catalytic activity',
                                                 'cell surface receptor signaling pathway')))
results$pval <- as.numeric(results$pval)
ggplot(results) +
  geom_point(aes(x = cluster, 
                 y = Term,
                 size = Fold_enrichment,
                 color = pval)) +
  theme_bw() + 
  scale_colour_gradient2(low="darkslateblue", mid="darkseagreen", high="khaki3", 
                         limits=c(0, 0.05),
                         midpoint = .025) +
  theme(axis.text.x = element_text(angle=90, vjust=.5, hjust=1))

# Cell-cell ixn ----

seu$tissue <- substr(seu$orig.ident, 1, 1) # I = ileum, J = jejunum

# Remove cells belonging to cell types not well represented in both tissues:
df <- as.data.frame.matrix(table(seu$celltype2, seu$tissue)) # create data frame of cell type quantities per tissue
order <- c('Activated B cells', 'Cycling B cells', 'Resting B cells', 'Transitioning B cells', 'Antibody-secreting cells',
           'Cycling CD4 ab T cells', 'Cycling CD8 ab T cells', 'Cycling gd T cells', 
           'Cytotoxic CD8 ab T cells', 'Cytotoxic gd T cells', 'Cytotoxic group 1 ILCs', 
           'Non-naive CD8 ab T cells', 'Non-naive gd T cells', 'Non-naive group 1 ILCs', 
           'SELLhi gd T cells', 'CD2neg GD T cells', 
           'Naive CD4/CD8 ab T cells', 'Non-naive CD4 ab T cells', 'Follicular CD4 ab T cells', 'Group 3 ILCs',
           'Dendritic cells', 'CD4pos macrophage', 'CD4neg macrophage', 'Mast cells',
           'Crypt cells', 'Enterocytes', 'BEST4 enterocytes', 'Goblet cells', 'NEUROD1lo EE cells', 'NEUROD1hi EE cells',
           'Endothelial cells', 'Fibroblasts', 'Muscle cells') # define order of cell type levels ot use
df <- df[match(order, rownames(df)), ]  # change to order we would like to present annotations in
rm <- rownames(df[df$I < 10 | df$J < 10,]) # identify cell types that don't have at least 10 cells recovered in each tissue & remove these from further analysis
df <- df[!(row.names(df) %in% rm),]
group.new <- rownames(df) # this is a vector of cell types 
Idents(seu) <- seu$celltype2
seu <- subset(seu, idents = c(group.new)) # subset to only cell types with >= 10 cells in each tissue
Idents(seu) <- seu$celltype2
levels(seu) <- c('Activated B cells', 'Cycling B cells', 'Resting B cells', 'Transitioning B cells', 'Antibody-secreting cells',
                 'Cycling CD4 ab T cells', 'Cycling CD8 ab T cells', 
                 'Cytotoxic CD8 ab T cells', 'Cytotoxic gd T cells', 'Cytotoxic group 1 ILCs', 
                 'Non-naive CD8 ab T cells', 'Non-naive gd T cells', 'Non-naive group 1 ILCs', 
                 'SELLhi gd T cells', 'CD2neg GD T cells', 
                 'Naive CD4/CD8 ab T cells', 'Non-naive CD4 ab T cells', 'Follicular CD4 ab T cells', 'Group 3 ILCs',
                 'Dendritic cells', 'CD4pos macrophage', 'CD4neg macrophage', 'Mast cells',
                 'Crypt cells', 'Enterocytes', 
                 'Endothelial cells', 'Fibroblasts')
seu$celltype2 <- Idents(seu)
DefaultAssay(seu) <- 'SCT'

# Convert to human ortholog genes:
# Load in ortho genes
orthoGenes <- read.delim("/home/Jayne.Wiarda/SI_PP_SC_ST/GeneAnnotationFiles/PigToHuman_GeneOrthos_v97.txt") # read in gene ortholog file
orthoGenes <- subset(orthoGenes, Human.homology.type == 'ortholog_one2one') # subset to only one to one orthologs

# Filter  data to include only one-to-one gene orthologs & convert to human gene symbols:
genes <- as.data.frame(rownames(seu[['SCT']]@data)) # extract pig gene names from dataset
colnames(genes) <- 'gene'
pigGenes <- read_excel('/home/Jayne.Wiarda/SI_PP_SC_ST/GeneAnnotationFiles/UpdatedGeneNameListForSus97GTF_06302021_JEW_SKS.xlsx') # read in file with an updated gene symbol annotation for Sus scrofa v97 annotation build
pigGenes$FinalList <-gsub("_", "-", pigGenes$FinalList) # replace all underscores with dashes since this occurred when processing data in a previous step
pigGenes <- pigGenes[pigGenes$FinalList %in% genes$gene, ] # slim down to only genes in our dataset
orthos <- intersect(pigGenes$ENSID, orthoGenes$Gene.stable.ID) # find which genes are one-to-one orthologs
#length(orthos) # how many genes are orthologs?
pigGenes <- pigGenes[pigGenes$ENSID %in% orthos, ]
#dim(pigGenes)
orthoGenes <- orthoGenes[orthoGenes$Gene.stable.ID %in% pigGenes$ENSID, ] # slim down to only ortho genes in our dataset
orthoGenes <- orthoGenes %>% distinct(orthoGenes$Gene.stable.ID, orthoGenes$Human.gene.stable.ID, .keep_all = TRUE)  # retain only unique combinations of pig & human Ensembl IDs, ignoring transcript IDs
#dim(orthoGenes) # should not have same number of rows as in pigGenes
counts <- seu[['SCT']]@data[rownames(seu[['SCT']]@data) %in% pigGenes$FinalList,]
pigGenes <- pigGenes %>% arrange(factor(FinalList, levels = rownames(counts))) # arrange pigGenes in same order as counts
orthoGenes <- orthoGenes %>% arrange(factor(Gene.stable.ID, levels = pigGenes$ENSID)) # arrange orthoGenes in same order as pigGenes (and consequently counts)
rownames(counts) <- orthoGenes$Human.gene.name # change pig genes to human gene names

# Do a custom bit here to incorporate more MHC-II SLA genes, since some are one-to-many orthologs but we are highly interested in MHC-II related pathways for this work...
add <- setdiff(rownames(seu[['SCT']]@data)[which(grepl("SLA-D", rownames(seu[['SCT']]@data)))],
               pigGenes$FinalList[which(grepl("SLA-D", pigGenes$FinalList))]) # ID MHC-II genes from Seurat object that weren't recognized as one-to-one orthologs
counts2 <- seu[['SCT']]@data[rownames(seu[['SCT']]@data) %in% add,]
rownames(counts2) <- gsub('SLA-','HLA-', rownames(counts2)) # change to human gene symbol IDs
counts <- rbind(counts, counts2)

# Subset to ileum and jejunum separately:
Meta <- seu@meta.data
meta <- t(Meta)
meta <- as.data.frame(meta)
jcells <- colnames(meta %>% dplyr::select(starts_with("J")))
icells <- colnames(meta %>% dplyr::select(starts_with("I")))
countsj <- as.matrix(counts[, jcells])
countsi <- as.matrix(counts[, icells])

Meta <- as.data.frame(t(Meta)) # transpose so we can use select() function on column names
metaj <- t(Meta %>% dplyr::select(starts_with("J"))) #need cells as rows in final format so use t()
metai <- t(Meta %>% dplyr::select(starts_with("I")))

# Cell-cell signaling ----
# Create CellChat object:
cellchatj <- createCellChat(object = countsj, # use new Seurat object with human gene names
                            meta = metaj,
                            group.by = "celltype2") # set cell identities to cell type annotations
cellchati <- createCellChat(object = countsi, # use new Seurat object with human gene names
                            meta = metai,
                            group.by = "celltype2") # set cell identities to cell type annotations

# reset order:
cellchatj@idents = factor(cellchatj@idents,
                          levels = c('Activated B cells', 'Cycling B cells', 'Resting B cells', 'Transitioning B cells', 'Antibody-secreting cells',
                                     'Cycling CD4 ab T cells', 'Cycling CD8 ab T cells', 
                                     'Cytotoxic CD8 ab T cells', 'Cytotoxic gd T cells', 'Cytotoxic group 1 ILCs', 
                                     'Non-naive CD8 ab T cells', 'Non-naive gd T cells', 'Non-naive group 1 ILCs', 
                                     'SELLhi gd T cells', 'CD2neg GD T cells', 
                                     'Naive CD4/CD8 ab T cells', 'Non-naive CD4 ab T cells', 'Follicular CD4 ab T cells', 'Group 3 ILCs',
                                     'Dendritic cells', 'CD4pos macrophage', 'CD4neg macrophage', 'Mast cells',
                                     'Crypt cells', 'Enterocytes', 
                                     'Endothelial cells', 'Fibroblasts'))
cellchati@idents = factor(cellchati@idents,
                          levels = c('Activated B cells', 'Cycling B cells', 'Resting B cells', 'Transitioning B cells', 'Antibody-secreting cells',
                                     'Cycling CD4 ab T cells', 'Cycling CD8 ab T cells', 
                                     'Cytotoxic CD8 ab T cells', 'Cytotoxic gd T cells', 'Cytotoxic group 1 ILCs', 
                                     'Non-naive CD8 ab T cells', 'Non-naive gd T cells', 'Non-naive group 1 ILCs', 
                                     'SELLhi gd T cells', 'CD2neg GD T cells', 
                                     'Naive CD4/CD8 ab T cells', 'Non-naive CD4 ab T cells', 'Follicular CD4 ab T cells', 'Group 3 ILCs',
                                     'Dendritic cells', 'CD4pos macrophage', 'CD4neg macrophage', 'Mast cells',
                                     'Crypt cells', 'Enterocytes', 
                                     'Endothelial cells', 'Fibroblasts'))

# Set cell interaction database:
CellChatDB <- CellChatDB.human # use CellChatDB.mouse if running on mouse data
#showDatabaseCategory(CellChatDB)

# use a subset of CellChatDB for cell-cell communication analysis
# CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") # use Secreted Signaling
# use all CellChatDB for cell-cell communication analysis
CellChatDB.use <- subsetDB(CellChatDB, search = "Cell-Cell Contact") # look at only cell-cell contact
#CellChatDB.use <- CellChatDB # simply use the default CellChatDB

# set the used database in the object
cellchatj@DB <- CellChatDB.use
cellchati@DB <- CellChatDB.use

# Preprocess expression data
# subset the expression data of signaling genes for saving computation cost
cellchatj <- subsetData(cellchatj) # This step is necessary even if using the whole database
cellchati <- subsetData(cellchati) # This step is necessary even if using the whole database
future::plan("multisession", workers = 4) # do parallel

cellchatj <- identifyOverExpressedGenes(cellchatj)
cellchatj <- identifyOverExpressedInteractions(cellchatj)
cellchati <- identifyOverExpressedGenes(cellchati)
cellchati <- identifyOverExpressedInteractions(cellchati)

# Compute communication probabilities and network inferences:
#cellchat <- computeCommunProb(cellchat)
cellchatj <- computeCommunProb(cellchatj, type =  "truncatedMean", trim = 0.1) # count as zero expression if in <10% of annotated cell type (default = 25%)
# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
#cellchatj <- filterCommunication(cellchatj, min.cells = 10)
cellchati <- computeCommunProb(cellchati, type =  "truncatedMean", trim = 0.1) # count as zero expression if in <10% of annotated cell type (default = 25%)
# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
#cellchati <- filterCommunication(cellchati, min.cells = 10)

# Save CellChat object: # save intermediate here...overwrite with final version later
saveRDS(cellchatj, file = '/home/Jayne.Wiarda/SI_PP_SC_ST/CellChat/CellChat_Jejunum_CellCellSignalOnly_CD4macrophages.rds')
saveRDS(cellchati, file = '/home/Jayne.Wiarda/SI_PP_SC_ST/CellChat/CellChat_Ileum_CellCellSignalOnly_CD4macrophages.rds')

# Secreted signaling ----

# Create CellChat object:
cellchatj <- createCellChat(object = countsj, # use new Seurat object with human gene names
                            meta = metaj,
                            group.by = "celltype2") # set cell identities to cell type annotations
cellchati <- createCellChat(object = countsi, # use new Seurat object with human gene names
                            meta = metai,
                            group.by = "celltype2") # set cell identities to cell type annotations

# reset order:
cellchatj@idents = factor(cellchatj@idents,
                          levels = c('Activated B cells', 'Cycling B cells', 'Resting B cells', 'Transitioning B cells', 'Antibody-secreting cells',
                                     'Cycling CD4 ab T cells', 'Cycling CD8 ab T cells', 
                                     'Cytotoxic CD8 ab T cells', 'Cytotoxic gd T cells', 'Cytotoxic group 1 ILCs', 
                                     'Non-naive CD8 ab T cells', 'Non-naive gd T cells', 'Non-naive group 1 ILCs', 
                                     'SELLhi gd T cells', 'CD2neg GD T cells', 
                                     'Naive CD4/CD8 ab T cells', 'Non-naive CD4 ab T cells', 'Follicular CD4 ab T cells', 'Group 3 ILCs',
                                     'Dendritic cells', 'CD4pos macrophage', 'CD4neg macrophage', 'Mast cells',
                                     'Crypt cells', 'Enterocytes', 
                                     'Endothelial cells', 'Fibroblasts'))
cellchati@idents = factor(cellchati@idents,
                          levels = c('Activated B cells', 'Cycling B cells', 'Resting B cells', 'Transitioning B cells', 'Antibody-secreting cells',
                                     'Cycling CD4 ab T cells', 'Cycling CD8 ab T cells', 
                                     'Cytotoxic CD8 ab T cells', 'Cytotoxic gd T cells', 'Cytotoxic group 1 ILCs', 
                                     'Non-naive CD8 ab T cells', 'Non-naive gd T cells', 'Non-naive group 1 ILCs', 
                                     'SELLhi gd T cells', 'CD2neg GD T cells', 
                                     'Naive CD4/CD8 ab T cells', 'Non-naive CD4 ab T cells', 'Follicular CD4 ab T cells', 'Group 3 ILCs',
                                     'Dendritic cells', 'CD4pos macrophage', 'CD4neg macrophage', 'Mast cells',
                                     'Crypt cells', 'Enterocytes', 
                                     'Endothelial cells', 'Fibroblasts'))

# Set cell interaction database:
CellChatDB <- CellChatDB.human # use CellChatDB.mouse if running on mouse data
#showDatabaseCategory(CellChatDB)

# use a subset of CellChatDB for cell-cell communication analysis
# CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") # use Secreted Signaling
# use all CellChatDB for cell-cell communication analysis
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") # look at only secreted signaling
#CellChatDB.use <- CellChatDB # simply use the default CellChatDB

# set the used database in the object
cellchatj@DB <- CellChatDB.use
cellchati@DB <- CellChatDB.use

# Preprocess expression data
# subset the expression data of signaling genes for saving computation cost
cellchatj <- subsetData(cellchatj) # This step is necessary even if using the whole database
cellchati <- subsetData(cellchati) # This step is necessary even if using the whole database
future::plan("multisession", workers = 4) # do parallel

cellchatj <- identifyOverExpressedGenes(cellchatj)
cellchatj <- identifyOverExpressedInteractions(cellchatj)
cellchati <- identifyOverExpressedGenes(cellchati)
cellchati <- identifyOverExpressedInteractions(cellchati)

# Compute communication probabilities and network inferences:
#cellchat <- computeCommunProb(cellchat)
cellchatj <- computeCommunProb(cellchatj, type =  "truncatedMean", trim = 0.1) # count as zero expression if in <10% of annotated cell type (default = 25%)
# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
#cellchatj <- filterCommunication(cellchatj, min.cells = 10)
cellchati <- computeCommunProb(cellchati, type =  "truncatedMean", trim = 0.1) # count as zero expression if in <10% of annotated cell type (default = 25%)
# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
#cellchati <- filterCommunication(cellchati, min.cells = 10)

# Save CellChat object: # save intermediate here...overwrite with final version later
saveRDS(cellchatj, file = '/home/Jayne.Wiarda/SI_PP_SC_ST/CellChat/CellChat_Jejunum_SecretedSignalOnly_CD4macrophages.rds')
saveRDS(cellchati, file = '/home/Jayne.Wiarda/SI_PP_SC_ST/CellChat/CellChat_Ileum_SecretedSignalOnly_CD4macrophages.rds')

# Further process & visualize individual CellChat objects ----

# Load all in:
cellj <- readRDS('/home/Jayne.Wiarda/SI_PP_SC_ST/CellChat/CellChat_Jejunum_CellCellSignalOnly_CD4macrophages.rds')
celli <- readRDS('/home/Jayne.Wiarda/SI_PP_SC_ST/CellChat/CellChat_Ileum_CellCellSignalOnly_CD4macrophages.rds')
secj <- readRDS("/home/Jayne.Wiarda/SI_PP_SC_ST/CellChat/CellChat_Jejunum_SecretedSignalOnly_CD4macrophages.rds")
seci <- readRDS("/home/Jayne.Wiarda/SI_PP_SC_ST/CellChat/CellChat_Ileum_SecretedSignalOnly_CD4macrophages.rds")

# Infer cell signaling pathway communication & calculate aggregated communication network:
cellj <- computeCommunProbPathway(cellj)
cellj <- aggregateNet(cellj)
celli <- computeCommunProbPathway(celli)
celli <- aggregateNet(celli)
secj <- computeCommunProbPathway(secj)
secj <- aggregateNet(secj)
seci <- computeCommunProbPathway(seci)
seci <- aggregateNet(seci)

# Plot interactions between cell types:
cols <- c('cyan4', 'gold3', 'chartreuse4', 'deeppink4', 'sandybrown',
          'cornflowerblue', 'navy', 'salmon', 
          'deepskyblue2', 'tan4', 'mediumpurple1', 
          'darkgreen', 'gray50', 'darkmagenta', 'red', 'hotpink', 'khaki', 
          'orange4', 'limegreen', 'cadetblue3', 'firebrick3', 'indianred4', 'deepskyblue4', 'darkseagreen', 'burlywood3', 
          'lightskyblue3', 'mistyrose3')

groupSize <- as.numeric(table(cellj@idents))
par(mfrow=c(1,1))
g1 <- netVisual_heatmap(cellj, measure = 'count', color.heatmap = "BuGn", color.use = cols)
g2 <- netVisual_heatmap(cellj, measure = 'weight', color.heatmap = "BuGn", color.use = cols)
g1 + g2

groupSize <- as.numeric(table(celli@idents))
par(mfrow=c(1,1))
g1 <- netVisual_heatmap(celli, measure = 'count', color.heatmap = "BuGn", color.use = cols)
g2 <- netVisual_heatmap(celli, measure = 'weight', color.heatmap = "BuGn", color.use = cols)
g1 + g2

groupSize <- as.numeric(table(secj@idents))
par(mfrow=c(1,1))
g1 <- netVisual_heatmap(secj, measure = 'count', color.heatmap = "BuGn", color.use = cols)
g2 <- netVisual_heatmap(secj, measure = 'weight', color.heatmap = "BuGn", color.use = cols)
g1 + g2

groupSize <- as.numeric(table(seci@idents))
par(mfrow=c(1,1))
g1 <- netVisual_heatmap(seci, measure = 'count', color.heatmap = "BuGn", color.use = cols)
g2 <- netVisual_heatmap(seci, measure = 'weight', color.heatmap = "BuGn", color.use = cols)
g1 + g2

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
#  [1] stats4    stats     graphics  grDevices utils     datasets  methods   base     

#other attached packages:
#  [1] ggnewscale_0.4.8      ggrepel_0.9.1         topGO_2.48.0          SparseM_1.81          GO.db_3.15.0          AnnotationDbi_1.58.0  IRanges_2.30.1        S4Vectors_0.34.0      Biobase_2.56.0        graph_1.74.0          BiocGenerics_0.42.0  
#[12] biomaRt_2.52.0        forcats_0.5.2         stringr_1.4.1         purrr_0.3.4           readr_2.1.2           tidyr_1.2.1           tibble_3.1.8          tidyverse_1.3.2       CellChat_1.6.1        igraph_1.3.4          writexl_1.4.0        
#[23] readxl_1.4.1          dplyr_1.0.10          ggplot2_3.3.6         SeuratDisk_0.0.0.9020 sp_1.5-0              SeuratObject_4.1.1    Seurat_4.1.1         

#loaded via a namespace (and not attached):
#  [1] rappdirs_0.3.3         scattermore_0.8        coda_0.19-4            pkgmaker_0.32.2        knitr_1.40             bit64_4.0.5            irlba_2.3.5            data.table_1.14.2      rpart_4.1.16           KEGGREST_1.36.3        RCurl_1.98-1.8        
#[12] doParallel_1.0.17      generics_0.1.3         cowplot_1.1.1          RSQLite_2.2.17         RANN_2.6.1             future_1.28.0          bit_4.0.4              tzdb_0.3.0             spatstat.data_2.2-0    xml2_1.3.3             lubridate_1.9.0       
#[23] httpuv_1.6.6           assertthat_0.2.1       gargle_1.2.1           xfun_0.33              hms_1.1.2              evaluate_0.16          promises_1.2.0.1       fansi_1.0.3            progress_1.2.2         dbplyr_2.2.1           DBI_1.1.3             
#[34] htmlwidgets_1.5.4      spatstat.geom_2.4-0    googledrive_2.0.0      ellipsis_0.3.2         RSpectra_0.16-1        ggpubr_0.4.0           backports_1.4.1        gridBase_0.4-7         deldir_1.0-6           vctrs_0.4.1            ggalluvial_0.12.3     
#[45] Cairo_1.6-0            ROCR_1.0-11            abind_1.4-5            cachem_1.0.6           withr_2.5.0            progressr_0.11.0       vroom_1.5.7            sctransform_0.3.4      sna_2.7                prettyunits_1.1.1      goftest_1.2-3         
#[56] svglite_2.1.0          cluster_2.1.4          lazyeval_0.2.2         crayon_1.5.1           hdf5r_1.3.5            labeling_0.4.2         pkgconfig_2.0.3        GenomeInfoDb_1.32.4    nlme_3.1-159           rlang_1.0.6            globals_0.16.1        
#[67] lifecycle_1.0.2        miniUI_0.1.1.1         registry_0.5-1         filelock_1.0.2         BiocFileCache_2.4.0    modelr_0.1.10          cellranger_1.1.0       polyclip_1.10-0        matrixStats_0.62.0     lmtest_0.9-40          rngtools_1.5.2        
#[78] Matrix_1.5-1           carData_3.0-5          zoo_1.8-10             reprex_2.0.2           ggridges_0.5.3         GlobalOptions_0.1.2    googlesheets4_1.0.1    png_0.1-7              viridisLite_0.4.1      rjson_0.2.21           bitops_1.0-7          
#[89] KernSmooth_2.23-20     ggnetwork_0.5.10       Biostrings_2.64.1      blob_1.2.3             shape_1.4.6            parallelly_1.32.1      spatstat.random_2.2-0  rstatix_0.7.0          ggsignif_0.6.3         scales_1.2.1           memoise_2.0.1         
#[100] magrittr_2.0.3         plyr_1.8.7             ica_1.0-3              zlibbioc_1.42.0        compiler_4.2.2         RColorBrewer_1.1-3     clue_0.3-61            fitdistrplus_1.1-8     cli_3.4.0              XVector_0.36.0         listenv_0.8.0         
#[111] patchwork_1.1.2        pbapply_1.5-0          MASS_7.3-58.1          mgcv_1.8-40            tidyselect_1.1.2       stringi_1.7.8          yaml_2.3.5             grid_4.2.2             tools_4.2.2            timechange_0.2.0       future.apply_1.9.1    
#[122] parallel_4.2.2         circlize_0.4.15        rstudioapi_0.14        foreach_1.5.2          gridExtra_2.3          farver_2.1.1           Rtsne_0.16             digest_0.6.29          rgeos_0.5-9            FNN_1.1.3.1            shiny_1.7.2           
#[133] Rcpp_1.0.9             car_3.1-0              broom_1.0.1            later_1.3.0            RcppAnnoy_0.0.19       httr_1.4.4             ComplexHeatmap_2.12.1  colorspace_2.0-3       rvest_1.0.3            XML_3.99-0.10          fs_1.5.2              
#[144] tensor_1.5             reticulate_1.26        splines_4.2.2          uwot_0.1.14            spatstat.utils_2.3-1   plotly_4.10.0          systemfonts_1.0.4      xtable_1.8-4           jsonlite_1.8.0         R6_2.5.1               pillar_1.8.1          
#[155] htmltools_0.5.3        mime_0.12              NMF_0.24.0             glue_1.6.2             fastmap_1.1.0          BiocParallel_1.30.3    BiocNeighbors_1.14.0   codetools_0.2-18       utf8_1.2.2             lattice_0.20-45        spatstat.sparse_2.1-1 
#[166] network_1.17.2         curl_4.3.2             leiden_0.4.3           limma_3.52.3           survival_3.4-0         rmarkdown_2.16         statnet.common_4.7.0   munsell_0.5.0          GetoptLong_1.0.5       GenomeInfoDbData_1.2.8 iterators_1.0.14      
#[177] haven_2.5.1            reshape2_1.4.4         gtable_0.3.1           spatstat.core_2.4-4   