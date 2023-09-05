library(Seurat)
library(SeuratDisk)
library(CellChat)
library(patchwork)
library(readxl)
library(future)
library(NMF) 
library(uwot)
library(dplyr)
library(ggplot2)
library(tidyr)
library(scales)
library(ComplexHeatmap)
library(viridis)

# Load Seurat object:
seu <- LoadH5Seurat('/home/Jayne.Wiarda/SI_PP_SC_ST/AnnotatedSeurat/AllSamples_annotated.h5seurat')
DefaultAssay(seu) <- 'SCT'
seu$tissue <- substr(seu$orig.ident, 1, 1) # I = ileum, J = jejunum

# Remove cells belonging to cell types not well represented in both tissues:
df <- as.data.frame.matrix(table(seu$celltype, seu$tissue)) # create data frame of cell type quantities per tissue
order <- c('Activated B cells', 'Cycling B cells', 'Resting B cells', 'Transitioning B cells', 'Antibody-secreting cells',
           'Cycling CD4 ab T cells', 'Cycling CD8 ab T cells', 'Cycling gd T cells', 
           'Cytotoxic CD8 ab T cells', 'Cytotoxic gd T cells', 'Cytotoxic group 1 ILCs', 
           'Non-naive CD8 ab T cells', 'Non-naive gd T cells', 'Non-naive group 1 ILCs', 
           'SELLhi gd T cells', 'CD2neg GD T cells', 
           'Naive CD4/CD8 ab T cells', 'Non-naive CD4 ab T cells', 'Follicular CD4 ab T cells', 'Group 3 ILCs',
           'Dendritic cells', 'Macrophages', 'Mast cells',
           'Crypt cells', 'Enterocytes', 'BEST4 enterocytes', 'Goblet cells', 'NEUROD1lo EE cells', 'NEUROD1hi EE cells',
           'Endothelial cells', 'Fibroblasts', 'Muscle cells') # define order of cell type levels ot use
df <- df[match(order, rownames(df)), ]  # change to order we would like to present annotations in
rm <- rownames(df[df$I < 10 | df$J < 10,]) # identify cell types that don't have at least 10 cells recovered in each tissue & remove these from further analysis
df <- df[!(row.names(df) %in% rm),]
group.new <- rownames(df) # this is a vector of cell types 
Idents(seu) <- seu$celltype
seu <- subset(seu, idents = c(group.new)) # subset to only cell types with >= 10 cells in each tissue
Idents(seu) <- seu$celltype
levels(seu) <- c('Activated B cells', 'Cycling B cells', 'Resting B cells', 'Transitioning B cells', 'Antibody-secreting cells',
                 'Cycling CD4 ab T cells', 'Cycling CD8 ab T cells', 
                 'Cytotoxic CD8 ab T cells', 'Cytotoxic gd T cells', 'Cytotoxic group 1 ILCs', 
                 'Non-naive CD8 ab T cells', 'Non-naive gd T cells', 'Non-naive group 1 ILCs', 
                 'SELLhi gd T cells', 'CD2neg GD T cells', 
                 'Naive CD4/CD8 ab T cells', 'Non-naive CD4 ab T cells', 'Follicular CD4 ab T cells', 'Group 3 ILCs',
                 'Dendritic cells', 'Macrophages', 'Mast cells',
                 'Crypt cells', 'Enterocytes', 
                 'Endothelial cells', 'Fibroblasts')
seu$celltype <- Idents(seu)
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
                            group.by = "celltype") # set cell identities to cell type annotations
cellchati <- createCellChat(object = countsi, # use new Seurat object with human gene names
                            meta = metai,
                            group.by = "celltype") # set cell identities to cell type annotations

# reset order:
cellchatj@idents = factor(cellchatj@idents,
                          levels = c('Activated B cells', 'Cycling B cells', 'Resting B cells', 'Transitioning B cells', 'Antibody-secreting cells',
                                     'Cycling CD4 ab T cells', 'Cycling CD8 ab T cells', 
                                     'Cytotoxic CD8 ab T cells', 'Cytotoxic gd T cells', 'Cytotoxic group 1 ILCs', 
                                     'Non-naive CD8 ab T cells', 'Non-naive gd T cells', 'Non-naive group 1 ILCs', 
                                     'SELLhi gd T cells', 'CD2neg GD T cells', 
                                     'Naive CD4/CD8 ab T cells', 'Non-naive CD4 ab T cells', 'Follicular CD4 ab T cells', 'Group 3 ILCs',
                                     'Dendritic cells', 'Macrophages', 'Mast cells',
                                     'Crypt cells', 'Enterocytes', 
                                     'Endothelial cells', 'Fibroblasts'))
cellchati@idents = factor(cellchati@idents,
                          levels = c('Activated B cells', 'Cycling B cells', 'Resting B cells', 'Transitioning B cells', 'Antibody-secreting cells',
                                     'Cycling CD4 ab T cells', 'Cycling CD8 ab T cells', 
                                     'Cytotoxic CD8 ab T cells', 'Cytotoxic gd T cells', 'Cytotoxic group 1 ILCs', 
                                     'Non-naive CD8 ab T cells', 'Non-naive gd T cells', 'Non-naive group 1 ILCs', 
                                     'SELLhi gd T cells', 'CD2neg GD T cells', 
                                     'Naive CD4/CD8 ab T cells', 'Non-naive CD4 ab T cells', 'Follicular CD4 ab T cells', 'Group 3 ILCs',
                                     'Dendritic cells', 'Macrophages', 'Mast cells',
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
saveRDS(cellchatj, file = '/home/Jayne.Wiarda/SI_PP_SC_ST/CellChat/CellChat_Jejunum_CellCellSignalOnly.rds')
saveRDS(cellchati, file = '/home/Jayne.Wiarda/SI_PP_SC_ST/CellChat/CellChat_Ileum_CellCellSignalOnly.rds')

# Secreted signaling ----

# Create CellChat object:
cellchatj <- createCellChat(object = countsj, # use new Seurat object with human gene names
                            meta = metaj,
                            group.by = "celltype") # set cell identities to cell type annotations
cellchati <- createCellChat(object = countsi, # use new Seurat object with human gene names
                            meta = metai,
                            group.by = "celltype") # set cell identities to cell type annotations

# reset order:
cellchatj@idents = factor(cellchatj@idents,
                          levels = c('Activated B cells', 'Cycling B cells', 'Resting B cells', 'Transitioning B cells', 'Antibody-secreting cells',
                                     'Cycling CD4 ab T cells', 'Cycling CD8 ab T cells', 
                                     'Cytotoxic CD8 ab T cells', 'Cytotoxic gd T cells', 'Cytotoxic group 1 ILCs', 
                                     'Non-naive CD8 ab T cells', 'Non-naive gd T cells', 'Non-naive group 1 ILCs', 
                                     'SELLhi gd T cells', 'CD2neg GD T cells', 
                                     'Naive CD4/CD8 ab T cells', 'Non-naive CD4 ab T cells', 'Follicular CD4 ab T cells', 'Group 3 ILCs',
                                     'Dendritic cells', 'Macrophages', 'Mast cells',
                                     'Crypt cells', 'Enterocytes', 
                                     'Endothelial cells', 'Fibroblasts'))
cellchati@idents = factor(cellchati@idents,
                          levels = c('Activated B cells', 'Cycling B cells', 'Resting B cells', 'Transitioning B cells', 'Antibody-secreting cells',
                                     'Cycling CD4 ab T cells', 'Cycling CD8 ab T cells', 
                                     'Cytotoxic CD8 ab T cells', 'Cytotoxic gd T cells', 'Cytotoxic group 1 ILCs', 
                                     'Non-naive CD8 ab T cells', 'Non-naive gd T cells', 'Non-naive group 1 ILCs', 
                                     'SELLhi gd T cells', 'CD2neg GD T cells', 
                                     'Naive CD4/CD8 ab T cells', 'Non-naive CD4 ab T cells', 'Follicular CD4 ab T cells', 'Group 3 ILCs',
                                     'Dendritic cells', 'Macrophages', 'Mast cells',
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
saveRDS(cellchatj, file = '/home/Jayne.Wiarda/SI_PP_SC_ST/CellChat/CellChat_Jejunum_SecretedSignalOnly.rds')
saveRDS(cellchati, file = '/home/Jayne.Wiarda/SI_PP_SC_ST/CellChat/CellChat_Ileum_SecretedSignalOnly.rds')

# ECM signaling ----

# Create CellChat object:
cellchatj <- createCellChat(object = countsj, # use new Seurat object with human gene names
                            meta = metaj,
                            group.by = "celltype") # set cell identities to cell type annotations
cellchati <- createCellChat(object = countsi, # use new Seurat object with human gene names
                            meta = metai,
                            group.by = "celltype") # set cell identities to cell type annotations

# reset order:
cellchatj@idents = factor(cellchatj@idents,
                          levels = c('Activated B cells', 'Cycling B cells', 'Resting B cells', 'Transitioning B cells', 'Antibody-secreting cells',
                                     'Cycling CD4 ab T cells', 'Cycling CD8 ab T cells', 
                                     'Cytotoxic CD8 ab T cells', 'Cytotoxic gd T cells', 'Cytotoxic group 1 ILCs', 
                                     'Non-naive CD8 ab T cells', 'Non-naive gd T cells', 'Non-naive group 1 ILCs', 
                                     'SELLhi gd T cells', 'CD2neg GD T cells', 
                                     'Naive CD4/CD8 ab T cells', 'Non-naive CD4 ab T cells', 'Follicular CD4 ab T cells', 'Group 3 ILCs',
                                     'Dendritic cells', 'Macrophages', 'Mast cells',
                                     'Crypt cells', 'Enterocytes', 
                                     'Endothelial cells', 'Fibroblasts'))
cellchati@idents = factor(cellchati@idents,
                          levels = c('Activated B cells', 'Cycling B cells', 'Resting B cells', 'Transitioning B cells', 'Antibody-secreting cells',
                                     'Cycling CD4 ab T cells', 'Cycling CD8 ab T cells', 
                                     'Cytotoxic CD8 ab T cells', 'Cytotoxic gd T cells', 'Cytotoxic group 1 ILCs', 
                                     'Non-naive CD8 ab T cells', 'Non-naive gd T cells', 'Non-naive group 1 ILCs', 
                                     'SELLhi gd T cells', 'CD2neg GD T cells', 
                                     'Naive CD4/CD8 ab T cells', 'Non-naive CD4 ab T cells', 'Follicular CD4 ab T cells', 'Group 3 ILCs',
                                     'Dendritic cells', 'Macrophages', 'Mast cells',
                                     'Crypt cells', 'Enterocytes', 
                                     'Endothelial cells', 'Fibroblasts'))

# Set cell interaction database:
CellChatDB <- CellChatDB.human # use CellChatDB.mouse if running on mouse data
#showDatabaseCategory(CellChatDB)

# use a subset of CellChatDB for cell-cell communication analysis
# CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") # use Secreted Signaling
# use all CellChatDB for cell-cell communication analysis
CellChatDB.use <- subsetDB(CellChatDB, search = "ECM-Receptor") # look at only ECM signaling
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
saveRDS(cellchatj, file = '/home/Jayne.Wiarda/SI_PP_SC_ST/CellChat/CellChat_Jejunum_ECMSignalOnly.rds')
saveRDS(cellchati, file = '/home/Jayne.Wiarda/SI_PP_SC_ST/CellChat/CellChat_Ileum_ECMSignalOnly.rds')

# Further process & visualize individual CellChat objects ----

# Load all in:
cellj <- readRDS("/home/Jayne.Wiarda/SI_PP_SC_ST/CellChat/CellChat_Jejunum_CellCellSignalOnly.rds")
celli <- readRDS("/home/Jayne.Wiarda/SI_PP_SC_ST/CellChat/CellChat_Ileum_CellCellSignalOnly.rds")
secj <- readRDS("/home/Jayne.Wiarda/SI_PP_SC_ST/CellChat/CellChat_Jejunum_SecretedSignalOnly.rds")
seci <- readRDS("/home/Jayne.Wiarda/SI_PP_SC_ST/CellChat/CellChat_Ileum_SecretedSignalOnly.rds")
ecmj <- readRDS("/home/Jayne.Wiarda/SI_PP_SC_ST/CellChat/CellChat_Jejunum_ECMSignalOnly.rds")
ecmi <- readRDS("/home/Jayne.Wiarda/SI_PP_SC_ST/CellChat/CellChat_Ileum_ECMSignalOnly.rds")

# Infer cell signaling pathway communication & calculate aggregated communication network:
cellj <- computeCommunProbPathway(cellj)
cellj <- aggregateNet(cellj)
celli <- computeCommunProbPathway(celli)
celli <- aggregateNet(celli)
secj <- computeCommunProbPathway(secj)
secj <- aggregateNet(secj)
seci <- computeCommunProbPathway(seci)
seci <- aggregateNet(seci)
ecmj <- computeCommunProbPathway(ecmj)
ecmj <- aggregateNet(ecmj)
ecmi <- computeCommunProbPathway(ecmi)
ecmi <- aggregateNet(ecmi)

# Plot interactions between cell types:
cols <- c('cyan4', 'gold3', 'chartreuse4', 'deeppink4', 'sandybrown',
         'cornflowerblue', 'navy', 'salmon', 
         'deepskyblue2', 'tan4', 'mediumpurple1', 
         'darkgreen', 'gray50', 'darkmagenta', 'red', 'hotpink', 'khaki', 
         'orange4', 'limegreen', 'cadetblue3', 'firebrick', 'deepskyblue4', 'darkseagreen', 'burlywood3', 
         'lightskyblue3', 'mistyrose3')

groupSize <- as.numeric(table(cellj@idents))
par(mfrow=c(1,1))
g1 <- netVisual_heatmap(cellj, measure = 'count', color.heatmap = "BuGn", color.use = cols)
g2 <- netVisual_heatmap(cellj, measure = 'weight', color.heatmap = "BuGn", color.use = cols)
g1 + g2
#levels(cellj@idents)
#netVisual_bubble(cellj, targets.use = c(1:5), remove.isolate = FALSE)
#netVisual_bubble(cellj, sources.use = c(1:5), remove.isolate = FALSE)

groupSize <- as.numeric(table(celli@idents))
par(mfrow=c(1,1))
g1 <- netVisual_heatmap(celli, measure = 'count', color.heatmap = "BuGn", color.use = cols)
g2 <- netVisual_heatmap(celli, measure = 'weight', color.heatmap = "BuGn", color.use = cols)
g1 + g2
#levels(celli@idents)
#netVisual_bubble(celli, targets.use = c(1:5), remove.isolate = FALSE)
#netVisual_bubble(celli, sources.use = c(1:5), remove.isolate = FALSE)

groupSize <- as.numeric(table(secj@idents))
par(mfrow=c(1,1))
g1 <- netVisual_heatmap(secj, measure = 'count', color.heatmap = "BuGn", color.use = cols)
g2 <- netVisual_heatmap(secj, measure = 'weight', color.heatmap = "BuGn", color.use = cols)
g1 + g2
#levels(secj@idents)
#netVisual_bubble(secj, targets.use = c(1:5), remove.isolate = FALSE)
#netVisual_bubble(secj, sources.use = c(1:5), remove.isolate = FALSE)

groupSize <- as.numeric(table(seci@idents))
par(mfrow=c(1,1))
g1 <- netVisual_heatmap(seci, measure = 'count', color.heatmap = "BuGn", color.use = cols)
g2 - netVisual_heatmap(seci, measure = 'weight', color.heatmap = "BuGn", color.use = cols)
g1 + g2
#levels(seci@idents)
#netVisual_bubble(seci, targets.use = c(1:5), remove.isolate = FALSE)
#netVisual_bubble(seci, sources.use = c(1:5), remove.isolate = FALSE)

groupSize <- as.numeric(table(ecmj@idents))
par(mfrow=c(1,1))
g1 <- netVisual_heatmap(ecmj, measure = 'count', color.heatmap = "BuGn", color.use = cols)
g2 <- netVisual_heatmap(ecmj, measure = 'weight', color.heatmap = "BuGn", color.use = cols)
g1 + g2
#levels(ecmj@idents)
#netVisual_bubble(ecmj, targets.use = c(1:5), remove.isolate = FALSE)
#netVisual_bubble(ecmj, sources.use = c(1:5), remove.isolate = FALSE)

groupSize <- as.numeric(table(ecmi@idents))
par(mfrow=c(1,1))
g1 <- netVisual_heatmap(ecmi, measure = 'count', color.heatmap = "BuGn", color.use = cols)
g2 <- netVisual_heatmap(ecmi, measure = 'weight', color.heatmap = "BuGn", color.use = cols)
g1 + g2
#levels(ecmi@idents)
#netVisual_bubble(ecmi, targets.use = c(1:5), remove.isolate = FALSE)
#netVisual_bubble(ecmi, sources.use = c(1:5), remove.isolate = FALSE)

# other plots to consider:
#par(mfrow = c(1,2), xpd=TRUE)
#netVisual_circle(cellj@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions - Jejunum") # hard to read these plots if we have lots of cell types
#netVisual_circle(cellj@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength - Jejunum")
#netVisual_chord_gene(cellj, targets.use = c(1), lab.cex = 0.5,legend.pos.y = 30, color.use = cols) # can run this for different cell types
#netVisual_chord_gene(cellj, sources.use = c(1), lab.cex = 0.5,legend.pos.y = 30, color.use = cols) # can run this for different cell types
#netVisual_chord_cell(cellj, lab.cex = 0.5, legend.pos.y = 30, slot.name = 'netP', signaling = 'MHC-II', color.use = cols) # can run this for different pathways

# Compute the network centrality scores
cellj <- netAnalysis_computeCentrality(cellj, slot.name = "netP") # the slot 'netP' means the inferred intercellular communication network of signaling pathways
celli <- netAnalysis_computeCentrality(celli, slot.name = "netP")
secj <- netAnalysis_computeCentrality(secj, slot.name = "netP")
seci <- netAnalysis_computeCentrality(seci, slot.name = "netP")
ecmj <- netAnalysis_computeCentrality(ecmj, slot.name = "netP") 
ecmi <- netAnalysis_computeCentrality(ecmi, slot.name = "netP")

netAnalysis_signalingRole_network(cellj, width = 11, height = 2.5, font.size = 10, color.heatmap = "BuGn", color.use = cols, signaling = 'MHC-II') # can change to look at other signaling pathways
#gg1 <- netAnalysis_signalingRole_scatter(cellj, color.use = cols)
#gg2 <- netAnalysis_signalingRole_scatter(cellj, color.use = cols)
#gg1 + gg2
#ht1 <- netAnalysis_signalingRole_heatmap(cellj, pattern = "outgoing", color.heatmap = "BuGn", color.use = cols, height = 14, width = 8)
#ht2 <- netAnalysis_signalingRole_heatmap(cellj, pattern = "incoming", color.heatmap = "BuGn", color.use = cols, height = 14, width = 8)
#ht1 + ht2

netAnalysis_signalingRole_network(celli, width = 11, height = 2.5, font.size = 10, color.heatmap = "BuGn", color.use = cols, signaling = 'MHC-II')
#gg1 <- netAnalysis_signalingRole_scatter(celli, color.use = cols)
#gg2 <- netAnalysis_signalingRole_scatter(celli, color.use = cols)
#gg1 + gg2
#ht1 <- netAnalysis_signalingRole_heatmap(celli, pattern = "outgoing", color.heatmap = "BuGn", color.use = cols, height = 14, width = 8)
#ht2 <- netAnalysis_signalingRole_heatmap(celli, pattern = "incoming", color.heatmap = "BuGn", color.use = cols, height = 14, width = 8)
#ht1 + ht2

netAnalysis_signalingRole_network(secj, width = 11, height = 2.5, font.size = 10, color.heatmap = "BuGn", color.use = cols)
#gg1 <- netAnalysis_signalingRole_scatter(secj, color.use = cols)
#gg2 <- netAnalysis_signalingRole_scatter(secj, color.use = cols)
#gg1 + gg2
#ht1 <- netAnalysis_signalingRole_heatmap(secj, pattern = "outgoing", color.heatmap = "BuGn", color.use = cols, height = 14, width = 8)
#ht2 <- netAnalysis_signalingRole_heatmap(secj, pattern = "incoming", color.heatmap = "BuGn", color.use = cols, height = 14, width = 8)
#ht1 + ht2

netAnalysis_signalingRole_network(seci, width = 11, height = 2.5, font.size = 10, color.heatmap = "BuGn", color.use = cols)
#gg1 <- netAnalysis_signalingRole_scatter(seci, color.use = cols)
#gg2 <- netAnalysis_signalingRole_scatter(seci, color.use = cols)
#gg1 + gg2
#ht1 <- netAnalysis_signalingRole_heatmap(seci, pattern = "outgoing", color.heatmap = "BuGn", color.use = cols, height = 14, width = 8)
#ht2 <- netAnalysis_signalingRole_heatmap(seci, pattern = "incoming", color.heatmap = "BuGn", color.use = cols, height = 14, width = 8)
#ht1 + ht2

netAnalysis_signalingRole_network(ecmj, width = 11, height = 2.5, font.size = 10, color.heatmap = "BuGn", color.use = cols)
#gg1 <- netAnalysis_signalingRole_scatter(ecmj, color.use = cols)
#gg2 <- netAnalysis_signalingRole_scatter(ecmj, color.use = cols)
#gg1 + gg2
#ht1 <- netAnalysis_signalingRole_heatmap(ecmj, pattern = "outgoing", color.heatmap = "BuGn", color.use = cols, height = 5, width = 8)
#ht2 <- netAnalysis_signalingRole_heatmap(ecmj, pattern = "incoming", color.heatmap = "BuGn", color.use = cols, height = 5, width = 8)
#ht1 + ht2

netAnalysis_signalingRole_network(ecmi, width = 11, height = 2.5, font.size = 10, color.heatmap = "BuGn", color.use = cols)
#gg1 <- netAnalysis_signalingRole_scatter(ecmi, color.use = cols)
#gg2 <- netAnalysis_signalingRole_scatter(ecmi, color.use = cols)
#gg1 + gg2
#ht1 <- netAnalysis_signalingRole_heatmap(ecmi, pattern = "outgoing", color.heatmap = "BuGn", color.use = cols, height = 5, width = 8)
#ht2 <- netAnalysis_signalingRole_heatmap(ecmi, pattern = "incoming", color.heatmap = "BuGn", color.use = cols, height = 5, width = 8)
#ht1 + ht2

# SOme additional pathways analyses you can do:
# Identify signaling groups based on their functional similarity
#cellj <- computeNetSimilarity(cellj, type = "functional")
#cellj <- netEmbedding(cellj, type = "functional", umap.method = 'uwot')
#cellj <- netClustering(cellj, type = "functional")
#netVisual_embedding(cellj, type = "functional", label.size = 3.5)
#netVisual_embeddingZoomIn(cellj, type = "functional", nCol = 2)

# Identify signaling groups based on structure similarity
#cellj <- computeNetSimilarity(cellj, type = "structural")
#cellj <- netEmbedding(cellj, type = "structural", umap.method = 'uwot')
#cellj <- netClustering(cellj, type = "structural")
#netVisual_embedding(cellj, type = "structural", label.size = 3.5)
#netVisual_embeddingZoomIn(cellj, type = "structural", nCol = 2)

# Save CellChat object:
saveRDS(cellj, file = '/home/Jayne.Wiarda/SI_PP_SC_ST/CellChat/CellChat_Jejunum_CellCellSignalOnly.rds')
saveRDS(celli, file = '/home/Jayne.Wiarda/SI_PP_SC_ST/CellChat/CellChat_Ileum_CellCellSignalOnly.rds')
saveRDS(secj, file = '/home/Jayne.Wiarda/SI_PP_SC_ST/CellChat/CellChat_Jejunum_SecretedSignalOnly.rds')
saveRDS(seci, file = '/home/Jayne.Wiarda/SI_PP_SC_ST/CellChat/CellChat_Ileum_SecretedSignalOnly.rds')
saveRDS(ecmj, file = '/home/Jayne.Wiarda/SI_PP_SC_ST/CellChat/CellChat_Jejunum_ECMSignalOnly.rds')
saveRDS(ecmi, file = '/home/Jayne.Wiarda/SI_PP_SC_ST/CellChat/CellChat_Ileum_ECMSignalOnly.rds')

# Merge CellChat objects ----
cellj <- updateCellChat(cellj)
celli <- updateCellChat(celli)
object.list <- list(jejunum = cellj, ileum = celli)
cell <- mergeCellChat(object.list, add.names = names(object.list)) # merge into a single CellChat object
# don't use liftCellChat() function since we filtered out cell types not represented in both tissues prior to making individual CellChat objects

secj <- updateCellChat(secj)
seci <- updateCellChat(seci)
object.list <- list(jejunum = secj, ileum = seci)
sec <- mergeCellChat(object.list, add.names = names(object.list))

ecmj <- updateCellChat(ecmj)
ecmi <- updateCellChat(ecmi)
object.list <- list(jejunum = ecmj, ileum = ecmi)
ecm <- mergeCellChat(object.list, add.names = names(object.list))

# Save CellChat object:
saveRDS(cell, file = '/home/Jayne.Wiarda/SI_PP_SC_ST/CellChat/CellChat_MergedIleumJejunum_CellCellSignalOnly.rds')
saveRDS(sec, file = '/home/Jayne.Wiarda/SI_PP_SC_ST/CellChat/CellChat_MergedIleumJejunum_SecretedSignalOnly.rds')
saveRDS(ecm, file = '/home/Jayne.Wiarda/SI_PP_SC_ST/CellChat/CellChat_MergedIleumJejunum_ECMSignalOnly.rds')

# Further process & visualize merged CellChat objects ----
# Compare number and strength of interactions between tissues:
gg1 <- compareInteractions(cell, show.legend = F, group = c(1,2), color.use = c('blue', 'red')) # compare number of interactions
gg2 <- compareInteractions(cell, show.legend = F, group = c(1,2), measure = "weight", color.use = c('blue', 'red')) # compare strength of interactions
gg1 + gg2

gg1 <- compareInteractions(sec, show.legend = F, group = c(1,2), color.use = c('blue', 'red')) # compare number of interactions
gg2 <- compareInteractions(sec, show.legend = F, group = c(1,2), measure = "weight", color.use = c('blue', 'red')) # compare strength of interactions
gg1 + gg2

gg1 <- compareInteractions(ecm, show.legend = F, group = c(1,2), color.use = c('blue', 'red')) # compare number of interactions
gg2 <- compareInteractions(ecm, show.legend = F, group = c(1,2), measure = "weight", color.use = c('blue', 'red')) # compare strength of interactions
gg1 + gg2

# Identify differential interactions: ## Too hard to read with complex data/many cell types
#netVisual_diffInteraction(cell, weight.scale = T) # number of interactions
#netVisual_diffInteraction(cell, weight.scale = T, measure = "weight") # strength of interactions
# red  (or blue) colored edges represent increased (or decreased) signaling in the ileum compared to the jejunum

# 2D plot of differential signaling by cell type
## Signaling number:
netAnalysis_diff_signalingRole_scatter(cell,
                                       x.measure = "outdeg_unweighted",
                                       y.measure = "indeg_unweighted",
                                       xlabel = "Outgoing interaction number",
                                       ylabel = "Incoming interaction number",
                                       color.use = cols)
## Signaling strength:
netAnalysis_diff_signalingRole_scatter(cell,
                                       x.measure = "outdeg",
                                       y.measure = "indeg",
                                       xlabel = "Outgoing interaction strength",
                                       ylabel = "Incoming interaction strength",
                                       color.use = cols)
# Positive values indicate the increase in the second dataset (ileum) while negative values indicate the increase in the first dataset (jejunum)

## Signaling number:
netAnalysis_diff_signalingRole_scatter(sec,
                                       x.measure = "outdeg_unweighted",
                                       y.measure = "indeg_unweighted",
                                       xlabel = "Outgoing interaction number",
                                       ylabel = "Incoming interaction number",
                                       color.use = cols)
## Signaling strength:
netAnalysis_diff_signalingRole_scatter(sec,
                                       x.measure = "outdeg",
                                       y.measure = "indeg",
                                       xlabel = "Outgoing interaction strength",
                                       ylabel = "Incoming interaction strength",
                                       color.use = cols)
# Positive values indicate the increase in the second dataset (ileum) while negative values indicate the increase in the first dataset (jejunum)

## Signaling number:
netAnalysis_diff_signalingRole_scatter(ecm,
                                       x.measure = "outdeg_unweighted",
                                       y.measure = "indeg_unweighted",
                                       xlabel = "Outgoing interaction number",
                                       ylabel = "Incoming interaction number",
                                       color.use = cols)
## Signaling strength:
netAnalysis_diff_signalingRole_scatter(ecm,
                                       x.measure = "outdeg",
                                       y.measure = "indeg",
                                       xlabel = "Outgoing interaction strength",
                                       ylabel = "Incoming interaction strength",
                                       color.use = cols)
# Positive values indicate the increase in the second dataset (ileum) while negative values indicate the increase in the first dataset (jejunum)

# View differential signaling in heatmap:
#red  (or blue) colored edges represent increased (or decreased) signaling in the ileum compared to the jejunum
# The top colored bar plot represents the sum of column of values displayed in the heatmap (incoming signaling). The right colored bar plot represents the sum of row of values (outgoing signaling). 
gg1 <- netVisual_heatmap(cell, color.use = cols) # number of interactions
gg2 <- netVisual_heatmap(cell, measure = "weight", color.use = cols) # strength of interactions
gg1 + gg2

gg1 <- netVisual_heatmap(sec, color.use = cols) # number of interactions
gg2 <- netVisual_heatmap(sec, measure = "weight", color.use = cols) # strength of interactions
gg1 + gg2

gg1 <- netVisual_heatmap(ecm, color.use = cols) # number of interactions
gg2 <- netVisual_heatmap(ecm, measure = "weight", color.use = cols) # strength of interactions
gg1 + gg2

# interaction number/strength based on broader cell lineage IDs
#levels(cell@idents)
group.cellType <- c(rep("B/ASC", 5), rep('T/ILC', 12), rep('Myeloid', 3), 
                    rep('Epithelial', 2), rep('Stromal', 2))
group.cellType <- factor(group.cellType, levels = c('B/ASC', 'T/ILC', 'Myeloid', 
                                                    'Epithelial', 'Stromal'))

object.list <- list(jejunum = cellj, ileum = celli)
object.list <- lapply(object.list, function(x) {mergeInteractions(x, group.cellType)})
cell <- mergeCellChat(object.list, add.names = names(object.list))
weight.max <- getMaxWeight(object.list, slot.name = c("idents", "net", "net"), attribute = c("idents","count", "count.merged"))
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_circle(object.list[[i]]@net$count.merged, weight.scale = T, label.edge= T, edge.weight.max = weight.max[3], edge.width.max = 12, title.name = paste0("Number of interactions - ", names(object.list)[i]))
}
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_circle(object.list[[i]]@net$weight.merged, weight.scale = T, label.edge= T, title.name = paste0("Strength of interactions - ", names(object.list)[i]))
}
par(mfrow = c(1,2), xpd=TRUE)
netVisual_diffInteraction(cell, weight.scale = T, measure = "count.merged", label.edge = T)
netVisual_diffInteraction(cell, weight.scale = T, measure = "weight.merged", label.edge = T)

object.list <- list(jejunum = secj, ileum = seci)
object.list <- lapply(object.list, function(x) {mergeInteractions(x, group.cellType)})
sec <- mergeCellChat(object.list, add.names = names(object.list))
weight.max <- getMaxWeight(object.list, slot.name = c("idents", "net", "net"), attribute = c("idents","count", "count.merged"))
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_circle(object.list[[i]]@net$count.merged, weight.scale = T, label.edge= T, edge.weight.max = weight.max[3], edge.width.max = 12, title.name = paste0("Number of interactions - ", names(object.list)[i]))
}
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_circle(object.list[[i]]@net$weight.merged, weight.scale = T, label.edge= T, title.name = paste0("Strength of interactions - ", names(object.list)[i]))
}
par(mfrow = c(1,2), xpd=TRUE)
netVisual_diffInteraction(sec, weight.scale = T, measure = "count.merged", label.edge = T)
netVisual_diffInteraction(sec, weight.scale = T, measure = "weight.merged", label.edge = T)

object.list <- list(jejunum = ecmj, ileum = ecmi)
object.list <- lapply(object.list, function(x) {mergeInteractions(x, group.cellType)})
ecm <- mergeCellChat(object.list, add.names = names(object.list))
weight.max <- getMaxWeight(object.list, slot.name = c("idents", "net", "net"), attribute = c("idents","count", "count.merged"))
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_circle(object.list[[i]]@net$count.merged, weight.scale = T, label.edge= T, edge.weight.max = weight.max[3], edge.width.max = 12, title.name = paste0("Number of interactions - ", names(object.list)[i]))
}
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_circle(object.list[[i]]@net$weight.merged, weight.scale = T, label.edge= T, title.name = paste0("Strength of interactions - ", names(object.list)[i]))
}
par(mfrow = c(1,2), xpd=TRUE)
netVisual_diffInteraction(ecm, weight.scale = T, measure = "count.merged", label.edge = T)
netVisual_diffInteraction(ecm, weight.scale = T, measure = "weight.merged", label.edge = T)

# Look at signaling changes in specific cell types:
# Positive values indicate the increase in the second dataset (ileum) while negative values indicate the increase in the first dataset (jejunum)
# Signal number:
gg1 <- netAnalysis_signalingChanges_scatter(cell, 
                                            idents.use = "Activated B cells",
                                            x.measure = "outdeg_unweighted",
                                            y.measure = "indeg_unweighted",
                                            xlabel = "Differential outgoing interaction number",
                                            ylabel = "Differential incoming interaction number",
                                            color.use = c('grey30', 'blue', 'red'))
# Signal strength:
gg2 <- netAnalysis_signalingChanges_scatter(cell, 
                                            idents.use = "Activated B cells",
                                            x.measure = "outdeg",
                                            y.measure = "indeg",
                                            xlabel = "Differential outgoing interaction strength",
                                            ylabel = "Differential incoming interaction strength",
                                            color.use = c('grey30', 'blue', 'red'))
patchwork::wrap_plots(plots = list(gg1,gg2))

gg1 <- netAnalysis_signalingChanges_scatter(sec, 
                                            idents.use = "Activated B cells",
                                            x.measure = "outdeg_unweighted",
                                            y.measure = "indeg_unweighted",
                                            xlabel = "Differential outgoing interaction number",
                                            ylabel = "Differential incoming interaction number",
                                            color.use = c('grey30', 'blue', 'red'))
gg2 <- netAnalysis_signalingChanges_scatter(sec, 
                                            idents.use = "Activated B cells",
                                            x.measure = "outdeg",
                                            y.measure = "indeg",
                                            xlabel = "Differential outgoing interaction strength",
                                            ylabel = "Differential incoming interaction strength",
                                            color.use = c('grey30', 'blue', 'red'))
patchwork::wrap_plots(plots = list(gg1,gg2))

gg1 <- netAnalysis_signalingChanges_scatter(ecm, 
                                            idents.use = "Activated B cells",
                                            x.measure = "outdeg_unweighted",
                                            y.measure = "indeg_unweighted",
                                            xlabel = "Differential outgoing interaction number",
                                            ylabel = "Differential incoming interaction number",
                                            color.use = c('grey30', 'blue', 'red'))
gg2 <- netAnalysis_signalingChanges_scatter(ecm, 
                                            idents.use = "Activated B cells",
                                            x.measure = "outdeg",
                                            y.measure = "indeg",
                                            xlabel = "Differential outgoing interaction strength",
                                            ylabel = "Differential incoming interaction strength",
                                            color.use = c('grey30', 'blue', 'red'))
patchwork::wrap_plots(plots = list(gg1,gg2))

gg1 <- netAnalysis_signalingChanges_scatter(cell, 
                                            idents.use = "Cycling B cells",
                                            x.measure = "outdeg_unweighted",
                                            y.measure = "indeg_unweighted",
                                            xlabel = "Differential outgoing interaction number",
                                            ylabel = "Differential incoming interaction number",
                                            color.use = c('grey30', 'blue', 'red'))
gg2 <- netAnalysis_signalingChanges_scatter(cell, 
                                            idents.use = "Cycling B cells",
                                            x.measure = "outdeg",
                                            y.measure = "indeg",
                                            xlabel = "Differential outgoing interaction strength",
                                            ylabel = "Differential incoming interaction strength",
                                            color.use = c('grey30', 'blue', 'red'))
patchwork::wrap_plots(plots = list(gg1,gg2))

gg1 <- netAnalysis_signalingChanges_scatter(sec, 
                                            idents.use = "Cycling B cells",
                                            x.measure = "outdeg_unweighted",
                                            y.measure = "indeg_unweighted",
                                            xlabel = "Differential outgoing interaction number",
                                            ylabel = "Differential incoming interaction number",
                                            color.use = c('grey30', 'blue', 'red'))
gg2 <- netAnalysis_signalingChanges_scatter(sec, 
                                            idents.use = "Cycling B cells",
                                            x.measure = "outdeg",
                                            y.measure = "indeg",
                                            xlabel = "Differential outgoing interaction strength",
                                            ylabel = "Differential incoming interaction strength",
                                            color.use = c('grey30', 'blue', 'red'))
patchwork::wrap_plots(plots = list(gg1,gg2))

gg1 <- netAnalysis_signalingChanges_scatter(ecm, 
                                            idents.use = "Cycling B cells",
                                            x.measure = "outdeg_unweighted",
                                            y.measure = "indeg_unweighted",
                                            xlabel = "Differential outgoing interaction number",
                                            ylabel = "Differential incoming interaction number",
                                            color.use = c('grey30', 'blue', 'red'))
gg2 <- netAnalysis_signalingChanges_scatter(ecm, 
                                            idents.use = "Cycling B cells",
                                            x.measure = "outdeg",
                                            y.measure = "indeg",
                                            xlabel = "Differential outgoing interaction strength",
                                            ylabel = "Differential incoming interaction strength",
                                            color.use = c('grey30', 'blue', 'red'))
patchwork::wrap_plots(plots = list(gg1,gg2))

gg1 <- netAnalysis_signalingChanges_scatter(cell, 
                                            idents.use = "Resting B cells",
                                            x.measure = "outdeg_unweighted",
                                            y.measure = "indeg_unweighted",
                                            xlabel = "Differential outgoing interaction number",
                                            ylabel = "Differential incoming interaction number",
                                            color.use = c('grey30', 'blue', 'red'))
gg2 <- netAnalysis_signalingChanges_scatter(cell, 
                                            idents.use = "Resting B cells",
                                            x.measure = "outdeg",
                                            y.measure = "indeg",
                                            xlabel = "Differential outgoing interaction strength",
                                            ylabel = "Differential incoming interaction strength",
                                            color.use = c('grey30', 'blue', 'red'))
patchwork::wrap_plots(plots = list(gg1,gg2))

gg1 <- netAnalysis_signalingChanges_scatter(sec, 
                                            idents.use = "Resting B cells",
                                            x.measure = "outdeg_unweighted",
                                            y.measure = "indeg_unweighted",
                                            xlabel = "Differential outgoing interaction number",
                                            ylabel = "Differential incoming interaction number",
                                            color.use = c('grey30', 'blue', 'red'))
gg2 <- netAnalysis_signalingChanges_scatter(sec, 
                                            idents.use = "Resting B cells",
                                            x.measure = "outdeg",
                                            y.measure = "indeg",
                                            xlabel = "Differential outgoing interaction strength",
                                            ylabel = "Differential incoming interaction strength",
                                            color.use = c('grey30', 'blue', 'red'))
patchwork::wrap_plots(plots = list(gg1,gg2))

gg1 <- netAnalysis_signalingChanges_scatter(ecm, 
                                            idents.use = "Resting B cells",
                                            x.measure = "outdeg_unweighted",
                                            y.measure = "indeg_unweighted",
                                            xlabel = "Differential outgoing interaction number",
                                            ylabel = "Differential incoming interaction number",
                                            color.use = c('grey30', 'blue', 'red'))
gg2 <- netAnalysis_signalingChanges_scatter(ecm, 
                                            idents.use = "Resting B cells",
                                            x.measure = "outdeg",
                                            y.measure = "indeg",
                                            xlabel = "Differential outgoing interaction strength",
                                            ylabel = "Differential incoming interaction strength",
                                            color.use = c('grey30', 'blue', 'red'))
patchwork::wrap_plots(plots = list(gg1,gg2))

gg1 <- netAnalysis_signalingChanges_scatter(cell, 
                                            idents.use = "Transitioning B cells",
                                            x.measure = "outdeg_unweighted",
                                            y.measure = "indeg_unweighted",
                                            xlabel = "Differential outgoing interaction number",
                                            ylabel = "Differential incoming interaction number",
                                            color.use = c('grey30', 'blue', 'red'))
gg2 <- netAnalysis_signalingChanges_scatter(cell, 
                                            idents.use = "Transitioning B cells",
                                            x.measure = "outdeg",
                                            y.measure = "indeg",
                                            xlabel = "Differential outgoing interaction strength",
                                            ylabel = "Differential incoming interaction strength",
                                            color.use = c('grey30', 'blue', 'red'))
patchwork::wrap_plots(plots = list(gg1,gg2))

gg1 <- netAnalysis_signalingChanges_scatter(sec, 
                                            idents.use = "Transitioning B cells",
                                            x.measure = "outdeg_unweighted",
                                            y.measure = "indeg_unweighted",
                                            xlabel = "Differential outgoing interaction number",
                                            ylabel = "Differential incoming interaction number",
                                            color.use = c('grey30', 'blue', 'red'))
gg2 <- netAnalysis_signalingChanges_scatter(sec, 
                                            idents.use = "Transitioning B cells",
                                            x.measure = "outdeg",
                                            y.measure = "indeg",
                                            xlabel = "Differential outgoing interaction strength",
                                            ylabel = "Differential incoming interaction strength",
                                            color.use = c('grey30', 'blue', 'red'))
patchwork::wrap_plots(plots = list(gg1,gg2))

gg1 <- netAnalysis_signalingChanges_scatter(ecm, 
                                            idents.use = "Transitioning B cells",
                                            x.measure = "outdeg_unweighted",
                                            y.measure = "indeg_unweighted",
                                            xlabel = "Differential outgoing interaction number",
                                            ylabel = "Differential incoming interaction number",
                                            color.use = c('grey30', 'blue', 'red'))
gg2 <- netAnalysis_signalingChanges_scatter(ecm, 
                                            idents.use = "Transitioning B cells",
                                            x.measure = "outdeg",
                                            y.measure = "indeg",
                                            xlabel = "Differential outgoing interaction strength",
                                            ylabel = "Differential incoming interaction strength",
                                            color.use = c('grey30', 'blue', 'red'))
patchwork::wrap_plots(plots = list(gg1,gg2))

gg1 <- netAnalysis_signalingChanges_scatter(cell, 
                                            idents.use = "Antibody-secreting cells",
                                            x.measure = "outdeg_unweighted",
                                            y.measure = "indeg_unweighted",
                                            xlabel = "Differential outgoing interaction number",
                                            ylabel = "Differential incoming interaction number",
                                            color.use = c('grey30', 'blue', 'red'))
gg2 <- netAnalysis_signalingChanges_scatter(cell, 
                                            idents.use = "Antibody-secreting cells",
                                            x.measure = "outdeg",
                                            y.measure = "indeg",
                                            xlabel = "Differential outgoing interaction strength",
                                            ylabel = "Differential incoming interaction strength",
                                            color.use = c('grey30', 'blue', 'red'))
patchwork::wrap_plots(plots = list(gg1,gg2))

gg1 <- netAnalysis_signalingChanges_scatter(sec, 
                                            idents.use = "Antibody-secreting cells",
                                            x.measure = "outdeg_unweighted",
                                            y.measure = "indeg_unweighted",
                                            xlabel = "Differential outgoing interaction number",
                                            ylabel = "Differential incoming interaction number",
                                            color.use = c('grey30', 'blue', 'red'))
gg2 <- netAnalysis_signalingChanges_scatter(sec, 
                                            idents.use = "Antibody-secreting cells",
                                            x.measure = "outdeg",
                                            y.measure = "indeg",
                                            xlabel = "Differential outgoing interaction strength",
                                            ylabel = "Differential incoming interaction strength",
                                            color.use = c('grey30', 'blue', 'red'))
patchwork::wrap_plots(plots = list(gg1,gg2))

gg1 <- netAnalysis_signalingChanges_scatter(ecm, 
                                            idents.use = "Antibody-secreting cells",
                                            x.measure = "outdeg_unweighted",
                                            y.measure = "indeg_unweighted",
                                            xlabel = "Differential outgoing interaction number",
                                            ylabel = "Differential incoming interaction number",
                                            color.use = c('grey30', 'blue', 'red'))
gg2 <- netAnalysis_signalingChanges_scatter(ecm, 
                                            idents.use = "Antibody-secreting cells",
                                            x.measure = "outdeg",
                                            y.measure = "indeg",
                                            xlabel = "Differential outgoing interaction strength",
                                            ylabel = "Differential incoming interaction strength",
                                            color.use = c('grey30', 'blue', 'red'))
patchwork::wrap_plots(plots = list(gg1,gg2))

# Identify enriched pathways (not cell type specific):
# red font = enriched in jejunum; blue font = enriched in ileum
gg1 <- rankNet(cell, mode = "comparison", stacked = T, do.stat = TRUE, color.use = c('blue', 'red'))
gg2 <- rankNet(cell, mode = "comparison", stacked = F, do.stat = TRUE, color.use = c('blue', 'red'))
gg1 + gg2

gg1 <- rankNet(sec, mode = "comparison", stacked = T, do.stat = TRUE, color.use = c('blue', 'red'))
gg2 <- rankNet(sec, mode = "comparison", stacked = F, do.stat = TRUE, color.use = c('blue', 'red'))
gg1 + gg2

gg1 <- rankNet(ecm, mode = "comparison", stacked = T, do.stat = TRUE, color.use = c('blue', 'red'))
gg2 <- rankNet(ecm, mode = "comparison", stacked = F, do.stat = TRUE, color.use = c('blue', 'red'))
gg1 + gg2

# Identify signalling affected by tissue, with B cells as targets
netVisual_bubble(cell, targets.use = 1,  comparison = c(1, 2), max.dataset = 2, title.name = "Increased signaling in Ileum", angle.x = 45, remove.isolate = T, color.heatmap = 'viridis', color.text = c('blue', 'red'))
netVisual_bubble(cell, targets.use = 1,  comparison = c(1, 2), max.dataset = 1, title.name = "Increased signaling in Jejunum", angle.x = 45, remove.isolate = T, color.heatmap = 'viridis', color.text = c('blue', 'red'))
netVisual_bubble(cell, targets.use = 2,  comparison = c(1, 2), max.dataset = 2, title.name = "Increased signaling in Ileum", angle.x = 45, remove.isolate = T, color.heatmap = 'viridis', color.text = c('blue', 'red'))
netVisual_bubble(cell, targets.use = 2,  comparison = c(1, 2), max.dataset = 1, title.name = "Increased signaling in Jejunum", angle.x = 45, remove.isolate = T, color.heatmap = 'viridis', color.text = c('blue', 'red'))
netVisual_bubble(cell, targets.use = 3,  comparison = c(1, 2), max.dataset = 2, title.name = "Increased signaling in Ileum", angle.x = 45, remove.isolate = T, color.heatmap = 'viridis', color.text = c('blue', 'red'))
netVisual_bubble(cell, targets.use = 3,  comparison = c(1, 2), max.dataset = 1, title.name = "Increased signaling in Jejunum", angle.x = 45, remove.isolate = T, color.heatmap = 'viridis', color.text = c('blue', 'red'))
netVisual_bubble(cell, targets.use = 4,  comparison = c(1, 2), max.dataset = 2, title.name = "Increased signaling in Ileum", angle.x = 45, remove.isolate = T, color.heatmap = 'viridis', color.text = c('blue', 'red'))
netVisual_bubble(cell, targets.use = 4,  comparison = c(1, 2), max.dataset = 1, title.name = "Increased signaling in Jejunum", angle.x = 45, remove.isolate = T, color.heatmap = 'viridis', color.text = c('blue', 'red'))
netVisual_bubble(cell, targets.use = 5,  comparison = c(1, 2), max.dataset = 2, title.name = "Increased signaling in Ileum", angle.x = 45, remove.isolate = T, color.heatmap = 'viridis', color.text = c('blue', 'red'))
netVisual_bubble(cell, targets.use = 5,  comparison = c(1, 2), max.dataset = 1, title.name = "Increased signaling in Jejunum", angle.x = 45, remove.isolate = T, color.heatmap = 'viridis', color.text = c('blue', 'red'))

netVisual_bubble(sec, targets.use = 1,  comparison = c(1, 2), max.dataset = 2, title.name = "Increased signaling in Ileum", angle.x = 45, remove.isolate = T, color.heatmap = 'viridis', color.text = c('blue', 'red'))
netVisual_bubble(sec, targets.use = 1,  comparison = c(1, 2), max.dataset = 1, title.name = "Increased signaling in Jejunum", angle.x = 45, remove.isolate = T, color.heatmap = 'viridis', color.text = c('blue', 'red'))
netVisual_bubble(sec, targets.use = 2,  comparison = c(1, 2), max.dataset = 2, title.name = "Increased signaling in Ileum", angle.x = 45, remove.isolate = T, color.heatmap = 'viridis', color.text = c('blue', 'red'))
netVisual_bubble(sec, targets.use = 2,  comparison = c(1, 2), max.dataset = 1, title.name = "Increased signaling in Jejunum", angle.x = 45, remove.isolate = T, color.heatmap = 'viridis', color.text = c('blue', 'red'))
netVisual_bubble(sec, targets.use = 3,  comparison = c(1, 2), max.dataset = 2, title.name = "Increased signaling in Ileum", angle.x = 45, remove.isolate = T, color.heatmap = 'viridis', color.text = c('blue', 'red'))
netVisual_bubble(sec, targets.use = 3,  comparison = c(1, 2), max.dataset = 1, title.name = "Increased signaling in Jejunum", angle.x = 45, remove.isolate = T, color.heatmap = 'viridis', color.text = c('blue', 'red'))
netVisual_bubble(sec, targets.use = 4,  comparison = c(1, 2), max.dataset = 2, title.name = "Increased signaling in Ileum", angle.x = 45, remove.isolate = T, color.heatmap = 'viridis', color.text = c('blue', 'red'))
netVisual_bubble(sec, targets.use = 4,  comparison = c(1, 2), max.dataset = 1, title.name = "Increased signaling in Jejunum", angle.x = 45, remove.isolate = T, color.heatmap = 'viridis', color.text = c('blue', 'red'))
netVisual_bubble(sec, targets.use = 5,  comparison = c(1, 2), max.dataset = 2, title.name = "Increased signaling in Ileum", angle.x = 45, remove.isolate = T, color.heatmap = 'viridis', color.text = c('blue', 'red'))
netVisual_bubble(sec, targets.use = 5,  comparison = c(1, 2), max.dataset = 1, title.name = "Increased signaling in Jejunum", angle.x = 45, remove.isolate = T, color.heatmap = 'viridis', color.text = c('blue', 'red'))

netVisual_bubble(ecm, targets.use = 1,  comparison = c(1, 2), max.dataset = 2, title.name = "Increased signaling in Ileum", angle.x = 45, remove.isolate = T, color.heatmap = 'viridis', color.text = c('blue', 'red'))
netVisual_bubble(ecm, targets.use = 1,  comparison = c(1, 2), max.dataset = 1, title.name = "Increased signaling in Jejunum", angle.x = 45, remove.isolate = T, color.heatmap = 'viridis', color.text = c('blue', 'red'))
netVisual_bubble(ecm, targets.use = 2,  comparison = c(1, 2), max.dataset = 2, title.name = "Increased signaling in Ileum", angle.x = 45, remove.isolate = T, color.heatmap = 'viridis', color.text = c('blue', 'red'))
netVisual_bubble(ecm, targets.use = 2,  comparison = c(1, 2), max.dataset = 1, title.name = "Increased signaling in Jejunum", angle.x = 45, remove.isolate = T, color.heatmap = 'viridis', color.text = c('blue', 'red'))
netVisual_bubble(ecm, targets.use = 3,  comparison = c(1, 2), max.dataset = 2, title.name = "Increased signaling in Ileum", angle.x = 45, remove.isolate = T, color.heatmap = 'viridis', color.text = c('blue', 'red'))
netVisual_bubble(ecm, targets.use = 3,  comparison = c(1, 2), max.dataset = 1, title.name = "Increased signaling in Jejunum", angle.x = 45, remove.isolate = T, color.heatmap = 'viridis', color.text = c('blue', 'red'))
netVisual_bubble(ecm, targets.use = 4,  comparison = c(1, 2), max.dataset = 2, title.name = "Increased signaling in Ileum", angle.x = 45, remove.isolate = T, color.heatmap = 'viridis', color.text = c('blue', 'red'))
netVisual_bubble(ecm, targets.use = 4,  comparison = c(1, 2), max.dataset = 1, title.name = "Increased signaling in Jejunum", angle.x = 45, remove.isolate = T, color.heatmap = 'viridis', color.text = c('blue', 'red'))
netVisual_bubble(ecm, targets.use = 5,  comparison = c(1, 2), max.dataset = 2, title.name = "Increased signaling in Ileum", angle.x = 45, remove.isolate = T, color.heatmap = 'viridis', color.text = c('blue', 'red'))
netVisual_bubble(ecm, targets.use = 5,  comparison = c(1, 2), max.dataset = 1, title.name = "Increased signaling in Jejunum", angle.x = 45, remove.isolate = T, color.heatmap = 'viridis', color.text = c('blue', 'red'))

# Identify signalling affected by tissue, with B cells as sources
netVisual_bubble(cell, sources.use = 1,  comparison = c(1, 2), max.dataset = 2, title.name = "Increased signaling in Ileum", angle.x = 45, remove.isolate = T, color.heatmap = 'viridis', color.text = c('blue', 'red'))
netVisual_bubble(cell, sources.use = 1,  comparison = c(1, 2), max.dataset = 1, title.name = "Increased signaling in Jejunum", angle.x = 45, remove.isolate = T, color.heatmap = 'viridis', color.text = c('blue', 'red'))
netVisual_bubble(cell, sources.use = 2,  comparison = c(1, 2), max.dataset = 2, title.name = "Increased signaling in Ileum", angle.x = 45, remove.isolate = T, color.heatmap = 'viridis', color.text = c('blue', 'red'))
netVisual_bubble(cell, sources.use = 2,  comparison = c(1, 2), max.dataset = 1, title.name = "Increased signaling in Jejunum", angle.x = 45, remove.isolate = T, color.heatmap = 'viridis', color.text = c('blue', 'red'))
netVisual_bubble(cell, sources.use = 3,  comparison = c(1, 2), max.dataset = 2, title.name = "Increased signaling in Ileum", angle.x = 45, remove.isolate = T, color.heatmap = 'viridis', color.text = c('blue', 'red'))
netVisual_bubble(cell, sources.use = 3,  comparison = c(1, 2), max.dataset = 1, title.name = "Increased signaling in Jejunum", angle.x = 45, remove.isolate = T, color.heatmap = 'viridis', color.text = c('blue', 'red'))
netVisual_bubble(cell, sources.use = 4,  comparison = c(1, 2), max.dataset = 2, title.name = "Increased signaling in Ileum", angle.x = 45, remove.isolate = T, color.heatmap = 'viridis', color.text = c('blue', 'red'))
netVisual_bubble(cell, sources.use = 4,  comparison = c(1, 2), max.dataset = 1, title.name = "Increased signaling in Jejunum", angle.x = 45, remove.isolate = T, color.heatmap = 'viridis', color.text = c('blue', 'red'))
netVisual_bubble(cell, sources.use = 5,  comparison = c(1, 2), max.dataset = 2, title.name = "Increased signaling in Ileum", angle.x = 45, remove.isolate = T, color.heatmap = 'viridis', color.text = c('blue', 'red'))
netVisual_bubble(cell, sources.use = 5,  comparison = c(1, 2), max.dataset = 1, title.name = "Increased signaling in Jejunum", angle.x = 45, remove.isolate = T, color.heatmap = 'viridis', color.text = c('blue', 'red'))

netVisual_bubble(sec, sources.use = 1,  comparison = c(1, 2), max.dataset = 2, title.name = "Increased signaling in Ileum", angle.x = 45, remove.isolate = T, color.heatmap = 'viridis', color.text = c('blue', 'red'))
netVisual_bubble(sec, sources.use = 1,  comparison = c(1, 2), max.dataset = 1, title.name = "Increased signaling in Jejunum", angle.x = 45, remove.isolate = T, color.heatmap = 'viridis', color.text = c('blue', 'red'))
netVisual_bubble(sec, sources.use = 2,  comparison = c(1, 2), max.dataset = 2, title.name = "Increased signaling in Ileum", angle.x = 45, remove.isolate = T, color.heatmap = 'viridis', color.text = c('blue', 'red'))
netVisual_bubble(sec, sources.use = 2,  comparison = c(1, 2), max.dataset = 1, title.name = "Increased signaling in Jejunum", angle.x = 45, remove.isolate = T, color.heatmap = 'viridis', color.text = c('blue', 'red'))
netVisual_bubble(sec, sources.use = 3,  comparison = c(1, 2), max.dataset = 2, title.name = "Increased signaling in Ileum", angle.x = 45, remove.isolate = T, color.heatmap = 'viridis', color.text = c('blue', 'red'))
netVisual_bubble(sec, sources.use = 3,  comparison = c(1, 2), max.dataset = 1, title.name = "Increased signaling in Jejunum", angle.x = 45, remove.isolate = T, color.heatmap = 'viridis', color.text = c('blue', 'red'))
netVisual_bubble(sec, sources.use = 4,  comparison = c(1, 2), max.dataset = 2, title.name = "Increased signaling in Ileum", angle.x = 45, remove.isolate = T, color.heatmap = 'viridis', color.text = c('blue', 'red'))
netVisual_bubble(sec, sources.use = 4,  comparison = c(1, 2), max.dataset = 1, title.name = "Increased signaling in Jejunum", angle.x = 45, remove.isolate = T, color.heatmap = 'viridis', color.text = c('blue', 'red'))
netVisual_bubble(sec, sources.use = 5,  comparison = c(1, 2), max.dataset = 2, title.name = "Increased signaling in Ileum", angle.x = 45, remove.isolate = T, color.heatmap = 'viridis', color.text = c('blue', 'red'))
netVisual_bubble(sec, sources.use = 5,  comparison = c(1, 2), max.dataset = 1, title.name = "Increased signaling in Jejunum", angle.x = 45, remove.isolate = T, color.heatmap = 'viridis', color.text = c('blue', 'red'))

netVisual_bubble(ecm, sources.use = 1,  comparison = c(1, 2), max.dataset = 2, title.name = "Increased signaling in Ileum", angle.x = 45, remove.isolate = T, color.heatmap = 'viridis', color.text = c('blue', 'red'))
netVisual_bubble(ecm, sources.use = 1,  comparison = c(1, 2), max.dataset = 1, title.name = "Increased signaling in Jejunum", angle.x = 45, remove.isolate = T, color.heatmap = 'viridis', color.text = c('blue', 'red'))
netVisual_bubble(ecm, sources.use = 2,  comparison = c(1, 2), max.dataset = 2, title.name = "Increased signaling in Ileum", angle.x = 45, remove.isolate = T, color.heatmap = 'viridis', color.text = c('blue', 'red'))
netVisual_bubble(ecm, sources.use = 2,  comparison = c(1, 2), max.dataset = 1, title.name = "Increased signaling in Jejunum", angle.x = 45, remove.isolate = T, color.heatmap = 'viridis', color.text = c('blue', 'red'))
netVisual_bubble(ecm, sources.use = 3,  comparison = c(1, 2), max.dataset = 2, title.name = "Increased signaling in Ileum", angle.x = 45, remove.isolate = T, color.heatmap = 'viridis', color.text = c('blue', 'red'))
netVisual_bubble(ecm, sources.use = 3,  comparison = c(1, 2), max.dataset = 1, title.name = "Increased signaling in Jejunum", angle.x = 45, remove.isolate = T, color.heatmap = 'viridis', color.text = c('blue', 'red'))
netVisual_bubble(ecm, sources.use = 4,  comparison = c(1, 2), max.dataset = 2, title.name = "Increased signaling in Ileum", angle.x = 45, remove.isolate = T, color.heatmap = 'viridis', color.text = c('blue', 'red'))
netVisual_bubble(ecm, sources.use = 4,  comparison = c(1, 2), max.dataset = 1, title.name = "Increased signaling in Jejunum", angle.x = 45, remove.isolate = T, color.heatmap = 'viridis', color.text = c('blue', 'red'))
netVisual_bubble(ecm, sources.use = 5,  comparison = c(1, 2), max.dataset = 2, title.name = "Increased signaling in Ileum", angle.x = 45, remove.isolate = T, color.heatmap = 'viridis', color.text = c('blue', 'red'))
netVisual_bubble(ecm, sources.use = 5,  comparison = c(1, 2), max.dataset = 1, title.name = "Increased signaling in Jejunum", angle.x = 45, remove.isolate = T, color.heatmap = 'viridis', color.text = c('blue', 'red'))

# More plots for MHC-II signaling----
netAnalysis_signalingRole_network(cellj, width = 8, height = 2.5, font.size = 10, color.heatmap = "BuGn", color.use = cols, signaling = 'MHC-II')
netAnalysis_signalingRole_network(celli, width = 8, height = 2.5, font.size = 10, color.heatmap = "BuGn", color.use = cols, signaling = 'MHC-II')

# extract signaling immportance using pieces of netAnalysis_signalingRole_network() source code:
# Plot difference in signaling role importance for jejunum vs ileum
centr0 <- cellj@netP$centr$`MHC-II`
mat <- matrix(unlist(centr0), ncol = length(centr0), byrow = FALSE)
mat <- t(mat)
rownames(mat) <- names(centr0); colnames(mat) <- names(centr0$outdeg)
mat <- sweep(mat, 1L, apply(mat, 1, max), '/', check.margin = FALSE)
matj <- t(mat)

centr0 <- celli@netP$centr$`MHC-II`
mat <- matrix(unlist(centr0), ncol = length(centr0), byrow = FALSE)
mat <- t(mat)
rownames(mat) <- names(centr0); colnames(mat) <- names(centr0$outdeg)
mat <- sweep(mat, 1L, apply(mat, 1, max), '/', check.margin = FALSE)
mati <- t(mat)

newmat <- as.data.frame(matj - mati) # calculate net difference in importance between jejunum and ileum
newmat$celltype <- rownames(newmat) 
newmat$Sender <- newmat$outdeg
newmat$Receiver <- newmat$indeg 
newmat$Mediator <- newmat$flowbet 
newmat$Influencer <- newmat$info 
newmat <- gather(newmat, Role, NetImportance, Sender:Influencer, factor_key=TRUE)
newmat$celltype <- factor(newmat$celltype, levels = rev(levels(cellj@idents)))

ggplot(newmat, aes(Role, celltype, fill = NetImportance)) +
  geom_tile() +
  scale_fill_gradientn(
    colors=c("red","white","blue"),
    values=rescale(c(-0.3,0,0.3)),
    limits=c(-0.3, 0.3), oob = squish) +
  theme_classic()
    
#Plot number and strength of MHC-II signaling in jejunum and ileum                                                
object.list <- list(jejunum = cellj, ileum = celli)
pathways.show <- c("MHC-II") 
par(mfrow = c(1,2), xpd=TRUE)
ht <- list()
for (i in 1:length(object.list)) {
  ht[[i]] <- netVisual_heatmap(object.list[[i]], signaling = pathways.show, color.use = cols,
                               color.heatmap = "BuGn",title.name = paste(pathways.show, "signaling ",names(object.list)[i]))
}
ComplexHeatmap::draw(ht[[1]] + ht[[2]], ht_gap = unit(0.5, "cm"))

# Plot MHC-II signaling from B cells
netVisual_aggregate(cellj, signaling = 'MHC-II', sources.use = c(1:4), color.use = cols, layout = 'chord')
netVisual_aggregate(celli, signaling = 'MHC-II', sources.use = c(1:4), color.use = cols, layout = 'chord')

# See receptor-ligand pairs for MHC-II singaling with B cells as senders
netVisual_bubble(cellj, sources.use = c(1:4), remove.isolate = FALSE, signaling = 'MHC-II', color.heatmap = 'viridis', min.quantile = 0, max.quantile = 1) + 
  scale_colour_viridis(option = 'viridis', limits=c(0, .4), oob = squish)
netVisual_bubble(celli, sources.use = c(1:4), remove.isolate = FALSE, signaling = 'MHC-II', color.heatmap = 'viridis', min.quantile = 0, max.quantile = 1) + 
  scale_colour_viridis(option = 'viridis', limits=c(0, .4), oob = squish)

# And show specific signaling between cell types for specified ligand receptor pairs
pairLR.MHC <- extractEnrichedLR(cellj, signaling = "MHC-II", geneLR.return = FALSE)
netVisual_chord_gene(cellj, sources.use = c(1:4), # only B cell sources
                     slot.name = 'net', pairLR.use = pairLR.MHC,
                     color.use = cols,
                     show.legend = TRUE, legend.pos.y = 5, legend.pos.x = 20)

pairLR.MHC <- extractEnrichedLR(celli, signaling = "MHC-II", geneLR.return = FALSE)
netVisual_chord_gene(celli, sources.use = c(1:4), # only B cell sources
                     slot.name = 'net', pairLR.use = pairLR.MHC,
                     color.use = cols,
                     show.legend = TRUE, legend.pos.y = 5, legend.pos.x = 20)

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
#  [1] grid      stats     graphics  grDevices utils     datasets  methods   base     

#other attached packages:
#  [1] viridis_0.6.2         viridisLite_0.4.1     ComplexHeatmap_2.12.1 scales_1.2.1          tidyr_1.2.1           uwot_0.1.14           Matrix_1.5-1          NMF_0.24.0            cluster_2.1.4         rngtools_1.5.2        pkgmaker_0.32.2      
#[12] registry_0.5-1        future_1.28.0         readxl_1.4.1          patchwork_1.1.2       SeuratDisk_0.0.0.9020 CellChat_1.6.1        Biobase_2.56.0        BiocGenerics_0.42.0   ggplot2_3.3.6         igraph_1.3.4          dplyr_1.0.10         
#[23] sp_1.5-0              SeuratObject_4.1.1    Seurat_4.1.1         

#loaded via a namespace (and not attached):
#  [1] utf8_1.2.2            reticulate_1.26       tidyselect_1.1.2      htmlwidgets_1.5.4     BiocParallel_1.30.3   Rtsne_0.16            munsell_0.5.0         codetools_0.2-18      ica_1.0-3             miniUI_0.1.1.1        withr_2.5.0          
#[12] spatstat.random_2.2-0 colorspace_2.0-3      progressr_0.11.0      knitr_1.40            ggalluvial_0.12.3     rstudioapi_0.14       stats4_4.2.2          ROCR_1.0-11           ggsignif_0.6.3        tensor_1.5            listenv_0.8.0        
#[23] polyclip_1.10-0       bit64_4.0.5           coda_0.19-4           parallelly_1.32.1     vctrs_0.4.1           generics_0.1.3        xfun_0.33             R6_2.5.1              doParallel_1.0.17     clue_0.3-61           hdf5r_1.3.5          
#[34] spatstat.utils_2.3-1  assertthat_0.2.1      promises_1.2.0.1      rgeos_0.5-9           gtable_0.3.1          globals_0.16.1        goftest_1.2-3         rlang_1.0.6           systemfonts_1.0.4     GlobalOptions_0.1.2   splines_4.2.2        
#[45] rstatix_0.7.0         lazyeval_0.2.2        spatstat.geom_2.4-0   broom_1.0.1           yaml_2.3.5            reshape2_1.4.4        abind_1.4-5           ggnetwork_0.5.10      backports_1.4.1       httpuv_1.6.6          tools_4.2.2          
#[56] gridBase_0.4-7        statnet.common_4.7.0  ellipsis_0.3.2        spatstat.core_2.4-4   RColorBrewer_1.1-3    ggridges_0.5.3        Rcpp_1.0.9            plyr_1.8.7            purrr_0.3.4           ggpubr_0.4.0          rpart_4.1.16         
#[67] deldir_1.0-6          pbapply_1.5-0         GetoptLong_1.0.5      cowplot_1.1.1         S4Vectors_0.34.0      zoo_1.8-10            ggrepel_0.9.1         magrittr_2.0.3        data.table_1.14.2     RSpectra_0.16-1       sna_2.7              
#[78] scattermore_0.8       circlize_0.4.15       lmtest_0.9-40         RANN_2.6.1            fitdistrplus_1.1-8    matrixStats_0.62.0    mime_0.12             evaluate_0.16         xtable_1.8-4          IRanges_2.30.1        gridExtra_2.3        
#[89] shape_1.4.6           compiler_4.2.2        tibble_3.1.8          KernSmooth_2.23-20    crayon_1.5.1          htmltools_0.5.3       mgcv_1.8-40           later_1.3.0           DBI_1.1.3             MASS_7.3-58.1         car_3.1-0            
#[100] cli_3.4.0             parallel_4.2.2        pkgconfig_2.0.3       plotly_4.10.0         spatstat.sparse_2.1-1 foreach_1.5.2         svglite_2.1.0         stringr_1.4.1         digest_0.6.29         sctransform_0.3.4     RcppAnnoy_0.0.19     
#[111] spatstat.data_2.2-0   rmarkdown_2.16        cellranger_1.1.0      leiden_0.4.3          shiny_1.7.2           rjson_0.2.21          lifecycle_1.0.2       nlme_3.1-159          jsonlite_1.8.0        carData_3.0-5         network_1.17.2       
#[122] BiocNeighbors_1.14.0  fansi_1.0.3           pillar_1.8.1          lattice_0.20-45       fastmap_1.1.0         httr_1.4.4            survival_3.4-0        glue_1.6.2            FNN_1.1.3.1           png_0.1-7             iterators_1.0.14     
#[133] bit_4.0.4             stringi_1.7.8         irlba_2.3.5           future.apply_1.9.1  