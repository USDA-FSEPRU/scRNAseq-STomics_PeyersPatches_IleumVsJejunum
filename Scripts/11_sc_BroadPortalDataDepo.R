library(DropletUtils) 
library(Seurat) 
library(dplyr)
library(readxl)
library(SeuratDisk)

# Make raw matrix (completely unfiltered and unnormalized data from all samples) ----
data_dir <- c(Ileum_pig2 = "/home/Jayne.Wiarda/SI_PP_SC_ST/CellRangerOutputs/I2/raw_feature_bc_matrix",
              Ileum_pig3 = "/home/Jayne.Wiarda/SI_PP_SC_ST/CellRangerOutputs/I3/raw_feature_bc_matrix", 
             Ileum_pig4 = "/home/Jayne.Wiarda/SI_PP_SC_ST/CellRangerOutputs/I4/raw_feature_bc_matrix", 
              Jejunum_pig2 = "/home/Jayne.Wiarda/SI_PP_SC_ST/CellRangerOutputs/J2/raw_feature_bc_matrix", 
              Jejunum_pig3= "/home/Jayne.Wiarda/SI_PP_SC_ST/CellRangerOutputs/J3/raw_feature_bc_matrix", 
              Jejunum_pig4= "/home/Jayne.Wiarda/SI_PP_SC_ST/CellRangerOutputs/J4/raw_feature_bc_matrix")
lapply(data_dir, dir) # Should show barcodes.tsv.gz, features.tsv.gz, and matrix.mtx.gz for each sample listed
sc <- Read10X(data.dir = data_dir) # read the 10X data from all samples into a data matrix
sc = CreateSeuratObject(counts = sc) # create a Seurat object of the data matrix

write10xCounts(x = sc@assays$RNA@counts, path = "/home/Jayne.Wiarda/SI_PP_SC_ST/DataDepo/Broad/RawMatrix_AllSamples", version = "3") 
try <- Read10X(data.dir = "/home/Jayne.Wiarda/SI_PP_SC_ST/DataDepo/Broad/RawMatrix_AllSamples") # test upload of newly created matrix files
try <- CreateSeuratObject(counts = try) # make matrix into Seurat object again
identical(try@assays$RNA@counts, sc@assays$RNA@counts) # should return TRUE, which shows we have properly preserved the data using write10xCounts()
rm(list = ls()) # clear all space

# Make processed matrix (the final matrix we used to analyze the data)
seu <- LoadH5Seurat('/home/Jayne.Wiarda/SI_PP_SC_ST/AnnotatedSeurat/AllSamples_annotated.h5seurat')
write10xCounts(x = seu@assays$SCT@data, # save the normalized data, not the raw or scaled data
               path = "/home/Jayne.Wiarda/SI_PP_SC_ST/DataDepo/Broad/ProcessedMatrix_AllSamples", 
               version = "3") # we use SCT-normalized data since that is what we used for most analyses.... could instead opt for log-transformed RNA data slot (seu@assays$RNA@data) 

# Make meta data file following SCP formatting ----
meta <- as.data.frame(seu@meta.data)
df <- data.frame(rownames(meta)) # make data frame with cell barcode IDs
colnames(df) <- 'NAME' # name all columns according to SCP conventions
df$biosample_id <- meta$orig.ident
df$donor_id <- paste0('pig', substr(df$biosample_id, nchar(df$biosample_id) - 1 + 1, nchar(df$biosample_id)))
df$species <- rep('NCBITaxon_9823', nrow(df))
df$species__ontology_label <- rep('Sus scrofa', nrow(df))
df$disease <- rep('PATO_0000461', nrow(df))
df$disease__ontology_label <- rep('normal', nrow(df))
df$organ <- rep('UBERON_0002116', nrow(df))
df$organ__ontology_label <- rep('ileum', nrow(df))
df$organ[df$biosample_id=='J2'] <- 'UBERON_0002115'
df$organ__ontology_label[df$biosample_id=='J2'] <- 'jejunum'
df$organ[df$biosample_id=='J3'] <- 'UBERON_0002115'
df$organ__ontology_label[df$biosample_id=='J3'] <- 'jejunum'
df$organ[df$biosample_id=='J4'] <- 'UBERON_0002115'
df$organ__ontology_label[df$biosample_id=='J4'] <- 'jejunum'
df$library_preparation_protocol <- rep('EFO_0009922', nrow(df))
df$library_preparation_protocol__ontology_label <- rep("10x 3' v3", nrow(df))
df$sex <- rep('female', nrow(df))
df <- rbind(c('TYPE', rep('group', 11)), df)
write.table(df, file='/home/Jayne.Wiarda/SI_PP_SC_ST/DataDepo/Broad/Meta_SCP.tsv', 
            quote=FALSE, 
            sep='\t', 
            row.names = FALSE,
            col.names = TRUE)
rm(df) # clear space

# Make cluster file for all cell combinations used in paper ----

##All cells ----
umap <- as.data.frame(Embeddings(seu@reductions$umap))
#meta <- as.data.frame(seu@meta.data)
clus <- data.frame(rownames(meta)) # make data frame with cell barcode IDs
colnames(clus) <- 'NAME' # name all columns according to SCP conventions
clus$X <- umap$UMAP_1
clus$Y <- umap$UMAP_2
clus <- cbind(clus, meta)
colnames(clus)
clus <- subset(clus, select = -c(8, 37, 39, 50)) # remove some meta data
colnames(clus)
clus <- clus[, c(1:3, 46, 4:45, 47)] # reorder meta data
colnames(clus)
clus$celltype <- as.character(clus$celltype)
id <- c('TYPE', 'numeric', 'numeric', 'group', rep('numeric', 4), 'group',
        rep('numeric', 28), 'group', rep('numeric', 8), 'group')
clus <- rbind(id, clus)
write.table(clus, 
            file='/home/Jayne.Wiarda/SI_PP_SC_ST/DataDepo/Broad/AllCells.tsv', 
            quote=FALSE, 
            sep='\t', 
            row.names = FALSE,
            col.names = TRUE)
rm(list = ls()) # clear all space

## T cells & ILCs----
seu <- LoadH5Seurat('/home/Jayne.Wiarda/SI_PP_SC_ST/AnnotatedSeurat/TILC.h5seurat')
meta <- as.data.frame(seu@meta.data)
umap <- as.data.frame(Embeddings(seu@reductions$umap))
clus <- data.frame(rownames(meta)) # make data frame with cell barcode IDs
colnames(clus) <- 'NAME' # name all columns according to SCP conventions
clus$X <- umap$UMAP_1
clus$Y <- umap$UMAP_2
clus <- cbind(clus, meta)
colnames(clus)
clus <- subset(clus, select = -c(8, 37, 39:48, 50:51)) # remove some meta data
colnames(clus)
clus <- clus[, c(1:3, 37, 4:36)] # reorder meta data
clus$celltype <- as.character(clus$celltype)
colnames(clus)
id <- c('TYPE', 'numeric', 'numeric', 'group', rep('numeric', 4), 'group',
        rep('numeric', 28))
clus <- rbind(id, clus)
write.table(clus, 
            file='/home/Jayne.Wiarda/SI_PP_SC_ST/DataDepo/Broad/Tcells_ILCs.tsv', 
            quote=FALSE, 
            sep='\t', 
            row.names = FALSE,
            col.names = TRUE)
rm(list = ls()) # clear all space

## B cells and ASCs----
seu <- LoadH5Seurat('/home/Jayne.Wiarda/SI_PP_SC_ST/PseudotimeTrajectory/B_traj.h5seurat')
meta <- as.data.frame(seu@meta.data)
umap <- as.data.frame(Embeddings(seu@reductions$umap))
clus <- data.frame(rownames(meta)) # make data frame with cell barcode IDs
colnames(clus) <- 'NAME' # name all columns according to SCP conventions
clus$X <- umap$UMAP_1
clus$Y <- umap$UMAP_2
clus <- cbind(clus, meta)
colnames(clus)
clus <- subset(clus, select = -c(8, 37, 39:48, 50:52, 55:56, 62)) # remove some meta data
colnames(clus)
clus <- clus[, c(1:3, 37, 4:36, 38:44)] # reorder meta data
clus$celltype <- as.character(clus$celltype)
colnames(clus)
id <- c('TYPE', 'numeric', 'numeric', 'group', rep('numeric', 4), 'group',
        rep('numeric', 30), 'group', 'group', 'numeric', 'numeric', 'group')
clus <- rbind(id, clus)
write.table(clus, 
            file='/home/Jayne.Wiarda/SI_PP_SC_ST/DataDepo/Broad/Bcells.tsv', 
            quote=FALSE, 
            sep='\t', 
            row.names = FALSE,
            col.names = TRUE)
rm(list = ls()) # clear all space

## Myeloid cells----
seu <- LoadH5Seurat('/home/Jayne.Wiarda/SI_PP_SC_ST/AnnotatedSeurat/Myeloid.h5seurat')
meta <- as.data.frame(seu@meta.data)
umap <- as.data.frame(Embeddings(seu@reductions$umap))
clus <- data.frame(rownames(meta)) # make data frame with cell barcode IDs
colnames(clus) <- 'NAME' # name all columns according to SCP conventions
clus$X <- umap$UMAP_1
clus$Y <- umap$UMAP_2
clus <- cbind(clus, meta)
colnames(clus)
clus <- subset(clus, select = -c(8, 37, 39:48, 50:51)) # remove some meta data
colnames(clus)
clus <- clus[, c(1:3, 37, 4:36)] # reorder meta data
clus$celltype <- as.character(clus$celltype)
colnames(clus)
id <- c('TYPE', 'numeric', 'numeric', 'group', rep('numeric', 4), 'group',
        rep('numeric', 28))
clus <- rbind(id, clus)
write.table(clus, 
            file='/home/Jayne.Wiarda/SI_PP_SC_ST/DataDepo/Broad/Myeloid.tsv', 
            quote=FALSE, 
            sep='\t', 
            row.names = FALSE,
            col.names = TRUE)
rm(list = ls()) # clear all space

## Epithelial cells----
seu <- LoadH5Seurat('/home/Jayne.Wiarda/SI_PP_SC_ST/AnnotatedSeurat/Epithelial.h5seurat')
meta <- as.data.frame(seu@meta.data)
umap <- as.data.frame(Embeddings(seu@reductions$umap))
clus <- data.frame(rownames(meta)) # make data frame with cell barcode IDs
colnames(clus) <- 'NAME' # name all columns according to SCP conventions
clus$X <- umap$UMAP_1
clus$Y <- umap$UMAP_2
clus <- cbind(clus, meta)
colnames(clus)
clus <- subset(clus, select = -c(8, 37, 39, 50:51)) # remove some meta data
colnames(clus)
clus <- clus[, c(1:3, 46, 4:45)] # reorder meta data
clus$celltype <- as.character(clus$celltype)
colnames(clus)
id <- c('TYPE', 'numeric', 'numeric', 'group', rep('numeric', 4), 'group',
        rep('numeric', 28), 'group', rep('numeric', 8))
clus <- rbind(id, clus)
write.table(clus, 
            file='/home/Jayne.Wiarda/SI_PP_SC_ST/DataDepo/Broad/Epithelial.tsv', 
            quote=FALSE, 
            sep='\t', 
            row.names = FALSE,
            col.names = TRUE)
rm(list = ls()) # clear all space

## Stromal cells----
seu <- LoadH5Seurat('/home/Jayne.Wiarda/SI_PP_SC_ST/AnnotatedSeurat/StromalCells.h5seurat')
meta <- as.data.frame(seu@meta.data)
umap <- as.data.frame(Embeddings(seu@reductions$umap))
clus <- data.frame(rownames(meta)) # make data frame with cell barcode IDs
colnames(clus) <- 'NAME' # name all columns according to SCP conventions
clus$X <- umap$UMAP_1
clus$Y <- umap$UMAP_2
clus <- cbind(clus, meta)
colnames(clus)
clus <- subset(clus, select = -c(8, 37, 39:48, 50, 52)) # remove some meta data
colnames(clus)
clus <- clus[, c(1:3, 37, 4:36, 38)] # reorder meta data
clus$celltype <- as.character(clus$celltype)
clus$seurat_clusters <- as.character(clus$seurat_clusters)
colnames(clus)
id <- c('TYPE', 'numeric', 'numeric', 'group', rep('numeric', 4), 'group',
        rep('numeric', 28), 'group')
clus <- rbind(id, clus)
write.table(clus, 
            file='/home/Jayne.Wiarda/SI_PP_SC_ST/DataDepo/Broad/Stromal.tsv', 
            quote=FALSE, 
            sep='\t', 
            row.names = FALSE,
            col.names = TRUE)
rm(list = ls()) # clear all space

##Ileum Only ----
clus <- read.delim('/home/Jayne.Wiarda/SI_PP_SC_ST/DataDepo/Broad/AllCells.tsv')
clus <- clus %>% 
  filter(startsWith(as.character(NAME), 'I'))
id <- c('TYPE', 'numeric', 'numeric', 'group', rep('numeric', 4), 'group',
        rep('numeric', 28), 'group', rep('numeric', 8), 'group')
clus <- rbind(id, clus)
write.table(clus, 
            file='/home/Jayne.Wiarda/SI_PP_SC_ST/DataDepo/Broad/IleumOnly.tsv', 
            quote=FALSE, 
            sep='\t', 
            row.names = FALSE,
            col.names = TRUE)
rm(list = ls()) # clear all space

##Jejunum Only ----
clus <- read.delim('/home/Jayne.Wiarda/SI_PP_SC_ST/DataDepo/Broad/AllCells.tsv')
clus <- clus %>% 
  filter(startsWith(as.character(NAME), 'J'))
id <- c('TYPE', 'numeric', 'numeric', 'group', rep('numeric', 4), 'group',
        rep('numeric', 28), 'group', rep('numeric', 8), 'group')
clus <- rbind(id, clus)
write.table(clus, 
            file='/home/Jayne.Wiarda/SI_PP_SC_ST/DataDepo/Broad/JejunumOnly.tsv', 
            quote=FALSE, 
            sep='\t', 
            row.names = FALSE,
            col.names = TRUE)
rm(list = ls()) # clear all space
