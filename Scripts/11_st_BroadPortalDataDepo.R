library(DropletUtils) 
library(Seurat) 
library(dplyr)
library(readxl)

# Make raw matrix (completely unfiltered and unnormalized data from all samples) ----
I2 <- Load10X_Spatial(data.dir = '/home/Jayne.Wiarda/SI_PP_SC_ST/SpaceRangerOutputs/I2', 
                      slice = 'I2')
I3 <- Load10X_Spatial(data.dir = '/home/Jayne.Wiarda/SI_PP_SC_ST/SpaceRangerOutputs/I3', 
                      slice = 'I3')
I4 <- Load10X_Spatial(data.dir = '/home/Jayne.Wiarda/SI_PP_SC_ST/SpaceRangerOutputs/I4', 
                      slice = 'I4')
J2 <- Load10X_Spatial(data.dir = '/home/Jayne.Wiarda/SI_PP_SC_ST/SpaceRangerOutputs/J2', 
                      slice = 'J2')
J3 <- Load10X_Spatial(data.dir = '/home/Jayne.Wiarda/SI_PP_SC_ST/SpaceRangerOutputs/J3', 
                      slice = 'J3')
J4 <- Load10X_Spatial(data.dir = '/home/Jayne.Wiarda/SI_PP_SC_ST/SpaceRangerOutputs/J4', 
                      slice = 'J4')
sc <- merge(I2, y = c(I3, I4, J2, J3, J4),
             add.cell.ids = c('I2', 'I3', 'I4', 'J2', 'J3', 'J4'))
write10xCounts(x = sc@assays$Spatial@counts, path = "/home/Jayne.Wiarda/SI_PP_SC_ST/DataDepo/Broad/RawMatrix_AllSamples_STomics", version = "3") 
rm(list = ls()) # clear all space

# Make processed matrix (the final matrix we used to analyze the data)
seu <- readRDS('/home/Jayne.Wiarda/SI_PP_SC_ST/AnnotatedSeurat/AllSamples_annotated_ST.rds')
write10xCounts(x = seu@assays$SCT@data, # save the normalized data, not the raw or scaled data
               path = "/home/Jayne.Wiarda/SI_PP_SC_ST/DataDepo/Broad/ProcessedMatrix_AllSamples_STomics", 
               version = "3") # we use SCT-normalized data since that is what we used for most analyses.... could instead opt for log-transformed Spatial data slot (seu@assays$Spatial@data) 

# Make meta data file following SCP formatting ----
meta <- as.data.frame(seu@meta.data)
df <- data.frame(rownames(meta)) # make data frame with cell barcode IDs
colnames(df) <- 'NAME' # name all columns according to SCP conventions
df$biosample_id <- meta$SampleID
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
df$library_preparation_protocol <- rep('EFO_0010961', nrow(df))
df$library_preparation_protocol__ontology_label <- rep("Visium Spatial Gene Expression", nrow(df))
df$sex <- rep('female', nrow(df))
df <- rbind(c('TYPE', rep('group', 11)), df)
write.table(df, file='/home/Jayne.Wiarda/SI_PP_SC_ST/DataDepo/Broad/Meta_SCP_STomics.tsv', 
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
clus <- subset(clus, select = -c(4,7,10,15)) # remove some meta data
colnames(clus)
clus <- clus[, c(1:3, 9, 10, 4:7, 11, 8)] # reorder meta data
colnames(clus)
clus$Region_ManAnn <- as.character(clus$Region_ManAnn)
clus$Region_Clust <- as.character(clus$Region_Clust)
clus$seurat_clusters <- as.character(clus$seurat_clusters)
id <- c('TYPE', 'numeric', 'numeric', 'group', 'group', rep('numeric', 4), 'group',
        'group')
clus <- rbind(id, clus)
colnames(clus) <- c('NAME', 'X', 'Y', 'Region_ManualAnnotation', 'Region_ClusteringAnnotation', 'nCount_Spatial', 
                    'nFeature_Spatial', 'nCount_SCT', 'nFeature_SCT', 'AnnotationConsensus', 'seurat_clusters')
write.table(clus, 
            file='/home/Jayne.Wiarda/SI_PP_SC_ST/DataDepo/Broad/AllSpots_STomics.tsv', 
            quote=FALSE, 
            sep='\t', 
            row.names = FALSE,
            col.names = TRUE)
rm(list = ls()) # clear all space

## Jejunum Only (clustering annotation):
seu <- readRDS('/home/Jayne.Wiarda/SI_PP_SC_ST/MappingPrediction/STomics_PredictedLocations_ClusteringAnnotation_Jejunum.rds')
umap <- as.data.frame(Embeddings(seu@reductions$umap))
meta <- as.data.frame(seu@meta.data)
clus <- data.frame(rownames(meta)) # make data frame with cell barcode IDs
colnames(clus) <- 'NAME' # name all columns according to SCP conventions
clus$X <- umap$UMAP_1
clus$Y <- umap$UMAP_2
clus <- cbind(clus, meta)
colnames(clus)
clus <- subset(clus, select = -c(4,7,10:12,14,15,50,51:52,57:73, 75:88,91:96,98:99,101:106,108:122)) # remove some meta data
colnames(clus)
clus <- clus[, c(1:3, 8, 4:7, 9:52)] # reorder meta data
colnames(clus)
colnames(clus) <- c(colnames(clus[,1:42]), 'Trajectory1_prediction.score.Increment2', 'Trajectory1_prediction.score.Increment3',
                    'Trajectory1_prediction.score.Increment5', 'Trajectory1_prediction.score.Increment4',
                    'Trajectory1_prediction.score.Increment1', 'Trajectory2_prediction.score.Increment2',
                    'Trajectory2_prediction.score.Increment3', 'Trajectory2_prediction.score.Increment4',
                    'Trajectory2_prediction.score.Increment1', 'Trajectory2_prediction.score.Increment5')
colnames(clus)
clus$Region_Clust <- as.character(clus$Region_Clust)
id <- c('TYPE', 'numeric', 'numeric', 'group', rep('numeric', 4), 'group',
        rep('numeric', 43))
clus <- rbind(id, clus)
write.table(clus, 
            file='/home/Jayne.Wiarda/SI_PP_SC_ST/DataDepo/Broad/Jejunum_ClusteringAnnotation_CellPredictions_STomics.tsv', 
            quote=FALSE, 
            sep='\t', 
            row.names = FALSE,
            col.names = TRUE)
rm(list = ls()) # clear all space

## Jejunum Only (manual annotation):
seu <- readRDS('/home/Jayne.Wiarda/SI_PP_SC_ST/MappingPrediction/STomics_PredictedLocations_ManualAnnotation_Jejunum.rds')
umap <- as.data.frame(Embeddings(seu@reductions$umap))
meta <- as.data.frame(seu@meta.data)
clus <- data.frame(rownames(meta)) # make data frame with cell barcode IDs
colnames(clus) <- 'NAME' # name all columns according to SCP conventions
clus$X <- umap$UMAP_1
clus$Y <- umap$UMAP_2
clus <- cbind(clus, meta)
colnames(clus)
clus <- subset(clus, select = -c(4,7,10, 11, 13,14,15,50,51:52,57:73, 75:88,91:96,98:99,101:106,108:122)) # remove some meta data
colnames(clus)
clus <- clus[, c(1:3, 8, 4:7, 9:52)] # reorder meta data
colnames(clus)
colnames(clus) <- c(colnames(clus[,1:42]), 'Trajectory1_prediction.score.Increment2', 'Trajectory1_prediction.score.Increment3',
                    'Trajectory1_prediction.score.Increment5', 'Trajectory1_prediction.score.Increment4',
                    'Trajectory1_prediction.score.Increment1', 'Trajectory2_prediction.score.Increment2',
                    'Trajectory2_prediction.score.Increment3', 'Trajectory2_prediction.score.Increment4',
                    'Trajectory2_prediction.score.Increment1', 'Trajectory2_prediction.score.Increment5')
colnames(clus)
clus$Region_ManAnn <- as.character(clus$Region_ManAnn)
id <- c('TYPE', 'numeric', 'numeric', 'group', rep('numeric', 4), 'group',
        rep('numeric', 43))
clus <- rbind(id, clus)
write.table(clus, 
            file='/home/Jayne.Wiarda/SI_PP_SC_ST/DataDepo/Broad/Jejunum_ManualAnnotation_CellPredictions_STomics.tsv', 
            quote=FALSE, 
            sep='\t', 
            row.names = FALSE,
            col.names = TRUE)
rm(list = ls()) # clear all space

## Ileum Only (clustering annotation):
seu <- readRDS('/home/Jayne.Wiarda/SI_PP_SC_ST/MappingPrediction/STomics_PredictedLocations_ClusteringAnnotation_Ileum.rds')
umap <- as.data.frame(Embeddings(seu@reductions$umap))
meta <- as.data.frame(seu@meta.data)
clus <- data.frame(rownames(meta)) # make data frame with cell barcode IDs
colnames(clus) <- 'NAME' # name all columns according to SCP conventions
clus$X <- umap$UMAP_1
clus$Y <- umap$UMAP_2
clus <- cbind(clus, meta)
colnames(clus)
clus <- subset(clus, select = -c(4,7,10:12,14,15,50,51, 54, 56, 58:63, 65:88, 90, 92:94, 96:98, 100:106, 108:122)) # remove some meta data
colnames(clus)
clus <- clus[, c(1:3, 8, 4:7, 9:52)] # reorder meta data
colnames(clus)
colnames(clus) <- c(colnames(clus[,1:42]), 'Trajectory1_prediction.score.Increment3', 'Trajectory1_prediction.score.Increment4',
                    'Trajectory1_prediction.score.Increment2', 'Trajectory1_prediction.score.Increment5',
                    'Trajectory1_prediction.score.Increment1', 'Trajectory2_prediction.score.Increment3',
                    'Trajectory2_prediction.score.Increment2', 'Trajectory2_prediction.score.Increment4',
                    'Trajectory2_prediction.score.Increment1', 'Trajectory2_prediction.score.Increment5')
colnames(clus)
clus$Region_Clust <- as.character(clus$Region_Clust)
id <- c('TYPE', 'numeric', 'numeric', 'group', rep('numeric', 4), 'group',
        rep('numeric', 43))
clus <- rbind(id, clus)
write.table(clus, 
            file='/home/Jayne.Wiarda/SI_PP_SC_ST/DataDepo/Broad/Ileum_ClusteringAnnotation_CellPredictions_STomics.tsv', 
            quote=FALSE, 
            sep='\t', 
            row.names = FALSE,
            col.names = TRUE)
rm(list = ls()) # clear all space

## Ileum Only (manual annotation):
seu <- readRDS('/home/Jayne.Wiarda/SI_PP_SC_ST/MappingPrediction/STomics_PredictedLocations_ManualAnnotation_Ileum.rds')
umap <- as.data.frame(Embeddings(seu@reductions$umap))
meta <- as.data.frame(seu@meta.data)
clus <- data.frame(rownames(meta)) # make data frame with cell barcode IDs
colnames(clus) <- 'NAME' # name all columns according to SCP conventions
clus$X <- umap$UMAP_1
clus$Y <- umap$UMAP_2
clus <- cbind(clus, meta)
colnames(clus)
clus <- subset(clus, select = -c(4,7,10:11, 13,14,15,50,51, 54, 56, 58:63, 65:88, 90, 92:94, 96:98, 100:106, 108:122)) # remove some meta data
colnames(clus)
clus <- clus[, c(1:3, 8, 4:7, 9:52)] # reorder meta data
colnames(clus)
colnames(clus) <- c(colnames(clus[,1:42]), 'Trajectory1_prediction.score.Increment3', 'Trajectory1_prediction.score.Increment4',
                    'Trajectory1_prediction.score.Increment2', 'Trajectory1_prediction.score.Increment5',
                    'Trajectory1_prediction.score.Increment1', 'Trajectory2_prediction.score.Increment3',
                    'Trajectory2_prediction.score.Increment2', 'Trajectory2_prediction.score.Increment4',
                    'Trajectory2_prediction.score.Increment1', 'Trajectory2_prediction.score.Increment5')
colnames(clus)
clus$Region_ManAnn <- as.character(clus$Region_ManAnn)
id <- c('TYPE', 'numeric', 'numeric', 'group', rep('numeric', 4), 'group',
        rep('numeric', 43))
clus <- rbind(id, clus)
write.table(clus, 
            file='/home/Jayne.Wiarda/SI_PP_SC_ST/DataDepo/Broad/Ileum_ManualAnnotation_CellPredictions_STomics.tsv', 
            quote=FALSE, 
            sep='\t', 
            row.names = FALSE,
            col.names = TRUE)
rm(list = ls()) # clear all space

## Create spatial files:
seu <- readRDS('/home/Jayne.Wiarda/SI_PP_SC_ST/AnnotatedSeurat/AllSamples_annotated_ST.rds')
meta <- as.data.frame(seu@meta.data)
colnames(meta)
meta <- subset(meta, select = -c(1, 7, 8, 12)) # remove some meta data

I2 <- data.frame(rownames(subset(meta, SampleID == 'I2')), seu@images$I2@coordinates$row, seu@images$I2@coordinates$col, subset(meta, SampleID == 'I2'))
colnames(I2)
colnames(I2) <- c('NAME', 'X', 'Y', colnames(I2[,4:11]))
colnames(I2)
I2 <- I2[, c(1:3, 9:10, 6, 4:5, 7:8, 11)] # reorder meta data
colnames(I2)
id <- c('TYPE', 'numeric', 'numeric', rep('group', 3), rep('numeric', 4), 'group')
I2$Region_Clust <- as.character(I2$Region_Clust)
I2$Region_ManAnn <- as.character(I2$Region_ManAnn)
I2 <- rbind(id, I2)
write.table(I2, 
            file='/home/Jayne.Wiarda/SI_PP_SC_ST/DataDepo/Broad/SpatialFile_I2.tsv', 
            quote=FALSE, 
            sep='\t', 
            row.names = FALSE,
            col.names = TRUE)

I3 <- data.frame(rownames(subset(meta, SampleID == 'I3')), seu@images$I3@coordinates$row, seu@images$I3@coordinates$col, subset(meta, SampleID == 'I3'))
colnames(I3)
colnames(I3) <- c('NAME', 'X', 'Y', colnames(I3[,4:11]))
colnames(I3)
I3 <- I3[, c(1:3, 9:10, 6, 4:5, 7:8, 11)] # reorder meta data
colnames(I3)
id <- c('TYPE', 'numeric', 'numeric', rep('group', 3), rep('numeric', 4), 'group')
I3$Region_Clust <- as.character(I3$Region_Clust)
I3$Region_ManAnn <- as.character(I3$Region_ManAnn)
I3 <- rbind(id, I3)
write.table(I3, 
            file='/home/Jayne.Wiarda/SI_PP_SC_ST/DataDepo/Broad/SpatialFile_I3.tsv', 
            quote=FALSE, 
            sep='\t', 
            row.names = FALSE,
            col.names = TRUE)

I4 <- data.frame(rownames(subset(meta, SampleID == 'I4')), seu@images$I4@coordinates$row, seu@images$I4@coordinates$col, subset(meta, SampleID == 'I4'))
colnames(I4)
colnames(I4) <- c('NAME', 'X', 'Y', colnames(I4[,4:11]))
colnames(I4)
I4 <- I4[, c(1:3, 9:10, 6, 4:5, 7:8, 11)] # reorder meta data
colnames(I4)
id <- c('TYPE', 'numeric', 'numeric', rep('group', 3), rep('numeric', 4), 'group')
I4$Region_Clust <- as.character(I4$Region_Clust)
I4$Region_ManAnn <- as.character(I4$Region_ManAnn)
I4 <- rbind(id, I4)
write.table(I4, 
            file='/home/Jayne.Wiarda/SI_PP_SC_ST/DataDepo/Broad/SpatialFile_I4.tsv', 
            quote=FALSE, 
            sep='\t', 
            row.names = FALSE,
            col.names = TRUE)

J2 <- data.frame(rownames(subset(meta, SampleID == 'J2')), seu@images$J2@coordinates$row, seu@images$J2@coordinates$col, subset(meta, SampleID == 'J2'))
colnames(J2)
colnames(J2) <- c('NAME', 'X', 'Y', colnames(J2[,4:11]))
colnames(J2)
J2 <- J2[, c(1:3, 9:10, 6, 4:5, 7:8, 11)] # reorder meta data
colnames(J2)
id <- c('TYPE', 'numeric', 'numeric', rep('group', 3), rep('numeric', 4), 'group')
J2$Region_Clust <- as.character(J2$Region_Clust)
J2$Region_ManAnn <- as.character(J2$Region_ManAnn)
J2 <- rbind(id, J2)
write.table(J2, 
            file='/home/Jayne.Wiarda/SI_PP_SC_ST/DataDepo/Broad/SpatialFile_J2.tsv', 
            quote=FALSE, 
            sep='\t', 
            row.names = FALSE,
            col.names = TRUE)

J3 <- data.frame(rownames(subset(meta, SampleID == 'J3')), seu@images$J3@coordinates$row, seu@images$J3@coordinates$col, subset(meta, SampleID == 'J3'))
colnames(J3)
colnames(J3) <- c('NAME', 'X', 'Y', colnames(J3[,4:11]))
colnames(J3)
J3 <- J3[, c(1:3, 9:10, 6, 4:5, 7:8, 11)] # reorder meta data
colnames(J3)
id <- c('TYPE', 'numeric', 'numeric', rep('group', 3), rep('numeric', 4), 'group')
J3$Region_Clust <- as.character(J3$Region_Clust)
J3$Region_ManAnn <- as.character(J3$Region_ManAnn)
J3 <- rbind(id, J3)
write.table(J3, 
            file='/home/Jayne.Wiarda/SI_PP_SC_ST/DataDepo/Broad/SpatialFile_J3.tsv', 
            quote=FALSE, 
            sep='\t', 
            row.names = FALSE,
            col.names = TRUE)

J4 <- data.frame(rownames(subset(meta, SampleID == 'J4')), seu@images$J4@coordinates$row, seu@images$J4@coordinates$col, subset(meta, SampleID == 'J4'))
colnames(J4)
colnames(J4) <- c('NAME', 'X', 'Y', colnames(J4[,4:11]))
colnames(J4)
J4 <- J4[, c(1:3, 9:10, 6, 4:5, 7:8, 11)] # reorder meta data
colnames(J4)
id <- c('TYPE', 'numeric', 'numeric', rep('group', 3), rep('numeric', 4), 'group')
J4$Region_Clust <- as.character(J4$Region_Clust)
J4$Region_ManAnn <- as.character(J4$Region_ManAnn)
J4 <- rbind(id, J4)
write.table(J4, 
            file='/home/Jayne.Wiarda/SI_PP_SC_ST/DataDepo/Broad/SpatialFile_J4.tsv', 
            quote=FALSE, 
            sep='\t', 
            row.names = FALSE,
            col.names = TRUE)
