### Load required software packages

library(Seurat)
library(SeuratObject)
library(SeuratDisk)
library(writexl)
library(dplyr)
library(ggplot2)
library(viridis)
library(scales)

# Run analysis on clustering annotations

## Load & process Seurat objects

#Load a processed Seurat object of spatial data (from both ileum + jejunum samples combined):
  
st <- readRDS('/home/Jayne.Wiarda/SI_PP_SC_ST/AnnotatedSeurat/AllSamples_annotated_ST.rds')
Idents(st) <- st$Region_Clust
st <- subset(st, idents = c('Follicle', 'Crypt', 'Villus', 'IFZ/PFZ')) # remove muscularis spots from spatial dataset since they aren't represented in scRNA-seq data (muscularis was removed from scRNA-seq samples)
st.list <- SplitObject(st, split.by = "SampleID") # split by sample IDs

st.list$I2@images$I3 <- NULL
st.list$I2@images$I4 <- NULL
st.list$I2@images$J2 <- NULL
st.list$I2@images$J3 <- NULL
st.list$I2@images$J4 <- NULL

st.list$I3@images$I2 <- NULL
st.list$I3@images$I4 <- NULL
st.list$I3@images$J2 <- NULL
st.list$I3@images$J3 <- NULL
st.list$I3@images$J4 <- NULL

st.list$I4@images$I2 <- NULL
st.list$I4@images$I3 <- NULL
st.list$I4@images$J2 <- NULL
st.list$I4@images$J3 <- NULL
st.list$I4@images$J4 <- NULL

st.list$J2@images$J3 <- NULL
st.list$J2@images$J4 <- NULL
st.list$J2@images$I2 <- NULL
st.list$J2@images$I3 <- NULL
st.list$J2@images$I4 <- NULL

st.list$J3@images$J2 <- NULL
st.list$J3@images$J4 <- NULL
st.list$J3@images$I2 <- NULL
st.list$J3@images$I3 <- NULL
st.list$J3@images$I4 <- NULL

st.list$J4@images$J2 <- NULL
st.list$J4@images$J3 <- NULL
st.list$J4@images$I2 <- NULL
st.list$J4@images$I3 <- NULL
st.list$J4@images$I4 <- NULL

st.list # see that each piece of list only has one associated image

#Normalize spatial data using SCT method and calculate PCs:
  
for (i in 1:length(st.list)) { # normalize data using SCTransform method
  st.list[[i]] <- SCTransform(st.list[[i]], 
                              return.only.var.genes = FALSE, 
                              verbose = TRUE,
                              assay = 'Spatial') 
} # use SCT normalization since this was also used to integrate scRNA-seq data

#Load scRNA-seq data:
#Data from ileum + jejunum: 
  
all <- LoadH5Seurat('/home/Jayne.Wiarda/SI_PP_SC_ST/AnnotatedSeurat/AllSamples_annotated.h5seurat')
Idents(all) <- all$celltype
all 

# Incorporate ranges for B cell trajectory
seu <- LoadH5Seurat('/home/Jayne.Wiarda/SI_PP_SC_ST/PseudotimeTrajectory/B_traj.h5seurat')

dif <- setdiff(colnames(all), colnames(seu)) # identify cell barcodes missing from the dataframe because they weren't predicted as B cells
psuedotime1_bin <- c(seu$psuedotime1_bin, rep('non-B', length(dif)))
psuedotime2_bin <- c(seu$psuedotime2_bin, rep('non-B', length(dif)))
CellBarcodes <- c(colnames(seu), dif)
meta <- data.frame(CellBarcodes, 
                   psuedotime1_bin,
                   psuedotime2_bin)
meta <- meta[ order(match(meta$CellBarcodes, colnames(all))), ]
identical(colnames(all), meta$CellBarcodes)
meta$celltype <- as.character(all$celltype)

meta <- meta %>% 
  mutate(psuedotime1_bin = ifelse(psuedotime1_bin == "non-B", celltype, psuedotime1_bin))
meta <- meta %>% 
  mutate(psuedotime2_bin = ifelse(psuedotime2_bin == "non-B", celltype, psuedotime2_bin))

all <- AddMetaData(all, metadata = c(meta))
rm(seu)

## Perform cell location prediction mapping in jejunum

#Load scRNA-seq data from only jejunum:

id <- data.frame(colnames(all), all$orig.ident, all$celltype, all$psuedotime1_bin, all$psuedotime2_bin)
colnames(id) <- c('barcode', 'SampleID', 'celltype', 'Btrajectory1', 'Btrajectory2')
id <- subset(id, SampleID == 'J2' | SampleID == 'J3' | SampleID == 'J4')

sc <- LoadH5Seurat('/home/Jayne.Wiarda/SI_PP_SC_ST/NormalizeIntegrateDimReduc_SConly/JejunumOnly.h5Seurat')
sc <- AddMetaData(sc, c(id))
DefaultAssay(sc) <- 'integrated' # use integrated since we have multiple samples, otherwise use SCT for a single sample

#Perform mapping prediction:

## Run this blurb for only one spatial slide (st)
# anchors <- FindTransferAnchors(reference = sc, query = st, reduction = 'cca', dims = 1:15,  normalization.method = 'SCT', recompute.residuals = TRUE)
# predictions <- TransferData(anchorset = anchors, refdata = list(cell_type = sc$celltype),  weight.reduction = 'cca', dims = 1:15)
# CellTypePredictions <- predictions

## Run this blurb for multiple spatial slides, as in st.list
CellTypePredictions <- list()
CellTrajectory1Predictions <- list()
CellTrajectory2Predictions <- list()
for(i in 4:6) { # 4:6 are jejunal slides
  anchors <- FindTransferAnchors(
    reference = sc,
    query = st.list[[i]],
    reduction = 'cca', # not using cca unless cross-modality comparison is being performed
    dims = 1:15, 
    normalization.method = 'SCT',
    recompute.residuals = TRUE) 
  predictions <- TransferData(anchorset = anchors, 
                              refdata = list(cell_type = sc$celltype, B_trajectory1 = sc$Btrajectory1, B_trajectory2 = sc$Btrajectory2), 
                              dims = 1:15,
                              weight.reduction = 'cca')
  CellTypePredictions[[i]] <- predictions$cell_type
  CellTrajectory1Predictions[[i]] <- predictions$B_trajectory1
  CellTrajectory2Predictions[[i]] <- predictions$B_trajectory2
}

CellTypePredictions <- do.call(rbind, CellTypePredictions)
CellTypePredictions <- as.data.frame(CellTypePredictions)
CellTypePredictions$CellBarcode <- rownames(CellTypePredictions)

CellTrajectory1Predictions <- do.call(rbind, CellTrajectory1Predictions)
CellTrajectory1Predictions <- as.data.frame(CellTrajectory1Predictions)
CellTrajectory1Predictions$CellBarcode <- rownames(CellTrajectory1Predictions)

CellTrajectory2Predictions <- do.call(rbind, CellTrajectory2Predictions)
CellTrajectory2Predictions <- as.data.frame(CellTrajectory2Predictions)
CellTrajectory2Predictions$CellBarcode <- rownames(CellTrajectory2Predictions)

colnames(CellTypePredictions) <- paste('CellType', colnames(CellTypePredictions), sep = "_")
colnames(CellTrajectory1Predictions) <- paste('Trajectory1', colnames(CellTrajectory1Predictions), sep = "_")
colnames(CellTrajectory2Predictions) <- paste('Trajectory2', colnames(CellTrajectory2Predictions), sep = "_")

write_xlsx(CellTypePredictions, '/home/Jayne.Wiarda/SI_PP_SC_ST/MappingPrediction/STomics_ClusteringAnnotation_CellTypePredictions_ClusteringAnnotation_Jejunum.xlsx')
write_xlsx(CellTrajectory1Predictions, '/home/Jayne.Wiarda/SI_PP_SC_ST/MappingPrediction/STomics_ClusteringAnnotation_CellTypePredictions_Btrajectory1Annotation_Jejunum.xlsx')
write_xlsx(CellTrajectory2Predictions, '/home/Jayne.Wiarda/SI_PP_SC_ST/MappingPrediction/STomics_ClusteringAnnotation_CellTypePredictions_Btrajectory2Annotation_Jejunum.xlsx')

#Incorporate predicted locations into spatial Seurat object:

Idents(st) <- st$SampleID
jej <- subset(st, idents = c('J2', 'J3', 'J4'))
jej <- AddMetaData(object = jej, 
                   metadata = c(CellTypePredictions, CellTrajectory1Predictions, CellTrajectory2Predictions))
saveRDS(jej, '/home/Jayne.Wiarda/SI_PP_SC_ST/MappingPrediction/STomics_PredictedLocations_ClusteringAnnotation_Jejunum.rds')

## Perform cell location prediction mapping in ileum

#Load scRNA-seq data from only ileum:

id <- data.frame(colnames(all), all$orig.ident, all$celltype, all$psuedotime1_bin, all$psuedotime2_bin)
colnames(id) <- c('barcode', 'SampleID', 'celltype', 'Btrajectory1', 'Btrajectory2')
id <- subset(id, SampleID == 'I2' | SampleID == 'I3' | SampleID == 'I4')

sc <- LoadH5Seurat('/home/Jayne.Wiarda/SI_PP_SC_ST/NormalizeIntegrateDimReduc_SConly/IleumOnly.h5Seurat')
sc <- AddMetaData(sc, c(id))
DefaultAssay(sc) <- 'integrated' # use integrated since we have multiple samples, otherwise use SCT for a single sample

#Perform mapping prediction:

## Run this blurb for only one spatial slide (st)
# anchors <- FindTransferAnchors(reference = sc, query = st, reduction = 'cca', dims = 1:15,  normalization.method = 'SCT', recompute.residuals = TRUE)
# predictions <- TransferData(anchorset = anchors, refdata = list(cell_type = sc$celltype),  weight.reduction = 'cca', dims = 1:15)
# CellTypePredictions <- predictions

## Run this blurb for multiple spatial slides, as in st.list
CellTypePredictions <- list()
CellTrajectory1Predictions <- list()
CellTrajectory2Predictions <- list()
for(i in 1:3) { # 1:3 are ileal slides
  anchors <- FindTransferAnchors(
    reference = sc,
    query = st.list[[i]],
    reduction = 'cca', # not using cca unless cross-modality comparison is being performed
    dims = 1:15, 
    normalization.method = 'SCT',
    recompute.residuals = TRUE) 
  predictions <- TransferData(anchorset = anchors, 
                              refdata = list(cell_type = sc$celltype, B_trajectory1 = sc$Btrajectory1, B_trajectory2 = sc$Btrajectory2), 
                              dims = 1:15,
                              weight.reduction = 'cca')
  CellTypePredictions[[i]] <- predictions$cell_type
  CellTrajectory1Predictions[[i]] <- predictions$B_trajectory1
  CellTrajectory2Predictions[[i]] <- predictions$B_trajectory2
}

CellTypePredictions <- do.call(rbind, CellTypePredictions)
CellTypePredictions <- as.data.frame(CellTypePredictions)
CellTypePredictions$CellBarcode <- rownames(CellTypePredictions)

CellTrajectory1Predictions <- do.call(rbind, CellTrajectory1Predictions)
CellTrajectory1Predictions <- as.data.frame(CellTrajectory1Predictions)
CellTrajectory1Predictions$CellBarcode <- rownames(CellTrajectory1Predictions)

CellTrajectory2Predictions <- do.call(rbind, CellTrajectory2Predictions)
CellTrajectory2Predictions <- as.data.frame(CellTrajectory2Predictions)
CellTrajectory2Predictions$CellBarcode <- rownames(CellTrajectory2Predictions)

colnames(CellTypePredictions) <- paste('CellType', colnames(CellTypePredictions), sep = "_")
colnames(CellTrajectory1Predictions) <- paste('Trajectory1', colnames(CellTrajectory1Predictions), sep = "_")
colnames(CellTrajectory2Predictions) <- paste('Trajectory2', colnames(CellTrajectory2Predictions), sep = "_")

write_xlsx(CellTypePredictions, '/home/Jayne.Wiarda/SI_PP_SC_ST/MappingPrediction/STomics_ClusteringAnnotation_CellTypePredictions_ClusteringAnnotation_Ileum.xlsx')
write_xlsx(CellTrajectory1Predictions, '/home/Jayne.Wiarda/SI_PP_SC_ST/MappingPrediction/STomics_ClusteringAnnotation_CellTypePredictions_Btrajectory1Annotation_Ileum.xlsx')
write_xlsx(CellTrajectory2Predictions, '/home/Jayne.Wiarda/SI_PP_SC_ST/MappingPrediction/STomics_ClusteringAnnotation_CellTypePredictions_Btrajectory2Annotation_Ileum.xlsx')

#Incorporate predicted locations into spatial Seurat object:

Idents(st) <- st$SampleID
il <- subset(st, idents = c('I2', 'I3', 'I4'))
il <- AddMetaData(object = il, 
                  metadata = c(CellTypePredictions, CellTrajectory1Predictions, CellTrajectory2Predictions))
saveRDS(il, '/home/Jayne.Wiarda/SI_PP_SC_ST/MappingPrediction/STomics_PredictedLocations_ClusteringAnnotation_Ileum.rds')

# Run analysis on manual annotations

## Load & process Seurat objects

#Load a processed Seurat object of spatial data (from both ileum + jejunum samples combined):

st <- readRDS('/home/Jayne.Wiarda/SI_PP_SC_ST/AnnotatedSeurat/AllSamples_annotated_ST.rds')
Idents(st) <- st$Region_ManAnn
st <- subset(st, idents = c('Follicle', 'Crypt', 'Villus', 'IFZ/PFZ')) # remove muscularis spots from spatial dataset since they aren't represented in scRNA-seq data (muscularis was removed from scRNA-seq samples)
st.list <- SplitObject(st, split.by = "SampleID") # split by sample IDs

st.list$I2@images$I3 <- NULL
st.list$I2@images$I4 <- NULL
st.list$I2@images$J2 <- NULL
st.list$I2@images$J3 <- NULL
st.list$I2@images$J4 <- NULL

st.list$I3@images$I2 <- NULL
st.list$I3@images$I4 <- NULL
st.list$I3@images$J2 <- NULL
st.list$I3@images$J3 <- NULL
st.list$I3@images$J4 <- NULL

st.list$I4@images$I2 <- NULL
st.list$I4@images$I3 <- NULL
st.list$I4@images$J2 <- NULL
st.list$I4@images$J3 <- NULL
st.list$I4@images$J4 <- NULL

st.list$J2@images$J3 <- NULL
st.list$J2@images$J4 <- NULL
st.list$J2@images$I2 <- NULL
st.list$J2@images$I3 <- NULL
st.list$J2@images$I4 <- NULL

st.list$J3@images$J2 <- NULL
st.list$J3@images$J4 <- NULL
st.list$J3@images$I2 <- NULL
st.list$J3@images$I3 <- NULL
st.list$J3@images$I4 <- NULL

st.list$J4@images$J2 <- NULL
st.list$J4@images$J3 <- NULL
st.list$J4@images$I2 <- NULL
st.list$J4@images$I3 <- NULL
st.list$J4@images$I4 <- NULL

st.list # see that each piece of list only has one associated image

#Normalize spatial data using SCT method and calculate PCs:

for (i in 1:length(st.list)) { # normalize data using SCTransform method
  st.list[[i]] <- SCTransform(st.list[[i]], 
                              return.only.var.genes = FALSE, 
                              verbose = TRUE,
                              assay = 'Spatial') 
} # use SCT normalization since this was also used to integrate scRNA-seq data

## Perform cell location prediction mapping in jejunum

#Load scRNA-seq data from only jejunum:

id <- data.frame(colnames(all), all$orig.ident, all$celltype, all$psuedotime1_bin, all$psuedotime2_bin)
colnames(id) <- c('barcode', 'SampleID', 'celltype', 'Btrajectory1', 'Btrajectory2')
id <- subset(id, SampleID == 'J2' | SampleID == 'J3' | SampleID == 'J4')

sc <- LoadH5Seurat('/home/Jayne.Wiarda/SI_PP_SC_ST/NormalizeIntegrateDimReduc_SConly/JejunumOnly.h5Seurat')
sc <- AddMetaData(sc, c(id))
DefaultAssay(sc) <- 'integrated' # use integrated since we have multiple samples, otherwise use SCT for a single sample

#Perform mapping prediction:

## Run this blurb for only one spatial slide (st)
# anchors <- FindTransferAnchors(reference = sc, query = st, reduction = 'cca', dims = 1:15,  normalization.method = 'SCT', recompute.residuals = TRUE)
# predictions <- TransferData(anchorset = anchors, refdata = list(cell_type = sc$celltype),  weight.reduction = 'cca', dims = 1:15)
# CellTypePredictions <- predictions

## Run this blurb for multiple spatial slides, as in st.list
CellTypePredictions <- list()
CellTrajectory1Predictions <- list()
CellTrajectory2Predictions <- list()
for(i in 4:6) { # 4:6 are jejunal slides
  anchors <- FindTransferAnchors(
    reference = sc,
    query = st.list[[i]],
    reduction = 'cca', # not using cca unless cross-modality comparison is being performed
    dims = 1:15, 
    normalization.method = 'SCT',
    recompute.residuals = TRUE) 
  predictions <- TransferData(anchorset = anchors, 
                              refdata = list(cell_type = sc$celltype, B_trajectory1 = sc$Btrajectory1, B_trajectory2 = sc$Btrajectory2), 
                              dims = 1:15,
                              weight.reduction = 'cca')
  CellTypePredictions[[i]] <- predictions$cell_type
  CellTrajectory1Predictions[[i]] <- predictions$B_trajectory1
  CellTrajectory2Predictions[[i]] <- predictions$B_trajectory2
}

CellTypePredictions <- do.call(rbind, CellTypePredictions)
CellTypePredictions <- as.data.frame(CellTypePredictions)
CellTypePredictions$CellBarcode <- rownames(CellTypePredictions)

CellTrajectory1Predictions <- do.call(rbind, CellTrajectory1Predictions)
CellTrajectory1Predictions <- as.data.frame(CellTrajectory1Predictions)
CellTrajectory1Predictions$CellBarcode <- rownames(CellTrajectory1Predictions)

CellTrajectory2Predictions <- do.call(rbind, CellTrajectory2Predictions)
CellTrajectory2Predictions <- as.data.frame(CellTrajectory2Predictions)
CellTrajectory2Predictions$CellBarcode <- rownames(CellTrajectory2Predictions)

colnames(CellTypePredictions) <- paste('CellType', colnames(CellTypePredictions), sep = "_")
colnames(CellTrajectory1Predictions) <- paste('Trajectory1', colnames(CellTrajectory1Predictions), sep = "_")
colnames(CellTrajectory2Predictions) <- paste('Trajectory2', colnames(CellTrajectory2Predictions), sep = "_")

write_xlsx(CellTypePredictions, '/home/Jayne.Wiarda/SI_PP_SC_ST/MappingPrediction/STomics_ManualAnnotation_CellTypePredictions_ClusteringAnnotation_Jejunum.xlsx')
write_xlsx(CellTrajectory1Predictions, '/home/Jayne.Wiarda/SI_PP_SC_ST/MappingPrediction/STomics_ManualAnnotation_CellTypePredictions_Btrajectory1Annotation_Jejunum.xlsx')
write_xlsx(CellTrajectory2Predictions, '/home/Jayne.Wiarda/SI_PP_SC_ST/MappingPrediction/STomics_ManualAnnotation_CellTypePredictions_Btrajectory2Annotation_Jejunum.xlsx')

#Incorporate predicted locations into spatial Seurat object:

Idents(st) <- st$SampleID
jej <- subset(st, idents = c('J2', 'J3', 'J4'))
jej <- AddMetaData(object = jej, 
                   metadata = c(CellTypePredictions, CellTrajectory1Predictions, CellTrajectory2Predictions))
saveRDS(jej, '/home/Jayne.Wiarda/SI_PP_SC_ST/MappingPrediction/STomics_PredictedLocations_ManualAnnotation_Jejunum.rds')

## Perform cell location prediction mapping in ileum

#Load scRNA-seq data from only ileum:

id <- data.frame(colnames(all), all$orig.ident, all$celltype, all$psuedotime1_bin, all$psuedotime2_bin)
colnames(id) <- c('barcode', 'SampleID', 'celltype', 'Btrajectory1', 'Btrajectory2')
id <- subset(id, SampleID == 'I2' | SampleID == 'I3' | SampleID == 'I4')

sc <- LoadH5Seurat('/home/Jayne.Wiarda/SI_PP_SC_ST/NormalizeIntegrateDimReduc_SConly/IleumOnly.h5Seurat')
sc <- AddMetaData(sc, c(id))
DefaultAssay(sc) <- 'integrated' # use integrated since we have multiple samples, otherwise use SCT for a single sample

#Perform mapping prediction:

## Run this blurb for only one spatial slide (st)
# anchors <- FindTransferAnchors(reference = sc, query = st, reduction = 'cca', dims = 1:15,  normalization.method = 'SCT', recompute.residuals = TRUE)
# predictions <- TransferData(anchorset = anchors, refdata = list(cell_type = sc$celltype),  weight.reduction = 'cca', dims = 1:15)
# CellTypePredictions <- predictions

## Run this blurb for multiple spatial slides, as in st.list
CellTypePredictions <- list()
CellTrajectory1Predictions <- list()
CellTrajectory2Predictions <- list()
for(i in 1:3) { # 1:3 are ileal slides
  anchors <- FindTransferAnchors(
    reference = sc,
    query = st.list[[i]],
    reduction = 'cca', # not using cca unless cross-modality comparison is being performed
    dims = 1:15, 
    normalization.method = 'SCT',
    recompute.residuals = TRUE) 
  predictions <- TransferData(anchorset = anchors, 
                              refdata = list(cell_type = sc$celltype, B_trajectory1 = sc$Btrajectory1, B_trajectory2 = sc$Btrajectory2), 
                              dims = 1:15,
                              weight.reduction = 'cca')
  CellTypePredictions[[i]] <- predictions$cell_type
  CellTrajectory1Predictions[[i]] <- predictions$B_trajectory1
  CellTrajectory2Predictions[[i]] <- predictions$B_trajectory2
}

CellTypePredictions <- do.call(rbind, CellTypePredictions)
CellTypePredictions <- as.data.frame(CellTypePredictions)
CellTypePredictions$CellBarcode <- rownames(CellTypePredictions)

CellTrajectory1Predictions <- do.call(rbind, CellTrajectory1Predictions)
CellTrajectory1Predictions <- as.data.frame(CellTrajectory1Predictions)
CellTrajectory1Predictions$CellBarcode <- rownames(CellTrajectory1Predictions)

CellTrajectory2Predictions <- do.call(rbind, CellTrajectory2Predictions)
CellTrajectory2Predictions <- as.data.frame(CellTrajectory2Predictions)
CellTrajectory2Predictions$CellBarcode <- rownames(CellTrajectory2Predictions)

colnames(CellTypePredictions) <- paste('CellType', colnames(CellTypePredictions), sep = "_")
colnames(CellTrajectory1Predictions) <- paste('Trajectory1', colnames(CellTrajectory1Predictions), sep = "_")
colnames(CellTrajectory2Predictions) <- paste('Trajectory2', colnames(CellTrajectory2Predictions), sep = "_")

write_xlsx(CellTypePredictions, '/home/Jayne.Wiarda/SI_PP_SC_ST/MappingPrediction/STomics_ManualAnnotation_CellTypePredictions_ClusteringAnnotation_Ileum.xlsx')
write_xlsx(CellTrajectory1Predictions, '/home/Jayne.Wiarda/SI_PP_SC_ST/MappingPrediction/STomics_ManualAnnotation_CellTypePredictions_Btrajectory1Annotation_Ileum.xlsx')
write_xlsx(CellTrajectory2Predictions, '/home/Jayne.Wiarda/SI_PP_SC_ST/MappingPrediction/STomics_ManualAnnotation_CellTypePredictions_Btrajectory2Annotation_Ileum.xlsx')

#Incorporate predicted locations into spatial Seurat object:

Idents(st) <- st$SampleID
il <- subset(st, idents = c('I2', 'I3', 'I4'))
il <- AddMetaData(object = il, 
                  metadata = c(CellTypePredictions, CellTrajectory1Predictions, CellTrajectory2Predictions))
saveRDS(il, '/home/Jayne.Wiarda/SI_PP_SC_ST/MappingPrediction/STomics_PredictedLocations_ManualAnnotation_Ileum.rds')

## Next, visualize results:

### Just showing results for the clustering annotation overlaid onto tissue sections but can do same with manual annotation as well

seu <- readRDS('/home/Jayne.Wiarda/SI_PP_SC_ST/MappingPrediction/STomics_PredictedLocations_ClusteringAnnotation_Jejunum.rds')
SpatialFeaturePlot(seu, alpha = c(.25, .8), crop = TRUE,                   
                   features = c('CellType_prediction.score.Activated.B.cells', 
                                'CellType_prediction.score.Resting.B.cells', 
                                'CellType_prediction.score.Cycling.B.cells', 
                                'CellType_prediction.score.Transitioning.B.cells', 
                                'CellType_prediction.score.Antibody.secreting.cells'), 
                   ncol = 3, 
                   images = 'J4') & NoAxes() & 
  scale_fill_viridis(                       
    begin = 0.3, end = 0.95,                     
    direction = 1,                       
    discrete = FALSE,                       
    option = "plasma", 
    #limits = c(0.0,1),                   
    oob = squish)

SpatialFeaturePlot(seu, alpha = c(.25, .8), crop = TRUE,                   
                   features = c('Trajectory1_prediction.score..0.33.5.', 
                                'Trajectory1_prediction.score..33.5.67.',
                                'Trajectory1_prediction.score..67.101.',
                                'Trajectory1_prediction.score..101.134.',
                                'Trajectory1_prediction.score..134.168.'), 
                   ncol = 3, 
                   images = 'J4') & NoAxes() & 
  scale_fill_viridis(                       
    begin = 0.3, end = 0.95,                       
    direction = 1,                       
    discrete = FALSE,                       
    option = "plasma", 
    limits = c(0.0,0.5),                         
    oob = squish)

SpatialFeaturePlot(seu, alpha = c(.25, .8), crop = TRUE,                   
                   features = c('Trajectory2_prediction.score..0.54.5.', 
                                'Trajectory2_prediction.score..54.5.109.', 
                                'Trajectory2_prediction.score..109.164.',
                                'Trajectory2_prediction.score..164.218.',
                                'Trajectory2_prediction.score..218.273.'), 
                   ncol = 3, 
                   images = 'J4') & NoAxes() & 
  scale_fill_viridis(                       
    begin = 0.3, end = 0.95,                           
    direction = 1,                       
    discrete = FALSE,                       
    option = "plasma", 
    limits = c(0.0,0.4),                      
    oob = squish)

seu <- readRDS('/home/Jayne.Wiarda/SI_PP_SC_ST/MappingPrediction/STomics_PredictedLocations_ClusteringAnnotation_Ileum.rds')
SpatialFeaturePlot(seu, alpha = c(.25, .8), crop = TRUE,                   
                   features = c('CellType_prediction.score.Activated.B.cells', 
                                'CellType_prediction.score.Resting.B.cells', 
                                'CellType_prediction.score.Cycling.B.cells', 
                                'CellType_prediction.score.Transitioning.B.cells', 
                                'CellType_prediction.score.Antibody.secreting.cells'), 
                   ncol = 3, 
                   images = 'I4') & NoAxes() & 
  scale_fill_viridis(                       
    begin = 0.3, end = 0.95,                              
    direction = 1,                       
    discrete = FALSE,                       
    option = "plasma", 
    #limits = c(0.0,0.4),                       
    oob = squish)

SpatialFeaturePlot(seu, alpha = c(.25, .8), crop = TRUE,                   
                   features = c('Trajectory1_prediction.score..0.33.5.', 
                                'Trajectory1_prediction.score..33.5.67.',
                                'Trajectory1_prediction.score..67.101.',
                                'Trajectory1_prediction.score..101.134.',
                                'Trajectory1_prediction.score..134.168.'), 
                   ncol = 3, 
                   images = 'I4') & NoAxes() & 
  scale_fill_viridis(                       
    begin = 0.3, end = 0.95,                            
    direction = 1,                       
    discrete = FALSE,                       
    option = "plasma", 
    limits = c(0.0,0.5),                     
    oob = squish)

SpatialFeaturePlot(seu, alpha = c(.25, .8), crop = TRUE,                   
                   features = c('Trajectory2_prediction.score..0.54.5.', 
                                'Trajectory2_prediction.score..54.5.109.', 
                                'Trajectory2_prediction.score..109.164.',
                                'Trajectory2_prediction.score..164.218.',
                                'Trajectory2_prediction.score..218.273.'), 
                   ncol = 3, 
                   images = 'I4') & NoAxes() & 
  scale_fill_viridis(                       
    begin = 0.3, end = 0.95,                             
    direction = 1,                       
    discrete = FALSE,                       
    option = "plasma", 
    limits = c(0.0,0.4),                      
    oob = squish)

### We can also make bar plots to show results ----

mapj <- read_xlsx('/home/Jayne.Wiarda/SI_PP_SC_ST/MappingPrediction/STomics_ClusteringAnnotation_CellTypePredictions_Btrajectory1Annotation_Jejunum.xlsx')
mapi <- read_xlsx('/home/Jayne.Wiarda/SI_PP_SC_ST/MappingPrediction/STomics_ClusteringAnnotation_CellTypePredictions_Btrajectory1Annotation_Ileum.xlsx')
mapClus <- as.data.frame(rbind(mapj, mapi))

mapj <- read_xlsx('/home/Jayne.Wiarda/SI_PP_SC_ST/MappingPrediction/STomics_ManualAnnotation_CellTypePredictions_Btrajectory1Annotation_Jejunum.xlsx')
mapi <- read_xlsx('/home/Jayne.Wiarda/SI_PP_SC_ST/MappingPrediction/STomics_ManualAnnotation_CellTypePredictions_Btrajectory1Annotation_Ileum.xlsx')
mapMan <- as.data.frame(rbind(mapj, mapi))

st <- readRDS('/home/Jayne.Wiarda/SI_PP_SC_ST/AnnotatedSeurat/AllSamples_annotated_ST.rds')
st$tissue <- substr(st$SampleID, 1, 1)   
Idents(st) <- st$Region_Clust
stClus <- subset(st, idents = c('Follicle', 'Crypt', 'Villus', 'IFZ/PFZ')) # remove muscularis spots from spatial dataset since they aren't represented in scRNA-seq data (muscularis was removed from scRNA-seq samples)
Idents(st) <- st$Region_ManAnn
stMan <- subset(st, idents = c('Follicle', 'Crypt', 'Villus', 'IFZ/PFZ')) # remove muscularis spots from spatial dataset since they aren't represented in scRNA-seq data (muscularis was removed from scRNA-seq samples)

order <- colnames(stClus)
mapClus <- mapClus[match(order, mapClus$Trajectory1_CellBarcode), ]
rownames(mapClus) <- mapClus$Trajectory1_CellBarcode
identical(colnames(stClus), rownames(mapClus))

order <- colnames(stMan)
mapMan <- mapMan[match(order, mapMan$Trajectory1_CellBarcode), ]
rownames(mapMan) <- mapMan$Trajectory1_CellBarcode
identical(colnames(stMan), rownames(mapMan))

stClus <- AddMetaData(stClus, mapClus)
stMan <- AddMetaData(stMan, mapMan)

Idents(stClus) <- stClus$tissue
stClusIl <- subset(stClus, idents = 'I')
stClusIl$Region_Clust <- factor(stClusIl$Region_Clust, levels = c('Villus', 'Crypt', 'IFZ/PFZ', 'Follicle'))
Idents(stClusIl) <- 'Region_Clust'

Idents(stClus) <- stClus$tissue
stClusJej <- subset(stClus, idents = 'J')
stClusJej$Region_Clust <- factor(stClusJej$Region_Clust, levels = c('Villus', 'Crypt', 'IFZ/PFZ', 'Follicle'))
Idents(stClusJej) <- 'Region_Clust'

Idents(stMan) <- stMan$tissue
stManIl <- subset(stMan, idents = 'I')
stManIl$Region_ManAnn <- factor(stManIl$Region_ManAnn, levels = c('Villus', 'Crypt', 'IFZ/PFZ', 'Follicle'))
Idents(stManIl) <- 'Region_ManAnn'

Idents(stMan) <- stMan$tissue
stManJej <- subset(stMan, idents = 'J')
stManJej$Region_ManAnn <- factor(stManJej$Region_ManAnn, levels = c('Villus', 'Crypt', 'IFZ/PFZ', 'Follicle'))
Idents(stManJej) <- 'Region_ManAnn'

sub1 <- as.data.frame(stClusJej@meta.data)
colnames(sub1)
sub1 <- sub1[,c(10,16:19,37)]
colnames(sub1) <- c('Region', 'Increment2', 'Increment3', 'Increment5', 'Increment4', 'Increment1')
sub1 <- gather(sub1, key="Increment", value="predictionScore", 2:6)
sub1$Annotation <- rep('Clustering', nrow(sub1))
sub1$Region <- factor(sub1$Region, levels = c('Villus', 'Crypt', 'IFZ/PFZ', 'Follicle'))

sub2 <- as.data.frame(stManJej@meta.data)
colnames(sub2)
sub2 <- sub2[,c(9,16:19,37)]
colnames(sub2) <- c('Region', 'Increment2', 'Increment3', 'Increment5', 'Increment4', 'Increment1')
sub2 <- gather(sub2, key="Increment", value="predictionScore", 2:6)
sub2$Annotation <- rep('Manual', nrow(sub2))
sub2$Region <- factor(sub2$Region, levels = c('Villus', 'Crypt', 'IFZ/PFZ', 'Follicle'))

sub <- rbind(sub1, sub2)

# Show both clustering and manual annotation side-by-side:
dat <- subset(sub, Increment == 'Increment1')
g1 <- ggplot(dat, aes(x=Region,
                y=predictionScore,
                fill = Annotation)) + 
  theme_bw() +
  geom_boxplot(aes(fill=Annotation), width=0.6, position = position_dodge(width=.8), 
               coef = 0, color = 'black', lwd = 0.6, outlier.shape = NA) +
  geom_point(shape=16, position = position_jitterdodge(jitter.width=.15, dodge.width = 0.8), size = .1, color = 'grey20') + 
  scale_fill_manual(breaks = dat$Annotation, values=c("grey85", "grey60")) + 
  stat_summary(fun.y="mean", geom="point", size=2.5,
               position=position_dodge(width=0.8), 
               color="red")+
  theme(axis.text.x=element_text(size=rel(1.75)), 
        axis.text.y=element_text(size=rel(1.75)))+ 
  labs(x="Region",y="Prediction score",
       title="Trajectory 1 - Increment 1 (Jejunum)") +
  theme(axis.title.x = element_text(size = 16),
        axis.text.x=element_text(size=12), 
        axis.title.y = element_text(size = 16),
        axis.text.y=element_text(size=12),
        legend.title = element_text (size=16),
        title = element_text(size = 16)) +
  theme(axis.text.x=element_text(size=rel(1.75)), axis.text.y=element_text(size=rel(1.75))) +
  scale_y_continuous(expand = c(0,0), lim = c(0,max(dat$predictionScore)))

dat <- subset(sub, Increment == 'Increment2')
g2 <- ggplot(dat, aes(x=Region,
                y=predictionScore,
                fill = Annotation)) + 
  theme_bw() +
  geom_boxplot(aes(fill=Annotation), width=0.6, position = position_dodge(width=.8), 
               coef = 0, color = 'black', lwd = 0.6, outlier.shape = NA) +
  geom_point(shape=16, position = position_jitterdodge(jitter.width=.15, dodge.width = 0.8), size = .1, color = 'grey20') + 
  scale_fill_manual(breaks = dat$Annotation, values=c("grey85", "grey60")) + 
  stat_summary(fun.y="mean", geom="point", size=2.5,
               position=position_dodge(width=0.8), 
               color="red")+
  theme(axis.text.x=element_text(size=rel(1.75)), 
        axis.text.y=element_text(size=rel(1.75)))+ 
  labs(x="Region",y="Prediction score",
       title="Trajectory 1 - Increment 2 (Jejunum)") +
  theme(axis.title.x = element_text(size = 16),
        axis.text.x=element_text(size=12), 
        axis.title.y = element_text(size = 16),
        axis.text.y=element_text(size=12),
        legend.title = element_text (size=16),
        title = element_text(size = 16)) +
  theme(axis.text.x=element_text(size=rel(1.75)), axis.text.y=element_text(size=rel(1.75))) +
  scale_y_continuous(expand = c(0,0), lim = c(0,max(dat$predictionScore)))

dat <- subset(sub, Increment == 'Increment3')
g3 <- ggplot(dat, aes(x=Region,
                y=predictionScore,
                fill = Annotation)) + 
  theme_bw() +
  geom_boxplot(aes(fill=Annotation), width=0.6, position = position_dodge(width=.8), 
               coef = 0, color = 'black', lwd = 0.6, outlier.shape = NA) +
  geom_point(shape=16, position = position_jitterdodge(jitter.width=.15, dodge.width = 0.8), size = .1, color = 'grey20') + 
  scale_fill_manual(breaks = dat$Annotation, values=c("grey85", "grey60")) + 
  stat_summary(fun.y="mean", geom="point", size=2.5,
               position=position_dodge(width=0.8), 
               color="red")+
  theme(axis.text.x=element_text(size=rel(1.75)), 
        axis.text.y=element_text(size=rel(1.75)))+ 
  labs(x="Region",y="Prediction score",
       title="Trajectory 1 - Increment 3 (Jejunum)") +
  theme(axis.title.x = element_text(size = 16),
        axis.text.x=element_text(size=12), 
        axis.title.y = element_text(size = 16),
        axis.text.y=element_text(size=12),
        legend.title = element_text (size=16),
        title = element_text(size = 16)) +
  theme(axis.text.x=element_text(size=rel(1.75)), axis.text.y=element_text(size=rel(1.75))) +
  scale_y_continuous(expand = c(0,0), lim = c(0,max(dat$predictionScore)))

dat <- subset(sub, Increment == 'Increment4')
g4 <- ggplot(dat, aes(x=Region,
                y=predictionScore,
                fill = Annotation)) + 
  theme_bw() +
  geom_boxplot(aes(fill=Annotation), width=0.6, position = position_dodge(width=.8), 
               coef = 0, color = 'black', lwd = 0.6, outlier.shape = NA) +
  geom_point(shape=16, position = position_jitterdodge(jitter.width=.15, dodge.width = 0.8), size = .1, color = 'grey20') + 
  scale_fill_manual(breaks = dat$Annotation, values=c("grey85", "grey60")) + 
  stat_summary(fun.y="mean", geom="point", size=2.5,
               position=position_dodge(width=0.8), 
               color="red")+
  theme(axis.text.x=element_text(size=rel(1.75)), 
        axis.text.y=element_text(size=rel(1.75)))+ 
  labs(x="Region",y="Prediction score",
       title="Trajectory 1 - Increment 4 (Jejunum)") +
  theme(axis.title.x = element_text(size = 16),
        axis.text.x=element_text(size=12), 
        axis.title.y = element_text(size = 16),
        axis.text.y=element_text(size=12),
        legend.title = element_text (size=16),
        title = element_text(size = 16)) +
  theme(axis.text.x=element_text(size=rel(1.75)), axis.text.y=element_text(size=rel(1.75))) +
  scale_y_continuous(expand = c(0,0), lim = c(0,max(dat$predictionScore)))

dat <- subset(sub, Increment == 'Increment5')
g5 <- ggplot(dat, aes(x=Region,
                y=predictionScore,
                fill = Annotation)) + 
  theme_bw() +
  geom_boxplot(aes(fill=Annotation), width=0.6, position = position_dodge(width=.8), 
               coef = 0, color = 'black', lwd = 0.6, outlier.shape = NA) +
  geom_point(shape=16, position = position_jitterdodge(jitter.width=.15, dodge.width = 0.8), size = .1, color = 'grey20') + 
  scale_fill_manual(breaks = dat$Annotation, values=c("grey85", "grey60")) + 
  stat_summary(fun.y="mean", geom="point", size=2.5,
               position=position_dodge(width=0.8), 
               color="red")+
  theme(axis.text.x=element_text(size=rel(1.75)), 
        axis.text.y=element_text(size=rel(1.75)))+ 
  labs(x="Region",y="Prediction score",
       title="Trajectory 1 - Increment 5 (Jejunum)") +
  theme(axis.title.x = element_text(size = 16),
        axis.text.x=element_text(size=12), 
        axis.title.y = element_text(size = 16),
        axis.text.y=element_text(size=12),
        legend.title = element_text (size=16),
        title = element_text(size = 16)) +
  theme(axis.text.x=element_text(size=rel(1.75)), axis.text.y=element_text(size=rel(1.75))) +
  scale_y_continuous(expand = c(0,0), lim = c(0,max(dat$predictionScore)))

g1+g2+g3+g4+g5

sub1 <- as.data.frame(stClusIl@meta.data)
colnames(sub1)
sub1 <- sub1[,c(10,16:19,37)]
colnames(sub1) <- c('Region', 'Increment2', 'Increment3', 'Increment5', 'Increment4', 'Increment1')
sub1 <- gather(sub1, key="Increment", value="predictionScore", 2:6)
sub1$Annotation <- rep('Clustering', nrow(sub1))
sub1$Region <- factor(sub1$Region, levels = c('Villus', 'Crypt', 'IFZ/PFZ', 'Follicle'))

sub2 <- as.data.frame(stManIl@meta.data)
colnames(sub2)
sub2 <- sub2[,c(9,16:19,37)]
colnames(sub2) <- c('Region', 'Increment2', 'Increment3', 'Increment5', 'Increment4', 'Increment1')
sub2 <- gather(sub2, key="Increment", value="predictionScore", 2:6)
sub2$Annotation <- rep('Manual', nrow(sub2))
sub2$Region <- factor(sub2$Region, levels = c('Villus', 'Crypt', 'IFZ/PFZ', 'Follicle'))

sub <- rbind(sub1, sub2)

dat <- subset(sub, Increment == 'Increment1')
g1 <- ggplot(dat, aes(x=Region,
                y=predictionScore,
                fill = Annotation)) + 
  theme_bw() +
  geom_boxplot(aes(fill=Annotation), width=0.6, position = position_dodge(width=.8), 
               coef = 0, color = 'black', lwd = 0.6, outlier.shape = NA) +
  geom_point(shape=16, position = position_jitterdodge(jitter.width=.15, dodge.width = 0.8), size = .1, color = 'grey20') + 
  scale_fill_manual(breaks = dat$Annotation, values=c("grey85", "grey60")) + 
  stat_summary(fun.y="mean", geom="point", size=2.5,
               position=position_dodge(width=0.8), 
               color="red")+
  theme(axis.text.x=element_text(size=rel(1.75)), 
        axis.text.y=element_text(size=rel(1.75)))+ 
  labs(x="Region",y="Prediction score",
       title="Trajectory 1 - Increment 1 (Ileum)") +
  theme(axis.title.x = element_text(size = 16),
        axis.text.x=element_text(size=12), 
        axis.title.y = element_text(size = 16),
        axis.text.y=element_text(size=12),
        legend.title = element_text (size=16),
        title = element_text(size = 16)) +
  theme(axis.text.x=element_text(size=rel(1.75)), axis.text.y=element_text(size=rel(1.75))) +
  scale_y_continuous(expand = c(0,0), lim = c(0,max(dat$predictionScore)))

dat <- subset(sub, Increment == 'Increment2')
g2 <- ggplot(dat, aes(x=Region,
                y=predictionScore,
                fill = Annotation)) + 
  theme_bw() +
  geom_boxplot(aes(fill=Annotation), width=0.6, position = position_dodge(width=.8), 
               coef = 0, color = 'black', lwd = 0.6, outlier.shape = NA) +
  geom_point(shape=16, position = position_jitterdodge(jitter.width=.15, dodge.width = 0.8), size = .1, color = 'grey20') + 
  scale_fill_manual(breaks = dat$Annotation, values=c("grey85", "grey60")) + 
  stat_summary(fun.y="mean", geom="point", size=2.5,
               position=position_dodge(width=0.8), 
               color="red")+
  theme(axis.text.x=element_text(size=rel(1.75)), 
        axis.text.y=element_text(size=rel(1.75)))+ 
  labs(x="Region",y="Prediction score",
       title="Trajectory 1 - Increment 2 (Ileum)") +
  theme(axis.title.x = element_text(size = 16),
        axis.text.x=element_text(size=12), 
        axis.title.y = element_text(size = 16),
        axis.text.y=element_text(size=12),
        legend.title = element_text (size=16),
        title = element_text(size = 16)) +
  theme(axis.text.x=element_text(size=rel(1.75)), axis.text.y=element_text(size=rel(1.75))) +
  scale_y_continuous(expand = c(0,0), lim = c(0,max(dat$predictionScore)))

dat <- subset(sub, Increment == 'Increment3')
g3 <- ggplot(dat, aes(x=Region,
                y=predictionScore,
                fill = Annotation)) + 
  theme_bw() +
  geom_boxplot(aes(fill=Annotation), width=0.6, position = position_dodge(width=.8), 
               coef = 0, color = 'black', lwd = 0.6, outlier.shape = NA) +
  geom_point(shape=16, position = position_jitterdodge(jitter.width=.15, dodge.width = 0.8), size = .1, color = 'grey20') + 
  scale_fill_manual(breaks = dat$Annotation, values=c("grey85", "grey60")) + 
  stat_summary(fun.y="mean", geom="point", size=2.5,
               position=position_dodge(width=0.8), 
               color="red")+
  theme(axis.text.x=element_text(size=rel(1.75)), 
        axis.text.y=element_text(size=rel(1.75)))+ 
  labs(x="Region",y="Prediction score",
       title="Trajectory 1 - Increment 3 (Ileum)") +
  theme(axis.title.x = element_text(size = 16),
        axis.text.x=element_text(size=12), 
        axis.title.y = element_text(size = 16),
        axis.text.y=element_text(size=12),
        legend.title = element_text (size=16),
        title = element_text(size = 16)) +
  theme(axis.text.x=element_text(size=rel(1.75)), axis.text.y=element_text(size=rel(1.75))) +
  scale_y_continuous(expand = c(0,0), lim = c(0,max(dat$predictionScore)))

dat <- subset(sub, Increment == 'Increment4')
g4 <- ggplot(dat, aes(x=Region,
                y=predictionScore,
                fill = Annotation)) + 
  theme_bw() +
  geom_boxplot(aes(fill=Annotation), width=0.6, position = position_dodge(width=.8), 
               coef = 0, color = 'black', lwd = 0.6, outlier.shape = NA) +
  geom_point(shape=16, position = position_jitterdodge(jitter.width=.15, dodge.width = 0.8), size = .1, color = 'grey20') + 
  scale_fill_manual(breaks = dat$Annotation, values=c("grey85", "grey60")) + 
  stat_summary(fun.y="mean", geom="point", size=2.5,
               position=position_dodge(width=0.8), 
               color="red")+
  theme(axis.text.x=element_text(size=rel(1.75)), 
        axis.text.y=element_text(size=rel(1.75)))+ 
  labs(x="Region",y="Prediction score",
       title="Trajectory 1 - Increment 4 (Ileum)") +
  theme(axis.title.x = element_text(size = 16),
        axis.text.x=element_text(size=12), 
        axis.title.y = element_text(size = 16),
        axis.text.y=element_text(size=12),
        legend.title = element_text (size=16),
        title = element_text(size = 16)) +
  theme(axis.text.x=element_text(size=rel(1.75)), axis.text.y=element_text(size=rel(1.75))) +
  scale_y_continuous(expand = c(0,0), lim = c(0,max(dat$predictionScore)))

dat <- subset(sub, Increment == 'Increment5')
g5 <- ggplot(dat, aes(x=Region,
                y=predictionScore,
                fill = Annotation)) + 
  theme_bw() +
  geom_boxplot(aes(fill=Annotation), width=0.6, position = position_dodge(width=.8), 
               coef = 0, color = 'black', lwd = 0.6, outlier.shape = NA) +
  geom_point(shape=16, position = position_jitterdodge(jitter.width=.15, dodge.width = 0.8), size = .1, color = 'grey20') + 
  scale_fill_manual(breaks = dat$Annotation, values=c("grey85", "grey60")) + 
  stat_summary(fun.y="mean", geom="point", size=2.5,
               position=position_dodge(width=0.8), 
               color="red")+
  theme(axis.text.x=element_text(size=rel(1.75)), 
        axis.text.y=element_text(size=rel(1.75)))+ 
  labs(x="Region",y="Prediction score",
       title="Trajectory 1 - Increment 5 (Ileum)") +
  theme(axis.title.x = element_text(size = 16),
        axis.text.x=element_text(size=12), 
        axis.title.y = element_text(size = 16),
        axis.text.y=element_text(size=12),
        legend.title = element_text (size=16),
        title = element_text(size = 16)) +
  theme(axis.text.x=element_text(size=rel(1.75)), axis.text.y=element_text(size=rel(1.75))) +
  scale_y_continuous(expand = c(0,0), lim = c(0,max(dat$predictionScore)))

g1+g2+g3+g4+g5

# Trajectory 2:

mapj <- read_xlsx('/home/Jayne.Wiarda/SI_PP_SC_ST/MappingPrediction/STomics_ClusteringAnnotation_CellTypePredictions_Btrajectory2Annotation_Jejunum.xlsx')
mapi <- read_xlsx('/home/Jayne.Wiarda/SI_PP_SC_ST/MappingPrediction/STomics_ClusteringAnnotation_CellTypePredictions_Btrajectory2Annotation_Ileum.xlsx')
mapClus <- as.data.frame(rbind(mapj, mapi))

mapj <- read_xlsx('/home/Jayne.Wiarda/SI_PP_SC_ST/MappingPrediction/STomics_ManualAnnotation_CellTypePredictions_Btrajectory2Annotation_Jejunum.xlsx')
mapi <- read_xlsx('/home/Jayne.Wiarda/SI_PP_SC_ST/MappingPrediction/STomics_ManualAnnotation_CellTypePredictions_Btrajectory2Annotation_Ileum.xlsx')
mapMan <- as.data.frame(rbind(mapj, mapi))

st <- readRDS('/home/Jayne.Wiarda/SI_PP_SC_ST/AnnotatedSeurat/AllSamples_annotated_ST.rds')
st$tissue <- substr(st$SampleID, 1, 1)   
Idents(st) <- st$Region_Clust
stClus <- subset(st, idents = c('Follicle', 'Crypt', 'Villus', 'IFZ/PFZ')) # remove muscularis spots from spatial dataset since they aren't represented in scRNA-seq data (muscularis was removed from scRNA-seq samples)
Idents(st) <- st$Region_ManAnn
stMan <- subset(st, idents = c('Follicle', 'Crypt', 'Villus', 'IFZ/PFZ')) # remove muscularis spots from spatial dataset since they aren't represented in scRNA-seq data (muscularis was removed from scRNA-seq samples)

order <- colnames(stClus)
mapClus <- mapClus[match(order, mapClus$Trajectory2_CellBarcode), ]
rownames(mapClus) <- mapClus$Trajectory2_CellBarcode
identical(colnames(stClus), rownames(mapClus))

order <- colnames(stMan)
mapMan <- mapMan[match(order, mapMan$Trajectory2_CellBarcode), ]
rownames(mapMan) <- mapMan$Trajectory2_CellBarcode
identical(colnames(stMan), rownames(mapMan))

stClus <- AddMetaData(stClus, mapClus)
stMan <- AddMetaData(stMan, mapMan)

Idents(stClus) <- stClus$tissue
stClusIl <- subset(stClus, idents = 'I')
stClusIl$Region_Clust <- factor(stClusIl$Region_Clust, levels = c('Villus', 'Crypt', 'IFZ/PFZ', 'Follicle'))
Idents(stClusIl) <- 'Region_Clust'

Idents(stClus) <- stClus$tissue
stClusJej <- subset(stClus, idents = 'J')
stClusJej$Region_Clust <- factor(stClusJej$Region_Clust, levels = c('Villus', 'Crypt', 'IFZ/PFZ', 'Follicle'))
Idents(stClusJej) <- 'Region_Clust'

Idents(stMan) <- stMan$tissue
stManIl <- subset(stMan, idents = 'I')
stManIl$Region_ManAnn <- factor(stManIl$Region_ManAnn, levels = c('Villus', 'Crypt', 'IFZ/PFZ', 'Follicle'))
Idents(stManIl) <- 'Region_ManAnn'

Idents(stMan) <- stMan$tissue
stManJej <- subset(stMan, idents = 'J')
stManJej$Region_ManAnn <- factor(stManJej$Region_ManAnn, levels = c('Villus', 'Crypt', 'IFZ/PFZ', 'Follicle'))
Idents(stManJej) <- 'Region_ManAnn'

sub1 <- as.data.frame(stClusJej@meta.data)
colnames(sub1)
sub1 <- sub1[,c(10, 16, 17, 24, 27, 34)]
colnames(sub1) <- c('Region', 'Increment2', 'Increment3', 'Increment4', 'Increment1', 'Increment5')
sub1 <- gather(sub1, key="Increment", value="predictionScore", 2:6)
sub1$Annotation <- rep('Clustering', nrow(sub1))
sub1$Region <- factor(sub1$Region, levels = c('Villus', 'Crypt', 'IFZ/PFZ', 'Follicle'))

sub2 <- as.data.frame(stManJej@meta.data)
colnames(sub2)
sub2 <- sub2[,c(10, 16, 17, 24, 27, 34)]
colnames(sub2) <- c('Region', 'Increment2', 'Increment3', 'Increment4', 'Increment1', 'Increment5')
sub2 <- gather(sub2, key="Increment", value="predictionScore", 2:6)
sub2$Annotation <- rep('Manual', nrow(sub2))
sub2$Region <- factor(sub2$Region, levels = c('Villus', 'Crypt', 'IFZ/PFZ', 'Follicle'))

sub <- rbind(sub1, sub2)

# Show both clustering and manual annotation side-by-side:
dat <- subset(sub, Increment == 'Increment1')
g1 <- ggplot(dat, aes(x=Region,
                      y=predictionScore,
                      fill = Annotation)) + 
  theme_bw() +
  geom_boxplot(aes(fill=Annotation), width=0.6, position = position_dodge(width=.8), 
               coef = 0, color = 'black', lwd = 0.6, outlier.shape = NA) +
  geom_point(shape=16, position = position_jitterdodge(jitter.width=.15, dodge.width = 0.8), size = .1, color = 'grey20') + 
  scale_fill_manual(breaks = dat$Annotation, values=c("grey85", "grey60")) + 
  stat_summary(fun.y="mean", geom="point", size=2.5,
               position=position_dodge(width=0.8), 
               color="red")+
  theme(axis.text.x=element_text(size=rel(1.75)), 
        axis.text.y=element_text(size=rel(1.75)))+ 
  labs(x="Region",y="Prediction score",
       title="Trajectory 2 - Increment 1 (Jejunum)") +
  theme(axis.title.x = element_text(size = 16),
        axis.text.x=element_text(size=12), 
        axis.title.y = element_text(size = 16),
        axis.text.y=element_text(size=12),
        legend.title = element_text (size=16),
        title = element_text(size = 16)) +
  theme(axis.text.x=element_text(size=rel(1.75)), axis.text.y=element_text(size=rel(1.75))) +
  scale_y_continuous(expand = c(0,0), lim = c(0,max(dat$predictionScore)))

dat <- subset(sub, Increment == 'Increment2')
g2 <- ggplot(dat, aes(x=Region,
                      y=predictionScore,
                      fill = Annotation)) + 
  theme_bw() +
  geom_boxplot(aes(fill=Annotation), width=0.6, position = position_dodge(width=.8), 
               coef = 0, color = 'black', lwd = 0.6, outlier.shape = NA) +
  geom_point(shape=16, position = position_jitterdodge(jitter.width=.15, dodge.width = 0.8), size = .1, color = 'grey20') + 
  scale_fill_manual(breaks = dat$Annotation, values=c("grey85", "grey60")) + 
  stat_summary(fun.y="mean", geom="point", size=2.5,
               position=position_dodge(width=0.8), 
               color="red")+
  theme(axis.text.x=element_text(size=rel(1.75)), 
        axis.text.y=element_text(size=rel(1.75)))+ 
  labs(x="Region",y="Prediction score",
       title="Trajectory 2 - Increment 2 (Jejunum)") +
  theme(axis.title.x = element_text(size = 16),
        axis.text.x=element_text(size=12), 
        axis.title.y = element_text(size = 16),
        axis.text.y=element_text(size=12),
        legend.title = element_text (size=16),
        title = element_text(size = 16)) +
  theme(axis.text.x=element_text(size=rel(1.75)), axis.text.y=element_text(size=rel(1.75))) +
  scale_y_continuous(expand = c(0,0), lim = c(0,max(dat$predictionScore)))

dat <- subset(sub, Increment == 'Increment3')
g3 <- ggplot(dat, aes(x=Region,
                      y=predictionScore,
                      fill = Annotation)) + 
  theme_bw() +
  geom_boxplot(aes(fill=Annotation), width=0.6, position = position_dodge(width=.8), 
               coef = 0, color = 'black', lwd = 0.6, outlier.shape = NA) +
  geom_point(shape=16, position = position_jitterdodge(jitter.width=.15, dodge.width = 0.8), size = .1, color = 'grey20') + 
  scale_fill_manual(breaks = dat$Annotation, values=c("grey85", "grey60")) + 
  stat_summary(fun.y="mean", geom="point", size=2.5,
               position=position_dodge(width=0.8), 
               color="red")+
  theme(axis.text.x=element_text(size=rel(1.75)), 
        axis.text.y=element_text(size=rel(1.75)))+ 
  labs(x="Region",y="Prediction score",
       title="Trajectory 2 - Increment 3 (Jejunum)") +
  theme(axis.title.x = element_text(size = 16),
        axis.text.x=element_text(size=12), 
        axis.title.y = element_text(size = 16),
        axis.text.y=element_text(size=12),
        legend.title = element_text (size=16),
        title = element_text(size = 16)) +
  theme(axis.text.x=element_text(size=rel(1.75)), axis.text.y=element_text(size=rel(1.75))) +
  scale_y_continuous(expand = c(0,0), lim = c(0,max(dat$predictionScore)))

dat <- subset(sub, Increment == 'Increment4')
g4 <- ggplot(dat, aes(x=Region,
                      y=predictionScore,
                      fill = Annotation)) + 
  theme_bw() +
  geom_boxplot(aes(fill=Annotation), width=0.6, position = position_dodge(width=.8), 
               coef = 0, color = 'black', lwd = 0.6, outlier.shape = NA) +
  geom_point(shape=16, position = position_jitterdodge(jitter.width=.15, dodge.width = 0.8), size = .1, color = 'grey20') + 
  scale_fill_manual(breaks = dat$Annotation, values=c("grey85", "grey60")) + 
  stat_summary(fun.y="mean", geom="point", size=2.5,
               position=position_dodge(width=0.8), 
               color="red")+
  theme(axis.text.x=element_text(size=rel(1.75)), 
        axis.text.y=element_text(size=rel(1.75)))+ 
  labs(x="Region",y="Prediction score",
       title="Trajectory 2 - Increment 4 (Jejunum)") +
  theme(axis.title.x = element_text(size = 16),
        axis.text.x=element_text(size=12), 
        axis.title.y = element_text(size = 16),
        axis.text.y=element_text(size=12),
        legend.title = element_text (size=16),
        title = element_text(size = 16)) +
  theme(axis.text.x=element_text(size=rel(1.75)), axis.text.y=element_text(size=rel(1.75))) +
  scale_y_continuous(expand = c(0,0), lim = c(0,max(dat$predictionScore)))

dat <- subset(sub, Increment == 'Increment5')
g5 <- ggplot(dat, aes(x=Region,
                      y=predictionScore,
                      fill = Annotation)) + 
  theme_bw() +
  geom_boxplot(aes(fill=Annotation), width=0.6, position = position_dodge(width=.8), 
               coef = 0, color = 'black', lwd = 0.6, outlier.shape = NA) +
  geom_point(shape=16, position = position_jitterdodge(jitter.width=.15, dodge.width = 0.8), size = .1, color = 'grey20') + 
  scale_fill_manual(breaks = dat$Annotation, values=c("grey85", "grey60")) + 
  stat_summary(fun.y="mean", geom="point", size=2.5,
               position=position_dodge(width=0.8), 
               color="red")+
  theme(axis.text.x=element_text(size=rel(1.75)), 
        axis.text.y=element_text(size=rel(1.75)))+ 
  labs(x="Region",y="Prediction score",
       title="Trajectory 2 - Increment 5 (Jejunum)") +
  theme(axis.title.x = element_text(size = 16),
        axis.text.x=element_text(size=12), 
        axis.title.y = element_text(size = 16),
        axis.text.y=element_text(size=12),
        legend.title = element_text (size=16),
        title = element_text(size = 16)) +
  theme(axis.text.x=element_text(size=rel(1.75)), axis.text.y=element_text(size=rel(1.75))) +
  scale_y_continuous(expand = c(0,0), lim = c(0,max(dat$predictionScore)))

g1+g2+g3+g4+g5

sub1 <- as.data.frame(stClusIl@meta.data)
colnames(sub1)
sub1 <- sub1[,c(10, 16, 17, 24, 27, 34)]
colnames(sub1) <- c('Region', 'Increment2', 'Increment3', 'Increment4', 'Increment1', 'Increment5')
sub1 <- gather(sub1, key="Increment", value="predictionScore", 2:6)
sub1$Annotation <- rep('Clustering', nrow(sub1))
sub1$Region <- factor(sub1$Region, levels = c('Villus', 'Crypt', 'IFZ/PFZ', 'Follicle'))

sub2 <- as.data.frame(stManIl@meta.data)
colnames(sub2)
sub2 <- sub2[,c(10, 16, 17, 24, 27, 34)]
colnames(sub2) <- c('Region', 'Increment2', 'Increment3', 'Increment4', 'Increment1', 'Increment5')
sub2 <- gather(sub2, key="Increment", value="predictionScore", 2:6)
sub2$Annotation <- rep('Manual', nrow(sub2))
sub2$Region <- factor(sub2$Region, levels = c('Villus', 'Crypt', 'IFZ/PFZ', 'Follicle'))

sub <- rbind(sub1, sub2)

dat <- subset(sub, Increment == 'Increment1')
g1 <- ggplot(dat, aes(x=Region,
                      y=predictionScore,
                      fill = Annotation)) + 
  theme_bw() +
  geom_boxplot(aes(fill=Annotation), width=0.6, position = position_dodge(width=.8), 
               coef = 0, color = 'black', lwd = 0.6, outlier.shape = NA) +
  geom_point(shape=16, position = position_jitterdodge(jitter.width=.15, dodge.width = 0.8), size = .1, color = 'grey20') + 
  scale_fill_manual(breaks = dat$Annotation, values=c("grey85", "grey60")) + 
  stat_summary(fun.y="mean", geom="point", size=2.5,
               position=position_dodge(width=0.8), 
               color="red")+
  theme(axis.text.x=element_text(size=rel(1.75)), 
        axis.text.y=element_text(size=rel(1.75)))+ 
  labs(x="Region",y="Prediction score",
       title="Trajectory 2 - Increment 1 (Ileum)") +
  theme(axis.title.x = element_text(size = 16),
        axis.text.x=element_text(size=12), 
        axis.title.y = element_text(size = 16),
        axis.text.y=element_text(size=12),
        legend.title = element_text (size=16),
        title = element_text(size = 16)) +
  theme(axis.text.x=element_text(size=rel(1.75)), axis.text.y=element_text(size=rel(1.75))) +
  scale_y_continuous(expand = c(0,0), lim = c(0,max(dat$predictionScore)))

dat <- subset(sub, Increment == 'Increment2')
g2 <- ggplot(dat, aes(x=Region,
                      y=predictionScore,
                      fill = Annotation)) + 
  theme_bw() +
  geom_boxplot(aes(fill=Annotation), width=0.6, position = position_dodge(width=.8), 
               coef = 0, color = 'black', lwd = 0.6, outlier.shape = NA) +
  geom_point(shape=16, position = position_jitterdodge(jitter.width=.15, dodge.width = 0.8), size = .1, color = 'grey20') + 
  scale_fill_manual(breaks = dat$Annotation, values=c("grey85", "grey60")) + 
  stat_summary(fun.y="mean", geom="point", size=2.5,
               position=position_dodge(width=0.8), 
               color="red")+
  theme(axis.text.x=element_text(size=rel(1.75)), 
        axis.text.y=element_text(size=rel(1.75)))+ 
  labs(x="Region",y="Prediction score",
       title="Trajectory 2 - Increment 2 (Ileum)") +
  theme(axis.title.x = element_text(size = 16),
        axis.text.x=element_text(size=12), 
        axis.title.y = element_text(size = 16),
        axis.text.y=element_text(size=12),
        legend.title = element_text (size=16),
        title = element_text(size = 16)) +
  theme(axis.text.x=element_text(size=rel(1.75)), axis.text.y=element_text(size=rel(1.75))) +
  scale_y_continuous(expand = c(0,0), lim = c(0,max(dat$predictionScore)))

dat <- subset(sub, Increment == 'Increment3')
g3 <- ggplot(dat, aes(x=Region,
                      y=predictionScore,
                      fill = Annotation)) + 
  theme_bw() +
  geom_boxplot(aes(fill=Annotation), width=0.6, position = position_dodge(width=.8), 
               coef = 0, color = 'black', lwd = 0.6, outlier.shape = NA) +
  geom_point(shape=16, position = position_jitterdodge(jitter.width=.15, dodge.width = 0.8), size = .1, color = 'grey20') + 
  scale_fill_manual(breaks = dat$Annotation, values=c("grey85", "grey60")) + 
  stat_summary(fun.y="mean", geom="point", size=2.5,
               position=position_dodge(width=0.8), 
               color="red")+
  theme(axis.text.x=element_text(size=rel(1.75)), 
        axis.text.y=element_text(size=rel(1.75)))+ 
  labs(x="Region",y="Prediction score",
       title="Trajectory 2 - Increment 3 (Ileum)") +
  theme(axis.title.x = element_text(size = 16),
        axis.text.x=element_text(size=12), 
        axis.title.y = element_text(size = 16),
        axis.text.y=element_text(size=12),
        legend.title = element_text (size=16),
        title = element_text(size = 16)) +
  theme(axis.text.x=element_text(size=rel(1.75)), axis.text.y=element_text(size=rel(1.75))) +
  scale_y_continuous(expand = c(0,0), lim = c(0,max(dat$predictionScore)))

dat <- subset(sub, Increment == 'Increment4')
g4 <- ggplot(dat, aes(x=Region,
                      y=predictionScore,
                      fill = Annotation)) + 
  theme_bw() +
  geom_boxplot(aes(fill=Annotation), width=0.6, position = position_dodge(width=.8), 
               coef = 0, color = 'black', lwd = 0.6, outlier.shape = NA) +
  geom_point(shape=16, position = position_jitterdodge(jitter.width=.15, dodge.width = 0.8), size = .1, color = 'grey20') + 
  scale_fill_manual(breaks = dat$Annotation, values=c("grey85", "grey60")) + 
  stat_summary(fun.y="mean", geom="point", size=2.5,
               position=position_dodge(width=0.8), 
               color="red")+
  theme(axis.text.x=element_text(size=rel(1.75)), 
        axis.text.y=element_text(size=rel(1.75)))+ 
  labs(x="Region",y="Prediction score",
       title="Trajectory 2 - Increment 4 (Ileum)") +
  theme(axis.title.x = element_text(size = 16),
        axis.text.x=element_text(size=12), 
        axis.title.y = element_text(size = 16),
        axis.text.y=element_text(size=12),
        legend.title = element_text (size=16),
        title = element_text(size = 16)) +
  theme(axis.text.x=element_text(size=rel(1.75)), axis.text.y=element_text(size=rel(1.75))) +
  scale_y_continuous(expand = c(0,0), lim = c(0,max(dat$predictionScore)))

dat <- subset(sub, Increment == 'Increment5')
g5 <- ggplot(dat, aes(x=Region,
                      y=predictionScore,
                      fill = Annotation)) + 
  theme_bw() +
  geom_boxplot(aes(fill=Annotation), width=0.6, position = position_dodge(width=.8), 
               coef = 0, color = 'black', lwd = 0.6, outlier.shape = NA) +
  geom_point(shape=16, position = position_jitterdodge(jitter.width=.15, dodge.width = 0.8), size = .1, color = 'grey20') + 
  scale_fill_manual(breaks = dat$Annotation, values=c("grey85", "grey60")) + 
  stat_summary(fun.y="mean", geom="point", size=2.5,
               position=position_dodge(width=0.8), 
               color="red")+
  theme(axis.text.x=element_text(size=rel(1.75)), 
        axis.text.y=element_text(size=rel(1.75)))+ 
  labs(x="Region",y="Prediction score",
       title="Trajectory 2 - Increment 5 (Ileum)") +
  theme(axis.title.x = element_text(size = 16),
        axis.text.x=element_text(size=12), 
        axis.title.y = element_text(size = 16),
        axis.text.y=element_text(size=12),
        legend.title = element_text (size=16),
        title = element_text(size = 16)) +
  theme(axis.text.x=element_text(size=rel(1.75)), axis.text.y=element_text(size=rel(1.75))) +
  scale_y_continuous(expand = c(0,0), lim = c(0,max(dat$predictionScore)))

g1+g2+g3+g4+g5

## Now let's just plot clustering results from trajectory 1 ----
mapj <- read_xlsx('/home/Jayne.Wiarda/SI_PP_SC_ST/MappingPrediction/STomics_ClusteringAnnotation_CellTypePredictions_Btrajectory1Annotation_Jejunum.xlsx')
mapi <- read_xlsx('/home/Jayne.Wiarda/SI_PP_SC_ST/MappingPrediction/STomics_ClusteringAnnotation_CellTypePredictions_Btrajectory1Annotation_Ileum.xlsx')
mapClus <- as.data.frame(rbind(mapj, mapi))

st <- readRDS('/home/Jayne.Wiarda/SI_PP_SC_ST/AnnotatedSeurat/AllSamples_annotated_ST.rds')
st$tissue <- substr(st$SampleID, 1, 1)   
Idents(st) <- st$Region_Clust
stClus <- subset(st, idents = c('Follicle', 'Crypt', 'Villus', 'IFZ/PFZ')) # remove muscularis spots from spatial dataset since they aren't represented in scRNA-seq data (muscularis was removed from scRNA-seq samples)

order <- colnames(stClus)
mapClus <- mapClus[match(order, mapClus$Trajectory1_CellBarcode), ]
rownames(mapClus) <- mapClus$Trajectory1_CellBarcode
identical(colnames(stClus), rownames(mapClus))

stClus <- AddMetaData(stClus, mapClus)

Idents(stClus) <- stClus$tissue
stClusIl <- subset(stClus, idents = 'I')
stClusIl$Region_Clust <- factor(stClusIl$Region_Clust, levels = c('Villus', 'Crypt', 'IFZ/PFZ', 'Follicle'))
Idents(stClusIl) <- 'Region_Clust'

Idents(stClus) <- stClus$tissue
stClusJej <- subset(stClus, idents = 'J')
stClusJej$Region_Clust <- factor(stClusJej$Region_Clust, levels = c('Villus', 'Crypt', 'IFZ/PFZ', 'Follicle'))
Idents(stClusJej) <- 'Region_Clust'

sub1 <- as.data.frame(stClusJej@meta.data)
colnames(sub1)
sub1 <- sub1[,c(10,16:19,37)]
colnames(sub1) <- c('Region', 'Increment2', 'Increment3', 'Increment5', 'Increment4', 'Increment1')
sub1 <- gather(sub1, key="Increment", value="predictionScore", 2:6)
sub1$Annotation <- rep('Clustering', nrow(sub1))
sub1$Region <- factor(sub1$Region, levels = c('Villus', 'Crypt', 'IFZ/PFZ', 'Follicle'))

sub2 <- as.data.frame(stClusIl@meta.data)
colnames(sub2)
sub2 <- sub2[,c(10,16:19,37)]
colnames(sub2) <- c('Region', 'Increment2', 'Increment3', 'Increment5', 'Increment4', 'Increment1')
sub2 <- gather(sub2, key="Increment", value="predictionScore", 2:6)
sub2$Annotation <- rep('Clustering', nrow(sub2))
sub2$Region <- factor(sub2$Region, levels = c('Villus', 'Crypt', 'IFZ/PFZ', 'Follicle'))

# Show just clustering annotation - Trajectory 1 jejunum ----
dat <- subset(sub1, Increment == 'Increment1')
g1 <- ggplot(dat, aes(x=Region,
                      y=predictionScore)) + 
  theme_bw() +
  geom_boxplot(fill = 'grey90', width=0.6, position = position_dodge(width=.8), 
               coef = 0, color = 'black', lwd = 0.6, outlier.shape = NA) +
  stat_summary(fun.y="mean", geom="point", size=3,   
               position=position_dodge(width=0.8),                                                 
               color="red") +
  theme(axis.text.x=element_text(size=rel(1.75)), 
        axis.text.y=element_text(size=rel(1.75)))+ 
  labs(x="Region",y="Prediction score",
       title="Trajectory 1 - Increment 1 (Jejunum)") +
  theme(axis.title.x = element_text(size = 16),
        axis.text.x=element_text(size=12), 
        axis.title.y = element_text(size = 16),
        axis.text.y=element_text(size=12),
        legend.title = element_text (size=16),
        title = element_text(size = 16)) +
  theme(axis.text.x=element_text(angle = 90, size=rel(1.75), hjust=0.9,vjust=0.2), axis.text.y=element_text(size=rel(1.75))) +
  scale_y_continuous(expand = c(0,0)) + coord_cartesian(ylim = c(0, 1))

dat <- subset(sub1, Increment == 'Increment2')
g2 <- ggplot(dat, aes(x=Region,
                      y=predictionScore)) + 
  theme_bw() +
  geom_boxplot(fill = 'grey90', width=0.6, position = position_dodge(width=.8), 
               coef = 0, color = 'black', lwd = 0.6, outlier.shape = NA) +
  stat_summary(fun.y="mean", geom="point", size=3,   
               position=position_dodge(width=0.8),                                                 
               color="red") +
  theme(axis.text.x=element_text(size=rel(1.75)), 
        axis.text.y=element_text(size=rel(1.75)))+ 
  labs(x="Region",y="Prediction score",
       title="Trajectory 1 - Increment 2 (Jejunum)") +
  theme(axis.title.x = element_text(size = 16),
        axis.text.x=element_text(size=12), 
        axis.title.y = element_text(size = 16),
        axis.text.y=element_text(size=12),
        legend.title = element_text (size=16),
        title = element_text(size = 16)) +
  theme(axis.text.x=element_text(angle = 90, size=rel(1.75), hjust=0.9,vjust=0.2), axis.text.y=element_text(size=rel(1.75))) +
  scale_y_continuous(expand = c(0,0)) + coord_cartesian(ylim = c(0, 1))

dat <- subset(sub1, Increment == 'Increment3')
g3 <- ggplot(dat, aes(x=Region,
                      y=predictionScore)) + 
  theme_bw() +
  geom_boxplot(fill = 'grey90', width=0.6, position = position_dodge(width=.8), 
               coef = 0, color = 'black', lwd = 0.6, outlier.shape = NA) +
  stat_summary(fun.y="mean", geom="point", size=3,   
               position=position_dodge(width=0.8),                                                 
               color="red") +
  theme(axis.text.x=element_text(size=rel(1.75)), 
        axis.text.y=element_text(size=rel(1.75)))+ 
  labs(x="Region",y="Prediction score",
       title="Trajectory 1 - Increment 3 (Jejunum)") +
  theme(axis.title.x = element_text(size = 16),
        axis.text.x=element_text(size=12), 
        axis.title.y = element_text(size = 16),
        axis.text.y=element_text(size=12),
        legend.title = element_text (size=16),
        title = element_text(size = 16)) +
  theme(axis.text.x=element_text(angle = 90, size=rel(1.75), hjust=0.9,vjust=0.2), axis.text.y=element_text(size=rel(1.75))) +
  scale_y_continuous(expand = c(0,0)) + coord_cartesian(ylim = c(0, 1))

dat <- subset(sub1, Increment == 'Increment4')
g4 <- ggplot(dat, aes(x=Region,
                      y=predictionScore)) + 
  theme_bw() +
  geom_boxplot(fill = 'grey90', width=0.6, position = position_dodge(width=.8), 
               coef = 0, color = 'black', lwd = 0.6, outlier.shape = NA) +
  stat_summary(fun.y="mean", geom="point", size=3,   
               position=position_dodge(width=0.8),                                                 
               color="red") +
  theme(axis.text.x=element_text(size=rel(1.75)), 
        axis.text.y=element_text(size=rel(1.75)))+ 
  labs(x="Region",y="Prediction score",
       title="Trajectory 1 - Increment 4 (Jejunum)") +
  theme(axis.title.x = element_text(size = 16),
        axis.text.x=element_text(size=12), 
        axis.title.y = element_text(size = 16),
        axis.text.y=element_text(size=12),
        legend.title = element_text (size=16),
        title = element_text(size = 16)) +
  theme(axis.text.x=element_text(angle = 90, size=rel(1.75), hjust=0.9,vjust=0.2), axis.text.y=element_text(size=rel(1.75))) +
  scale_y_continuous(expand = c(0,0)) + coord_cartesian(ylim = c(0, 1))

dat <- subset(sub1, Increment == 'Increment5')
g5 <- ggplot(dat, aes(x=Region,
                      y=predictionScore)) + 
  theme_bw() +
  geom_boxplot(fill = 'grey90', width=0.6, position = position_dodge(width=.8), 
               coef = 0, color = 'black', lwd = 0.6, outlier.shape = NA) +
  stat_summary(fun.y="mean", geom="point", size=3,   
               position=position_dodge(width=0.8),                                                 
               color="red") +
  theme(axis.text.x=element_text(size=rel(1.75)), 
        axis.text.y=element_text(size=rel(1.75)))+ 
  labs(x="Region",y="Prediction score",
       title="Trajectory 1 - Increment 5 (Jejunum)") +
  theme(axis.title.x = element_text(size = 16),
        axis.text.x=element_text(size=12), 
        axis.title.y = element_text(size = 16),
        axis.text.y=element_text(size=12),
        legend.title = element_text (size=16),
        title = element_text(size = 16)) +
  theme(axis.text.x=element_text(angle = 90, size=rel(1.75), hjust=0.9,vjust=0.2), axis.text.y=element_text(size=rel(1.75))) +
  scale_y_continuous(expand = c(0,0)) + coord_cartesian(ylim = c(0, 1))

g1+g2+g3+g4+g5

# Show just clustering annotation - Trajectory 1 Ileum ----
dat <- subset(sub2, Increment == 'Increment1')
g1 <- ggplot(dat, aes(x=Region,
                      y=predictionScore)) + 
  theme_bw() +
  geom_boxplot(fill = 'grey90', width=0.6, position = position_dodge(width=.8), 
               coef = 0, color = 'black', lwd = 0.6, outlier.shape = NA) +
  stat_summary(fun.y="mean", geom="point", size=3,   
               position=position_dodge(width=0.8),                                                 
               color="red") +
  theme(axis.text.x=element_text(size=rel(1.75)), 
        axis.text.y=element_text(size=rel(1.75)))+ 
  labs(x="Region",y="Prediction score",
       title="Trajectory 1 - Increment 1 (Ileum)") +
  theme(axis.title.x = element_text(size = 16),
        axis.text.x=element_text(size=12), 
        axis.title.y = element_text(size = 16),
        axis.text.y=element_text(size=12),
        legend.title = element_text (size=16),
        title = element_text(size = 16)) +
  theme(axis.text.x=element_text(angle = 90, size=rel(1.75), hjust=0.9,vjust=0.2), axis.text.y=element_text(size=rel(1.75))) +
  scale_y_continuous(expand = c(0,0)) + coord_cartesian(ylim = c(0, 1))

dat <- subset(sub2, Increment == 'Increment2')
g2 <- ggplot(dat, aes(x=Region,
                      y=predictionScore)) + 
  theme_bw() +
  geom_boxplot(fill = 'grey90', width=0.6, position = position_dodge(width=.8), 
               coef = 0, color = 'black', lwd = 0.6, outlier.shape = NA) +
  stat_summary(fun.y="mean", geom="point", size=3,   
               position=position_dodge(width=0.8),                                                 
               color="red") +
  theme(axis.text.x=element_text(size=rel(1.75)), 
        axis.text.y=element_text(size=rel(1.75)))+ 
  labs(x="Region",y="Prediction score",
       title="Trajectory 1 - Increment 2 (Ileum)") +
  theme(axis.title.x = element_text(size = 16),
        axis.text.x=element_text(size=12), 
        axis.title.y = element_text(size = 16),
        axis.text.y=element_text(size=12),
        legend.title = element_text (size=16),
        title = element_text(size = 16)) +
  theme(axis.text.x=element_text(angle = 90, size=rel(1.75), hjust=0.9,vjust=0.2), axis.text.y=element_text(size=rel(1.75))) +
  scale_y_continuous(expand = c(0,0)) + coord_cartesian(ylim = c(0, 1))

dat <- subset(sub2, Increment == 'Increment3')
g3 <- ggplot(dat, aes(x=Region,
                      y=predictionScore)) + 
  theme_bw() +
  geom_boxplot(fill = 'grey90', width=0.6, position = position_dodge(width=.8), 
               coef = 0, color = 'black', lwd = 0.6, outlier.shape = NA) +
  stat_summary(fun.y="mean", geom="point", size=3,   
               position=position_dodge(width=0.8),                                                 
               color="red") +
  theme(axis.text.x=element_text(size=rel(1.75)), 
        axis.text.y=element_text(size=rel(1.75)))+ 
  labs(x="Region",y="Prediction score",
       title="Trajectory 1 - Increment 3 (Ileum)") +
  theme(axis.title.x = element_text(size = 16),
        axis.text.x=element_text(size=12), 
        axis.title.y = element_text(size = 16),
        axis.text.y=element_text(size=12),
        legend.title = element_text (size=16),
        title = element_text(size = 16)) +
  theme(axis.text.x=element_text(angle = 90, size=rel(1.75), hjust=0.9,vjust=0.2), axis.text.y=element_text(size=rel(1.75))) +
  scale_y_continuous(expand = c(0,0)) + coord_cartesian(ylim = c(0, 1))

dat <- subset(sub2, Increment == 'Increment4')
g4 <- ggplot(dat, aes(x=Region,
                      y=predictionScore)) + 
  theme_bw() +
  geom_boxplot(fill = 'grey90', width=0.6, position = position_dodge(width=.8), 
               coef = 0, color = 'black', lwd = 0.6, outlier.shape = NA) +
  stat_summary(fun.y="mean", geom="point", size=3,   
               position=position_dodge(width=0.8),                                                 
               color="red") +
  theme(axis.text.x=element_text(size=rel(1.75)), 
        axis.text.y=element_text(size=rel(1.75)))+ 
  labs(x="Region",y="Prediction score",
       title="Trajectory 1 - Increment 4 (Ileum)") +
  theme(axis.title.x = element_text(size = 16),
        axis.text.x=element_text(size=12), 
        axis.title.y = element_text(size = 16),
        axis.text.y=element_text(size=12),
        legend.title = element_text (size=16),
        title = element_text(size = 16)) +
  theme(axis.text.x=element_text(angle = 90, size=rel(1.75), hjust=0.9,vjust=0.2), axis.text.y=element_text(size=rel(1.75))) +
  scale_y_continuous(expand = c(0,0)) + coord_cartesian(ylim = c(0, 1))

dat <- subset(sub2, Increment == 'Increment5')
g5 <- ggplot(dat, aes(x=Region,
                      y=predictionScore)) + 
  theme_bw() +
  geom_boxplot(fill = 'grey90', width=0.6, position = position_dodge(width=.8), 
               coef = 0, color = 'black', lwd = 0.6, outlier.shape = NA) +
  stat_summary(fun.y="mean", geom="point", size=3,   
               position=position_dodge(width=0.8),                                                 
               color="red") +
  theme(axis.text.x=element_text(size=rel(1.75)), 
        axis.text.y=element_text(size=rel(1.75)))+ 
  labs(x="Region",y="Prediction score",
       title="Trajectory 1 - Increment 5 (Ileum)") +
  theme(axis.title.x = element_text(size = 16),
        axis.text.x=element_text(size=12), 
        axis.title.y = element_text(size = 16),
        axis.text.y=element_text(size=12),
        legend.title = element_text (size=16),
        title = element_text(size = 16)) +
  theme(axis.text.x=element_text(angle = 90, size=rel(1.75), hjust=0.9,vjust=0.2), axis.text.y=element_text(size=rel(1.75))) +
  scale_y_continuous(expand = c(0,0)) + coord_cartesian(ylim = c(0, 1))

g1+g2+g3+g4+g5

## Now let's just plot clustering results from trajectory 2 ----
mapj <- read_xlsx('/home/Jayne.Wiarda/SI_PP_SC_ST/MappingPrediction/STomics_ClusteringAnnotation_CellTypePredictions_Btrajectory2Annotation_Jejunum.xlsx')
mapi <- read_xlsx('/home/Jayne.Wiarda/SI_PP_SC_ST/MappingPrediction/STomics_ClusteringAnnotation_CellTypePredictions_Btrajectory2Annotation_Ileum.xlsx')
mapClus <- as.data.frame(rbind(mapj, mapi))

st <- readRDS('/home/Jayne.Wiarda/SI_PP_SC_ST/AnnotatedSeurat/AllSamples_annotated_ST.rds')
st$tissue <- substr(st$SampleID, 1, 1)   
Idents(st) <- st$Region_Clust
stClus <- subset(st, idents = c('Follicle', 'Crypt', 'Villus', 'IFZ/PFZ')) # remove muscularis spots from spatial dataset since they aren't represented in scRNA-seq data (muscularis was removed from scRNA-seq samples)

order <- colnames(stClus)
mapClus <- mapClus[match(order, mapClus$Trajectory2_CellBarcode), ]
rownames(mapClus) <- mapClus$Trajectory2_CellBarcode
identical(colnames(stClus), rownames(mapClus))

stClus <- AddMetaData(stClus, mapClus)

Idents(stClus) <- stClus$tissue
stClusIl <- subset(stClus, idents = 'I')
stClusIl$Region_Clust <- factor(stClusIl$Region_Clust, levels = c('Villus', 'Crypt', 'IFZ/PFZ', 'Follicle'))
Idents(stClusIl) <- 'Region_Clust'

Idents(stClus) <- stClus$tissue
stClusJej <- subset(stClus, idents = 'J')
stClusJej$Region_Clust <- factor(stClusJej$Region_Clust, levels = c('Villus', 'Crypt', 'IFZ/PFZ', 'Follicle'))
Idents(stClusJej) <- 'Region_Clust'

sub1 <- as.data.frame(stClusJej@meta.data)
colnames(sub1)
sub1 <- sub1[,c(10,16:17, 24, 27, 34)]
colnames(sub1) <- c('Region', 'Increment2', 'Increment3', 'Increment4', 'Increment1', 'Increment5')
sub1 <- gather(sub1, key="Increment", value="predictionScore", 2:6)
sub1$Annotation <- rep('Clustering', nrow(sub1))
sub1$Region <- factor(sub1$Region, levels = c('Villus', 'Crypt', 'IFZ/PFZ', 'Follicle'))

sub2 <- as.data.frame(stClusIl@meta.data)
colnames(sub2)
sub2 <- sub2[,c(10,16:17, 24, 27, 34)]
colnames(sub2) <- c('Region', 'Increment2', 'Increment3', 'Increment4', 'Increment1', 'Increment5')
sub2 <- gather(sub2, key="Increment", value="predictionScore", 2:6)
sub2$Annotation <- rep('Clustering', nrow(sub2))
sub2$Region <- factor(sub2$Region, levels = c('Villus', 'Crypt', 'IFZ/PFZ', 'Follicle'))

# Show just clustering annotation - Trajectory 2 jejunum ----
dat <- subset(sub1, Increment == 'Increment1')
g1 <- ggplot(dat, aes(x=Region,
                      y=predictionScore)) + 
  theme_bw() +
  geom_boxplot(fill = 'grey90', width=0.6, position = position_dodge(width=.8), 
               coef = 0, color = 'black', lwd = 0.6, outlier.shape = NA) +
  stat_summary(fun.y="mean", geom="point", size=3,   
               position=position_dodge(width=0.8),                                                 
               color="red") +
  theme(axis.text.x=element_text(size=rel(1.75)), 
        axis.text.y=element_text(size=rel(1.75)))+ 
  labs(x="Region",y="Prediction score",
       title="Trajectory 2 - Increment 1 (Jejunum)") +
  theme(axis.title.x = element_text(size = 16),
        axis.text.x=element_text(size=12), 
        axis.title.y = element_text(size = 16),
        axis.text.y=element_text(size=12),
        legend.title = element_text (size=16),
        title = element_text(size = 16)) +
  theme(axis.text.x=element_text(angle = 90, size=rel(1.75), hjust=0.9,vjust=0.2), axis.text.y=element_text(size=rel(1.75))) +
  scale_y_continuous(expand = c(0,0)) + coord_cartesian(ylim = c(0, 0.25))

dat <- subset(sub1, Increment == 'Increment2')
g2 <- ggplot(dat, aes(x=Region,
                      y=predictionScore)) + 
  theme_bw() +
  geom_boxplot(fill = 'grey90', width=0.6, position = position_dodge(width=.8), 
               coef = 0, color = 'black', lwd = 0.6, outlier.shape = NA) +
  stat_summary(fun.y="mean", geom="point", size=3,   
               position=position_dodge(width=0.8),                                                 
               color="red") +
  theme(axis.text.x=element_text(size=rel(1.75)), 
        axis.text.y=element_text(size=rel(1.75)))+ 
  labs(x="Region",y="Prediction score",
       title="Trajectory 2 - Increment 2 (Jejunum)") +
  theme(axis.title.x = element_text(size = 16),
        axis.text.x=element_text(size=12), 
        axis.title.y = element_text(size = 16),
        axis.text.y=element_text(size=12),
        legend.title = element_text (size=16),
        title = element_text(size = 16)) +
  theme(axis.text.x=element_text(angle = 90, size=rel(1.75), hjust=0.9,vjust=0.2), axis.text.y=element_text(size=rel(1.75))) +
  scale_y_continuous(expand = c(0,0)) + coord_cartesian(ylim = c(0, 0.25))

dat <- subset(sub1, Increment == 'Increment3')
g3 <- ggplot(dat, aes(x=Region,
                      y=predictionScore)) + 
  theme_bw() +
  geom_boxplot(fill = 'grey90', width=0.6, position = position_dodge(width=.8), 
               coef = 0, color = 'black', lwd = 0.6, outlier.shape = NA) +
  stat_summary(fun.y="mean", geom="point", size=3,   
               position=position_dodge(width=0.8),                                                 
               color="red") +
  theme(axis.text.x=element_text(size=rel(1.75)), 
        axis.text.y=element_text(size=rel(1.75)))+ 
  labs(x="Region",y="Prediction score",
       title="Trajectory 2 - Increment 3 (Jejunum)") +
  theme(axis.title.x = element_text(size = 16),
        axis.text.x=element_text(size=12), 
        axis.title.y = element_text(size = 16),
        axis.text.y=element_text(size=12),
        legend.title = element_text (size=16),
        title = element_text(size = 16)) +
  theme(axis.text.x=element_text(angle = 90, size=rel(1.75), hjust=0.9,vjust=0.2), axis.text.y=element_text(size=rel(1.75))) +
  scale_y_continuous(expand = c(0,0)) + coord_cartesian(ylim = c(0, 0.25))

dat <- subset(sub1, Increment == 'Increment4')
g4 <- ggplot(dat, aes(x=Region,
                      y=predictionScore)) + 
  theme_bw() +
  geom_boxplot(fill = 'grey90', width=0.6, position = position_dodge(width=.8), 
               coef = 0, color = 'black', lwd = 0.6, outlier.shape = NA) +
  stat_summary(fun.y="mean", geom="point", size=3,   
               position=position_dodge(width=0.8),                                                 
               color="red") +
  theme(axis.text.x=element_text(size=rel(1.75)), 
        axis.text.y=element_text(size=rel(1.75)))+ 
  labs(x="Region",y="Prediction score",
       title="Trajectory 2 - Increment 4 (Jejunum)") +
  theme(axis.title.x = element_text(size = 16),
        axis.text.x=element_text(size=12), 
        axis.title.y = element_text(size = 16),
        axis.text.y=element_text(size=12),
        legend.title = element_text (size=16),
        title = element_text(size = 16)) +
  theme(axis.text.x=element_text(angle = 90, size=rel(1.75), hjust=0.9,vjust=0.2), axis.text.y=element_text(size=rel(1.75))) +
  scale_y_continuous(expand = c(0,0)) + coord_cartesian(ylim = c(0, 0.25))

dat <- subset(sub1, Increment == 'Increment5')
g5 <- ggplot(dat, aes(x=Region,
                      y=predictionScore)) + 
  theme_bw() +
  geom_boxplot(fill = 'grey90', width=0.6, position = position_dodge(width=.8), 
               coef = 0, color = 'black', lwd = 0.6, outlier.shape = NA) +
  stat_summary(fun.y="mean", geom="point", size=3,   
               position=position_dodge(width=0.8),                                                 
               color="red") +
  theme(axis.text.x=element_text(size=rel(1.75)), 
        axis.text.y=element_text(size=rel(1.75)))+ 
  labs(x="Region",y="Prediction score",
       title="Trajectory 2 - Increment 5 (Jejunum)") +
  theme(axis.title.x = element_text(size = 16),
        axis.text.x=element_text(size=12), 
        axis.title.y = element_text(size = 16),
        axis.text.y=element_text(size=12),
        legend.title = element_text (size=16),
        title = element_text(size = 16)) +
  theme(axis.text.x=element_text(angle = 90, size=rel(1.75), hjust=0.9,vjust=0.2), axis.text.y=element_text(size=rel(1.75))) +
  scale_y_continuous(expand = c(0,0)) + coord_cartesian(ylim = c(0, 0.25))

g1+g2+g3+g4+g5

# Show just clustering annotation - Trajectory 2 Ileum ----
dat <- subset(sub2, Increment == 'Increment1')
g1 <- ggplot(dat, aes(x=Region,
                      y=predictionScore)) + 
  theme_bw() +
  geom_boxplot(fill = 'grey90', width=0.6, position = position_dodge(width=.8), 
               coef = 0, color = 'black', lwd = 0.6, outlier.shape = NA) +
  stat_summary(fun.y="mean", geom="point", size=3,   
               position=position_dodge(width=0.8),                                                 
               color="red") +
  theme(axis.text.x=element_text(size=rel(1.75)), 
        axis.text.y=element_text(size=rel(1.75)))+ 
  labs(x="Region",y="Prediction score",
       title="Trajectory 2 - Increment 1 (Ileum)") +
  theme(axis.title.x = element_text(size = 16),
        axis.text.x=element_text(size=12), 
        axis.title.y = element_text(size = 16),
        axis.text.y=element_text(size=12),
        legend.title = element_text (size=16),
        title = element_text(size = 16)) +
  theme(axis.text.x=element_text(angle = 90, size=rel(1.75), hjust=0.9,vjust=0.2), axis.text.y=element_text(size=rel(1.75))) +
  scale_y_continuous(expand = c(0,0)) + coord_cartesian(ylim = c(0, 0.25))

dat <- subset(sub2, Increment == 'Increment2')
g2 <- ggplot(dat, aes(x=Region,
                      y=predictionScore)) + 
  theme_bw() +
  geom_boxplot(fill = 'grey90', width=0.6, position = position_dodge(width=.8), 
               coef = 0, color = 'black', lwd = 0.6, outlier.shape = NA) +
  stat_summary(fun.y="mean", geom="point", size=3,   
               position=position_dodge(width=0.8),                                                 
               color="red") +
  theme(axis.text.x=element_text(size=rel(1.75)), 
        axis.text.y=element_text(size=rel(1.75)))+ 
  labs(x="Region",y="Prediction score",
       title="Trajectory 2 - Increment 2 (Ileum)") +
  theme(axis.title.x = element_text(size = 16),
        axis.text.x=element_text(size=12), 
        axis.title.y = element_text(size = 16),
        axis.text.y=element_text(size=12),
        legend.title = element_text (size=16),
        title = element_text(size = 16)) +
  theme(axis.text.x=element_text(angle = 90, size=rel(1.75), hjust=0.9,vjust=0.2), axis.text.y=element_text(size=rel(1.75))) +
  scale_y_continuous(expand = c(0,0)) + coord_cartesian(ylim = c(0, 0.25))

dat <- subset(sub2, Increment == 'Increment3')
g3 <- ggplot(dat, aes(x=Region,
                      y=predictionScore)) + 
  theme_bw() +
  geom_boxplot(fill = 'grey90', width=0.6, position = position_dodge(width=.8), 
               coef = 0, color = 'black', lwd = 0.6, outlier.shape = NA) +
  stat_summary(fun.y="mean", geom="point", size=3,   
               position=position_dodge(width=0.8),                                                 
               color="red") +
  theme(axis.text.x=element_text(size=rel(1.75)), 
        axis.text.y=element_text(size=rel(1.75)))+ 
  labs(x="Region",y="Prediction score",
       title="Trajectory 2 - Increment 3 (Ileum)") +
  theme(axis.title.x = element_text(size = 16),
        axis.text.x=element_text(size=12), 
        axis.title.y = element_text(size = 16),
        axis.text.y=element_text(size=12),
        legend.title = element_text (size=16),
        title = element_text(size = 16)) +
  theme(axis.text.x=element_text(angle = 90, size=rel(1.75), hjust=0.9,vjust=0.2), axis.text.y=element_text(size=rel(1.75))) +
  scale_y_continuous(expand = c(0,0)) + coord_cartesian(ylim = c(0, 0.25))

dat <- subset(sub2, Increment == 'Increment4')
g4 <- ggplot(dat, aes(x=Region,
                      y=predictionScore)) + 
  theme_bw() +
  geom_boxplot(fill = 'grey90', width=0.6, position = position_dodge(width=.8), 
               coef = 0, color = 'black', lwd = 0.6, outlier.shape = NA) +
  stat_summary(fun.y="mean", geom="point", size=3,   
               position=position_dodge(width=0.8),                                                 
               color="red") +
  theme(axis.text.x=element_text(size=rel(1.75)), 
        axis.text.y=element_text(size=rel(1.75)))+ 
  labs(x="Region",y="Prediction score",
       title="Trajectory 2 - Increment 4 (Ileum)") +
  theme(axis.title.x = element_text(size = 16),
        axis.text.x=element_text(size=12), 
        axis.title.y = element_text(size = 16),
        axis.text.y=element_text(size=12),
        legend.title = element_text (size=16),
        title = element_text(size = 16)) +
  theme(axis.text.x=element_text(angle = 90, size=rel(1.75), hjust=0.9,vjust=0.2), axis.text.y=element_text(size=rel(1.75))) +
  scale_y_continuous(expand = c(0,0)) + coord_cartesian(ylim = c(0, 0.25))

dat <- subset(sub2, Increment == 'Increment5')
g5 <- ggplot(dat, aes(x=Region,
                      y=predictionScore)) + 
  theme_bw() +
  geom_boxplot(fill = 'grey90', width=0.6, position = position_dodge(width=.8), 
               coef = 0, color = 'black', lwd = 0.6, outlier.shape = NA) +
  stat_summary(fun.y="mean", geom="point", size=3,   
               position=position_dodge(width=0.8),                                                 
               color="red") +
  theme(axis.text.x=element_text(size=rel(1.75)), 
        axis.text.y=element_text(size=rel(1.75)))+ 
  labs(x="Region",y="Prediction score",
       title="Trajectory 2 - Increment 5 (Ileum)") +
  theme(axis.title.x = element_text(size = 16),
        axis.text.x=element_text(size=12), 
        axis.title.y = element_text(size = 16),
        axis.text.y=element_text(size=12),
        legend.title = element_text (size=16),
        title = element_text(size = 16)) +
  theme(axis.text.x=element_text(angle = 90, size=rel(1.75), hjust=0.9,vjust=0.2), axis.text.y=element_text(size=rel(1.75))) +
  scale_y_continuous(expand = c(0,0)) + coord_cartesian(ylim = c(0, 0.25))

g1+g2+g3+g4+g5

# Make bar plots for prediction scores only in follicle regions ----
## Comparing annotation methods ----
# Split barplot of mapping location predictions by cell type in follicles
mapj <- read_xlsx('/home/Jayne.Wiarda/SI_PP_SC_ST/MappingPrediction/STomics_CellTypePredictions_ClusteringAnnotation_Jejunum.xlsx')
mapi <- read_xlsx('/home/Jayne.Wiarda/SI_PP_SC_ST/MappingPrediction/STomics_CellTypePredictions_ClusteringAnnotation_Ileum.xlsx')
mapClus <- as.data.frame(rbind(mapj, mapi))
rownames(mapClus) <- mapClus$CellBarcode
mapClus <- subset(mapClus, select=-c(CellBarcode))

mapj <- read_xlsx('/home/Jayne.Wiarda/SI_PP_SC_ST/MappingPrediction/STomics_CellTypePredictions_ManualAnnotation_Jejunum.xlsx')
mapi <- read_xlsx('/home/Jayne.Wiarda/SI_PP_SC_ST/MappingPrediction/STomics_CellTypePredictions_ManualAnnotation_Ileum.xlsx')
mapMan <- as.data.frame(rbind(mapj, mapi))
rownames(mapMan) <- mapMan$CellBarcode
mapMan <- subset(mapMan, select=-c(CellBarcode))

st <- readRDS('/home/Jayne.Wiarda/SI_PP_SC_ST/AnnotatedSeurat/AllSamples_annotated_ST.rds')
st$tissue <- substr(st$SampleID, 1, 1)   
Idents(st) <- st$Region_Clust
stClus <- subset(st, idents = c('Follicle', 'Crypt', 'Villus', 'IFZ/PFZ')) # remove muscularis spots from spatial dataset since they aren't represented in scRNA-seq data (muscularis was removed from scRNA-seq samples)
Idents(st) <- st$Region_ManAnn
stMan <- subset(st, idents = c('Follicle', 'Crypt', 'Villus', 'IFZ/PFZ')) # remove muscularis spots from spatial dataset since they aren't represented in scRNA-seq data (muscularis was removed from scRNA-seq samples)

stClus <- AddMetaData(stClus, mapClus)
stMan <- AddMetaData(stMan, mapMan)

Idents(stClus) <- stClus$tissue
stClusIl <- subset(stClus, idents = 'I')
stClusIl$Region_Clust <- factor(stClusIl$Region_Clust, levels = c('Villus', 'Crypt', 'IFZ/PFZ', 'Follicle'))

Idents(stMan) <- stMan$tissue
stManIl <- subset(stMan, idents = 'I')
stManIl$Region_ManAnn <- factor(stManIl$Region_ManAnn, levels = c('Villus', 'Crypt', 'IFZ/PFZ', 'Follicle'))

Idents(stClusIl) <- 'Region_Clust'
sub1 <- subset(stClusIl, idents = 'Follicle')
sub1 <- as.data.frame(sub1@meta.data)
colnames(sub1)
sub1 <- sub1[,c(15:46 )]
colnames(sub1) <- c('Non.naive.CD8.ab.T.cells', 'Resting.B.cells', 
                    "Cycling.B.cells", "Activated.B.cells", "Non.naive.CD4.ab.T.cells",
                    "Naive.CD4.CD8.ab.T.cells", "Follicular.CD4.ab.T.cells", "CD2neg.GD.T.cells",
                    "Non.naive.group.1.ILCs",  "Antibody.secreting.cells", 
                    "Macrophages", "Group.3.ILCs","Transitioning.B.cells", "Cytotoxic.gd.T.cells",     
                    "Enterocytes", "Cytotoxic.CD8.ab.T.cells",  "Crypt.cells", "Cytotoxic.group.1.ILCs",   
                    "Non.naive.gd.T.cells", "Fibroblasts", "Dendritic.cells","Cycling.CD4.ab.T.cells",   
                    "SELLhi.gd.T.cells", "Cycling.CD8.ab.T.cells",  "Mast.cells", "NEUROD1hi.EE.cells",       
                    "Endothelial.cells", "Cycling.gd.T.cells",  "Goblet.cells",  "Muscle.cells",             
                    "NEUROD1lo.EE.cells", "BEST4.enterocytes")
sub1 <- gather(sub1, key="celltype", value="predictionScore", 1:32)
sub1$Annotation <- rep('Clustering', nrow(sub1))
sub1$celltype <- factor(sub1$celltype, levels = c("Activated.B.cells", "Cycling.B.cells", "Resting.B.cells",
                                                  "Transitioning.B.cells", "Antibody.secreting.cells",
                                                  "Cycling.CD4.ab.T.cells", "Cycling.CD8.ab.T.cells",
                                                  "Cycling.gd.T.cells", "Cytotoxic.CD8.ab.T.cells",
                                                  "Cytotoxic.gd.T.cells", "Cytotoxic.group.1.ILCs",
                                                  "Non.naive.CD8.ab.T.cells",
                                                  "Non.naive.gd.T.cells", "Non.naive.group.1.ILCs",
                                                  "SELLhi.gd.T.cells", "CD2neg.GD.T.cells", 
                                                  "Naive.CD4.CD8.ab.T.cells", "Non.naive.CD4.ab.T.cells",
                                                  "Follicular.CD4.ab.T.cells", "Group.3.ILCs",
                                                  "Dendritic.cells", "Macrophages", "Mast.cells",
                                                  "Crypt.cells", "Enterocytes", "BEST4.enterocytes",
                                                  "Goblet.cells",
                                                  "NEUROD1lo.EE.cells", "NEUROD1hi.EE.cells",
                                                  "Endothelial.cells", "Fibroblasts", "Muscle.cells"))

Idents(stManIl) <- 'Region_ManAnn'
sub2 <- subset(stManIl, idents = 'Follicle')
sub2 <- as.data.frame(sub2@meta.data)
colnames(sub2)
sub2 <- sub2[,c(15:46 )]
colnames(sub2) <- c('Non.naive.CD8.ab.T.cells', 'Resting.B.cells', 
                    "Cycling.B.cells", "Activated.B.cells", "Non.naive.CD4.ab.T.cells",
                    "Naive.CD4.CD8.ab.T.cells", "Follicular.CD4.ab.T.cells", "CD2neg.GD.T.cells",
                    "Non.naive.group.1.ILCs",  "Antibody.secreting.cells", 
                    "Macrophages", "Group.3.ILCs","Transitioning.B.cells", "Cytotoxic.gd.T.cells",     
                    "Enterocytes", "Cytotoxic.CD8.ab.T.cells",  "Crypt.cells", "Cytotoxic.group.1.ILCs",   
                    "Non.naive.gd.T.cells", "Fibroblasts", "Dendritic.cells","Cycling.CD4.ab.T.cells",   
                    "SELLhi.gd.T.cells", "Cycling.CD8.ab.T.cells",  "Mast.cells", "NEUROD1hi.EE.cells",       
                    "Endothelial.cells", "Cycling.gd.T.cells",  "Goblet.cells",  "Muscle.cells",             
                    "NEUROD1lo.EE.cells", "BEST4.enterocytes")
sub2 <- gather(sub2, key="celltype", value="predictionScore", 1:32)
sub2$Annotation <- rep('Manual', nrow(sub2))
sub2$celltype <- factor(sub2$celltype, levels = c("Activated.B.cells", "Cycling.B.cells", "Resting.B.cells",
                                                  "Transitioning.B.cells", "Antibody.secreting.cells",
                                                  "Cycling.CD4.ab.T.cells", "Cycling.CD8.ab.T.cells",
                                                  "Cycling.gd.T.cells", "Cytotoxic.CD8.ab.T.cells",
                                                  "Cytotoxic.gd.T.cells", "Cytotoxic.group.1.ILCs",
                                                  "Non.naive.CD8.ab.T.cells",
                                                  "Non.naive.gd.T.cells", "Non.naive.group.1.ILCs",
                                                  "SELLhi.gd.T.cells", "CD2neg.GD.T.cells", 
                                                  "Naive.CD4.CD8.ab.T.cells", "Non.naive.CD4.ab.T.cells",
                                                  "Follicular.CD4.ab.T.cells", "Group.3.ILCs",
                                                  "Dendritic.cells", "Macrophages", "Mast.cells",
                                                  "Crypt.cells", "Enterocytes", "BEST4.enterocytes",
                                                  "Goblet.cells",
                                                  "NEUROD1lo.EE.cells", "NEUROD1hi.EE.cells",
                                                  "Endothelial.cells", "Fibroblasts", "Muscle.cells"))

dat <- rbind(sub1, sub2)
dat$Annotation <- factor(dat$Annotation, levels = c('Manual', 'Clustering'))
ggplot(dat, aes(x=celltype,
                y=predictionScore,
                fill = Annotation)) + 
  theme_bw() +
  geom_boxplot(aes(fill=Annotation), width=0.6, position = position_dodge(width=.8), 
               coef = 0, color = 'black', lwd = 0.6, outlier.shape = NA) +
  geom_point(shape=16, position = position_jitterdodge(jitter.width=.3, dodge.width = 0.8), size = .1, color = 'grey30') + 
  scale_fill_manual(values=c('grey90', 'grey65')) + 
  stat_summary(fun.y="mean", geom="point", size=2.5,
               position=position_dodge(width=0.8), 
               color="red")+
  theme(axis.text.x=element_text(size=rel(1.75)), 
        axis.text.y=element_text(size=rel(1.75)))+ 
  labs(x="Cell type",y="Prediction score",
       title="Prediction scores to follicular regions (Ileum)") +
  theme(axis.title.x = element_text(size = 16),
        axis.text.x=element_text(size=12, angle = 45, hjust=1), 
        axis.title.y = element_text(size = 16),
        axis.text.y=element_text(size=12),
        legend.title = element_text (size=16),
        title = element_text(size = 16)) +
  theme(axis.text.x=element_text(size=rel(1.75)), axis.text.y=element_text(size=rel(1.75))) +
  scale_y_continuous(expand = c(0,0), lim = c(0,max(dat$predictionScore)))

## Only for clustering annotation: 
dat <- subset(dat, Annotation == 'Clustering')
ggplot(dat, aes(x=celltype,
                y=predictionScore)) + 
  theme_bw() +
  geom_boxplot(fill = 'grey90', width=0.6, position = position_dodge(width=.8), 
               coef = 0, color = 'black', lwd = 0.6, outlier.shape = NA) +
  stat_summary(fun.y="mean", geom="point", size=3,   
               position=position_dodge(width=0.8),                                                 
               color="red") +
  theme(axis.text.x=element_text(size=rel(1.75)), 
        axis.text.y=element_text(size=rel(1.75)))+ 
  labs(x="Region",y="Prediction score",
       title="Cell type prediction scores in follicles (Ileum)") +
  theme(axis.title.x = element_text(size = 16),
        axis.text.x=element_text(size=12), 
        axis.title.y = element_text(size = 16),
        axis.text.y=element_text(size=12),
        legend.title = element_text (size=16),
        title = element_text(size = 16)) +
  theme(axis.text.x=element_text(angle = 90, size=rel(1.75), hjust=0.9,vjust=0.2), axis.text.y=element_text(size=rel(1.75))) +
  scale_y_continuous(expand = c(0,0)) + coord_cartesian(ylim = c(0, 1))

Idents(stClus) <- stClus$tissue
stClusJej <- subset(stClus, idents = 'J')
stClusJej$Region_Clust <- factor(stClusJej$Region_Clust, levels = c('Villus', 'Crypt', 'IFZ/PFZ', 'Follicle'))

Idents(stMan) <- stMan$tissue
stManJej <- subset(stMan, idents = 'J')
stManJej$Region_ManAnn <- factor(stManJej$Region_ManAnn, levels = c('Villus', 'Crypt', 'IFZ/PFZ', 'Follicle'))

Idents(stClusJej) <- 'Region_Clust'
sub1 <- subset(stClusJej, idents = 'Follicle')
sub1 <- as.data.frame(sub1@meta.data)
colnames(sub1)
sub1 <- sub1[,c(15:46 )]
colnames(sub1) <- c('Non.naive.CD8.ab.T.cells', 'Resting.B.cells', 
                    "Cycling.B.cells", "Activated.B.cells", "Non.naive.CD4.ab.T.cells",
                    "Naive.CD4.CD8.ab.T.cells", "Follicular.CD4.ab.T.cells", "CD2neg.GD.T.cells",
                    "Non.naive.group.1.ILCs",  "Antibody.secreting.cells", 
                    "Macrophages", "Group.3.ILCs","Transitioning.B.cells", "Cytotoxic.gd.T.cells",     
                    "Enterocytes", "Cytotoxic.CD8.ab.T.cells",  "Crypt.cells", "Cytotoxic.group.1.ILCs",   
                    "Non.naive.gd.T.cells", "Fibroblasts", "Dendritic.cells","Cycling.CD4.ab.T.cells",   
                    "SELLhi.gd.T.cells", "Cycling.CD8.ab.T.cells",  "Mast.cells", "NEUROD1hi.EE.cells",       
                    "Endothelial.cells", "Cycling.gd.T.cells",  "Goblet.cells",  "Muscle.cells",             
                    "NEUROD1lo.EE.cells", "BEST4.enterocytes")
sub1 <- gather(sub1, key="celltype", value="predictionScore", 1:32)
sub1$Annotation <- rep('Clustering', nrow(sub1))
sub1$celltype <- factor(sub1$celltype, levels = c("Activated.B.cells", "Cycling.B.cells", "Resting.B.cells",
                                                  "Transitioning.B.cells", "Antibody.secreting.cells",
                                                  "Cycling.CD4.ab.T.cells", "Cycling.CD8.ab.T.cells",
                                                  "Cycling.gd.T.cells", "Cytotoxic.CD8.ab.T.cells",
                                                  "Cytotoxic.gd.T.cells", "Cytotoxic.group.1.ILCs",
                                                  "Non.naive.CD8.ab.T.cells",
                                                  "Non.naive.gd.T.cells", "Non.naive.group.1.ILCs",
                                                  "SELLhi.gd.T.cells", "CD2neg.GD.T.cells", 
                                                  "Naive.CD4.CD8.ab.T.cells", "Non.naive.CD4.ab.T.cells",
                                                  "Follicular.CD4.ab.T.cells", "Group.3.ILCs",
                                                  "Dendritic.cells", "Macrophages", "Mast.cells",
                                                  "Crypt.cells", "Enterocytes", "BEST4.enterocytes",
                                                  "Goblet.cells",
                                                  "NEUROD1lo.EE.cells", "NEUROD1hi.EE.cells",
                                                  "Endothelial.cells", "Fibroblasts", "Muscle.cells"))

Idents(stManJej) <- 'Region_ManAnn'
sub2 <- subset(stManJej, idents = 'Follicle')
sub2 <- as.data.frame(sub2@meta.data)
colnames(sub2)
sub2 <- sub2[,c(15:46 )]
colnames(sub2) <- c('Non.naive.CD8.ab.T.cells', 'Resting.B.cells', 
                    "Cycling.B.cells", "Activated.B.cells", "Non.naive.CD4.ab.T.cells",
                    "Naive.CD4.CD8.ab.T.cells", "Follicular.CD4.ab.T.cells", "CD2neg.GD.T.cells",
                    "Non.naive.group.1.ILCs",  "Antibody.secreting.cells", 
                    "Macrophages", "Group.3.ILCs","Transitioning.B.cells", "Cytotoxic.gd.T.cells",     
                    "Enterocytes", "Cytotoxic.CD8.ab.T.cells",  "Crypt.cells", "Cytotoxic.group.1.ILCs",   
                    "Non.naive.gd.T.cells", "Fibroblasts", "Dendritic.cells","Cycling.CD4.ab.T.cells",   
                    "SELLhi.gd.T.cells", "Cycling.CD8.ab.T.cells",  "Mast.cells", "NEUROD1hi.EE.cells",       
                    "Endothelial.cells", "Cycling.gd.T.cells",  "Goblet.cells",  "Muscle.cells",             
                    "NEUROD1lo.EE.cells", "BEST4.enterocytes")
sub2 <- gather(sub2, key="celltype", value="predictionScore", 1:32)
sub2$Annotation <- rep('Manual', nrow(sub2))
sub2$celltype <- factor(sub2$celltype, levels = c("Activated.B.cells", "Cycling.B.cells", "Resting.B.cells",
                                                  "Transitioning.B.cells", "Antibody.secreting.cells",
                                                  "Cycling.CD4.ab.T.cells", "Cycling.CD8.ab.T.cells",
                                                  "Cycling.gd.T.cells", "Cytotoxic.CD8.ab.T.cells",
                                                  "Cytotoxic.gd.T.cells", "Cytotoxic.group.1.ILCs",
                                                  "Non.naive.CD8.ab.T.cells",
                                                  "Non.naive.gd.T.cells", "Non.naive.group.1.ILCs",
                                                  "SELLhi.gd.T.cells", "CD2neg.GD.T.cells", 
                                                  "Naive.CD4.CD8.ab.T.cells", "Non.naive.CD4.ab.T.cells",
                                                  "Follicular.CD4.ab.T.cells", "Group.3.ILCs",
                                                  "Dendritic.cells", "Macrophages", "Mast.cells",
                                                  "Crypt.cells", "Enterocytes", "BEST4.enterocytes",
                                                  "Goblet.cells",
                                                  "NEUROD1lo.EE.cells", "NEUROD1hi.EE.cells",
                                                  "Endothelial.cells", "Fibroblasts", "Muscle.cells"))

dat <- rbind(sub1, sub2)
dat$Annotation <- factor(dat$Annotation, levels = c('Manual', 'Clustering'))
ggplot(dat, aes(x=celltype,
                y=predictionScore,
                fill = Annotation)) + 
  theme_bw() +
  geom_boxplot(aes(fill=Annotation), width=0.6, position = position_dodge(width=.8), 
               coef = 0, color = 'black', lwd = 0.6, outlier.shape = NA) +
  geom_point(shape=16, position = position_jitterdodge(jitter.width=.3, dodge.width = 0.8), size = .1, color = 'grey30') + 
  scale_fill_manual(values=c('grey90', 'grey65')) + 
  stat_summary(fun.y="mean", geom="point", size=2.5,
               position=position_dodge(width=0.8), 
               color="red")+
  theme(axis.text.x=element_text(size=rel(1.75)), 
        axis.text.y=element_text(size=rel(1.75)))+ 
  labs(x="Cell type",y="Prediction score",
       title="Prediction scores to follicular regions (Jejunum)") +
  theme(axis.title.x = element_text(size = 16),
        axis.text.x=element_text(size=12, angle = 45, hjust=1), 
        axis.title.y = element_text(size = 16),
        axis.text.y=element_text(size=12),
        legend.title = element_text (size=16),
        title = element_text(size = 16)) +
  theme(axis.text.x=element_text(size=rel(1.75)), axis.text.y=element_text(size=rel(1.75))) +
  scale_y_continuous(expand = c(0,0), lim = c(0,max(dat$predictionScore)))

## Only for clustering annotation:
dat <- subset(dat, Annotation == 'Clustering')
ggplot(dat, aes(x=celltype,
                y=predictionScore)) + 
  theme_bw() +
  geom_boxplot(fill = 'grey90', width=0.6, position = position_dodge(width=.8), 
               coef = 0, color = 'black', lwd = 0.6, outlier.shape = NA) +
  stat_summary(fun.y="mean", geom="point", size=3,   
               position=position_dodge(width=0.8),                                                 
               color="red") +
  theme(axis.text.x=element_text(size=rel(1.75)), 
        axis.text.y=element_text(size=rel(1.75)))+ 
  labs(x="Region",y="Prediction score",
       title="Cell type prediction scores in follicles (Jejunum)") +
  theme(axis.title.x = element_text(size = 16),
        axis.text.x=element_text(size=12), 
        axis.title.y = element_text(size = 16),
        axis.text.y=element_text(size=12),
        legend.title = element_text (size=16),
        title = element_text(size = 16)) +
  theme(axis.text.x=element_text(angle = 90, size=rel(1.75), hjust=0.9,vjust=0.2), axis.text.y=element_text(size=rel(1.75))) +
  scale_y_continuous(expand = c(0,0)) + coord_cartesian(ylim = c(0, 1))

# Make example plot of cycling B cell prediction scores ----
mapj <- read_xlsx('/home/Jayne.Wiarda/SI_PP_SC_ST/MappingPrediction/STomics_CellTypePredictions_ClusteringAnnotation_Jejunum.xlsx')
mapi <- read_xlsx('/home/Jayne.Wiarda/SI_PP_SC_ST/MappingPrediction/STomics_CellTypePredictions_ClusteringAnnotation_Ileum.xlsx')
mapClus <- as.data.frame(rbind(mapj, mapi))
rownames(mapClus) <- mapClus$CellBarcode
mapClus <- subset(mapClus, select=-c(CellBarcode))

mapj <- read_xlsx('/home/Jayne.Wiarda/SI_PP_SC_ST/MappingPrediction/STomics_CellTypePredictions_ManualAnnotation_Jejunum.xlsx')
mapi <- read_xlsx('/home/Jayne.Wiarda/SI_PP_SC_ST/MappingPrediction/STomics_CellTypePredictions_ManualAnnotation_Ileum.xlsx')
mapMan <- as.data.frame(rbind(mapj, mapi))
rownames(mapMan) <- mapMan$CellBarcode
mapMan <- subset(mapMan, select=-c(CellBarcode))

Idents(st) <- st$Region_Clust
stClus <- subset(st, idents = c('Follicle', 'Crypt', 'Villus', 'IFZ/PFZ')) # remove muscularis spots from spatial dataset since they aren't represented in scRNA-seq data (muscularis was removed from scRNA-seq samples)
Idents(seu) <- seu$Region_ManAnn
stMan <- subset(seu, idents = c('Follicle', 'Crypt', 'Villus', 'IFZ/PFZ')) # remove muscularis spots from spatial dataset since they aren't represented in scRNA-seq data (muscularis was removed from scRNA-seq samples)

stClus <- AddMetaData(stClus, mapClus)
stMan <- AddMetaData(stMan, mapMan)

sub1 <- as.data.frame(stClus@meta.data)
colnames(sub1)
sub1 <- sub1[,c(10,17)]
colnames(sub1) <- c('Region', 'Cycling B cells')
sub1 <- gather(sub1, key="Increment", value="predictionScore", 2)
sub1$Annotation <- rep('Clustering', nrow(sub1))
sub1$Region <- factor(sub1$Region, levels = c('Crypt', 'Villus', 'IFZ/PFZ', 'Follicle'))

sub2 <- as.data.frame(stMan@meta.data)
colnames(sub2)
sub2 <- sub2[,c(9,16)]
colnames(sub2) <- c('Region', 'Cycling B cells')
sub2 <- gather(sub2, key="Increment", value="predictionScore", 2)
sub2$Annotation <- rep('Manual', nrow(sub2))
sub2$Region <- factor(sub2$Region, levels = c('Crypt', 'Villus', 'IFZ/PFZ', 'Follicle'))

sub <- rbind(sub1, sub2)
sub$Annotation <- factor(sub$Annotation, levels = c('Manual', 'Clustering'))

dat <- subset(sub, Increment == 'Cycling B cells')
ggplot(dat, aes(x=Region,
                y=predictionScore,
                fill = Annotation)) + 
  theme_bw() +
  geom_boxplot(aes(fill=Annotation), width=0.6, position = position_dodge(width=.8), 
               coef = 0, color = 'black', lwd = 0.6, outlier.shape = NA) +
  geom_point(shape=16, position = position_jitterdodge(jitter.width=.3, dodge.width = 0.8), size = .1, color = 'grey30') + 
  scale_fill_manual(values=c('grey90', 'grey65')) + 
  stat_summary(fun.y="mean", geom="point", size=4,
               position=position_dodge(width=0.8), 
               color="red")+
  theme(axis.text.x=element_text(size=rel(1.75)), 
        axis.text.y=element_text(size=rel(1.75)))+ 
  labs(x="Region",y="Prediction score",
       title="Cycling B cell Prediction Scores") +
  theme(axis.title.x = element_text(size = 16),
        axis.text.x=element_text(size=12), 
        axis.title.y = element_text(size = 16),
        axis.text.y=element_text(size=12),
        legend.title = element_text (size=16),
        title = element_text(size = 16)) +
  theme(axis.text.x=element_text(size=rel(1.75)), axis.text.y=element_text(size=rel(1.75))) +
  scale_y_continuous(expand = c(0,0), lim = c(0,1))

ggplot(dat, aes(x=Region,
                y=predictionScore)) + 
  theme_bw() +
  geom_boxplot(fill = 'grey90', width=0.6, position = position_dodge(width=.8), 
               coef = 0, color = 'black', lwd = 0.6, outlier.shape = NA) +
  stat_summary(fun.y="mean", geom="point", size=3,   
               position=position_dodge(width=0.8),                                                 
               color="red") +
  theme(axis.text.x=element_text(size=rel(1.75)), 
        axis.text.y=element_text(size=rel(1.75)))+ 
  labs(x="Region",y="Prediction score",
       title="Cell type prediction scores in follicles (Jejunum)") +
  theme(axis.title.x = element_text(size = 16),
        axis.text.x=element_text(size=12), 
        axis.title.y = element_text(size = 16),
        axis.text.y=element_text(size=12),
        legend.title = element_text (size=16),
        title = element_text(size = 16)) +
  theme(axis.text.x=element_text(angle = 90, size=rel(1.75), hjust=0.9,vjust=0.2), axis.text.y=element_text(size=rel(1.75))) +
  scale_y_continuous(expand = c(0,0)) + coord_cartesian(ylim = c(0, 1))

sessionInfo()

#R version 4.2.2 Patched (2022-11-10 r83330)
#Platform: x86_64-pc-linux-gnu (64-bit)
#Running under: Ubuntu 22.04.1 LTS

#Matrix products: default
#BLAS:   /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3
#LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/libopenblasp-r0.3.20.so

#Random number generation:
#  RNG:     Mersenne-Twister 
#Normal:  Inversion 
#Sample:  Rounding 

#locale:
#  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C               LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8     LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8    LC_PAPER=en_US.UTF-8       LC_NAME=C                  LC_ADDRESS=C               LC_TELEPHONE=C             LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       

#attached base packages:
#  [1] grid      stats4    stats     graphics  grDevices utils     datasets  methods   base     

#other attached packages:
#  [1] scales_1.2.1          viridis_0.6.2         viridisLite_0.4.1     ggplot2_3.3.6         dplyr_1.0.10          writexl_1.4.0         SeuratDisk_0.0.0.9020 sp_1.5-0              SeuratObject_4.1.1    Seurat_4.1.1         

#loaded via a namespace (and not attached):
#  [1] scattermore_0.8             ragg_1.2.2                  pkgmaker_0.32.2             tidyr_1.2.1                 knitr_1.40                  bit64_4.0.5                 irlba_2.3.5                 DelayedArray_0.22.0         data.table_1.14.2           rpart_4.1.16                KEGGREST_1.36.3            
#[12] RCurl_1.98-1.8              doParallel_1.0.17           generics_0.1.3              BiocGenerics_0.42.0         ScaledMatrix_1.4.1          cowplot_1.1.1               RSQLite_2.2.17              RANN_2.6.1                  miloR_1.4.0                 future_1.28.0               bit_4.0.4                  
#[23] phylobase_0.8.10            spatstat.data_2.2-0         xml2_1.3.3                  httpuv_1.6.6                SummarizedExperiment_1.26.1 assertthat_0.2.1            xfun_0.33                   hms_1.1.2                   evaluate_0.16               promises_1.2.0.1            fansi_1.0.3                
#[34] progress_1.2.2              igraph_1.3.4                DBI_1.1.3                   htmlwidgets_1.5.4           spatstat.geom_2.4-0         purrr_0.3.4                 ellipsis_0.3.2              annotate_1.74.0             gridBase_0.4-7              locfdr_1.1-8                deldir_1.0-6               
#[45] MatrixGenerics_1.8.1        vctrs_0.4.1                 SingleCellExperiment_1.18.0 Biobase_2.56.0              ROCR_1.0-11                 abind_1.4-5                 cachem_1.0.6                withr_2.5.0                 ggforce_0.3.4               progressr_0.11.0            sctransform_0.3.4          
#[56] prettyunits_1.1.1           goftest_1.2-3               softImpute_1.4-1            cluster_2.1.4               ape_5.6-2                   lazyeval_0.2.2              crayon_1.5.1                genefilter_1.78.0           hdf5r_1.3.5                 labeling_0.4.2              edgeR_3.38.4               
#[67] pkgconfig_2.0.3             tweenr_2.0.2                GenomeInfoDb_1.32.4         nlme_3.1-159                vipor_0.4.5                 rlang_1.0.6                 globals_0.16.1              lifecycle_1.0.2             miniUI_0.1.1.1              registry_0.5-1              rsvd_1.0.5                 
#[78] cellranger_1.1.0            polyclip_1.10-0             matrixStats_0.62.0          lmtest_0.9-40               rngtools_1.5.2              Matrix_1.5-1                Rhdf5lib_1.18.2             zoo_1.8-10                  beeswarm_0.4.0              ggridges_0.5.3              png_0.1-7                  
#[89] bitops_1.0-7                rncl_0.8.6                  KernSmooth_2.23-20          rhdf5filters_1.8.0          Biostrings_2.64.1           blob_1.2.3                  stringr_1.4.1               zinbwave_1.18.0             parallelly_1.32.1           spatstat.random_2.2-0       S4Vectors_0.34.0           
#[100] beachmat_2.12.0             memoise_2.0.1               magrittr_2.0.3              plyr_1.8.7                  ica_1.0-3                   howmany_0.3-1               zlibbioc_1.42.0             compiler_4.2.2              RColorBrewer_1.1-3          fitdistrplus_1.1-8          cli_3.4.0                  
#[111] ade4_1.7-19                 XVector_0.36.0              listenv_0.8.0               patchwork_1.1.2             pbapply_1.5-0               MASS_7.3-58.1               mgcv_1.8-40                 tidyselect_1.1.2            stringi_1.7.8               textshaping_0.3.6           yaml_2.3.5                 
#[122] BiocSingular_1.12.0         locfit_1.5-9.6              ggrepel_0.9.1               tools_4.2.2                 future.apply_1.9.1          parallel_4.2.2              rstudioapi_0.14             uuid_1.1-0                  foreach_1.5.2               RNeXML_2.4.7                gridExtra_2.3              
#[133] farver_2.1.1                Rtsne_0.16                  ggraph_2.0.6                digest_0.6.29               rgeos_0.5-9                 shiny_1.7.2                 Rcpp_1.0.9                  GenomicRanges_1.48.0        later_1.3.0                 RcppAnnoy_0.0.19            httr_1.4.4                 
#[144] AnnotationDbi_1.58.0        kernlab_0.9-31              colorspace_2.0-3            XML_3.99-0.10               tensor_1.5                  reticulate_1.26             clusterExperiment_2.16.0    IRanges_2.30.1              splines_4.2.2               uwot_0.1.14                 spatstat.utils_2.3-1       
#[155] graphlayouts_0.8.1          plotly_4.10.0               systemfonts_1.0.4           xtable_1.8-4                jsonlite_1.8.0              tidygraph_1.2.2             R6_2.5.1                    pillar_1.8.1                htmltools_0.5.3             mime_0.12                   NMF_0.24.0                 
#[166] glue_1.6.2                  fastmap_1.1.0               BiocParallel_1.30.3         BiocNeighbors_1.14.0        codetools_0.2-18            utf8_1.2.2                  lattice_0.20-45             spatstat.sparse_2.1-1       tibble_3.1.8                ggbeeswarm_0.6.0            leiden_0.4.3               
#[177] gtools_3.9.3                survival_3.4-0              limma_3.52.3                rmarkdown_2.16              munsell_0.5.0               rhdf5_2.40.0                GenomeInfoDbData_1.2.8      iterators_1.0.14            HDF5Array_1.24.2            reshape2_1.4.4              gtable_0.3.1               
#[188] spatstat.core_2.4-4   