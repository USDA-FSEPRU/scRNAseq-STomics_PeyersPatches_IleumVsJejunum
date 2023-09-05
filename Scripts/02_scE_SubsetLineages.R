library(Seurat)
library(SeuratObject)
library(ggplot2)
library(SeuratDisk)

## B lineage

#Subset data:

seu <- LoadH5Seurat('/home/Jayne.Wiarda/SI_PP_SC_ST/AnnotatedSeurat/AllSamples_annotated.h5seurat')
DefaultAssay(seu) <- 'RNA'
Idents(seu) <- seu$celltype
seu <- subset(seu, idents = c('Cycling B cells', "Activated B cells",
                              "Transitioning B cells",
                              "Resting B cells",
                              "Antibody-secreting cells"))
DefaultAssay(seu) <- 'RNA'
counts <- as.data.frame(seu[['RNA']]@counts)
keep <- rowSums(counts) > 0
keep <- rownames(counts[keep,])
seu <- DietSeurat(seu, 
                  counts = TRUE,
                  data = TRUE,
                  scale.data = FALSE, # remove the scaled data
                  dimreducs = NULL,
                  features = keep, # keep only genes with non-zero counts across all cells
                  assays = 'RNA') # keep only RNA assay and remove SCT and integrated

#Integrate samples:

seu.list <- SplitObject(seu, split.by = "orig.ident") # split by sample IDs
for (i in 1:length(seu.list)) { # normalize data using SCTransform method
  seu.list[[i]] <- SCTransform(seu.list[[i]], 
                               return.only.var.genes = FALSE, 
                               verbose = TRUE) 
}
seu.features <- SelectIntegrationFeatures(seu.list, # select the genes to use for integration
                                          verbose = TRUE) 
seu.list <- PrepSCTIntegration(seu.list, 
                               anchor.features = seu.features,
                               verbose = TRUE)
seu.anchors <- FindIntegrationAnchors(seu.list, # identify anchors for integration from top 30 data dimensions
                                      normalization.method = "SCT", 
                                      anchor.features = seu.features, 
                                      dims = 1:30)
seu.integrated <- IntegrateData(seu.anchors, # integrate data
                                normalization.method = "SCT", 
                                dims = 1:30)

#Dimensionality reduction:

seu.integrated <- RunPCA(seu.integrated, # run PCA analysis for 50 dimensions of the data
                         npcs = 50, 
                         verbose = TRUE) 
pct <- seu.integrated[["pca"]]@stdev / sum(seu.integrated[["pca"]]@stdev) * 100 # find standard deviation for each PC
cumu <- cumsum(pct) # find cumulative percentages for PCs
co1 <- which(cumu > 90 & pct < 5)[1] # find PC representing cumulative percent >90% and less than 5% associated with the single PC
co1 # list PC
co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1 # find last PC where change in percent variation is more than 0.1%
co2 # list PC
pcs <- min(co1, co2) # find the minimum PC from the 2 methods used above
pcs # list PC
plot_df <- data.frame(pct = pct, # put PC values into dataframe for plotting
                      cumu = cumu, 
                      rank = 1:length(pct))
#ggplot(plot_df, aes(cumu, pct, label = rank, color = rank > pcs)) + # visualize PCs to use in elbow plot
#  geom_text() + 
#  geom_vline(xintercept = 90, color = "grey") + 
#  geom_hline(yintercept = min(pct[pct > 5]), color = "grey") +
#  theme_bw()
PCdims <- 1:pcs # use the minimum PC from the quantitative method above to set the PCs for subsequent steps
rm(pct, cumu, co1, co2, pcs)
seu.integrated <- RunTSNE(seu.integrated, 
                          dims = PCdims, # use our calculated number of PCs
                          reduction = "pca", 
                          dim_embed = 3, # calculate 3 plot dimensions in case we want to try 3D plotting later
                          assay = "SCT") # create tSNE plot 
seu.integrated <- RunUMAP(seu.integrated, 
                          dims = PCdims, # use our calculated number of PCs
                          reduction = "pca", 
                          dim_embed = 3, # calculate 3 plot dimensions in case we want to try 3D plotting later
                          min.dist = 0.5,
                          spread = 0.2,
                          assay = "SCT") # create UMAP

#Normalize and scale counts:

seu.integrated <- NormalizeData(seu.integrated,  # normalize the RNA counts data per cell
                                normalization.method = "LogNormalize", 
                                scale.factor = 10000, 
                                assay = "RNA")
seu.integrated <- ScaleData(seu.integrated, # scale the RNA counts data relative to other cells
                            assay = "RNA")
seu.integrated <- ScaleData(seu.integrated, # scale the SCT counts data relative to other cells
                            assay = "SCT")
#dim(seu.integrated[["RNA"]]@scale.data) # see that all genes are scaled in RNA assay now
#dim(seu.integrated[["SCT"]]@scale.data) # see that all genes are scaled in SCT assay now

#Rearrange identity levels:

Idents(seu.integrated) <- seu.integrated$celltype
levels(seu.integrated) <- c('Activated B cells', 'Cycling B cells', 'Resting B cells', 'Transitioning B cells', 'Antibody-secreting cells')
seu.integrated$celltype <- Idents(seu.integrated)

#Plotting:

DefaultAssay(seu.integrated) <- 'SCT'
DimPlot(seu.integrated, reduction = 'umap', group.by = 'celltype', cols = c('cyan4', 'gold3', 'chartreuse4', 'deeppink4', 'sandybrown'))
DotPlot(seu.integrated,
        features = c('CD79A', 'CD79B', 'CD19', 'MS4A1', 'PAX5', # B cells
                     'JCHAIN', 'XBP1', 'PRDM1', 'IRF4', # ASCs
                     'CD69', 'CD83', 'SLA-DQB1', 'SLA-DRA', # B activation
                     'GPR183', 'CCR7', 'KLF2', 'SELL', 'FCER2', 'CD40', # resting B
                     'AICDA', 'CD86', 'BCL6', # follicular B
                     'PCLAF', 'BIRC5', 'TOP2A', 'STMN1'), # cycling B
        col.min = -.5, col.max = 1
) + RotatedAxis() + 
  scale_colour_gradient2(low="khaki1", mid="darkseagreen", high="darkslategray4")

#Save data:

SaveH5Seurat(seu.integrated, '/home/Jayne.Wiarda/SI_PP_SC_ST/AnnotatedSeurat/B.h5seurat', overwrite = TRUE)

## T/ILC lineage

#Subset data:

seu <- LoadH5Seurat('/home/Jayne.Wiarda/SI_PP_SC_ST/AnnotatedSeurat/AllSamples_annotated.h5seurat')
DefaultAssay(seu) <- 'RNA'
Idents(seu) <- seu$celltype
seu <- subset(seu, idents = c('Cycling CD4 ab T cells', 'Cycling CD8 ab T cells', 'Cycling gd T cells',
                              'Cytotoxic CD8 ab T cells', 'Cytotoxic gd T cells', 'Cytotoxic group 1 ILCs', 
                              'Non-naive CD8 ab T cells', 'Non-naive gd T cells', 'Non-naive group 1 ILCs', 
                              'SELLhi gd T cells', 'CD2neg GD T cells', 
                              'Naive CD4/CD8 ab T cells', 'Non-naive CD4 ab T cells', 'Follicular CD4 ab T cells', 'Group 3 ILCs'))
DefaultAssay(seu) <- 'RNA'
counts <- as.data.frame(seu[['RNA']]@counts)
keep <- rowSums(counts) > 0
keep <- rownames(counts[keep,])
seu <- DietSeurat(seu, 
                  counts = TRUE,
                  data = TRUE,
                  scale.data = FALSE, # remove the scaled data
                  dimreducs = NULL,
                  features = keep, # keep only genes with non-zero counts across all cells
                  assays = 'RNA') # keep only RNA assay and remove SCT and integrated

#Integrate samples:

seu.list <- SplitObject(seu, split.by = "orig.ident") # split by sample IDs
for (i in 1:length(seu.list)) { # normalize data using SCTransform method
  seu.list[[i]] <- SCTransform(seu.list[[i]], 
                               return.only.var.genes = FALSE, 
                               verbose = TRUE) 
}
seu.features <- SelectIntegrationFeatures(seu.list, # select the genes to use for integration
                                          verbose = TRUE) 
seu.list <- PrepSCTIntegration(seu.list, 
                               anchor.features = seu.features,
                               verbose = TRUE)
seu.anchors <- FindIntegrationAnchors(seu.list, # identify anchors for integration from top 30 data dimensions
                                      normalization.method = "SCT", 
                                      anchor.features = seu.features, 
                                      dims = 1:30)
seu.integrated <- IntegrateData(seu.anchors, # integrate data
                                normalization.method = "SCT", 
                                dims = 1:30)

#Dimensionality reduction:

seu.integrated <- RunPCA(seu.integrated, # run PCA analysis for 50 dimensions of the data
                         npcs = 50, 
                         verbose = TRUE) 
pct <- seu.integrated[["pca"]]@stdev / sum(seu.integrated[["pca"]]@stdev) * 100 # find standard deviation for each PC
cumu <- cumsum(pct) # find cumulative percentages for PCs
co1 <- which(cumu > 90 & pct < 5)[1] # find PC representing cumulative percent >90% and less than 5% associated with the single PC
co1 # list PC
co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1 # find last PC where change in percent variation is more than 0.1%
co2 # list PC
pcs <- min(co1, co2) # find the minimum PC from the 2 methods used above
pcs # list PC
plot_df <- data.frame(pct = pct, # put PC values into dataframe for plotting
                      cumu = cumu, 
                      rank = 1:length(pct))
#ggplot(plot_df, aes(cumu, pct, label = rank, color = rank > pcs)) + # visualize PCs to use in elbow plot
#  geom_text() + 
#  geom_vline(xintercept = 90, color = "grey") + 
#  geom_hline(yintercept = min(pct[pct > 5]), color = "grey") +
#  theme_bw()
PCdims <- 1:pcs # use the minimum PC from the quantitative method above to set the PCs for subsequent steps
rm(pct, cumu, co1, co2, pcs)
seu.integrated <- RunTSNE(seu.integrated, 
                          dims = PCdims, # use our calculated number of PCs
                          reduction = "pca", 
                          dim_embed = 3, # calculate 3 plot dimensions in case we want to try 3D plotting later
                          assay = "SCT") # create tSNE plot 
seu.integrated <- RunUMAP(seu.integrated, 
                          dims = PCdims, # use our calculated number of PCs
                          reduction = "pca", 
                          dim_embed = 3, # calculate 3 plot dimensions in case we want to try 3D plotting later
                          min.dist = 0.5,
                          spread = 0.2,
                          assay = "SCT") # create UMAP

#Normalize and scale counts:

seu.integrated <- NormalizeData(seu.integrated,  # normalize the RNA counts data per cell
                                normalization.method = "LogNormalize", 
                                scale.factor = 10000, 
                                assay = "RNA")
seu.integrated <- ScaleData(seu.integrated, # scale the RNA counts data relative to other cells
                            assay = "RNA") # scales all genes instead of just highly variable
seu.integrated <- ScaleData(seu.integrated, # scale the SCT counts data relative to other cells
                            assay = "SCT")
#dim(seu.integrated[["RNA"]]@scale.data) # see that all genes are scaled in RNA assay now
#dim(seu.integrated[["SCT"]]@scale.data) # see that all genes are scaled in SCT assay now

#Rearrange identity levels:

Idents(seu.integrated) <- seu.integrated$celltype
levels(seu.integrated) <- c('Cycling CD4 ab T cells', 'Cycling CD8 ab T cells', 'Cycling gd T cells',
                            'Cytotoxic CD8 ab T cells', 'Cytotoxic gd T cells', 'Cytotoxic group 1 ILCs', 
                            'Non-naive CD8 ab T cells', 'Non-naive gd T cells', 'Non-naive group 1 ILCs', 
                            'SELLhi gd T cells', 'CD2neg GD T cells', 
                            'Naive CD4/CD8 ab T cells', 'Non-naive CD4 ab T cells', 'Follicular CD4 ab T cells', 'Group 3 ILCs')
seu.integrated$celltype <- Idents(seu.integrated)

#Plotting:

DefaultAssay(seu.integrated) <- 'SCT'
DimPlot(seu.integrated, reduction = 'umap', group.by = 'celltype', cols = c('cornflowerblue', 'navy', 'lightpink', 'gold3', 'salmon', 'deepskyblue2', 'chartreuse4', 'deeppink4', 'tan4', 'mediumpurple1', 'darkgreen', 'gray50', 'cyan4', 'sandybrown', 'darkmagenta'))
DotPlot(seu.integrated,
        features = c('CD3E', 'CD3G', 'CD247', # T cell
                     'CD4', 'CD8B', 'TRDC', 'CD2', 'CD8A', # T cell subsets
                     'PCLAF', 'BIRC5', 'TOP2A', 'STMN1', # cycling
                     'CCL5', 'ITGAE', # effector/resident
                     'GZMA-16903', 'GZMB', 'GNLY', # cytotoxic
                     'CTSW', 'XCL1', 'SLA-DRA', 'SLA-DQB1', 'CCR9', 'KLRK1', # activation
                     'FCER1G', 'KLRG1', 'ITGB1', 'ITGB7', 'SELL', # SELLhi gd
                     'ID3', 'RHEX', 'BLK', 'SAMSN1', 'IL26', # CD2- gd
                     'CCR7', 'S1PR1', 'LEF1', 'KLF2', # naive ab
                     'ICOS', 'CTLA4', 'CD40LG', 'IL10', # Non-naive + follicular CD4
                     'CD52', 'IFITM3', 'GPR183', # non-naive CD4
                     'PDCD1', 'CXCR4', 'CD69', # follicular CD4
                     'LTB', 'ID2', 'KIT', 'IL7R', 'IL22', 'KLRB1', 'RORC', 'CXCL8'), # group 3 ILC
        col.min = -.5, col.max = 1
) + RotatedAxis() + 
  scale_colour_gradient2(low="khaki1", mid="darkseagreen", high="darkslategray4")

#Save data:

SaveH5Seurat(seu.integrated, '/home/Jayne.Wiarda/SI_PP_SC_ST/AnnotatedSeurat/TILC.h5seurat', overwrite = TRUE)

## Myeloid lineage

#Subset data:

seu <- LoadH5Seurat('/home/Jayne.Wiarda/SI_PP_SC_ST/AnnotatedSeurat/AllSamples_annotated.h5seurat')
DefaultAssay(seu) <- 'RNA'
Idents(seu) <- seu$celltype
seu <- subset(seu, idents = c('Dendritic cells', 'Macrophages', 'Mast cells'))
DefaultAssay(seu) <- 'RNA'
counts <- as.data.frame(seu[['RNA']]@counts)
keep <- rowSums(counts) > 0
keep <- rownames(counts[keep,])
seu <- DietSeurat(seu, 
                  counts = TRUE,
                  data = TRUE,
                  scale.data = FALSE, # remove the scaled data
                  dimreducs = NULL,
                  features = keep, # keep only genes with non-zero counts across all cells
                  assays = 'RNA') # keep only RNA assay and remove SCT and integrated

#Integrate samples:

seu.list <- SplitObject(seu, split.by = "orig.ident") # split by sample IDs
for (i in 1:length(seu.list)) { # normalize data using SCTransform method
  seu.list[[i]] <- SCTransform(seu.list[[i]], 
                               return.only.var.genes = FALSE, 
                               verbose = TRUE) 
}
seu.features <- SelectIntegrationFeatures(seu.list, # select the genes to use for integration
                                          verbose = TRUE) 
seu.list <- PrepSCTIntegration(seu.list, 
                               anchor.features = seu.features,
                               verbose = TRUE)
seu.anchors <- FindIntegrationAnchors(seu.list, # identify anchors for integration from top 30 data dimensions
                                      normalization.method = "SCT", 
                                      anchor.features = seu.features, 
                                      dims = 1:30)
seu.integrated <- IntegrateData(seu.anchors, # integrate data
                                normalization.method = "SCT", 
                                dims = 1:30,
                                k.weight = 55) # reduce to number of cells in smallest sample

#Dimensionality reduction:

seu.integrated <- RunPCA(seu.integrated, # run PCA analysis for 50 dimensions of the data
                         npcs = 50, 
                         verbose = TRUE) 
pct <- seu.integrated[["pca"]]@stdev / sum(seu.integrated[["pca"]]@stdev) * 100 # find standard deviation for each PC
cumu <- cumsum(pct) # find cumulative percentages for PCs
co1 <- which(cumu > 90 & pct < 5)[1] # find PC representing cumulative percent >90% and less than 5% associated with the single PC
co1 # list PC
co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1 # find last PC where change in percent variation is more than 0.1%
co2 # list PC
pcs <- min(co1, co2) # find the minimum PC from the 2 methods used above
pcs # list PC
plot_df <- data.frame(pct = pct, # put PC values into dataframe for plotting
                      cumu = cumu, 
                      rank = 1:length(pct))
#ggplot(plot_df, aes(cumu, pct, label = rank, color = rank > pcs)) + # visualize PCs to use in elbow plot
#  geom_text() + 
#  geom_vline(xintercept = 90, color = "grey") + 
#  geom_hline(yintercept = min(pct[pct > 5]), color = "grey") +
#  theme_bw()
PCdims <- 1:pcs # use the minimum PC from the quantitative method above to set the PCs for subsequent steps
rm(pct, cumu, co1, co2, pcs)
seu.integrated <- RunTSNE(seu.integrated, 
                          dims = PCdims, # use our calculated number of PCs
                          reduction = "pca", 
                          dim_embed = 3, # calculate 3 plot dimensions in case we want to try 3D plotting later
                          assay = "SCT") # create tSNE plot 
seu.integrated <- RunUMAP(seu.integrated, 
                          dims = PCdims, # use our calculated number of PCs
                          reduction = "pca", 
                          dim_embed = 3, # calculate 3 plot dimensions in case we want to try 3D plotting later
                          min.dist = 0.5,
                          spread = 0.2,
                          assay = "SCT") # create UMAP

#Normalize and scale counts:

seu.integrated <- NormalizeData(seu.integrated,  # normalize the RNA counts data per cell
                                normalization.method = "LogNormalize", 
                                scale.factor = 10000, 
                                assay = "RNA")
seu.integrated <- ScaleData(seu.integrated, # scale the RNA counts data relative to other cells
                            assay = "RNA") # scales all genes instead of just highly variable
seu.integrated <- ScaleData(seu.integrated, # scale the SCT counts data relative to other cells
                            assay = "SCT")
#dim(seu.integrated[["RNA"]]@scale.data) # see that all genes are scaled in RNA assay now
#dim(seu.integrated[["SCT"]]@scale.data) # see that all genes are scaled in SCT assay now

#Rearrange identity levels:

Idents(seu.integrated) <- seu.integrated$celltype
levels(seu.integrated) <- c('Dendritic cells', 'Macrophages', 'Mast cells')
seu.integrated$celltype <- Idents(seu.integrated)

#Plotting:

DefaultAssay(seu.integrated) <- 'SCT'
DimPlot(seu.integrated, reduction = 'umap', group.by = 'celltype', cols = c('cyan4', 'gold3', 'chartreuse4'))
DotPlot(seu.integrated,
        features = c('FLT3', 'SLA-DRA', 'SLA-DQB1', # DC
                     'SIRPA', 'CD68', 'CXCL2', 'C1QA', 'C1QB', 'C1QC', # macrophage
                     'ICAM1', 'CSF2RB', 'MS4A2', 'FCER1A'), # mast cell
        col.min = -.5, col.max = 1
) + RotatedAxis() + 
  scale_colour_gradient2(low="khaki1", mid="darkseagreen", high="darkslategray4")

#Save data:

SaveH5Seurat(seu.integrated, '/home/Jayne.Wiarda/SI_PP_SC_ST/AnnotatedSeurat/Myeloid.h5seurat', overwrite = TRUE)

## Epithelial lineage

#Subset data:

seu <- LoadH5Seurat('/home/Jayne.Wiarda/SI_PP_SC_ST/AnnotatedSeurat/AllSamples_annotated.h5seurat')
DefaultAssay(seu) <- 'RNA'
Idents(seu) <- seu$celltype
seu <- subset(seu, idents = c('Crypt cells', 'Enterocytes', 'BEST4 enterocytes', 'Goblet cells', 'NEUROD1lo EE cells', 'NEUROD1hi EE cells'))
DefaultAssay(seu) <- 'RNA'
counts <- as.data.frame(seu[['RNA']]@counts)
keep <- rowSums(counts) > 0
keep <- rownames(counts[keep,])
seu <- DietSeurat(seu, 
                  counts = TRUE,
                  data = TRUE,
                  scale.data = FALSE, # remove the scaled data
                  dimreducs = NULL,
                  features = keep, # keep only genes with non-zero counts across all cells
                  assays = 'RNA') # keep only RNA assay and remove SCT and integrated

#Integrate samples:

seu.list <- SplitObject(seu, split.by = "orig.ident") # split by sample IDs
for (i in 1:length(seu.list)) { # normalize data using SCTransform method
  seu.list[[i]] <- SCTransform(seu.list[[i]], 
                               return.only.var.genes = FALSE, 
                               verbose = TRUE) 
}
seu.features <- SelectIntegrationFeatures(seu.list, # select the genes to use for integration
                                          verbose = TRUE) 
seu.list <- PrepSCTIntegration(seu.list, 
                               anchor.features = seu.features,
                               verbose = TRUE)
seu.anchors <- FindIntegrationAnchors(seu.list, # identify anchors for integration from top 30 data dimensions
                                      normalization.method = "SCT", 
                                      anchor.features = seu.features, 
                                      k.filter = 11, # reduce to 11 since smallest sample has 12 cells
                                      k.score = 11, # reduce to 11 since smallest sample has 12 cells
                                      dims = 1:11) # reduce to 11 since smallest sample has 12 cells
seu.integrated <- IntegrateData(seu.anchors, # integrate data
                                normalization.method = "SCT", 
                                k.weight = 11, # reduce to 11 since smallest sample has 12 cells
                                dims = 1:11) # reduce to 11 since smallest sample has 12 cells

#Dimensionality reduction:

seu.integrated <- RunPCA(seu.integrated, # run PCA analysis for 50 dimensions of the data
                         npcs = 50, 
                         verbose = TRUE) 
pct <- seu.integrated[["pca"]]@stdev / sum(seu.integrated[["pca"]]@stdev) * 100 # find standard deviation for each PC
cumu <- cumsum(pct) # find cumulative percentages for PCs
co1 <- which(cumu > 90 & pct < 5)[1] # find PC representing cumulative percent >90% and less than 5% associated with the single PC
co1 # list PC
co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1 # find last PC where change in percent variation is more than 0.1%
co2 # list PC
pcs <- min(co1, co2) # find the minimum PC from the 2 methods used above
pcs # list PC
plot_df <- data.frame(pct = pct, # put PC values into dataframe for plotting
                      cumu = cumu, 
                      rank = 1:length(pct))
#ggplot(plot_df, aes(cumu, pct, label = rank, color = rank > pcs)) + # visualize PCs to use in elbow plot
#  geom_text() + 
#  geom_vline(xintercept = 90, color = "grey") + 
#  geom_hline(yintercept = min(pct[pct > 5]), color = "grey") +
#  theme_bw()
PCdims <- 1:pcs # use the minimum PC from the quantitative method above to set the PCs for subsequent steps
rm(pct, cumu, co1, co2, pcs)
seu.integrated <- RunTSNE(seu.integrated, 
                          dims = PCdims, # use our calculated number of PCs
                          reduction = "pca", 
                          dim_embed = 3, # calculate 3 plot dimensions in case we want to try 3D plotting later
                          assay = "SCT") # create tSNE plot 
seu.integrated <- RunUMAP(seu.integrated, 
                          dims = PCdims, # use our calculated number of PCs
                          reduction = "pca", 
                          dim_embed = 3, # calculate 3 plot dimensions in case we want to try 3D plotting later
                          min.dist = 0.5,
                          spread = 0.2,
                          assay = "SCT") # create UMAP

#Normalize and scale counts:

seu.integrated <- NormalizeData(seu.integrated,  # normalize the RNA counts data per cell
                                normalization.method = "LogNormalize", 
                                scale.factor = 10000, 
                                assay = "RNA")
seu.integrated <- ScaleData(seu.integrated, # scale the RNA counts data relative to other cells
                            assay = "RNA") # scales all genes instead of just highly variable
seu.integrated <- ScaleData(seu.integrated, # scale the SCT counts data relative to other cells
                            assay = "SCT")
#dim(seu.integrated[["RNA"]]@scale.data) # see that all genes are scaled in RNA assay now
#dim(seu.integrated[["SCT"]]@scale.data) # see that all genes are scaled in SCT assay now

#Rearrange identity levels:

Idents(seu.integrated) <- seu.integrated$celltype
levels(seu.integrated) <- c('Crypt cells', 'Enterocytes', 'BEST4 enterocytes', 'Goblet cells', 'NEUROD1lo EE cells', 'NEUROD1hi EE cells')
seu.integrated$celltype <- Idents(seu.integrated)

#Plotting:

DefaultAssay(seu.integrated) <- 'SCT'
DimPlot(seu.integrated, reduction = 'umap', group.by = 'celltype', cols = c('cyan4', 'gold3', 'chartreuse4', 'deeppink4', 'sandybrown', 'navy'))
DotPlot(seu.integrated,
        features = c('RPL5', 'RPS6', 'EEF1B2', 'OLFM4', 'PIGR', 'LYZ', # crypts
                     'FABP2', 'FABP1', 'CLCA4', 'SLC5A1', 'SI', 'ACE2',  # enterocyte
                     'GUCA2A', 'GUCA2B', 'BEST4', 'CFTR', 'NOTCH2', 'OTOP2', # best4 enterocyte
                     'TFF3', 'REG4', 'CLCA1', 'SPINK4', 'MUC2', 'CXCL8', # goblet
                     'PYY', 'GAST', 'SST', 'CCK', 'TTR', 'NTS', # NEUROD1lo EE
                     'NEUROD1', 'CHGA', 'CHGB', 'KRT7', 'SCT', 'PENK'), # NEUROD1hi EE
        col.min = -.5, col.max = 1
) + RotatedAxis() + 
  scale_colour_gradient2(low="khaki1", mid="darkseagreen", high="darkslategray4")

#Save data:

SaveH5Seurat(seu.integrated, '/home/Jayne.Wiarda/SI_PP_SC_ST/AnnotatedSeurat/Epithelial.h5seurat', overwrite = TRUE)

## Stromal lineage

#We can skip straight to plotting since we already subsetted and re-processed stromal cells

seu.integrated <- LoadH5Seurat('/home/Jayne.Wiarda/SI_PP_SC_ST/AnnotatedSeurat/StromalCells.h5seurat')
DefaultAssay(seu.integrated) <- 'SCT'
Idents(seu.integrated) <- seu.integrated$celltype

#Plotting:

DefaultAssay(seu.integrated) <- 'SCT'
DimPlot(seu.integrated, reduction = 'umap', group.by = 'celltype', cols = c('cyan4', 'gold3', 'chartreuse4'))
DotPlot(seu.integrated,
        features = c('PECAM1', 'CDH5', # endothelial
                     'ECM1', 'COL1A1', 'COL1A2', # fibroblast
                     'TAGLN', 'MYH11', 'ACTG2'), # muscle
        col.min = -.5, col.max = 1
) + RotatedAxis() + 
  scale_colour_gradient2(low="khaki1", mid="darkseagreen", high="darkslategray4")

### View session information

sessionInfo()
#R version 4.2.2 Patched (2022-11-10 r83330)
#Platform: x86_64-pc-linux-gnu (64-bit)
#Running under: Ubuntu 22.04.1 LTS

#Matrix products: default
#BLAS:   /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3
#LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/libopenblasp-r0.3.20.so

#locale:
#  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C               LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8     LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8    LC_PAPER=en_US.UTF-8       LC_NAME=C                 
#[9] LC_ADDRESS=C               LC_TELEPHONE=C             LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       

#attached base packages:
#  [1] stats4    stats     graphics  grDevices utils     datasets  methods   base     

#other attached packages:
#  [1] SeuratDisk_0.0.0.9020 ggplot2_3.3.6         sp_1.5-0              SeuratObject_4.1.1    Seurat_4.1.1         

#loaded via a namespace (and not attached):
#  [1] scattermore_0.8             R.methodsS3_1.8.2           pkgmaker_0.32.2             tidyr_1.2.1                 bit64_4.0.5                 knitr_1.40                  irlba_2.3.5                 DelayedArray_0.22.0        
#[9] R.utils_2.12.0              data.table_1.14.2           rpart_4.1.16                KEGGREST_1.36.3             RCurl_1.98-1.8              doParallel_1.0.17           generics_0.1.3              BiocGenerics_0.42.0        
#[17] ScaledMatrix_1.4.1          cowplot_1.1.1               RSQLite_2.2.17              RANN_2.6.1                  future_1.28.0               bit_4.0.4                   phylobase_0.8.10            spatstat.data_2.2-0        
#[25] xml2_1.3.3                  httpuv_1.6.6                SummarizedExperiment_1.26.1 assertthat_0.2.1            viridis_0.6.2               xfun_0.33                   hms_1.1.2                   evaluate_0.16              
#[33] promises_1.2.0.1            fansi_1.0.3                 progress_1.2.2              igraph_1.3.4                DBI_1.1.3                   htmlwidgets_1.5.4           spatstat.geom_2.4-0         purrr_0.3.4                
#[41] ellipsis_0.3.2              dplyr_1.0.10                V8_4.2.2                    annotate_1.74.0             gridBase_0.4-7              locfdr_1.1-8                deldir_1.0-6                sparseMatrixStats_1.8.0    
#[49] MatrixGenerics_1.8.1        vctrs_0.4.1                 SingleCellExperiment_1.18.0 Biobase_2.56.0              ROCR_1.0-11                 abind_1.4-5                 cachem_1.0.6                withr_2.5.0                
#[57] ggforce_0.3.4               progressr_0.11.0            sctransform_0.3.4           prettyunits_1.1.1           goftest_1.2-3               softImpute_1.4-1            cluster_2.1.4               ape_5.6-2                  
#[65] lazyeval_0.2.2              crayon_1.5.1                genefilter_1.78.0           hdf5r_1.3.5                 edgeR_3.38.4                pkgconfig_2.0.3             tweenr_2.0.2                GenomeInfoDb_1.32.4        
#[73] nlme_3.1-159                vipor_0.4.5                 rlang_1.0.5                 globals_0.16.1              lifecycle_1.0.2             miniUI_0.1.1.1              registry_0.5-1              rsvd_1.0.5                 
#[81] dichromat_2.0-0.1           cellranger_1.1.0            polyclip_1.10-0             matrixStats_0.62.0          lmtest_0.9-40               rngtools_1.5.2              Matrix_1.5-1                Rhdf5lib_1.18.2            
#[89] zoo_1.8-10                  beeswarm_0.4.0              ggridges_0.5.3              png_0.1-7                   viridisLite_0.4.1           bitops_1.0-7                R.oo_1.25.0                 rncl_0.8.6                 
#[97] KernSmooth_2.23-20          rhdf5filters_1.8.0          Biostrings_2.64.1           blob_1.2.3                  DelayedMatrixStats_1.18.0   stringr_1.4.1               zinbwave_1.18.0             parallelly_1.32.1          
#[105] spatstat.random_2.2-0       S4Vectors_0.34.0            beachmat_2.12.0             scales_1.2.1                memoise_2.0.1               magrittr_2.0.3              plyr_1.8.7                  ica_1.0-3                  
#[113] howmany_0.3-1               zlibbioc_1.42.0             compiler_4.2.2              dqrng_0.3.0                 RColorBrewer_1.1-3          fitdistrplus_1.1-8          cli_3.4.0                   ade4_1.7-19                
#[121] XVector_0.36.0              listenv_0.8.0               patchwork_1.1.2             pbapply_1.5-0               MASS_7.3-58.1               mgcv_1.8-40                 tidyselect_1.1.2            stringi_1.7.8              
#[129] yaml_2.3.5                  BiocSingular_1.12.0         locfit_1.5-9.6              ggrepel_0.9.1               grid_4.2.2                  tools_4.2.2                 future.apply_1.9.1          parallel_4.2.2             
#[137] rstudioapi_0.14             uuid_1.1-0                  foreach_1.5.2               RNeXML_2.4.7                gridExtra_2.3               farver_2.1.1                Rtsne_0.16                  ggraph_2.0.6               
#[145] digest_0.6.29               rgeos_0.5-9                 shiny_1.7.2                 Rcpp_1.0.9                  GenomicRanges_1.48.0        scuttle_1.6.3               later_1.3.0                 RcppAnnoy_0.0.19           
#[153] AnnotationDbi_1.58.0        httr_1.4.4                  kernlab_0.9-31              colorspace_2.0-3            XML_3.99-0.10               tensor_1.5                  reticulate_1.26             clusterExperiment_2.16.0   
#[161] IRanges_2.30.1              splines_4.2.2               uwot_0.1.14                 spatstat.utils_2.3-1        graphlayouts_0.8.1          mapproj_1.2.8               plotly_4.10.0               xtable_1.8-4               
#[169] jsonlite_1.8.0              tidygraph_1.2.2             R6_2.5.1                    pillar_1.8.1                htmltools_0.5.3             mime_0.12                   NMF_0.24.0                  glue_1.6.2                 
#[177] fastmap_1.1.0               BiocParallel_1.30.3         BiocNeighbors_1.14.0        codetools_0.2-18            maps_3.4.0                  utf8_1.2.2                  lattice_0.20-45             spatstat.sparse_2.1-1      
#[185] tibble_3.1.8                curl_4.3.2                  ggbeeswarm_0.6.0            leiden_0.4.3                gtools_3.9.3                survival_3.4-0              limma_3.52.3                rmarkdown_2.16             
#[193] munsell_0.5.0               rhdf5_2.40.0                GenomeInfoDbData_1.2.8      iterators_1.0.14            HDF5Array_1.24.2            reshape2_1.4.4              gtable_0.3.1                spatstat.core_2.4-4  