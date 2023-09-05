library(tidyverse)
library(biomaRt)  # used to map GOterms to Ensembl IDs
library(topGO)
library(Seurat)
library(SeuratDisk)
library(ggplot2)
library(readxl)
library(ggnewscale)

seu <- readRDS('/home/Jayne.Wiarda/SI_PP_SC_ST/AnnotatedSeurat/AllSamples_annotated_ST.rds')

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
    
    all_GO <- read_tsv('/home/Jayne.Wiarda/SI_PP_SC_ST/GeneOntology/AllSpatialDots_gene_to_GO.tsv')
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

# spatial transcriptomics regions ----

# to map gene 'FinalList' annotation to Ensembl gene_IDs
gene_IDs <- read_excel("/home/Jayne.Wiarda/SI_PP_SC_ST/GeneAnnotationFiles/UpdatedGeneNameListForSus97GTF_06302021_JEW_SKS.xlsx")
gene_IDs$EnsemblID <- gene_IDs$ENSID # make identical to column name in EG2GO

### Manual annotation

background <- read_excel('/home/Jayne.Wiarda/SI_PP_SC_ST/DGE/SpatialRegions_DEGs_BackgroundGeneList_AllSamples_ManualAnnotation.xlsx')
detected_genes <- background$FinalAnnot

# all GO terms detected in this dataset

GO_gene_universe <- 
  gene_IDs %>% 
  filter(FinalList %in% detected_genes) %>% 
  left_join(EG2GO) %>% 
  filter(!is.na(GO))

GO_gene_universe %>%
  dplyr::select(FinalList, GO) %>% 
  write_tsv('/home/Jayne.Wiarda/SI_PP_SC_ST/GeneOntology/AllSpatialDots_gene_to_GO.tsv')

SpatialManDGE <- read_xlsx('/home/Jayne.Wiarda/SI_PP_SC_ST/DGE/SpatialRegions_DEGs_AllSamples_ManualAnnotation.xlsx')

SpatialMan_results <- 
  SpatialManDGE %>% 
  filter(avg_log2FC > 0) %>% 
  group_by(cluster) %>% 
  nest() %>% 
  mutate(enriched_genes=map(data, ~.x %>% pull(gene))) %>% 
  mutate(GO_results=
           map(enriched_genes, ~topGO_wrapper(myInterestingGenes = .x,
                                              mapping_file = '/home/Jayne.Wiarda/SI_PP_SC_ST/GeneOntology/AllSpatialDots_gene_to_GO.tsv'))) %>% 
  mutate(filt_results=purrr::map(.x = GO_results, .f = ~filter(.x, pval < 0.05))) %>% 
  dplyr::select(cluster, filt_results) %>% 
  unnest(cols=filt_results)

SpatialMan_results$Fold_enrichment <- SpatialMan_results$Significant/SpatialMan_results$Expected # add in stat for fold enrichment
SpatialMan_results <- subset(SpatialMan_results, Significant > 1) # eliminate terms with only 1 gene present

SpatialMan_results %>% write_tsv('/home/Jayne.Wiarda/SI_PP_SC_ST/GeneOntology/AllSpatialDots_ManualAnnotation_GOresults.tsv')

# Create dot plot of selected GO processes for spatial regions
terms <- c('intestinal epithelial cell differentiation',
           'transepithelial transport',
           'carbohydrate transport',
           'fatty acid transport',
           'lipid transport',
           'positive regulation of interferon-gamma production',
           'positive regulation of tumor necrosis factor production',
           'leukocyte mediated cytotoxicity',
           'T cell costimulation',
           'T cell differentiation',
           'positive regulation of regulatory T cell differentiation',
           'positive regulation of alpha-beta T cell proliferation',
           'immunological synapse formation',
           'B cell receptor signaling pathway',
           'positive regulation of B cell activation',
           'mature B cell differentiation involved in immune response',
           'B cell differentiation',
           'B cell apoptotic process',
           'immunoglobulin production',
           'isotype switching',
           'cell cycle process',
           'DNA damage checkpoint signaling',
           'cellular response to DNA damage stimulus',
           'collagen biosynthetic process',
           'basement membrane assembly',
           'basement membrane organization')
SpatialMan_results <- SpatialMan_results[SpatialMan_results$Term %in% terms, ]
SpatialMan_results$cluster <- factor(SpatialMan_results$cluster, 
                                     levels=c(levels(seu$Region_Clust)))
SpatialMan_results$Term <- factor(SpatialMan_results$Term,
                          levels = c('intestinal epithelial cell differentiation',
                                     'transepithelial transport',
                                     'carbohydrate transport',
                                     'fatty acid transport',
                                     'lipid transport',
                                     'positive regulation of interferon-gamma production',
                                     'positive regulation of tumor necrosis factor production',
                                     'leukocyte mediated cytotoxicity',
                                     'T cell costimulation',
                                     'T cell differentiation',
                                     'positive regulation of regulatory T cell differentiation',
                                     'positive regulation of alpha-beta T cell proliferation',
                                     'immunological synapse formation',
                                     'B cell receptor signaling pathway',
                                     'positive regulation of B cell activation',
                                     'mature B cell differentiation involved in immune response',
                                     'B cell differentiation',
                                     'B cell apoptotic process',
                                     'immunoglobulin production',
                                     'isotype switching',
                                     'cell cycle process',
                                     'DNA damage checkpoint signaling',
                                     'cellular response to DNA damage stimulus',
                                     'collagen biosynthetic process',
                                     'basement membrane assembly',
                                     'basement membrane organization'))
SpatialMan_results$pval <- as.numeric(SpatialMan_results$pval)
ggplot(SpatialMan_results) +
  geom_point(aes(x = cluster, 
                 y = Term,
                 size = Fold_enrichment,
                 color = pval)) +
  theme_bw() + 
  scale_colour_gradient2(low="darkslateblue", mid="darkseagreen", high="khaki3", 
                         limits=c(0, 0.05),
                         midpoint = .025) +
  theme(axis.text.x = element_text(angle=90, vjust=.5, hjust=1))

### Clustering annotation

background <- read_excel('/home/Jayne.Wiarda/SI_PP_SC_ST/DGE/SpatialRegions_DEGs_BackgroundGeneList_AllSamples_ClusteringAnnotation.xlsx')
detected_genes <- background$FinalAnnot

# all GO terms detected in this dataset

GO_gene_universe <- 
  gene_IDs %>% 
  filter(FinalList %in% detected_genes) %>% 
  left_join(EG2GO) %>% 
  filter(!is.na(GO))

GO_gene_universe %>%
  dplyr::select(FinalList, GO) %>% 
  write_tsv('/home/Jayne.Wiarda/SI_PP_SC_ST/GeneOntology/AllSpatialDots_gene_to_GO.tsv')

SpatialClusDGE <- read_xlsx('/home/Jayne.Wiarda/SI_PP_SC_ST/DGE/SpatialRegions_DEGs_AllSamples_ClusteringAnnotation.xlsx')

SpatialClus_results <- 
  SpatialClusDGE %>% 
  filter(avg_log2FC > 0) %>% 
  group_by(cluster) %>% 
  nest() %>% 
  mutate(enriched_genes=map(data, ~.x %>% pull(gene))) %>% 
  mutate(GO_results=
           map(enriched_genes, ~topGO_wrapper(myInterestingGenes = .x,
                                              mapping_file = '/home/Jayne.Wiarda/SI_PP_SC_ST/GeneOntology/AllSpatialDots_gene_to_GO.tsv'))) %>% 
  mutate(filt_results=purrr::map(.x = GO_results, .f = ~filter(.x, pval < 0.05))) %>% 
  dplyr::select(cluster, filt_results) %>% 
  unnest(cols=filt_results)

SpatialClus_results$Fold_enrichment <- SpatialClus_results$Significant/SpatialClus_results$Expected # add in stat for fold enrichment
SpatialClus_results <- subset(SpatialClus_results, Significant > 1) # eliminate terms with only 1 gene present

SpatialClus_results %>% write_tsv('/home/Jayne.Wiarda/SI_PP_SC_ST/GeneOntology/AllSpatialDots_ClusteringAnnotation_GOresults.tsv')

# Create dot plot of selected GO processes for spatial regions
terms <- c('intestinal epithelial cell differentiation',
           'transepithelial transport',
           'carbohydrate transport',
           'fatty acid transport',
           'lipid transport',
           'positive regulation of interferon-gamma production',
           'positive regulation of tumor necrosis factor production',
           'leukocyte mediated cytotoxicity',
           'T cell costimulation',
           'T cell differentiation',
           'positive regulation of regulatory T cell differentiation',
           'positive regulation of alpha-beta T cell proliferation',
           'immunological synapse formation',
           'B cell receptor signaling pathway',
           'positive regulation of B cell activation',
           'mature B cell differentiation involved in immune response',
           'B cell differentiation',
           'B cell apoptotic process',
           'immunoglobulin production',
           'isotype switching',
           'cell cycle process',
           'DNA damage checkpoint signaling',
           'cellular response to DNA damage stimulus',
           'collagen biosynthetic process',
           'basement membrane assembly',
           'basement membrane organization')
SpatialClus_results <- SpatialClus_results[SpatialClus_results$Term %in% terms, ]
SpatialClus_results$cluster <- factor(SpatialClus_results$cluster, 
                                     levels=c(levels(seu$Region_Clust)))
SpatialClus_results$Term <- factor(SpatialClus_results$Term,
                                  levels = c('intestinal epithelial cell differentiation',
                                             'transepithelial transport',
                                             'carbohydrate transport',
                                             'fatty acid transport',
                                             'lipid transport',
                                             'positive regulation of interferon-gamma production',
                                             'positive regulation of tumor necrosis factor production',
                                             'leukocyte mediated cytotoxicity',
                                             'T cell costimulation',
                                             'T cell differentiation',
                                             'positive regulation of regulatory T cell differentiation',
                                             'positive regulation of alpha-beta T cell proliferation',
                                             'immunological synapse formation',
                                             'B cell receptor signaling pathway',
                                             'positive regulation of B cell activation',
                                             'mature B cell differentiation involved in immune response',
                                             'B cell differentiation',
                                             'B cell apoptotic process',
                                             'immunoglobulin production',
                                             'isotype switching',
                                             'cell cycle process',
                                             'DNA damage checkpoint signaling',
                                             'cellular response to DNA damage stimulus',
                                             'collagen biosynthetic process',
                                             'basement membrane assembly',
                                             'basement membrane organization'))
SpatialClus_results$pval <- as.numeric(SpatialClus_results$pval)
ggplot(SpatialClus_results) +
  geom_point(aes(x = cluster, 
                 y = Term,
                 size = Fold_enrichment,
                 color = pval)) +
  theme_bw() + 
  scale_colour_gradient2(low="darkslateblue", mid="darkseagreen", high="khaki3", 
                         limits=c(0, 0.05),
                         midpoint = .025) +
  theme(axis.text.x = element_text(angle=90, vjust=.5, hjust=1))

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
#  [1] grid      stats4    stats     graphics  grDevices utils     datasets  methods   base     

#other attached packages:
#  [1] ggnewscale_0.4.8      readxl_1.4.1          SeuratDisk_0.0.0.9020 sp_1.5-0              SeuratObject_4.1.1    Seurat_4.1.1          topGO_2.48.0          SparseM_1.81          GO.db_3.15.0          AnnotationDbi_1.58.0  IRanges_2.30.1       
#[12] S4Vectors_0.34.0      Biobase_2.56.0        graph_1.74.0          BiocGenerics_0.42.0   biomaRt_2.52.0        forcats_0.5.2         stringr_1.4.1         dplyr_1.0.10          purrr_0.3.4           readr_2.1.2           tidyr_1.2.1          
#[23] tibble_3.1.8          ggplot2_3.3.6         tidyverse_1.3.2      

#loaded via a namespace (and not attached):
#  [1] utf8_1.2.2             reticulate_1.26        tidyselect_1.1.2       RSQLite_2.2.17         htmlwidgets_1.5.4      Rtsne_0.16             munsell_0.5.0          codetools_0.2-18       ica_1.0-3              future_1.28.0         
#[11] miniUI_0.1.1.1         withr_2.5.0            spatstat.random_2.2-0  colorspace_2.0-3       progressr_0.11.0       filelock_1.0.2         knitr_1.40             rstudioapi_0.14        ROCR_1.0-11            tensor_1.5            
#[21] listenv_0.8.0          labeling_0.4.2         GenomeInfoDbData_1.2.8 polyclip_1.10-0        farver_2.1.1           bit64_4.0.5            parallelly_1.32.1      vctrs_0.4.1            generics_0.1.3         xfun_0.33             
#[31] timechange_0.2.0       BiocFileCache_2.4.0    R6_2.5.1               GenomeInfoDb_1.32.4    hdf5r_1.3.5            bitops_1.0-7           spatstat.utils_2.3-1   cachem_1.0.6           assertthat_0.2.1       vroom_1.5.7           
#[41] promises_1.2.0.1       scales_1.2.1           googlesheets4_1.0.1    rgeos_0.5-9            gtable_0.3.1           globals_0.16.1         goftest_1.2-3          rlang_1.0.6            splines_4.2.2          lazyeval_0.2.2        
#[51] gargle_1.2.1           spatstat.geom_2.4-0    broom_1.0.1            yaml_2.3.5             reshape2_1.4.4         abind_1.4-5            modelr_0.1.10          backports_1.4.1        httpuv_1.6.6           tools_4.2.2           
#[61] ellipsis_0.3.2         spatstat.core_2.4-4    RColorBrewer_1.1-3     ggridges_0.5.3         Rcpp_1.0.9             plyr_1.8.7             progress_1.2.2         zlibbioc_1.42.0        RCurl_1.98-1.8         prettyunits_1.1.1     
#[71] rpart_4.1.16           deldir_1.0-6           pbapply_1.5-0          cowplot_1.1.1          zoo_1.8-10             haven_2.5.1            ggrepel_0.9.1          cluster_2.1.4          fs_1.5.2               magrittr_2.0.3        
#[81] data.table_1.14.2      scattermore_0.8        lmtest_0.9-40          reprex_2.0.2           RANN_2.6.1             googledrive_2.0.0      fitdistrplus_1.1-8     matrixStats_0.62.0     hms_1.1.2              patchwork_1.1.2       
#[91] mime_0.12              evaluate_0.16          xtable_1.8-4           XML_3.99-0.10          gridExtra_2.3          compiler_4.2.2         KernSmooth_2.23-20     crayon_1.5.1           htmltools_0.5.3        mgcv_1.8-40           
#[101] later_1.3.0            tzdb_0.3.0             lubridate_1.9.0        DBI_1.1.3              dbplyr_2.2.1           MASS_7.3-58.1          rappdirs_0.3.3         Matrix_1.5-1           cli_3.4.0              parallel_4.2.2        
#[111] igraph_1.3.4           pkgconfig_2.0.3        plotly_4.10.0          spatstat.sparse_2.1-1  xml2_1.3.3             XVector_0.36.0         rvest_1.0.3            digest_0.6.29          sctransform_0.3.4      RcppAnnoy_0.0.19      
#[121] spatstat.data_2.2-0    Biostrings_2.64.1      rmarkdown_2.16         cellranger_1.1.0       leiden_0.4.3           uwot_0.1.14            curl_4.3.2             shiny_1.7.2            lifecycle_1.0.2        nlme_3.1-159          
#[131] jsonlite_1.8.0         viridisLite_0.4.1      fansi_1.0.3            pillar_1.8.1           lattice_0.20-45        KEGGREST_1.36.3        fastmap_1.1.0          httr_1.4.4             survival_3.4-0         glue_1.6.2            
#[141] png_0.1-7              bit_4.0.4              stringi_1.7.8          blob_1.2.3             memoise_2.0.1          irlba_2.3.5            future.apply_1.9.1  