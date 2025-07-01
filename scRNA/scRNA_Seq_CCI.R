##########################PHYSICAL CELL-CELL INTERACTIONS###########################
############################ Load Functions#################################################
####generating gene signature#########
library(Seurat)
library(tidyverse)
library(decoupleR)
GetSignature = function(seurat_obj, ident_col = NULL, n = 100, p_val = 0.05){
  if (class(seurat_obj) != 'Seurat') {
    stop('Input data is not a Seurat object, please provide a valid Seurat object')
  }
  if(!is.null(ident_col)){
    print(paste('using the specified seurat ident to generate signatures'))
    Idents(seurat_obj) <- ident_col
  }
  else{
    print('using the default seurat ident to generate signatures')
  }
  Markers <- FindAllMarkers(seurat_obj, only.pos = T)
  Markers <- Markers %>% filter(p_val_adj <0.05) %>% group_by(cluster) %>% 
    slice_max(., order_by = avg_log2FC, n = n) %>% dplyr:: select(c(cluster, gene)) %>%
    mutate(mor= 1, cluster= as.character(cluster))
  colnames(Markers) <- c('source', 'target', 'mor')
  return(Markers)
  
}

#######getting cell scores###############
GetCellScores <- function(seurat_obj, signatures,  assay = 'RNA', slot = 'counts'){
  mat <- GetAssayData(seurat_obj, assay = assay, slot = slot )
  acts <- run_ulm(mat=mat, net= signatures, .source='source', .target='target',
                  .mor='mor', minsize = 5)
  colnames(acts)[2:3] <- c('celltype', 'barcode')
  acts <- acts[, c('barcode', 'celltype', 'score', 'p_value', 'statistic')]
  return(acts)
}

#####getting final cell assignments############
GetCellAssignments <- function(score_data, p_val = 0.05, cut_off = 1){
  acts_filt <- score_data %>% filter(p_value <= p_val & score > cut_off) %>% 
    mutate(count_ulm = 1,
           celltype = str_replace_all(celltype, '_', ' '))
  cell_class <- acts_filt %>% 
    group_by(barcode) %>% 
    mutate(count_ulm= sum(count_ulm),
           celltype_ulm= paste(celltype, collapse = '_'),
           avg_pvalue = mean(p_value),
           avg_score = mean(score)) %>% 
    dplyr::select(-c(celltype, p_value, score)) %>% distinct()
  
  return(cell_class)
}

######adding cell assignments to seurat object########

AddMetaObject <- function(seurat_obj, cell_class_df) {
  new_meta <- seurat_obj@meta.data %>% rownames_to_column('barcode_ulm') %>% 
    left_join(cell_class_df, by= c('barcode_ulm'= 'barcode')) %>% column_to_rownames('barcode_ulm')
  seurat_obj@meta.data <- new_meta
  return(seurat_obj)
}

#####obtaining multiplet data########

GetMultiplet <- function(seurat_obj, minCells = 2) {
  multObj <- subset(seurat_obj, subset = count_ulm >= minCells)
  multSummary <- as.data.frame(table(multObj$celltype_ulm))
  colnames(multSummary) <- c('multiplet type', 'frequency')
  return(list(multSummary =multSummary, multObj = multObj))
  
}

####filtering multiplets########

FilterMultiplet <- function(seurat_obj, minCells = 2, minFreq = 10){
  multObj <- subset(seurat_obj, subset = count_ulm >= minCells)
  multSummary <- as.data.frame(table(multObj$celltype_ulm))
  colnames(multSummary) <- c('multiplet type', 'number')
  multSummary <- multSummary[multSummary$number >= minFreq,]
  Idents(multObj) <- multObj$celltype_ulm
  celltypes <- unique(multSummary$`multiplet type`)
  multObj <- subset(multObj, idents = celltypes)
  return(list(multSummaryFilt = multSummary, multObjFilt = multObj))
}

######getting nodes and edges##########
GetNodeDF <- function(mat){
  mat <- as.data.frame(mat)
  my_network <- apply(mat, 1, function (x){
    vec_split <- str_split(x[1], '_', simplify =F)
    vec_df <- as.data.frame(vec_split)
    colnames(vec_df)[1] <- 'cells'
    vec_comb <- t(combn(vec_df$cells, 2)) %>% 
      as.data.frame() %>% 
      mutate(n_cells = as.numeric(x[2])
      )
  })
  network_df <- do.call(rbind, my_network)
  colnames(network_df) <- c('Cell1', 'Cell2', 'n_cells')
  network_df <- aggregate(n_cells~Cell1+Cell2, data= network_df, FUN = sum)
  return(network_df)
}

#########plotting network###############
library(igraph)
library(ggraph)
library(tidyverse)
library(tidygraph)

PlotNetwork <- function(network_df, node_size = 20, node_color = 'blue', 
                        node_text_size = 4, edge_width_factor = 50, 
                        edge_color = 'red', network_layout = 'fr', 
                        legend_title = 'scaled counts', legend_position = 'bottom',
                        min_edge_width = 0.5, max_edge_width = 3, 
                        main = 'Network plot', main_size = 15,
                        hjust = 0.5, legend_text_size = 12, legend_title_size = 14) {
  g <- graph_from_data_frame(network_df, directed = FALSE)
  
  tg <- as_tbl_graph(g)
  
  # Set node and edge aesthetics
  tg <- tg %>%
    activate(nodes) %>%
    mutate(node_size = node_size) %>%
    activate(edges) %>%
    mutate(edge_width = n_cells / edge_width_factor, 
           edge_color = edge_color)
  
  # Create the network plot using ggraph
  net_plot <- ggraph(tg, layout = network_layout) +  # You can use other layouts like "kk", "lgl", etc.
    geom_edge_link(aes(width = edge_width), color = edge_color, show.legend = TRUE) +
    geom_node_point(aes(size = node_size), color = node_color, show.legend = FALSE) +
    geom_node_text(aes(label = name), vjust = 1.8, size = node_text_size, show.legend = FALSE) +
    scale_edge_width_continuous(range = c(min_edge_width, max_edge_width), name = legend_title) +
    theme_void() +
    ggtitle(main) +
    theme(
      legend.position = legend_position,
      legend.text = element_text(size = legend_text_size),  
      legend.title = element_text(size = legend_title_size),
      plot.title = element_text(size = main_size, hjust = hjust)
    ) +
    guides(size = "none")
  print(net_plot)
}



########################INFER CCC NETWORK#####################################
seuData <- readRDS("/mnt/8TB/users/shameed/shameed/Olbretch/2095-Olbrecht_Seurat_Celltype.rds")
head(seuData@meta.data)
table(seuData$celltype_cluster)

################generate signature
set.seed(1610242)
OC_sig <- GetSignature(seuData, ident_col = seuData$celltype_cluster)
saveRDS(OC_sig, 'OC_sig.rds')
###########score cells
my_scores <- GetCellScores(seurat_obj = seuData, signatures = OC_sig, assay = 'RNA', slot = 'data')
#########assign cell labels
my_ass <- GetCellAssignments(score_data = my_scores, cut_off = 1)
############add labels to the metadata 
seuData <- AddMetaObject(seuData, cell_class_df = my_ass)
colnames(seuData@meta.data)
################subset multiplets 
my_mult <- GetMultiplet(seuData)
####################filter to retain frequent multiplets
my_mult_filt <- FilterMultiplet(seuData, minFreq = 5)
###############plot network
my_network_df <- GetNodeDF(mat = my_mult_filt$multSummaryFilt)
p1<-PlotNetwork(network_df = my_network_df, edge_width_factor = 10, node_text_size = 5,
                legend_text_size = 15, legend_title_size = 15, main_size =20 )


########################DEG AND PATHWAY ANALYSIS (T cells vs Macrophages)###############################
table(seuData$celltype_ulm)
Idents(seuData) <- seuData$celltype_ulm
deg <- FindMarkers(seuData, ident.1 = 'T cells', ident.2 = "Macrophage",min.pct = 0.25)
EnhancedVolcano::EnhancedVolcano(deg, lab = rownames(deg), 
                                 x='avg_log2FC',y='p_val_adj',
                                 title ='Differential genes',
                                 pCutoff = 0.05,
                                 FCcutoff = 1.5,
                                 #pointSize = 3.0,
                                 #labSize = 6.0,
                                 #col=c('black', 'black', 'black', 'red3'),
                                 colAlpha = 1)

#######################Enriched pathways
##get significant genes and order by fold change
gsea_gene<- deg %>% rownames_to_column('Gene') %>% 
  filter(p_val <=0.05) %>% dplyr::select(c(Gene, avg_log2FC))
gsea_gene<- gsea_gene[order(-gsea_gene$avg_log2FC),] 

###make a named list
gene_list<- gsea_gene$avg_log2FC
#gene_list<- jitter(gsea_gene$avg_log2FC, factor = 0.01) ##if there is error due to many matches in avg_log2FC
names(gene_list)<- gsea_gene$Gene

#######obtain enriched pathways
gsea_path<- fgseaSimple(pathways = H.symbol, stats = gene_list, nperm = 1000)
gsea_path$pathway<- str_replace(gsea_path$pathway, 'REACTOME_', '')
gsea_path <- gsea_path%>% filter(pval <0.05)

######visualise pathways
gsea_path %>% slice_max(., n = 50, order_by = abs(NES)) %>% 
  ggplot(., aes(x= reorder(pathway, NES), y= NES, fill= pval)) + 
  geom_col() + coord_flip() + labs(y='normalised enrichment scores',
                                   x= 'Reactome pathways') + 
  ggtitle('Enriched pathways (T cells vs Macrophage)')
