library(Seurat)
library(tidyverse)
library(harmony)
library(gprofiler2)
library(SCPA)
library(decoupleR)
library(EnhancedVolcano)
library(caret)
library(RColorBrewer)

#### Input and processing ####
Szabo <- readRDS("Szabo_Seurat_Not_Processed.rds") # path to your file
Szabo_seurat <- Szabo

for(i in 1:length(Szabo_seurat)){
  Szabo_seurat[[i]] <- Szabo_seurat[[i]] %>%
    Seurat::NormalizeData(verbose = FALSE) %>%
    FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>% 
    ScaleData(verbose = FALSE) %>% 
    RunPCA(npcs = 20, verbose = FALSE) %>%
    RunHarmony("donor", plot_convergence = F) %>% # can skip if you want to not run batch process
    RunUMAP(reduction = "harmony", dims = 1:20) %>% 
    FindNeighbors(reduction = "harmony", dims = 1:20) %>% 
    FindClusters(resolution = 0.7) %>% 
    identity()
}

# RESEARCH QUESTIONS
# 1. Describe the difference between TCR-stimulated and unstimulated samples at the mRNA (gene) level.
# 2. What transcription factors are regulating TCR activation in the dataset? 
# 3. Is TCR stimulation different in different tissues, or does it give rise to the same transcriptional profile?
# 4. Develop a machine learning model to predict TCR and unstimulated cells.


#### TCR vs Resting ####
DE_markers <- Szabo_seurat
for(i in 1:length(Szabo_seurat)){
  Idents(Szabo_seurat[[i]]) <- "State"
  DE_markers[[i]] <- FindMarkers(Szabo_seurat[[i]], ident.1 = "Activated", ident.2 = "Resting")
}

DE_markers_signif <- DE_markers
for(i in 1:length(DE_markers_signif)){
  DE_markers_signif[[i]] <- DE_markers_signif[[i]] %>%
    dplyr::filter(p_val_adj < 0.05 & abs(avg_log2FC) > log2(1.5))
}
ora_result <- DE_markers_signif
for( i in 1:length(DE_markers_signif)){
  ora_result[[i]] <- gost(query = unique( rownames(DE_markers_signif[[i]])),
                          organism = "hsapiens", sources = "KEGG" )
  ora_result[[i]] <- ora_result[[i]]$result
}
FeaturePlot(Szabo_seurat$Blood, features = "IL2RA")

#### Difference of TCR in different tissues? ####
DE_markers_signif_vec <- DE_markers_signif
for(i in 1:length(DE_markers_signif_vec)){
  DE_markers_signif_vec[[i]]$gene <- rownames(DE_markers_signif_vec[[i]])
  DE_markers_signif_vec[[i]] <- DE_markers_signif_vec[[i]] %>%
    dplyr::pull(gene)
}

library(UpSetR)
input_data <- fromList(DE_markers_signif_vec)
upset(input_data, 
      order.by = "freq", 
      main.bar.color = "steelblue", 
      sets.bar.color = "darkgray",
      matrix.color = "black",
      point.size = 3.5, 
      line.size = 1)

#### TF activities ####
load("TFN.RData") # load the TF-targets database

for(i in 1:length(Szabo_seurat)){
  temp <- Szabo_seurat[[i]]@assays$RNA$counts
  sample_acts <- decoupleR::run_ulm(mat = temp,
                                    net = TFN,
                                    .source = 'source',
                                    .target = 'target',
                                    .mor = 'mor',
                                    minsize = 5)

  Szabo_seurat[[i]][['tfsulm']] <- sample_acts %>%
    tidyr::pivot_wider(id_cols = 'source',
                       names_from = 'condition',
                       values_from = 'score') %>%
    tibble::column_to_rownames('source') %>%
    Seurat::CreateAssayObject(.)

  # Change assay
  DefaultAssay(object = Szabo_seurat[[i]]) <- "tfsulm"
}
saveRDS(Szabo_seurat, "szabo_seurat_ULM.rds")
Szabo_seurat_ULM <- readRDS("szabo_seurat_ULM.rds")
p1 <- Seurat::FeaturePlot(Szabo_seurat_ULM$Blood, features = "EGR1") 
p2 <-  DimPlot(Szabo_seurat_ULM$Blood, group.by = "State")
cowplot::plot_grid(p1, p2)

#### ML model ####
Seurat_Blood <- Szabo_seurat$Blood
DefaultAssay(object = Seurat_Blood ) <- "RNA"

features <- VariableFeatures(Seurat_Blood)
data_matrix <- as.data.frame(t(as.matrix(GetAssayData(Seurat_Blood, layer = "data")[features, ])))

data_matrix$State <- Seurat_Blood$State 
data_matrix$State <- as.factor(data_matrix$State)

# Split into Training (70%) and Testing (30%)
set.seed(111111)
trainIndex <- createDataPartition(data_matrix$State, p = 0.7, list = FALSE)
train_set <- data_matrix[trainIndex, ]
test_set  <- data_matrix[-trainIndex, ]

# 5 cv
train_control <- trainControl(method = "cv", number = 5)
model <- train(State ~ .,
               data = train_set,
               method = "ranger",
               trControl = train_control,
               importance = 'impurity')
saveRDS(model, "RF_Model_Gene_Prediction.rds")
model <- readRDS("RF_Model_Gene_Prediction.rds")

predictions <- predict(model, newdata = test_set)

# Check accuracy
conf_matrix <- confusionMatrix(predictions, test_set$State)
print(conf_matrix)
importance <- varImp(model)

# top 30 features
plot(importance, top = 30)
