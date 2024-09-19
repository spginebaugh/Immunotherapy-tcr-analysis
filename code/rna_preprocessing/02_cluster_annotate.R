# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                              Load Libraries                              ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
library(Seurat)
library(qs)
library(ggplot2)
library(magrittr)
library(dplyr)
library(harmony)
library(stringr)
set.seed(42555)

source("code/utils/rna_utils.R")
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                                Import data                               ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
seurat <- qread("data/processed/filtered_seurat.qs")

load("/media/largedata/universal_files/cycle.rda") # cell cycle

## split into CD3 and CD45 groups ---------------------------------------------
# sc_cd3 <- subset(seurat, flow_cell == "CD3")
# sc_cd45 <- subset(seurat, flow_cell == "CD45")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                            Normalize and Cluster                         ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
seurat <- seurat %>% NormalizeData()
seurat <- seurat %>% CellCycleScoring(s.features = s_genes, g2m.features = g2m_genes)
seurat <- seurat %>% FindVariableFeatures()
seurat <- seurat %>% ScaleData()
seurat <- seurat %>% RunPCA()
seurat <- seurat %>% RunUMAP(reduction = "pca", dims = 1:40)

seurat <- seurat %>% FindNeighbors(reduction = "pca", dims = 1:40)
seurat <- seurat %>% FindClusters(resolution = c(0.4))

### basic plotting ------------------------------------------------------------
DimPlot(seurat, group.by = c("flow_cell"))
DimPlot(seurat, group.by = c("patient_group"))
DimPlot(seurat, group.by = c("sample", "Phase"))
FeaturePlot(seurat, min.cutoff = "q1", max.cutoff = "q99", 
            features = c("nCount_RNA","nFeature_RNA","mitoRatio","log10GenesPerUMI"))

### broad celltype marker plotting -------------------------------------------
DimPlot(seurat, group.by = c("RNA_snn_res.0.4"), label = TRUE)

panglao_scores <- score_panglao(seurat)

manuscript_markers <- list()
manuscript_markers[["plasma"]] <- c("IGHA1","IGHG1")
manuscript_markers[["mast"]] <- c("TPSAB1")
manuscript_markers[["myeloid"]] <- c("S100A9","C1QC")
manuscript_markers[["B_cell"]] <- c("MS4A1","FCER2","MEF2B")
manuscript_markers[["ILC"]] <- c("CSF2","IL7R","GZMB")

FeaturePlot(seurat, min.cutoff = "q1", max.cutoff = "q99", 
            features = manuscript_markers[["ILC"]])

seurat$annotation_level1 <- revalue(seurat$SCT_snn_res.1, 
                                    c('0' = "",
                                      '1' = "",
                                      '2' = "",
                                      '3' = "",
                                      '4' = "",
                                      '5' = "",
                                      '6' = "",
                                      '7' = "",
                                      '8' = "",
                                      '9' = "",
                                      '10' = "",
                                      '11' = "",
                                      '12' = "",
                                      '13' = "",
                                      '14' = "",
                                      '15' = "", 
                                      '16' = "",
                                      '17' = "",
                                      '18' = "",
                                      '19' = "",
                                      '20' = "",
                                      '21' = "",
                                      '22' = ""



