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


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                                  Annotate                                ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
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
            features = manuscript_markers[[3]])
FeaturePlot(seurat, min.cutoff = "q1", max.cutoff = "q99", 
            features = c("XCL1","XCL2","GNLY","NKG7"))
FeaturePlot(seurat, min.cutoff = "q1", max.cutoff = "q99", 
            features = c("CD3E","CD3D","CD3G"))
FeaturePlot(seurat, min.cutoff = "q1", max.cutoff = "q99", 
            features = c("PTPRC"))
### add annotation ----------------------------------------------------------
seurat$annotation_level1 <- plyr::revalue(seurat$RNA_snn_res.0.4, 
                                    c('0' = "T_cell",
                                      '1' = "T_cell",
                                      '2' = "T_cell",
                                      '3' = "T_cell",
                                      '4' = "T_cell",
                                      '5' = "T_cell",
                                      '6' = "T_cell",
                                      '7' = "Plasma",
                                      '8' = "Plasma",
                                      '9' = "Plasma",
                                      '10' = "B_cell",
                                      '11' = "T_cell",
                                      '12' = "B_cell",
                                      '13' = "Myeloid",
                                      '14' = "T_cell",
                                      '15' = "Mast", 
                                      '16' = "T_cell",
                                      '17' = "T_cell",
                                      '18' = "ILC",
                                      '19' = "B_cell",
                                      '20' = "T_cell",
                                      '21' = "Other",
                                      '22' = "Other"))

DimPlot(seurat, group.by = "annotation_level1", label = TRUE)


