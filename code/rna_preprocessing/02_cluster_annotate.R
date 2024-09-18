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
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                                Import data                               ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
seurat <- qread("data/processed/filtered_seurat.qs")

load("data/metadata/cycle.rda") # cell cycle

## split into CD3 and CD45 groups ---------------------------------------------
sc_cd3 <- subset(seurat, flow_cell == "CD3")
sc_cd45 <- subset(seurat, flow_cell == "CD45")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                            Normalize and Cluster                         ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## CD3 clustering --------------------------------------------------------------
sc_cd3 <- sc_cd3 %>% NormalizeData()
sc_cd3 <- sc_cd3 %>% CellCycleScoring(s.features = s_genes, g2m.features = g2m_genes)
sc_cd3 <- sc_cd3 %>% FindVariableFeatures()
sc_cd3 <- sc_cd3 %>% ScaleData()
sc_cd3 <- sc_cd3 %>% RunPCA()

sc_cd3 <- sc_cd3 %>% RunUMAP(reduction = "pca", dims = 1:40)
DimPlot(sc_cd3, group.by = c("patient_group"))
DimPlot(sc_cd3, group.by = c("sample", "Phase"))
DimPlot(sc_cd3, group.by = c("scDbl_class","doubletfind_class"))

FeaturePlot(sc_cd3, min.cutoff = "q1", max.cutoff = "q99", features = c("CD3E","CD3D","CD3G"))
FeaturePlot(sc_cd3, min.cutoff = "q1", max.cutoff = "q99", 
            features = c("nCount_RNA","nFeature_RNA","mitoRatio","log10GenesPerUMI"))

### remove doublets and 