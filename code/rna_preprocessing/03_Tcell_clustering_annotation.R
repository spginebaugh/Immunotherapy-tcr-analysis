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

library(SeuratData)
library(Azimuth)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                                Import data                               ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
sc_all <- qread("data/processed/annotated_seurat.qs")

seurat <- sc_all[,sc_all$annotation_level1 == "T_cell"]


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                                  Recluster                               ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
seurat <- seurat %>% FindVariableFeatures(nfeatures = 2200)
exclude_genes <- c(grep("TR[A-Z]V",rownames(seurat), value = TRUE),
                   grep("IG[A-Z]V",rownames(seurat), value = TRUE))
var_feat <- VariableFeatures(seurat)
var_feat <- var_feat[!(var_feat %in% exclude_genes)]
c("CD4","CD8A","CD8B") %in% var_feat
var_feat <- c(var_feat,"CD4") #CD4 was not considered highly variable, so add it in
VariableFeatures(seurat) <- var_feat

seurat <- seurat %>% ScaleData(vars.to.regress = c("S.Score","G2M.Score","nCount_RNA"))
seurat <- seurat %>% RunPCA()
seurat <- seurat %>% RunUMAP(reduction = "pca", dims = 1:30)


seurat <- seurat %>% RunHarmony(group.by.vars = "sample", dims.use = 1:30)
seurat <- seurat %>% RunUMAP(reduction = "harmony", dims = 1:30)

seurat <- seurat %>% FindNeighbors(reduction = "harmony", dims = 1:30)
seurat <- seurat %>% FindClusters(resolution = c(1))

## explicit labeling of subtype based on CD4 and CD8A expression ----------------
cd_expr <- GetAssayData(seurat, layer = "RNA", slot = "counts")[c("CD8A","CD4"),] %>% 
  data.matrix() %>%
  t()
cd8_cell <- rownames(cd_expr)[(cd_expr[,"CD8A"] > 1) & (cd_expr[,"CD4"] < 1)]
cd4_cell <- rownames(cd_expr)[(cd_expr[,"CD8A"] < 1) & (cd_expr[,"CD4"] > 1)]

seurat$CD8_expr <- FALSE
seurat$CD8_expr[colnames(seurat) %in% cd8_cell] <- TRUE

seurat$CD4_expr <- FALSE
seurat$CD4_expr[colnames(seurat) %in% cd4_cell] <- TRUE

DimPlot(seurat, group.by = c("flow_cell"))
DimPlot(seurat, group.by = c("patient_group"))
DimPlot(seurat, group.by = c("sample", "Phase"))
DimPlot(seurat, group.by = c("CD8_expr", "CD4_expr"))
FeaturePlot(seurat, min.cutoff = "q1", max.cutoff = "q99", 
            features = c("nCount_RNA","nFeature_RNA","mitoRatio","log10GenesPerUMI"))
FeaturePlot(seurat, min.cutoff = "q1", max.cutoff = "q99", 
            features = c("CD4","CD8A","CD8B"))
DimPlot(seurat, group.by = c("RNA_snn_res.1"), label = TRUE)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                             Azimuth annotation                           ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
sc_azimuth <- seurat
sc_azimuth <- RunAzimuth(sc_azimuth, reference = "pbmcref")

DimPlot(sc_azimuth, group.by = c("predicted.celltype.l1","predicted.celltype.l2","predicted.celltype.l3"))

seurat$azimuth_celltype_1 <- sc_azimuth$predicted.celltype.l1
seurat$azimuth_celltype_2 <- sc_azimuth$predicted.celltype.l2
seurat$azimuth_celltype_3 <- sc_azimuth$predicted.celltype.l3

## will use azimuth celltype level1 to select CD8 and CD4

qsave(seurat, "data/processed/clustered_tcell.qs")
