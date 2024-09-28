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


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                               subcluster CD4                             ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
cd4 <- seurat[,seurat$azimuth_celltype_1 == "CD4 T"]
cd4 <- cd4 %>% FindVariableFeatures()
cd4 <- cd4 %>% ScaleData()
cd4 <- cd4 %>% RunPCA()
cd4 <- cd4 %>% RunHarmony(group.by.vars = "sample", dims.use = 1:30)
cd4 <- cd4 %>% RunUMAP(reduction = "harmony", dims = 1:30)
cd4 <- cd4 %>% FindNeighbors(reduction = "harmony", dims = 1:30)
cd4 <- cd4 %>% FindClusters(resolution = c(0.6))


DimPlot(cd4, group.by = c("flow_cell"))
DimPlot(cd4, group.by = c("patient_group"))
DimPlot(cd4, group.by = c("sample", "Phase"))
FeaturePlot(cd4, min.cutoff = "q1", max.cutoff = "q99", 
            features = c("nCount_RNA","nFeature_RNA","mitoRatio","log10GenesPerUMI"))
DimPlot(cd4, group.by = c("RNA_snn_res.0.6", "azimuth_celltype_2", "azimuth_celltype_3"), label = TRUE)

cd4_subtype_markers <- c('IL4I1','IL23R','NR4A1','EGR1','CCL5','ANXA1','S1PR1','KLF2',
                         'CCR7','TOX2','DRAIC','IKZF2','FOXP3','GZMH','GZMK','GBP5',
                         'GPR25','MKI67','TYMS')

DimPlot(cd4, group.by = c("RNA_snn_res.0.6"), label = TRUE)
VlnPlot(cd4, features = cd4_subtype_markers, split.by = "RNA_snn_res.0.6", stack = TRUE, pt.size = 0)


cd4$annotation_level3 <- plyr::revalue(cd4$RNA_snn_res.0.6, 
                                    c('0' = "Th1_Trm_R",
                                      '1' = "Th1_Trm_A",
                                      '2' = "Th1_effector",
                                      '3' = "Th17_Trm",
                                      '4' = "Tfh",
                                      '5' = "Central_memory",
                                      '6' = "Treg",
                                      '7' = "Naive",
                                      '8' = "Naive",
                                      '9' = "Cytotoxic",
                                      '10' = "Cycling",
                                      '11' = "Th1_Trm_R",
                                      '12' = "Tcm_1",
                                      '13' = "Tcm_2",
                                      '14' = "Th1_effector",
                                      '15' = "Treg"))
DimPlot(cd4, group.by = "annotation_level3")
qsave(cd4, "data/processed/annotated_cd4.qs")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                               subcluster CD8                             ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
cd8 <- seurat[,seurat$azimuth_celltype_1 == "CD8 T"]
cd8 <- cd8 %>% FindVariableFeatures()
cd8 <- cd8 %>% ScaleData()
cd8 <- cd8 %>% RunPCA()
cd8 <- cd8 %>% RunHarmony(group.by.vars = "sample", dims.use = 1:30)
cd8 <- cd8 %>% RunUMAP(reduction = "harmony", dims = 1:30)
cd8 <- cd8 %>% FindNeighbors(reduction = "harmony", dims = 1:30)
cd8 <- cd8 %>% FindClusters(resolution = c(0.6))


DimPlot(cd8, group.by = c("flow_cell"))
DimPlot(cd8, group.by = c("patient_group"))
DimPlot(cd8, group.by = c("sample", "Phase"))
FeaturePlot(cd8, min.cutoff = "q1", max.cutoff = "q99", 
            features = c("nCount_RNA","nFeature_RNA","mitoRatio","log10GenesPerUMI"))
DimPlot(cd8, group.by = c("RNA_snn_res.1","RNA_snn_res.0.6"), label = TRUE)
DimPlot(cd8, group.by = c("RNA_snn_res.0.6", "azimuth_celltype_2", "azimuth_celltype_3"), label = TRUE)

cd8_subtype_markers <- c('KIR2DL4','TRDC','ANXA1','IL7R','GZMK','EOMES',
                         'SLC4A10','KLRB1','MIAT','RNF213','CCR7','KLF2',
                         'UCP2','GZMB','TYMS','MKI67')

DimPlot(cd8, group.by = c("RNA_snn_res.0.6"), label = TRUE)
VlnPlot(cd8, features = cd8_subtype_markers, split.by = "RNA_snn_res.0.6", stack = TRUE, pt.size = 0)

cd8$annotation_level3 <- plyr::revalue(cd8$RNA_snn_res.0.6, 
                                    c('0' = "Trm_LP2",
                                      '1' = "Trm_IEL",
                                      '2' = "Trm_LP1",
                                      '3' = "Trm_LP1",
                                      '4' = "Cytotoxic_effector",
                                      '5' = "Cytotoxic_effector",
                                      '6' = "CM_Naive",
                                      '7' = "Cytotoxic_effector",
                                      '8' = "Cycling",
                                      '9' = "Trm_LP1",
                                      '10' = "MAIT",
                                      '11' = "Cycling",
                                      '12' = "Term_effector"))
DimPlot(cd8, group.by = "annotation_level3", label = TRUE)
qsave(cd8, "data/processed/annotated_cd8.qs")



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                        add annotation to larger obj                      ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
metadata <- seurat@meta.data
metadata$annotation_level2 <- NA
metadata$annotation_level2[metadata$barcode %in% cd4$barcode] <- "CD4_Tcell"
metadata$annotation_level2[metadata$barcode %in% cd8$barcode] <- "CD8_Tcell"

cd_meta <- rbind(cd4@meta.data[,c("barcode","annotation_level3")],
                 cd8@meta.data[,c("barcode","annotation_level3")])

metadata <- left_join(metadata, cd_meta, by = "barcode")
rownames(metadata) <- metadata$barcode

seurat@meta.data <- metadata
DimPlot(seurat, group.by = c("annotation_level2","annotation_level3"))

qsave(seurat, "data/processed/annotated_Tcell.qs")
