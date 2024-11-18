# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                              Load Libraries                              ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
library(Seurat)
library(qs)
library(magrittr)
library(dplyr)
library(harmony)
library(stringr)
library(tidyverse)
library(tibble)
set.seed(42555)

library(scRepertoire)
library(ggraph)

library(ggplot2)
library(ggpmisc)
library(ggprism)
library(ggrepel)
library(knitr)
library(ggpubr)
library(cowplot)
library(patchwork)

library(clusterProfiler)
library(msigdbr)
library(circlize)
library(scales)

source("code/utils/TCR_utils.R")


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                                  Functions                               ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
prep_data <- function(seurat_obj) {
  seurat_obj <- add_contig_barcodes(seurat_obj)
  seurat_obj$patient_group <- factor(seurat_obj$patient_group,
    levels = c("Control", "CPI no colitis", "CPI colitis"),
    ordered = FALSE
  )
  donor_ids <- c(
    paste0("CT", 1:8),
    paste0("NC", 1:6),
    paste0("C", 1:8)
  )
  donor_ids <- donor_ids[donor_ids %in% seurat_obj$patient]

  seurat_obj$donor_id <- factor(seurat_obj$donor_id,
    levels = donor_ids,
    ordered = FALSE
  )
  return(seurat_obj)
}



get_donor_colors <- function(meta_data) {
  axis_colors <- sapply(levels(meta_data$donor_id), function(x) {
    if (grepl("^CT", x)) {
      return(ctr_col)
    } else if (grepl("^C[0-9]", x)) {
      return(cpi_col)
    } else if (grepl("^NC", x)) {
      return(cpi_nc_col)
    }
  })
  return(axis_colors)
}


colitis_pseudobulk <- function(seurat_obj) {
  pb_dat <- AggregateExpression(seurat_obj, assay = "RNA", return.seurat = TRUE, group.by = "donor_id")
  pb_dat$patient_group <- sapply(pb_dat$donor_id, function(x) {
    if (grepl("^CT", x)) {
      return("Control")
    } else if (grepl("^C[0-9]", x)) {
      return("CPI colitis")
    } else if (grepl("^NC", x)) {
      return("CPI no colitis")
    }
  })
  Idents(pb_dat) <- pb_dat$patient_group
  colitis_deg <- FindMarkers(pb_dat,
    ident.1 = "CPI colitis",
    ident.2 = c("CPI no colitis", "Control"),
    test.use = "DESeq2"
  )
  colitis_deg$gene_name <- rownames(colitis_deg)
  return(colitis_deg)
}



add_summary_metadata <- function(metadata) {
  metadata$sample <- word(metadata$patient_celltype, 1, 1, "_")
  metadata$celltype <- word(metadata$patient_celltype, 2, -1, "_")
  metadata$patient_group <- sapply(metadata$sample, function(x) {
    if (grepl("^CT", x)) {
      return("Control")
    } else if (grepl("^C[0-9]", x)) {
      return("CPI colitis")
    } else if (grepl("^NC", x)) {
      return("CPI no colitis")
    }
  })
  return(metadata)
}



get_patient_group_factor <- function(patient_group) {
  patient_group <- factor(patient_group, levels = c("Control", "CPI no colitis", "CPI colitis"), ordered = FALSE)
  return(patient_group)
}



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                                Import Data                               ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
output_dir <- file.path("output_figures/")

sc_all <- qread("data/processed/annotated_seurat.qs")
seurat <- qread("data/processed/annotated_Tcell.qs")
cd4 <- qread("data/processed/annotated_cd4.qs")
cd8 <- qread("data/processed/annotated_cd8.qs")

sc_all <- prep_data(sc_all)
seurat <- prep_data(seurat)
cd4 <- prep_data(cd4)
cd8 <- prep_data(cd8)

combined_tcr <- read.csv("data/processed/combined_tcr.csv")

sc_meta <- seurat@meta.data
sc_meta <- sc_meta[, c(
  "barcode_contig", "patient_group", "treatment", "age", "sex", "malignancy",
  "recist_response", "overal_surival_months_from_cpi_initiation",
  "diarrhea_grade_CTCAE", "mayo_endoscopic_score_MES",
  "time_from_first_treatment_to_symptom_onset_days", "infliximab_required_for_treatment",
  "azimuth_celltype_1", "azimuth_celltype_2", "azimuth_celltype_3",
  "annotation_level1", "annotation_level2", "annotation_level3"
)]
colnames(sc_meta)[1] <- "barcode"

combined_tcr_merge <- left_join(combined_tcr, sc_meta, by = "barcode")

split_tcr <- split(combined_tcr, f = combined_tcr_merge$donor_id)

# set colors for plotting
ctr_col <- "#00bfa0"
cpi_nc_col <- "#ffa300"
cpi_col <- "#9b19f5"

colors <- c(
  "Control" = ctr_col,
  "CPI no colitis" = cpi_nc_col,
  "CPI colitis" = cpi_col
)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                              Dataset overview                            ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## overview umap plots --------------------------------
p1 <- DimPlot(sc_all, group.by = c("annotation_level1"), label = TRUE, raster = TRUE) + ggtitle("Cell Type")
p2 <- DimPlot(sc_all, group.by = c("patient_group"), label = FALSE, raster = TRUE, cols = colors) + ggtitle("Patient Group")
p3 <- DimPlot(sc_all, group.by = c("flow_cell"), label = FALSE, raster = TRUE)
p4 <- DimPlot(sc_all, group.by = c("sex"), label = FALSE, raster = TRUE)

png(filename = paste0(output_dir,"overview_umaps.png"), width =720, height = 720, unit = "px")
(p1 + p2)/(p3 + p4)
dev.off()
## cd45 cell prop ---------------------
meta_45 <- sc_all@meta.data[sc_all@meta.data$flow_cell == "CD45", ]
meta_45$donor_id <- factor(meta_45$donor_id)

axis_colors <- get_donor_colors(meta_45)


p1 <- ggplot(meta_45, aes(x = donor_id, fill = annotation_level1)) +
  geom_bar(position = "fill") +
  theme_prism() +
  theme(axis.text.x = element_text(color = axis_colors, angle = 90, vjust = 0.5, hjust = 1)) +
  ylab("Cell Proportion in CD45+ Separated Cells") +
  xlab("Patient ID")

png(filename = paste0(output_dir,"cd45_cell_prop.png"), width =420, height = 360, unit = "px")
p1
dev.off()

## cluster split by patient group umap -------------
p1 <- DimPlot(sc_all[, sc_all@meta.data$flow_cell == "CD45"], 
        group.by = "RNA_snn_res.0.4", 
        split.by = "patient_group") +
  ggtitle("")

png(filename = paste0(output_dir,"cluster_split_patient_group.png"), width =720, height = 420, unit = "px")
p1
dev.off()

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                               Tcell Overview                             ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## T cell marker umap ----------------
p1 <- DimPlot(seurat, group.by = c("patient_group"), label = FALSE, raster = FALSE, cols = colors)
p2 <- DimPlot(seurat, group.by = c("flow_cell"), label = FALSE, raster = FALSE)
p3 <- DimPlot(seurat, group.by = c("sex"), label = FALSE, raster = FALSE)
p4 <- DimPlot(seurat, group.by = c("RNA_snn_res.1"), label = TRUE, raster = FALSE)

png(filename = paste0(output_dir,"tcell_overview_umaps.png"), width =720, height = 720, unit = "px")
(p1 + p2)/(p3 + p4)
dev.off()


p1 <- FeaturePlot(seurat, min.cutoff = "q1", max.cutoff = "q99", features = c("CD4", "CD8A", "CD8B"))
png(filename = paste0(output_dir,"tcell_cd_markers.png"), width =720, height = 720, unit = "px")
p1
dev.off()

# FeaturePlot(seurat, min.cutoff = "q1", max.cutoff = "q99", features = c("CD3E","CD3D","CD3G"))

## Azimuth umap -----------------
seurat$annotation_level2[is.na(seurat$annotation_level2)] <- "Other"
p1 <- DimPlot(seurat, 
              group.by = c("azimuth_celltype_1", "azimuth_celltype_2", "annotation_level2"), 
              label = FALSE, raster = FALSE)

png(filename = paste0(output_dir,"azimuth_umaps.png"), width =840, height = 360, unit = "px")
p1
dev.off()


## Azimuth cell prop -------------
axis_colors <- get_donor_colors(seurat@meta.data)
p1 <- ggplot(seurat@meta.data, aes(x = donor_id, fill = annotation_level2)) +
  geom_bar(position = "fill") +
  theme_prism() +
  theme(axis.text.x = element_text(color = axis_colors, angle = 90, vjust = 0.5, hjust = 1))

png(filename = paste0(output_dir,"azimuth_cell_prop.png"), width =420, height = 360, unit = "px")
p1
dev.off()


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                            Tcell subclustering                           ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## CD8 umap ------------------
DimPlot(cd8, group.by = c("patient_group"), label = TRUE, cols = colors) +
  ggtitle("Patient Group")
DimPlot(cd8, group.by = c("annotation_level3"), label = TRUE) +
  ggtitle("CD8+ T-cell sub-clustering")

## CD8 cell prop -------------
axis_colors <- get_donor_colors(cd8@meta.data)
ggplot(cd8@meta.data, aes(x = donor_id, fill = annotation_level3)) +
  geom_bar(position = "fill") +
  theme_prism() +
  theme(axis.text.x = element_text(color = axis_colors, angle = 90, vjust = 0.5, hjust = 1))

cd8_prop <- prop.table(table(cd8$patient_group, cd8$annotation_level3), margin = 1) %>%
  data.frame()
cd8_prop$Var1 <- factor(cd8_prop$Var1, levels = c("Control", "CPI no colitis", "CPI colitis"))
ggplot(cd8_prop, aes(x = Var2, y = Freq, fill = Var1)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  theme_prism() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
  scale_fill_manual(values = colors) +
  ylab("Cell proportion") +
  xlab("CD8 cell subtype")

## CD4 umap -----------------
DimPlot(cd4, group.by = c("patient_group"), label = TRUE, cols = colors) +
  ggtitle("Patient Group")
DimPlot(cd4, group.by = c("annotation_level3"), label = TRUE) +
  ggtitle("CD4+ T-cell sub-clustering")

## CD4 cell prop ---------------
axis_colors <- get_donor_colors(cd4@meta.data)
ggplot(cd4@meta.data, aes(x = donor_id, fill = annotation_level3)) +
  geom_bar(position = "fill") +
  theme_prism() +
  theme(axis.text.x = element_text(color = axis_colors, angle = 90, vjust = 0.5, hjust = 1)) +
  ylab("Cell Proportions")

cd4_prop <- prop.table(table(cd4$patient_group, cd4$annotation_level3), margin = 1) %>%
  data.frame()
cd4_prop$Var1 <- factor(cd4_prop$Var1, levels = c("Control", "CPI no colitis", "CPI colitis"))
ggplot(cd4_prop, aes(x = Var2, y = Freq, fill = Var1)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  theme_prism() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
  scale_fill_manual(values = colors) +
  ylab("Cell Proportions") +
  xlab("CD4 cell subtype")


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                         Checkpoint gene expression                       ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## cd8 checkpoint expr plots --------------------
# having a metadata column named CTLA4 causes errors when trying to plot CTLA4 gene expression
cd8_plot <- cd8
cd8_plot$CTLA4 <- NULL

cd4_plot <- cd4
cd4_plot$CTLA4 <- NULL

checkpoint_genes <- c("CTLA4", "PDCD1", "HAVCR2", "LAG3", "TIGIT", "CD38")

VlnPlot(cd8_plot,
  features = checkpoint_genes, pt.size = 0,
  group.by = "patient_group", cols = colors, adjust = 2
)

VlnPlot(cd4_plot,
  features = checkpoint_genes, pt.size = 0,
  group.by = "patient_group", cols = colors, adjust = 2
)

# VlnPlot(cd8_plot, features = checkpoint_genes, pt.size = 0,
#         group.by = "annotation_level3")
#
# VlnPlot(cd4_plot, features = checkpoint_genes, pt.size = 0,
#         group.by = "annotation_level3")


VlnPlot(cd8_plot,
  features = checkpoint_genes, pt.size = 0,
  group.by = "annotation_level3", split.by = "patient_group",
  stack = TRUE, flip = TRUE, cols = colors
)

VlnPlot(cd4_plot,
  features = checkpoint_genes, pt.size = 0,
  group.by = "annotation_level3", split.by = "patient_group",
  stack = TRUE, flip = TRUE, cols = colors
)

FeaturePlot(cd8_plot, min.cutoff = "q1", max.cutoff = "q99", features = checkpoint_genes)
FeaturePlot(cd4_plot, min.cutoff = "q1", max.cutoff = "q99", features = checkpoint_genes)

## pseudobulk prep ----------------------
cd8_pb <- colitis_pseudobulk(cd8)
cd4_pb <- colitis_pseudobulk(cd4)

cd8_pb[checkpoint_genes, ]
cd4_pb[checkpoint_genes, ]

de_merge <- inner_join(cd4_pb, cd8_pb, by = "gene_name", suffix = c(".CD4", ".CD8"))
de_merge <- de_merge[!is.na(de_merge$p_val_adj.CD4) & !is.na(de_merge$p_val_adj.CD8), ]

## pseudobulk comparison plot ------------------
ggplot(de_merge, aes(x = avg_log2FC.CD4, y = avg_log2FC.CD8)) +
  geom_point() +
  theme_prism() +
  stat_poly_line() +
  stat_poly_eq() +
  geom_text_repel(
    data = subset(
      de_merge,
      avg_log2FC.CD4 > 2 |
        avg_log2FC.CD4 < -0.5 |
        avg_log2FC.CD8 > 2 |
        avg_log2FC.CD8 < -1.2
    ),
    mapping = aes(x = avg_log2FC.CD4, y = avg_log2FC.CD8, label = gene_name)
  ) +
  xlab("CD4+ log2FC") +
  ylab("CD8+ log2FC")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                                  TCR stats                               ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## TCR prep ------------------
cd8_tcr <- combined_tcr_merge[combined_tcr_merge$annotation_level2 %in% "CD8_Tcell", ]
cd8_split_tcr <- split(cd8_tcr, f = combined_tcr_merge$donor_id)
cd8_merge <- cd8[, cd8$flow_cell == "CD3"]
cd8_merge <- combineExpression(cd8_split_tcr,
  cd8_merge,
  cloneCall = "gene",
  proportion = TRUE,
  group.by = "patient_group",
  filterNA = TRUE
)
cd8_merge$cloneSize <- as.character(cd8_merge$cloneSize)
cd8_merge$cloneSize[cd8_merge$clonalFrequency == 1] <- "Single"
cd8_merge$expanded <- ifelse(cd8_merge$cloneSize == "Single", "Single", "Expanded")
Idents(cd8_merge) <- cd8_merge$annotation_level3

cd4_tcr <- combined_tcr_merge[combined_tcr_merge$annotation_level2 %in% "CD4_Tcell", ]
cd4_split_tcr <- split(cd4_tcr, f = combined_tcr_merge$donor_id)
cd4_merge <- cd4[, cd4$flow_cell == "CD3"]
cd4_merge <- combineExpression(cd4_split_tcr,
  cd4_merge,
  cloneCall = "gene",
  proportion = TRUE,
  group.by = "patient_group",
  filterNA = TRUE
)
cd4_merge$cloneSize <- as.character(cd4_merge$cloneSize)
cd4_merge$cloneSize[cd4_merge$clonalFrequency == 1] <- "Single"
cd4_merge$expanded <- ifelse(cd4_merge$cloneSize == "Single", "Single", "Expanded")
Idents(cd4_merge) <- cd4_merge$annotation_level3

## TCR stats prep --------------------
cd_all_meta <- rbind(cd8_merge@meta.data, cd4_merge@meta.data)
cd_all_meta$patient_celltype <- paste0(cd_all_meta$patient, "_", cd_all_meta$azimuth_celltype_1)
cd_all_meta$patient_subcelltype <- paste0(cd_all_meta$patient, "_", cd_all_meta$annotation_level3)

unique_clones <- cd_all_meta %>%
  group_by(patient_celltype) %>%
  summarize(unique_clones = n_distinct(CTaa))
unique_clones <- add_summary_metadata(unique_clones)
unique_clones$patient_group <- get_patient_group_factor(unique_clones$patient_group)


expansion_ratio <- cd_all_meta %>%
  group_by(patient_celltype) %>%
  summarize(
    expanded_unique = n_distinct(CTaa[clonalFrequency > 1]),
    single_count = sum(clonalFrequency == 1),
    ratio = expanded_unique / (single_count + expanded_unique)
  )
expansion_ratio <- add_summary_metadata(expansion_ratio)
expansion_ratio$patient_group <- get_patient_group_factor(expansion_ratio$patient_group)

cd8_expanded <- cd_all_meta %>%
  dplyr::filter(azimuth_celltype_1 == "CD8 T") %>%
  group_by(patient_subcelltype) %>%
  summarize(
    expanded_unique = n_distinct(CTaa[clonalFrequency > 1]),
    single_count = sum(clonalFrequency == 1),
    ratio = expanded_unique / (single_count + expanded_unique)
  )
colnames(cd8_expanded)[1] <- "patient_celltype"
cd8_expanded <- add_summary_metadata(cd8_expanded)
cd8_expanded$patient_group <- get_patient_group_factor(cd8_expanded$patient_group)

cd4_expanded <- cd_all_meta %>%
  dplyr::filter(azimuth_celltype_1 == "CD4 T") %>%
  group_by(patient_subcelltype) %>%
  summarize(
    expanded_unique = n_distinct(CTaa[clonalFrequency > 1]),
    single_count = sum(clonalFrequency == 1),
    ratio = expanded_unique / (single_count + expanded_unique)
  )
colnames(cd4_expanded)[1] <- "patient_celltype"
cd4_expanded <- add_summary_metadata(cd4_expanded)
cd4_expanded$patient_group <- get_patient_group_factor(cd4_expanded$patient_group)

## TCR stats plots --------------------
ggplot(unique_clones, aes(x = patient_group, y = unique_clones, fill = celltype)) +
  geom_boxplot() +
  theme_prism() +
  theme(axis.text.x = element_text(color = colors, angle = 90, vjust = 0.5, hjust = 1)) +
  ylab("# unique clonotypes per patient") +
  xlab("")

ggplot(expansion_ratio, aes(x = patient_group, y = expanded_unique, fill = celltype)) +
  geom_boxplot() +
  theme_prism() +
  theme(axis.text.x = element_text(color = colors, angle = 90, vjust = 0.5, hjust = 1)) +
  ylab("# unique expanded clonotypes\nper patient") +
  xlab("")

ggplot(expansion_ratio, aes(x = patient_group, y = ratio, fill = celltype)) +
  geom_boxplot() +
  theme_prism() +
  theme(axis.text.x = element_text(color = colors, angle = 90, vjust = 0.5, hjust = 1)) +
  ylab("Expanded clonotypes\n(% of total) per patient") +
  xlab("")


ggplot(cd8_expanded, aes(x = celltype, y = expanded_unique, fill = patient_group)) +
  geom_boxplot() +
  theme_prism() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  ylab("# unique expanded clonotypes\nper patient per celltype") +
  xlab("") +
  scale_fill_manual(values = colors)

ggplot(cd4_expanded, aes(x = celltype, y = expanded_unique, fill = patient_group)) +
  geom_boxplot() +
  theme_prism() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  ylab("# unique expanded clonotypes\nper patient per celltype") +
  xlab("") +
  scale_fill_manual(values = colors)

## circle plot ---------------------
circles <- getCirclize(cd8_merge, group.by = "annotation_level3")
grid.cols <- hue_pal()(length(unique(cd8_merge$annotation_level3)))
names(grid.cols) <- unique(cd8_merge$annotation_level3)

# Graphing the chord diagram
chordDiagram(circles, self.link = 1, grid.col = grid.cols)


## clonal diversity and overlap -------------------
clonalDiversity(cd8_split_tcr, cloneCall = "aa", group.by = "patient_group")
clonalOverlap(cd8_split_tcr, cloneCall = "aa", group.by = "annotation_level3", method = "raw") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

## top AA plot ------------------------
top_aa <- cd8_merge@meta.data %>%
  dplyr::filter(annotation_level3 == "Cytotoxic_effector") %>%
  pull(CTaa) %>%
  table()
top_aa <- top_aa[order(top_aa, decreasing = TRUE)]

cd8_merge$top_aa <- NA
cd8_merge$top_aa[cd8_merge$CTaa %in% names(top_aa)[1:10]] <- cd8_merge$CTaa[cd8_merge$CTaa %in% names(top_aa)[1:10]]

cells_list <- list()
for (i in 1:20) {
  cells_list[[i]] <- which(cd8_merge$CTaa == names(top_aa)[i])
}
names(cells_list) <- names(top_aa)[1:20]


DimPlot(cd8_merge, group.by = "top_aa", raster = FALSE, pt.size = 4, na.value = "black")
DimPlot(cd8_merge,
  cells.highlight = cells_list, cols.highlight = hue_pal()(length(cells_list)),
  raster = FALSE, sizes.highlight = 3
) +
  theme(
    legend.position = "bottom",
    legend.text = element_text(size = 9),
    legend.key.size = unit(0.3, "cm"),
    legend.spacing.x = unit(0.2, "cm"),
    legend.margin = margin(0, 0, 0, 0),
    legend.box = "horizontal",
    legend.text.align = 0
  ) +
  guides(color = guide_legend(nrow = 7))
