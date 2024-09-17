# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                              Load Libraries                              ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
library(Seurat)
library(qs)
library(ggplot2)
library(magrittr)
library(dplyr)
library(stringr)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                          Import Cellbender Outputs                       ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
seurat <- qread("data/processed/doublets_detected.qs")

seurat$patient <- word(seurat$sample,2,2,"_") %>% word(.,1,1,"-")
seurat$flow_cell <- word(seurat$sample,2,2,"-")
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                              Add in Metadata                             ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
meta_1 <- read.csv("data/metadata/patient_metadata.csv")
meta_1[] <- lapply(meta_1, function(x) word(x,1,1, fixed(" (")))
meta_1[] <- lapply(meta_1, function(x) word(x,1,1, fixed("(")))

meta_2 <- read.csv("data/metadata/supp_table_2.csv", na.strings = c("","NA"))

meta_3 <- read.csv("data/metadata/supp_table_3.csv")

meta_4 <- read.csv("data/metadata/supp_table_4.csv")

## merge all -------------------------------------------------------------------
metadata <- seurat@meta.data
metadata <- left_join(metadata, meta_1, by = "patient")
metadata <- left_join(metadata, meta_2, by = "patient")
metadata <- left_join(metadata, meta_3, by = "patient")
metadata <- left_join(metadata, meta_4, by = "patient")
rownames(metadata) <- metadata$barcode
seurat@meta.data <- metadata


## add QC to metadata ---------------------------------------------------------
seurat$mitoRatio <- PercentageFeatureSet(object = seurat, pattern = "^MT-", assay = "RNA") / 100
seurat$riboRatio <- PercentageFeatureSet(object = seurat, 
                                                pattern = "^RP[SL][[:digit:]]|^RPLP[[:digit:]]|^RPSA", 
                                                assay = "RNA") / 100
seurat$log10GenesPerUMI <- log10(seurat$nFeature_RNA) / log10(seurat$nCount_RNA)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                                     QC                                   ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
metadata <- seurat@meta.data
metadata %>%
  ggplot(aes(x = sample, fill = sample)) +
  geom_bar() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold")) +
  ggtitle("NCells")

metadata %>%
  ggplot(aes(color = sample, x = nFeature_RNA, fill = sample)) +
  geom_density(alpha = 0.2) +
  scale_x_log10() +
  theme_classic() +
  ylab("Cell density") +
  geom_vline(xintercept = 800)

metadata %>%
  ggplot(aes(color = sample, x = nCount_RNA, fill = sample)) +
  geom_density(alpha = 0.2) +
  scale_x_log10() +
  theme_classic() +
  ylab("Cell density") +
  geom_vline(xintercept = 1400)

metadata %>%
  ggplot(aes(color = sample, x = mitoRatio, fill = sample)) +
  geom_density(alpha = 0.2) +
  scale_x_log10() +
  theme_classic() +
  geom_vline(xintercept = 0.15)

metadata %>%
  ggplot(aes(x = sample, log10GenesPerUMI, fill = sample)) +
  geom_violin() +
  geom_boxplot(width = 0.1, fill = alpha("white", 0.7)) +
  theme_classic() +
  geom_hline(yintercept = c(0.80)) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold")) +
  NoLegend()

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                                  Filtering                               ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
seurat <- subset(seurat, nFeature_RNA > 800 & nCount_RNA > 1400 & mitoRatio < 0.15)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                                  Save Data                               ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
qsave(seurat, "data/processed/filtered_seurat.qs")
