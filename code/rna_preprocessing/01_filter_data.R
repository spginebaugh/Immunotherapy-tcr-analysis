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