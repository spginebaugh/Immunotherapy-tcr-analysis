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

library(scRepertoire)
library(ggraph)

source("code/utils/TCR_utils.R")
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                                Import data                               ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
seurat <- qread("data/processed/annotated_Tcell.qs")
cd4 <- qread("data/processed/annotated_cd4.qs")
cd8 <- qread("data/processed/annotated_cd8.qs")

seurat <- add_contig_barcodes(seurat)

filt_contig <- read.csv("data/GSE144469/GSE144469_TCR_filtered_contig_annotations_all.csv", row.names = "X")
filt_contig <- filt_contig[filt_contig$barcode %in% seurat$barcode_contig,]

## read in additional contig annotations
contig_path <- "data/GSE144469/GSE144469_RAW/tcr_csvs/"
contig_files <- list.files(contig_path)
contig_names <- contig_files %>% word(.,2,2,"_") %>% word(.,1,1,"-")
contig_list <- list()
for (i in 1:length(contig_files)){
  contig_tmp <- read.csv(paste0(contig_path, contig_files[i]))
  contig_tmp$barcode <- paste0(word(contig_tmp$barcode,1,1,"-"), "-", contig_names[i])
  contig_tmp$contig_id <- paste0(contig_tmp$barcode, "_", word(contig_tmp$contig_id,2,3,"_"))
  contig_list[[contig_names[i]]] <- contig_tmp
}
filt_contig_2 <- do.call(rbind, contig_list)
filt_contig_2 <- filt_contig_2[filt_contig_2$barcode %in% seurat$barcode_contig,]
filt_contig_2 <- filt_contig_2[!(filt_contig_2$contig_id %in% filt_contig$contig_id),]

filt_contig <- rbind(filt_contig, filt_contig_2)
rm(filt_contig_2)


combined_tcr <- combineTCR(filt_contig)[["S1"]]
combined_tcr$donor_id <- word(combined_tcr$barcode,2,2,"-")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                                  save tcr                                ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
write.csv(combined_tcr, "data/processed/combined_tcr.csv", row.names = FALSE)







