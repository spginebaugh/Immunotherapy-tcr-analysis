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
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                                Import data                               ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
seurat <- qread("data/processed/clustered_tcell.qs")
seurat$donor_id <- word(seurat$sample,2,2,"_") %>% word(.,1,1,"-")
barcode_stripped <- word(seurat$barcode,-1,-1,"_") %>% word(.,1,1,"-")
seurat$barcode_contig <- paste0(barcode_stripped,"-",seurat$donor_id)

# only CD3 filtered cells have TCR data
sc_cd3 <- seurat[,seurat$flow_cell == "CD3" & seurat$azimuth_celltype_1 %in% c("CD4 T", "CD8 T")]


# contig_files <- list.files("data/GSE144469/GSE144469_RAW/tcr_csvs/")
# contig_names <- word(contig_files,1,1,"-")
# 
# contig_list <- list()
# for (i in 1:length(contig_files)){
#   contig_list[[contig_names[i]]] <- read.csv(paste0("data/GSE144469/GSE144469_RAW/tcr_csvs/",contig_files[i]))
# }
# all_contig <- do.call(rbind, contig_list)



filt_contig <- read.csv("data/GSE144469/GSE144469_TCR_filtered_contig_annotations_all.csv", row.names = "X")
filt_contig <- filt_contig[filt_contig$barcode %in% sc_cd3$barcode_contig,]

sc_meta <- sc_cd3@meta.data
sc_meta <- sc_meta[,c("barcode_contig","patient_group","treatment","age","sex","malignancy",
                      "recist_response","overal_surival_months_from_cpi_initiation",
                      "diarrhea_grade_CTCAE","mayo_endoscopic_score_MES",
                      "time_from_first_treatment_to_symptom_onset_days","infliximab_required_for_treatment",
                      "azimuth_celltype_1","azimuth_celltype_2","azimuth_celltype_3")]
colnames(sc_meta)[1] <- "barcode"

combined_tcr <- combineTCR(filt_contig)[["S1"]]
combined_tcr$donor_id <- word(combined_tcr$barcode,2,2,"-")
combined_tcr_merge <- left_join(combined_tcr, sc_meta, by = "barcode")

split_tcr <- split(combined_tcr_merge, f = combined_tcr_merge$donor_id)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                         Basic Clonal Visualization                       ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
clonalQuant(split_tcr, cloneCall = "strict", chain = "both", scale = TRUE)
clonalQuant(split_tcr, cloneCall = "strict", chain = "both", scale = TRUE, group.by = "patient_group")
clonalQuant(split_tcr, cloneCall = "strict", chain = "both", scale = TRUE, group.by = "azimuth_celltype_1")

clonalAbundance(split_tcr, cloneCall = "strict", scale = TRUE)
clonalAbundance(split_tcr, cloneCall = "strict", scale = TRUE, group.by = "patient_group")
clonalAbundance(split_tcr, cloneCall = "strict", scale = TRUE, group.by = "azimuth_celltype_1")

clonalLength(split_tcr, cloneCall = "aa", chain = "both", scale = TRUE, group.by = "patient_group")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                        combine clones with scRNA obj                     ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
sc_merge <- sc_cd3
colnames(sc_merge) <- sc_merge$barcode_contig
# sc_merge <- combineExpression(split_tcr,
#                               sc_merge,
#                               cloneCall = "strict",
#                               proportion = TRUE,
#                               group.by = "donor_id")

sc_merge <- combineExpression(split_tcr,
                         sc_merge,
                         cloneCall = "gene",
                         proportion = FALSE,
                         group.by = "donor_id",
                         cloneSize = c(Single=1, Small=5, Medium=20, Large=100, Hyperexpanded=500))
DimPlot(sc_merge, group.by = "cloneSize")

clonalOverlay(sc_merge,
              reduction = "umap", 
              cutpoint = 1, 
              bins = 10, 
              facet.by = "patient_group") + 
  guides(color = "none")

clonalNetwork(sc_merge, 
              reduction = "umap", 
              group.by = "azimuth_celltype_3",
              filter.clones = NULL,
              filter.identity = "CD8 Naive_2",
              cloneCall = "aa")


clonalOccupy(sc_merge, 
             x.axis = "azimuth_celltype_3")














