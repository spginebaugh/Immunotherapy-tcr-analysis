# pulls information from seurat_obj barcode name and adds it to seurat_obj@meta.data
add_contig_barcodes <- function(seurat_obj){
  seurat_obj$donor_id <- word(seurat_obj$sample,2,2,"_") %>% word(.,1,1,"-")
  barcode_stripped <- word(seurat_obj$barcode,-1,-1,"_") %>% word(.,1,1,"-")
  barcode_contig <- paste0(barcode_stripped,"-",seurat_obj$donor_id)
  barcode_contig[seurat_obj$flow_cell == "CD45"] <- paste0(barcode_contig[seurat_obj$flow_cell == "CD45"],"_CD45")
  seurat_obj$barcode_contig <- barcode_contig
  colnames(seurat_obj) <- seurat_obj$barcode_contig
  return(seurat_obj)
}