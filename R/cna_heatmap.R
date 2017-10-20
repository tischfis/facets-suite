cna_heatmap <- function(gene_level_calls){
  mat <- as.matrix(
    dcast.data.table(
      gene_level_calls,
      Tumor_Sample_Barcode ~ Hugo_Symbol,
      value.var = "tcn"
    )[,-1]
  )
  mat <- log(mat/ 2 + 0.01)
  gplots::heatmap.2(mat, trace = "none")
}
