#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
})

args <- commandArgs(trailingOnly = TRUE)
if (!length(args) %in% c(9, 10)) {
  stop(
    "Usage: pyscenic_export.R <seurat_rds> <expr_csv> <metadata_tsv> ",
    "<assay> <layer> <min_counts_per_gene> <min_cells_per_gene> ",
    "<downsample_cells> <export_chunk_cells> [cell_list_tsv]"
  )
}

seurat_rds <- args[[1]]
expr_csv <- args[[2]]
metadata_tsv <- args[[3]]
assay_name <- args[[4]]
layer_name <- args[[5]]
min_counts_per_gene <- as.numeric(args[[6]])
min_cells_per_gene <- as.numeric(args[[7]])
downsample_cells <- as.integer(args[[8]])
export_chunk_cells <- as.integer(args[[9]])
cell_list_tsv <- if (length(args) >= 10) args[[10]] else ""

if (is.na(export_chunk_cells) || export_chunk_cells < 1) {
  stop("export_chunk_cells must be a positive integer")
}

seurat_obj <- readRDS(seurat_rds)
if (!assay_name %in% names(seurat_obj@assays)) {
  stop("Assay not found in Seurat object: ", assay_name)
}

expr_matrix <- tryCatch(
  GetAssayData(seurat_obj, assay = assay_name, layer = layer_name),
  error = function(e) GetAssayData(seurat_obj, assay = assay_name, slot = layer_name)
)

if (!is.na(cell_list_tsv) && nzchar(cell_list_tsv)) {
  if (!file.exists(cell_list_tsv)) {
    stop("Cell list does not exist: ", cell_list_tsv)
  }
  selected_cells <- read.table(
    cell_list_tsv,
    header = TRUE,
    sep = "\t",
    stringsAsFactors = FALSE,
    check.names = FALSE
  )[[1]]
  selected_cells <- intersect(selected_cells, colnames(expr_matrix))
  if (length(selected_cells) == 0) {
    stop("No overlap between cell list and Seurat object cells: ", cell_list_tsv)
  }
  expr_matrix <- expr_matrix[, selected_cells, drop = FALSE]
  seurat_obj <- subset(seurat_obj, cells = selected_cells)
  message("Using selected cell list: ", cell_list_tsv)
}

if (downsample_cells > 0 && ncol(expr_matrix) > downsample_cells) {
  set.seed(1234)
  cells_to_keep <- sample(colnames(expr_matrix), downsample_cells)
  expr_matrix <- expr_matrix[, cells_to_keep, drop = FALSE]
  seurat_obj <- subset(seurat_obj, cells = cells_to_keep)
}

genes_to_keep <- Matrix::rowSums(expr_matrix) >= min_counts_per_gene &
  Matrix::rowSums(expr_matrix > 0) >= min_cells_per_gene
expr_matrix <- expr_matrix[genes_to_keep, , drop = FALSE]

message("Exporting ", nrow(expr_matrix), " genes x ", ncol(expr_matrix), " cells")

write_csv_field <- function(values, connection) {
  write.table(
    values,
    file = connection,
    sep = ",",
    quote = FALSE,
    row.names = FALSE,
    col.names = FALSE,
    append = TRUE
  )
}

if (file.exists(expr_csv)) {
  file.remove(expr_csv)
}
con <- file(expr_csv, open = "wt")
writeLines(paste(c("", rownames(expr_matrix)), collapse = ","), con = con)

cell_indices <- seq_len(ncol(expr_matrix))
chunk_starts <- seq(1, length(cell_indices), by = export_chunk_cells)
for (start_idx in chunk_starts) {
  end_idx <- min(start_idx + export_chunk_cells - 1, length(cell_indices))
  chunk_idx <- cell_indices[start_idx:end_idx]
  chunk_matrix <- as.matrix(t(expr_matrix[, chunk_idx, drop = FALSE]))
  chunk_df <- data.frame(cell_id = rownames(chunk_matrix), chunk_matrix, check.names = FALSE)
  write_csv_field(chunk_df, con)
  rm(chunk_matrix, chunk_df)
  gc(verbose = FALSE)
}
close(con)

metadata_columns <- intersect(
  c("barcode", "orig.ident", "group", "sample_name", "cell_type", "nCount_RNA", "nFeature_RNA"),
  colnames(seurat_obj@meta.data)
)
metadata_df <- seurat_obj@meta.data[, metadata_columns, drop = FALSE]
metadata_df <- data.frame(cell_id = rownames(metadata_df), metadata_df, check.names = FALSE)
write.table(
  metadata_df,
  file = metadata_tsv,
  quote = FALSE,
  sep = "\t",
  row.names = FALSE
)
