#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 6) {
  stop(
    "Usage: pyscenic_select_cells.R <seurat_rds> <out_tsv> <group_by> ",
    "<target_cells> <min_cells_per_group> <seed>"
  )
}

seurat_rds <- args[[1]]
out_tsv <- args[[2]]
group_by <- args[[3]]
target_cells <- as.integer(args[[4]])
min_cells_per_group <- as.integer(args[[5]])
seed <- as.integer(args[[6]])

if (is.na(target_cells) || target_cells < 1) {
  stop("target_cells must be a positive integer")
}
if (is.na(min_cells_per_group) || min_cells_per_group < 0) {
  stop("min_cells_per_group must be a non-negative integer")
}
if (is.na(seed)) {
  stop("seed must be an integer")
}

seurat_obj <- readRDS(seurat_rds)
if (!group_by %in% colnames(seurat_obj@meta.data)) {
  stop("Missing metadata column for balanced pySCENIC sampling: ", group_by)
}

metadata <- data.frame(
  cell_id = colnames(seurat_obj),
  group = as.character(seurat_obj@meta.data[[group_by]]),
  stringsAsFactors = FALSE
)
metadata <- metadata[!is.na(metadata$group) & metadata$group != "", , drop = FALSE]
if (nrow(metadata) == 0) {
  stop("No cells with non-empty metadata values in column: ", group_by)
}

target_cells <- min(target_cells, nrow(metadata))
group_sizes <- sort(table(metadata$group), decreasing = TRUE)

allocate_cells <- function(group_sizes, target_cells, min_cells_per_group) {
  sizes <- as.integer(group_sizes)
  names(sizes) <- names(group_sizes)

  if (target_cells >= sum(sizes)) {
    return(sizes)
  }

  allocation <- pmin(sizes, min_cells_per_group)
  if (sum(allocation) > target_cells) {
    raw_allocation <- target_cells * sizes / sum(sizes)
    allocation <- pmax(1L, floor(raw_allocation))
    while (sum(allocation) > target_cells) {
      candidates <- names(allocation)[allocation > 1L]
      remove_from <- candidates[which.min(raw_allocation[candidates] - floor(raw_allocation[candidates]))]
      allocation[remove_from] <- allocation[remove_from] - 1L
    }
    while (sum(allocation) < target_cells) {
      capacity <- sizes - allocation
      add_to <- names(capacity)[which.max(capacity)]
      if (capacity[[add_to]] <= 0L) {
        break
      }
      allocation[add_to] <- allocation[add_to] + 1L
    }
    return(allocation)
  }

  remaining <- target_cells - sum(allocation)
  capacity <- sizes - allocation
  while (remaining > 0 && any(capacity > 0)) {
    raw_extra <- remaining * capacity / sum(capacity)
    extra <- pmin(capacity, floor(raw_extra))
    if (sum(extra) == 0L) {
      add_to <- names(capacity)[which.max(capacity)]
      extra[add_to] <- 1L
    }
    allocation <- allocation + extra
    remaining <- target_cells - sum(allocation)
    capacity <- sizes - allocation
  }

  allocation
}

allocation <- allocate_cells(group_sizes, target_cells, min_cells_per_group)

set.seed(seed)
selected <- unlist(
  lapply(names(allocation), function(group_name) {
    group_cells <- metadata$cell_id[metadata$group == group_name]
    sample(group_cells, allocation[[group_name]])
  }),
  use.names = FALSE
)

selected_metadata <- metadata[match(selected, metadata$cell_id), , drop = FALSE]
dir.create(dirname(out_tsv), recursive = TRUE, showWarnings = FALSE)
write.table(
  selected_metadata,
  file = out_tsv,
  quote = FALSE,
  sep = "\t",
  row.names = FALSE
)

message("Selected ", nrow(selected_metadata), " cells for pySCENIC from ", nrow(metadata), " available cells.")
message("Balanced by metadata column: ", group_by)
print(table(selected_metadata$group))
