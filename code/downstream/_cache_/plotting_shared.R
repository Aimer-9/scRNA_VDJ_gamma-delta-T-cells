group_levels <- c("Naive", "Blank", "AB3", "ZOL", "PAN", "MSH2")

sample_name_levels <- c(
  "Naive_4", "Naive_5", "Naive_7", "Blank_3",
  "Blank_4", "Blank_5", "AB3_1", "AB3_2", "AB3_3",
  "ZOL_6", "ZOL_7", "ZOL_8", "PAN_1", "PAN_2", "PAN_3",
  "MSH2_3", "MSH2_4", "MSH2_5"
)

cell_type_levels <- c(
  "Effector Memory Vd2",
  "Pre-activated Vd2",
  "ZOL Effector Vd2",
  "ZOL FOXP3+ Vd2",
  "PAN Effector Vd2",
  "Naive Vd1",
  "Pre-activated Vd1",
  "Effector Vd1"
)

color_group <- c(
  "Naive" = "#1f77b4",
  "Blank" = "#ff7f0e",
  "AB3" = "#2ca02c",
  "ZOL" = "#17becf",
  "PAN" = "#d62728",
  "MSH2" = "#9467bd"
)

color_sample_name <- c(
  "Naive_4" = "#6BAED6",
  "Naive_5" = "#3182BD",
  "Naive_7" = "#08519C",
  "Blank_3" = "#FDBF6F",
  "Blank_4" = "#FF7F00",
  "Blank_5" = "#B15928",
  "AB3_1" = "#A1D99B",
  "AB3_2" = "#31A354",
  "AB3_3" = "#006D2C",
  "ZOL_6" = "#9EDAE5",
  "ZOL_7" = "#17BECF",
  "ZOL_8" = "#087E8B",
  "PAN_1" = "#FB9A99",
  "PAN_2" = "#E31A1C",
  "PAN_3" = "#A50F15",
  "MSH2_3" = "#BCBDDC",
  "MSH2_4" = "#756BB1",
  "MSH2_5" = "#54278F"
)

color_celltype <- c(
  "Effector Memory Vd2" = "#D55E00",
  "Pre-activated Vd2" = "#E69F00",
  "ZOL Effector Vd2" = "#CC79A7",
  "ZOL FOXP3+ Vd2" = "#7F3C8D",
  "PAN Effector Vd2" = "#E15759",
  "Naive Vd1" = "#009E73",
  "Pre-activated Vd1" = "#56B4E9",
  "Effector Vd1" = "#0072B2"
)

color_cellcycle <- c(
  "G1" = "#4E79A7",
  "S" = "#F28E2B",
  "G2M" = "#59A14F"
)

# Shared script setup and data-access helpers used by multiple downstream
# analyses. Keep analysis-specific statistics and plotting in the numbered
# scripts so the workflow remains easy to audit.
read_rds_checked <- function(path, label) {
  if (!file.exists(path)) {
    stop("Missing ", label, ": ", normalizePath(path, mustWork = FALSE), call. = FALSE)
  }
  if (file.info(path)$size == 0) {
    stop("Empty ", label, ": ", normalizePath(path, mustWork = FALSE), call. = FALSE)
  }
  readRDS(path)
}

read_optional_rds <- function(path, label) {
  if (!file.exists(path) || file.info(path)$size == 0) {
    message("[SKIP] Missing optional ", label, ": ", normalizePath(path, mustWork = FALSE))
    return(NULL)
  }
  readRDS(path)
}

normalise_metadata_levels <- function(seurat_obj) {
  seurat_obj$group <- factor(seurat_obj$group, levels = group_levels)
  seurat_obj$sample_name <- factor(seurat_obj$sample_name, levels = sample_name_levels)
  seurat_obj$cell_type <- factor(seurat_obj$cell_type, levels = cell_type_levels)
  seurat_obj
}

normalise_annotation_levels <- function(annotation) {
  annotation %>%
    dplyr::mutate(
      group = factor(group, levels = group_levels),
      sample_name = factor(sample_name, levels = sample_name_levels),
      cell_type = factor(cell_type, levels = cell_type_levels)
    )
}

get_assay_names <- function(seurat_obj) {
  assay_names <- names(seurat_obj@assays)
  if (is.null(assay_names)) {
    assay_names <- character()
  }
  as.character(assay_names)
}

get_reduction_names <- function(seurat_obj) {
  reduction_names <- names(seurat_obj@reductions)
  if (is.null(reduction_names)) {
    reduction_names <- character()
  }
  as.character(reduction_names)
}

join_assay_layers_if_needed <- function(seurat_obj, assay = SeuratObject::DefaultAssay(seurat_obj), label = NULL) {
  if (!assay %in% get_assay_names(seurat_obj)) {
    return(seurat_obj)
  }
  assay_obj <- seurat_obj[[assay]]
  assay_layers <- tryCatch(SeuratObject::Layers(assay_obj), error = function(e) character())
  layer_prefixes <- sub("\\..*$", "", assay_layers)
  has_split_layers <- length(unique(assay_layers)) > length(unique(layer_prefixes))
  has_multiple_same_type_layers <- any(table(layer_prefixes) > 1)

  if (inherits(assay_obj, "Assay5") && (has_split_layers || has_multiple_same_type_layers)) {
    suffix <- if (!is.null(label)) paste0(" for ", label) else ""
    message("[RUN] Joining Seurat v5 ", assay, " assay layers", suffix, ".")
    seurat_obj <- JoinLayers(seurat_obj, assay = assay)
  }
  seurat_obj
}

prepare_expression_object <- function(seurat_obj, assay = "RNA", label = NULL) {
  if (assay %in% get_assay_names(seurat_obj)) {
    DefaultAssay(seurat_obj) <- assay
  }
  join_assay_layers_if_needed(seurat_obj, assay = DefaultAssay(seurat_obj), label = label)
}

get_umap_reduction <- function(seurat_obj) {
  if ("umap.unintegrated" %in% get_reduction_names(seurat_obj)) {
    return("umap.unintegrated")
  }
  if ("umap" %in% get_reduction_names(seurat_obj)) {
    return("umap")
  }
  stop("No UMAP reduction found. Expected `umap.unintegrated` or `umap`.", call. = FALSE)
}

present_genes <- function(seurat_obj, genes) {
  genes[genes %in% rownames(seurat_obj)]
}

check_features <- function(seurat_obj, features) {
  missing_features <- setdiff(features, rownames(seurat_obj))
  if (length(missing_features) > 0) {
    stop("Missing feature(s) in Seurat object: ", paste(missing_features, collapse = ", "), call. = FALSE)
  }
  features
}

plot_empty <- function(title, subtitle = NULL) {
  ggplot2::ggplot() +
    ggplot2::theme_void() +
    ggplot2::labs(title = title, subtitle = subtitle)
}

output_files_exist <- function(paths) {
  all(file.exists(paths) & file.info(paths)$size > 0)
}

skip_existing_output <- function(paths, label) {
  if (output_files_exist(paths)) {
    message("[SKIP] Existing ", label, ": ", paste(paths, collapse = ", "))
    return(TRUE)
  }
  FALSE
}

save_plot <- function(plot, filename, width, height, dpi = 300, overwrite = FALSE) {
  output_paths <- file.path(figure_dir, paste0(filename, c(".png", ".pdf")))
  if (!overwrite && skip_existing_output(output_paths, "plot")) {
    return(invisible(NULL))
  }

  ggsave(output_paths[1],
    plot = plot,
    width = width,
    height = height,
    dpi = dpi
  )
  ggsave(output_paths[2],
    plot = plot,
    width = width,
    height = height
  )
}

save_heatmap <- function(heatmap, filename, width, height, overwrite = FALSE) {
  output_paths <- file.path(figure_dir, paste0(filename, c(".png", ".pdf")))
  if (!overwrite && skip_existing_output(output_paths, "heatmap")) {
    return(invisible(NULL))
  }

  if (inherits(heatmap, "ggplot")) {
    save_plot(heatmap, filename, width, height, overwrite = overwrite)
    return(invisible(NULL))
  }

  pdf(output_paths[2],
    width = width,
    height = height
  )
  ComplexHeatmap::draw(heatmap)
  dev.off()

  png(output_paths[1],
    width = width,
    height = height,
    units = "in",
    res = 300
  )
  ComplexHeatmap::draw(heatmap)
  dev.off()
}

save_rds_if_missing <- function(object, path, label = "RDS", overwrite = FALSE) {
  if (!overwrite && skip_existing_output(path, label)) {
    return(invisible(FALSE))
  }
  saveRDS(object, file = path)
  invisible(TRUE)
}

write_csv_if_missing <- function(data, path, label = "CSV", overwrite = FALSE) {
  if (!overwrite && skip_existing_output(path, label)) {
    return(invisible(FALSE))
  }
  readr::write_csv(data, path)
  invisible(TRUE)
}

add_fixed_umap_coordinates <- function(plot) {
  # UMAP panels must keep identical x/y scaling. Patchwork grids need the
  # coordinate system applied to every child plot, not just the last panel.
  if (inherits(plot, "patchwork")) {
    return(plot & coord_fixed() & theme_umap_arrows())
  }
  plot + coord_fixed() + theme_umap_arrows()
}

make_umap_plot <- function(seurat_obj, group_by, cols = NULL, ncol = NULL) {
  umap_plot <- DimPlot2(seurat_obj,
    group.by = group_by,
    reduction = "umap.unintegrated",
    theme = NoAxes(),
    label = TRUE,
    box = TRUE,
    label.color = "black",
    repel = TRUE,
    cols = cols,
    ncol = ncol
  )
  add_fixed_umap_coordinates(umap_plot)
}
