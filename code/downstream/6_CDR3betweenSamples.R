# R version 4.5.2 (2025-10-31)
rm(list = ls())
setwd("/path/to/project")
library(Seurat)
library(SeuratExtend)
library(tidyverse)
library(gplots)
library(ggpubr)
library(patchwork)
library(circlize)
library(ComplexHeatmap)
library(grid)
options(
  tibble.width = Inf,
  print.width = Inf,
  max.print = 200
)

rds_dir <- "rds"
figure_dir <- file.path("figures", "6_CDR3betweenSamples")

seurat_celltype_rds <- file.path(rds_dir, "all_seurat_celltype.rds")
included_annotation_rds <- file.path(rds_dir, "all_annotation_included.rds")

top_cdr3_n <- 100
shared_heatmap_gap_cells <- 2
shared_heatmap_width_margin_in <- 5
shared_heatmap_height_margin_in <- 3
force_cdr3_between_samples <- FALSE
force_cdr3_between_samples_plot <- TRUE

source_plotting_shared <- function() {
  candidates <- c(
    "code/downstream/_cache_/plotting_shared.R",
    "_cache_/plotting_shared.R",
    "cache/plotting_shared.R"
  )
  plotting_shared <- candidates[file.exists(candidates)][1]
  if (is.na(plotting_shared)) {
    stop("Missing plotting_shared.R. Checked: ", paste(candidates, collapse = ", "))
  }
  source(plotting_shared)
}
source_plotting_shared()

dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

normalise_metadata_levels <- function(seurat_obj, annotation) {
  seurat_obj$group <- factor(seurat_obj$group, levels = group_levels)
  seurat_obj$sample_name <- factor(seurat_obj$sample_name, levels = sample_name_levels)
  seurat_obj$cell_type <- factor(seurat_obj$cell_type, levels = cell_type_levels)

  annotation <- annotation %>%
    mutate(
      group = factor(group, levels = group_levels),
      sample_name = factor(sample_name, levels = sample_name_levels),
      cell_type = factor(cell_type, levels = cell_type_levels)
    )

  list(seurat = seurat_obj, annotation = annotation)
}

get_present_levels <- function(values, levels) {
  levels[levels %in% as.character(values)]
}

summarise_top_cdr3 <- function(annotation, group_column, group_levels, top_n = top_cdr3_n) {
  annotation %>%
    group_by(.data[[group_column]], chain, cdr3) %>%
    summarise(n = n(), .groups = "drop") %>%
    group_by(.data[[group_column]], chain) %>%
    arrange(desc(n), .by_group = TRUE) %>%
    slice_head(n = top_n) %>%
    ungroup() %>%
    mutate(
      group_value = factor(
        as.character(.data[[group_column]]),
        levels = get_present_levels(annotation[[group_column]], group_levels)
      )
    )
}

count_shared_cdr3 <- function(top_cdr3, group_values, chain_name, group_a, group_b) {
  set_a <- top_cdr3 %>%
    filter(chain == chain_name, group_value == group_a) %>%
    pull(cdr3)
  set_b <- top_cdr3 %>%
    filter(chain == chain_name, group_value == group_b) %>%
    pull(cdr3)
  length(intersect(set_a, set_b))
}

make_shared_cdr3_table <- function(annotation, group_column, group_levels) {
  group_values <- get_present_levels(annotation[[group_column]], group_levels)
  top_cdr3 <- summarise_top_cdr3(annotation, group_column, group_levels)

  expand_grid(Var1 = group_values, Var2 = group_values) %>%
    mutate(
      inter_trd = map2_int(
        Var1,
        Var2,
        ~ count_shared_cdr3(top_cdr3, group_values, "TRD", .x, .y)
      ),
      inter_trg = map2_int(
        Var1,
        Var2,
        ~ count_shared_cdr3(top_cdr3, group_values, "TRG", .x, .y)
      ),
      Var1 = factor(Var1, levels = group_values),
      Var2 = factor(Var2, levels = group_values)
    )
}

make_shared_matrix <- function(shared_table, value_column) {
  shared_table %>%
    select(Var1, Var2, all_of(value_column)) %>%
    pivot_wider(names_from = Var2, values_from = all_of(value_column)) %>%
    column_to_rownames("Var1") %>%
    as.matrix()
}

make_heatmap_colors <- function(matrix) {
  colorRamp2(
    c(0, max(1, max(matrix, na.rm = TRUE))),
    hcl_palette = "Blues 3",
    reverse = TRUE
  )
}

make_shared_cdr3_heatmap <- function(matrix,
                                     chain_label,
                                     triangle,
                                     cell_size_mm,
                                     cell_number_font_size,
                                     row_name_font_size,
                                     title_font_size) {
  ComplexHeatmap::Heatmap(
    matrix,
    rect_gp = gpar(type = "none"),
    width = ncol(matrix) * unit(cell_size_mm, "mm"),
    height = nrow(matrix) * unit(cell_size_mm, "mm"),
    show_row_names = triangle == "lower",
    show_column_names = FALSE,
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    row_names_side = "left",
    row_names_gp = gpar(fontsize = row_name_font_size),
    show_heatmap_legend = FALSE,
    column_title = paste("Number of shared", chain_label, "in top100"),
    column_title_side = ifelse(triangle == "lower", "bottom", "top"),
    column_title_gp = gpar(fontsize = title_font_size),
    cell_fun = function(j, i, x, y, w, h, fill) {
      should_draw <- if (triangle == "lower") i >= j else i <= j
      if (should_draw) {
        grid.rect(x, y, w, h, gp = gpar(fill = fill, col = fill))
        grid.text(sprintf("%d", matrix[i, j]),
          x,
          y,
          gp = gpar(fontsize = cell_number_font_size)
        )
      }
    },
    col = make_heatmap_colors(matrix)
  )
}

save_heatmap_pair <- function(trd_heatmap, trg_heatmap, filename, width, height, ht_gap_mm, overwrite = FALSE) {
  output_paths <- file.path(figure_dir, paste0(filename, c(".png", ".pdf")))
  if (output_files_exist(output_paths) && !overwrite) {
    message("[SKIP] Existing heatmap pair: ", paste(output_paths, collapse = ", "))
    return(invisible(NULL))
  } else if (output_files_exist(output_paths) && overwrite) {
    message("[FORCE] Rebuilding heatmap pair: ", paste(output_paths, collapse = ", "))
  }

  pdf(output_paths[2],
    width = width,
    height = height
  )
  ComplexHeatmap::draw(
    trd_heatmap + trg_heatmap,
    ht_gap = unit(ht_gap_mm, "mm")
  )
  dev.off()

  png(output_paths[1],
    width = width,
    height = height,
    units = "in",
    res = 300
  )
  ComplexHeatmap::draw(
    trd_heatmap + trg_heatmap,
    ht_gap = unit(ht_gap_mm, "mm")
  )
  dev.off()
}

get_heatmap_pair_gap_mm <- function(matrix, cell_size_mm, gap_cells = shared_heatmap_gap_cells) {
  gap_cells <- min(gap_cells, ncol(matrix))
  -(ncol(matrix) - gap_cells) * cell_size_mm
}

get_heatmap_pair_size_in <- function(matrix,
                                     cell_size_mm,
                                     gap_cells = shared_heatmap_gap_cells,
                                     width_margin_in = shared_heatmap_width_margin_in,
                                     height_margin_in = shared_heatmap_height_margin_in) {
  gap_cells <- min(gap_cells, ncol(matrix))
  effective_width_mm <- (ncol(matrix) + gap_cells) * cell_size_mm
  height_mm <- nrow(matrix) * cell_size_mm

  list(
    width = effective_width_mm / 25.4 + width_margin_in,
    height = height_mm / 25.4 + height_margin_in
  )
}

run_shared_cdr3_heatmaps <- function(annotation,
                                     group_column,
                                     group_levels,
                                     filename,
                                     cell_size_mm,
                                     width = NULL,
                                     height = NULL,
                                     ht_gap_mm = NULL,
                                     cell_number_font_size = 10,
                                     row_name_font_size = 10,
                                     title_font_size = 12) {
  shared_table <- make_shared_cdr3_table(annotation, group_column, group_levels)
  trd_matrix <- make_shared_matrix(shared_table, "inter_trd")
  trg_matrix <- make_shared_matrix(shared_table, "inter_trg")

  trd_heatmap <- make_shared_cdr3_heatmap(
    trd_matrix,
    "CDR3d",
    triangle = "lower",
    cell_size_mm = cell_size_mm,
    cell_number_font_size = cell_number_font_size,
    row_name_font_size = row_name_font_size,
    title_font_size = title_font_size
  )
  trg_heatmap <- make_shared_cdr3_heatmap(
    trg_matrix,
    "CDR3g",
    triangle = "upper",
    cell_size_mm = cell_size_mm,
    cell_number_font_size = cell_number_font_size,
    row_name_font_size = row_name_font_size,
    title_font_size = title_font_size
  )

  if (is.null(ht_gap_mm)) {
    ht_gap_mm <- get_heatmap_pair_gap_mm(trd_matrix, cell_size_mm)
  }
  if (is.null(width) || is.null(height)) {
    plot_size <- get_heatmap_pair_size_in(trd_matrix, cell_size_mm)
    if (is.null(width)) {
      width <- plot_size$width
    }
    if (is.null(height)) {
      height <- plot_size$height
    }
  }

  save_heatmap_pair(
    trd_heatmap,
    trg_heatmap,
    filename,
    width,
    height,
    ht_gap_mm,
    overwrite = force_cdr3_between_samples_plot
  )

  list(
    shared_table = shared_table,
    trd_matrix = trd_matrix,
    trg_matrix = trg_matrix
  )
}

# 1. Read final cell-type object and included productive VDJ annotations.
all_seurat_celltype <- readRDS(seurat_celltype_rds)
all_annotation_included <- readRDS(included_annotation_rds)
normalised <- normalise_metadata_levels(all_seurat_celltype, all_annotation_included)
all_seurat_celltype <- normalised$seurat
all_annotation_included <- normalised$annotation
rm(normalised)

# 2. Plot shared top100 CDR3s between cell types.
shared_cdr3_cell_type <- run_shared_cdr3_heatmaps(
  all_annotation_included,
  group_column = "cell_type",
  group_levels = cell_type_levels,
  filename = "shared_cdr3_between_cell_type",
  cell_size_mm = 30,
  cell_number_font_size = 30,
  row_name_font_size = 30,
  title_font_size = 50
)


# 3. Plot shared top100 CDR3s between samples.
shared_cdr3_samples <- run_shared_cdr3_heatmaps(
  all_annotation_included,
  group_column = "sample_name",
  group_levels = sample_name_levels,
  filename = "shared_cdr3_between_samples",
  cell_size_mm = 20,
  cell_number_font_size = 30,
  row_name_font_size = 30,
  title_font_size = 50
)
