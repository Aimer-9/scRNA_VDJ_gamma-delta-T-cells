# R version 4.5.2 (2025-10-31)
rm(list = ls())
setwd("/path/to/project")
library(Seurat)
library(SeuratExtend)
library(tidyverse)
library(ggpubr)
library(patchwork)
library(circlize)
library(ComplexHeatmap)
library(grid)
options(
  max.print = 100,
  tibble.width = Inf,
  spe = "human"
)

rds_dir <- "rds"
figure_dir <- file.path("figures", "8_Vd1vs2")
table_dir <- "table"

seurat_celltype_rds <- file.path(rds_dir, "all_seurat_celltype.rds")
hallmark_stats_csv <- file.path(table_dir, "vd1_vd2_hallmark_celltype_stats.csv")

waterfall_top_n <- 8
hallmark_heatmap_width <- 8
hallmark_heatmap_height <- 5
waterfall_width <- 10
waterfall_height <- 6
force_vd1_vs_vd2 <- FALSE
force_vd1_vs_vd2_plot <- TRUE

vd1_vd2_comparisons <- list(
  naive_vd1_vs_effector_memory_vd2 = c("Naive Vd1", "Effector Memory Vd2"),
  effector_vd1_vs_naive_vd1 = c("Effector Vd1", "Naive Vd1"),
  effector_vd1_vs_zol_effector_vd2 = c("Effector Vd1", "ZOL Effector Vd2"),
  effector_vd1_vs_pan_effector_vd2 = c("Effector Vd1", "PAN Effector Vd2"),
  pan_effector_vd2_vs_zol_effector_vd2 = c("PAN Effector Vd2", "ZOL Effector Vd2"),
  pan_effector_vd2_vs_effector_memory_vd2 = c("PAN Effector Vd2", "Effector Memory Vd2"),
  zol_effector_vd2_vs_effector_memory_vd2 = c("ZOL Effector Vd2", "Effector Memory Vd2")
)

# Load shared palettes plus common IO, metadata, assay, and plotting helpers.
source_plotting_shared <- function() {
  candidates <- c(
    "code/downstream/lib/plotting_shared.R",
    "lib/plotting_shared.R",
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
dir.create(table_dir, showWarnings = FALSE)

run_hallmark_scoring <- function(seurat_obj) {
  GeneSetAnalysis(seurat_obj, genesets = hall50$human)
}

get_hallmark_scores <- function(seurat_obj) {
  scores <- seurat_obj@misc$AUCell$genesets
  if (is.null(scores)) {
    stop("Missing AUCell hallmark scores after GeneSetAnalysis.", call. = FALSE)
  }
  scores
}

make_hallmark_celltype_heatmap <- function(scores, seurat_obj) {
  celltype_stats <- CalcStats(scores, f = seurat_obj$cell_type)
  write_csv_if_missing(
    as.data.frame(celltype_stats) %>% rownames_to_column("pathway"),
    hallmark_stats_csv,
    "Vd1/Vd2 Hallmark cell-type stats CSV",
    overwrite = force_vd1_vs_vd2
  )
  SeuratExtend::Heatmap(celltype_stats, lab_fill = "zscore")
}

get_present_comparisons <- function(seurat_obj, comparisons) {
  present_cell_types <- unique(as.character(seurat_obj$cell_type))
  keep <- vapply(
    comparisons,
    function(pair) all(pair %in% present_cell_types),
    logical(1)
  )
  missing_comparisons <- names(comparisons)[!keep]
  if (length(missing_comparisons) > 0) {
    message("Skipping missing cell-type comparison(s): ", paste(missing_comparisons, collapse = ", "))
  }
  comparisons[keep]
}

call_seuratextend_waterfall <- function(seurat_obj, scores, comparison, comparison_name) {
  if (!exists("WaterfallPlot", envir = asNamespace("SeuratExtend"), inherits = FALSE)) {
    stop("SeuratExtend::WaterfallPlot() is not available in the loaded SeuratExtend package.", call. = FALSE)
  }

  metadata <- seurat_obj@meta.data %>%
    rownames_to_column("cell_id") %>%
    select(cell_id, cell_type) %>%
    mutate(cell_type = as.character(cell_type))
  common_cells <- intersect(colnames(scores), metadata$cell_id)
  if (length(common_cells) == 0) {
    stop("No overlap between Hallmark score columns and Seurat cell IDs.", call. = FALSE)
  }

  metadata <- metadata %>%
    filter(cell_id %in% common_cells) %>%
    arrange(match(cell_id, common_cells))
  scores <- scores[, metadata$cell_id, drop = FALSE]

  SeuratExtend::WaterfallPlot(
    scores,
    f = metadata$cell_type,
    ident.1 = comparison[1],
    ident.2 = comparison[2],
    top.n = waterfall_top_n,
    style = "bar",
    length = "tscore",
    color = "p",
    title = paste(comparison[1], "vs", comparison[2])
  )
}

save_waterfall_plots <- function(scores, seurat_obj, comparisons) {
  present_comparisons <- get_present_comparisons(seurat_obj, comparisons)
  if (length(present_comparisons) == 0) {
    stop("None of the configured Vd1/Vd2 comparisons are present in cell_type metadata.", call. = FALSE)
  }

  plots <- imap(present_comparisons, function(comparison, comparison_name) {
    plot <- call_seuratextend_waterfall(
      seurat_obj = seurat_obj,
      scores = scores,
      comparison = comparison,
      comparison_name = comparison_name
    )
    save_plot(plot, paste0("waterfall_", comparison_name), waterfall_width, waterfall_height, overwrite = force_vd1_vs_vd2_plot)
    plot
  })

  combined_plot <- wrap_plots(plots, ncol = 1)
  save_plot(
    combined_plot,
    "waterfall_vd1_vd2_hallmark_comparisons",
    waterfall_width,
    waterfall_height * length(plots),
    overwrite = force_vd1_vs_vd2_plot
  )
  plots
}

# 1. Read annotated Seurat object and normalize metadata levels.
all_seurat_celltype <- read_rds_checked(seurat_celltype_rds, "annotated Seurat RDS")
all_seurat_celltype <- normalise_metadata_levels(all_seurat_celltype)

# 2. Score Hallmark pathways and summarize by current Vd1/Vd2 cell types.
all_seurat_celltype <- run_hallmark_scoring(all_seurat_celltype)
hallmark_scores <- get_hallmark_scores(all_seurat_celltype)
hallmark_celltype_heatmap <- make_hallmark_celltype_heatmap(
  hallmark_scores,
  all_seurat_celltype
)
save_heatmap(
  hallmark_celltype_heatmap,
  "vd1_vd2_hallmark_celltype_heatmap",
  hallmark_heatmap_width,
  hallmark_heatmap_height,
  overwrite = force_vd1_vs_vd2_plot
)

# 3. Plot selected Vd1 versus Vd2 waterfall comparisons with SeuratExtend.
waterfall_plots <- save_waterfall_plots(
  hallmark_scores,
  all_seurat_celltype,
  vd1_vd2_comparisons
)
