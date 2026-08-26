# R version 4.5.2 (2025-10-31)
rm(list = ls())

project_dir <- "/path/to/project"
if (dir.exists(project_dir)) {
  setwd(project_dir)
}

library(Seurat)
library(SeuratExtend)
library(tidyverse)
library(patchwork)
options(
  tibble.width = Inf,
  print.width = Inf,
  max.print = 200
)

rds_dir <- "rds"
figure_dir <- file.path("figures", "15_CD80_CD86_expression")

seurat_celltype_rds <- file.path(rds_dir, "all_seurat_celltype.rds")
genes_to_plot <- c("CD80", "CD86")
force_cd80_cd86_plot <- TRUE

# Load shared palettes plus common IO, metadata, assay, and plotting helpers.
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

plot_cd80_cd86_umap <- function(seurat_obj, features) {
  reduction_name <- get_umap_reduction(seurat_obj)
  umap_plot <- DimPlot2(
    seurat_obj,
    features = features,
    reduction = reduction_name,
    theme = NoAxes(),
    ncol = length(features)
  )
  add_fixed_umap_coordinates(umap_plot)
}

plot_cd80_cd86_group_dotplot <- function(seurat_obj, features) {
  DotPlot2(
    seurat_obj,
    features = features,
    group.by = "group",
    cols = c("lightgrey", "#B2182B"),
    show_grid = FALSE
  ) +
    xlab(NULL) +
    ylab(NULL) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

# 1. Read the final annotated Seurat object.
all_seurat_celltype <- read_rds_checked(seurat_celltype_rds, "annotated Seurat RDS")
DefaultAssay(all_seurat_celltype) <- "RNA"
all_seurat_celltype <- join_assay_layers_if_needed(all_seurat_celltype, assay = "RNA")
all_seurat_celltype$group <- factor(all_seurat_celltype$group, levels = group_levels)
all_seurat_celltype$sample_name <- factor(all_seurat_celltype$sample_name, levels = sample_name_levels)
all_seurat_celltype$cell_type <- factor(all_seurat_celltype$cell_type, levels = cell_type_levels)
genes_to_plot <- check_features(all_seurat_celltype, genes_to_plot)

# 2. Plot CD80 and CD86 expression on UMAP.
cd80_cd86_umap <- plot_cd80_cd86_umap(all_seurat_celltype, genes_to_plot)
save_plot(
  cd80_cd86_umap,
  "CD80_CD86_expression_umap",
  10,
  5,
  overwrite = force_cd80_cd86_plot
)

# 3. Plot CD80 and CD86 expression by sample group.
cd80_cd86_group_dotplot <- plot_cd80_cd86_group_dotplot(all_seurat_celltype, genes_to_plot)
save_plot(
  cd80_cd86_group_dotplot,
  "CD80_CD86_expression_group_dotplot",
  6,
  4,
  overwrite = force_cd80_cd86_plot
)
