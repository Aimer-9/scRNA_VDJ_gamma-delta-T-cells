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

# Retained Vd1/Vd2 panels consolidated from former steps 14 and 16.
retained_pairwise_comparisons <- list(
  naive_vd1_vs_effector_memory_vd2 = c("Naive Vd1", "Effector Memory Vd2"),
  effector_vd1_vs_pan_effector_vd2 = c("Effector Vd1", "PAN Effector Vd2")
)
retained_vd1_vd2_cell_types <- cell_type_levels
retained_marker_genes <- c(
  "TRDV1", "TRDV2", "TRGV9", "CCR7", "SELL", "IL7R", "TCF7", "LEF1",
  "LTB", "MAL", "BACH2", "SOX4", "NKG7", "GNLY", "PRF1", "GZMA",
  "GZMB", "GZMH", "IFNG", "TNF", "CCL3", "CCL4", "CX3CR1", "FGFBP2",
  "TBX21", "EOMES", "ZEB2", "XBP1", "CEBPB", "CD40LG", "CD70", "ICOS",
  "IL2RA", "PDCD1", "CTLA4", "LAG3", "TIGIT", "HAVCR2", "CD80", "CD86"
)
retained_module_gene_sets <- list(
  Naive_memory = c("CCR7", "SELL", "IL7R", "TCF7", "LEF1", "LTB"),
  Cytotoxicity = c("NKG7", "GNLY", "PRF1", "GZMA", "GZMB", "GZMH", "KLRD1"),
  Inflammation_activation = c("IFNG", "TNF", "CCL3", "CCL4", "CD40LG", "CD70", "IL2RA", "ICOS"),
  Checkpoint = c("PDCD1", "CTLA4", "LAG3", "TIGIT", "HAVCR2"),
  Tissue_migration = c("CX3CR1", "CXCR3", "KLRG1", "FGFBP2")
)

subset_retained_vd1_vd2 <- function(seurat_obj) {
  cells <- rownames(seurat_obj@meta.data)[as.character(seurat_obj$cell_type) %in% retained_vd1_vd2_cell_types]
  if (!length(cells)) stop("No configured Vd1/Vd2 cell types found in `cell_type` metadata.", call. = FALSE)
  subset(seurat_obj, cells = cells)
}

plot_retained_pairwise_marker_dotplot <- function(seurat_obj) {
  cell_types <- unique(unlist(get_present_comparisons(seurat_obj, retained_pairwise_comparisons)))
  plot_obj <- subset(seurat_obj, cells = rownames(seurat_obj@meta.data)[as.character(seurat_obj$cell_type) %in% cell_types])
  genes <- present_genes(plot_obj, retained_marker_genes)
  if (!length(genes)) return(plot_empty("No marker genes available in object"))
  SeuratExtend::DotPlot2(plot_obj, features = genes, group.by = "cell_type", show_grid = FALSE, flip = TRUE) +
    theme(axis.title = element_blank())
}

plot_retained_trdv_violin <- function(seurat_obj) {
  plot_obj <- subset_retained_vd1_vd2(seurat_obj)
  genes <- present_genes(plot_obj, c("TRDV1", "TRDV2"))
  if (!length(genes)) return(plot_empty("TRDV1/TRDV2 genes not found in object"))
  FetchData(plot_obj, vars = genes) %>% rownames_to_column("cell_id") %>% as_tibble() %>%
    left_join(plot_obj@meta.data %>% rownames_to_column("cell_id") %>% select(cell_id, cell_type), by = "cell_id") %>%
    pivot_longer(cols = all_of(genes), names_to = "gene", values_to = "expression") %>%
    mutate(cell_type = factor(as.character(cell_type), levels = retained_vd1_vd2_cell_types)) %>%
    ggplot(aes(cell_type, expression, fill = cell_type)) +
    geom_violin(scale = "width", trim = TRUE) + geom_boxplot(width = 0.08, outlier.shape = NA) +
    facet_wrap(~gene, scales = "free_y", ncol = 1) + scale_fill_manual(values = color_celltype, guide = "none") +
    theme_test() + theme(axis.title.x = element_blank(), axis.text.x = element_text(angle = 35, hjust = 1))
}

add_retained_module_scores <- function(seurat_obj) {
  for (module_name in names(retained_module_gene_sets)) {
    genes <- present_genes(seurat_obj, retained_module_gene_sets[[module_name]])
    score_col <- paste0(module_name, "_score")
    if (!length(genes)) {
      seurat_obj@meta.data[[score_col]] <- NA_real_
    } else {
      seurat_obj <- AddModuleScore(seurat_obj, features = list(genes), name = score_col)
      seurat_obj@meta.data[[score_col]] <- seurat_obj@meta.data[[paste0(score_col, "1")]]
    }
  }
  seurat_obj
}

plot_retained_module_violin <- function(seurat_obj) {
  scores <- seurat_obj@meta.data %>% rownames_to_column("cell_id") %>% as_tibble() %>%
    select(cell_id, cell_type, all_of(paste0(names(retained_module_gene_sets), "_score"))) %>%
    pivot_longer(cols = -c(cell_id, cell_type), names_to = "module", values_to = "score") %>%
    mutate(module = str_remove(module, "_score$"), cell_type = factor(as.character(cell_type), levels = retained_vd1_vd2_cell_types))
  ggplot(scores, aes(cell_type, score, fill = cell_type)) +
    geom_violin(scale = "width", trim = TRUE) + geom_boxplot(width = 0.08, outlier.shape = NA) +
    facet_wrap(~module, scales = "free_y", ncol = 2) + scale_fill_manual(values = color_celltype, guide = "none") +
    theme_test() + theme(axis.title.x = element_blank(), axis.text.x = element_text(angle = 35, hjust = 1)) +
    labs(y = "Module score", title = "Functional module scores across Vd1/Vd2 states")
}

retained_vd1_vd2_obj <- subset_retained_vd1_vd2(all_seurat_celltype)
save_plot(plot_retained_pairwise_marker_dotplot(all_seurat_celltype), "vd1_vd2_pairwise_marker_dotplot", 10, 8, overwrite = force_vd1_vs_vd2_plot)
save_plot(SeuratExtend::DotPlot2(retained_vd1_vd2_obj, features = present_genes(retained_vd1_vd2_obj, retained_marker_genes), group.by = "cell_type", show_grid = FALSE, flip = TRUE), "vd1_vd2_curated_marker_dotplot", 12, 9, overwrite = force_vd1_vs_vd2_plot)
save_plot(plot_retained_trdv_violin(all_seurat_celltype), "vd1_vd2_TRDV1_TRDV2_expression_violin", 10, 6, overwrite = force_vd1_vs_vd2_plot)
retained_vd1_vd2_obj <- add_retained_module_scores(retained_vd1_vd2_obj)
save_plot(plot_retained_module_violin(retained_vd1_vd2_obj), "vd1_vd2_module_score_violin", 13, 8, overwrite = force_vd1_vs_vd2_plot)

# Retained ZOL/PAN effector-Vd2 module-score violin consolidated from former
# step 13.
zol_pan_module_cell_types <- c("ZOL Effector Vd2", "ZOL FOXP3+ Vd2", "PAN Effector Vd2")
zol_pan_module_comparisons <- list(
  c("ZOL Effector Vd2", "ZOL FOXP3+ Vd2"),
  c("ZOL Effector Vd2", "PAN Effector Vd2"),
  c("ZOL FOXP3+ Vd2", "PAN Effector Vd2")
)
zol_pan_module_gene_sets <- list(
  Cytotoxicity = c("NKG7", "GNLY", "PRF1", "GZMA", "GZMB", "GZMH", "KLRD1"),
  Cytokine_inflammatory = c("IFNG", "TNF", "CCL3", "CCL4", "CXCR3"),
  Costimulation = c("CD40LG", "CD70", "ICOS", "IL2RA", "CD80", "CD86"),
  Checkpoint = c("PDCD1", "CTLA4", "LAG3", "TIGIT", "HAVCR2"),
  Proliferation = c("MKI67", "TOP2A", "STMN1", "TYMS", "HMGB2"),
  Antigen_presentation = c("HLA-DRA", "HLA-DRB1", "HLA-DPA1", "HLA-DPB1", "CD74"),
  Stress_AP1 = c("JUN", "JUNB", "FOS", "FOSB", "DUSP1", "IER2")
)

add_zol_pan_module_scores <- function(seurat_obj) {
  for (module_name in names(zol_pan_module_gene_sets)) {
    genes <- present_genes(seurat_obj, zol_pan_module_gene_sets[[module_name]])
    score_col <- paste0("zol_pan_", module_name, "_score")
    if (!length(genes)) {
      seurat_obj@meta.data[[score_col]] <- NA_real_
    } else {
      seurat_obj <- AddModuleScore(seurat_obj, features = list(genes), name = score_col)
      seurat_obj@meta.data[[score_col]] <- seurat_obj@meta.data[[paste0(score_col, "1")]]
    }
  }
  seurat_obj
}

plot_zol_pan_module_violin <- function(seurat_obj) {
  score_columns <- paste0("zol_pan_", names(zol_pan_module_gene_sets), "_score")
  scores <- seurat_obj@meta.data %>% rownames_to_column("cell_id") %>% as_tibble() %>%
    select(cell_id, cell_type, all_of(score_columns)) %>%
    pivot_longer(cols = all_of(score_columns), names_to = "module", values_to = "score") %>%
    mutate(
      module = str_remove(str_remove(module, "^zol_pan_"), "_score$"),
      cell_type = factor(as.character(cell_type), levels = zol_pan_module_cell_types)
    )
  comparisons <- Filter(function(pair) all(pair %in% unique(as.character(scores$cell_type))), zol_pan_module_comparisons)
  ggplot(scores, aes(cell_type, score, fill = cell_type)) +
    geom_violin(scale = "width", trim = TRUE) + geom_boxplot(width = 0.08, outlier.shape = NA) +
    stat_compare_means(comparisons = comparisons, method = "wilcox.test", label = "p.signif", hide.ns = FALSE, size = 3) +
    facet_wrap(~module, scales = "free_y", ncol = 2) + scale_fill_manual(values = color_celltype, guide = "none") +
    theme_test() + theme(axis.title.x = element_blank(), axis.text.x = element_text(angle = 25, hjust = 1)) +
    labs(y = "Module score", title = "ZOL/PAN effector Vd2 functional module scores")
}

zol_pan_cells <- rownames(all_seurat_celltype@meta.data)[as.character(all_seurat_celltype$cell_type) %in% zol_pan_module_cell_types]
if (length(zol_pan_cells)) {
  zol_pan_obj <- add_zol_pan_module_scores(subset(all_seurat_celltype, cells = zol_pan_cells))
  save_plot(plot_zol_pan_module_violin(zol_pan_obj), "zol_pan_effector_vd2_module_score_violin", 12, 7, overwrite = force_vd1_vs_vd2_plot)
} else {
  message("[SKIP] No configured ZOL/PAN Vd2 states found for module-score violin.")
}
