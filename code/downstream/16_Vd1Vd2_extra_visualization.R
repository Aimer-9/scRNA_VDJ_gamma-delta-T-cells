# R version 4.5.2 (2025-10-31)
rm(list = ls())
setwd("/path/to/project")

library(Seurat)
library(SeuratExtend)
library(tidyverse)
library(patchwork)
options(
  tibble.width = Inf,
  print.width = Inf,
  max.print = 200,
  spe = "human"
)

rds_dir <- "rds"
figure_dir <- file.path("figures", "16_Vd1Vd2_extra_visualization")
table_dir <- "table"

seurat_celltype_rds <- file.path(rds_dir, "all_seurat_celltype.rds")
included_annotation_rds <- file.path(rds_dir, "all_annotation_included.rds")
paired_cdr3_rds <- file.path(rds_dir, "barcode_trgd_paired.rds")

cell_fraction_csv <- file.path(table_dir, "vd1_vd2_extra_cell_fraction_by_sample.csv")
marker_summary_csv <- file.path(table_dir, "vd1_vd2_extra_marker_summary.csv")
module_score_summary_csv <- file.path(table_dir, "vd1_vd2_extra_module_score_summary.csv")
de_markers_csv <- file.path(table_dir, "vd1_vd2_extra_de_markers.csv")
clone_size_csv <- file.path(table_dir, "vd1_vd2_extra_clone_size_by_celltype.csv")
top_cdr3_csv <- file.path(table_dir, "vd1_vd2_extra_top_cdr3_by_celltype.csv")
paired_clone_summary_csv <- file.path(table_dir, "vd1_vd2_extra_paired_clone_summary.csv")
paired_clone_sharing_csv <- file.path(table_dir, "vd1_vd2_extra_trdg_pair_celltype_sharing.csv")
trdv_chain_expression_csv <- file.path(table_dir, "vd1_vd2_TRDV1_TRDV2_expression_by_trd_chain_type.csv")

force_vd1_vd2_extra <- FALSE
force_vd1_vd2_extra_plot <- TRUE

vd1_vd2_cell_types <- c(
  "Naive Vd1",
  "Pre-activated Vd1",
  "Effector Vd1",
  "Effector Memory Vd2",
  "Pre-activated Vd2",
  "ZOL Effector Vd2",
  "ZOL FOXP3+ Vd2",
  "PAN Effector Vd2"
)

trd_chain_type_levels <- c("TRDV1", "TRDV2", "TRDV1+TRDV2", "Other/unknown")

vd1_vd2_comparisons <- list(
  naive_vd1_vs_effector_memory_vd2 = c("Naive Vd1", "Effector Memory Vd2"),
  effector_vd1_vs_naive_vd1 = c("Effector Vd1", "Naive Vd1"),
  effector_vd1_vs_zol_effector_vd2 = c("Effector Vd1", "ZOL Effector Vd2"),
  effector_vd1_vs_pan_effector_vd2 = c("Effector Vd1", "PAN Effector Vd2"),
  preactivated_vd1_vs_preactivated_vd2 = c("Pre-activated Vd1", "Pre-activated Vd2"),
  pan_effector_vd2_vs_zol_effector_vd2 = c("PAN Effector Vd2", "ZOL Effector Vd2"),
  pan_effector_vd2_vs_effector_memory_vd2 = c("PAN Effector Vd2", "Effector Memory Vd2"),
  zol_effector_vd2_vs_effector_memory_vd2 = c("ZOL Effector Vd2", "Effector Memory Vd2")
)

marker_gene_classes <- tribble(
  ~gene_class, ~gene,
  "TCR identity", "TRDV1",
  "TCR identity", "TRDV2",
  "TCR identity", "TRGV9",
  "Naive memory", "CCR7",
  "Naive memory", "SELL",
  "Naive memory", "IL7R",
  "Naive memory", "TCF7",
  "Naive memory", "LEF1",
  "Naive memory", "LTB",
  "Naive memory", "MAL",
  "Developmental transcription", "BACH2",
  "Developmental transcription", "SOX4",
  "Cytotoxicity", "NKG7",
  "Cytotoxicity", "GNLY",
  "Cytotoxicity", "PRF1",
  "Cytotoxicity", "GZMA",
  "Cytotoxicity", "GZMB",
  "Cytotoxicity", "GZMH",
  "Inflammatory cytokine", "IFNG",
  "Inflammatory cytokine", "TNF",
  "Inflammatory cytokine", "CCL3",
  "Inflammatory cytokine", "CCL4",
  "Migration NK effector", "KLRD1",
  "Migration NK effector", "KLRG1",
  "Migration NK effector", "CX3CR1",
  "Migration NK effector", "FGFBP2",
  "Effector transcription", "TBX21",
  "Effector transcription", "EOMES",
  "Effector transcription", "ZEB2",
  "Effector transcription", "XBP1",
  "Effector transcription", "CEBPB",
  "Costimulation activation", "CD40LG",
  "Costimulation activation", "CD70",
  "Costimulation activation", "ICOS",
  "Costimulation activation", "IL2RA",
  "Checkpoint exhaustion", "PDCD1",
  "Checkpoint exhaustion", "CTLA4",
  "Checkpoint exhaustion", "LAG3",
  "Checkpoint exhaustion", "TIGIT",
  "Checkpoint exhaustion", "HAVCR2",
  "APC costimulation", "CD80",
  "APC costimulation", "CD86"
)
marker_gene_classes$gene_class <- factor(
  marker_gene_classes$gene_class,
  levels = unique(marker_gene_classes$gene_class)
)
marker_gene_panel <- unique(marker_gene_classes$gene)

module_gene_sets <- list(
  Naive_memory = c("CCR7", "SELL", "IL7R", "TCF7", "LEF1", "LTB"),
  Cytotoxicity = c("NKG7", "GNLY", "PRF1", "GZMA", "GZMB", "GZMH", "KLRD1"),
  Inflammation_activation = c("IFNG", "TNF", "CCL3", "CCL4", "CD40LG", "CD70", "IL2RA", "ICOS"),
  Checkpoint = c("PDCD1", "CTLA4", "LAG3", "TIGIT", "HAVCR2"),
  Tissue_migration = c("CX3CR1", "CXCR3", "KLRG1", "FGFBP2")
)

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
dir.create(table_dir, showWarnings = FALSE)

get_vd1_vd2_cells <- function(seurat_obj) {
  rownames(seurat_obj@meta.data)[as.character(seurat_obj$cell_type) %in% vd1_vd2_cell_types]
}

subset_vd1_vd2 <- function(seurat_obj) {
  selected_cells <- get_vd1_vd2_cells(seurat_obj)
  if (length(selected_cells) == 0) {
    stop("No configured Vd1/Vd2 cell types found in `cell_type` metadata.", call. = FALSE)
  }
  subset(seurat_obj, cells = selected_cells)
}

get_present_comparisons <- function(seurat_obj) {
  present_cell_types <- unique(as.character(seurat_obj$cell_type))
  keep <- vapply(vd1_vd2_comparisons, function(pair) all(pair %in% present_cell_types), logical(1))
  if (any(!keep)) {
    message("[SKIP] Missing Vd1/Vd2 comparison(s): ", paste(names(vd1_vd2_comparisons)[!keep], collapse = ", "))
  }
  vd1_vd2_comparisons[keep]
}

make_vd1_vd2_umap_highlight <- function(seurat_obj) {
  plot_obj <- seurat_obj
  plot_obj$vd1_vd2_cell_type <- ifelse(
    as.character(plot_obj$cell_type) %in% vd1_vd2_cell_types,
    as.character(plot_obj$cell_type),
    "Other"
  )
  plot_obj$vd1_vd2_cell_type <- factor(plot_obj$vd1_vd2_cell_type, levels = c(vd1_vd2_cell_types, "Other"))
  plot_colors <- c(color_celltype[vd1_vd2_cell_types], "Other" = "grey88")
  umap_plot <- DimPlot2(
    plot_obj,
    group.by = "vd1_vd2_cell_type",
    reduction = get_umap_reduction(plot_obj),
    theme = NoAxes(),
    cols = plot_colors,
    label = TRUE,
    box = TRUE,
    label.color = "black",
    repel = TRUE
  )
  add_fixed_umap_coordinates(umap_plot) +
    labs(title = "Vd1 and Vd2 cell states on UMAP", color = "Cell type")
}

make_vd1_vd2_subset_umap <- function(seurat_obj) {
  plot_obj <- subset_vd1_vd2(seurat_obj)
  plot_obj$cell_type <- factor(as.character(plot_obj$cell_type), levels = vd1_vd2_cell_types)
  umap_plot <- DimPlot2(
    plot_obj,
    group.by = "cell_type",
    reduction = get_umap_reduction(plot_obj),
    theme = NoAxes(),
    cols = color_celltype[vd1_vd2_cell_types],
    label = TRUE,
    box = TRUE,
    label.color = "black",
    repel = TRUE
  )
  add_fixed_umap_coordinates(umap_plot) +
    labs(title = "Vd1+Vd2 cells only", color = "Cell type")
}

calculate_cell_fraction_by_sample <- function(seurat_obj) {
  metadata <- seurat_obj@meta.data %>%
    as_tibble() %>%
    filter(!is.na(group), !is.na(sample_name), !is.na(cell_type)) %>%
    mutate(
      cell_type = factor(as.character(cell_type), levels = vd1_vd2_cell_types),
      sample_name = factor(as.character(sample_name), levels = sample_name_levels),
      group = factor(as.character(group), levels = group_levels)
    )
  sample_totals <- metadata %>%
    dplyr::count(group, sample_name, name = "sample_total_cells")
  metadata %>%
    filter(!is.na(cell_type)) %>%
    dplyr::count(group, sample_name, cell_type, name = "n_cells") %>%
    left_join(sample_totals, by = c("group", "sample_name")) %>%
    mutate(percent_of_sample = n_cells / sample_total_cells * 100) %>%
    arrange(group, sample_name, cell_type)
}

plot_cell_fraction_stacked <- function(cell_fraction) {
  if (nrow(cell_fraction) == 0) {
    return(plot_empty("No Vd1/Vd2 cell fraction data"))
  }
  ggplot(cell_fraction, aes(x = sample_name, y = percent_of_sample, fill = cell_type)) +
    geom_col(width = 0.75, color = "white", linewidth = 0.2) +
    scale_fill_manual(values = color_celltype, drop = FALSE) +
    theme_test() +
    theme(axis.title.x = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(y = "Cells in sample (%)", fill = "Cell type", title = "Vd1/Vd2 composition by sample")
}

plot_cell_fraction_points <- function(cell_fraction) {
  if (nrow(cell_fraction) == 0) {
    return(plot_empty("No Vd1/Vd2 cell fraction data"))
  }
  ggplot(cell_fraction, aes(x = group, y = percent_of_sample, color = cell_type)) +
    geom_boxplot(outlier.shape = NA, width = 0.55, alpha = 0.25) +
    geom_point(aes(shape = sample_name), position = position_jitter(width = 0.12), size = 2, alpha = 0.9) +
    facet_wrap(~cell_type, scales = "free_y", ncol = 4) +
    scale_color_manual(values = color_celltype, drop = FALSE, guide = "none") +
    theme_test() +
    theme(axis.title.x = element_blank(), axis.text.x = element_text(angle = 35, hjust = 1), strip.background = element_blank()) +
    labs(y = "Cells in sample (%)", shape = "Sample", title = "Sample-level Vd1/Vd2 abundance")
}

make_marker_summary <- function(seurat_obj, genes) {
  if (length(genes) == 0) {
    return(tibble())
  }
  expression_data <- FetchData(seurat_obj, vars = genes) %>%
    rownames_to_column("cell_id") %>%
    as_tibble() %>%
    left_join(
      seurat_obj@meta.data %>%
        rownames_to_column("cell_id") %>%
        select(cell_id, group, sample_name, cell_type),
      by = "cell_id"
    ) %>%
    filter(as.character(cell_type) %in% vd1_vd2_cell_types) %>%
    pivot_longer(cols = all_of(genes), names_to = "gene", values_to = "expression")
  expression_data %>%
    group_by(gene, cell_type) %>%
    dplyr::summarise(
      mean_expression = mean(expression, na.rm = TRUE),
      median_expression = median(expression, na.rm = TRUE),
      percent_expressing = mean(expression > 0, na.rm = TRUE) * 100,
      n_cells = n(),
      .groups = "drop"
    ) %>%
    left_join(marker_gene_classes, by = "gene") %>%
    mutate(
      gene_class = ifelse(is.na(gene_class), "Other", as.character(gene_class)),
      gene_class = factor(gene_class, levels = c(levels(marker_gene_classes$gene_class), "Other"))
    )
}

plot_marker_dotplot <- function(seurat_obj, genes) {
  plot_genes <- present_genes(seurat_obj, genes)
  if (length(plot_genes) == 0) {
    return(plot_empty("No marker genes available in object"))
  }
  plot_obj <- subset_vd1_vd2(seurat_obj)
  expression_data <- FetchData(plot_obj, vars = plot_genes) %>%
    rownames_to_column("cell_id") %>%
    as_tibble() %>%
    left_join(
      plot_obj@meta.data %>%
        rownames_to_column("cell_id") %>%
        select(cell_id, cell_type),
      by = "cell_id"
    ) %>%
    pivot_longer(cols = all_of(plot_genes), names_to = "gene", values_to = "expression")
  plot_data <- expression_data %>%
    group_by(cell_type, gene) %>%
    dplyr::summarise(
      mean_expression = mean(expression, na.rm = TRUE),
      percent_expressing = mean(expression > 0, na.rm = TRUE) * 100,
      .groups = "drop"
    ) %>%
    group_by(gene) %>%
    mutate(
      scaled_expression = as.numeric(scale(mean_expression)),
      scaled_expression = ifelse(is.na(scaled_expression), 0, scaled_expression)
    ) %>%
    ungroup() %>%
    left_join(marker_gene_classes, by = "gene") %>%
    mutate(
      gene_class = ifelse(is.na(gene_class), "Other", as.character(gene_class)),
      gene_class = factor(gene_class, levels = c(levels(marker_gene_classes$gene_class), "Other")),
      cell_type = factor(as.character(cell_type), levels = vd1_vd2_cell_types),
      gene = factor(gene, levels = rev(marker_gene_classes$gene[marker_gene_classes$gene %in% plot_genes]))
    )
  ggplot(plot_data, aes(x = cell_type, y = gene)) +
    geom_point(aes(size = percent_expressing, color = scaled_expression), alpha = 0.9) +
    facet_grid(gene_class ~ ., scales = "free_y", space = "free_y") +
    scale_color_gradient2(low = "#2166AC", mid = "lightgrey", high = "#B2182B", midpoint = 0, name = "Scaled\nmean") +
    scale_size(range = c(0.4, 4.2), name = "% expressing") +
    theme_test() +
    theme(
      axis.title = element_blank(),
      axis.text.x = element_text(angle = 35, hjust = 1),
      strip.background = element_blank(),
      panel.spacing.y = grid::unit(0.08, "lines")
    ) +
    labs(title = "Curated Vd1/Vd2 marker program by gene class")
}

plot_trdv_expression_violin <- function(seurat_obj, genes = c("TRDV1", "TRDV2")) {
  plot_genes <- present_genes(seurat_obj, genes)
  if (length(plot_genes) == 0) {
    return(plot_empty("TRDV1/TRDV2 genes not found in object"))
  }
  plot_obj <- subset_vd1_vd2(seurat_obj)
  expression_data <- FetchData(plot_obj, vars = plot_genes) %>%
    rownames_to_column("cell_id") %>%
    as_tibble() %>%
    left_join(
      plot_obj@meta.data %>%
        rownames_to_column("cell_id") %>%
        select(cell_id, cell_type),
      by = "cell_id"
    ) %>%
    pivot_longer(cols = all_of(plot_genes), names_to = "gene", values_to = "expression") %>%
    mutate(
      cell_type = factor(as.character(cell_type), levels = vd1_vd2_cell_types),
      gene = factor(gene, levels = genes[genes %in% plot_genes])
    )
  ggplot(expression_data, aes(x = cell_type, y = expression, fill = cell_type)) +
    geom_violin(scale = "width", trim = TRUE, width = 0.72, linewidth = 0.18) +
    geom_boxplot(width = 0.08, outlier.shape = NA, color = "grey25", alpha = 0.7, linewidth = 0.18) +
    facet_wrap(~gene, scales = "free_y", ncol = 1, drop = FALSE) +
    scale_fill_manual(values = color_celltype, drop = FALSE, guide = "none") +
    scale_x_discrete(drop = FALSE, expand = expansion(mult = c(0.02, 0.02))) +
    coord_cartesian(clip = "off") +
    theme_test() +
    theme(
      axis.title.x = element_blank(),
      axis.text.x = element_text(angle = 35, hjust = 1, vjust = 1),
      strip.background = element_blank(),
      strip.text = element_text(margin = margin(1.5, 0, 1.5, 0)),
      panel.spacing.y = grid::unit(0.25, "lines"),
      plot.margin = margin(4, 6, 4, 4)
    ) +
    labs(y = "Expression", title = "TRDV1 and TRDV2 expression across Vd1/Vd2 states")
}

infer_trd_chain_type <- function(annotation) {
  if (is.null(annotation) || nrow(annotation) == 0 || !"barcode" %in% colnames(annotation)) {
    return(tibble())
  }
  has_vd1 <- "Vd1" %in% colnames(annotation)
  has_vd2 <- "Vd2" %in% colnames(annotation)
  if (!has_vd1 && !has_vd2) {
    return(tibble())
  }
  annotation %>%
    mutate(
      has_TRDV1 = if (has_vd1) !is.na(Vd1) & as.character(Vd1) != "" & as.character(Vd1) != "0" else FALSE,
      has_TRDV2 = if (has_vd2) !is.na(Vd2) & as.character(Vd2) != "" & as.character(Vd2) != "0" else FALSE,
      trd_chain_type = case_when(
        has_TRDV1 & has_TRDV2 ~ "TRDV1+TRDV2",
        has_TRDV1 ~ "TRDV1",
        has_TRDV2 ~ "TRDV2",
        TRUE ~ "Other/unknown"
      )
    ) %>%
    distinct(barcode, trd_chain_type) %>%
    mutate(trd_chain_type = factor(trd_chain_type, levels = trd_chain_type_levels))
}

make_trdv_expression_by_chain_type <- function(seurat_obj, annotation, genes = c("TRDV1", "TRDV2")) {
  plot_genes <- present_genes(seurat_obj, genes)
  chain_type <- infer_trd_chain_type(annotation)
  if (length(plot_genes) == 0 || nrow(chain_type) == 0) {
    return(tibble())
  }
  plot_obj <- subset_vd1_vd2(seurat_obj)
  FetchData(plot_obj, vars = plot_genes) %>%
    rownames_to_column("cell_id") %>%
    as_tibble() %>%
    left_join(
      plot_obj@meta.data %>%
        rownames_to_column("cell_id") %>%
        select(cell_id, barcode, group, sample_name, cell_type),
      by = "cell_id"
    ) %>%
    inner_join(chain_type, by = "barcode") %>%
    pivot_longer(cols = all_of(plot_genes), names_to = "gene", values_to = "expression") %>%
    mutate(
      gene = factor(gene, levels = genes[genes %in% plot_genes]),
      trd_chain_type = factor(as.character(trd_chain_type), levels = trd_chain_type_levels),
      cell_type = factor(as.character(cell_type), levels = vd1_vd2_cell_types),
      group = factor(as.character(group), levels = group_levels),
      sample_name = factor(as.character(sample_name), levels = sample_name_levels)
    )
}

plot_trdv_expression_by_chain_type <- function(expression_data) {
  if (nrow(expression_data) == 0) {
    return(plot_empty("No TRDV expression and TRD chain type data available"))
  }
  ggplot(expression_data, aes(x = trd_chain_type, y = expression, fill = trd_chain_type)) +
    geom_violin(scale = "width", trim = TRUE, width = 0.72, linewidth = 0.18) +
    geom_boxplot(width = 0.08, outlier.shape = NA, color = "grey25", alpha = 0.7, linewidth = 0.18) +
    facet_wrap(~gene, scales = "free_y", ncol = 1, drop = FALSE) +
    scale_fill_manual(
      values = c(
        "TRDV1" = "#1B9E77",
        "TRDV2" = "#D95F02",
        "TRDV1+TRDV2" = "#7570B3",
        "Other/unknown" = "grey70"
      ),
      drop = FALSE,
      guide = "none"
    ) +
    scale_x_discrete(drop = FALSE, expand = expansion(mult = c(0.02, 0.02))) +
    coord_cartesian(clip = "off") +
    theme_test() +
    theme(
      axis.title.x = element_blank(),
      axis.text.x = element_text(angle = 25, hjust = 1, vjust = 1),
      strip.background = element_blank(),
      strip.text = element_text(margin = margin(1.5, 0, 1.5, 0)),
      panel.spacing.y = grid::unit(0.25, "lines"),
      plot.margin = margin(4, 6, 4, 4)
    ) +
    labs(y = "Expression", title = "TRDV1/TRDV2 expression by called TRD chain type")
}

plot_marker_heatmap <- function(marker_summary) {
  if (nrow(marker_summary) == 0) {
    return(plot_empty("No marker expression summary available"))
  }
  plot_data <- marker_summary %>%
    group_by(gene) %>%
    mutate(
      scaled_expression = as.numeric(scale(mean_expression)),
      scaled_expression = ifelse(is.na(scaled_expression), 0, scaled_expression)
    ) %>%
    ungroup() %>%
    mutate(
      cell_type = factor(as.character(cell_type), levels = vd1_vd2_cell_types),
      gene_class = ifelse(is.na(gene_class), "Other", as.character(gene_class)),
      gene_class = factor(gene_class, levels = c(levels(marker_gene_classes$gene_class), "Other")),
      gene = factor(gene, levels = rev(unique(gene)))
    )
  ggplot(plot_data, aes(x = cell_type, y = gene, fill = scaled_expression)) +
    geom_tile(color = "white", linewidth = 0.25) +
    facet_grid(gene_class ~ ., scales = "free_y", space = "free_y") +
    scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0) +
    theme_test() +
    theme(axis.title = element_blank(), axis.text.x = element_text(angle = 35, hjust = 1), strip.background = element_blank()) +
    labs(fill = "Row-scaled\nmean", title = "Average marker expression by Vd1/Vd2 state")
}

add_module_scores <- function(seurat_obj, gene_sets) {
  for (module_name in names(gene_sets)) {
    module_genes <- present_genes(seurat_obj, gene_sets[[module_name]])
    score_col <- paste0(module_name, "_score")
    if (length(module_genes) == 0) {
      seurat_obj[[score_col]] <- NA_real_
      message("[SKIP] No genes found for module score: ", module_name)
      next
    }
    seurat_obj <- AddModuleScore(
      seurat_obj,
      features = list(module_genes),
      name = score_col
    )
    seurat_obj@meta.data[[score_col]] <- seurat_obj@meta.data[[paste0(score_col, "1")]]
  }
  seurat_obj
}

make_module_score_long <- function(seurat_obj) {
  score_columns <- paste0(names(module_gene_sets), "_score")
  seurat_obj@meta.data %>%
    rownames_to_column("cell_id") %>%
    as_tibble() %>%
    filter(as.character(cell_type) %in% vd1_vd2_cell_types) %>%
    select(cell_id, group, sample_name, cell_type, all_of(score_columns)) %>%
    pivot_longer(cols = all_of(score_columns), names_to = "module", values_to = "score") %>%
    mutate(
      module = str_replace(module, "_score$", ""),
      cell_type = factor(as.character(cell_type), levels = vd1_vd2_cell_types),
      group = factor(as.character(group), levels = group_levels)
    )
}

summarise_module_scores <- function(module_scores) {
  if (nrow(module_scores) == 0) {
    return(tibble())
  }
  module_scores %>%
    group_by(module, group, sample_name, cell_type) %>%
    dplyr::summarise(
      n_cells = n(),
      mean_score = mean(score, na.rm = TRUE),
      median_score = median(score, na.rm = TRUE),
      .groups = "drop"
    )
}

plot_module_score_violin <- function(module_scores) {
  if (nrow(module_scores) == 0) {
    return(plot_empty("No module score data available"))
  }
  plot_data <- module_scores %>%
    mutate(
      cell_type = factor(as.character(cell_type), levels = vd1_vd2_cell_types),
      module = factor(as.character(module), levels = names(module_gene_sets))
    )
  ggplot(plot_data, aes(x = cell_type, y = score, fill = cell_type)) +
    geom_violin(scale = "width", trim = TRUE, width = 0.72, linewidth = 0.18) +
    geom_boxplot(width = 0.08, outlier.shape = NA, color = "grey25", alpha = 0.75, linewidth = 0.18) +
    facet_wrap(~module, scales = "free_y", ncol = 2, drop = FALSE) +
    scale_fill_manual(values = color_celltype, drop = FALSE, guide = "none") +
    scale_x_discrete(drop = FALSE, expand = expansion(mult = c(0.02, 0.02))) +
    coord_cartesian(clip = "off") +
    theme_test() +
    theme(
      axis.title.x = element_blank(),
      axis.text.x = element_text(angle = 35, hjust = 1, vjust = 1),
      strip.background = element_blank(),
      strip.text = element_text(margin = margin(1.5, 0, 1.5, 0)),
      panel.spacing.x = grid::unit(0.25, "lines"),
      panel.spacing.y = grid::unit(0.25, "lines"),
      plot.margin = margin(4, 6, 4, 4)
    ) +
    labs(y = "Module score", title = "Functional module scores across Vd1/Vd2 states")
}

plot_module_score_group_dotplot <- function(module_summary) {
  if (nrow(module_summary) == 0) {
    return(plot_empty("No module score summary available"))
  }
  plot_data <- module_summary %>%
    group_by(module, group, cell_type) %>%
    dplyr::summarise(mean_score = mean(mean_score, na.rm = TRUE), .groups = "drop") %>%
    mutate(
      group_cell_type = interaction(group, cell_type, sep = " | ", drop = TRUE),
      group_cell_type = factor(group_cell_type, levels = unique(group_cell_type))
    )
  ggplot(plot_data, aes(x = module, y = group_cell_type, color = mean_score, size = abs(mean_score))) +
    geom_point(alpha = 0.9) +
    scale_color_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0) +
    scale_size(range = c(1.2, 5.5), guide = "none") +
    theme_test() +
    theme(axis.title = element_blank(), axis.text.x = element_text(angle = 35, hjust = 1)) +
    labs(color = "Mean score", title = "Functional module score summary by group and cell state")
}

run_pairwise_de <- function(seurat_obj, comparisons) {
  if (length(comparisons) == 0) {
    return(tibble())
  }
  imap_dfr(comparisons, function(comparison, comparison_name) {
    selected_cells <- rownames(seurat_obj@meta.data)[as.character(seurat_obj$cell_type) %in% comparison]
    comparison_obj <- subset(seurat_obj, cells = selected_cells)
    Idents(comparison_obj) <- factor(as.character(comparison_obj$cell_type), levels = comparison)
    if (min(table(Idents(comparison_obj))) < 3) {
      message("[SKIP] Too few cells for DE: ", comparison_name)
      return(tibble())
    }
    FindMarkers(
      comparison_obj,
      ident.1 = comparison[1],
      ident.2 = comparison[2],
      only.pos = FALSE,
      logfc.threshold = 0,
      min.pct = 0.05
    ) %>%
      rownames_to_column("gene") %>%
      as_tibble() %>%
      mutate(
        comparison = comparison_name,
        ident_1 = comparison[1],
        ident_2 = comparison[2],
        higher_in = ifelse(avg_log2FC >= 0, ident_1, ident_2),
        abs_log2FC = abs(avg_log2FC),
        .before = 1
      )
  })
}

plot_de_volcano <- function(de_markers, comparison_name, top_n = 12) {
  plot_data <- de_markers %>%
    filter(comparison == comparison_name) %>%
    mutate(
      neg_log10_padj = -log10(pmax(p_val_adj, .Machine$double.xmin)),
      significant = p_val_adj < 0.05 & abs_log2FC >= 0.25
    )
  if (nrow(plot_data) == 0) {
    return(NULL)
  }
  label_data <- plot_data %>%
    filter(significant) %>%
    slice_max(order_by = abs_log2FC, n = top_n, with_ties = FALSE)
  volcano_plot <- ggplot(plot_data, aes(x = avg_log2FC, y = neg_log10_padj, color = higher_in)) +
    geom_point(aes(alpha = significant), size = 0.9) +
    geom_vline(xintercept = 0, color = "grey45", linewidth = 0.35) +
    scale_alpha_manual(values = c(`TRUE` = 0.85, `FALSE` = 0.25), guide = "none") +
    scale_color_manual(values = color_celltype, drop = FALSE) +
    theme_test() +
    labs(x = "avg_log2FC", y = "-log10 adjusted P", color = "Higher in", title = comparison_name)
  if (nrow(label_data) > 0 && requireNamespace("ggrepel", quietly = TRUE)) {
    volcano_plot <- volcano_plot +
      ggrepel::geom_text_repel(
        data = label_data,
        aes(label = gene),
        size = 3,
        max.overlaps = Inf,
        show.legend = FALSE
      )
  } else if (nrow(label_data) > 0) {
    volcano_plot <- volcano_plot +
      geom_text(data = label_data, aes(label = gene), size = 3, vjust = -0.5, show.legend = FALSE)
  }
  volcano_plot
}

plot_de_overlap_matrix <- function(de_markers, top_n = 40) {
  if (nrow(de_markers) == 0) {
    return(plot_empty("No DE marker data available"))
  }
  sig_genes <- de_markers %>%
    filter(!is.na(p_val_adj), p_val_adj < 0.05, abs_log2FC >= 0.25) %>%
    group_by(gene) %>%
    dplyr::summarise(
      comparison_n = n_distinct(comparison),
      max_abs_log2FC = max(abs_log2FC, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(desc(comparison_n), desc(max_abs_log2FC), gene) %>%
    slice_head(n = top_n)
  if (nrow(sig_genes) == 0) {
    return(plot_empty("No significant DE genes for overlap matrix"))
  }
  plot_data <- expand_grid(
    gene = sig_genes$gene,
    comparison = unique(de_markers$comparison)
  ) %>%
    left_join(
      de_markers %>%
        filter(gene %in% sig_genes$gene) %>%
        select(gene, comparison, avg_log2FC, p_val_adj),
      by = c("gene", "comparison")
    ) %>%
    mutate(
      significant = !is.na(p_val_adj) & p_val_adj < 0.05 & abs(avg_log2FC) >= 0.25,
      signed_presence = case_when(
        significant & avg_log2FC > 0 ~ 1,
        significant & avg_log2FC < 0 ~ -1,
        TRUE ~ 0
      ),
      gene = factor(gene, levels = rev(sig_genes$gene)),
      comparison = factor(comparison, levels = unique(de_markers$comparison))
    )
  ggplot(plot_data, aes(x = comparison, y = gene, fill = signed_presence)) +
    geom_tile(color = "white", linewidth = 0.2) +
    scale_fill_gradient2(
      low = "#2166AC",
      mid = "white",
      high = "#B2182B",
      midpoint = 0,
      breaks = c(-1, 0, 1),
      labels = c("Lower in ident.1", "Not significant", "Higher in ident.1")
    ) +
    theme_test() +
    theme(axis.title = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(fill = "DE direction", title = "Shared differential genes across Vd1/Vd2 comparisons")
}

plot_top_de_heatmap <- function(de_markers, seurat_obj, top_n_per_comparison = 8) {
  if (nrow(de_markers) == 0) {
    return(plot_empty("No DE marker data available"))
  }
  top_genes <- de_markers %>%
    filter(!is.na(p_val_adj), p_val_adj < 0.05) %>%
    group_by(comparison) %>%
    slice_max(order_by = abs_log2FC, n = top_n_per_comparison, with_ties = FALSE) %>%
    ungroup() %>%
    pull(gene) %>%
    unique()
  top_genes <- present_genes(seurat_obj, top_genes)
  if (length(top_genes) == 0) {
    return(plot_empty("No significant DE genes found in object"))
  }
  marker_summary <- make_marker_summary(seurat_obj, top_genes)
  plot_marker_heatmap(marker_summary) +
    labs(title = "Top DE gene expression across Vd1/Vd2 states")
}

calculate_clone_size_by_celltype <- function(annotation) {
  if (is.null(annotation) || nrow(annotation) == 0) {
    return(tibble())
  }
  annotation %>%
    filter(as.character(cell_type) %in% vd1_vd2_cell_types, !is.na(cdr3), cdr3 != "") %>%
    distinct(barcode, cell_type, chain, cdr3) %>%
    dplyr::count(cell_type, chain, cdr3, name = "clone_cells") %>%
    group_by(cell_type, chain) %>%
    mutate(
      total_cells = sum(clone_cells),
      clone_fraction = clone_cells / total_cells,
      clone_size_class = case_when(
        clone_fraction < 0.001 ~ "Rare (<0.1%)",
        clone_fraction < 0.01 ~ "Large (0.1-1%)",
        clone_fraction < 0.1 ~ "Expanded (1-10%)",
        TRUE ~ "Hyperexpanded (>10%)"
      ),
      clone_size_class = factor(
        clone_size_class,
        levels = c("Rare (<0.1%)", "Large (0.1-1%)", "Expanded (1-10%)", "Hyperexpanded (>10%)")
      )
    ) %>%
    ungroup()
}

plot_clone_size_stacked <- function(clone_size) {
  if (nrow(clone_size) == 0) {
    return(plot_empty("No clone-size data available"))
  }
  clone_size_colors <- c(
    "Rare (<0.1%)" = "#D9D9D9",
    "Large (0.1-1%)" = "#9ECAE1",
    "Expanded (1-10%)" = "#3182BD",
    "Hyperexpanded (>10%)" = "#08519C"
  )
  plot_data <- clone_size %>%
    group_by(cell_type, chain, clone_size_class) %>%
    dplyr::summarise(clone_cells = sum(clone_cells), .groups = "drop_last") %>%
    mutate(percent = clone_cells / sum(clone_cells) * 100) %>%
    ungroup() %>%
    mutate(cell_type = factor(as.character(cell_type), levels = vd1_vd2_cell_types))
  ggplot(plot_data, aes(x = cell_type, y = percent, fill = clone_size_class)) +
    geom_col(width = 0.72, color = "white", linewidth = 0.2) +
    facet_wrap(~chain, ncol = 1) +
    scale_fill_manual(values = clone_size_colors, drop = FALSE) +
    theme_test() +
    theme(axis.title.x = element_blank(), axis.text.x = element_text(angle = 35, hjust = 1), strip.background = element_blank()) +
    labs(y = "Cells in clone-size class (%)", fill = "Clone-size class", title = "TCR clone-size classes across Vd1/Vd2 states")
}

calculate_top_cdr3_by_celltype <- function(annotation, top_n = 10) {
  if (is.null(annotation) || nrow(annotation) == 0) {
    return(tibble())
  }
  annotation %>%
    filter(as.character(cell_type) %in% vd1_vd2_cell_types, !is.na(cdr3), cdr3 != "") %>%
    distinct(barcode, cell_type, chain, cdr3) %>%
    dplyr::count(cell_type, chain, cdr3, name = "n_cells") %>%
    group_by(cell_type, chain) %>%
    mutate(percent = n_cells / sum(n_cells) * 100) %>%
    slice_max(order_by = n_cells, n = top_n, with_ties = FALSE) %>%
    arrange(cell_type, chain, desc(n_cells)) %>%
    ungroup()
}

plot_top_cdr3_composition <- function(top_cdr3) {
  if (nrow(top_cdr3) == 0) {
    return(plot_empty("No top CDR3 data available"))
  }
  plot_data <- top_cdr3 %>%
    mutate(
      cell_type = factor(as.character(cell_type), levels = vd1_vd2_cell_types),
      cdr3_label = str_trunc(cdr3, width = 24)
    )
  ggplot(plot_data, aes(x = reorder(cdr3_label, percent), y = percent, fill = cell_type)) +
    geom_col(width = 0.75) +
    coord_flip() +
    facet_grid(chain ~ cell_type, scales = "free_y", space = "free_y") +
    scale_fill_manual(values = color_celltype, drop = FALSE, guide = "none") +
    theme_test() +
    theme(axis.title.y = element_blank(), strip.background = element_blank()) +
    labs(y = "Cells with CDR3 (%)", title = "Top CDR3 clones by Vd1/Vd2 state")
}

prepare_paired_clone_metadata <- function(paired_cdr3, seurat_obj) {
  if (is.null(paired_cdr3) || nrow(paired_cdr3) == 0) {
    return(tibble())
  }
  cell_metadata <- seurat_obj@meta.data %>%
    rownames_to_column("cell_id") %>%
    as_tibble() %>%
    select(cell_id, barcode, group, sample_name, cell_type) %>%
    filter(as.character(cell_type) %in% vd1_vd2_cell_types)
  paired_cdr3 %>%
    select(barcode, TRD, TRG) %>%
    distinct() %>%
    inner_join(cell_metadata, by = "barcode") %>%
    filter(!is.na(TRD), !is.na(TRG), TRD != "", TRG != "") %>%
    mutate(
      TRDG_pair = paste(TRD, TRG, sep = "||"),
      cell_type = factor(as.character(cell_type), levels = vd1_vd2_cell_types),
      group = factor(as.character(group), levels = group_levels)
    ) %>%
    distinct()
}

summarise_paired_clones <- function(paired_metadata, top_n = 12) {
  if (nrow(paired_metadata) == 0) {
    return(tibble())
  }
  top_pairs <- paired_metadata %>%
    dplyr::count(TRDG_pair, TRD, TRG, name = "total_cells") %>%
    arrange(desc(total_cells), TRDG_pair) %>%
    slice_head(n = top_n) %>%
    mutate(pair_rank = paste0("pair_", row_number()))
  paired_metadata %>%
    inner_join(top_pairs, by = c("TRDG_pair", "TRD", "TRG")) %>%
    dplyr::count(pair_rank, TRDG_pair, TRD, TRG, group, sample_name, cell_type, total_cells, name = "n_cells") %>%
    group_by(pair_rank) %>%
    mutate(percent_of_pair = n_cells / sum(n_cells) * 100) %>%
    ungroup()
}

plot_clone_alluvial <- function(paired_clone_summary) {
  if (nrow(paired_clone_summary) == 0) {
    return(plot_empty("No paired TRD+TRG clone data available"))
  }
  if (!requireNamespace("ggalluvial", quietly = TRUE)) {
    return(plot_empty("ggalluvial is not installed", "Install ggalluvial to draw group -> cell type -> paired clone alluvial plots."))
  }
  plot_data <- paired_clone_summary %>%
    mutate(
      clone_label = factor(pair_rank, levels = paste0("pair_", seq_len(length(unique(pair_rank))))),
      cell_type = factor(as.character(cell_type), levels = vd1_vd2_cell_types)
    )
  ggplot(plot_data, aes(y = n_cells, axis1 = group, axis2 = cell_type, axis3 = clone_label)) +
    ggalluvial::geom_alluvium(aes(fill = cell_type), alpha = 0.65, width = 1 / 12) +
    ggalluvial::geom_stratum(width = 1 / 6, fill = "grey95", color = "grey55") +
    ggalluvial::stat_stratum(geom = "text", aes(label = after_stat(stratum)), size = 3) +
    scale_x_discrete(limits = c("Group", "Cell type", "Top TRD+TRG pair"), expand = c(0.08, 0.08)) +
    scale_fill_manual(values = color_celltype, drop = FALSE) +
    theme_test() +
    theme(axis.title = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
    labs(fill = "Cell type", title = "Group to Vd1/Vd2 state to top paired TRD+TRG clone")
}

calculate_paired_clone_sharing <- function(paired_metadata) {
  if (nrow(paired_metadata) == 0) {
    return(tibble())
  }
  pair_celltype <- paired_metadata %>%
    distinct(TRDG_pair, cell_type) %>%
    dplyr::count(TRDG_pair, name = "celltype_n") %>%
    filter(celltype_n >= 2)
  paired_metadata %>%
    semi_join(pair_celltype, by = "TRDG_pair") %>%
    dplyr::count(cell_type, TRDG_pair, name = "n_cells") %>%
    group_by(TRDG_pair) %>%
    mutate(total_pair_cells = sum(n_cells)) %>%
    ungroup() %>%
    arrange(desc(total_pair_cells), TRDG_pair, cell_type)
}

plot_paired_clone_sharing_heatmap <- function(sharing_table, top_n = 40) {
  if (nrow(sharing_table) == 0) {
    return(plot_empty("No shared paired TRD+TRG clones across Vd1/Vd2 states"))
  }
  top_pairs <- sharing_table %>%
    distinct(TRDG_pair, total_pair_cells) %>%
    arrange(desc(total_pair_cells), TRDG_pair) %>%
    slice_head(n = top_n) %>%
    pull(TRDG_pair)
  plot_data <- sharing_table %>%
    filter(TRDG_pair %in% top_pairs) %>%
    mutate(
      pair_label = factor(str_trunc(TRDG_pair, width = 36), levels = rev(str_trunc(top_pairs, width = 36))),
      cell_type = factor(as.character(cell_type), levels = vd1_vd2_cell_types)
    )
  ggplot(plot_data, aes(x = cell_type, y = pair_label, fill = n_cells)) +
    geom_tile(color = "white", linewidth = 0.2) +
    scale_fill_gradient(low = "white", high = "#B2182B") +
    theme_test() +
    theme(axis.title = element_blank(), axis.text.x = element_text(angle = 35, hjust = 1)) +
    labs(fill = "Cells", title = "Paired TRD+TRG clone sharing across Vd1/Vd2 states")
}

# 1. Load final object and prepare expression assay.
all_seurat_celltype <- read_rds_checked(seurat_celltype_rds, "annotated Seurat RDS")
all_seurat_celltype <- normalise_metadata_levels(all_seurat_celltype)
all_seurat_celltype <- prepare_expression_object(all_seurat_celltype)

vd1_vd2_obj <- subset_vd1_vd2(all_seurat_celltype)
present_comparisons <- get_present_comparisons(vd1_vd2_obj)
marker_genes <- present_genes(vd1_vd2_obj, marker_gene_panel)

# 2. Identity overview and abundance.
save_plot(
  make_vd1_vd2_umap_highlight(all_seurat_celltype),
  "vd1_vd2_umap_highlight_all_states",
  9,
  7,
  overwrite = force_vd1_vd2_extra_plot
)
save_plot(
  make_vd1_vd2_subset_umap(vd1_vd2_obj),
  "vd1_vd2_umap_subset_only",
  8,
  7,
  overwrite = force_vd1_vd2_extra_plot
)

cell_fraction <- calculate_cell_fraction_by_sample(all_seurat_celltype)
write_csv_if_missing(cell_fraction, cell_fraction_csv, "Vd1/Vd2 cell fraction by sample CSV", overwrite = force_vd1_vd2_extra)
save_plot(plot_cell_fraction_stacked(cell_fraction), "vd1_vd2_cell_fraction_by_sample", 11, 5, overwrite = force_vd1_vd2_extra_plot)
save_plot(plot_cell_fraction_points(cell_fraction), "vd1_vd2_cell_fraction_group_points", 12, 7, overwrite = force_vd1_vd2_extra_plot)

# 3. Curated marker program plots.
marker_summary <- make_marker_summary(vd1_vd2_obj, marker_genes)
write_csv_if_missing(marker_summary, marker_summary_csv, "Vd1/Vd2 marker expression summary CSV", overwrite = force_vd1_vd2_extra)
save_plot(plot_marker_dotplot(vd1_vd2_obj, marker_genes), "vd1_vd2_curated_marker_dotplot", 12, 9, overwrite = force_vd1_vd2_extra_plot)
save_plot(plot_trdv_expression_violin(vd1_vd2_obj), "vd1_vd2_TRDV1_TRDV2_expression_violin", 10, 6, overwrite = force_vd1_vd2_extra_plot)
save_plot(plot_marker_heatmap(marker_summary), "vd1_vd2_curated_marker_heatmap", 10, max(6, length(marker_genes) * 0.22), overwrite = force_vd1_vd2_extra_plot)

# 4. Functional module scores.
vd1_vd2_obj <- add_module_scores(vd1_vd2_obj, module_gene_sets)
module_scores <- make_module_score_long(vd1_vd2_obj)
module_score_summary <- summarise_module_scores(module_scores)
write_csv_if_missing(module_score_summary, module_score_summary_csv, "Vd1/Vd2 module score summary CSV", overwrite = force_vd1_vd2_extra)
save_plot(plot_module_score_violin(module_scores), "vd1_vd2_module_score_violin", 13, 8, overwrite = force_vd1_vd2_extra_plot)
save_plot(plot_module_score_group_dotplot(module_score_summary), "vd1_vd2_module_score_group_dotplot", 11, 8, overwrite = force_vd1_vd2_extra_plot)

# 5. Expanded differential-expression views across all configured comparisons.
de_markers <- run_pairwise_de(vd1_vd2_obj, present_comparisons)
write_csv_if_missing(de_markers, de_markers_csv, "Vd1/Vd2 extra DE markers CSV", overwrite = force_vd1_vd2_extra)
iwalk(present_comparisons, function(comparison, comparison_name) {
  volcano_plot <- plot_de_volcano(de_markers, comparison_name)
  if (!is.null(volcano_plot)) {
    save_plot(volcano_plot, paste0("volcano_", comparison_name), 7, 6, overwrite = force_vd1_vd2_extra_plot)
  }
})
save_plot(plot_de_overlap_matrix(de_markers), "vd1_vd2_de_gene_overlap_matrix", 12, 10, overwrite = force_vd1_vd2_extra_plot)
save_plot(plot_top_de_heatmap(de_markers, vd1_vd2_obj), "vd1_vd2_top_de_gene_heatmap", 10, 10, overwrite = force_vd1_vd2_extra_plot)

# 6. Repertoire-aware plots from productive VDJ annotation and paired TRD+TRG clones.
all_annotation_included <- read_optional_rds(included_annotation_rds, "included VDJ annotation RDS")
if (!is.null(all_annotation_included)) {
  all_annotation_included <- all_annotation_included %>%
    mutate(
      group = factor(as.character(group), levels = group_levels),
      sample_name = factor(as.character(sample_name), levels = sample_name_levels),
      cell_type = factor(as.character(cell_type), levels = cell_type_levels)
    )
  clone_size <- calculate_clone_size_by_celltype(all_annotation_included)
  top_cdr3 <- calculate_top_cdr3_by_celltype(all_annotation_included)
  trdv_chain_expression <- make_trdv_expression_by_chain_type(vd1_vd2_obj, all_annotation_included)
  write_csv_if_missing(clone_size, clone_size_csv, "Vd1/Vd2 clone-size class CSV", overwrite = force_vd1_vd2_extra)
  write_csv_if_missing(top_cdr3, top_cdr3_csv, "Vd1/Vd2 top CDR3 CSV", overwrite = force_vd1_vd2_extra)
  write_csv_if_missing(trdv_chain_expression, trdv_chain_expression_csv, "TRDV expression by TRD chain type CSV", overwrite = force_vd1_vd2_extra)
  save_plot(plot_clone_size_stacked(clone_size), "vd1_vd2_clone_size_stacked_bar", 11, 7, overwrite = force_vd1_vd2_extra_plot)
  save_plot(plot_top_cdr3_composition(top_cdr3), "vd1_vd2_top_cdr3_composition", 14, 10, overwrite = force_vd1_vd2_extra_plot)
  save_plot(plot_trdv_expression_by_chain_type(trdv_chain_expression), "vd1_vd2_TRDV1_TRDV2_expression_by_TRD_chain_type", 8, 6, overwrite = force_vd1_vd2_extra_plot)
}

barcode_trgd_paired <- read_optional_rds(paired_cdr3_rds, "paired TRD/TRG CDR3 RDS")
if (!is.null(barcode_trgd_paired)) {
  paired_metadata <- prepare_paired_clone_metadata(barcode_trgd_paired, all_seurat_celltype)
  paired_clone_summary <- summarise_paired_clones(paired_metadata)
  paired_clone_sharing <- calculate_paired_clone_sharing(paired_metadata)
  write_csv_if_missing(paired_clone_summary, paired_clone_summary_csv, "Vd1/Vd2 paired clone summary CSV", overwrite = force_vd1_vd2_extra)
  write_csv_if_missing(paired_clone_sharing, paired_clone_sharing_csv, "Vd1/Vd2 paired clone sharing CSV", overwrite = force_vd1_vd2_extra)
  save_plot(plot_clone_alluvial(paired_clone_summary), "vd1_vd2_group_celltype_clone_alluvial", 13, 7, overwrite = force_vd1_vd2_extra_plot)
  save_plot(plot_paired_clone_sharing_heatmap(paired_clone_sharing), "vd1_vd2_trdg_pair_sharing_heatmap", 10, 10, overwrite = force_vd1_vd2_extra_plot)
}
