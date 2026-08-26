# R version 4.5.2 (2025-10-31)
rm(list = ls())
setwd("/path/to/project")

library(Seurat)
library(SeuratExtend)
library(tidyverse)
library(patchwork)
library(ggpubr)
options(
  tibble.width = Inf,
  print.width = Inf,
  max.print = 200,
  spe = "human"
)

rds_dir <- "rds"
figure_dir <- file.path("figures", "17_ZOL_PAN_effector_Vd2_comparison")
table_dir <- "table"

seurat_celltype_rds <- file.path(rds_dir, "all_seurat_celltype.rds")
included_annotation_rds <- file.path(rds_dir, "all_annotation_included.rds")
paired_cdr3_rds <- file.path(rds_dir, "barcode_trgd_paired.rds")

cell_fraction_csv <- file.path(table_dir, "zol_pan_effector_vd2_cell_fraction_by_sample.csv")
marker_summary_csv <- file.path(table_dir, "zol_pan_effector_vd2_marker_summary.csv")
module_score_summary_csv <- file.path(table_dir, "zol_pan_effector_vd2_module_score_summary.csv")
de_markers_csv <- file.path(table_dir, "zol_pan_effector_vd2_de_markers.csv")
hallmark_delta_csv <- file.path(table_dir, "zol_pan_effector_vd2_hallmark_delta.csv")
hallmark_stats_csv <- file.path(table_dir, "zol_pan_effector_vd2_hallmark_stats.csv")
clone_size_csv <- file.path(table_dir, "zol_pan_effector_vd2_clone_size.csv")
top_cdr3_csv <- file.path(table_dir, "zol_pan_effector_vd2_top_cdr3.csv")
paired_clone_summary_csv <- file.path(table_dir, "zol_pan_effector_vd2_paired_clone_summary.csv")
paired_clone_overlap_csv <- file.path(table_dir, "zol_pan_effector_vd2_paired_clone_overlap.csv")

force_zol_pan_effector_vd2 <- FALSE
force_zol_pan_effector_vd2_plot <- TRUE

target_cell_types <- c("ZOL Effector Vd2", "PAN Effector Vd2")
comparison_name <- "pan_effector_vd2_vs_zol_effector_vd2"
comparison_pair <- c("PAN Effector Vd2", "ZOL Effector Vd2")

marker_gene_panel <- c(
  "TRDV2", "TRGV9",
  "NKG7", "GNLY", "PRF1", "GZMA", "GZMB", "GZMH",
  "IFNG", "TNF", "CCL3", "CCL4",
  "KLRD1", "KLRG1", "CX3CR1", "FGFBP2",
  "CD40LG", "CD70", "ICOS", "IL2RA",
  "PDCD1", "CTLA4", "LAG3", "TIGIT", "HAVCR2",
  "MKI67", "TOP2A", "STMN1",
  "XBP1", "CEBPB", "SOX4", "TBX21", "EOMES",
  "CD80", "CD86", "HLA-DRA", "HLA-DRB1", "CD74"
)

module_gene_sets <- list(
  Cytotoxicity = c("NKG7", "GNLY", "PRF1", "GZMA", "GZMB", "GZMH", "KLRD1"),
  Cytokine_inflammatory = c("IFNG", "TNF", "CCL3", "CCL4", "CXCR3"),
  Costimulation = c("CD40LG", "CD70", "ICOS", "IL2RA", "CD80", "CD86"),
  Checkpoint = c("PDCD1", "CTLA4", "LAG3", "TIGIT", "HAVCR2"),
  Proliferation = c("MKI67", "TOP2A", "STMN1", "TYMS", "HMGB2"),
  Antigen_presentation = c("HLA-DRA", "HLA-DRB1", "HLA-DPA1", "HLA-DPB1", "CD74"),
  Stress_AP1 = c("JUN", "JUNB", "FOS", "FOSB", "DUSP1", "IER2")
)

selected_hallmark_pathways <- c(
  "HALLMARK_INTERFERON_GAMMA_RESPONSE",
  "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
  "HALLMARK_INFLAMMATORY_RESPONSE",
  "HALLMARK_IL6_JAK_STAT3_SIGNALING",
  "HALLMARK_MTORC1_SIGNALING",
  "HALLMARK_OXIDATIVE_PHOSPHORYLATION",
  "HALLMARK_APOPTOSIS",
  "HALLMARK_DNA_REPAIR"
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

subset_target_cells <- function(seurat_obj) {
  selected_cells <- rownames(seurat_obj@meta.data)[as.character(seurat_obj$cell_type) %in% target_cell_types]
  if (length(selected_cells) == 0) {
    stop("No ZOL/PAN effector Vd2 cells found in `cell_type` metadata.", call. = FALSE)
  }
  subset(seurat_obj, cells = selected_cells)
}

make_umap_highlight <- function(seurat_obj) {
  plot_obj <- seurat_obj
  plot_obj$zol_pan_effector_vd2 <- ifelse(
    as.character(plot_obj$cell_type) %in% target_cell_types,
    as.character(plot_obj$cell_type),
    "Other"
  )
  plot_obj$zol_pan_effector_vd2 <- factor(plot_obj$zol_pan_effector_vd2, levels = c(target_cell_types, "Other"))
  plot_colors <- c(color_celltype[target_cell_types], "Other" = "grey88")
  umap_plot <- DimPlot2(
    plot_obj,
    group.by = "zol_pan_effector_vd2",
    reduction = get_umap_reduction(plot_obj),
    theme = NoAxes(),
    cols = plot_colors,
    label = TRUE,
    box = TRUE,
    label.color = "black",
    repel = TRUE
  )
  add_fixed_umap_coordinates(umap_plot) +
    labs(title = "ZOL and PAN effector Vd2 cells on UMAP", color = "Cell type")
}

make_sample_umap <- function(seurat_obj) {
  plot_obj <- subset_target_cells(seurat_obj)
  umap_plot <- DimPlot2(
    plot_obj,
    group.by = "sample_name",
    reduction = get_umap_reduction(plot_obj),
    theme = NoAxes(),
    cols = color_sample_name
  )
  add_fixed_umap_coordinates(umap_plot) +
    labs(title = "ZOL/PAN effector Vd2 cells colored by sample", color = "Sample")
}

make_umap_density <- function(seurat_obj) {
  reduction_name <- get_umap_reduction(seurat_obj)
  plot_obj <- subset_target_cells(seurat_obj)
  umap_data <- Embeddings(plot_obj, reduction_name) %>%
    as.data.frame() %>%
    rownames_to_column("cell_id") %>%
    as_tibble()
  colnames(umap_data)[2:3] <- c("UMAP_1", "UMAP_2")
  plot_data <- plot_obj@meta.data %>%
    rownames_to_column("cell_id") %>%
    as_tibble() %>%
    left_join(umap_data, by = "cell_id") %>%
    filter(is.finite(UMAP_1), is.finite(UMAP_2))
  ggplot(plot_data, aes(x = UMAP_1, y = UMAP_2, color = cell_type)) +
    geom_point(size = 0.25, alpha = 0.35) +
    geom_density_2d(linewidth = 0.45, alpha = 0.85) +
    facet_wrap(~cell_type) +
    scale_color_manual(values = color_celltype, drop = FALSE, guide = "none") +
    NoAxes() +
    coord_fixed() +
    theme(strip.background = element_blank()) +
    labs(title = "UMAP density of ZOL/PAN effector Vd2 states")
}

calculate_cell_fraction <- function(seurat_obj) {
  metadata <- seurat_obj@meta.data %>%
    as_tibble() %>%
    filter(!is.na(group), !is.na(sample_name), !is.na(cell_type)) %>%
    mutate(
      cell_type = factor(as.character(cell_type), levels = target_cell_types),
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
    return(plot_empty("No ZOL/PAN effector Vd2 abundance data"))
  }
  ggplot(cell_fraction, aes(x = sample_name, y = percent_of_sample, fill = cell_type)) +
    geom_col(width = 0.72, color = "white", linewidth = 0.2) +
    scale_fill_manual(values = color_celltype, drop = FALSE) +
    theme_test() +
    theme(axis.title.x = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(y = "Cells in sample (%)", fill = "Cell type", title = "ZOL/PAN effector Vd2 abundance by sample")
}

plot_cell_fraction_points <- function(cell_fraction) {
  if (nrow(cell_fraction) == 0) {
    return(plot_empty("No ZOL/PAN effector Vd2 abundance data"))
  }
  ggplot(cell_fraction, aes(x = cell_type, y = percent_of_sample, color = cell_type)) +
    geom_boxplot(outlier.shape = NA, width = 0.55) +
    geom_point(aes(shape = sample_name), position = position_jitter(width = 0.08), size = 2, alpha = 0.9) +
    scale_color_manual(values = color_celltype, drop = FALSE, guide = "none") +
    theme_test() +
    theme(axis.title.x = element_blank(), axis.text.x = element_text(angle = 20, hjust = 1)) +
    labs(y = "Cells in sample (%)", shape = "Sample", title = "Sample-level abundance")
}

plot_cell_fraction_ratio <- function(cell_fraction) {
  ratio_data <- cell_fraction %>%
    select(group, sample_name, cell_type, percent_of_sample) %>%
    mutate(cell_type = as.character(cell_type)) %>%
    pivot_wider(names_from = cell_type, values_from = percent_of_sample, values_fill = 0)
  if (!"PAN Effector Vd2" %in% colnames(ratio_data)) {
    ratio_data[["PAN Effector Vd2"]] <- 0
  }
  if (!"ZOL Effector Vd2" %in% colnames(ratio_data)) {
    ratio_data[["ZOL Effector Vd2"]] <- 0
  }
  ratio_data <- ratio_data %>%
    mutate(
      pan_to_zol_ratio = (`PAN Effector Vd2` + 1e-6) / (`ZOL Effector Vd2` + 1e-6),
      log2_pan_to_zol_ratio = log2(pan_to_zol_ratio)
    )
  if (nrow(ratio_data) == 0) {
    return(plot_empty("No ZOL/PAN effector Vd2 ratio data"))
  }
  ggplot(ratio_data, aes(x = group, y = log2_pan_to_zol_ratio, color = group)) +
    geom_hline(yintercept = 0, color = "grey50", linewidth = 0.35) +
    geom_boxplot(outlier.shape = NA, width = 0.5) +
    geom_point(aes(shape = sample_name), position = position_jitter(width = 0.08), size = 2, alpha = 0.9) +
    scale_color_manual(values = color_group, drop = FALSE, guide = "none") +
    theme_test() +
    theme(axis.title.x = element_blank()) +
    labs(y = "log2(PAN effector Vd2 / ZOL effector Vd2)", shape = "Sample", title = "PAN-to-ZOL effector Vd2 abundance ratio")
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
    filter(as.character(cell_type) %in% target_cell_types) %>%
    pivot_longer(cols = all_of(genes), names_to = "gene", values_to = "expression")
  expression_data %>%
    group_by(gene, group, sample_name, cell_type) %>%
    dplyr::summarise(
      mean_expression = mean(expression, na.rm = TRUE),
      median_expression = median(expression, na.rm = TRUE),
      percent_expressing = mean(expression > 0, na.rm = TRUE) * 100,
      n_cells = n(),
      .groups = "drop"
    )
}

plot_marker_dotplot <- function(seurat_obj, genes) {
  plot_genes <- present_genes(seurat_obj, genes)
  if (length(plot_genes) == 0) {
    return(plot_empty("No marker genes available in object"))
  }
  DotPlot2(
    subset_target_cells(seurat_obj),
    features = plot_genes,
    group.by = "cell_type",
    cols = c("lightgrey", "#B2182B"),
    show_grid = FALSE,
    flip = TRUE
  ) +
    theme(axis.title = element_blank(), axis.text.x = element_text(angle = 35, hjust = 1)) +
    labs(title = "ZOL/PAN effector Vd2 marker program")
}

plot_marker_heatmap <- function(marker_summary) {
  if (nrow(marker_summary) == 0) {
    return(plot_empty("No marker expression summary available"))
  }
  plot_data <- marker_summary %>%
    group_by(gene, cell_type) %>%
    dplyr::summarise(mean_expression = mean(mean_expression, na.rm = TRUE), .groups = "drop") %>%
    group_by(gene) %>%
    mutate(
      scaled_expression = as.numeric(scale(mean_expression)),
      scaled_expression = ifelse(is.na(scaled_expression), 0, scaled_expression)
    ) %>%
    ungroup() %>%
    mutate(
      cell_type = factor(as.character(cell_type), levels = target_cell_types),
      gene = factor(gene, levels = rev(unique(gene)))
    )
  ggplot(plot_data, aes(x = cell_type, y = gene, fill = scaled_expression)) +
    geom_tile(color = "white", linewidth = 0.25) +
    scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0) +
    theme_test() +
    theme(axis.title = element_blank(), axis.text.x = element_text(angle = 25, hjust = 1)) +
    labs(fill = "Row-scaled\nmean", title = "Average marker expression")
}

plot_marker_by_sample <- function(marker_summary, genes) {
  plot_data <- marker_summary %>%
    filter(gene %in% genes) %>%
    mutate(cell_type = factor(as.character(cell_type), levels = target_cell_types))
  if (nrow(plot_data) == 0) {
    return(plot_empty("No sample-level marker expression data"))
  }
  ggplot(plot_data, aes(x = cell_type, y = mean_expression, color = cell_type)) +
    geom_boxplot(outlier.shape = NA, width = 0.55) +
    geom_point(aes(shape = sample_name), position = position_jitter(width = 0.08), size = 1.8, alpha = 0.9) +
    facet_wrap(~gene, scales = "free_y", ncol = 4) +
    scale_color_manual(values = color_celltype, drop = FALSE, guide = "none") +
    theme_test() +
    theme(axis.title.x = element_blank(), axis.text.x = element_text(angle = 25, hjust = 1), strip.background = element_blank()) +
    labs(y = "Mean expression by sample", shape = "Sample", title = "Sample-level marker expression")
}

add_module_scores <- function(seurat_obj, gene_sets) {
  for (module_name in names(gene_sets)) {
    module_genes <- present_genes(seurat_obj, gene_sets[[module_name]])
    score_col <- paste0(module_name, "_score")
    if (length(module_genes) == 0) {
      seurat_obj@meta.data[[score_col]] <- NA_real_
      message("[SKIP] No genes found for module score: ", module_name)
      next
    }
    seurat_obj <- AddModuleScore(seurat_obj, features = list(module_genes), name = score_col)
    seurat_obj@meta.data[[score_col]] <- seurat_obj@meta.data[[paste0(score_col, "1")]]
  }
  seurat_obj
}

make_module_score_long <- function(seurat_obj) {
  score_columns <- paste0(names(module_gene_sets), "_score")
  seurat_obj@meta.data %>%
    rownames_to_column("cell_id") %>%
    as_tibble() %>%
    filter(as.character(cell_type) %in% target_cell_types) %>%
    select(cell_id, group, sample_name, cell_type, all_of(score_columns)) %>%
    pivot_longer(cols = all_of(score_columns), names_to = "module", values_to = "score") %>%
    mutate(
      module = str_replace(module, "_score$", ""),
      cell_type = factor(as.character(cell_type), levels = target_cell_types),
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
  comparisons <- list(target_cell_types)
  plot_data <- module_scores %>%
    mutate(
      cell_type = factor(as.character(cell_type), levels = target_cell_types),
      module = factor(as.character(module), levels = names(module_gene_sets))
    )
  ggplot(plot_data, aes(x = cell_type, y = score, fill = cell_type)) +
    geom_violin(scale = "width", trim = TRUE, width = 0.72, alpha = 0.82, linewidth = 0.18) +
    geom_boxplot(width = 0.08, outlier.shape = NA, color = "grey25", alpha = 0.7, linewidth = 0.18) +
    stat_compare_means(
      comparisons = comparisons,
      method = "wilcox.test",
      label = "p.signif",
      hide.ns = FALSE,
      size = 3
    ) +
    facet_wrap(~module, scales = "free_y", ncol = 3, drop = FALSE) +
    scale_fill_manual(values = color_celltype, drop = FALSE, guide = "none") +
    scale_x_discrete(drop = FALSE, expand = expansion(mult = c(0.02, 0.02))) +
    coord_cartesian(clip = "off") +
    theme_test() +
    theme(
      axis.title.x = element_blank(),
      axis.text.x = element_text(angle = 25, hjust = 1, vjust = 1),
      strip.background = element_blank(),
      strip.text = element_text(margin = margin(1.5, 0, 1.5, 0)),
      panel.spacing.x = grid::unit(0.25, "lines"),
      panel.spacing.y = grid::unit(0.25, "lines"),
      plot.margin = margin(4, 6, 4, 4)
    ) +
    labs(y = "Module score", title = "Functional module scores with Wilcoxon significance")
}

plot_module_score_by_sample <- function(module_summary) {
  if (nrow(module_summary) == 0) {
    return(plot_empty("No module score summary available"))
  }
  ggplot(module_summary, aes(x = cell_type, y = mean_score, color = cell_type)) +
    geom_boxplot(outlier.shape = NA, width = 0.55) +
    geom_point(aes(shape = sample_name), position = position_jitter(width = 0.08), size = 1.9, alpha = 0.9) +
    facet_wrap(~module, scales = "free_y", ncol = 3) +
    scale_color_manual(values = color_celltype, drop = FALSE, guide = "none") +
    theme_test() +
    theme(axis.title.x = element_blank(), axis.text.x = element_text(angle = 25, hjust = 1), strip.background = element_blank()) +
    labs(y = "Mean module score by sample", shape = "Sample", title = "Sample-level functional module scores")
}

run_de <- function(seurat_obj) {
  comparison_obj <- subset_target_cells(seurat_obj)
  Idents(comparison_obj) <- factor(as.character(comparison_obj$cell_type), levels = comparison_pair)
  if (!all(comparison_pair %in% levels(Idents(comparison_obj))) || min(table(Idents(comparison_obj))) < 3) {
    message("[SKIP] Too few cells for DE: ", comparison_name)
    return(tibble())
  }
  FindMarkers(
    comparison_obj,
    ident.1 = comparison_pair[1],
    ident.2 = comparison_pair[2],
    only.pos = FALSE,
    logfc.threshold = 0,
    min.pct = 0.05
  ) %>%
    rownames_to_column("gene") %>%
    as_tibble() %>%
    mutate(
      comparison = comparison_name,
      ident_1 = comparison_pair[1],
      ident_2 = comparison_pair[2],
      higher_in = ifelse(avg_log2FC >= 0, ident_1, ident_2),
      abs_log2FC = abs(avg_log2FC),
      .before = 1
    )
}

plot_de_volcano <- function(de_markers, top_n = 12) {
  if (nrow(de_markers) == 0) {
    return(plot_empty("No DE marker data available"))
  }
  plot_data <- de_markers %>%
    mutate(
      neg_log10_padj = -log10(pmax(p_val_adj, .Machine$double.xmin)),
      significant = p_val_adj < 0.05 & abs_log2FC >= 0.25
    )
  label_data <- plot_data %>%
    filter(significant) %>%
    slice_max(order_by = abs_log2FC, n = top_n, with_ties = FALSE)
  volcano_plot <- ggplot(plot_data, aes(x = avg_log2FC, y = neg_log10_padj, color = higher_in)) +
    geom_point(aes(alpha = significant), size = 0.9) +
    geom_vline(xintercept = 0, color = "grey45", linewidth = 0.35) +
    scale_alpha_manual(values = c(`TRUE` = 0.85, `FALSE` = 0.25), guide = "none") +
    scale_color_manual(values = color_celltype, drop = FALSE) +
    theme_test() +
    labs(x = "avg_log2FC: PAN Effector Vd2 - ZOL Effector Vd2", y = "-log10 adjusted P", color = "Higher in", title = "PAN vs ZOL effector Vd2 DE")
  if (nrow(label_data) > 0 && requireNamespace("ggrepel", quietly = TRUE)) {
    volcano_plot <- volcano_plot +
      ggrepel::geom_text_repel(data = label_data, aes(label = gene), size = 3, max.overlaps = Inf, show.legend = FALSE)
  } else if (nrow(label_data) > 0) {
    volcano_plot <- volcano_plot +
      geom_text(data = label_data, aes(label = gene), size = 3, vjust = -0.5, show.legend = FALSE)
  }
  volcano_plot
}

plot_de_lollipop <- function(de_markers, top_n = 15) {
  if (nrow(de_markers) == 0) {
    return(plot_empty("No DE marker data available"))
  }
  plot_data <- de_markers %>%
    filter(!is.na(p_val_adj), p_val_adj < 0.05) %>%
    group_by(higher_in) %>%
    slice_max(order_by = abs_log2FC, n = top_n, with_ties = FALSE) %>%
    ungroup() %>%
    mutate(gene = factor(gene, levels = gene[order(avg_log2FC)]))
  if (nrow(plot_data) == 0) {
    return(plot_empty("No significant DE genes for lollipop plot"))
  }
  ggplot(plot_data, aes(x = avg_log2FC, y = gene, color = higher_in)) +
    geom_vline(xintercept = 0, color = "grey55", linewidth = 0.35) +
    geom_segment(aes(x = 0, xend = avg_log2FC, yend = gene), linewidth = 0.8) +
    geom_point(size = 2.4) +
    scale_color_manual(values = color_celltype, drop = FALSE) +
    theme_test() +
    theme(axis.title.y = element_blank()) +
    labs(x = "avg_log2FC: PAN Effector Vd2 - ZOL Effector Vd2", color = "Higher in", title = "Top DE genes")
}

plot_top_de_heatmap <- function(de_markers, seurat_obj, top_n = 30) {
  if (nrow(de_markers) == 0) {
    return(plot_empty("No DE marker data available"))
  }
  top_genes <- de_markers %>%
    filter(!is.na(p_val_adj), p_val_adj < 0.05) %>%
    slice_max(order_by = abs_log2FC, n = top_n, with_ties = FALSE) %>%
    pull(gene) %>%
    unique()
  top_genes <- present_genes(seurat_obj, top_genes)
  if (length(top_genes) == 0) {
    return(plot_empty("No significant DE genes found in object"))
  }
  marker_summary <- make_marker_summary(seurat_obj, top_genes)
  plot_marker_heatmap(marker_summary) +
    labs(title = "Top DE gene expression")
}

calculate_hallmark_delta <- function(seurat_obj) {
  analysis_obj <- subset_target_cells(seurat_obj)
  analysis_obj <- GeneSetAnalysis(analysis_obj, genesets = hall50$human)
  scores <- analysis_obj@misc$AUCell$genesets
  if (is.null(scores)) {
    stop("Missing AUCell Hallmark scores after GeneSetAnalysis().", call. = FALSE)
  }
  metadata <- analysis_obj@meta.data %>%
    rownames_to_column("cell_id") %>%
    as_tibble() %>%
    select(cell_id, cell_type)
  common_cells <- intersect(colnames(scores), metadata$cell_id)
  scores <- scores[, common_cells, drop = FALSE]
  metadata <- metadata %>% filter(cell_id %in% common_cells)
  celltype_stats <- CalcStats(scores, f = metadata$cell_type)
  delta_table <- tibble(
    comparison = comparison_name,
    ident_1 = comparison_pair[1],
    ident_2 = comparison_pair[2],
    pathway = rownames(scores),
    mean_ident_1 = rowMeans(scores[, metadata$cell_id[as.character(metadata$cell_type) == comparison_pair[1]], drop = FALSE], na.rm = TRUE),
    mean_ident_2 = rowMeans(scores[, metadata$cell_id[as.character(metadata$cell_type) == comparison_pair[2]], drop = FALSE], na.rm = TRUE)
  ) %>%
    mutate(
      delta = mean_ident_1 - mean_ident_2,
      abs_delta = abs(delta),
      higher_in = ifelse(delta >= 0, ident_1, ident_2)
    ) %>%
    arrange(desc(abs_delta), pathway)
  list(scores = scores, stats = celltype_stats, delta = delta_table)
}

plot_hallmark_delta <- function(delta_table, top_n = 12) {
  if (nrow(delta_table) == 0) {
    return(plot_empty("No Hallmark delta data available"))
  }
  plot_data <- delta_table %>%
    slice_max(order_by = abs_delta, n = top_n, with_ties = FALSE) %>%
    mutate(
      pathway_label = str_replace(pathway, "^HALLMARK_", ""),
      pathway_label = str_replace_all(pathway_label, "_", " "),
      pathway_label = factor(pathway_label, levels = pathway_label[order(delta)])
    )
  ggplot(plot_data, aes(x = delta, y = pathway_label, color = higher_in)) +
    geom_segment(aes(x = 0, xend = delta, yend = pathway_label), linewidth = 0.8) +
    geom_point(size = 2.5) +
    geom_vline(xintercept = 0, color = "grey45", linewidth = 0.35) +
    scale_color_manual(values = color_celltype, drop = FALSE) +
    theme_test() +
    theme(axis.title.y = element_blank()) +
    labs(x = "Mean AUCell delta: PAN - ZOL", color = "Higher in", title = "Hallmark pathway delta")
}

plot_selected_hallmark_heatmap <- function(stats_matrix) {
  if (is.null(stats_matrix) || nrow(stats_matrix) == 0) {
    return(plot_empty("No Hallmark stats available"))
  }
  pathways <- intersect(selected_hallmark_pathways, rownames(stats_matrix))
  if (length(pathways) == 0) {
    return(plot_empty("Selected Hallmark pathways not available"))
  }
  plot_data <- as.data.frame(stats_matrix[pathways, , drop = FALSE]) %>%
    rownames_to_column("pathway") %>%
    as_tibble() %>%
    pivot_longer(cols = -pathway, names_to = "cell_type", values_to = "score") %>%
    mutate(
      pathway = str_replace(pathway, "^HALLMARK_", ""),
      pathway = str_replace_all(pathway, "_", " "),
      pathway = factor(pathway, levels = rev(str_replace_all(str_replace(selected_hallmark_pathways, "^HALLMARK_", ""), "_", " "))),
      cell_type = factor(cell_type, levels = target_cell_types)
    )
  ggplot(plot_data, aes(x = cell_type, y = pathway, fill = score)) +
    geom_tile(color = "white", linewidth = 0.25) +
    scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0) +
    theme_test() +
    theme(axis.title = element_blank(), axis.text.x = element_text(angle = 25, hjust = 1)) +
    labs(fill = "AUCell z-score", title = "Selected Hallmark programs")
}

calculate_clone_size <- function(annotation) {
  if (is.null(annotation) || nrow(annotation) == 0) {
    return(tibble())
  }
  annotation %>%
    filter(as.character(cell_type) %in% target_cell_types, !is.na(cdr3), cdr3 != "") %>%
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
      clone_size_class = factor(clone_size_class, levels = c("Rare (<0.1%)", "Large (0.1-1%)", "Expanded (1-10%)", "Hyperexpanded (>10%)"))
    ) %>%
    ungroup()
}

plot_clone_size <- function(clone_size) {
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
    mutate(cell_type = factor(as.character(cell_type), levels = target_cell_types))
  ggplot(plot_data, aes(x = cell_type, y = percent, fill = clone_size_class)) +
    geom_col(width = 0.72, color = "white", linewidth = 0.2) +
    facet_wrap(~chain, ncol = 1) +
    scale_fill_manual(values = clone_size_colors, drop = FALSE) +
    theme_test() +
    theme(axis.title.x = element_blank(), axis.text.x = element_text(angle = 25, hjust = 1), strip.background = element_blank()) +
    labs(y = "Cells in clone-size class (%)", fill = "Clone-size class", title = "Clone-size classes in ZOL/PAN effector Vd2")
}

calculate_top_cdr3 <- function(annotation, top_n = 10) {
  if (is.null(annotation) || nrow(annotation) == 0) {
    return(tibble())
  }
  annotation %>%
    filter(as.character(cell_type) %in% target_cell_types, !is.na(cdr3), cdr3 != "") %>%
    distinct(barcode, cell_type, chain, cdr3) %>%
    dplyr::count(cell_type, chain, cdr3, name = "n_cells") %>%
    group_by(cell_type, chain) %>%
    mutate(percent = n_cells / sum(n_cells) * 100) %>%
    slice_max(order_by = n_cells, n = top_n, with_ties = FALSE) %>%
    arrange(cell_type, chain, desc(n_cells)) %>%
    ungroup()
}

plot_top_cdr3 <- function(top_cdr3) {
  if (nrow(top_cdr3) == 0) {
    return(plot_empty("No top CDR3 data available"))
  }
  plot_data <- top_cdr3 %>%
    mutate(
      cell_type = factor(as.character(cell_type), levels = target_cell_types),
      cdr3_label = str_trunc(cdr3, width = 24)
    )
  ggplot(plot_data, aes(x = reorder(cdr3_label, percent), y = percent, fill = cell_type)) +
    geom_col(width = 0.75) +
    coord_flip() +
    facet_grid(chain ~ cell_type, scales = "free_y", space = "free_y") +
    scale_fill_manual(values = color_celltype, drop = FALSE, guide = "none") +
    theme_test() +
    theme(axis.title.y = element_blank(), strip.background = element_blank()) +
    labs(y = "Cells with CDR3 (%)", title = "Top CDR3 clones in ZOL/PAN effector Vd2")
}

prepare_paired_clone_metadata <- function(paired_cdr3, seurat_obj) {
  if (is.null(paired_cdr3) || nrow(paired_cdr3) == 0) {
    return(tibble())
  }
  cell_metadata <- seurat_obj@meta.data %>%
    rownames_to_column("cell_id") %>%
    as_tibble() %>%
    select(cell_id, barcode, group, sample_name, cell_type) %>%
    filter(as.character(cell_type) %in% target_cell_types)
  paired_cdr3 %>%
    select(barcode, TRD, TRG) %>%
    distinct() %>%
    inner_join(cell_metadata, by = "barcode") %>%
    filter(!is.na(TRD), !is.na(TRG), TRD != "", TRG != "") %>%
    mutate(
      TRDG_pair = paste(TRD, TRG, sep = "||"),
      cell_type = factor(as.character(cell_type), levels = target_cell_types),
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

plot_paired_clone_composition <- function(paired_clone_summary) {
  if (nrow(paired_clone_summary) == 0) {
    return(plot_empty("No paired TRD+TRG clone data available"))
  }
  plot_data <- paired_clone_summary %>%
    mutate(
      pair_rank = factor(pair_rank, levels = paste0("pair_", seq_len(length(unique(pair_rank))))),
      cell_type = factor(as.character(cell_type), levels = target_cell_types)
    ) %>%
    group_by(pair_rank, cell_type) %>%
    dplyr::summarise(n_cells = sum(n_cells), .groups = "drop_last") %>%
    mutate(percent = n_cells / sum(n_cells) * 100) %>%
    ungroup()
  ggplot(plot_data, aes(x = pair_rank, y = percent, fill = cell_type)) +
    geom_col(width = 0.72, color = "white", linewidth = 0.2) +
    scale_fill_manual(values = color_celltype, drop = FALSE) +
    theme_test() +
    theme(axis.title.x = element_blank(), axis.text.x = element_text(angle = 35, hjust = 1)) +
    labs(y = "Cells in paired clone (%)", fill = "Cell type", title = "Top paired TRD+TRG clone composition")
}

calculate_paired_clone_overlap <- function(paired_metadata) {
  if (nrow(paired_metadata) == 0) {
    return(tibble())
  }
  paired_metadata %>%
    distinct(TRDG_pair, cell_type) %>%
    dplyr::count(TRDG_pair, name = "celltype_n") %>%
    right_join(
      paired_metadata %>%
        dplyr::count(TRDG_pair, cell_type, name = "n_cells") %>%
        group_by(TRDG_pair) %>%
        mutate(total_pair_cells = sum(n_cells)) %>%
        ungroup(),
      by = "TRDG_pair"
    ) %>%
    arrange(desc(total_pair_cells), TRDG_pair, cell_type)
}

plot_paired_clone_overlap_heatmap <- function(overlap_table, top_n = 40) {
  if (nrow(overlap_table) == 0) {
    return(plot_empty("No paired TRD+TRG clone overlap data available"))
  }
  top_pairs <- overlap_table %>%
    distinct(TRDG_pair, total_pair_cells) %>%
    arrange(desc(total_pair_cells), TRDG_pair) %>%
    slice_head(n = top_n) %>%
    pull(TRDG_pair)
  plot_data <- overlap_table %>%
    filter(TRDG_pair %in% top_pairs) %>%
    mutate(
      pair_label = factor(str_trunc(TRDG_pair, width = 36), levels = rev(str_trunc(top_pairs, width = 36))),
      cell_type = factor(as.character(cell_type), levels = target_cell_types)
    )
  ggplot(plot_data, aes(x = cell_type, y = pair_label, fill = n_cells)) +
    geom_tile(color = "white", linewidth = 0.2) +
    scale_fill_gradient(low = "white", high = "#B2182B") +
    theme_test() +
    theme(axis.title = element_blank(), axis.text.x = element_text(angle = 25, hjust = 1)) +
    labs(fill = "Cells", title = "Paired TRD+TRG clone overlap between ZOL/PAN effector Vd2")
}

plot_clone_alluvial <- function(paired_clone_summary) {
  if (nrow(paired_clone_summary) == 0) {
    return(plot_empty("No paired TRD+TRG clone data available"))
  }
  if (!requireNamespace("ggalluvial", quietly = TRUE)) {
    return(plot_empty("ggalluvial is not installed", "Install ggalluvial to draw sample -> cell type -> paired clone alluvial plots."))
  }
  plot_data <- paired_clone_summary %>%
    mutate(
      pair_rank = factor(pair_rank, levels = paste0("pair_", seq_len(length(unique(pair_rank))))),
      cell_type = factor(as.character(cell_type), levels = target_cell_types),
      sample_name = as.character(sample_name)
    )
  ggplot(plot_data, aes(y = n_cells, axis1 = sample_name, axis2 = cell_type, axis3 = pair_rank)) +
    ggalluvial::geom_alluvium(aes(fill = cell_type), alpha = 0.65, width = 1 / 12) +
    ggalluvial::geom_stratum(width = 1 / 6, fill = "grey95", color = "grey55") +
    ggalluvial::stat_stratum(geom = "text", aes(label = after_stat(stratum)), size = 3) +
    scale_x_discrete(limits = c("Sample", "Cell type", "Top TRD+TRG pair"), expand = c(0.08, 0.08)) +
    scale_fill_manual(values = color_celltype, drop = FALSE) +
    theme_test() +
    theme(axis.title = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
    labs(fill = "Cell type", title = "Sample to ZOL/PAN effector Vd2 state to top paired clone")
}

# 1. Load final annotated Seurat object.
all_seurat_celltype <- read_rds_checked(seurat_celltype_rds, "annotated Seurat RDS")
all_seurat_celltype <- normalise_metadata_levels(all_seurat_celltype)
all_seurat_celltype <- prepare_expression_object(all_seurat_celltype)
target_obj <- subset_target_cells(all_seurat_celltype)
marker_genes <- present_genes(target_obj, marker_gene_panel)

# 2. UMAP context and abundance.
save_plot(make_umap_highlight(all_seurat_celltype), "zol_pan_effector_vd2_umap_highlight", 8, 6, overwrite = force_zol_pan_effector_vd2_plot)
save_plot(make_sample_umap(all_seurat_celltype), "zol_pan_effector_vd2_umap_by_sample", 8, 6, overwrite = force_zol_pan_effector_vd2_plot)
save_plot(make_umap_density(all_seurat_celltype), "zol_pan_effector_vd2_umap_density", 9, 5, overwrite = force_zol_pan_effector_vd2_plot)

cell_fraction <- calculate_cell_fraction(all_seurat_celltype)
write_csv_if_missing(cell_fraction, cell_fraction_csv, "ZOL/PAN effector Vd2 cell fraction CSV", overwrite = force_zol_pan_effector_vd2)
save_plot(plot_cell_fraction_stacked(cell_fraction), "zol_pan_effector_vd2_cell_fraction_by_sample", 8, 5, overwrite = force_zol_pan_effector_vd2_plot)
save_plot(plot_cell_fraction_points(cell_fraction), "zol_pan_effector_vd2_cell_fraction_points", 6, 5, overwrite = force_zol_pan_effector_vd2_plot)
save_plot(plot_cell_fraction_ratio(cell_fraction), "zol_pan_effector_vd2_pan_to_zol_ratio", 7, 5, overwrite = force_zol_pan_effector_vd2_plot)

# 3. Marker program and module-score comparisons.
marker_summary <- make_marker_summary(target_obj, marker_genes)
write_csv_if_missing(marker_summary, marker_summary_csv, "ZOL/PAN effector Vd2 marker summary CSV", overwrite = force_zol_pan_effector_vd2)
save_plot(plot_marker_dotplot(target_obj, marker_genes), "zol_pan_effector_vd2_marker_dotplot", 11, 8, overwrite = force_zol_pan_effector_vd2_plot)
save_plot(plot_marker_heatmap(marker_summary), "zol_pan_effector_vd2_marker_heatmap", 7, max(6, length(marker_genes) * 0.25), overwrite = force_zol_pan_effector_vd2_plot)
save_plot(plot_marker_by_sample(marker_summary, head(marker_genes, 24)), "zol_pan_effector_vd2_marker_by_sample", 12, 7, overwrite = force_zol_pan_effector_vd2_plot)

target_obj <- add_module_scores(target_obj, module_gene_sets)
module_scores <- make_module_score_long(target_obj)
module_score_summary <- summarise_module_scores(module_scores)
write_csv_if_missing(module_score_summary, module_score_summary_csv, "ZOL/PAN effector Vd2 module score summary CSV", overwrite = force_zol_pan_effector_vd2)
save_plot(plot_module_score_violin(module_scores), "zol_pan_effector_vd2_module_score_violin", 12, 7, overwrite = force_zol_pan_effector_vd2_plot)
save_plot(plot_module_score_by_sample(module_score_summary), "zol_pan_effector_vd2_module_score_by_sample", 12, 7, overwrite = force_zol_pan_effector_vd2_plot)

# 4. Differential expression and Hallmark pathway comparison.
de_markers <- run_de(target_obj)
write_csv_if_missing(de_markers, de_markers_csv, "ZOL/PAN effector Vd2 DE marker CSV", overwrite = force_zol_pan_effector_vd2)
save_plot(plot_de_volcano(de_markers), "zol_pan_effector_vd2_de_volcano", 7, 6, overwrite = force_zol_pan_effector_vd2_plot)
save_plot(plot_de_lollipop(de_markers), "zol_pan_effector_vd2_top_de_lollipop", 8, 7, overwrite = force_zol_pan_effector_vd2_plot)
save_plot(plot_top_de_heatmap(de_markers, target_obj), "zol_pan_effector_vd2_top_de_heatmap", 7, 8, overwrite = force_zol_pan_effector_vd2_plot)

hallmark_result <- calculate_hallmark_delta(target_obj)
write_csv_if_missing(hallmark_result$delta, hallmark_delta_csv, "ZOL/PAN effector Vd2 Hallmark delta CSV", overwrite = force_zol_pan_effector_vd2)
write_csv_if_missing(as.data.frame(hallmark_result$stats) %>% rownames_to_column("pathway"), hallmark_stats_csv, "ZOL/PAN effector Vd2 Hallmark stats CSV", overwrite = force_zol_pan_effector_vd2)
save_plot(plot_hallmark_delta(hallmark_result$delta), "zol_pan_effector_vd2_hallmark_delta", 8, 6, overwrite = force_zol_pan_effector_vd2_plot)
save_plot(plot_selected_hallmark_heatmap(hallmark_result$stats), "zol_pan_effector_vd2_selected_hallmark_heatmap", 7, 5, overwrite = force_zol_pan_effector_vd2_plot)

# 5. Repertoire-aware comparison.
all_annotation_included <- read_optional_rds(included_annotation_rds, "included VDJ annotation RDS")
if (!is.null(all_annotation_included)) {
  all_annotation_included <- all_annotation_included %>%
    mutate(
      group = factor(as.character(group), levels = group_levels),
      sample_name = factor(as.character(sample_name), levels = sample_name_levels),
      cell_type = factor(as.character(cell_type), levels = cell_type_levels)
    )
  clone_size <- calculate_clone_size(all_annotation_included)
  top_cdr3 <- calculate_top_cdr3(all_annotation_included)
  write_csv_if_missing(clone_size, clone_size_csv, "ZOL/PAN effector Vd2 clone-size CSV", overwrite = force_zol_pan_effector_vd2)
  write_csv_if_missing(top_cdr3, top_cdr3_csv, "ZOL/PAN effector Vd2 top CDR3 CSV", overwrite = force_zol_pan_effector_vd2)
  save_plot(plot_clone_size(clone_size), "zol_pan_effector_vd2_clone_size", 8, 6, overwrite = force_zol_pan_effector_vd2_plot)
  save_plot(plot_top_cdr3(top_cdr3), "zol_pan_effector_vd2_top_cdr3", 10, 8, overwrite = force_zol_pan_effector_vd2_plot)
}

barcode_trgd_paired <- read_optional_rds(paired_cdr3_rds, "paired TRD/TRG CDR3 RDS")
if (!is.null(barcode_trgd_paired)) {
  paired_metadata <- prepare_paired_clone_metadata(barcode_trgd_paired, all_seurat_celltype)
  paired_clone_summary <- summarise_paired_clones(paired_metadata)
  paired_clone_overlap <- calculate_paired_clone_overlap(paired_metadata)
  write_csv_if_missing(paired_clone_summary, paired_clone_summary_csv, "ZOL/PAN effector Vd2 paired clone summary CSV", overwrite = force_zol_pan_effector_vd2)
  write_csv_if_missing(paired_clone_overlap, paired_clone_overlap_csv, "ZOL/PAN effector Vd2 paired clone overlap CSV", overwrite = force_zol_pan_effector_vd2)
  save_plot(plot_paired_clone_composition(paired_clone_summary), "zol_pan_effector_vd2_paired_clone_composition", 9, 5, overwrite = force_zol_pan_effector_vd2_plot)
  save_plot(plot_paired_clone_overlap_heatmap(paired_clone_overlap), "zol_pan_effector_vd2_paired_clone_overlap_heatmap", 8, 10, overwrite = force_zol_pan_effector_vd2_plot)
  save_plot(plot_clone_alluvial(paired_clone_summary), "zol_pan_effector_vd2_sample_celltype_clone_alluvial", 12, 6, overwrite = force_zol_pan_effector_vd2_plot)
}
