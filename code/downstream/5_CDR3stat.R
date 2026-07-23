# R version 4.5.2 (2025-10-31)
rm(list = ls())
setwd("/data/huotong/project_tcr/2026May")
library(Seurat)
library(SeuratExtend)
library(tidyverse)
library(ggseqlogo)
library(gplots)
library(ggpubr)
library(patchwork)
library(ggalluvial)
options(
  tibble.width = Inf,
  print.width = Inf,
  max.print = 200
)

rds_dir <- "rds"
figure_dir <- file.path("figures", "5_CDR3stat")
table_dir <- "table"

seurat_celltype_rds <- file.path(rds_dir, "all_seurat_celltype.rds")
included_annotation_rds <- file.path(rds_dir, "all_annotation_included.rds")
tcr_diversity_metrics_csv <- file.path(table_dir, "tcr_diversity_clonality_metrics.csv")
tcr_clone_size_classes_csv <- file.path(table_dir, "tcr_clone_size_classes.csv")

top_cdr3_n <- 100
top_cdr3_plot_limit <- 10
logo_grid_width <- 8
logo_height_per_group <- 2.2
logo_min_height <- 4
cdr3_length_a4_width <- 11.69
cdr3_length_a4_height <- 8.27
expanded_cdr3_to_exclude <- c("CALDTTFPIGDRGYTDKLIF")
expanded_cdr3_label <- "CALDTTFPIGDRGYTDKLIF"
force_cdr3_stat <- FALSE
force_cdr3_stat_plot <- TRUE

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
vd2_cell_types <- c(
  "Effector Memory Vd2",
  "Pre-activated Vd2",
  "ZOL Effector Vd2",
  "ZOL FOXP3+ Vd2",
  "PAN Effector Vd2"
)
vd1_cell_types <- c(
  "Naive Vd1",
  "Pre-activated Vd1",
  "Effector Vd1"
)
cdr3_logo_sets <- list(
  list(
    chain = "TRD",
    length = 19,
    prefix = "VDJ_CDR3d_length_19_TRDV1_sequence_logo",
    title = "TRDV1 CDR3d Length 19 AA",
    groups = list(
      "Naive Vd1" = "Naive Vd1",
      "Pre-activated Vd1" = "Pre-activated Vd1",
      "Effector Vd1" = "Effector Vd1"
    )
  ),
  list(
    chain = "TRD",
    length = 16,
    prefix = "VDJ_CDR3d_length_16_TRDV2_sequence_logo",
    title = "TRDV2 CDR3d Length 16 AA",
    groups = list(
      "Effector Memory Vd2" = "Effector Memory Vd2",
      "Pre-activated Vd2" = "Pre-activated Vd2",
      "ZOL Effector Vd2" = "ZOL Effector Vd2",
      "ZOL FOXP3+ Vd2" = "ZOL FOXP3+ Vd2",
      "PAN Effector Vd2" = "PAN Effector Vd2"
    )
  ),
  list(
    chain = "TRD",
    length = 18,
    prefix = "VDJ_CDR3d_length_18_TRDV2_sequence_logo",
    title = "TRDV2 CDR3d Length 18 AA",
    groups = list(
      "Effector Memory Vd2" = "Effector Memory Vd2",
      "Pre-activated Vd2" = "Pre-activated Vd2",
      "ZOL Effector Vd2" = "ZOL Effector Vd2",
      "ZOL FOXP3+ Vd2" = "ZOL FOXP3+ Vd2",
      "PAN Effector Vd2" = "PAN Effector Vd2"
    )
  ),
  list(
    chain = "TRG",
    length = 13,
    prefix = "VDJ_CDR3g_length_13_TRDV1_sequence_logo",
    title = "TRDV1 CDR3g Length 13 AA",
    groups = list(
      "Naive Vd1" = "Naive Vd1",
      "Pre-activated Vd1" = "Pre-activated Vd1",
      "Effector Vd1" = "Effector Vd1"
    )
  ),
  list(
    chain = "TRG",
    length = 16,
    prefix = "VDJ_CDR3g_length_16_TRDV2_sequence_logo",
    title = "TRDV2 CDR3g Length 16 AA",
    groups = list(
      "Effector Memory Vd2" = "Effector Memory Vd2",
      "Pre-activated Vd2" = "Pre-activated Vd2",
      "ZOL Effector Vd2" = "ZOL Effector Vd2",
      "ZOL FOXP3+ Vd2" = "ZOL FOXP3+ Vd2",
      "PAN Effector Vd2" = "PAN Effector Vd2"
    )
  )
)

dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, showWarnings = FALSE)

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

summarise_cdr3_frequency <- function(annotation, group_columns, top_n = top_cdr3_n) {
  annotation %>%
    group_by(across(all_of(c(group_columns, "chain", "cdr3")))) %>%
    summarise(sum_umis = sum(umis), n = n(), .groups = "drop_last") %>%
    group_by(across(all_of(c(group_columns, "chain")))) %>%
    mutate(percent = n / sum(n) * 100) %>%
    arrange(desc(percent), .by_group = TRUE) %>%
    slice_head(n = top_n) %>%
    mutate(rank = row_number()) %>%
    ungroup()
}

plot_top_cdr3_by_sample <- function(annotation, title = NULL) {
  summarise_cdr3_frequency(annotation, c("sample_name", "group")) %>%
    ggplot(aes(x = rank, y = percent, color = group)) +
    geom_line() +
    geom_point() +
    facet_grid(sample_name ~ chain) +
    theme_test() +
    theme(
      legend.position = "bottom",
      strip.background = element_rect(fill = "white")
    ) +
    scale_color_manual(values = color_group) +
    scale_x_continuous(
      breaks = seq(0, top_cdr3_plot_limit, by = 1),
      limits = c(1, top_cdr3_plot_limit)
    ) +
    labs(
      title = title,
      x = "Rank of CDR3 within each sample and chain",
      y = "Percentage of cells with CDR3"
    )
}

plot_top_cdr3_by_cell_type <- function(annotation) {
  summarise_cdr3_frequency(annotation, c("cell_type")) %>%
    ggplot(aes(x = rank, y = percent, color = cell_type)) +
    geom_line() +
    geom_point() +
    facet_grid(cell_type ~ chain) +
    theme_test() +
    theme(
      legend.position = "bottom",
      strip.background = element_blank()
    ) +
    scale_color_manual(values = color_celltype) +
    scale_x_continuous(
      breaks = seq(0, top_cdr3_plot_limit, by = 1),
      limits = c(1, top_cdr3_plot_limit)
    ) +
    labs(
      x = "Rank of CDR3 within each cell type and chain",
      y = "Percentage of cells with CDR3"
    )
}

plot_marked_cdr3_umap <- function(seurat_obj, marked_barcodes) {
  seurat_obj$mark <- ifelse(
    seurat_obj$barcode %in% marked_barcodes,
    expanded_cdr3_label,
    "other"
  )
  umap_plot <- DimPlot2(seurat_obj,
    features = "mark",
    reduction = "umap.unintegrated",
    theme = NoAxes(),
    cols = setNames(c("red", "grey"), c(expanded_cdr3_label, "other"))
  )
  add_fixed_umap_coordinates(umap_plot)
}

plot_marked_cdr3_distribution <- function(seurat_obj, marked_barcodes) {
  mark <- ifelse(
    seurat_obj$barcode %in% marked_barcodes,
    expanded_cdr3_label,
    "other"
  )
  ClusterDistrBar(
    origin = mark,
    cluster = seurat_obj$cell_type,
    cols = color_celltype,
    width = 0.3
  )
}

add_trd_type <- function(annotation) {
  annotation %>%
    mutate(
      cdr3_length = nchar(cdr3),
      trd_type = case_when(
        cell_type %in% vd2_cell_types ~ "Vd2",
        cell_type %in% vd1_cell_types ~ "Vd1",
        TRUE ~ "Other"
      )
    )
}

plot_cdr3_length_distribution <- function(annotation, chain_name) {
  chain_label <- ifelse(chain_name == "TRD", "CDR3d", "CDR3g")
  zoom_xlim <- if (chain_name == "TRD") c(14, 20) else c(12, 18)
  dashed_lines <- if (chain_name == "TRD") c(16, 18) else 16
  solid_lines <- if (chain_name == "TRD") 19 else 13

  plot_data <- add_trd_type(annotation) %>%
    group_by(trd_type, cell_type, chain, cdr3_length) %>%
    tally() %>%
    group_by(trd_type, cell_type, chain) %>%
    mutate(prop = n / sum(n) * 100) %>%
    ungroup() %>%
    filter(chain == chain_name)

  if (nrow(plot_data) == 0) {
    return(ggplot() +
      theme_void() +
      labs(title = paste(chain_label, "length distribution: no data")))
  }

  peak_data <- plot_data %>%
    group_by(trd_type, cell_type) %>%
    slice_max(order_by = prop, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    mutate(
      peak_label = paste0(cdr3_length, " aa\n", sprintf("%.1f%%", prop))
    )

  base_plot <- ggplot(plot_data, aes(cdr3_length, prop, color = cell_type)) +
    geom_line(linewidth = 0.5, aes(linetype = trd_type)) +
    geom_point(
      data = peak_data,
      aes(cdr3_length, prop, fill = cell_type),
      inherit.aes = FALSE,
      shape = 21,
      color = "black",
      size = 2.4,
      stroke = 0.25
    ) +
    theme_bw() +
    theme(
      legend.position = "bottom",
      strip.background = element_rect(fill = "white")
    ) +
    scale_color_manual(values = color_celltype, drop = FALSE) +
    scale_fill_manual(values = color_celltype, drop = FALSE, guide = "none") +
    geom_vline(xintercept = dashed_lines, linetype = "dashed") +
    geom_vline(xintercept = solid_lines, linetype = "solid") +
    labs(
      title = paste(chain_label, "length distribution across cell types"),
      x = "CDR3 length (aa)",
      y = "Cell proportion (%)"
    )

  full_plot <- base_plot +
    facet_wrap(~trd_type, ncol = 1, scales = "free_y") +
    labs(subtitle = "Full length range")

  zoom_data <- filter(plot_data, cdr3_length >= zoom_xlim[1], cdr3_length <= zoom_xlim[2])
  zoom_peak_data <- peak_data %>%
    filter(cdr3_length >= zoom_xlim[1], cdr3_length <= zoom_xlim[2])
  zoom_plot <- ggplot(zoom_data, aes(cdr3_length, prop, color = cell_type)) +
    geom_line(linewidth = 0.7, aes(linetype = trd_type)) +
    geom_point(
      data = zoom_peak_data,
      aes(cdr3_length, prop, fill = cell_type),
      inherit.aes = FALSE,
      shape = 21,
      color = "black",
      size = 3,
      stroke = 0.3
    ) +
    geom_text(
      data = zoom_peak_data,
      aes(cdr3_length, prop, label = peak_label),
      inherit.aes = FALSE,
      size = 2.8,
      vjust = -0.45,
      check_overlap = TRUE
    ) +
    geom_vline(xintercept = dashed_lines, linetype = "dashed") +
    geom_vline(xintercept = solid_lines, linetype = "solid") +
    scale_color_manual(values = color_celltype, drop = FALSE) +
    scale_fill_manual(values = color_celltype, drop = FALSE, guide = "none") +
    facet_wrap(~trd_type, ncol = 1, scales = "free_y") +
    coord_cartesian(xlim = zoom_xlim) +
    theme_bw() +
    theme(
      legend.position = "bottom",
      strip.background = element_rect(fill = "white")
    ) +
    labs(
      subtitle = paste0("Peak-focused window: ", zoom_xlim[1], "-", zoom_xlim[2], " aa"),
      x = "CDR3 length (aa)",
      y = "Cell proportion (%)"
    )

  full_plot / zoom_plot +
    plot_layout(heights = c(2, 1), guides = "collect") &
    theme(legend.position = "bottom")
}

calculate_tcr_diversity_metrics <- function(annotation) {
  clone_counts <- annotation %>%
    filter(!is.na(cdr3), cdr3 != "") %>%
    distinct(barcode, group, sample_name, chain, cdr3) %>%
    count(group, sample_name, chain, cdr3, name = "clone_cells")

  if (nrow(clone_counts) == 0) {
    return(tibble())
  }

  clone_counts %>%
    group_by(group, sample_name, chain) %>%
    mutate(
      total_cells = sum(clone_cells),
      clone_fraction = clone_cells / total_cells
    ) %>%
    summarise(
      total_cells = dplyr::first(total_cells),
      richness = n(),
      shannon = -sum(clone_fraction * log(clone_fraction)),
      simpson = 1 - sum(clone_fraction^2),
      inverse_simpson = 1 / sum(clone_fraction^2),
      gini = calculate_gini(clone_cells),
      max_clone_fraction = max(clone_fraction),
      .groups = "drop"
    ) %>%
    mutate(
      group = factor(group, levels = group_levels),
      sample_name = factor(sample_name, levels = sample_name_levels),
      metric_label = paste(sample_name, chain, sep = "_")
    )
}

calculate_gini <- function(x) {
  x <- sort(as.numeric(x))
  n <- length(x)
  total <- sum(x)
  if (n == 0 || total == 0) {
    return(NA_real_)
  }
  if (n == 1) {
    return(0)
  }
  sum((2 * seq_len(n) - n - 1) * x) / (n * total)
}

plot_tcr_diversity_clonality <- function(metrics, show_stats = TRUE) {
  if (nrow(metrics) == 0) {
    return(ggplot() +
      theme_void() +
      labs(title = "TCR diversity and clonality: no data"))
  }

  plot_data <- metrics %>%
    select(
      group, sample_name, chain, total_cells, shannon, gini
    ) %>%
    pivot_longer(
      cols = c(shannon, gini),
      names_to = "metric",
      values_to = "value"
    ) %>%
    mutate(
      metric = factor(
        metric,
        levels = c("shannon", "gini"),
        labels = c("Diversity: Shannon index", "Clonality: Gini coefficient")
      )
    )

  stat_data <- tibble()

  if (show_stats && "Naive" %in% as.character(plot_data$group)) {
    comparison_groups <- setdiff(as.character(na.omit(unique(plot_data$group))), "Naive")
    naive_comparisons <- map(comparison_groups, ~c("Naive", .x))

    facet_ranges <- plot_data %>%
      group_by(metric, chain) %>%
      summarise(
        y_max = max(value, na.rm = TRUE),
        y_min = min(value, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      mutate(
        y_step = case_when(
          is.finite(y_max - y_min) & y_max > y_min ~ (y_max - y_min) * 0.12,
          is.finite(y_max) & y_max != 0 ~ abs(y_max) * 0.12,
          TRUE ~ 0.1
        )
      )

    if (length(naive_comparisons) > 0) {
      stat_data <- compare_means(
        value ~ group,
        data = plot_data,
        group.by = c("metric", "chain"),
        comparisons = naive_comparisons,
        method = "wilcox.test",
        method.args = list(exact = FALSE)
      ) %>%
        left_join(facet_ranges, by = c("metric", "chain")) %>%
        group_by(metric, chain) %>%
        arrange(match(group2, group_levels), .by_group = TRUE) %>%
        mutate(
          y.position = y_max + y_step * row_number(),
          p_label = if_else(p.signif == "ns", paste0("p=", p.format), p.signif)
        ) %>%
        ungroup()
    }
  }

  plot <- ggplot(plot_data, aes(x = group, y = value, fill = group)) +
    geom_boxplot(width = 0.65, outlier.shape = NA, alpha = 0.75, color = "grey25") +
    geom_point(
      aes(color = group),
      position = position_jitter(width = 0.18, height = 0),
      size = 1.8,
      alpha = 0.9
    ) +
    facet_grid(metric ~ chain, scales = "free_y") +
    scale_fill_manual(values = color_group, drop = FALSE) +
    scale_color_manual(values = color_group, drop = FALSE) +
    theme_bw() +
    theme(
      legend.position = "bottom",
      axis.text.x = element_text(angle = 45, hjust = 1),
      strip.background = element_rect(fill = "white"),
      panel.grid.minor = element_blank()
    ) +
    labs(
      title = "TCR diversity and clonality across condition groups",
      subtitle = "Each point is one sample within a condition group and chain",
      x = NULL,
      y = NULL,
      fill = "Group",
      color = "Group"
    )

  if (nrow(stat_data) > 0) {
    plot <- plot +
      stat_pvalue_manual(
        stat_data,
        label = "p_label",
        tip.length = 0.01,
        size = 3
      )
  }

  plot
}

calculate_clone_size_classes <- function(annotation) {
  clone_counts <- annotation %>%
    filter(!is.na(cdr3), cdr3 != "") %>%
    distinct(barcode, group, sample_name, chain, cdr3) %>%
    count(group, sample_name, chain, cdr3, name = "clone_cells")

  if (nrow(clone_counts) == 0) {
    return(tibble())
  }

  clone_counts %>%
    group_by(group, sample_name, chain) %>%
    mutate(
      total_cells = sum(clone_cells),
      clone_fraction = clone_cells / total_cells,
      clone_size_class = case_when(
        clone_fraction < 0.001 ~ "Rare (<0.1%)",
        clone_fraction < 0.01 ~ "Large (0.1-1%)",
        TRUE ~ "Hyperexpanded (>1%)"
      )
    ) %>%
    ungroup() %>%
    mutate(
      group = factor(group, levels = group_levels),
      sample_name = factor(sample_name, levels = sample_name_levels),
      clone_size_class = factor(
        clone_size_class,
        levels = c("Rare (<0.1%)", "Large (0.1-1%)", "Hyperexpanded (>1%)")
      )
    )
}

plot_clone_size_stacked_bar <- function(clone_classes) {
  if (nrow(clone_classes) == 0) {
    return(ggplot() +
      theme_void() +
      labs(title = "TCR clone-size classes: no data"))
  }

  plot_data <- clone_classes %>%
    group_by(group, chain, clone_size_class) %>%
    summarise(clone_cells = sum(clone_cells), .groups = "drop_last") %>%
    mutate(percent = clone_cells / sum(clone_cells) * 100) %>%
    ungroup()

  clone_size_colors <- c(
    "Rare (<0.1%)" = "#A6CEE3",
    "Large (0.1-1%)" = "#FDBF6F",
    "Hyperexpanded (>1%)" = "#E31A1C"
  )

  ggplot(plot_data, aes(x = group, y = percent, fill = clone_size_class)) +
    geom_col(width = 0.72, color = "white", linewidth = 0.25) +
    facet_wrap(~chain, nrow = 1) +
    scale_fill_manual(values = clone_size_colors, drop = FALSE) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.03))) +
    theme_bw() +
    theme(
      legend.position = "bottom",
      strip.background = element_rect(fill = "white"),
      panel.grid.minor = element_blank()
    ) +
    labs(
      title = "TCR clone-size distribution across condition groups",
      x = NULL,
      y = "Cells in clone-size class (%)",
      fill = "Clone-size class"
    )
}

make_sequence_logo <- function(annotation, chain_name, cdr3_length, cell_types, title) {
  cdr3_values <- annotation %>%
    mutate(cdr3_length = nchar(cdr3)) %>%
    filter(
      chain == chain_name,
      cdr3_length == !!cdr3_length,
      cell_type %in% cell_types
    ) %>%
    pull(cdr3)

  if (length(cdr3_values) == 0) {
    return(ggplot() +
      theme_void() +
      labs(title = paste(title, "(no sequences)")))
  }

  ggseqlogo(cdr3_values,
    method = "bits",
    seq_type = "aa"
  ) +
    ggtitle(title)
}

make_logo_grid <- function(annotation, logo_set) {
  logo_plots <- imap(logo_set$groups, function(cell_types, group_name) {
    cell_type_title <- paste(cell_types, collapse = ", ")
    make_sequence_logo(
      annotation,
      logo_set$chain,
      logo_set$length,
      cell_types,
      paste0(logo_set$title, " Sequence Logo\nCell type: ", cell_type_title)
    )
  })
  wrap_plots(logo_plots, ncol = 1) +
    plot_layout(guides = "collect") &
    theme(legend.position = "bottom")
}

get_logo_grid_height <- function(logo_set) {
  max(logo_min_height, length(logo_set$groups) * logo_height_per_group)
}

# 1. Read final cell-type object and included productive VDJ annotations.
all_seurat_celltype <- readRDS(seurat_celltype_rds)
all_annotation_included <- readRDS(included_annotation_rds)
normalised <- normalise_metadata_levels(all_seurat_celltype, all_annotation_included)
all_seurat_celltype <- normalised$seurat
all_annotation_included <- normalised$annotation
rm(normalised)

# 2. Inspect duplicated TRD CDR3 calls within single cells.
duplicated_trd_cdr3 <- all_annotation_included %>%
  filter(chain == "TRD") %>%
  group_by(barcode, cdr3) %>%
  tally() %>%
  arrange(desc(n))
duplicated_trd_cdr3 %>%
  filter(n > 1)

# 3. Plot top CDR3 clone frequencies by sample.
plot_top_cdr3_samples <- plot_top_cdr3_by_sample(all_annotation_included)
save_plot(plot_top_cdr3_samples, "top200_cdr3_samples", 12, 8, overwrite = force_cdr3_stat_plot)

# 4. Mark and inspect the extreme MSH2-expanded TRD CDR3 before exclusion.
expanded_cdr3_barcodes <- all_annotation_included %>%
  filter(cdr3 %in% expanded_cdr3_to_exclude) %>%
  pull(barcode)
umap_mark <- plot_marked_cdr3_umap(all_seurat_celltype, expanded_cdr3_barcodes)
cluster_distr_mark <- plot_marked_cdr3_distribution(all_seurat_celltype, expanded_cdr3_barcodes)
save_plot(umap_mark, "umap_mark", 8, 6, overwrite = force_cdr3_stat_plot)
save_plot(cluster_distr_mark, "cluster_distr_mark", 8, 6, overwrite = force_cdr3_stat_plot)

# 5. Exclude the extreme expanded clone and replot CDR3 frequencies.
all_annotation_included_excluded <- all_annotation_included %>%
  filter(!cdr3 %in% expanded_cdr3_to_exclude)
plot_top_cdr3_samples_excluded <- plot_top_cdr3_by_sample(
  all_annotation_included_excluded,
  title = "Top CDR3 clone frequencies after expanded-clone exclusion"
)
save_plot(plot_top_cdr3_samples_excluded, "top200_cdr3_samples_excluded", 12, 8, overwrite = force_cdr3_stat_plot)

plot_top_cdr3_cell_type <- plot_top_cdr3_by_cell_type(all_annotation_included_excluded)
save_plot(plot_top_cdr3_cell_type, "top200_cdr3_cell_type_excluded", 12, 8, overwrite = force_cdr3_stat_plot)

# 6. Summarize TCR diversity, clonality, and clone-size classes by condition group.
tcr_diversity_metrics <- calculate_tcr_diversity_metrics(all_annotation_included_excluded)
write_csv_if_missing(
  tcr_diversity_metrics,
  tcr_diversity_metrics_csv,
  "TCR diversity and clonality metrics CSV",
  overwrite = force_cdr3_stat
)
tcr_diversity_clonality_boxplot <- plot_tcr_diversity_clonality(
  tcr_diversity_metrics,
  show_stats = FALSE
)
save_plot(
  tcr_diversity_clonality_boxplot,
  "TCR_diversity_clonality_boxplot",
  13,
  7,
  overwrite = force_cdr3_stat_plot
)
tcr_diversity_clonality_boxplot_with_stats <- plot_tcr_diversity_clonality(
  tcr_diversity_metrics,
  show_stats = TRUE
)
save_plot(
  tcr_diversity_clonality_boxplot_with_stats,
  "TCR_diversity_clonality_boxplot_with_stats",
  13,
  7,
  overwrite = force_cdr3_stat_plot
)

tcr_clone_size_classes <- calculate_clone_size_classes(all_annotation_included_excluded)
write_csv_if_missing(
  tcr_clone_size_classes,
  tcr_clone_size_classes_csv,
  "TCR clone-size class CSV",
  overwrite = force_cdr3_stat
)
tcr_clone_size_stacked_bar <- plot_clone_size_stacked_bar(tcr_clone_size_classes)
save_plot(
  tcr_clone_size_stacked_bar,
  "TCR_clone_size_class_stacked_bar",
  10,
  5.5,
  overwrite = force_cdr3_stat_plot
)

# 7. Plot CDR3 length distributions by chain and cell type.
cdr3_length_trg_lineplot <- plot_cdr3_length_distribution(
  all_annotation_included_excluded,
  "TRG"
)
save_plot(
  cdr3_length_trg_lineplot,
  "VDJ_CDR3_length_TRG_lineplot",
  cdr3_length_a4_width,
  cdr3_length_a4_height,
  overwrite = force_cdr3_stat_plot
)

cdr3_length_trd_lineplot <- plot_cdr3_length_distribution(
  all_annotation_included_excluded,
  "TRD"
)
save_plot(
  cdr3_length_trd_lineplot,
  "VDJ_CDR3_length_TRD_lineplot",
  cdr3_length_a4_width,
  cdr3_length_a4_height,
  overwrite = force_cdr3_stat_plot
)

# 8. Generate sequence logos for selected CDR3 length peaks.
for (logo_set in cdr3_logo_sets) {
  logo_grid <- make_logo_grid(all_annotation_included_excluded, logo_set)
  save_plot(
    logo_grid,
    logo_set$prefix,
    logo_grid_width,
    get_logo_grid_height(logo_set),
    overwrite = force_cdr3_stat_plot
  )
}
