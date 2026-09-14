# R version 4.5.2 (2025-10-31)
rm(list = ls())
setwd("/path/to/project")
library(Seurat)
library(SeuratExtend)
library(SeuratWrappers)
library(monocle3)
library(tidyverse)
library(patchwork)
options(
  tibble.width = Inf,
  print.width = Inf,
  max.print = 200
)

rds_dir <- "rds"
figure_dir <- file.path("figures", "13_Vd2_pseudotime")
table_dir <- "table"

seurat_celltype_rds <- file.path(rds_dir, "all_seurat_celltype.rds")
paired_cdr3_rds <- file.path(rds_dir, "barcode_trgd_paired.rds")
pair_rank_metadata_rds <- file.path(rds_dir, "trdg_pair_rank_metadata.rds")
toprank_cells_rds <- file.path(rds_dir, "all_seurat_celltype_toprank_cells.rds")
toprank_vd2_pairs_csv <- file.path(table_dir, "toprank_vd2_trdg_pairs_for_monocle3.csv")
vd2_pseudotime_metadata_csv <- file.path(table_dir, "vd2_state_pseudotime_metadata.csv")
vd2_pseudotime_summary_csv <- file.path(table_dir, "vd2_state_pseudotime_summary.csv")
vd2_gene_expression_summary_csv <- file.path(table_dir, "vd2_state_gene_expression_summary.csv")
vd2_gene_expression_by_sample_csv <- file.path(table_dir, "vd2_state_gene_expression_by_sample.csv")
vd2_pseudotime_gene_correlation_csv <- file.path(table_dir, "vd2_state_pseudotime_gene_correlation.csv")

vd2_trajectory_states <- c(
  "Effector Memory Vd2",
  "Pre-activated Vd2",
  "ZOL Effector Vd2",
  "ZOL FOXP3+ Vd2",
  "PAN Effector Vd2"
)
trajectory_root_cell_type <- "Effector Memory Vd2"
top_rank_n <- 10
toprank_target_cells <- 5000
excluded_shared_pair_groups <- c("AB3")
required_top_clone_groups <- c("Naive", "ZOL", "PAN")
marker_gene_panel <- c(
  "TRDV2", "TRGV9", "NKG7", "GNLY", "PRF1", "GZMA", "GZMB", "GZMH",
  "IFNG", "TNF", "CCL3", "CCL4", "CD40LG", "CD70", "IL2RA", "ICOS",
  "FOXP3", "IL2RB", "PDCD1", "CTLA4", "LAG3", "HAVCR2", "TIGIT",
  "IRF1", "NR3C1", "GABPB1", "NFATC3", "XBP1", "CEBPB", "TFDP1"
)
force_vd2_pseudotime <- FALSE
force_vd2_pseudotime_plot <- TRUE

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

get_single_pair_records <- function(paired_cdr3) {
  single_barcodes <- paired_cdr3 %>%
    distinct(barcode, TRD, TRG) %>%
    group_by(barcode) %>%
    tally() %>%
    filter(n == 1) %>%
    pull(barcode)
  paired_cdr3 %>%
    filter(barcode %in% single_barcodes) %>%
    distinct(barcode, TRD, TRG, .keep_all = TRUE)
}

rank_paired_clones <- function(paired_single) {
  pair_rank <- paired_single %>%
    group_by(TRD, TRG) %>%
    summarise(
      n = n_distinct(barcode),
      clone_groups = list(unique(as.character(group))),
      .groups = "drop"
    ) %>%
    filter(map_lgl(clone_groups, ~all(required_top_clone_groups %in% .x))) %>%
    select(-clone_groups) %>%
    arrange(dplyr::desc(n), TRD, TRG) %>%
    mutate(rank = row_number())
  paired_single %>%
    left_join(pair_rank, by = c("TRD", "TRG"))
}

make_pair_rank_metadata <- function(paired_ranked) {
  paired_ranked %>%
    transmute(
      barcode,
      TRD,
      TRG,
      rank,
      TRDG_pair_n = n,
      rank_top10 = ifelse(
        is.na(rank),
        "other_pair",
        ifelse(rank <= top_rank_n, paste0("rank_", rank, "_pair"), "other_pair")
      )
    ) %>%
    distinct(barcode, .keep_all = TRUE) %>%
    mutate(
      rank_top10 = factor(
        rank_top10,
        levels = c(paste0("rank_", seq_len(top_rank_n), "_pair"), "other_pair")
      )
    )
}

read_or_make_pair_rank_metadata <- function(paired_ranked) {
  if (output_files_exist(pair_rank_metadata_rds) && !force_vd2_pseudotime) {
    message("[SKIP] Existing paired clone rank metadata RDS: ", pair_rank_metadata_rds)
    return(read_rds_checked(pair_rank_metadata_rds, "paired clone rank metadata RDS"))
  }
  if (output_files_exist(pair_rank_metadata_rds) && force_vd2_pseudotime) {
    message("[FORCE] Rebuilding paired clone rank metadata: ", pair_rank_metadata_rds)
  }
  pair_rank_metadata <- make_pair_rank_metadata(paired_ranked)
  save_rds_if_missing(
    pair_rank_metadata,
    pair_rank_metadata_rds,
    "paired clone rank metadata RDS",
    overwrite = force_vd2_pseudotime
  )
  pair_rank_metadata
}

add_pair_rank_metadata <- function(seurat_obj, pair_rank_metadata) {
  columns_to_remove <- intersect(
    c(
      "TRD", "TRG", "rank", "rank_top10", "TRDG_pair_n",
      "TRD_rank", "TRG_rank",
      "TRD_clone_n", "TRG_clone_n",
      "TRD_top10_clone", "TRG_top10_clone"
    ),
    colnames(seurat_obj@meta.data)
  )
  if (length(columns_to_remove) > 0) {
    seurat_obj@meta.data[, columns_to_remove] <- NULL
  }
  pair_metadata <- seurat_obj@meta.data %>%
    rownames_to_column("cell_id") %>%
    left_join(
      pair_rank_metadata %>% select(barcode, TRD, TRG, rank, TRDG_pair_n, rank_top10),
      by = "barcode"
    ) %>%
    column_to_rownames("cell_id")
  seurat_obj@meta.data <- pair_metadata
  seurat_obj$rank_top10 <- as.character(seurat_obj$rank_top10)
  seurat_obj$rank_top10[is.na(seurat_obj$rank_top10)] <- "other_pair"
  seurat_obj$rank_top10 <- factor(
    seurat_obj$rank_top10,
    levels = c(paste0("rank_", seq_len(top_rank_n), "_pair"), "other_pair")
  )
  seurat_obj
}

choose_top_vd2_pairs <- function(paired_ranked, target_cells = toprank_target_cells) {
  top_pairs <- paired_ranked %>%
    filter(
      !is.na(rank),
      !as.character(group) %in% excluded_shared_pair_groups,
      str_detect(as.character(cell_type), "Vd2")
    ) %>%
    group_by(TRD, TRG) %>%
    summarise(
      n_cells = n_distinct(barcode),
      n_group = n_distinct(group),
      n_sample = n_distinct(sample_name),
      rank = min(rank),
      .groups = "drop"
    ) %>%
    arrange(dplyr::desc(n_cells), rank, TRD, TRG) %>%
    mutate(
      top_pair_rank = row_number(),
      cumulative_cells = cumsum(n_cells)
    )
  if (nrow(top_pairs) == 0) {
    stop("No eligible single-pair Vd2 TRD/TRG clones found for pseudotime analysis.", call. = FALSE)
  }
  selected_n <- which(top_pairs$cumulative_cells >= target_cells)[1]
  if (is.na(selected_n)) {
    selected_n <- nrow(top_pairs)
  }
  top_pairs %>% slice_head(n = selected_n)
}

make_toprank_cell_metadata <- function(seurat_obj, top_pairs) {
  selected_pair_ids <- paste(top_pairs$TRD, top_pairs$TRG, sep = "||")
  seurat_obj@meta.data %>%
    rownames_to_column("cell_id") %>%
    mutate(pair_id = paste(TRD, TRG, sep = "||")) %>%
    filter(
      !as.character(group) %in% excluded_shared_pair_groups,
      str_detect(as.character(cell_type), "Vd2"),
      as.character(cell_type) %in% vd2_trajectory_states,
      pair_id %in% selected_pair_ids
    )
}

read_or_make_toprank_cell_metadata <- function(seurat_obj, paired_ranked) {
  if (output_files_exist(toprank_cells_rds) && !force_vd2_pseudotime) {
    message("[SKIP] Existing top-pair Vd2 cell metadata RDS: ", toprank_cells_rds)
    return(read_rds_checked(toprank_cells_rds, "top-pair Vd2 cell metadata RDS"))
  }
  if (output_files_exist(toprank_cells_rds) && force_vd2_pseudotime) {
    message("[FORCE] Rebuilding top-pair Vd2 cell metadata: ", toprank_cells_rds)
  } else {
    message("[RUN] Building top-pair Vd2 cell metadata: ", toprank_cells_rds)
  }
  toprank_vd2_pairs <- choose_top_vd2_pairs(paired_ranked, target_cells = toprank_target_cells)
  write_csv_if_missing(
    toprank_vd2_pairs,
    toprank_vd2_pairs_csv,
    "top-pair Vd2 TRD/TRG table",
    overwrite = force_vd2_pseudotime
  )
  message(
    "Selected ", nrow(toprank_vd2_pairs), " Vd2 TRD/TRG pairs covering ",
    max(toprank_vd2_pairs$cumulative_cells), " single-pair cells for pseudotime."
  )
  toprank_cell_metadata <- make_toprank_cell_metadata(seurat_obj, toprank_vd2_pairs)
  save_rds_if_missing(
    toprank_cell_metadata,
    toprank_cells_rds,
    "top-pair Vd2 cell metadata RDS",
    overwrite = force_vd2_pseudotime
  )
  toprank_cell_metadata
}

subset_from_cell_metadata <- function(seurat_obj, cell_metadata) {
  selected_cells <- intersect(as.character(cell_metadata$cell_id), colnames(seurat_obj))
  if (length(selected_cells) == 0) {
    stop("No cells from top-pair Vd2 cell metadata are present in the Seurat object.", call. = FALSE)
  }
  subset(seurat_obj, cells = selected_cells)
}

prepare_monocle_object <- function(seurat_obj) {
  if ("umap.unintegrated" %in% get_reduction_names(seurat_obj)) {
    seurat_obj@reductions$umap <- seurat_obj@reductions$umap.unintegrated
  }
  seurat_obj$sample <- seurat_obj$sample_name
  as.cell_data_set(seurat_obj)
}

get_earliest_principal_node <- function(cds, cell_type) {
  cell_ids <- which(colData(cds)$cell_type == cell_type)
  if (length(cell_ids) == 0) {
    fallback_cell_type <- names(sort(table(colData(cds)$cell_type), decreasing = TRUE))[1]
    message(
      "Root cell type '", cell_type, "' not found in Monocle3 subset; using '",
      fallback_cell_type, "' instead."
    )
    cell_ids <- which(colData(cds)$cell_type == fallback_cell_type)
  }
  closest_vertex <- cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
  igraph::V(principal_graph(cds)[["UMAP"]])$name[
    as.numeric(names(which.max(table(closest_vertex[cell_ids, ]))))
  ]
}

run_vd2_pseudotime <- function(seurat_obj) {
  cds <- prepare_monocle_object(seurat_obj)
  cds <- cluster_cells(cds, reduction_method = "UMAP")
  cds <- learn_graph(cds, use_partition = TRUE)
  root_node <- get_earliest_principal_node(cds, trajectory_root_cell_type)
  order_cells(cds, root_pr_nodes = root_node)
}

make_pseudotime_metadata <- function(cds) {
  umap <- reducedDims(cds)$UMAP
  if (is.null(umap)) {
    stop("Missing UMAP reducedDims(cds)$UMAP for pseudotime metadata.", call. = FALSE)
  }
  colData(cds) %>%
    as.data.frame() %>%
    rownames_to_column("cell_id") %>%
    as_tibble() %>%
    mutate(
      UMAP_1 = umap[cell_id, 1],
      UMAP_2 = umap[cell_id, 2],
      pseudotime = monocle3::pseudotime(cds)[cell_id],
      group = factor(as.character(group), levels = group_levels),
      sample_name = factor(as.character(sample_name), levels = sample_name_levels),
      cell_type = factor(as.character(cell_type), levels = vd2_trajectory_states)
    ) %>%
    filter(is.finite(pseudotime))
}

summarise_pseudotime <- function(pseudotime_metadata) {
  pseudotime_metadata %>%
    group_by(group, sample_name, cell_type) %>%
    summarise(
      n_cells = n(),
      mean_pseudotime = mean(pseudotime, na.rm = TRUE),
      median_pseudotime = median(pseudotime, na.rm = TRUE),
      q25_pseudotime = quantile(pseudotime, 0.25, na.rm = TRUE),
      q75_pseudotime = quantile(pseudotime, 0.75, na.rm = TRUE),
      .groups = "drop"
    )
}

plot_pseudotime_umap <- function(pseudotime_metadata) {
  umap_plot <- ggplot(pseudotime_metadata, aes(x = UMAP_1, y = UMAP_2, color = pseudotime)) +
    geom_point(size = 0.45, alpha = 0.85) +
    scale_color_viridis_c(option = "plasma", name = "Pseudotime") +
    NoAxes()
  add_fixed_umap_coordinates(umap_plot)
}

plot_celltype_umap <- function(pseudotime_metadata) {
  umap_plot <- ggplot(pseudotime_metadata, aes(x = UMAP_1, y = UMAP_2, color = cell_type)) +
    geom_point(size = 0.45, alpha = 0.85) +
    scale_color_manual(values = color_celltype, drop = FALSE) +
    NoAxes()
  add_fixed_umap_coordinates(umap_plot) +
    labs(color = "Vd2 state")
}

plot_group_umap <- function(pseudotime_metadata) {
  umap_plot <- ggplot(pseudotime_metadata, aes(x = UMAP_1, y = UMAP_2, color = group)) +
    geom_point(size = 0.45, alpha = 0.85) +
    scale_color_manual(values = color_group, drop = FALSE) +
    NoAxes()
  add_fixed_umap_coordinates(umap_plot) +
    labs(color = "Group")
}

plot_pseudotime_by_state <- function(pseudotime_metadata) {
  ggplot(pseudotime_metadata, aes(x = cell_type, y = pseudotime, fill = cell_type)) +
    geom_violin(scale = "width", trim = TRUE, alpha = 0.8) +
    geom_boxplot(width = 0.16, outlier.shape = NA, color = "grey25") +
    scale_fill_manual(values = color_celltype, drop = FALSE) +
    theme_test() +
    theme(
      axis.title.x = element_blank(),
      axis.text.x = element_text(angle = 35, hjust = 1),
      legend.position = "none"
    ) +
    labs(y = "Pseudotime")
}

plot_pseudotime_by_sample <- function(pseudotime_summary) {
  ggplot(pseudotime_summary, aes(x = sample_name, y = median_pseudotime, fill = cell_type)) +
    geom_col(position = position_dodge(width = 0.8), width = 0.75) +
    facet_grid(. ~ group, scales = "free_x", space = "free_x") +
    scale_fill_manual(values = color_celltype, drop = FALSE) +
    theme_test() +
    theme(
      axis.title.x = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1),
      strip.background = element_blank()
    ) +
    labs(y = "Median pseudotime", fill = "Vd2 state")
}

plot_pseudotime_density <- function(pseudotime_metadata) {
  ggplot(pseudotime_metadata, aes(x = pseudotime, color = cell_type, fill = cell_type)) +
    geom_density(alpha = 0.18, linewidth = 0.8) +
    facet_wrap(~cell_type, ncol = 1, scales = "free_y") +
    scale_color_manual(values = color_celltype, drop = FALSE) +
    scale_fill_manual(values = color_celltype, drop = FALSE) +
    theme_test() +
    theme(
      legend.position = "none",
      strip.background = element_blank()
    ) +
    labs(x = "Pseudotime", y = "Density")
}

plot_gene_pseudotime_trends <- function(seurat_obj, pseudotime_metadata, genes) {
  seurat_obj <- join_assay_layers_if_needed(seurat_obj)
  genes <- genes[genes %in% rownames(seurat_obj)]
  if (length(genes) == 0) {
    message("[SKIP] No configured marker genes available for pseudotime trend plot.")
    return(NULL)
  }
  genes <- head(genes, 18)
  expression_data <- FetchData(seurat_obj, vars = genes) %>%
    rownames_to_column("cell_id") %>%
    as_tibble() %>%
    filter(cell_id %in% pseudotime_metadata$cell_id) %>%
    pivot_longer(
      cols = all_of(genes),
      names_to = "gene",
      values_to = "expression"
    ) %>%
    left_join(
      pseudotime_metadata %>% select(cell_id, pseudotime, cell_type),
      by = "cell_id"
    )
  ggplot(expression_data, aes(x = pseudotime, y = expression, color = cell_type)) +
    geom_smooth(se = FALSE, method = "loess", span = 0.55, linewidth = 0.75) +
    facet_wrap(~gene, scales = "free_y", ncol = 3) +
    scale_color_manual(values = color_celltype, drop = FALSE) +
    theme_test() +
    theme(strip.background = element_blank()) +
    labs(x = "Pseudotime", y = "Expression", color = "Vd2 state")
}

make_gene_expression_tables <- function(seurat_obj, pseudotime_metadata, genes) {
  seurat_obj <- join_assay_layers_if_needed(seurat_obj)
  genes <- genes[genes %in% rownames(seurat_obj)]
  if (length(genes) == 0) {
    return(list(group_summary = tibble(), sample_summary = tibble(), cell_expression = tibble()))
  }
  expression_data <- FetchData(seurat_obj, vars = genes) %>%
    rownames_to_column("cell_id") %>%
    as_tibble() %>%
    filter(cell_id %in% pseudotime_metadata$cell_id) %>%
    pivot_longer(
      cols = all_of(genes),
      names_to = "gene",
      values_to = "expression"
    ) %>%
    left_join(
      pseudotime_metadata %>% select(cell_id, group, sample_name, cell_type, pseudotime),
      by = "cell_id"
    )
  group_summary <- expression_data %>%
    group_by(gene, cell_type) %>%
    summarise(
      mean_expression = mean(expression, na.rm = TRUE),
      median_expression = median(expression, na.rm = TRUE),
      percent_expressing = mean(expression > 0, na.rm = TRUE) * 100,
      mean_pseudotime = mean(pseudotime, na.rm = TRUE),
      n_cells = n(),
      .groups = "drop"
    )
  sample_summary <- expression_data %>%
    group_by(gene, group, sample_name, cell_type) %>%
    summarise(
      mean_expression = mean(expression, na.rm = TRUE),
      median_expression = median(expression, na.rm = TRUE),
      percent_expressing = mean(expression > 0, na.rm = TRUE) * 100,
      mean_pseudotime = mean(pseudotime, na.rm = TRUE),
      n_cells = n(),
      .groups = "drop"
    )
  list(
    group_summary = group_summary,
    sample_summary = sample_summary,
    cell_expression = expression_data
  )
}

calculate_pseudotime_gene_correlations <- function(expression_data) {
  if (nrow(expression_data) == 0) {
    return(tibble())
  }
  expression_data %>%
    group_by(gene) %>%
    summarise(
      pseudotime_spearman = suppressWarnings(cor(expression, pseudotime, method = "spearman", use = "pairwise.complete.obs")),
      pseudotime_pearson = suppressWarnings(cor(expression, pseudotime, method = "pearson", use = "pairwise.complete.obs")),
      n_cells = n(),
      .groups = "drop"
    ) %>%
    arrange(dplyr::desc(abs(pseudotime_spearman)), gene)
}

plot_gene_expression_heatmap <- function(expression_summary) {
  if (nrow(expression_summary) == 0) {
    message("[SKIP] No gene expression summary for heatmap.")
    return(NULL)
  }
  plot_data <- expression_summary %>%
    group_by(gene) %>%
    mutate(
      scaled_expression = as.numeric(scale(mean_expression)),
      scaled_expression = ifelse(is.na(scaled_expression), 0, scaled_expression)
    ) %>%
    ungroup() %>%
    mutate(
      cell_type = factor(as.character(cell_type), levels = vd2_trajectory_states),
      gene = factor(gene, levels = unique(gene[order(-abs(scaled_expression))]))
    )
  ggplot(plot_data, aes(x = cell_type, y = gene, fill = scaled_expression)) +
    geom_tile(color = "white", linewidth = 0.25) +
    scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0) +
    theme_test() +
    theme(
      axis.title = element_blank(),
      axis.text.x = element_text(angle = 35, hjust = 1)
    ) +
    labs(fill = "Row-scaled\nmean expression")
}

plot_gene_expression_by_sample <- function(sample_summary, correlation_table, top_n = 16) {
  if (nrow(sample_summary) == 0 || nrow(correlation_table) == 0) {
    message("[SKIP] No sample gene expression data for sample plot.")
    return(NULL)
  }
  selected_genes <- correlation_table %>%
    slice_max(order_by = abs(pseudotime_spearman), n = top_n, with_ties = FALSE) %>%
    pull(gene)
  plot_data <- sample_summary %>%
    filter(gene %in% selected_genes) %>%
    mutate(cell_type = factor(as.character(cell_type), levels = vd2_trajectory_states))
  ggplot(plot_data, aes(x = cell_type, y = mean_expression, color = cell_type)) +
    geom_boxplot(outlier.shape = NA, width = 0.55) +
    geom_point(aes(shape = group), position = position_jitter(width = 0.08), size = 1.7, alpha = 0.85) +
    facet_wrap(~gene, scales = "free_y", ncol = 4) +
    scale_color_manual(values = color_celltype, drop = FALSE) +
    theme_test() +
    theme(
      axis.title.x = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      strip.background = element_blank()
    ) +
    labs(y = "Mean expression by sample", color = "Vd2 state", shape = "Group")
}

plot_pseudotime_gene_correlation <- function(correlation_table, top_n = 20) {
  if (nrow(correlation_table) == 0) {
    message("[SKIP] No gene-pseudotime correlations to plot.")
    return(NULL)
  }
  plot_data <- correlation_table %>%
    slice_max(order_by = abs(pseudotime_spearman), n = top_n, with_ties = FALSE) %>%
    mutate(
      direction = ifelse(pseudotime_spearman >= 0, "Increases with pseudotime", "Decreases with pseudotime"),
      gene = factor(gene, levels = gene[order(pseudotime_spearman)])
    )
  ggplot(plot_data, aes(x = pseudotime_spearman, y = gene, fill = direction)) +
    geom_col(width = 0.75) +
    geom_vline(xintercept = 0, color = "grey40", linewidth = 0.4) +
    scale_fill_manual(values = c(
      "Increases with pseudotime" = "#B2182B",
      "Decreases with pseudotime" = "#2166AC"
    )) +
    theme_test() +
    theme(axis.title.y = element_blank()) +
    labs(x = "Spearman correlation with pseudotime", fill = "")
}

plot_pseudotime_expression_heatmap <- function(expression_data, correlation_table, top_n = 24, n_bins = 50) {
  if (nrow(expression_data) == 0 || nrow(correlation_table) == 0) {
    message("[SKIP] No expression/correlation data for pseudotime heatmap.")
    return(NULL)
  }
  selected_genes <- correlation_table %>%
    filter(is.finite(pseudotime_spearman)) %>%
    slice_max(order_by = abs(pseudotime_spearman), n = top_n, with_ties = FALSE) %>%
    pull(gene)
  if (length(selected_genes) == 0) {
    message("[SKIP] No finite gene-pseudotime correlations for pseudotime heatmap.")
    return(NULL)
  }

  cell_bins <- expression_data %>%
    distinct(cell_id, pseudotime, cell_type) %>%
    filter(is.finite(pseudotime)) %>%
    arrange(pseudotime) %>%
    mutate(
      pseudotime_bin = dplyr::ntile(pseudotime, n_bins)
    )

  binned_expression <- expression_data %>%
    filter(gene %in% selected_genes) %>%
    select(cell_id, gene, expression) %>%
    inner_join(
      cell_bins %>% select(cell_id, pseudotime, cell_type, pseudotime_bin),
      by = "cell_id"
    ) %>%
    mutate(
      pseudotime_bin = factor(pseudotime_bin, levels = seq_len(n_bins))
    )

  if (nrow(binned_expression) == 0) {
    message("[SKIP] No binned expression data for pseudotime heatmap.")
    return(NULL)
  }

  gene_order <- correlation_table %>%
    filter(gene %in% selected_genes) %>%
    arrange(pseudotime_spearman) %>%
    pull(gene)

  plot_data <- binned_expression %>%
    group_by(gene, pseudotime_bin) %>%
    summarise(
      mean_expression = mean(expression, na.rm = TRUE),
      mean_pseudotime = mean(pseudotime, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    group_by(gene) %>%
    mutate(
      scaled_expression = as.numeric(scale(mean_expression)),
      scaled_expression = ifelse(is.na(scaled_expression), 0, scaled_expression)
    ) %>%
    ungroup() %>%
    mutate(
      gene = factor(gene, levels = gene_order),
      pseudotime_bin = factor(pseudotime_bin, levels = seq_len(n_bins))
    ) %>%
    complete(
      gene = factor(gene_order, levels = gene_order),
      pseudotime_bin = factor(seq_len(n_bins), levels = seq_len(n_bins)),
      fill = list(mean_expression = 0, mean_pseudotime = NA_real_, scaled_expression = 0)
    )

  heatmap_matrix <- plot_data %>%
    select(gene, pseudotime_bin, scaled_expression) %>%
    pivot_wider(
      names_from = pseudotime_bin,
      values_from = scaled_expression,
      values_fn = list(scaled_expression = mean),
      values_fill = 0
    ) %>%
    arrange(gene) %>%
    as.data.frame()
  heatmap_gene_ids <- as.character(heatmap_matrix$gene)
  heatmap_matrix$gene <- NULL
  heatmap_matrix <- data.frame(lapply(heatmap_matrix, as.numeric), check.names = FALSE)
  rownames(heatmap_matrix) <- heatmap_gene_ids
  heatmap_matrix <- as.matrix(heatmap_matrix)
  missing_genes <- setdiff(gene_order, rownames(heatmap_matrix))
  if (length(missing_genes) > 0) {
    missing_gene_matrix <- matrix(
      0,
      nrow = length(missing_genes),
      ncol = ncol(heatmap_matrix),
      dimnames = list(missing_genes, colnames(heatmap_matrix))
    )
    heatmap_matrix <- rbind(heatmap_matrix, missing_gene_matrix)
  }
  missing_bins <- setdiff(as.character(seq_len(n_bins)), colnames(heatmap_matrix))
  if (length(missing_bins) > 0) {
    missing_bin_matrix <- matrix(
      0,
      nrow = nrow(heatmap_matrix),
      ncol = length(missing_bins),
      dimnames = list(rownames(heatmap_matrix), missing_bins)
    )
    heatmap_matrix <- cbind(heatmap_matrix, missing_bin_matrix)
  }
  heatmap_matrix <- heatmap_matrix[gene_order, as.character(seq_len(n_bins)), drop = FALSE]
  heatmap_matrix[heatmap_matrix > 2.5] <- 2.5
  heatmap_matrix[heatmap_matrix < -2.5] <- -2.5

  plot_data <- as.data.frame(heatmap_matrix) %>%
    rownames_to_column("gene") %>%
    as_tibble() %>%
    pivot_longer(
      cols = -gene,
      names_to = "pseudotime_bin",
      values_to = "scaled_expression"
    ) %>%
    mutate(
      gene = factor(gene, levels = gene_order),
      pseudotime_bin = as.integer(pseudotime_bin)
    )

  ggplot(plot_data, aes(x = pseudotime_bin, y = gene, fill = scaled_expression)) +
    geom_tile(color = NA) +
    scale_fill_gradient2(
      low = "#2166AC",
      mid = "white",
      high = "#B2182B",
      midpoint = 0,
      limits = c(-2.5, 2.5),
      oob = scales::squish
    ) +
    scale_x_continuous(expand = c(0, 0)) +
    theme_test() +
    theme(
      axis.title.y = element_blank(),
      axis.text.y = element_text(size = 7),
      panel.grid = element_blank()
    ) +
    labs(
      x = "Pseudotime bin",
      fill = "Row-scaled\nmean expression"
    )
}

all_seurat_celltype <- read_rds_checked(seurat_celltype_rds, "annotated Seurat RDS")
all_seurat_celltype <- normalise_metadata_levels(all_seurat_celltype)

barcode_trgd_paired <- read_rds_checked(paired_cdr3_rds, "paired TRD/TRG CDR3 RDS") %>%
  normalise_annotation_levels()
barcode_paired_single_all <- get_single_pair_records(barcode_trgd_paired)
barcode_paired_ranked <- rank_paired_clones(barcode_paired_single_all)
pair_rank_metadata <- read_or_make_pair_rank_metadata(barcode_paired_ranked)
all_seurat_celltype <- add_pair_rank_metadata(all_seurat_celltype, pair_rank_metadata)

toprank_cell_metadata <- read_or_make_toprank_cell_metadata(all_seurat_celltype, barcode_paired_ranked)
vd2_obj <- subset_from_cell_metadata(all_seurat_celltype, toprank_cell_metadata)
vd2_obj <- add_pair_rank_metadata(vd2_obj, pair_rank_metadata)

cell_type_values <- as.character(vd2_obj[[]]$cell_type)
available_states <- intersect(vd2_trajectory_states, unique(cell_type_values))
missing_states <- setdiff(vd2_trajectory_states, available_states)
if (length(missing_states) > 0) {
  message("[SKIP] Missing configured Vd2 state(s) in selected top-pair clone cells: ", paste(missing_states, collapse = ", "))
}
if (length(available_states) < 2) {
  stop("Need at least two configured Vd2 states in selected top-pair clone cells for pseudotime analysis.", call. = FALSE)
}

vd2_cells <- colnames(vd2_obj)[cell_type_values %in% available_states]
vd2_obj <- subset(vd2_obj, cells = vd2_cells)
if ("RNA" %in% get_assay_names(vd2_obj)) {
  DefaultAssay(vd2_obj) <- "RNA"
}
vd2_obj <- join_assay_layers_if_needed(vd2_obj)
vd2_obj$cell_type <- factor(as.character(vd2_obj$cell_type), levels = vd2_trajectory_states)
message("Vd2 pseudotime cells from selected paired clones: ", ncol(vd2_obj))

cds_vd2 <- run_vd2_pseudotime(vd2_obj)
pseudotime_metadata <- make_pseudotime_metadata(cds_vd2)
write_csv_if_missing(
  pseudotime_metadata,
  vd2_pseudotime_metadata_csv,
  "Vd2 state pseudotime metadata CSV",
  overwrite = force_vd2_pseudotime
)

pseudotime_summary <- summarise_pseudotime(pseudotime_metadata)
write_csv_if_missing(
  pseudotime_summary,
  vd2_pseudotime_summary_csv,
  "Vd2 state pseudotime summary CSV",
  overwrite = force_vd2_pseudotime
)

save_plot(
  plot_pseudotime_umap(pseudotime_metadata),
  "vd2_state_pseudotime_umap",
  7,
  5,
  overwrite = force_vd2_pseudotime_plot
)
save_plot(
  plot_celltype_umap(pseudotime_metadata),
  "vd2_state_celltype_umap",
  8,
  5,
  overwrite = force_vd2_pseudotime_plot
)
save_plot(
  plot_group_umap(pseudotime_metadata),
  "vd2_state_group_umap",
  8,
  5,
  overwrite = force_vd2_pseudotime_plot
)
save_plot(
  plot_pseudotime_by_state(pseudotime_metadata),
  "vd2_state_pseudotime_violin",
  8,
  5,
  overwrite = force_vd2_pseudotime_plot
)
save_plot(
  plot_pseudotime_by_sample(pseudotime_summary),
  "vd2_state_pseudotime_by_sample",
  12,
  5,
  overwrite = force_vd2_pseudotime_plot
)
save_plot(
  plot_pseudotime_density(pseudotime_metadata),
  "vd2_state_pseudotime_density",
  8,
  10,
  overwrite = force_vd2_pseudotime_plot
)

gene_trend_plot <- plot_gene_pseudotime_trends(vd2_obj, pseudotime_metadata, marker_gene_panel)
if (!is.null(gene_trend_plot)) {
  save_plot(
    gene_trend_plot,
    "vd2_state_gene_pseudotime_trends",
    11,
    12,
    overwrite = force_vd2_pseudotime_plot
  )
}

gene_expression_tables <- make_gene_expression_tables(vd2_obj, pseudotime_metadata, marker_gene_panel)
write_csv_if_missing(
  gene_expression_tables$group_summary,
  vd2_gene_expression_summary_csv,
  "Vd2 state gene expression summary CSV",
  overwrite = force_vd2_pseudotime
)
write_csv_if_missing(
  gene_expression_tables$sample_summary,
  vd2_gene_expression_by_sample_csv,
  "Vd2 state gene expression by sample CSV",
  overwrite = force_vd2_pseudotime
)

gene_correlations <- calculate_pseudotime_gene_correlations(gene_expression_tables$cell_expression)
write_csv_if_missing(
  gene_correlations,
  vd2_pseudotime_gene_correlation_csv,
  "Vd2 state pseudotime gene correlation CSV",
  overwrite = force_vd2_pseudotime
)

gene_expression_heatmap <- plot_gene_expression_heatmap(gene_expression_tables$group_summary)
if (!is.null(gene_expression_heatmap)) {
  save_plot(
    gene_expression_heatmap,
    "vd2_state_gene_expression_heatmap",
    8,
    max(5, length(unique(gene_expression_tables$group_summary$gene)) * 0.24),
    overwrite = force_vd2_pseudotime_plot
  )
}

gene_expression_sample_plot <- plot_gene_expression_by_sample(
  gene_expression_tables$sample_summary,
  gene_correlations,
  top_n = 16
)
if (!is.null(gene_expression_sample_plot)) {
  save_plot(
    gene_expression_sample_plot,
    "vd2_state_gene_expression_by_sample",
    12,
    8,
    overwrite = force_vd2_pseudotime_plot
  )
}

gene_correlation_plot <- plot_pseudotime_gene_correlation(gene_correlations, top_n = 20)
if (!is.null(gene_correlation_plot)) {
  save_plot(
    gene_correlation_plot,
    "vd2_state_pseudotime_gene_correlation",
    8,
    6,
    overwrite = force_vd2_pseudotime_plot
  )
}

pseudotime_expression_heatmap <- plot_pseudotime_expression_heatmap(
  gene_expression_tables$cell_expression,
  gene_correlations,
  top_n = 24,
  n_bins = 50
)
if (!is.null(pseudotime_expression_heatmap)) {
  save_plot(
    pseudotime_expression_heatmap,
    "vd2_state_pseudotime_expression_heatmap",
    9,
    7,
    overwrite = force_vd2_pseudotime_plot
  )
}
