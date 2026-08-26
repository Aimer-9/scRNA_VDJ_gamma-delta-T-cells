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
figure_dir <- file.path("figures", "14_Vd1Vd2_pairwise")
table_dir <- "table"

seurat_celltype_rds <- file.path(rds_dir, "all_seurat_celltype.rds")
de_markers_csv <- file.path(table_dir, "vd1_vd2_pairwise_de_markers.csv")
marker_summary_csv <- file.path(table_dir, "vd1_vd2_pairwise_marker_summary.csv")
gene_expression_by_sample_csv <- file.path(table_dir, "vd1_vd2_pairwise_gene_expression_by_sample.csv")
hallmark_delta_csv <- file.path(table_dir, "vd1_vd2_pairwise_hallmark_delta.csv")
hallmark_celltype_stats_csv <- file.path(table_dir, "vd1_vd2_pairwise_hallmark_celltype_stats.csv")

force_vd1_vd2_pairwise <- FALSE
force_vd1_vd2_pairwise_plot <- TRUE

pairwise_comparisons <- list(
  naive_vd1_vs_effector_memory_vd2 = c("Naive Vd1", "Effector Memory Vd2"),
  effector_vd1_vs_pan_effector_vd2 = c("Effector Vd1", "PAN Effector Vd2")
)

marker_gene_panel <- c(
  "TRDV1", "TRDV2", "TRGV9", "CCR7", "SELL", "IL7R", "TCF7", "LEF1",
  "LTB", "MAL", "BACH2", "SOX4", "NKG7", "GNLY", "PRF1", "GZMA",
  "GZMB", "GZMH", "IFNG", "TNF", "CCL3", "CCL4", "CX3CR1", "FGFBP2",
  "TBX21", "EOMES", "ZEB2", "XBP1", "CEBPB", "CXCR3", "KLRD1", "KLRG1"
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

get_present_comparisons <- function(seurat_obj) {
  present_cell_types <- unique(as.character(seurat_obj$cell_type))
  keep <- vapply(pairwise_comparisons, function(pair) all(pair %in% present_cell_types), logical(1))
  if (any(!keep)) {
    message("[SKIP] Missing comparison(s): ", paste(names(pairwise_comparisons)[!keep], collapse = ", "))
  }
  pairwise_comparisons[keep]
}

make_comparison_subset <- function(seurat_obj, comparison) {
  selected_cells <- rownames(seurat_obj@meta.data)[as.character(seurat_obj$cell_type) %in% comparison]
  subset_obj <- subset(seurat_obj, cells = selected_cells)
  Idents(subset_obj) <- factor(as.character(subset_obj$cell_type), levels = comparison)
  subset_obj
}

run_pairwise_de <- function(seurat_obj, comparisons) {
  imap_dfr(comparisons, function(comparison, comparison_name) {
    comparison_obj <- make_comparison_subset(seurat_obj, comparison)
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
        abs_log2FC = abs(avg_log2FC)
      ) %>%
      relocate(comparison, ident_1, ident_2, gene)
  })
}

get_selected_genes <- function(de_markers, seurat_obj, top_n = 12) {
  if (nrow(de_markers) == 0 || !"p_val_adj" %in% colnames(de_markers)) {
    return(intersect(marker_gene_panel, rownames(seurat_obj)))
  }
  top_de_genes <- de_markers %>%
    filter(!is.na(p_val_adj)) %>%
    group_by(comparison) %>%
    slice_max(order_by = abs_log2FC, n = top_n, with_ties = FALSE) %>%
    ungroup() %>%
    pull(gene)
  unique(c(marker_gene_panel, top_de_genes)) %>%
    intersect(rownames(seurat_obj))
}

make_gene_expression_tables <- function(seurat_obj, genes, comparisons) {
  if (length(genes) == 0) {
    return(list(marker_summary = tibble(), sample_summary = tibble()))
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
    pivot_longer(
      cols = all_of(genes),
      names_to = "gene",
      values_to = "expression"
    ) %>%
    filter(as.character(cell_type) %in% unique(unlist(comparisons)))

  marker_summary <- expression_data %>%
    group_by(gene, cell_type) %>%
    summarise(
      mean_expression = mean(expression, na.rm = TRUE),
      median_expression = median(expression, na.rm = TRUE),
      percent_expressing = mean(expression > 0, na.rm = TRUE) * 100,
      n_cells = n(),
      .groups = "drop"
    )

  sample_summary <- expression_data %>%
    group_by(gene, group, sample_name, cell_type) %>%
    summarise(
      mean_expression = mean(expression, na.rm = TRUE),
      median_expression = median(expression, na.rm = TRUE),
      percent_expressing = mean(expression > 0, na.rm = TRUE) * 100,
      n_cells = n(),
      .groups = "drop"
    )

  list(marker_summary = marker_summary, sample_summary = sample_summary)
}

calculate_hallmark_delta <- function(seurat_obj, comparisons) {
  seurat_obj <- GeneSetAnalysis(seurat_obj, genesets = hall50$human)
  scores <- seurat_obj@misc$AUCell$genesets
  if (is.null(scores)) {
    stop("Missing AUCell Hallmark scores after GeneSetAnalysis().", call. = FALSE)
  }

  metadata <- seurat_obj@meta.data %>%
    rownames_to_column("cell_id") %>%
    select(cell_id, group, sample_name, cell_type)
  common_cells <- intersect(colnames(scores), metadata$cell_id)
  scores <- scores[, common_cells, drop = FALSE]
  metadata <- metadata %>% filter(cell_id %in% common_cells)

  celltype_stats <- CalcStats(scores, f = metadata$cell_type)
  write_csv_if_missing(
    as.data.frame(celltype_stats) %>% rownames_to_column("pathway"),
    hallmark_celltype_stats_csv,
    "Vd1/Vd2 pairwise Hallmark cell-type stats CSV",
    overwrite = force_vd1_vd2_pairwise
  )

  delta_table <- imap_dfr(comparisons, function(comparison, comparison_name) {
    ident_1_cells <- metadata$cell_id[as.character(metadata$cell_type) == comparison[1]]
    ident_2_cells <- metadata$cell_id[as.character(metadata$cell_type) == comparison[2]]
    if (length(ident_1_cells) == 0 || length(ident_2_cells) == 0) {
      return(tibble())
    }
    tibble(
      comparison = comparison_name,
      ident_1 = comparison[1],
      ident_2 = comparison[2],
      pathway = rownames(scores),
      mean_ident_1 = rowMeans(scores[, ident_1_cells, drop = FALSE], na.rm = TRUE),
      mean_ident_2 = rowMeans(scores[, ident_2_cells, drop = FALSE], na.rm = TRUE)
    ) %>%
      mutate(
        delta = mean_ident_1 - mean_ident_2,
        abs_delta = abs(delta),
        higher_in = ifelse(delta >= 0, ident_1, ident_2)
      ) %>%
      arrange(comparison, dplyr::desc(abs_delta), pathway)
  })

  list(scores = scores, celltype_stats = celltype_stats, delta = delta_table)
}

plot_umap_highlight <- function(seurat_obj, comparisons) {
  plot_data <- seurat_obj@meta.data %>%
    rownames_to_column("cell_id") %>%
    mutate(
      comparison_cell_type = ifelse(
        as.character(cell_type) %in% unique(unlist(comparisons)),
        as.character(cell_type),
        "Other"
      ),
      comparison_cell_type = factor(comparison_cell_type, levels = c(unique(unlist(comparisons)), "Other"))
    )
  reduction_name <- if ("umap.unintegrated" %in% get_reduction_names(seurat_obj)) {
    "umap.unintegrated"
  } else {
    "umap"
  }
  umap_data <- Embeddings(seurat_obj, reduction_name) %>%
    as.data.frame() %>%
    rownames_to_column("cell_id") %>%
    as_tibble()
  colnames(umap_data)[2:3] <- c("UMAP_1", "UMAP_2")
  plot_data <- plot_data %>%
    left_join(umap_data, by = "cell_id")
  plot_cols <- c(color_celltype[unique(unlist(comparisons))], "Other" = "grey88")

  umap_plot <- ggplot(plot_data, aes(x = UMAP_1, y = UMAP_2, color = comparison_cell_type)) +
    geom_point(size = 0.35, alpha = 0.75) +
    scale_color_manual(values = plot_cols, drop = FALSE) +
    NoAxes()
  add_fixed_umap_coordinates(umap_plot) +
    labs(color = "Cell type")
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
      geom_text(
        data = label_data,
        aes(label = gene),
        size = 3,
        vjust = -0.5,
        show.legend = FALSE
      )
  }
  volcano_plot
}

plot_marker_dotplot <- function(seurat_obj, genes, comparisons) {
  plot_genes <- genes[genes %in% rownames(seurat_obj)]
  if (length(plot_genes) == 0) {
    return(NULL)
  }
  selected_cells <- rownames(seurat_obj@meta.data)[as.character(seurat_obj$cell_type) %in% unique(unlist(comparisons))]
  plot_obj <- subset(seurat_obj, cells = selected_cells)
  Idents(plot_obj) <- factor(as.character(plot_obj$cell_type), levels = unique(unlist(comparisons)))
  DotPlot(plot_obj, features = plot_genes, group.by = "cell_type") +
    coord_flip() +
    theme_test() +
    theme(axis.title = element_blank(), axis.text.x = element_text(angle = 35, hjust = 1))
}

plot_gene_expression_heatmap <- function(marker_summary, genes) {
  plot_data <- marker_summary %>%
    filter(gene %in% genes) %>%
    group_by(gene) %>%
    mutate(
      scaled_expression = as.numeric(scale(mean_expression)),
      scaled_expression = ifelse(is.na(scaled_expression), 0, scaled_expression)
    ) %>%
    ungroup() %>%
    mutate(
      cell_type = factor(as.character(cell_type), levels = unique(unlist(pairwise_comparisons))),
      gene = factor(gene, levels = unique(gene[order(-abs(scaled_expression))]))
    )
  if (nrow(plot_data) == 0) {
    return(NULL)
  }
  ggplot(plot_data, aes(x = cell_type, y = gene, fill = scaled_expression)) +
    geom_tile(color = "white", linewidth = 0.25) +
    scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0) +
    theme_test() +
    theme(axis.title = element_blank(), axis.text.x = element_text(angle = 35, hjust = 1)) +
    labs(fill = "Row-scaled\nmean expression")
}

plot_gene_expression_by_sample <- function(sample_summary, genes) {
  plot_data <- sample_summary %>%
    filter(gene %in% genes) %>%
    mutate(cell_type = factor(as.character(cell_type), levels = unique(unlist(pairwise_comparisons))))
  if (nrow(plot_data) == 0) {
    return(NULL)
  }
  ggplot(plot_data, aes(x = cell_type, y = mean_expression, color = cell_type)) +
    geom_boxplot(outlier.shape = NA, width = 0.55) +
    geom_point(aes(shape = group), position = position_jitter(width = 0.08), size = 1.6, alpha = 0.85) +
    facet_wrap(~gene, scales = "free_y", ncol = 4) +
    scale_color_manual(values = color_celltype, drop = FALSE) +
    theme_test() +
    theme(
      axis.title.x = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      strip.background = element_blank()
    ) +
    labs(y = "Mean expression by sample", color = "Cell type", shape = "Group")
}

plot_hallmark_delta <- function(delta_table, comparison_name, top_n = 12) {
  plot_data <- delta_table %>%
    filter(comparison == comparison_name) %>%
    slice_max(order_by = abs_delta, n = top_n, with_ties = FALSE) %>%
    mutate(
      pathway_label = str_replace(pathway, "^HALLMARK_", ""),
      pathway_label = str_replace_all(pathway_label, "_", " "),
      pathway_label = factor(pathway_label, levels = pathway_label[order(delta)])
    )
  if (nrow(plot_data) == 0) {
    return(NULL)
  }
  ggplot(plot_data, aes(x = delta, y = pathway_label, color = higher_in)) +
    geom_segment(aes(x = 0, xend = delta, yend = pathway_label), linewidth = 0.8) +
    geom_point(size = 2.4) +
    geom_vline(xintercept = 0, color = "grey45", linewidth = 0.35) +
    scale_color_manual(values = color_celltype, drop = FALSE) +
    theme_test() +
    theme(axis.title.y = element_blank()) +
    labs(x = "Mean AUCell delta", color = "Higher in", title = comparison_name)
}

all_seurat_celltype <- read_rds_checked(seurat_celltype_rds, "annotated Seurat RDS")
all_seurat_celltype <- normalise_metadata_levels(all_seurat_celltype)
all_seurat_celltype <- prepare_expression_object(all_seurat_celltype)

present_comparisons <- get_present_comparisons(all_seurat_celltype)
if (length(present_comparisons) == 0) {
  stop("None of the configured Vd1/Vd2 pairwise comparisons are present.", call. = FALSE)
}

de_markers <- run_pairwise_de(all_seurat_celltype, present_comparisons)
write_csv_if_missing(
  de_markers,
  de_markers_csv,
  "Vd1/Vd2 pairwise DE marker CSV",
  overwrite = force_vd1_vd2_pairwise
)

selected_genes <- get_selected_genes(de_markers, all_seurat_celltype, top_n = 10)
gene_tables <- make_gene_expression_tables(all_seurat_celltype, selected_genes, present_comparisons)
write_csv_if_missing(
  gene_tables$marker_summary,
  marker_summary_csv,
  "Vd1/Vd2 pairwise marker summary CSV",
  overwrite = force_vd1_vd2_pairwise
)
write_csv_if_missing(
  gene_tables$sample_summary,
  gene_expression_by_sample_csv,
  "Vd1/Vd2 pairwise gene expression by sample CSV",
  overwrite = force_vd1_vd2_pairwise
)

hallmark_result <- calculate_hallmark_delta(all_seurat_celltype, present_comparisons)
write_csv_if_missing(
  hallmark_result$delta,
  hallmark_delta_csv,
  "Vd1/Vd2 pairwise Hallmark delta CSV",
  overwrite = force_vd1_vd2_pairwise
)

save_plot(
  plot_umap_highlight(all_seurat_celltype, present_comparisons),
  "vd1_vd2_pairwise_umap_highlight",
  8,
  6,
  overwrite = force_vd1_vd2_pairwise_plot
)

marker_dotplot <- plot_marker_dotplot(all_seurat_celltype, selected_genes, present_comparisons)
if (!is.null(marker_dotplot)) {
  save_plot(marker_dotplot, "vd1_vd2_pairwise_marker_dotplot", 10, 8, overwrite = force_vd1_vd2_pairwise_plot)
}

gene_expression_heatmap <- plot_gene_expression_heatmap(gene_tables$marker_summary, selected_genes)
if (!is.null(gene_expression_heatmap)) {
  save_plot(
    gene_expression_heatmap,
    "vd1_vd2_pairwise_gene_expression_heatmap",
    8,
    max(5, length(unique(gene_tables$marker_summary$gene)) * 0.24),
    overwrite = force_vd1_vd2_pairwise_plot
  )
}

gene_expression_sample_plot <- plot_gene_expression_by_sample(gene_tables$sample_summary, head(selected_genes, 24))
if (!is.null(gene_expression_sample_plot)) {
  save_plot(
    gene_expression_sample_plot,
    "vd1_vd2_pairwise_gene_expression_by_sample",
    12,
    8,
    overwrite = force_vd1_vd2_pairwise_plot
  )
}

iwalk(present_comparisons, function(comparison, comparison_name) {
  volcano_plot <- plot_de_volcano(de_markers, comparison_name)
  if (!is.null(volcano_plot)) {
    save_plot(volcano_plot, paste0("volcano_", comparison_name), 7, 6, overwrite = force_vd1_vd2_pairwise_plot)
  }
  hallmark_delta_plot <- plot_hallmark_delta(hallmark_result$delta, comparison_name)
  if (!is.null(hallmark_delta_plot)) {
    save_plot(hallmark_delta_plot, paste0("hallmark_delta_", comparison_name), 8, 6, overwrite = force_vd1_vd2_pairwise_plot)
  }
})
