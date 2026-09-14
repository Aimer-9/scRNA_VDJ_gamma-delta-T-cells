# R version 4.5.2 (2025-10-31)
rm(list = ls())
setwd("/path/to/project")
library(Seurat)
library(SeuratExtend)
library(tidyverse)
library(ggpubr)
library(patchwork)
library(ComplexHeatmap)
library(circlize)
library(grid)
options(
  max.print = 100,
  tibble.width = Inf,
  spe = "human"
)

rds_dir <- "rds"
figure_dir <- file.path("figures", "10_MSH2")
table_dir <- "table"

seurat_celltype_rds <- file.path(rds_dir, "all_seurat_celltype.rds")
included_annotation_rds <- file.path(rds_dir, "all_annotation_included.rds")
marker_cell_table <- file.path(table_dir, "msh2_marker_cdr3_cells.csv")
marker_summary_table <- file.path(table_dir, "msh2_marker_group_summary.csv")
hallmark_stats_table <- file.path(table_dir, "msh2_marker_hallmark_stats.csv")
hallmark_delta_table <- file.path(table_dir, "msh2_marker_hallmark_delta.csv")
effector_score_table <- file.path(table_dir, "msh2_marker_effector_score_summary.csv")
centroid_distance_table <- file.path(table_dir, "msh2_marker_umap_centroid_distance.csv")
centroid_distance_summary_table <- file.path(table_dir, "msh2_marker_umap_centroid_distance_summary.csv")
de_marker_table <- file.path(table_dir, "msh2_marker_de_markers.csv")
cdr3g_pair_summary_table <- file.path(table_dir, "msh2_marker_cdr3g_to_caldttfpigdrgytdklif_top10.csv")
cdr3d_pair_summary_table <- file.path(table_dir, "msh2_marker_top_cdr3g_to_cdr3d_top10.csv")

msh2_marker_cdr3 <- "CALDTTFPIGDRGYTDKLIF"
top_msh2_cdr3g_n <- 10
nonmarker_control_groups <- c("Blank", "MSH2")
marker_levels <- c("Naive", "Non-Marker", "Marker", "Effector")
marker_colors <- c(
  "Naive" = "#009E73",
  "Non-Marker" = "grey75",
  "Marker" = "#D62728",
  "Effector" = "#0072B2"
)
vd1_marker_gene_classes <- tribble(
  ~gene_class, ~gene, ~marker_dotplot, ~effector_score,
  "Cytotoxicity", "NKG7", FALSE, TRUE,
  "Cytotoxicity", "GNLY", FALSE, TRUE,
  "Cytotoxicity", "PRF1", FALSE, TRUE,
  "Cytotoxicity", "GZMA", TRUE, TRUE,
  "Cytotoxicity", "GZMB", TRUE, TRUE,
  "Cytotoxicity", "GZMH", TRUE, TRUE,
  "Inflammatory cytokine", "IFNG", TRUE, TRUE,
  "Inflammatory cytokine", "TNF", TRUE, TRUE,
  "Inflammatory cytokine", "CCL3", FALSE, TRUE,
  "Inflammatory cytokine", "CCL4", FALSE, TRUE,
  "Costimulation activation", "CD40LG", TRUE, TRUE,
  "Costimulation activation", "CD70", TRUE, TRUE,
  "Costimulation activation", "IL2RA", TRUE, TRUE,
  "Costimulation activation", "ICOS", FALSE, TRUE,
  "Checkpoint exhaustion", "CTLA4", TRUE, FALSE,
  "Checkpoint exhaustion", "PDCD1", TRUE, FALSE,
  "Checkpoint exhaustion", "LAG3", TRUE, FALSE,
  "Checkpoint exhaustion", "HAVCR2", TRUE, FALSE,
  "Checkpoint exhaustion", "TIGIT", TRUE, FALSE,
  "Migration NK effector", "CX3CR1", FALSE, TRUE,
  "Migration NK effector", "KLRD1", FALSE, TRUE,
  "Migration NK effector", "KLRG1", FALSE, TRUE
)
vd1_marker_gene_classes$gene_class <- factor(
  vd1_marker_gene_classes$gene_class,
  levels = unique(vd1_marker_gene_classes$gene_class)
)
vd1_marker_genes <- vd1_marker_gene_classes %>%
  filter(marker_dotplot) %>%
  pull(gene)
vd1_effector_activation_genes <- vd1_marker_gene_classes %>%
  filter(effector_score) %>%
  pull(gene)
vd1_hallmark_pathways <- c(
  "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
  "HALLMARK_P53_PATHWAY",
  "HALLMARK_WNT_BETA_CATENIN_SIGNALING",
  "HALLMARK_INTERFERON_ALPHA_RESPONSE",
  "HALLMARK_INTERFERON_GAMMA_RESPONSE",
  "HALLMARK_IL6_JAK_STAT3_SIGNALING",
  "HALLMARK_INFLAMMATORY_RESPONSE",
  "HALLMARK_MTORC1_SIGNALING",
  "HALLMARK_FATTY_ACID_METABOLISM",
  "HALLMARK_TGF_BETA_SIGNALING",
  "HALLMARK_APOPTOSIS",
  "HALLMARK_OXIDATIVE_PHOSPHORYLATION",
  "HALLMARK_MITOTIC_SPINDLE",
  "HALLMARK_DNA_REPAIR"
)
force_msh2 <- FALSE
force_msh2_plot <- TRUE

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

normalise_msh2_metadata_levels <- function(seurat_obj, annotation) {
  list(
    seurat = normalise_metadata_levels(seurat_obj),
    annotation = normalise_annotation_levels(annotation)
  )
}

get_marker_barcodes <- function(annotation, marker_cdr3) {
  annotation %>%
    dplyr::filter(cdr3 == marker_cdr3) %>%
    dplyr::distinct(barcode, sample_name, group, cell_type, chain, cdr3)
}

make_marker_cdr3g_pairs <- function(annotation, marker_cells, marker_cdr3) {
  marker_barcodes <- unique(marker_cells$barcode)
  if (length(marker_barcodes) == 0) {
    return(tibble(
      barcode = character(),
      CDR3d = character(),
      CDR3g = character(),
      sample_name = character(),
      group = character(),
      cell_type = character()
    ))
  }
  pair_records <- annotation %>%
    filter(barcode %in% marker_barcodes, chain %in% c("TRD", "TRG")) %>%
    distinct(barcode, chain, cdr3) %>%
    pivot_wider(names_from = chain, values_from = cdr3, values_fn = list) %>%
    split(.$barcode) %>%
    lapply(function(record) {
      trd <- unlist(record$TRD)
      trg <- unlist(record$TRG)
      if (length(trd) == 0 || length(trg) == 0 || !marker_cdr3 %in% trd) {
        return(NULL)
      }
      tibble(
        barcode = unique(record$barcode),
        CDR3d = marker_cdr3,
        CDR3g = trg
      )
    }) %>%
    bind_rows()

  if (nrow(pair_records) == 0) {
    return(tibble(
      barcode = character(),
      CDR3d = character(),
      CDR3g = character(),
      sample_name = character(),
      group = character(),
      cell_type = character()
    ))
  }

  pair_records %>%
    left_join(
      marker_cells %>% distinct(barcode, sample_name, group, cell_type),
      by = "barcode"
    ) %>%
    distinct()
}

summarise_marker_cdr3g_pairs <- function(marker_pairs, marker_cdr3) {
  if (nrow(marker_pairs) == 0) {
    return(tibble(
      rank = integer(),
      CDR3g = character(),
      CDR3d = character(),
      n = integer(),
      percent = numeric()
    ))
  }
  marker_pairs %>%
    dplyr::count(CDR3g, name = "n") %>%
    arrange(desc(n), CDR3g) %>%
    mutate(
      rank = row_number(),
      percent = n / sum(n) * 100,
      CDR3d = marker_cdr3
    ) %>%
    select(rank, CDR3g, CDR3d, n, percent)
}

spread_sankey_label_positions <- function(label_y, total, min_gap = NULL, pad = NULL) {
  if (length(label_y) <= 1 || total <= 0) {
    return(label_y)
  }
  if (is.null(min_gap)) {
    min_gap <- max(total * 0.045, 0.8)
  }
  if (is.null(pad)) {
    pad <- min(total * 0.04, 1)
  }
  lower <- pad
  upper <- total - pad
  if ((length(label_y) - 1) * min_gap > (upper - lower)) {
    return(seq(from = upper, to = lower, length.out = length(label_y))[rank(-label_y, ties.method = "first")])
  }

  ord <- order(label_y)
  y <- label_y[ord]
  for (i in seq_along(y)[-1]) {
    y[i] <- max(y[i], y[i - 1] + min_gap)
  }
  overflow <- max(y) - upper
  if (overflow > 0) {
    y <- y - overflow
  }
  for (i in rev(seq_along(y)[-length(y)])) {
    y[i] <- min(y[i], y[i + 1] - min_gap)
  }
  underflow <- lower - min(y)
  if (underflow > 0) {
    y <- y + underflow
  }
  out <- numeric(length(label_y))
  out[ord] <- y
  out
}

plot_marker_cdr3g_sankey <- function(marker_cdr3g_summary, top_n = 10) {
  if (nrow(marker_cdr3g_summary) == 0) {
    return(NULL)
  }
  target_cdr3d <- unique(marker_cdr3g_summary$CDR3d)
  top_cdr3g <- marker_cdr3g_summary %>%
    slice_head(n = top_n) %>%
    pull(CDR3g)

  plot_counts <- marker_cdr3g_summary %>%
    mutate(CDR3g_plot = if_else(CDR3g %in% top_cdr3g, CDR3g, "Other CDR3g")) %>%
    group_by(CDR3g_plot) %>%
    summarise(n = sum(n), .groups = "drop") %>%
    mutate(
      CDR3g_plot = factor(CDR3g_plot, levels = c(top_cdr3g, "Other CDR3g")),
      percent = n / sum(n) * 100
    ) %>%
    arrange(CDR3g_plot) %>%
    filter(n > 0) %>%
    mutate(
      left_ymin = lag(cumsum(n), default = 0),
      left_ymax = cumsum(n),
      right_ymin = left_ymin,
      right_ymax = left_ymax,
      label_y = (left_ymin + left_ymax) / 2,
      label = sprintf("%s\nn=%d (%.1f%%)", as.character(CDR3g_plot), n, percent)
    )

  label_data <- plot_counts %>%
    mutate(label_plot_y = spread_sankey_label_positions(label_y, sum(plot_counts$n)))

  smoothstep <- function(x) {
    3 * x^2 - 2 * x^3
  }
  make_ribbon <- function(row, steps = 80) {
    x <- seq(0, 1, length.out = steps)
    s <- smoothstep(x)
    upper <- tibble(x = x, y = (1 - s) * row$left_ymax + s * row$right_ymax)
    lower <- tibble(x = rev(x), y = rev((1 - s) * row$left_ymin + s * row$right_ymin))
    bind_rows(upper, lower) %>%
      mutate(CDR3g_plot = row$CDR3g_plot)
  }

  ribbon_data <- bind_rows(lapply(seq_len(nrow(plot_counts)), function(i) make_ribbon(plot_counts[i, ])))
  node_width <- 0.045
  left_nodes <- plot_counts %>%
    transmute(
      xmin = -node_width,
      xmax = node_width,
      ymin = left_ymin,
      ymax = left_ymax,
      CDR3g_plot
    )
  right_node <- tibble(
    xmin = 1 - node_width,
    xmax = 1 + node_width,
    ymin = 0,
    ymax = sum(plot_counts$n),
    label_y = sum(plot_counts$n) / 2,
    label = sprintf("%s\nn=%d", target_cdr3d, sum(plot_counts$n))
  )
  fill_values <- c(
    "#4E79A7", "#F28E2B", "#59A14F", "#E15759", "#76B7B2",
    "#EDC948", "#B07AA1", "#FF9DA7", "#9C755F", "#BAB0AC",
    "grey80"
  )
  names(fill_values) <- levels(plot_counts$CDR3g_plot)

  ggplot() +
    geom_polygon(
      data = ribbon_data,
      aes(x = x, y = y, group = CDR3g_plot, fill = CDR3g_plot),
      alpha = 0.72,
      color = NA
    ) +
    geom_rect(
      data = left_nodes,
      aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = CDR3g_plot),
      color = "white",
      linewidth = 0.35
    ) +
    geom_rect(
      data = right_node,
      aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
      fill = "#3B3B3B",
      color = "white",
      linewidth = 0.35
    ) +
    geom_segment(
      data = label_data,
      aes(x = -0.05, xend = -0.32, y = label_y, yend = label_plot_y),
      linewidth = 0.25,
      color = "grey55"
    ) +
    geom_text(
      data = label_data,
      aes(x = -0.34, y = label_plot_y, label = label),
      hjust = 1,
      size = 2.35,
      lineheight = 0.9
    ) +
    geom_text(
      data = right_node,
      aes(x = 1.08, y = label_y, label = label),
      hjust = 0,
      size = 3.4,
      lineheight = 0.9
    ) +
    annotate("text", x = 0, y = sum(plot_counts$n) * 1.04, label = "CDR3g", fontface = "bold", size = 4) +
    annotate("text", x = 1, y = sum(plot_counts$n) * 1.04, label = "CDR3d", fontface = "bold", size = 4) +
    scale_fill_manual(values = fill_values, drop = FALSE) +
    scale_x_continuous(limits = c(-1.05, 1.78), expand = expansion(mult = c(0.01, 0.01))) +
    scale_y_continuous(expand = expansion(mult = c(0.01, 0.08))) +
    labs(
      title = "Paired CDR3g sequences to dominant MSH2 CDR3d",
      subtitle = sprintf("Top %d CDR3g sequences paired with %s in MSH2 marker cells", top_n, target_cdr3d)
    ) +
    theme_void(base_size = 11) +
    theme(
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA),
      legend.position = "none",
      plot.title = element_text(face = "bold", hjust = 0.5, margin = margin(b = 4)),
      plot.subtitle = element_text(hjust = 0.5, color = "grey30", margin = margin(b = 10)),
      plot.margin = margin(12, 18, 12, 18)
    )
}

summarise_top_cdr3g_cdr3d_pairs <- function(annotation, marker_cdr3g_summary, top_n = 10) {
  if (nrow(marker_cdr3g_summary) == 0) {
    return(tibble(
      rank = integer(),
      CDR3g = character(),
      CDR3d = character(),
      n = integer(),
      percent = numeric()
    ))
  }
  top_cdr3g <- marker_cdr3g_summary %>%
    arrange(rank) %>%
    slice_head(n = 1) %>%
    pull(CDR3g)

  pair_records <- annotation %>%
    filter(chain %in% c("TRD", "TRG"), !is.na(cdr3), cdr3 != "") %>%
    distinct(barcode, chain, cdr3) %>%
    pivot_wider(names_from = chain, values_from = cdr3, values_fn = list) %>%
    split(.$barcode) %>%
    lapply(function(record) {
      trd <- unlist(record$TRD)
      trg <- unlist(record$TRG)
      if (length(trd) == 0 || length(trg) == 0 || !top_cdr3g %in% trg) {
        return(NULL)
      }
      tibble(
        barcode = unique(record$barcode),
        CDR3g = top_cdr3g,
        CDR3d = trd
      )
    }) %>%
    bind_rows()

  if (nrow(pair_records) == 0) {
    return(tibble(
      rank = integer(),
      CDR3g = character(),
      CDR3d = character(),
      n = integer(),
      percent = numeric()
    ))
  }

  pair_records %>%
    distinct(barcode, CDR3g, CDR3d) %>%
    dplyr::count(CDR3g, CDR3d, name = "n") %>%
    arrange(desc(n), CDR3d) %>%
    mutate(
      rank = row_number(),
      percent = n / sum(n) * 100
    ) %>%
    slice_head(n = top_n)
}

plot_top_cdr3g_to_cdr3d_sankey <- function(top_cdr3g_cdr3d_summary, marker_cdr3, top_n = 10) {
  if (nrow(top_cdr3g_cdr3d_summary) == 0) {
    return(NULL)
  }
  source_cdr3g <- unique(top_cdr3g_cdr3d_summary$CDR3g)
  plot_counts <- top_cdr3g_cdr3d_summary %>%
    mutate(
      CDR3d_plot = if_else(rank <= top_n, CDR3d, "Other CDR3d"),
      CDR3d_plot = factor(CDR3d_plot, levels = unique(CDR3d_plot))
    ) %>%
    group_by(CDR3d_plot) %>%
    summarise(n = sum(n), .groups = "drop") %>%
    mutate(
      percent = n / sum(n) * 100,
      left_ymin = lag(cumsum(n), default = 0),
      left_ymax = cumsum(n),
      right_ymin = left_ymin,
      right_ymax = left_ymax,
      label_y = (right_ymin + right_ymax) / 2,
      label = sprintf("%s\nn=%d (%.1f%%)", as.character(CDR3d_plot), n, percent),
      is_marker_cdr3 = as.character(CDR3d_plot) == marker_cdr3
    )

  label_data <- plot_counts %>%
    mutate(label_plot_y = spread_sankey_label_positions(label_y, sum(plot_counts$n)))
  smoothstep <- function(x) {
    3 * x^2 - 2 * x^3
  }
  make_ribbon <- function(row, steps = 80) {
    x <- seq(0, 1, length.out = steps)
    s <- smoothstep(x)
    upper <- tibble(x = x, y = (1 - s) * row$left_ymax + s * row$right_ymax)
    lower <- tibble(x = rev(x), y = rev((1 - s) * row$left_ymin + s * row$right_ymin))
    bind_rows(upper, lower) %>%
      mutate(CDR3d_plot = row$CDR3d_plot)
  }
  ribbon_data <- bind_rows(lapply(seq_len(nrow(plot_counts)), function(i) make_ribbon(plot_counts[i, ])))
  node_width <- 0.045
  left_node <- tibble(
    xmin = -node_width,
    xmax = node_width,
    ymin = 0,
    ymax = sum(plot_counts$n),
    label_y = sum(plot_counts$n) / 2,
    label = sprintf("%s\nn=%d", source_cdr3g, sum(plot_counts$n))
  )
  right_nodes <- plot_counts %>%
    transmute(
      xmin = 1 - node_width,
      xmax = 1 + node_width,
      ymin = right_ymin,
      ymax = right_ymax,
      CDR3d_plot,
      is_marker_cdr3
    )
  fill_values <- c(
    "#4E79A7", "#F28E2B", "#59A14F", "#E15759", "#76B7B2",
    "#EDC948", "#B07AA1", "#FF9DA7", "#9C755F", "#BAB0AC",
    "grey80"
  )
  names(fill_values) <- levels(plot_counts$CDR3d_plot)

  ggplot() +
    geom_polygon(
      data = ribbon_data,
      aes(x = x, y = y, group = CDR3d_plot, fill = CDR3d_plot),
      alpha = 0.72,
      color = NA
    ) +
    geom_rect(
      data = left_node,
      aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
      fill = "#3B3B3B",
      color = "white",
      linewidth = 0.35
    ) +
    geom_rect(
      data = right_nodes,
      aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = CDR3d_plot),
      color = "white",
      linewidth = 0.35
    ) +
    geom_segment(
      data = label_data,
      aes(x = 1.05, xend = 1.32, y = label_y, yend = label_plot_y),
      linewidth = 0.25,
      color = "grey55"
    ) +
    geom_text(
      data = label_data,
      aes(x = 1.34, y = label_plot_y, label = label, fontface = if_else(is_marker_cdr3, "bold", "plain")),
      hjust = 0,
      size = 2.35,
      lineheight = 0.9
    ) +
    geom_text(
      data = left_node,
      aes(x = -0.08, y = label_y, label = label),
      hjust = 1,
      size = 3.2,
      lineheight = 0.9
    ) +
    annotate("text", x = 0, y = sum(plot_counts$n) * 1.04, label = "CDR3g", fontface = "bold", size = 4) +
    annotate("text", x = 1, y = sum(plot_counts$n) * 1.04, label = "CDR3d", fontface = "bold", size = 4) +
    scale_fill_manual(values = fill_values, drop = FALSE) +
    scale_x_continuous(limits = c(-0.78, 2.05), expand = expansion(mult = c(0.01, 0.01))) +
    scale_y_continuous(expand = expansion(mult = c(0.01, 0.08))) +
    labs(
      title = "CDR3d sequences paired with the top MSH2-associated CDR3g",
      subtitle = sprintf("Top %d CDR3d partners of %s; %s is highlighted", top_n, source_cdr3g, marker_cdr3)
    ) +
    theme_void(base_size = 11) +
    theme(
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA),
      legend.position = "none",
      plot.title = element_text(face = "bold", hjust = 0.5, margin = margin(b = 4)),
      plot.subtitle = element_text(hjust = 0.5, color = "grey30", margin = margin(b = 10)),
      plot.margin = margin(12, 18, 12, 18)
    )
}

plot_combined_marker_pair_sankey <- function(marker_cdr3g_summary, top_cdr3g_cdr3d_summary, marker_cdr3, top_n = 10) {
  left_plot <- plot_marker_cdr3g_sankey(marker_cdr3g_summary, top_n = top_n)
  right_plot <- plot_top_cdr3g_to_cdr3d_sankey(top_cdr3g_cdr3d_summary, marker_cdr3 = marker_cdr3, top_n = top_n)
  if (is.null(left_plot) || is.null(right_plot)) {
    return(NULL)
  }

  top_cdr3g <- marker_cdr3g_summary %>%
    arrange(rank) %>%
    slice_head(n = 1) %>%
    pull(CDR3g)
  top_pair <- marker_cdr3g_summary %>%
    filter(CDR3g == top_cdr3g) %>%
    slice_head(n = 1)

  left_plot <- left_plot +
    labs(title = NULL, subtitle = NULL) +
    theme(
      plot.title = element_blank(),
      plot.subtitle = element_blank(),
      plot.margin = margin(8, 8, 8, 8)
    )
  right_plot <- right_plot +
    labs(title = NULL, subtitle = NULL) +
    theme(
      plot.title = element_blank(),
      plot.subtitle = element_blank(),
      plot.margin = margin(8, 8, 8, 8)
    )

  subtitle <- sprintf(
    "Left: CDR3g partners of %s. Right: CDR3d partners of the dominant CDR3g %s. Dominant pair: n=%d (%.1f%%).",
    marker_cdr3,
    top_cdr3g,
    top_pair$n,
    top_pair$percent
  )

  left_plot + right_plot +
    plot_layout(widths = c(1, 1), guides = "collect") +
    plot_annotation(
      title = "MSH2-associated gamma-delta CDR3 pairing",
      subtitle = subtitle
    ) &
    theme(
      plot.background = element_rect(fill = "white", color = NA),
      plot.title = element_text(face = "bold", hjust = 0.5, size = 15, margin = margin(b = 4)),
      plot.subtitle = element_text(hjust = 0.5, color = "grey30", size = 10, margin = margin(b = 8))
    )
}

is_vd1_cell <- function(cell_type) {
  str_detect(as.character(cell_type), "Vd1")
}

is_nonmarker_control_cell <- function(metadata, marker_barcodes) {
  is_vd1_cell(metadata$cell_type) &
    !(metadata$barcode %in% marker_barcodes$barcode) &
    tolower(as.character(metadata$group)) %in% tolower(nonmarker_control_groups)
}

add_msh2_marker_metadata <- function(seurat_obj, marker_barcodes) {
  # Non-Marker is intentionally a matched control: Vd1 cells from Blank/MSH2
  # samples that do not carry the dominant MSH2 TRD CDR3. Other cells are left
  # unlabelled so they cannot enter marker-vs-control summaries by accident.
  nonmarker_control <- is_nonmarker_control_cell(seurat_obj@meta.data, marker_barcodes)
  seurat_obj$MSH2_marker_cdr3 <- case_when(
    seurat_obj$barcode %in% marker_barcodes$barcode ~ "Marker",
    nonmarker_control ~ "Non-Marker",
    TRUE ~ NA_character_
  )
  seurat_obj$MSH2_marker_cdr3 <- factor(seurat_obj$MSH2_marker_cdr3, levels = c("Non-Marker", "Marker"))
  seurat_obj
}

subset_vd1_and_label_marker_groups <- function(seurat_obj, marker_barcodes) {
  vd1_obj <- subset(seurat_obj, is_vd1_cell(cell_type))
  nonmarker_control <- is_nonmarker_control_cell(vd1_obj@meta.data, marker_barcodes)
  vd1_obj$MSH2_marker_group <- case_when(
    vd1_obj$barcode %in% marker_barcodes$barcode ~ "Marker",
    vd1_obj$cell_type == "Naive Vd1" ~ "Naive",
    vd1_obj$cell_type == "Effector Vd1" ~ "Effector",
    nonmarker_control ~ "Non-Marker",
    TRUE ~ NA_character_
  )
  vd1_obj$MSH2_marker_group <- factor(vd1_obj$MSH2_marker_group, levels = marker_levels)
  vd1_obj
}

plot_marker_umap <- function(seurat_obj) {
  marker_cells <- rownames(seurat_obj@meta.data)[!is.na(seurat_obj$MSH2_marker_cdr3)]
  if (length(marker_cells) == 0) {
    return(NULL)
  }
  seurat_obj <- subset(seurat_obj, cells = marker_cells)
  reduction_name <- if ("umap.unintegrated" %in% names(seurat_obj@reductions)) {
    "umap.unintegrated"
  } else {
    "umap"
  }
  umap_plot <- DimPlot2(seurat_obj,
    group.by = "MSH2_marker_cdr3",
    reduction = reduction_name,
    cols = c("Non-Marker" = "lightgray", "Marker" = "#D62728"),
    theme = NoAxes(),
    label = TRUE,
    box = TRUE,
    label.color = "black",
    repel = TRUE
  )
  add_fixed_umap_coordinates(umap_plot) +
    labs(title = "Dominant MSH2 CDR3 cells and Blank/MSH2 Vd1 non-marker control")
}

plot_marker_group_umap <- function(seurat_obj) {
  reduction_name <- if ("umap.unintegrated" %in% names(seurat_obj@reductions)) {
    "umap.unintegrated"
  } else {
    "umap"
  }
  umap_plot <- DimPlot2(seurat_obj,
    group.by = "MSH2_marker_group",
    reduction = reduction_name,
    cols = marker_colors,
    theme = NoAxes(),
    label = TRUE,
    box = TRUE,
    label.color = "black",
    repel = TRUE
  )
  add_fixed_umap_coordinates(umap_plot) +
    labs(title = "Dominant MSH2 CDR3-bearing Vd1 cells relative to Vd1 states")
}

make_marker_summary <- function(seurat_obj) {
  seurat_obj@meta.data %>%
    filter(!is.na(MSH2_marker_cdr3)) %>%
    dplyr::count(group, sample_name, cell_type, MSH2_marker_cdr3, name = "n") %>%
    dplyr::group_by(group, sample_name, cell_type) %>%
    dplyr::mutate(percent = n / sum(n) * 100) %>%
    dplyr::ungroup()
}

plot_marker_celltype_composition <- function(marker_summary) {
  plot_data <- marker_summary %>%
    filter(
      MSH2_marker_cdr3 == "Marker" |
        (
          MSH2_marker_cdr3 == "Non-Marker" &
            tolower(as.character(group)) %in% tolower(nonmarker_control_groups) &
            is_vd1_cell(cell_type)
        )
    ) %>%
    group_by(MSH2_marker_cdr3, cell_type) %>%
    summarise(n = sum(n), .groups = "drop") %>%
    group_by(MSH2_marker_cdr3) %>%
    mutate(
      percent = n / sum(n) * 100,
      MSH2_marker_cdr3 = factor(as.character(MSH2_marker_cdr3), levels = c("Non-Marker", "Marker")),
      cell_type = factor(as.character(cell_type), levels = cell_type_levels)
    ) %>%
    ungroup() %>%
    filter(!is.na(cell_type), n > 0)
  if (nrow(plot_data) == 0) {
    return(NULL)
  }

  ggplot(plot_data, aes(x = MSH2_marker_cdr3, y = percent, fill = cell_type)) +
    geom_col(width = 0.65, color = "white") +
    geom_text(aes(label = ifelse(percent >= 4, paste0(n, " (", sprintf("%.1f", percent), "%)"), "")),
      position = position_stack(vjust = 0.5),
      size = 3
    ) +
    scale_fill_manual(values = color_celltype, drop = FALSE) +
    theme_test() +
    theme(axis.title.x = element_blank()) +
    ylab("Percent within CDR3 group") +
    labs(fill = "Cell type", title = "Cell-type composition of dominant MSH2 CDR3 cells and Blank/MSH2 non-marker control")
}

add_effector_activation_score <- function(seurat_obj, features) {
  present_features <- intersect(features, rownames(seurat_obj))
  if (length(present_features) == 0) {
    message("[SKIP] No effector activation score genes found in object.")
    seurat_obj$effector_activation_score <- NA_real_
    return(seurat_obj)
  }

  seurat_obj <- AddModuleScore(
    seurat_obj,
    features = list(present_features),
    name = "effector_activation_score"
  )
  seurat_obj$effector_activation_score <- seurat_obj$effector_activation_score1
  seurat_obj
}

summarise_effector_activation_score <- function(seurat_obj) {
  seurat_obj@meta.data %>%
    filter(!is.na(MSH2_marker_group), !is.na(effector_activation_score)) %>%
    group_by(MSH2_marker_group) %>%
    summarise(
      n = n(),
      mean_score = mean(effector_activation_score, na.rm = TRUE),
      median_score = median(effector_activation_score, na.rm = TRUE),
      sd_score = sd(effector_activation_score, na.rm = TRUE),
      .groups = "drop"
    )
}

plot_effector_activation_score <- function(seurat_obj) {
  plot_data <- seurat_obj@meta.data %>%
    filter(!is.na(MSH2_marker_group), !is.na(effector_activation_score)) %>%
    mutate(MSH2_marker_group = factor(as.character(MSH2_marker_group), levels = marker_levels))
  if (nrow(plot_data) == 0) {
    return(NULL)
  }
  present_levels <- levels(droplevels(plot_data$MSH2_marker_group))
  comparisons <- list(
    c("Marker", "Non-Marker"),
    c("Marker", "Naive"),
    c("Marker", "Effector")
  )
  comparisons <- comparisons[vapply(comparisons, function(x) all(x %in% present_levels), logical(1))]

  score_plot <- ggplot(plot_data, aes(x = MSH2_marker_group, y = effector_activation_score, fill = MSH2_marker_group)) +
    geom_boxplot(width = 0.55, outlier.shape = NA, alpha = 0.85) +
    geom_jitter(width = 0.18, size = 0.25, alpha = 0.25) +
    scale_fill_manual(values = marker_colors, drop = FALSE) +
    theme_test() +
    theme(legend.position = "none") +
    xlab(NULL) +
    ylab("Effector activation module score") +
    labs(title = "Dominant CDR3-bearing Vd1 cells show effector activation score")
  if (length(comparisons) > 0) {
    score_plot <- score_plot +
      stat_compare_means(comparisons = comparisons, method = "wilcox.test", label = "p.format", hide.ns = FALSE)
  }
  score_plot
}

make_hallmark_heatmap <- function(scores, groups, pathways = NULL) {
  if (!is.null(pathways)) {
    pathways <- pathways[pathways %in% rownames(scores)]
    scores <- scores[pathways, , drop = FALSE]
  }
  stats <- CalcStats(scores, f = groups)
  SeuratExtend::Heatmap(stats, lab_fill = "zscore")
}

calculate_hallmark_delta <- function(scores, groups, ident_1, ident_2, comparison_name) {
  names(groups) <- colnames(scores)
  common_cells <- intersect(colnames(scores), names(groups))
  scores <- scores[, common_cells, drop = FALSE]
  groups <- as.character(groups[common_cells])
  ident_1_cells <- names(groups)[groups == ident_1]
  ident_2_cells <- names(groups)[groups == ident_2]
  if (length(ident_1_cells) == 0 || length(ident_2_cells) == 0) {
    message(
      "[SKIP] Missing cells for comparison: ", ident_1, " vs ", ident_2,
      " (", ident_1, "=", length(ident_1_cells), ", ",
      ident_2, "=", length(ident_2_cells), ")"
    )
    return(tibble(
      comparison = character(),
      pathway = character(),
      ident_1 = character(),
      ident_2 = character(),
      mean_ident_1 = numeric(),
      mean_ident_2 = numeric(),
      n_ident_1 = integer(),
      n_ident_2 = integer(),
      delta = numeric(),
      abs_delta = numeric(),
      higher_in = character()
    ))
  }

  tibble(
    comparison = comparison_name,
    pathway = rownames(scores),
    ident_1 = ident_1,
    ident_2 = ident_2,
    mean_ident_1 = rowMeans(scores[, ident_1_cells, drop = FALSE]),
    mean_ident_2 = rowMeans(scores[, ident_2_cells, drop = FALSE]),
    n_ident_1 = length(ident_1_cells),
    n_ident_2 = length(ident_2_cells)
  ) %>%
    mutate(
      delta = mean_ident_1 - mean_ident_2,
      abs_delta = abs(delta),
      higher_in = ifelse(delta >= 0, ident_1, ident_2)
    ) %>%
    arrange(desc(abs_delta), pathway)
}

make_delta_plot <- function(delta_table, top_n = 10) {
  if (nrow(delta_table) == 0) {
    return(NULL)
  }
  plot_data <- delta_table %>%
    slice_max(order_by = abs_delta, n = top_n, with_ties = FALSE) %>%
    mutate(
      pathway_label = str_replace(pathway, "^HALLMARK_", ""),
      pathway_label = factor(pathway_label, levels = pathway_label[order(delta)])
    )

  ggplot(plot_data, aes(x = delta, y = pathway_label, color = higher_in)) +
    geom_vline(xintercept = 0, color = "grey65", linewidth = 0.4) +
    geom_segment(aes(x = 0, xend = delta, yend = pathway_label), linewidth = 0.8) +
    geom_point(size = 2.6) +
    scale_color_manual(values = marker_colors, drop = FALSE) +
    theme_test() +
    theme(axis.title.y = element_blank()) +
    labs(
      x = paste0("Mean AUCell delta: ", unique(delta_table$ident_1), " - ", unique(delta_table$ident_2)),
      color = "Higher in",
      title = paste(unique(delta_table$ident_1), "vs", unique(delta_table$ident_2))
    )
}

save_delta_plot_if_available <- function(delta_table, filename, width, height) {
  plot <- make_delta_plot(delta_table, top_n = 10)
  if (is.null(plot)) {
    message("[SKIP] No delta plot data for: ", filename)
    return(invisible(NULL))
  }
  save_plot(plot, filename, width, height, overwrite = force_msh2_plot)
}

calculate_umap_centroid_distances <- function(seurat_obj) {
  reduction_name <- if ("umap.unintegrated" %in% names(seurat_obj@reductions)) {
    "umap.unintegrated"
  } else if ("umap" %in% names(seurat_obj@reductions)) {
    "umap"
  } else {
    NA_character_
  }
  if (is.na(reduction_name)) {
    message("[SKIP] No UMAP reduction available for MSH2 centroid-distance analysis.")
    return(tibble())
  }

  umap_data <- Embeddings(seurat_obj, reduction_name) %>%
    as.data.frame() %>%
    rownames_to_column("cell_id") %>%
    as_tibble()
  colnames(umap_data)[2:3] <- c("UMAP_1", "UMAP_2")
  cell_data <- seurat_obj@meta.data %>%
    rownames_to_column("cell_id") %>%
    as_tibble() %>%
    left_join(umap_data, by = "cell_id") %>%
    filter(!is.na(MSH2_marker_group), is.finite(UMAP_1), is.finite(UMAP_2))

  reference_centroids <- cell_data %>%
    filter(MSH2_marker_group %in% c("Naive", "Non-Marker", "Effector")) %>%
    group_by(reference_group = as.character(MSH2_marker_group)) %>%
    summarise(
      reference_centroid_1 = mean(UMAP_1, na.rm = TRUE),
      reference_centroid_2 = mean(UMAP_2, na.rm = TRUE),
      n_reference_cells = n(),
      .groups = "drop"
    )
  marker_cells <- cell_data %>%
    filter(MSH2_marker_group == "Marker") %>%
    select(cell_id, barcode, group, sample_name, cell_type, UMAP_1, UMAP_2)
  if (nrow(reference_centroids) == 0 || nrow(marker_cells) == 0) {
    return(tibble())
  }

  marker_cells %>%
    mutate(.join_key = 1L) %>%
    full_join(reference_centroids %>% mutate(.join_key = 1L), by = ".join_key") %>%
    select(-.join_key) %>%
    mutate(
      reduction = reduction_name,
      distance_to_reference_centroid = sqrt(
        (UMAP_1 - reference_centroid_1)^2 +
          (UMAP_2 - reference_centroid_2)^2
      )
    ) %>%
    arrange(cell_id, distance_to_reference_centroid)
}

summarise_umap_centroid_distances <- function(distance_table) {
  if (nrow(distance_table) == 0) {
    return(tibble())
  }
  distance_table %>%
    group_by(reference_group) %>%
    summarise(
      n_marker_cells = n_distinct(cell_id),
      median_distance = median(distance_to_reference_centroid, na.rm = TRUE),
      mean_distance = mean(distance_to_reference_centroid, na.rm = TRUE),
      sd_distance = sd(distance_to_reference_centroid, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(mean_distance)
}

plot_umap_centroid_distances <- function(distance_table) {
  if (nrow(distance_table) == 0) {
    return(NULL)
  }
  plot_data <- distance_table %>%
    mutate(reference_group = factor(reference_group, levels = marker_levels))
  ggplot(plot_data, aes(x = reference_group, y = distance_to_reference_centroid, fill = reference_group)) +
    geom_boxplot(width = 0.55, outlier.shape = NA, alpha = 0.85) +
    geom_jitter(width = 0.18, size = 0.35, alpha = 0.35) +
    scale_fill_manual(values = marker_colors, drop = FALSE) +
    theme_test() +
    theme(legend.position = "none") +
    xlab("Reference Vd1 state centroid") +
    ylab("Marker-cell UMAP distance to centroid") +
    labs(title = "Dominant CDR3 cells are compared to Vd1 state centroids")
}

run_marker_group_de <- function(seurat_obj, comparisons) {
  seurat_obj <- join_assay_layers_if_needed(seurat_obj)
  Idents(seurat_obj) <- "MSH2_marker_group"
  present_groups <- levels(droplevels(seurat_obj$MSH2_marker_group))
  group_counts <- table(seurat_obj$MSH2_marker_group)
  purrr::imap_dfr(comparisons, function(comparison, comparison_name) {
    ident_1 <- comparison[[1]]
    ident_2 <- comparison[[2]]
    if (!all(c(ident_1, ident_2) %in% present_groups)) {
      message("[SKIP] Missing cells for DE comparison: ", ident_1, " vs ", ident_2)
      return(tibble())
    }
    if (unname(group_counts[[ident_1]]) < 3 || unname(group_counts[[ident_2]]) < 3) {
      message("[SKIP] Too few cells for DE comparison: ", ident_1, " vs ", ident_2)
      return(tibble())
    }
    markers <- FindMarkers(
      seurat_obj,
      ident.1 = ident_1,
      ident.2 = ident_2,
      logfc.threshold = 0,
      min.pct = 0.05
    )
    markers %>%
      rownames_to_column("gene") %>%
      mutate(
        comparison = comparison_name,
        ident_1 = ident_1,
        ident_2 = ident_2,
        .before = 1
      )
  })
}

plot_marker_effector_gene_heatmap <- function(seurat_obj, features) {
  features <- intersect(features, rownames(seurat_obj))
  if (length(features) == 0) {
    return(NULL)
  }
  expression_data <- FetchData(seurat_obj, vars = c("MSH2_marker_group", features)) %>%
    as_tibble() %>%
    filter(!is.na(MSH2_marker_group)) %>%
    mutate(MSH2_marker_group = factor(as.character(MSH2_marker_group), levels = marker_levels))
  if (nrow(expression_data) == 0) {
    return(NULL)
  }
  heatmap_data <- expression_data %>%
    pivot_longer(cols = all_of(features), names_to = "gene", values_to = "expression") %>%
    group_by(MSH2_marker_group, gene) %>%
    dplyr::summarise(mean_expression = mean(expression, na.rm = TRUE), .groups = "drop") %>%
    group_by(gene) %>%
    mutate(scaled_expression = as.numeric(scale(mean_expression))) %>%
    ungroup() %>%
    left_join(vd1_marker_gene_classes %>% select(gene_class, gene), by = "gene") %>%
    mutate(
      scaled_expression = ifelse(is.na(scaled_expression), 0, scaled_expression),
      gene_class = ifelse(is.na(gene_class), "Other", as.character(gene_class)),
      gene_class = factor(gene_class, levels = c(levels(vd1_marker_gene_classes$gene_class), "Other")),
      gene = factor(gene, levels = rev(features))
    )

  ggplot(heatmap_data, aes(x = MSH2_marker_group, y = gene, fill = scaled_expression)) +
    geom_tile(color = "white", linewidth = 0.4) +
    facet_grid(gene_class ~ ., scales = "free_y", space = "free_y") +
    scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0, name = "Scaled\nmean") +
    theme_test() +
    theme(axis.title = element_blank(), strip.background = element_blank()) +
    labs(title = "Grouped effector-gene expression pattern across Vd1 marker groups")
}

plot_marker_gene_dotplot <- function(seurat_obj, features) {
  features <- features[features %in% rownames(seurat_obj)]
  if (length(features) == 0) {
    message("No configured Vd1 marker genes found for MSH2 dotplot.")
    return(NULL)
  }
  expression_data <- FetchData(seurat_obj, vars = c("MSH2_marker_group", features)) %>%
    as_tibble() %>%
    filter(!is.na(MSH2_marker_group)) %>%
    mutate(MSH2_marker_group = factor(as.character(MSH2_marker_group), levels = marker_levels)) %>%
    pivot_longer(cols = all_of(features), names_to = "gene", values_to = "expression")
  plot_data <- expression_data %>%
    group_by(MSH2_marker_group, gene) %>%
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
    left_join(vd1_marker_gene_classes %>% select(gene_class, gene), by = "gene") %>%
    mutate(
      gene_class = ifelse(is.na(gene_class), "Other", as.character(gene_class)),
      gene_class = factor(gene_class, levels = c(levels(vd1_marker_gene_classes$gene_class), "Other")),
      gene = factor(gene, levels = rev(vd1_marker_gene_classes$gene[vd1_marker_gene_classes$gene %in% features]))
    )
  ggplot(plot_data, aes(x = MSH2_marker_group, y = gene)) +
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
    labs(title = "Grouped Vd1 marker genes across MSH2 marker groups")
}

# 1. Read annotated cells and productive VDJ annotations.
all_seurat_celltype <- read_rds_checked(seurat_celltype_rds, "annotated Seurat RDS")
all_annotation_included <- read_rds_checked(included_annotation_rds, "included VDJ annotation RDS")
normalised <- normalise_msh2_metadata_levels(all_seurat_celltype, all_annotation_included)
all_seurat_celltype <- normalised$seurat
all_annotation_included <- normalised$annotation

# 2. Mark cells carrying the MSH2-expanded TRD CDR3 sequence.
msh2_marker_cells <- get_marker_barcodes(all_annotation_included, msh2_marker_cdr3)
write_csv_if_missing(msh2_marker_cells, marker_cell_table, "MSH2 marker CDR3 cell table", overwrite = force_msh2)
all_seurat_celltype <- add_msh2_marker_metadata(all_seurat_celltype, msh2_marker_cells)

marker_summary <- make_marker_summary(all_seurat_celltype)
write_csv_if_missing(marker_summary, marker_summary_table, "MSH2 marker group summary table", overwrite = force_msh2)

marker_cdr3g_pairs <- make_marker_cdr3g_pairs(all_annotation_included, msh2_marker_cells, msh2_marker_cdr3)
marker_cdr3g_summary <- summarise_marker_cdr3g_pairs(marker_cdr3g_pairs, msh2_marker_cdr3)
write_csv_if_missing(
  marker_cdr3g_summary,
  cdr3g_pair_summary_table,
  "MSH2 marker paired CDR3g summary table",
  overwrite = force_msh2
)
marker_cdr3g_sankey <- plot_marker_cdr3g_sankey(marker_cdr3g_summary, top_n = top_msh2_cdr3g_n)
if (!is.null(marker_cdr3g_sankey)) {
  save_plot(
    marker_cdr3g_sankey,
    "msh2_marker_cdr3g_to_caldttfpigdrgytdklif_top10_sankey",
    12,
    8,
    overwrite = force_msh2_plot
  )
} else {
  message("[SKIP] No paired CDR3g records available for MSH2 marker Sankey.")
}

top_cdr3g_cdr3d_summary <- summarise_top_cdr3g_cdr3d_pairs(
  all_annotation_included,
  marker_cdr3g_summary,
  top_n = top_msh2_cdr3g_n
)
write_csv_if_missing(
  top_cdr3g_cdr3d_summary,
  cdr3d_pair_summary_table,
  "MSH2 marker top CDR3g paired CDR3d summary table",
  overwrite = force_msh2
)
top_cdr3g_cdr3d_sankey <- plot_top_cdr3g_to_cdr3d_sankey(
  top_cdr3g_cdr3d_summary,
  marker_cdr3 = msh2_marker_cdr3,
  top_n = top_msh2_cdr3g_n
)
if (!is.null(top_cdr3g_cdr3d_sankey)) {
  save_plot(
    top_cdr3g_cdr3d_sankey,
    "msh2_marker_top_cdr3g_to_cdr3d_top10_sankey",
    12,
    8,
    overwrite = force_msh2_plot
  )
} else {
  message("[SKIP] No paired CDR3d records available for top MSH2-associated CDR3g Sankey.")
}

combined_pair_sankey <- plot_combined_marker_pair_sankey(
  marker_cdr3g_summary,
  top_cdr3g_cdr3d_summary,
  marker_cdr3 = msh2_marker_cdr3,
  top_n = top_msh2_cdr3g_n
)
if (!is.null(combined_pair_sankey)) {
  save_plot(
    combined_pair_sankey,
    "msh2_marker_cdr3g_cdr3d_pairing_combined_sankey",
    18,
    8,
    overwrite = force_msh2_plot
  )
} else {
  message("[SKIP] No complete CDR3g/CDR3d records available for combined MSH2 pairing Sankey.")
}

marker_umap <- plot_marker_umap(all_seurat_celltype)
if (!is.null(marker_umap)) {
  save_plot(marker_umap, "msh2_marker_cdr3_umap", 10, 8, overwrite = force_msh2_plot)
} else {
  message("[SKIP] No marker/non-marker cells available for MSH2 CDR3 UMAP.")
}

marker_celltype_composition <- plot_marker_celltype_composition(marker_summary)
if (!is.null(marker_celltype_composition)) {
  save_plot(marker_celltype_composition, "msh2_marker_celltype_composition", 7, 5, overwrite = force_msh2_plot)
}

# 3. Compare Vd1 marker, non-marker, naive, and effector groups by Hallmark activity.
all_seurat_celltype_vd1 <- subset_vd1_and_label_marker_groups(all_seurat_celltype, msh2_marker_cells)
all_seurat_celltype_vd1 <- join_assay_layers_if_needed(all_seurat_celltype_vd1)

marker_group_umap <- plot_marker_group_umap(all_seurat_celltype_vd1)
save_plot(marker_group_umap, "msh2_marker_group_vd1_umap", 10, 8, overwrite = force_msh2_plot)

all_seurat_celltype_vd1 <- add_effector_activation_score(all_seurat_celltype_vd1, vd1_effector_activation_genes)
effector_score_summary <- summarise_effector_activation_score(all_seurat_celltype_vd1)
write_csv_if_missing(
  effector_score_summary,
  effector_score_table,
  "MSH2 marker effector activation score summary",
  overwrite = force_msh2
)
effector_score_plot <- plot_effector_activation_score(all_seurat_celltype_vd1)
if (!is.null(effector_score_plot)) {
  save_plot(effector_score_plot, "msh2_marker_effector_activation_score", 7, 5, overwrite = force_msh2_plot)
}

centroid_distances <- calculate_umap_centroid_distances(all_seurat_celltype_vd1)
if (nrow(centroid_distances) > 0) {
  write_csv_if_missing(
    centroid_distances,
    centroid_distance_table,
    "MSH2 marker UMAP centroid-distance table",
    overwrite = force_msh2
  )
  centroid_distance_summary <- summarise_umap_centroid_distances(centroid_distances)
  write_csv_if_missing(
    centroid_distance_summary,
    centroid_distance_summary_table,
    "MSH2 marker UMAP centroid-distance summary",
    overwrite = force_msh2
  )
  centroid_distance_plot <- plot_umap_centroid_distances(centroid_distances)
  if (!is.null(centroid_distance_plot)) {
    save_plot(centroid_distance_plot, "msh2_marker_umap_centroid_distance", 7, 5, overwrite = force_msh2_plot)
  }
} else {
  message("[SKIP] No MSH2 marker UMAP centroid-distance output available.")
}

marker_group_de <- run_marker_group_de(
  all_seurat_celltype_vd1,
  comparisons = list(
    marker_vs_nonmarker = c("Marker", "Non-Marker"),
    marker_vs_naive = c("Marker", "Naive"),
    marker_vs_effector = c("Marker", "Effector")
  )
)
if (nrow(marker_group_de) > 0) {
  write_csv_if_missing(marker_group_de, de_marker_table, "MSH2 marker differential-expression table", overwrite = force_msh2)
} else {
  message("[SKIP] No MSH2 marker differential-expression comparisons available to save.")
}

all_seurat_celltype_vd1 <- GeneSetAnalysis(all_seurat_celltype_vd1, genesets = hall50$human)
hallmark_scores <- all_seurat_celltype_vd1@misc$AUCell$genesets

hallmark_stats <- CalcStats(hallmark_scores, f = all_seurat_celltype_vd1$MSH2_marker_group)
write_csv_if_missing(
  as.data.frame(hallmark_stats) %>% rownames_to_column("pathway"),
  hallmark_stats_table,
  "MSH2 marker Hallmark stats table",
  overwrite = force_msh2
)

hallmark_heatmap_marker <- make_hallmark_heatmap(
  hallmark_scores,
  all_seurat_celltype_vd1$MSH2_marker_group
)
save_heatmap(hallmark_heatmap_marker, "msh2_marker_hallmark_heatmap", 10, 8, overwrite = force_msh2_plot)

hallmark_heatmap_marker_sub <- make_hallmark_heatmap(
  hallmark_scores,
  all_seurat_celltype_vd1$MSH2_marker_group,
  pathways = vd1_hallmark_pathways
)
save_heatmap(hallmark_heatmap_marker_sub, "msh2_marker_hallmark_heatmap_selected", 10, 8, overwrite = force_msh2_plot)

delta_marker_vs_nonmarker <- calculate_hallmark_delta(
  hallmark_scores,
  all_seurat_celltype_vd1$MSH2_marker_group,
  "Marker",
  "Non-Marker",
  "marker_vs_nonmarker"
)
delta_marker_vs_effector <- calculate_hallmark_delta(
  hallmark_scores,
  all_seurat_celltype_vd1$MSH2_marker_group,
  "Marker",
  "Effector",
  "marker_vs_effector"
)
hallmark_delta <- bind_rows(delta_marker_vs_nonmarker, delta_marker_vs_effector)
if (nrow(hallmark_delta) > 0) {
  write_csv_if_missing(hallmark_delta, hallmark_delta_table, "MSH2 marker Hallmark delta table", overwrite = force_msh2)
} else {
  message("[SKIP] No MSH2 marker Hallmark delta comparisons available to save.")
}

save_delta_plot_if_available(
  delta_marker_vs_nonmarker,
  "msh2_marker_vs_nonmarker_hallmark_delta",
  10,
  6
)
save_delta_plot_if_available(
  delta_marker_vs_effector,
  "msh2_marker_vs_effector_hallmark_delta",
  10,
  6
)

# 4. Dotplot selected co-stimulatory, effector, and checkpoint genes in Vd1 marker groups.
marker_gene_dotplot <- plot_marker_gene_dotplot(all_seurat_celltype_vd1, vd1_marker_genes)
if (!is.null(marker_gene_dotplot)) {
  save_plot(marker_gene_dotplot, "msh2_marker_vd1_gene_dotplot", 10, 5, overwrite = force_msh2_plot)
}

effector_gene_heatmap <- plot_marker_effector_gene_heatmap(all_seurat_celltype_vd1, vd1_effector_activation_genes)
if (!is.null(effector_gene_heatmap)) {
  save_plot(effector_gene_heatmap, "msh2_marker_effector_gene_heatmap", 8, 6, overwrite = force_msh2_plot)
}
