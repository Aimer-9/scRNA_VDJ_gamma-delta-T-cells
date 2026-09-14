# R version 4.5.2 (2025-10-31)
rm(list = ls())
setwd("/path/to/project")
library(Seurat)
library(SeuratExtend)
library(clusterProfiler)
library(tidyverse)
library(gplots)
library(ggpubr)
library(patchwork)
library(RColorBrewer)
library(ggsci)
library(monocle3)
library(SeuratWrappers)
library(ComplexHeatmap)
library(circlize)
library(msigdbr)
library(grid)
options(
  max.print = 100,
  tibble.width = Inf,
  spe = "human"
)

rds_dir <- "rds"
figure_dir <- file.path("figures", "9_CDR3paired")
table_dir <- "table"

seurat_celltype_rds <- file.path(rds_dir, "all_seurat_celltype.rds")
paired_cdr3_rds <- file.path(rds_dir, "barcode_trgd_paired.rds")
pair_rank_metadata_rds <- file.path(rds_dir, "trdg_pair_rank_metadata.rds")
toprank_cells_rds <- file.path(rds_dir, "all_seurat_celltype_toprank_cells.rds")
graph_test_rds <- file.path(rds_dir, "pr_graph_test_res_toprank.rds")
toprank_vd2_pairs_csv <- file.path(table_dir, "toprank_vd2_trdg_pairs_for_monocle3.csv")
trdg_clone_dispersion_csv <- file.path(table_dir, "trdg_clone_dispersion_metrics.csv")
trdg_pair_rank_metadata_csv <- file.path(table_dir, "trdg_pair_rank_metadata.csv")

top_rank_n <- 10
toprank_target_cells <- 5000
# Keep FALSE for normal checkpointed runs. Set TRUE after changing ranking
# rules, top-pair selection, or figure definitions so metadata caches and plots
# are rebuilt from the current paired CDR3 table.
force_ranked_rds <- FALSE
# Redraw figures/heatmaps without forcing rank metadata, selected-cell
# metadata, or graph-test RDS caches to be rebuilt.
force_ranked_plot <- TRUE
excluded_shared_pair_groups <- c("AB3")
# Top clone labels are exact paired TRD+TRG clonotypes. A ranked clone must be
# observed in all main comparison groups. Blank and MSH2 may also contain that
# same pair, but pairs found only in Blank or MSH2 are treated as background.
required_top_clone_groups <- c("Naive", "ZOL", "PAN")
trajectory_root_cell_type <- "Effector Memory Vd2"
min_cells_per_gene_for_graph_test <- 10
graph_test_core_fraction <- 0.7

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
colors_rank_pair <- c(brewer.pal(12, "Set3"), brewer.pal(8, "Set2"))

vd2_hallmark_pathways <- c(
  "HALLMARK_P53_PATHWAY",
  "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
  "HALLMARK_HYPOXIA",
  "HALLMARK_INTERFERON_ALPHA_RESPONSE",
  "HALLMARK_INTERFERON_GAMMA_RESPONSE",
  "HALLMARK_WNT_BETA_CATENIN_SIGNALING",
  "HALLMARK_IL6_JAK_STAT3_SIGNALING",
  "HALLMARK_INFLAMMATORY_RESPONSE",
  "HALLMARK_MTORC1_SIGNALING",
  "HALLMARK_FATTY_ACID_METABOLISM",
  "HALLMARK_OXIDATIVE_PHOSPHORYLATION",
  "HALLMARK_MITOTIC_SPINDLE"
)
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
  "HALLMARK_OXIDATIVE_PHOSPHORYLATION",
  "HALLMARK_MITOTIC_SPINDLE"
)

dir.create(rds_dir, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, showWarnings = FALSE)

get_single_pair_records <- function(paired_cdr3) {
  # Keep only cells with exactly one productive TRD/TRG pair. Cells with more
  # than one pair can reflect multiple contigs, doublets, or ambiguous pairing;
  # excluding them keeps pair-rank and clone-rank interpretation single-cell
  # specific.
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
  # Rank exact TRD+TRG pairs by the number of cells carrying that same paired
  # receptor. A pair is eligible for top-rank display only if the exact same
  # TRD+TRG pair appears in Naive, ZOL, and PAN. Blank and MSH2 are optional
  # extra groups, but cannot make a pair eligible by themselves.
  pair_rank <- paired_single %>%
    group_by(TRD, TRG) %>%
    summarise(
      n = n_distinct(barcode),
      clone_groups = list(unique(as.character(group))),
      .groups = "drop"
    ) %>%
    filter(map_lgl(clone_groups, ~all(required_top_clone_groups %in% .x))) %>%
    select(-clone_groups) %>%
    arrange(desc(n), TRD, TRG) %>%
    mutate(rank = row_number())
  paired_single %>%
    left_join(pair_rank, by = c("TRD", "TRG"))
}

make_pair_rank_metadata <- function(paired_ranked) {
  # Convert the full paired CDR3 table into a small barcode-level metadata
  # table. This is the only state needed to reattach pair rank labels to a
  # Seurat object; saving it avoids repeated writes of full Seurat objects.
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

save_pair_rank_metadata <- function(pair_rank_metadata) {
  # Save both RDS and CSV: RDS is used by the script for exact types/factors,
  # while CSV gives an inspectable table for troubleshooting clone-rank joins.
  if (force_ranked_rds) {
    saveRDS(pair_rank_metadata, pair_rank_metadata_rds)
    readr::write_csv(pair_rank_metadata, trdg_pair_rank_metadata_csv)
    message("[FORCE] Saved paired clone rank metadata RDS: ", pair_rank_metadata_rds)
    message("[FORCE] Saved paired clone rank metadata CSV: ", trdg_pair_rank_metadata_csv)
  } else {
    save_rds_if_missing(
      pair_rank_metadata,
      pair_rank_metadata_rds,
      "paired clone rank metadata RDS"
    )
    write_csv_if_missing(
      pair_rank_metadata,
      trdg_pair_rank_metadata_csv,
      "paired clone rank metadata CSV"
    )
  }
}

read_or_make_pair_rank_metadata <- function(paired_ranked) {
  # Pair ranks are deterministic for a fixed paired CDR3 input and eligibility
  # rule. Reuse the cache unless force mode explicitly requests a rebuild.
  if (output_files_exist(pair_rank_metadata_rds) && !force_ranked_rds) {
    message("[SKIP] Existing paired clone rank metadata RDS: ", pair_rank_metadata_rds)
    return(read_rds_checked(pair_rank_metadata_rds, "paired clone rank metadata RDS"))
  }
  if (output_files_exist(pair_rank_metadata_rds) && force_ranked_rds) {
    message("[FORCE] Rebuilding paired clone rank metadata: ", pair_rank_metadata_rds)
  }
  pair_rank_metadata <- make_pair_rank_metadata(paired_ranked)
  save_pair_rank_metadata(pair_rank_metadata)
  pair_rank_metadata
}

add_pair_rank_metadata <- function(seurat_obj, pair_rank_metadata) {
  # Remove older pair-rank columns before joining. Several users rerun this
  # script with cached RDS files, and stale metadata would otherwise produce
  # duplicated columns such as TRD.x/TRD.y or rank.x/rank.y. Old chain-only rank
  # columns are also removed because top-rank clone now means exact TRD+TRG pair.
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

calculate_entropy <- function(values) {
  counts <- table(values, useNA = "no")
  if (length(counts) == 0 || sum(counts) == 0) {
    return(NA_real_)
  }
  proportions <- as.numeric(counts) / sum(counts)
  -sum(proportions * log(proportions))
}

calculate_evenness <- function(values) {
  counts <- table(values, useNA = "no")
  if (length(counts) <= 1 || sum(counts) == 0) {
    return(0)
  }
  calculate_entropy(values) / log(length(counts))
}

calculate_simpson_dispersion <- function(values) {
  counts <- table(values, useNA = "no")
  if (length(counts) == 0 || sum(counts) == 0) {
    return(NA_real_)
  }
  proportions <- as.numeric(counts) / sum(counts)
  1 - sum(proportions^2)
}

summarise_rank_value <- function(rank_values) {
  rank_values <- suppressWarnings(as.numeric(rank_values))
  if (all(is.na(rank_values))) {
    return(NA_integer_)
  }
  as.integer(min(rank_values, na.rm = TRUE))
}

summarise_top10_label <- function(labels) {
  labels <- as.character(labels)
  labels <- labels[!is.na(labels) & labels != "other_pair"]
  if (length(labels) == 0) {
    return("other_pair")
  }
  labels[1]
}

rescale_dispersion_metric <- function(values) {
  values <- as.numeric(values)
  finite_values <- values[is.finite(values)]
  if (length(finite_values) == 0) {
    return(rep(NA_real_, length(values)))
  }
  value_range <- range(finite_values)
  if (value_range[1] == value_range[2]) {
    return(ifelse(is.finite(values), 0, NA_real_))
  }
  (values - value_range[1]) / (value_range[2] - value_range[1])
}

safe_mean <- function(values) {
  values <- as.numeric(values)
  values <- values[is.finite(values)]
  if (length(values) == 0) {
    return(NA_real_)
  }
  mean(values)
}

safe_median <- function(values) {
  values <- as.numeric(values)
  values <- values[is.finite(values)]
  if (length(values) == 0) {
    return(NA_real_)
  }
  median(values)
}

safe_max <- function(values) {
  values <- as.numeric(values)
  values <- values[is.finite(values)]
  if (length(values) == 0) {
    return(NA_real_)
  }
  max(values)
}

safe_quantile <- function(values, probability) {
  values <- as.numeric(values)
  values <- values[is.finite(values)]
  if (length(values) == 0) {
    return(NA_real_)
  }
  as.numeric(quantile(values, probability, names = FALSE))
}

calculate_trdg_clone_dispersion <- function(seurat_obj) {
  # This table is a clone-level dictionary for measuring how dispersed each
  # exact TRD+TRG clonotype is. Dispersion is intentionally based on cell
  # position and functional pathway activity only:
  #   1. UMAP spread: mean distance of clone cells from their UMAP centroid.
  #   2. Hallmark spread: mean per-pathway AUCell standard deviation across
  #      cells in the clone.
  # Group/sample/cell-type columns are kept as descriptive annotations, but
  # they do not contribute to dispersion_score.
  reduction_name <- if ("umap.unintegrated" %in% names(seurat_obj@reductions)) {
    "umap.unintegrated"
  } else if ("umap" %in% names(seurat_obj@reductions)) {
    "umap"
  } else {
    NA_character_
  }

  clone_metadata <- seurat_obj@meta.data %>%
    rownames_to_column("cell_id") %>%
    filter(!is.na(TRD), TRD != "", !is.na(TRG), TRG != "") %>%
    mutate(
      TRDG_pair = paste(TRD, TRG, sep = "||"),
      rank_top10 = as.character(rank_top10)
    )

  if (nrow(clone_metadata) == 0) {
    return(tibble())
  }

  # Hallmark AUCell scores provide a functional-space companion to UMAP space.
  # A clone with similar UMAP positions but heterogeneous pathway activity will
  # still receive a higher functional dispersion value.
  seurat_obj <- GeneSetAnalysis(seurat_obj, genesets = hall50$human)
  hallmark_scores <- seurat_obj@misc$AUCell$genesets
  hallmark_scores <- hallmark_scores[, intersect(colnames(hallmark_scores), clone_metadata$cell_id), drop = FALSE]

  if (!is.na(reduction_name)) {
    umap_data <- Embeddings(seurat_obj, reduction_name) %>%
      as.data.frame() %>%
      rownames_to_column("cell_id") %>%
      as_tibble()
    colnames(umap_data)[2:3] <- c("UMAP_1", "UMAP_2")
    clone_metadata <- clone_metadata %>%
      left_join(umap_data %>% select(cell_id, UMAP_1, UMAP_2), by = "cell_id")
  } else {
    clone_metadata <- clone_metadata %>%
      mutate(UMAP_1 = NA_real_, UMAP_2 = NA_real_)
  }

  hallmark_dispersion <- clone_metadata %>%
    group_by(TRDG_pair) %>%
    summarise(
      hallmark_mean_pathway_sd = {
        clone_cells <- intersect(cell_id, colnames(hallmark_scores))
        if (length(clone_cells) <= 1) {
          0
        } else {
          pathway_sd <- apply(hallmark_scores[, clone_cells, drop = FALSE], 1, sd, na.rm = TRUE)
          safe_mean(pathway_sd)
        }
      },
      hallmark_median_pathway_sd = {
        clone_cells <- intersect(cell_id, colnames(hallmark_scores))
        if (length(clone_cells) <= 1) {
          0
        } else {
          pathway_sd <- apply(hallmark_scores[, clone_cells, drop = FALSE], 1, sd, na.rm = TRUE)
          safe_median(pathway_sd)
        }
      },
      hallmark_max_pathway_sd = {
        clone_cells <- intersect(cell_id, colnames(hallmark_scores))
        if (length(clone_cells) <= 1) {
          0
        } else {
          pathway_sd <- apply(hallmark_scores[, clone_cells, drop = FALSE], 1, sd, na.rm = TRUE)
          safe_max(pathway_sd)
        }
      },
      .groups = "drop"
    )

  clone_metadata %>%
    group_by(TRD, TRG, TRDG_pair) %>%
    mutate(
      clone_umap_centroid_1 = mean(UMAP_1, na.rm = TRUE),
      clone_umap_centroid_2 = mean(UMAP_2, na.rm = TRUE),
      clone_umap_distance = sqrt(
        (UMAP_1 - clone_umap_centroid_1)^2 +
          (UMAP_2 - clone_umap_centroid_2)^2
      )
    ) %>%
    summarise(
      pair_rank = summarise_rank_value(rank),
      rank_top10 = summarise_top10_label(rank_top10),
      n_cells = n_distinct(cell_id),
      n_groups = n_distinct(group),
      n_samples = n_distinct(sample_name),
      n_cell_types = n_distinct(cell_type),
      groups_present = paste(sort(unique(as.character(group))), collapse = ";"),
      samples_present = paste(sort(unique(as.character(sample_name))), collapse = ";"),
      cell_types_present = paste(sort(unique(as.character(cell_type))), collapse = ";"),
      group_entropy = calculate_entropy(group),
      group_evenness = calculate_evenness(group),
      group_simpson_dispersion = calculate_simpson_dispersion(group),
      sample_entropy = calculate_entropy(sample_name),
      sample_evenness = calculate_evenness(sample_name),
      sample_simpson_dispersion = calculate_simpson_dispersion(sample_name),
      cell_type_entropy = calculate_entropy(cell_type),
      cell_type_evenness = calculate_evenness(cell_type),
      cell_type_simpson_dispersion = calculate_simpson_dispersion(cell_type),
      umap_reduction = reduction_name,
      umap_centroid_1 = dplyr::first(clone_umap_centroid_1),
      umap_centroid_2 = dplyr::first(clone_umap_centroid_2),
      umap_mean_distance_to_centroid = safe_mean(clone_umap_distance),
      umap_median_distance_to_centroid = safe_median(clone_umap_distance),
      umap_max_distance_to_centroid = safe_max(clone_umap_distance),
      umap_radius_90 = safe_quantile(clone_umap_distance, 0.9),
      .groups = "drop"
    ) %>%
    left_join(hallmark_dispersion, by = "TRDG_pair") %>%
    mutate(
      pair_rank = ifelse(is.infinite(pair_rank), NA_integer_, pair_rank),
      is_ranked_pair = !is.na(pair_rank),
      umap_mean_distance_scaled = rescale_dispersion_metric(umap_mean_distance_to_centroid),
      hallmark_mean_pathway_sd_scaled = rescale_dispersion_metric(hallmark_mean_pathway_sd),
      dispersion_score = rowMeans(
        cbind(umap_mean_distance_scaled, hallmark_mean_pathway_sd_scaled),
        na.rm = TRUE
      ),
      dispersion_score = ifelse(is.nan(dispersion_score), NA_real_, dispersion_score)
    ) %>%
    arrange(desc(is_ranked_pair), pair_rank, desc(n_cells), TRD, TRG)
}

prepare_ranked_clone_dispersion_plot_data <- function(dispersion_table, top_n = top_rank_n) {
  # Restrict visualization to ranked paired clones. Unranked pairs can be very
  # numerous and are still available in the CSV dictionary, but plotting them
  # together would obscure the biologically prioritized top clones.
  if (nrow(dispersion_table) == 0) {
    return(tibble())
  }

  dispersion_table %>%
    filter(is_ranked_pair, !is.na(pair_rank), pair_rank <= top_n) %>%
    arrange(pair_rank) %>%
    mutate(
      clone_label = paste0("R", pair_rank),
      clone_label = factor(clone_label, levels = paste0("R", seq_len(top_n)))
    )
}

plot_clone_dispersion_score <- function(dispersion_table, top_n = top_rank_n) {
  # Overall dispersion_score averages two dimensions: scaled UMAP spread and
  # scaled Hallmark pathway-score spread. Higher values indicate a clone whose
  # cells are more spatially and functionally dispersed.
  plot_data <- prepare_ranked_clone_dispersion_plot_data(dispersion_table, top_n)
  if (nrow(plot_data) == 0) {
    return(ggplot() +
      theme_void() +
      labs(title = "No ranked paired clones available for dispersion plot"))
  }
  clone_label_colors <- setNames(
    colors_rank_pair[seq_len(top_n)],
    paste0("R", seq_len(top_n))
  )

  ggplot(plot_data, aes(x = clone_label, y = dispersion_score, fill = clone_label)) +
    geom_col(color = "white", width = 0.8) +
    geom_text(aes(label = paste0("n=", n_cells)), vjust = -0.35, size = 3) +
    scale_fill_manual(values = clone_label_colors, drop = FALSE, guide = "none") +
    scale_y_continuous(limits = c(0, 1), expand = expansion(mult = c(0, 0.12))) +
    theme_test() +
    xlab("Top paired TRD+TRG clone rank") +
    ylab("Dispersion Score") +
    labs(title = "UMAP and Hallmark dispersion of top paired TRD+TRG clones")
}

plot_clone_dispersion_components <- function(dispersion_table, top_n = top_rank_n) {
  # Component heatmap for interpreting whether clone dispersion comes from UMAP
  # position, Hallmark pathway activity, or both. Both values are scaled across
  # clones in the dictionary.
  plot_data <- prepare_ranked_clone_dispersion_plot_data(dispersion_table, top_n)
  if (nrow(plot_data) == 0) {
    return(ggplot() +
      theme_void() +
      labs(title = "No ranked paired clones available for dispersion components"))
  }

  component_data <- plot_data %>%
    select(
      clone_label,
      umap_mean_distance_scaled,
      hallmark_mean_pathway_sd_scaled
    ) %>%
    pivot_longer(
      cols = -clone_label,
      names_to = "component",
      values_to = "value"
    ) %>%
    mutate(
      component = recode(
        component,
        umap_mean_distance_scaled = "UMAP position",
        hallmark_mean_pathway_sd_scaled = "Hallmark score"
      ),
      component = factor(component, levels = c("UMAP position", "Hallmark score"))
    )

  ggplot(component_data, aes(x = clone_label, y = component, fill = value)) +
    geom_tile(color = "white", linewidth = 0.5) +
    geom_text(aes(label = ifelse(is.na(value), "NA", sprintf("%.2f", value))), size = 3) +
    scale_fill_gradient2(
      low = "#2166AC",
      mid = "#F7F7F7",
      high = "#B2182B",
      midpoint = 0.5,
      limits = c(0, 1),
      na.value = "grey85",
      name = "Dispersion"
    ) +
    theme_test() +
    theme(axis.title = element_blank()) +
    labs(title = "UMAP and Hallmark dispersion components of top paired TRD+TRG clones")
}

plot_top10_pair_proportion <- function(seurat_obj) {
  # Show how much of each condition group is occupied by the globally top
  # paired TRD+TRG clones. AB3 is excluded here because it is not part of the
  # intended shared-clone comparison for this plot.
  seurat_obj@meta.data %>%
    group_by(rank_top10, group) %>%
    tally() %>%
    group_by(group) %>%
    mutate(rank_percent = n / sum(n)) %>%
    filter(rank_top10 != "other_pair", !group %in% excluded_shared_pair_groups) %>%
    mutate(rank_top10 = factor(
      rank_top10,
      levels = paste0("rank_", seq_len(top_rank_n), "_pair")
    )) %>%
    ggplot(aes(x = group, y = rank_percent * 100, fill = rank_top10)) +
    geom_col(color = "white", width = 0.8) +
    coord_flip() +
    scale_fill_manual(values = colors_rank_pair) +
    theme_test() +
    ylab("Percentage of Cells (%)") +
    xlab("") +
    labs(fill = "Top 10 TRDG clone")
}

plot_top10_pair_celltype_proportion <- function(seurat_obj) {
  # For each top paired TRD+TRG clone, show which annotated cell types contain
  # its cells. This answers whether an expanded clone is concentrated in one
  # state or spread across multiple Vd1/Vd2 states. The first bar is a reference
  # distribution from all paired TRD+TRG Vd2 clone cells, so ranked clones can
  # be compared against the overall Vd2-clone background.
  plot_metadata <- seurat_obj@meta.data %>%
    mutate(
      rank_top10 = factor(
        as.character(rank_top10),
        levels = paste0("rank_", seq_len(top_rank_n), "_pair")
      ),
      cell_type = factor(cell_type, levels = cell_type_levels)
    )

  ranked_clone_data <- plot_metadata %>%
    filter(!is.na(rank_top10), !is.na(cell_type)) %>%
    group_by(rank_top10, cell_type) %>%
    summarise(n = dplyr::n(), .groups = "drop") %>%
    group_by(rank_top10) %>%
    mutate(celltype_percent = n / sum(n) * 100) %>%
    ungroup() %>%
    transmute(
      clone_group = as.character(rank_top10),
      cell_type,
      n,
      celltype_percent
    )

  all_vd2_clone_data <- plot_metadata %>%
    filter(
      !is.na(cell_type),
      str_detect(as.character(cell_type), "Vd2"),
      !is.na(TRD), TRD != "",
      !is.na(TRG), TRG != ""
    ) %>%
    group_by(cell_type) %>%
    summarise(n = dplyr::n(), .groups = "drop") %>%
    mutate(
      clone_group = "All Vd2 clones",
      celltype_percent = n / sum(n) * 100
    ) %>%
    select(clone_group, cell_type, n, celltype_percent)

  bind_rows(all_vd2_clone_data, ranked_clone_data) %>%
    mutate(
      clone_group = factor(
        clone_group,
        levels = c("All Vd2 clones", paste0("rank_", seq_len(top_rank_n), "_pair"))
      )
    ) %>%
    ungroup() %>%
    ggplot(aes(x = clone_group, y = celltype_percent, fill = cell_type)) +
    geom_col(color = "white", width = 0.8) +
    scale_fill_manual(values = color_celltype, drop = FALSE) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.03))) +
    theme_test() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    xlab("Paired TRD+TRG clone group") +
    ylab("Percentage of Clone Cells (%)") +
    labs(fill = "Cell Type")
}

plot_top_pair_clone_umap <- function(seurat_obj) {
  # Color cells by paired TRD+TRG clone rank, not by cell type. This is the
  # direct UMAP companion to plot_rank_top10_proportion().
  pair_levels <- c(paste0("rank_", seq_len(top_rank_n), "_pair"), "other_pair")
  pair_colors <- setNames(
    c(colors_rank_pair[seq_len(top_rank_n)], "grey85"),
    pair_levels
  )

  umap_plot <- DimPlot2(seurat_obj,
    group.by = "rank_top10",
    reduction = "umap.unintegrated",
    theme = NoAxes(),
    label = FALSE,
    cols = pair_colors
  )
  add_fixed_umap_coordinates(umap_plot) +
    labs(
      title = paste0("Top ", top_rank_n, " paired TRD+TRG clones"),
      color = "Paired clone rank"
    )
}

choose_top_vd2_pairs <- function(paired_ranked, target_cells = toprank_target_cells) {
  # Select high-frequency Vd2 paired clones until the cumulative selected cell
  # count reaches the target size. This keeps the Monocle3 analysis focused on
  # recurrent Vd2 clonotypes while avoiding the full 100k-cell object.
  top_pairs <- paired_ranked %>%
    filter(
      !is.na(rank),
      !group %in% excluded_shared_pair_groups,
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
    arrange(desc(n_cells), rank, TRD, TRG) %>%
    mutate(
      top_pair_rank = row_number(),
      cumulative_cells = cumsum(n_cells)
    )

  if (nrow(top_pairs) == 0) {
    stop("No single-pair Vd2 TRD/TRG records available for top-pair Monocle3 subset.")
  }

  selected_n <- which(top_pairs$cumulative_cells >= target_cells)[1]
  if (is.na(selected_n)) {
    selected_n <- nrow(top_pairs)
  }
  top_pairs %>% slice_head(n = selected_n)
}

subset_toprank_vd2 <- function(seurat_obj, top_pairs) {
  # Subset by exact paired sequence ID. This avoids accidentally selecting cells
  # that share only TRD or only TRG with a selected pair.
  selected_pair_ids <- paste(top_pairs$TRD, top_pairs$TRG, sep = "||")
  cell_metadata <- seurat_obj@meta.data %>%
    rownames_to_column("cell_id") %>%
    mutate(pair_id = paste(TRD, TRG, sep = "||"))
  selected_cells <- cell_metadata %>%
    filter(
      !group %in% excluded_shared_pair_groups,
      str_detect(as.character(cell_type), "Vd2"),
      pair_id %in% selected_pair_ids
    ) %>%
    pull(cell_id)
  subset(seurat_obj, cells = selected_cells)
}

make_toprank_cell_metadata <- function(seurat_obj, top_pairs) {
  # Cache selected cells as metadata only. The downstream trajectory subset is
  # rebuilt from all_seurat_celltype.rds on each run, which avoids writing a
  # second large Seurat object for a selection that is just a cell list.
  selected_pair_ids <- paste(top_pairs$TRD, top_pairs$TRG, sep = "||")
  seurat_obj@meta.data %>%
    rownames_to_column("cell_id") %>%
    mutate(pair_id = paste(TRD, TRG, sep = "||")) %>%
    filter(
      !group %in% excluded_shared_pair_groups,
      str_detect(as.character(cell_type), "Vd2"),
      pair_id %in% selected_pair_ids
    )
}

subset_from_cell_metadata <- function(seurat_obj, cell_metadata) {
  # Reconstruct the top-rank Seurat subset from cached cell IDs. If no cells
  # overlap, the cache belongs to a different upstream object and should be
  # regenerated with force_ranked_rds <- TRUE.
  selected_cells <- intersect(cell_metadata$cell_id, colnames(seurat_obj))
  if (length(selected_cells) == 0) {
    stop("No cells from top-rank cell metadata are present in the Seurat object.", call. = FALSE)
  }
  subset(seurat_obj, cells = selected_cells)
}

plot_umap_cell_type <- function(seurat_obj) {
  umap_plot <- DimPlot2(seurat_obj,
    group.by = "cell_type",
    theme = NoAxes(),
    label = TRUE,
    box = TRUE,
    label.color = "black",
    repel = TRUE,
    cols = color_celltype
  )
  add_fixed_umap_coordinates(umap_plot)
}

make_hallmark_heatmap <- function(seurat_obj, pathways, group_by = "cell_type") {
  # Generic 2D Hallmark heatmap: rows are pathways and columns are categories
  # from group_by, usually cell_type. The clone-specific heatmap below uses a
  # custom ComplexHeatmap because it needs a third dimension as column splits.
  seurat_obj <- GeneSetAnalysis(seurat_obj, genesets = hall50$human)
  scores <- seurat_obj@misc$AUCell$genesets
  SeuratExtend::Heatmap(
    CalcStats(scores[pathways, ], f = seurat_obj@meta.data[[group_by]]),
    lab_fill = "zscore"
  )
}

make_top_pair_hallmark_heatmap <- function(seurat_obj, pathways) {
  # Build a three-dimensional clone heatmap:
  #   1. rows: Hallmark pathways,
  #   2. column_split: cell types,
  #   3. columns inside each split: top-ranked exact TRD+TRG pair labels.
  # Each column is the mean AUCell score for cells sharing both the paired clone
  # rank and cell type. The matrix is row z-scored so each pathway highlights
  # relative activity differences across pair/cell-type combinations.
  if (!"rank_top10" %in% colnames(seurat_obj@meta.data)) {
    stop("Missing paired clone metadata column: rank_top10", call. = FALSE)
  }

  # Only exact paired-clone rank labels are shown. Background cells are excluded
  # rather than plotted as an "other" column because that column would dominate
  # many heatmaps and obscure clone-specific patterns.
  clone_levels <- paste0("rank_", seq_len(top_rank_n), "_pair")
  clone_metadata <- seurat_obj$rank_top10
  selected_cells <- rownames(seurat_obj@meta.data)[
    !is.na(clone_metadata) & as.character(clone_metadata) %in% clone_levels
  ]

  if (length(selected_cells) == 0) {
    return(ggplot() +
      theme_void() +
      labs(title = "No top paired TRD+TRG clone cells available"))
  }

  clone_obj <- subset(seurat_obj, cells = selected_cells)
  clone_obj$rank_top10 <- factor(
    as.character(clone_obj$rank_top10),
    levels = clone_levels
  )
  # The interaction label defines one heatmap column per observed clone x
  # cell-type combination. Cell type and clone rank stay separate in
  # group_metadata so the visual can show cell types as column splits and clone
  # ranks as compact column labels within each cell-type block.
  clone_obj$clone_celltype_group <- interaction(
    clone_obj$rank_top10,
    clone_obj$cell_type,
    sep = " | ",
    drop = TRUE
  )

  # GeneSetAnalysis stores AUCell pathway scores in object@misc$AUCell$genesets.
  # The requested pathway vector is intersected with available rows to protect
  # against Hallmark set name changes or missing pathway scores.
  clone_obj <- GeneSetAnalysis(clone_obj, genesets = hall50$human)
  scores <- clone_obj@misc$AUCell$genesets
  present_pathways <- intersect(pathways, rownames(scores))
  if (length(present_pathways) == 0) {
    return(ggplot() +
      theme_void() +
      labs(title = "No Hallmark pathways available for top paired clones"))
  }

  group_metadata <- clone_obj@meta.data %>%
    rownames_to_column("cell_id") %>%
    transmute(
      cell_id,
      clone_group = factor(as.character(rank_top10), levels = clone_levels),
      cell_type = factor(cell_type, levels = cell_type_levels),
      clone_celltype_group = as.character(clone_celltype_group)
    )

  score_matrix <- as.matrix(scores[present_pathways, group_metadata$cell_id, drop = FALSE])
  group_stats <- group_metadata %>%
    distinct(clone_group, cell_type, clone_celltype_group) %>%
    arrange(cell_type, clone_group)

  # Average AUCell scores across cells for every clone x cell-type combination.
  # Missing combinations are not inserted; only observed clone/cell-type columns
  # appear in the heatmap.
  heatmap_matrix <- map_dfc(group_stats$clone_celltype_group, function(group_name) {
    group_cells <- group_metadata %>%
      filter(clone_celltype_group == group_name) %>%
      pull(cell_id)
    tibble(value = rowMeans(score_matrix[, group_cells, drop = FALSE]))
  }) %>%
    as.matrix()
  rownames(heatmap_matrix) <- present_pathways
  colnames(heatmap_matrix) <- str_replace(
    as.character(group_stats$clone_group),
    "^rank_(\\d+)_pair$",
    "R\\1"
  )

  # Row z-score scaling emphasizes pathway-specific relative differences across
  # clone/cell-type columns. Constant rows become NA after scaling, so they are
  # set to 0 to render as the neutral color.
  zscore_matrix <- t(scale(t(heatmap_matrix)))
  zscore_matrix[!is.finite(zscore_matrix)] <- 0

  Heatmap(
    zscore_matrix,
    name = "z-score",
    col = colorRamp2(c(-2, 0, 2), c("#053061", "#F7F7F7", "#67001F")),
    column_split = group_stats$cell_type,
    cluster_rows = TRUE,
    cluster_columns = FALSE,
    show_column_names = TRUE,
    column_names_rot = 0,
    row_names_gp = gpar(fontsize = 7),
    column_names_gp = gpar(fontsize = 6),
    column_title_gp = gpar(fontsize = 8, fontface = "bold"),
    column_gap = unit(2, "mm"),
    column_title = "%s",
    heatmap_legend_param = list(title = "z-score")
  )
}

prepare_monocle_object <- function(seurat_obj) {
  if ("umap.unintegrated" %in% names(seurat_obj@reductions)) {
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

run_pseudotime <- function(seurat_obj) {
  cds <- prepare_monocle_object(seurat_obj)
  cds <- cluster_cells(cds, reduction_method = "UMAP")
  cds <- learn_graph(cds, use_partition = TRUE)
  root_node <- get_earliest_principal_node(cds, trajectory_root_cell_type)
  order_cells(cds, root_pr_nodes = root_node)
}

plot_pseudotime_umap <- function(cds) {
  umap <- reducedDims(cds)$UMAP
  if (is.null(umap)) {
    stop("Missing UMAP reducedDims(cds)$UMAP for pseudotime plotting.", call. = FALSE)
  }
  plot_data <- tibble(
    cell_id = colnames(cds),
    UMAP_1 = umap[, 1],
    UMAP_2 = umap[, 2],
    pseudotime = pseudotime(cds)
  ) %>%
    filter(is.finite(pseudotime))

  umap_plot <- ggplot(plot_data, aes(x = UMAP_1, y = UMAP_2, color = pseudotime)) +
    geom_point(size = 0.8, alpha = 0.85) +
    scale_color_viridis_c(option = "plasma", name = "Pseudotime") +
    NoAxes()
  add_fixed_umap_coordinates(umap_plot)
}

plot_cell_cycle_umap <- function(cds) {
  umap <- reducedDims(cds)$UMAP
  if (is.null(umap)) {
    stop("Missing UMAP reducedDims(cds)$UMAP for cell-cycle plotting.", call. = FALSE)
  }
  if (!"Phase" %in% colnames(colData(cds))) {
    stop("Missing Phase column in Monocle3 colData for cell-cycle plotting.", call. = FALSE)
  }
  plot_data <- tibble(
    cell_id = colnames(cds),
    UMAP_1 = umap[, 1],
    UMAP_2 = umap[, 2],
    Phase = factor(as.character(colData(cds)$Phase), levels = names(color_cellcycle))
  )

  umap_plot <- ggplot(plot_data, aes(x = UMAP_1, y = UMAP_2, color = Phase)) +
    geom_point(size = 0.8, alpha = 0.85) +
    scale_color_manual(values = color_cellcycle, na.value = "grey80", drop = FALSE) +
    theme(legend.position = "right") +
    NoAxes()
  add_fixed_umap_coordinates(umap_plot)
}

run_graph_test <- function(cds) {
  sf_available <- requireNamespace("sf", quietly = TRUE) && tryCatch({
    library(sf, quietly = TRUE, character.only = TRUE)
    TRUE
  }, error = function(e) {
    message("[SKIP] Monocle3 graph_test unavailable because sf cannot load: ", conditionMessage(e))
    FALSE
  })
  if (!sf_available) {
    return(NULL)
  }

  cores <- max(1, round(parallel::detectCores() * graph_test_core_fraction))
  cells_per_gene <- rowSums(exprs(cds) > 0)
  genes_to_test <- names(cells_per_gene[cells_per_gene > min_cells_per_gene_for_graph_test])
  subset_cds <- cds[genes_to_test, ]
  tryCatch(
    graph_test(
      subset_cds,
      neighbor_graph = "principal_graph",
      cores = cores
    ),
    error = function(e) {
      message("[SKIP] Monocle3 graph_test failed: ", conditionMessage(e))
      NULL
    }
  )
}

make_pseudotime_heatmap_matrix <- function(cds, graph_test_result) {
  if (is.null(graph_test_result) || nrow(graph_test_result) == 0) {
    return(NULL)
  }
  sig_genes <- graph_test_result %>%
    filter(q_value == 0 & morans_test_statistic > 0.3) %>%
    arrange(desc(morans_test_statistic)) %>%
    rownames()
  if (length(sig_genes) == 0) {
    message("[SKIP] No significant graph_test genes for pseudotime heatmap.")
    return(NULL)
  }
  pt_matrix <- as.matrix(exprs(cds)[
    match(sig_genes, rownames(rowData(cds))),
    order(pseudotime(cds))
  ])
  pt_matrix <- t(apply(pt_matrix, 1, function(x) smooth.spline(x, df = 3)$y))
  t(apply(pt_matrix, 1, function(x) (x - mean(x)) / sd(x)))
}

make_pseudotime_heatmap <- function(pt_matrix) {
  Heatmap(
    pt_matrix,
    name = "z-score",
    col = colorRamp2(c(-2, 0, 2), c("#053061", "#F7F7F7", "#67001F")),
    show_row_names = FALSE,
    show_column_names = FALSE,
    row_names_gp = gpar(fontsize = 6),
    km = 5,
    row_title_rot = 0,
    cluster_rows = TRUE,
    cluster_row_slices = FALSE,
    cluster_columns = FALSE
  )
}

make_hallmark_term2gene <- function() {
  msigdbr(species = "Homo sapiens", category = "H") %>%
    dplyr::select(gs_name, gene_symbol)
}

run_hallmark_enrich <- function(gene_list, term2gene) {
  enricher(
    gene = gene_list,
    TERM2GENE = term2gene,
    pvalueCutoff = 0.1
  )
}

plot_hallmark_enrichment <- function(hallmark_enrichment, top_n = 10) {
  if (is.null(hallmark_enrichment) || nrow(as.data.frame(hallmark_enrichment)) == 0) {
    return(NULL)
  }
  parse_ratio <- function(ratio_str) {
    sapply(
      strsplit(as.character(ratio_str), "/"),
      function(x) as.numeric(x[1]) / as.numeric(x[2])
    )
  }
  plot_data <- as.data.frame(hallmark_enrichment) %>%
    mutate(
      Ratio = parse_ratio(GeneRatio),
      Description = str_replace(ID, "HALLMARK_", "")
    ) %>%
    slice_min(p.adjust, n = top_n)

  ggplot(plot_data, aes(x = Ratio, y = reorder(Description, Ratio), fill = p.adjust)) +
    geom_col(width = 0.75, color = "white", linewidth = 0.3) +
    scale_fill_gradient2(
      low = "#8B0000",
      mid = "#F2F2F2",
      high = "#003366",
      midpoint = 0.05
    ) +
    scale_x_continuous(labels = abs, expand = expansion(mult = c(0.3, 0.3))) +
    theme_test() +
    theme(axis.title.y = element_blank()) +
    labs(x = "Gene Ratio", title = "MSigDB Hallmark Enrichment") +
    geom_vline(xintercept = 0, color = "black")
}

# 1. Read annotated cells and paired productive TRD/TRG CDR3 records.
all_seurat_celltype <- read_rds_checked(seurat_celltype_rds, "annotated Seurat RDS")
# Normalize metadata levels early so joins, plots, and heatmap split order are
# consistent with the rest of the downstream analysis.
all_seurat_celltype$group <- factor(all_seurat_celltype$group, levels = group_levels)
all_seurat_celltype$sample_name <- factor(all_seurat_celltype$sample_name, levels = sample_name_levels)
all_seurat_celltype$cell_type <- factor(all_seurat_celltype$cell_type, levels = cell_type_levels)

barcode_trgd_paired <- read_rds_checked(paired_cdr3_rds, "paired TRD/TRG CDR3 RDS") %>%
  mutate(
    group = factor(group, levels = group_levels),
    sample_name = factor(sample_name, levels = sample_name_levels),
    cell_type = factor(cell_type, levels = cell_type_levels)
  )

# 2. Rank single-cell exact TRD+TRG pairs and add paired-clone rank metadata.
barcode_paired_single_all <- get_single_pair_records(barcode_trgd_paired)
barcode_paired_ranked <- rank_paired_clones(barcode_paired_single_all)
# Rank metadata is saved separately from the Seurat object. Joining these few
# columns in memory is much faster than repeatedly writing a full ranked Seurat
# RDS after every ranking-rule change.
pair_rank_metadata <- read_or_make_pair_rank_metadata(barcode_paired_ranked)
all_seurat_celltype <- add_pair_rank_metadata(all_seurat_celltype, pair_rank_metadata)

trdg_clone_dispersion <- calculate_trdg_clone_dispersion(all_seurat_celltype)
write_csv_if_missing(
  trdg_clone_dispersion,
  trdg_clone_dispersion_csv,
  "TRD+TRG clone dispersion dictionary",
  overwrite = force_ranked_rds
)

plot_trdg_clone_dispersion_score <- plot_clone_dispersion_score(trdg_clone_dispersion)
save_plot(
  plot_trdg_clone_dispersion_score,
  "trdg_clone_dispersion_score",
  8,
  5,
  overwrite = force_ranked_plot
)

plot_trdg_clone_dispersion_components <- plot_clone_dispersion_components(trdg_clone_dispersion)
save_plot(
  plot_trdg_clone_dispersion_components,
  "trdg_clone_dispersion_components",
  8,
  4,
  overwrite = force_ranked_plot
)

# Pair-level clone burden by condition group and UMAP distribution of exact
# paired TRD+TRG clone ranks.
plot_rank_top10_proportion <- plot_top10_pair_proportion(all_seurat_celltype)
save_plot(plot_rank_top10_proportion, "plot_rank_top10_proportion", 10, 2, overwrite = force_ranked_plot)

# Complement the condition-level clone-burden plot with a clone-intrinsic view:
# each bar is one ranked paired clone, and the stack shows which cell states
# make up that clone.
plot_rank_top10_celltype_proportion <- plot_top10_pair_celltype_proportion(all_seurat_celltype)
save_plot(
  plot_rank_top10_celltype_proportion,
  "plot_rank_top10_celltype_proportion",
  10,
  5,
  overwrite = force_ranked_plot
)

plot_toprank_pair_clone_umap <- plot_top_pair_clone_umap(all_seurat_celltype)
save_plot(plot_toprank_pair_clone_umap, "toprank_umap_TRDG_pair_clone", 8, 6, overwrite = force_ranked_plot)

# 3. Select top-frequency Vd2 TRD/TRG pairs to about 5k cells and score cell
# cycle. Only the selected cell metadata is cached; the Seurat subset is rebuilt
# in memory from the main annotated object so expression/reduction data are not
# written repeatedly.
if (output_files_exist(toprank_cells_rds) && !force_ranked_rds) {
  message("[SKIP] Existing top-pair Vd2 cell metadata RDS: ", toprank_cells_rds)
  toprank_cell_metadata <- read_rds_checked(toprank_cells_rds, "top-pair Vd2 cell metadata RDS")
} else {
  if (output_files_exist(toprank_cells_rds) && force_ranked_rds) {
    message("[FORCE] Rebuilding top-pair Vd2 cell metadata from current paired clone ranks: ", toprank_cells_rds)
  } else {
    message("[RUN] Building top-pair Vd2 cell metadata: ", toprank_cells_rds)
  }
  toprank_vd2_pairs <- choose_top_vd2_pairs(
    barcode_paired_ranked,
    target_cells = toprank_target_cells
  )
  write_csv_if_missing(
    toprank_vd2_pairs,
    toprank_vd2_pairs_csv,
    "top-pair Vd2 TRD/TRG table",
    overwrite = force_ranked_rds
  )
  message(
    "Selected ", nrow(toprank_vd2_pairs), " Vd2 TRD/TRG pairs covering ",
    max(toprank_vd2_pairs$cumulative_cells), " single-pair cells for Monocle3."
  )
  toprank_cell_metadata <- make_toprank_cell_metadata(all_seurat_celltype, toprank_vd2_pairs)
  save_rds_if_missing(
    toprank_cell_metadata,
    toprank_cells_rds,
    "top-pair Vd2 cell metadata RDS",
    overwrite = force_ranked_rds
  )
}

all_seurat_celltype_toprank <- subset_from_cell_metadata(all_seurat_celltype, toprank_cell_metadata)
message("Monocle3 Vd2 top-pair Seurat subset cells: ", ncol(all_seurat_celltype_toprank))

# The top-rank subset is rebuilt from cached cell metadata, so refresh clone
# rank metadata after subsetting. This guarantees downstream clone heatmaps use
# the current ranking rules.
all_seurat_celltype_toprank <- add_pair_rank_metadata(all_seurat_celltype_toprank, pair_rank_metadata)

# Cell-cycle scores are added after clone metadata refresh because Monocle3 and
# later UMAP plots use the same top-rank object.
all_seurat_celltype_toprank <- CellCycleScoring(all_seurat_celltype_toprank,
  s.features = cc.genes$s.genes,
  g2m.features = cc.genes$g2m.genes,
  set.ident = TRUE
)
plot_toprank_celltype <- plot_umap_cell_type(all_seurat_celltype_toprank)
save_plot(plot_toprank_celltype, "toprank_umap_celltype", 8, 6, overwrite = force_ranked_plot)

plot_toprank_subset_pair_clone_umap <- plot_top_pair_clone_umap(all_seurat_celltype_toprank)
save_plot(
  plot_toprank_subset_pair_clone_umap,
  "toprank_subset_umap_TRDG_pair_clone",
  8,
  6,
  overwrite = force_ranked_plot
)

# 4. Hallmark heatmaps for top-pair Vd2, top clone groups, all Vd2, and all Vd1 cells.
# This first heatmap keeps the older view: Hallmark pathways summarized by cell
# type within the selected top-pair Vd2 subset.
hallmark_heatmap_top_vd2_pairs <- make_hallmark_heatmap(
  all_seurat_celltype_toprank,
  vd2_hallmark_pathways
)
save_plot(hallmark_heatmap_top_vd2_pairs, "hallmark_heatmap_top_vd2_pairs", 8, 5, overwrite = force_ranked_plot)

# This heatmap is the paired-clone-focused view: cell type is the main
# ComplexHeatmap column split/title, and compact R1-R10 clone labels are the
# small columns inside each cell-type block.
hallmark_heatmap_toprank_pair_clones <- make_top_pair_hallmark_heatmap(
  all_seurat_celltype_toprank,
  vd2_hallmark_pathways
)
save_heatmap(
  hallmark_heatmap_toprank_pair_clones,
  "hallmark_heatmap_toprank_TRDG_pair_clones",
  10,
  6,
  overwrite = force_ranked_plot
)

all_seurat_celltype_vd2 <- subset(all_seurat_celltype, str_detect(cell_type, "Vd2"))
hallmark_heatmap_all_vd2 <- make_hallmark_heatmap(
  all_seurat_celltype_vd2,
  vd2_hallmark_pathways
)
save_plot(hallmark_heatmap_all_vd2, "hallmark_heatmap_all_vd2", 8, 5, overwrite = force_ranked_plot)

all_seurat_celltype_vd1 <- subset(all_seurat_celltype, str_detect(cell_type, "Vd1"))
hallmark_heatmap_all_vd1 <- make_hallmark_heatmap(
  all_seurat_celltype_vd1,
  vd1_hallmark_pathways
)
save_plot(hallmark_heatmap_all_vd1, "hallmark_heatmap_all_vd1", 8, 5, overwrite = force_ranked_plot)

# 5. Monocle3 trajectory and cell-cycle UMAPs for top-frequency Vd2 TRD/TRG-pair cells.
# Monocle3 uses the top-pair Vd2 subset rather than all cells to keep the
# trajectory focused on recurrent Vd2 clonotypes and to keep graph learning
# tractable.
cds_toprank <- run_pseudotime(all_seurat_celltype_toprank)
plot_umap_pseudotime <- plot_pseudotime_umap(cds_toprank)
save_plot(plot_umap_pseudotime, "plot_umap_pseudotime", 7, 5, overwrite = force_ranked_plot)

plot_umap_cellcycle <- plot_cell_cycle_umap(cds_toprank)
save_plot(plot_umap_cellcycle, "plot_umap_cellcycle", 10, 8, overwrite = force_ranked_plot)

# 8. Graph test, pseudotime heatmap, and Hallmark enrichment of gene modules.
# graph_test depends on sf and system PROJ libraries. run_graph_test() returns
# NULL instead of stopping when sf is unavailable, so upstream clone plots and
# heatmaps are still produced in imperfect R environments.
if (output_files_exist(graph_test_rds) && !force_ranked_rds) {
  message("[SKIP] Existing Monocle3 graph test RDS: ", graph_test_rds)
  pr_graph_test_res_toprank <- read_rds_checked(graph_test_rds, "Monocle3 graph test RDS")
} else {
  if (output_files_exist(graph_test_rds) && force_ranked_rds) {
    message("[FORCE] Rebuilding Monocle3 graph test RDS: ", graph_test_rds)
  }
  pr_graph_test_res_toprank <- run_graph_test(cds_toprank)
  if (!is.null(pr_graph_test_res_toprank)) {
    save_rds_if_missing(
      pr_graph_test_res_toprank,
      graph_test_rds,
      "Monocle3 graph test RDS",
      overwrite = force_ranked_rds
    )
  }
}

pt_matrix <- make_pseudotime_heatmap_matrix(cds_toprank, pr_graph_test_res_toprank)
if (is.null(pt_matrix)) {
  message("[SKIP] Pseudotime heatmap and Hallmark enrichment because graph_test output is unavailable.")
} else {
  # The pseudotime heatmap clusters graph-test genes, then each gene module is
  # interpreted with Hallmark enrichment. Currently only module 3 is plotted
  # because that was the biologically inspected module in the original analysis.
  pseudotime_heatmap <- make_pseudotime_heatmap(pt_matrix)
  save_heatmap(pseudotime_heatmap, "pseudotime_heatmap", 10, 8, overwrite = force_ranked_plot)

  drawn_pseudotime_heatmap <- draw(pseudotime_heatmap)
  row_dendrogram <- row_dend(drawn_pseudotime_heatmap)
  hallmark_term2gene <- make_hallmark_term2gene()
  hallmark_enrichment_list <- lapply(row_dendrogram, function(dendrogram) {
    dendrogram %>%
      as.hclust() %>%
      cutree(h = 1) %>%
      names() %>%
      run_hallmark_enrich(hallmark_term2gene)
  })
  plot_hallmark_module_3 <- plot_hallmark_enrichment(hallmark_enrichment_list[[3]])
  if (is.null(plot_hallmark_module_3)) {
    message("[SKIP] Hallmark enrichment plot for cluster 3 has no enriched terms.")
  } else {
    plot_hallmark_module_3 <- plot_hallmark_module_3 +
      ggtitle("Hallmark Enrichment for Cluster 3")
    save_plot(plot_hallmark_module_3, "hallmark_enrichment_cluster_3", 8, 5, overwrite = force_ranked_plot)
  }
}

# 6. Hallmark heatmap across top-ranked cell types.
# Final compact summary of Hallmark programs by cell type in the same top-pair
# Vd2 subset. This complements the clone x cell-type split heatmaps above.
all_seurat_celltype_toprank <- GeneSetAnalysis(all_seurat_celltype_toprank,
  genesets = hall50$human
)
toprank_stats <- CalcStats(
  all_seurat_celltype_toprank@misc$AUCell$genesets,
  f = all_seurat_celltype_toprank$cell_type
)
all_seurat_celltype_toprank_heatmap_hallmark <- SeuratExtend::Heatmap(
  toprank_stats[vd2_hallmark_pathways, ],
  lab_fill = "zscore"
)
save_plot(
  all_seurat_celltype_toprank_heatmap_hallmark,
  "all_seurat_celltype_toprank_heatmap_hallmark",
  10,
  10,
  overwrite = force_ranked_plot
)
