#!/usr/bin/env Rscript

# Visualize pySCENIC regulon activity from auc_mtx.loom.
#
# This script is intentionally separate from 11_pySCENIC.sh because the
# pySCENIC analysis runs in a pinned Python environment, while visualization
# uses the project's R/Seurat environment.
#
# Required inputs:
#   1. pyscenic_10k/auc_mtx.loom
#      Produced by the pySCENIC AUCell step. The loom must contain:
#        /col_attrs/CellID
#        /col_attrs/RegulonsAUC
#   2. rds/all_seurat_celltype.rds
#      Supplies metadata groups and the existing UMAP coordinates.
#
# Default command:
#   Rscript code/downstream/12_pySCENIC_visualization.R
#
# Explicit command:
#   Rscript code/downstream/12_pySCENIC_visualization.R \
#     --project-dir /path/to/project \
#     --auc-loom pyscenic_10k/auc_mtx.loom \
#     --seurat-rds rds/all_seurat_celltype.rds \
#     --group-by cell_type \
#     --reduction umap.unintegrated \
#     --n-top-regulons 30 \
#
# Analysis performed:
#   1. Validate and read the AUCell loom into a regulons-by-cells matrix.
#   2. Align loom CellID values with Seurat cell names. Cells absent from
#      either input are dropped; no positional matching is used.
#   3. Calculate mean regulon AUC for each metadata group.
#   4. Standardize each regulon across groups and rank regulons by the range
#      of their group z-scores. A large range identifies regulons that vary
#      strongly between the selected cell groups.
#   5. Draw:
#        - grouped z-score heatmap for the top regulons;
#        - dot plot where color is mean AUC and size is the percentage of
#          cells above that regulon's global median AUC;
#
# Outputs:
#   figures/12_pySCENIC_visualization/pyscenic_regulon_activity_heatmap_by_<group>.{png,pdf}
#   figures/12_pySCENIC_visualization/pyscenic_regulon_activity_dotplot_by_<group>.{png,pdf}
#   figures/12_pySCENIC_visualization/pyscenic_curated_tf_heatmap.{png,pdf}
#   figures/12_pySCENIC_visualization/pyscenic_vd2_effector_regulon_dotplot.{png,pdf}
#   figures/12_pySCENIC_visualization/pyscenic_zol_foxp3_vd2_function_gene_heatmap.{png,pdf}
#   figures/12_pySCENIC_visualization/pyscenic_zol_foxp3_vd2_function_gene_dotplot.{png,pdf}
#   figures/12_pySCENIC_visualization/pyscenic_sox4_target_expression_heatmap.{png,pdf}
#   table/pyscenic/pyscenic_regulon_activity_by_group.csv
#   table/pyscenic/pyscenic_selected_regulons.txt
#   table/pyscenic/pyscenic_curated_tf_interpretation.csv
#   table/pyscenic/pyscenic_vd2_effector_regulon_delta.csv
#   table/pyscenic/pyscenic_zol_foxp3_vd2_function_gene_summary.csv
#   table/pyscenic/pyscenic_zol_foxp3_vd2_function_gene_delta.csv
#   table/pyscenic/pyscenic_sox4_target_associations.csv
#   rds/pyscenic_visualization_result.rds
#
# Existing PNG/PDF/RDS outputs are reused by default. Pass --overwrite after
# changing grouping, regulon counts, or source inputs. Pass --overwrite-plots
# to redraw figures without overwriting the cached result RDS.

suppressPackageStartupMessages({
  library(hdf5r)
  library(Seurat)
  library(Matrix)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(ggplot2)
  library(patchwork)
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
})

options(stringsAsFactors = FALSE)
force_pyscenic_visualization <- FALSE
force_pyscenic_visualization_plot <- TRUE

source_plotting_shared <- function() {
  candidates <- c(
    "code/downstream/lib/plotting_shared.R",
    "lib/plotting_shared.R",
    "_cache_/plotting_shared.R",
    "cache/plotting_shared.R"
  )
  plotting_shared <- candidates[file.exists(candidates)][1]
  if (!is.na(plotting_shared)) {
    source(plotting_shared)
  }
}

infer_project_dir <- function() {
  command_args <- commandArgs(trailingOnly = FALSE)
  file_arg <- "--file="
  script_path <- sub(file_arg, "", command_args[startsWith(command_args, file_arg)][1])
  if (!is.na(script_path) && nzchar(script_path)) {
    script_path <- normalizePath(script_path, mustWork = FALSE)
    script_dir <- dirname(script_path)
    if (basename(script_dir) == "downstream_scripts") {
      return(dirname(script_dir))
    }
  }
  getwd()
}

get_group_palette <- function(group_by, group_values) {
  group_values <- unique(as.character(group_values))
  if (group_by == "cell_type" && exists("color_celltype", inherits = TRUE)) {
    palette <- get("color_celltype", inherits = TRUE)
  } else if (group_by == "group" && exists("color_group", inherits = TRUE)) {
    palette <- get("color_group", inherits = TRUE)
  } else if (group_by == "sample_name" && exists("color_sample_name", inherits = TRUE)) {
    palette <- get("color_sample_name", inherits = TRUE)
  } else {
    return(NULL)
  }
  palette[intersect(names(palette), group_values)]
}

add_group_color_scale <- function(plot, group_palette) {
  if (is.null(group_palette) || length(group_palette) == 0L) {
    return(plot)
  }
  plot + scale_color_manual(values = group_palette, drop = FALSE)
}

# Parse command-line options without requiring an additional CLI package.
# Paths are interpreted relative to project_dir after parsing.
parse_args <- function(args) {
  config <- list(
    project_dir = infer_project_dir(),
    auc_loom = "pyscenic_10k/auc_mtx.loom",
    regulons_csv = "pyscenic_10k/regulons.csv",
    adj_tsv = "pyscenic_10k/adjacencies.tsv",
    seurat_rds = "rds/all_seurat_celltype.rds",
    assay = "RNA",
    expr_layer = "data",
    group_by = "cell_type",
    reduction = "umap.unintegrated",
    n_top_regulons = 30L,
    sox4_top_targets = 30L,
    figure_dir = file.path("figures", "12_pySCENIC_visualization"),
    table_dir = "table/pyscenic",
    result_rds = "rds/pyscenic_visualization_result.rds",
    overwrite = force_pyscenic_visualization,
    overwrite_plots = force_pyscenic_visualization_plot
  )

  i <- 1L
  while (i <= length(args)) {
    flag <- args[[i]]
    if (flag %in% c("-h", "--help")) {
      cat(
        "Usage: Rscript code/downstream/12_pySCENIC_visualization.R [OPTIONS]\n\n",
        "  --project-dir DIR\n",
        "  --auc-loom PATH\n",
        "  --regulons-csv PATH\n",
        "  --adj-tsv PATH\n",
        "  --seurat-rds PATH\n",
        "  --assay NAME\n",
        "  --expr-layer NAME\n",
        "  --group-by COLUMN\n",
        "  --reduction NAME\n",
        "  --n-top-regulons INT\n",
        "  --sox4-top-targets INT\n",
        "  --figure-dir DIR\n",
        "  --table-dir DIR\n",
        "  --result-rds PATH\n",
        "  --overwrite\n",
        "  --overwrite-plots\n",
        sep = ""
      )
      quit(status = 0)
    }
    if (flag == "--overwrite") {
      config$overwrite <- TRUE
      config$overwrite_plots <- TRUE
      i <- i + 1L
      next
    }
    if (flag == "--overwrite-plots") {
      config$overwrite_plots <- TRUE
      i <- i + 1L
      next
    }
    if (i == length(args)) {
      stop("Missing value for ", flag, call. = FALSE)
    }
    value <- args[[i + 1L]]
    key <- switch(flag,
      "--project-dir" = "project_dir",
      "--auc-loom" = "auc_loom",
      "--regulons-csv" = "regulons_csv",
      "--adj-tsv" = "adj_tsv",
      "--seurat-rds" = "seurat_rds",
      "--assay" = "assay",
      "--expr-layer" = "expr_layer",
      "--group-by" = "group_by",
      "--reduction" = "reduction",
      "--n-top-regulons" = "n_top_regulons",
      "--sox4-top-targets" = "sox4_top_targets",
      "--figure-dir" = "figure_dir",
      "--table-dir" = "table_dir",
      "--result-rds" = "result_rds",
      stop("Unknown argument: ", flag, call. = FALSE)
    )
    config[[key]] <- value
    i <- i + 2L
  }

  config$n_top_regulons <- as.integer(config$n_top_regulons)
  config$sox4_top_targets <- as.integer(config$sox4_top_targets)
  if (is.na(config$n_top_regulons) || config$n_top_regulons < 1L) {
    stop("--n-top-regulons must be a positive integer.", call. = FALSE)
  }
  if (is.na(config$sox4_top_targets) || config$sox4_top_targets < 1L) {
    stop("--sox4-top-targets must be a positive integer.", call. = FALSE)
  }
  config$project_dir <- normalizePath(config$project_dir, mustWork = FALSE)
  for (path_key in c("auc_loom", "regulons_csv", "adj_tsv", "seurat_rds", "figure_dir", "table_dir", "result_rds")) {
    if (!grepl("^/", config[[path_key]])) {
      config[[path_key]] <- file.path(config$project_dir, config[[path_key]])
    }
  }
  config
}

# Reject missing and zero-byte inputs before loading large R packages/objects.
check_file <- function(path, label) {
  if (!file.exists(path) || file.info(path)$size == 0) {
    extra <- ""
    if (grepl("auc_mtx[.]loom$", path)) {
      extra <- paste0(
        "\nRun pySCENIC first, for example:\n",
        "  bash scripts/run_downstream.sh --existing-run <run_id> --step 11 ",
        "--tf-list <TF_LIST> --ranking-db <RANKING_DB> --motif-annotations <MOTIF_ANNOTATIONS>\n",
        "Then rerun step 12."
      )
    }
    stop("Missing or empty ", label, ": ", normalizePath(path, mustWork = FALSE), extra, call. = FALSE)
  }
}

# HDF5 files begin with a fixed eight-byte signature. Checking it here gives
# a useful message when a failed AUCell command leaves text, gzip data, or a
# truncated file named auc_mtx.loom.
validate_hdf5_signature <- function(path) {
  connection <- file(path, open = "rb")
  on.exit(close(connection), add = TRUE)
  signature <- readBin(connection, what = "raw", n = 8L)
  hdf5_signature <- as.raw(c(0x89, 0x48, 0x44, 0x46, 0x0d, 0x0a, 0x1a, 0x0a))
  if (identical(signature, hdf5_signature)) {
    return(invisible(TRUE))
  }

  if (length(signature) >= 2L && identical(signature[1:2], as.raw(c(0x1f, 0x8b)))) {
    stop(
      "AUCell output is gzip-compressed, not a directly readable loom file: ", path, "\n",
      "Decompress it first or rerun pySCENIC AUCell with an uncompressed .loom output.",
      call. = FALSE
    )
  }

  preview_connection <- file(path, open = "rb")
  on.exit(close(preview_connection), add = TRUE)
  preview_raw <- readBin(preview_connection, what = "raw", n = 160L)
  preview_bytes <- as.integer(preview_raw)
  preview <- paste0(
    vapply(
      preview_bytes,
      function(byte) if (byte == 9L || (byte >= 32L && byte <= 126L)) intToUtf8(byte) else " ",
      character(1)
    ),
    collapse = ""
  )
  preview <- trimws(preview)
  preview_message <- if (nzchar(preview)) paste0("\nFile begins with: ", preview) else ""
  stop(
    "AUCell output is not an HDF5/loom file: ", path, preview_message, "\n",
    "Inspect pyscenic_10k/logs/pyscenic_aucell.log, remove the invalid ",
    "auc_mtx.loom, and rerun code/downstream/11_pySCENIC.sh.",
    call. = FALSE
  )
}

# hdf5r may return fixed-length HDF5 strings with padding; convert them to
# ordinary trimmed R character vectors before matching cell IDs.
decode_hdf5_strings <- function(values) {
  trimws(as.character(as.vector(values)))
}

# Read the two pySCENIC loom attributes needed for downstream visualization.
# pySCENIC stores RegulonsAUC as a compound HDF5 field in common loom outputs,
# but hdf5r can expose it as a data frame, matrix, or named list depending on
# package/HDF5 versions. All supported representations are normalized here.
#
# Returned orientation is always:
#   rows    = regulons
#   columns = cells
read_regulon_auc <- function(loom_path) {
  validate_hdf5_signature(loom_path)
  loom <- H5File$new(loom_path, mode = "r")
  on.exit(loom$close_all(), add = TRUE)

  if (!loom$exists("col_attrs/CellID")) {
    stop("Loom is missing /col_attrs/CellID: ", loom_path, call. = FALSE)
  }
  if (!loom$exists("col_attrs/RegulonsAUC")) {
    stop("Loom is missing /col_attrs/RegulonsAUC: ", loom_path, call. = FALSE)
  }

  cell_ids <- decode_hdf5_strings(loom[["col_attrs/CellID"]][])
  auc_raw <- loom[["col_attrs/RegulonsAUC"]][]

  if (is.data.frame(auc_raw)) {
    auc_cells_by_regulon <- as.matrix(auc_raw)
  } else if (is.matrix(auc_raw)) {
    auc_cells_by_regulon <- auc_raw
  } else if (is.list(auc_raw) && !is.null(names(auc_raw))) {
    auc_cells_by_regulon <- do.call(cbind, auc_raw)
    colnames(auc_cells_by_regulon) <- names(auc_raw)
  } else {
    stop(
      "Unsupported /col_attrs/RegulonsAUC representation: ",
      paste(class(auc_raw), collapse = ", "),
      call. = FALSE
    )
  }

  if (nrow(auc_cells_by_regulon) != length(cell_ids) &&
    ncol(auc_cells_by_regulon) == length(cell_ids)) {
    auc_cells_by_regulon <- t(auc_cells_by_regulon)
  }
  if (nrow(auc_cells_by_regulon) != length(cell_ids)) {
    stop(
      "RegulonsAUC cell dimension does not match CellID: ",
      nrow(auc_cells_by_regulon), " versus ", length(cell_ids),
      call. = FALSE
    )
  }
  if (is.null(colnames(auc_cells_by_regulon))) {
    stop("RegulonsAUC has no regulon names in the loom compound fields.", call. = FALSE)
  }

  rownames(auc_cells_by_regulon) <- cell_ids
  auc_matrix <- t(auc_cells_by_regulon)
  storage.mode(auc_matrix) <- "numeric"
  auc_matrix
}

# Match by cell names, never by column position. The original Seurat UMAP and
# metadata are retained only for cells present in the AUCell result.
align_auc_and_seurat <- function(auc_matrix, seurat_obj, group_by, reduction) {
  if (!group_by %in% colnames(seurat_obj@meta.data)) {
    stop("Seurat metadata column not found: ", group_by, call. = FALSE)
  }
  if (!reduction %in% names(seurat_obj@reductions)) {
    stop(
      "Seurat reduction not found: ", reduction,
      ". Available: ", paste(names(seurat_obj@reductions), collapse = ", "),
      call. = FALSE
    )
  }

  common_cells <- intersect(colnames(auc_matrix), colnames(seurat_obj))
  if (length(common_cells) == 0L) {
    stop("No common cell IDs between AUCell loom and Seurat object.", call. = FALSE)
  }

  auc_matrix <- auc_matrix[, common_cells, drop = FALSE]
  seurat_obj <- subset(seurat_obj, cells = common_cells)
  seurat_obj <- seurat_obj[, common_cells]

  groups <- as.character(seurat_obj@meta.data[[group_by]])
  names(groups) <- colnames(seurat_obj)
  keep <- !is.na(groups) & groups != ""
  if (!any(keep)) {
    stop("No cells have a non-empty value in metadata column: ", group_by, call. = FALSE)
  }

  cells <- names(groups)[keep]
  list(
    auc = auc_matrix[, cells, drop = FALSE],
    seurat = seurat_obj[, cells],
    groups = groups[cells]
  )
}

# Summarize activity by the selected metadata group.
#
# mean_matrix:
#   Mean AUCell AUC for each regulon in each group.
# z_matrix:
#   Per-regulon z-score across group means. This highlights relative group
#   enrichment and should not be interpreted as an absolute activity value.
# specificity:
#   max(group z-score) - min(group z-score). Regulons with the largest range
#   are selected for plotting because they best distinguish the groups.
summarise_regulons <- function(auc_matrix, groups) {
  group_levels <- unique(groups)
  mean_matrix <- vapply(
    group_levels,
    function(group_name) rowMeans(auc_matrix[, groups == group_name, drop = FALSE]),
    numeric(nrow(auc_matrix))
  )
  if (is.null(dim(mean_matrix))) {
    mean_matrix <- matrix(mean_matrix, ncol = 1L)
  }
  rownames(mean_matrix) <- rownames(auc_matrix)
  colnames(mean_matrix) <- group_levels

  row_sd <- apply(mean_matrix, 1, sd)
  row_sd[is.na(row_sd) | row_sd == 0] <- 1
  z_matrix <- sweep(mean_matrix, 1, rowMeans(mean_matrix), "-")
  z_matrix <- sweep(z_matrix, 1, row_sd, "/")

  specificity <- apply(z_matrix, 1, function(x) max(x) - min(x))
  top_regulons <- names(sort(specificity, decreasing = TRUE))

  mean_stats <- as.data.frame(mean_matrix) %>%
    tibble::rownames_to_column("regulon") %>%
    pivot_longer(-regulon, names_to = "group", values_to = "mean_auc")
  z_stats <- as.data.frame(z_matrix) %>%
    tibble::rownames_to_column("regulon") %>%
    pivot_longer(-regulon, names_to = "group", values_to = "zscore")
  stats <- mean_stats %>%
    left_join(z_stats, by = c("regulon", "group")) %>%
    mutate(specificity = specificity[regulon])

  list(
    mean = mean_matrix,
    zscore = z_matrix,
    specificity = specificity,
    top_regulons = top_regulons,
    stats = stats
  )
}

# Prepare the dot-plot summary. "Active" is defined descriptively as a cell
# whose AUC is above the global median for that regulon. This is not a formal
# AUCell binary threshold and is used only to provide a prevalence measure.
make_dotplot_data <- function(auc_matrix, groups, regulons) {
  thresholds <- apply(auc_matrix[regulons, , drop = FALSE], 1, median)
  bind_rows(lapply(regulons, function(regulon) {
    bind_rows(lapply(unique(groups), function(group_name) {
      values <- auc_matrix[regulon, groups == group_name]
      tibble(
        regulon = regulon,
        group = group_name,
        mean_auc = mean(values),
        percent_active = mean(values > thresholds[[regulon]]) * 100
      )
    }))
  }))
}

vd2_effector_compare_states <- c(
  "ZOL Effector Vd2",
  "ZOL FOXP3+ Vd2",
  "PAN Effector Vd2"
)

vd2_effector_curated_tfs <- c(
  "IRF1", "FOXP3", "NR3C1", "NFATC3", "GABPB1",
  "TFEB", "XBP1", "CEBPB", "TFDP1", "STAT3"
)

make_vd2_effector_comparison <- function(
  auc_matrix,
  seurat_obj,
  states = vd2_effector_compare_states,
  curated_tfs = vd2_effector_curated_tfs,
  top_n = 30L
) {
  if (!"cell_type" %in% colnames(seurat_obj@meta.data)) {
    message("[SKIP] Seurat metadata has no cell_type column for Vd2 effector pySCENIC comparison.")
    return(NULL)
  }

  cell_types <- as.character(seurat_obj@meta.data[colnames(auc_matrix), "cell_type", drop = TRUE])
  names(cell_types) <- colnames(auc_matrix)
  present_states <- states[states %in% unique(cell_types)]
  if (length(present_states) < 2L) {
    message(
      "[SKIP] Need at least two target Vd2 states for pySCENIC comparison. Present: ",
      paste(present_states, collapse = ", ")
    )
    return(NULL)
  }

  target_cells <- names(cell_types)[cell_types %in% present_states]
  target_auc <- auc_matrix[, target_cells, drop = FALSE]
  target_groups <- factor(cell_types[target_cells], levels = states)
  target_groups <- droplevels(target_groups)
  names(target_groups) <- target_cells

  summary_result <- summarise_regulons(target_auc, as.character(target_groups))
  curated_table <- make_curated_tf_annotations() %>%
    filter(tf %in% curated_tfs) %>%
    match_curated_regulons(rownames(target_auc)) %>%
    filter(detected)
  curated_regulons <- curated_table$matched_regulon
  variable_regulons <- summary_result$top_regulons[
    seq_len(min(top_n, length(summary_result$top_regulons)))
  ]
  selected_regulons <- unique(c(curated_regulons, variable_regulons))
  selected_regulons <- selected_regulons[selected_regulons %in% rownames(target_auc)]
  if (length(selected_regulons) == 0L) {
    message("[SKIP] No regulons available for Vd2 effector pySCENIC comparison.")
    return(NULL)
  }

  labels <- tibble(
    regulon = selected_regulons,
    tf = get_regulon_tf(selected_regulons)
  ) %>%
    mutate(
      display_label = ifelse(tf %in% curated_table$tf, paste0(tf, " *"), tf),
      display_label = make.unique(display_label)
    )
  heatmap_matrix <- summary_result$zscore[selected_regulons, levels(target_groups), drop = FALSE]
  rownames(heatmap_matrix) <- labels$display_label

  dotplot_data <- make_dotplot_data(target_auc, as.character(target_groups), selected_regulons) %>%
    left_join(labels, by = "regulon") %>%
    mutate(
      group = factor(group, levels = states),
      tf = factor(tf, levels = rev(unique(labels$tf)))
    )

  pairwise_comparisons <- list(
    zol_effector_vs_pan_effector = c("ZOL Effector Vd2", "PAN Effector Vd2"),
    zol_foxp3_vs_pan_effector = c("ZOL FOXP3+ Vd2", "PAN Effector Vd2"),
    zol_foxp3_vs_zol_effector = c("ZOL FOXP3+ Vd2", "ZOL Effector Vd2")
  )
  delta_table <- bind_rows(lapply(names(pairwise_comparisons), function(comparison_name) {
    comparison <- pairwise_comparisons[[comparison_name]]
    ident_1 <- comparison[[1]]
    ident_2 <- comparison[[2]]
    if (!all(c(ident_1, ident_2) %in% as.character(target_groups))) {
      return(tibble())
    }
    ident_1_cells <- names(target_groups)[as.character(target_groups) == ident_1]
    ident_2_cells <- names(target_groups)[as.character(target_groups) == ident_2]
    tibble(
      comparison = comparison_name,
      regulon = rownames(target_auc),
      tf = get_regulon_tf(rownames(target_auc)),
      ident_1 = ident_1,
      ident_2 = ident_2,
      mean_ident_1 = rowMeans(target_auc[, ident_1_cells, drop = FALSE]),
      mean_ident_2 = rowMeans(target_auc[, ident_2_cells, drop = FALSE]),
      n_ident_1 = length(ident_1_cells),
      n_ident_2 = length(ident_2_cells)
    ) %>%
      mutate(
        delta = mean_ident_1 - mean_ident_2,
        abs_delta = abs(delta),
        higher_in = ifelse(delta >= 0, ident_1, ident_2),
        curated = tf %in% curated_tfs
      )
  })) %>%
    arrange(comparison, desc(curated), desc(abs_delta), regulon)

  top_delta_regulons <- delta_table %>%
    group_by(comparison) %>%
    slice_max(order_by = abs_delta, n = 12, with_ties = FALSE) %>%
    ungroup() %>%
    distinct(regulon)

  list(
    states = levels(target_groups),
    cells = target_cells,
    selected_regulons = selected_regulons,
    curated_table = curated_table,
    heatmap_matrix = heatmap_matrix,
    dotplot_data = dotplot_data,
    delta_table = delta_table,
    delta_plot_data = delta_table %>% semi_join(top_delta_regulons, by = "regulon")
  )
}

plot_vd2_effector_dotplot <- function(dotplot_data) {
  ggplot(dotplot_data, aes(group, tf)) +
    geom_point(aes(size = percent_active, color = mean_auc)) +
    scale_color_viridis_c(option = "plasma") +
    scale_size(range = c(1.5, 8)) +
    theme_bw() +
    theme(
      axis.title = element_blank(),
      axis.text.x = element_text(angle = 35, hjust = 1),
      axis.text.y = element_text(face = "bold"),
      panel.grid = element_line(linewidth = 0.2, color = "grey90")
    ) +
    labs(
      title = "Regulon activity in ZOL/PAN effector Vd2 states",
      subtitle = "Curated TFs are prioritized, with top variable regulons added from the all-cell pySCENIC result",
      color = "Mean AUC",
      size = "% above median"
    )
}

zol_foxp3_vd2_states <- c(
  "Effector Memory Vd2",
  "Pre-activated Vd2",
  "ZOL FOXP3+ Vd2",
  "ZOL Effector Vd2",
  "PAN Effector Vd2"
)

get_all_vd2_states <- function(cell_types, preferred_states = zol_foxp3_vd2_states) {
  all_states <- unique(as.character(cell_types))
  all_states <- all_states[!is.na(all_states) & grepl("Vd2", all_states)]
  if (length(all_states) == 0L) {
    return(character())
  }
  if (exists("cell_type_levels", inherits = TRUE)) {
    level_order <- get("cell_type_levels", inherits = TRUE)
    ordered_states <- level_order[level_order %in% all_states]
    ordered_states <- c(ordered_states, setdiff(all_states, ordered_states))
  } else {
    ordered_states <- sort(all_states)
  }
  c(preferred_states[preferred_states %in% ordered_states], setdiff(ordered_states, preferred_states))
}

zol_foxp3_function_gene_sets <- list(
  "FOXP3 regulatory program" = c(
    "FOXP3", "IL2RA", "CTLA4", "TIGIT", "IKZF2", "IKZF4",
    "TNFRSF18", "TNFRSF4", "CCR8", "ENTPD1", "IL10", "TGFB1"
  ),
  "Activation checkpoint" = c(
    "ICOS", "PDCD1", "LAG3", "HAVCR2", "TOX", "BATF", "IRF1",
    "NR3C1", "CXCR3", "CXCR4"
  ),
  "Effector cytotoxic context" = c(
    "NKG7", "GNLY", "PRF1", "GZMA", "GZMB", "GZMH",
    "IFNG", "TNF", "CCL3", "CCL4", "FASLG"
  ),
  "Stress secretory fitness" = c(
    "XBP1", "CEBPB", "HIF1A", "SLC2A1", "MKI67", "TOP2A"
  )
)

get_gene_function_map <- function(gene_sets) {
  bind_rows(lapply(names(gene_sets), function(function_group) {
    tibble(gene = gene_sets[[function_group]], function_group = function_group)
  })) %>%
    distinct(gene, .keep_all = TRUE)
}

summarise_function_gene_expression <- function(expression_matrix, groups, gene_function_map) {
  present_genes <- intersect(gene_function_map$gene, rownames(expression_matrix))
  if (length(present_genes) == 0L) {
    return(tibble())
  }
  bind_rows(lapply(present_genes, function(gene) {
    expression <- dense_feature_vector(expression_matrix, gene)
    bind_rows(lapply(unique(groups), function(group_name) {
      cells <- names(groups)[groups == group_name]
      values <- expression[cells]
      tibble(
        gene = gene,
        group = group_name,
        n_cells = length(values),
        mean_expression = mean(values, na.rm = TRUE),
        median_expression = median(values, na.rm = TRUE),
        percent_expressing = mean(values > 0, na.rm = TRUE) * 100
      )
    }))
  })) %>%
    left_join(gene_function_map, by = "gene") %>%
    relocate(function_group, .after = gene)
}

calculate_function_gene_delta <- function(expression_summary, ident_1 = "ZOL FOXP3+ Vd2") {
  if (nrow(expression_summary) == 0L) {
    return(tibble())
  }
  references <- setdiff(unique(expression_summary$group), ident_1)
  bind_rows(lapply(references, function(ident_2) {
    ident_1_data <- expression_summary %>%
      filter(group == ident_1) %>%
      select(gene, function_group, mean_ident_1 = mean_expression, pct_ident_1 = percent_expressing, n_ident_1 = n_cells)
    ident_2_data <- expression_summary %>%
      filter(group == ident_2) %>%
      select(gene, mean_ident_2 = mean_expression, pct_ident_2 = percent_expressing, n_ident_2 = n_cells)
    ident_1_data %>%
      inner_join(ident_2_data, by = "gene") %>%
      mutate(
        comparison = paste0(make.names(ident_1), "_vs_", make.names(ident_2)),
        ident_1 = ident_1,
        ident_2 = ident_2,
        delta_mean_expression = mean_ident_1 - mean_ident_2,
        delta_percent_expressing = pct_ident_1 - pct_ident_2,
        abs_delta_mean_expression = abs(delta_mean_expression),
        higher_in = ifelse(delta_mean_expression >= 0, ident_1, ident_2),
        .before = 1
      )
  })) %>%
    arrange(comparison, desc(abs_delta_mean_expression), gene)
}

make_function_gene_heatmap <- function(expression_summary, states = zol_foxp3_vd2_states) {
  if (nrow(expression_summary) == 0L) {
    return(NULL)
  }
  function_levels <- names(zol_foxp3_function_gene_sets)
  mean_matrix <- expression_summary %>%
    mutate(group = factor(group, levels = states)) %>%
    select(gene, function_group, group, mean_expression) %>%
    tidyr::pivot_wider(names_from = group, values_from = mean_expression, values_fill = 0) %>%
    arrange(function_group, gene)
  row_info <- mean_matrix %>% select(gene, function_group)
  matrix_data <- mean_matrix %>%
    select(-gene, -function_group) %>%
    as.matrix()
  rownames(matrix_data) <- row_info$gene
  matrix_data <- matrix_data[, intersect(states, colnames(matrix_data)), drop = FALSE]
  row_sd <- apply(matrix_data, 1, sd)
  row_sd[is.na(row_sd) | row_sd == 0] <- 1
  z_matrix <- sweep(matrix_data, 1, rowMeans(matrix_data), "-")
  z_matrix <- sweep(z_matrix, 1, row_sd, "/")
  row_function <- factor(row_info$function_group, levels = function_levels)
  names(row_function) <- row_info$gene
  module_colors <- setNames(
    c("#B2182B", "#E69F00", "#2166AC", "#4D9221")[seq_along(function_levels)],
    function_levels
  )
  row_ha <- rowAnnotation(
    module = row_function,
    col = list(module = module_colors),
    show_annotation_name = TRUE
  )

  Heatmap(
    z_matrix,
    name = "expr z-score",
    col = colorRamp2(c(-2, 0, 2), c("#2166AC", "#F7F7F7", "#B2182B")),
    cluster_rows = TRUE,
    cluster_columns = FALSE,
    show_row_names = TRUE,
    row_names_gp = gpar(fontsize = 7),
    column_names_gp = gpar(fontsize = 9),
    left_annotation = row_ha,
    row_title = "FOXP3-related genes",
    column_title = "Mean FOXP3-related gene expression by Vd2 state"
  )
}

plot_function_gene_dotplot <- function(
  expression_summary,
  seurat_obj,
  cells,
  states = zol_foxp3_vd2_states,
  gene_sets = zol_foxp3_function_gene_sets
) {
  if (nrow(expression_summary) == 0L) {
    return(NULL)
  }
  function_levels <- names(gene_sets)
  gene_levels <- rev(unlist(gene_sets[function_levels], use.names = FALSE))
  gene_levels <- gene_levels[gene_levels %in% expression_summary$gene]
  plot_data <- expression_summary %>%
    mutate(
      group = factor(group, levels = states),
      function_group = factor(function_group, levels = function_levels),
      gene = factor(gene, levels = gene_levels)
    )
  ggplot(plot_data, aes(group, gene)) +
    geom_point(aes(size = percent_expressing, color = mean_expression)) +
    facet_grid(function_group ~ ., scales = "free_y", space = "free_y") +
    scale_color_viridis_c(option = "plasma") +
    scale_size(range = c(0.8, 6.5), limits = c(0, 100)) +
    theme_bw() +
    theme(
      axis.title = element_blank(),
      axis.text.x = element_text(angle = 35, hjust = 1),
      strip.text.y = element_text(angle = 0, hjust = 0),
      panel.spacing.y = unit(0.08, "in"),
      panel.grid = element_line(linewidth = 0.2, color = "grey90")
    ) +
    labs(
      title = "FOXP3-related function-gene expression across Vd2 states",
      color = "Mean expression",
      size = "% expressing"
    )
}

make_zol_foxp3_function_gene_analysis <- function(
  auc_matrix,
  seurat_obj,
  assay,
  layer,
  states = zol_foxp3_vd2_states,
  gene_sets = zol_foxp3_function_gene_sets
) {
  if (!"cell_type" %in% colnames(seurat_obj@meta.data)) {
    message("[SKIP] Seurat metadata has no cell_type column for ZOL FOXP3+ function-gene analysis.")
    return(NULL)
  }
  cell_types <- as.character(seurat_obj@meta.data[colnames(auc_matrix), "cell_type", drop = TRUE])
  names(cell_types) <- colnames(auc_matrix)
  states <- get_all_vd2_states(cell_types, preferred_states = states)
  present_states <- states[states %in% unique(cell_types)]
  if (!"ZOL FOXP3+ Vd2" %in% present_states || length(present_states) < 2L) {
    message("[SKIP] FOXP3-related Vd2 function-gene analysis needs ZOL FOXP3+ Vd2 and at least one other Vd2 state.")
    return(NULL)
  }
  target_cells <- names(cell_types)[cell_types %in% present_states]
  target_groups <- factor(cell_types[target_cells], levels = states)
  names(target_groups) <- target_cells
  target_group_values <- as.character(target_groups)
  names(target_group_values) <- target_cells
  gene_function_map <- get_gene_function_map(gene_sets)
  expression_matrix <- get_expression_matrix(
    seurat_obj = seurat_obj,
    assay = assay,
    layer = layer,
    features = gene_function_map$gene,
    cells = target_cells
  )
  expression_summary <- summarise_function_gene_expression(
    expression_matrix,
    target_group_values,
    gene_function_map
  )
  delta_table <- calculate_function_gene_delta(expression_summary, ident_1 = "ZOL FOXP3+ Vd2")

  list(
    cells = target_cells,
    groups = target_groups,
    states = present_states,
    expression_matrix = expression_matrix,
    expression_summary = expression_summary,
    delta_table = delta_table,
    heatmap = make_function_gene_heatmap(expression_summary, states = states),
    dotplot = plot_function_gene_dotplot(
      expression_summary,
      seurat_obj = seurat_obj,
      cells = target_cells,
      states = states,
      gene_sets = gene_sets
    )
  )
}

# Curated TF regulons selected for biological interpretation. These
# descriptions are project-specific annotations supplied for the Vd1/Vd2
# analysis and are not inferred by pySCENIC itself.
make_curated_tf_annotations <- function() {
  tibble::tribble(
    ~tf, ~highest_associated_cell_types, ~established_t_cell_function, ~predicted_interpretation,
    "BACH2", "Naive Vd1",
    "Maintains naive/quiescent states; restrains terminal effector differentiation",
    "Preservation of an undifferentiated Vd1 program",
    "SOX4", "Naive Vd1",
    "Supports thymocyte development and early T-cell differentiation",
    "Developmental or naive-cell identity",
    "NFATC1", "Naive Vd1",
    "Mediates TCR-dependent activation; can support activation, tolerance, or exhaustion depending on partners",
    "TCR-responsive but non-terminal activation state",
    "MTF1", "Naive Vd1",
    "Regulates metal-ion and oxidative-stress responses; emerging role in immune-cell fitness",
    "Cellular stress adaptation in naive cells",
    "TFEB", "Effector Memory Vd2",
    "Controls lysosomal biogenesis, autophagy, and metabolic adaptation",
    "Enhanced lysosomal and memory-associated metabolic fitness",
    "NR3C1", "ZOL Effector Vd2; Pre-activated Vd2",
    "Glucocorticoid receptor that suppresses inflammatory and TCR-driven responses",
    "Feedback control of activation and inflammatory stress",
    "GABPB1", "ZOL Effector Vd2",
    "Supports mitochondrial gene expression, proliferation, and lymphocyte fitness",
    "Increased energetic requirements during effector activation",
    "NFATC3", "ZOL Effector Vd2",
    "Participates in TCR/calcium-dependent transcription and cytokine responses",
    "Strong antigen-receptor-dependent effector activation",
    "IRF1", "ZOL Effector Vd2; ZOL FOXP3+ Vd2; Pre-activated Vd2",
    "Regulates interferon responses, antigen presentation, and inflammatory T-cell differentiation",
    "Interferon-responsive activation induced by ZOL treatment",
    "FOXP3", "ZOL FOXP3+ Vd2",
    "Master regulator of regulatory T-cell differentiation and suppressive function",
    "Potential regulatory or activation-limiting Vd2 state",
    "STAT3", "Effector Vd1",
    "Transduces IL-6, IL-21, and IL-23 signals; supports survival and inflammatory differentiation",
    "Cytokine-driven inflammatory effector Vd1 program",
    "XBP1", "Effector Vd1; PAN Effector Vd2",
    "Coordinates the unfolded-protein response, metabolism, and sustained effector function",
    "Adaptation to high biosynthetic and secretory demands",
    "CEBPB", "Effector Vd1; PAN Effector Vd2",
    "Regulates inflammatory, metabolic, and stress-response genes in activated immune cells",
    "Terminal inflammatory-effector differentiation",
    "TFDP1", "PAN Effector Vd2",
    "E2F partner controlling cell-cycle progression; not T-cell-specific",
    "Increased proliferation within PAN effector Vd2 cells",
    "HMGA1", "Pre-activated Vd1",
    "Chromatin regulator associated with lymphocyte activation and proliferation",
    "Early chromatin remodeling during Vd1 activation"
  )
}

# pySCENIC commonly names regulons as TF(+), TF(-), or TF_extended(+).
# Match annotations by the TF prefix and prefer a positive regulon when both
# positive and negative regulons are available.
match_curated_regulons <- function(annotations, regulon_names) {
  regulon_index <- tibble(
    regulon = regulon_names,
    tf = sub("_extended$", "", sub("\\(.*$", "", regulon_names)),
    positive = grepl("\\(\\+\\)$", regulon_names),
    extended = grepl("_extended", regulon_names)
  ) %>%
    arrange(tf, desc(positive), extended, regulon)

  matched <- regulon_index %>%
    filter(tf %in% annotations$tf) %>%
    distinct(tf, .keep_all = TRUE) %>%
    transmute(
      tf,
      matched_regulon = regulon,
      regulon_sign = case_when(
        grepl("\\(\\+\\)$", regulon) ~ "+",
        grepl("\\(-\\)$", regulon) ~ "-",
        TRUE ~ NA_character_
      )
    )

  annotations %>%
    left_join(matched, by = "tf") %>%
    mutate(detected = !is.na(matched_regulon))
}

# Save both raster and vector forms. Existing non-empty pairs are checkpoints
# unless --overwrite is supplied.
save_ggplot <- function(plot, stem, width, height, overwrite) {
  paths <- paste0(stem, c(".png", ".pdf"))
  if (!overwrite && all(file.exists(paths) & file.info(paths)$size > 0)) {
    message("[SKIP] Existing plot: ", paste(paths, collapse = ", "))
    return(invisible(FALSE))
  }
  ggsave(paths[[1]], plot, width = width, height = height, dpi = 300)
  ggsave(paths[[2]], plot, width = width, height = height)
  invisible(TRUE)
}

# ComplexHeatmap objects are drawn directly to graphics devices rather than
# passed to ggsave. The same checkpoint behavior is used as for ggplot output.
save_complex_heatmap <- function(heatmap, stem, width, height, overwrite) {
  paths <- paste0(stem, c(".png", ".pdf"))
  if (!overwrite && all(file.exists(paths) & file.info(paths)$size > 0)) {
    message("[SKIP] Existing heatmap: ", paste(paths, collapse = ", "))
    return(invisible(FALSE))
  }
  png(paths[[1]], width = width, height = height, units = "in", res = 300)
  ComplexHeatmap::draw(heatmap)
  dev.off()
  pdf(paths[[2]], width = width, height = height)
  ComplexHeatmap::draw(heatmap)
  dev.off()
  invisible(TRUE)
}

# ---------------------------------------------------------------------------
# SOX4 target relationship visualizations
# ---------------------------------------------------------------------------

# Extract the TF symbol prefix from pySCENIC regulon names such as SOX4(+),
# SOX4(-), or SOX4_extended(+).
get_regulon_tf <- function(regulons) {
  sub("_extended$", "", sub("\\(.*$", "", regulons))
}

choose_tf_regulon <- function(auc_matrix, tf = "SOX4") {
  regulon_index <- tibble(
    regulon = rownames(auc_matrix),
    tf_symbol = get_regulon_tf(rownames(auc_matrix)),
    positive = grepl("\\(\\+\\)$", rownames(auc_matrix)),
    extended = grepl("_extended", rownames(auc_matrix))
  ) %>%
    filter(tf_symbol == tf) %>%
    arrange(desc(positive), extended, regulon)

  if (nrow(regulon_index) == 0L) {
    return(NA_character_)
  }
  regulon_index$regulon[[1]]
}

# pySCENIC ctx output differs slightly by version. TargetGenes is often a
# serialized Python list of target tuples, but it may also be a simpler string.
# This parser intentionally keeps only plausible gene symbols and removes
# common serialization words/numbers.
extract_target_tokens <- function(values) {
  text <- paste(values, collapse = " ")
  hits <- gregexpr("[A-Za-z][A-Za-z0-9._-]*", text, perl = TRUE)
  tokens <- regmatches(text, hits)[[1]]
  tokens <- unique(tokens)
  tokens <- tokens[!tokens %in% c(
    "NA", "NaN", "None", "True", "False", "TargetGenes", "Regulon",
    "MotifID", "Annotation", "NES", "AUC", "Context", "Enrichment",
    "MotifSimilarityQvalue", "OrthologousIdentity", "RankAtMax", "MotifID",
    "FeatureID", "TF", "gene", "genes", "target", "targets"
  )]
  tokens <- tokens[!grepl("^[0-9.]+$", tokens)]
  tokens <- tokens[!grepl("^(cisbp|transfac|jaspar|hocomoco|taipale|swissregulon|metacluster)", tokens, ignore.case = TRUE)]
  tokens <- tokens[!grepl("^[A-Za-z]+[.][.][.][0-9]+$", tokens)]
  tokens
}

read_sox4_regulon_targets <- function(regulons_csv, sox4_regulon) {
  if (is.na(sox4_regulon) || !file.exists(regulons_csv) || file.info(regulons_csv)$size == 0) {
    return(tibble(target = character(), source = character()))
  }

  regulons <- suppressMessages(readr::read_csv(
    regulons_csv,
    show_col_types = FALSE,
    col_types = readr::cols(.default = readr::col_character())
  ))
  regulon_col <- intersect(
    c(
      "Regulons", "Regulon", "regulon", "regulons", "TF", "tf",
      "TranscriptionFactor", "transcription_factor", "module", "Module"
    ),
    colnames(regulons)
  )[1]
  target_col <- intersect(
    c(
      "TargetGenes", "Target_Genes", "target_genes", "targetGenes",
      "targets", "Targets", "genes", "Genes", "gene", "Gene",
      "EnrichedGenes", "enriched_genes"
    ),
    colnames(regulons)
  )[1]
  enrichment_cols <- grep("^Enrichment", colnames(regulons), value = TRUE)

  if (is.na(regulon_col)) {
    sox4_rows <- regulons %>%
      filter(if_any(everything(), ~ grepl("SOX4", .x, fixed = TRUE)))
  } else {
    sox4_rows <- regulons %>%
      filter(.data[[regulon_col]] == sox4_regulon | get_regulon_tf(.data[[regulon_col]]) == "SOX4")
  }

  if (nrow(sox4_rows) == 0L) {
    message(
      "[SKIP] No SOX4 rows found in ", regulons_csv,
      ". Columns: ", paste(colnames(regulons), collapse = ", ")
    )
    return(tibble(target = character(), source = character()))
  }

  if (is.na(target_col)) {
    if (length(enrichment_cols) > 0L) {
      message(
        "[INFO] Using pySCENIC enrichment columns as target-gene source in ",
        regulons_csv, ": ", paste(enrichment_cols, collapse = ", ")
      )
      target_values <- unlist(sox4_rows[, enrichment_cols, drop = FALSE], use.names = FALSE)
    } else {
      message(
        "[WARN] Could not identify target-gene column in ", regulons_csv,
        ". Columns: ", paste(colnames(regulons), collapse = ", "),
        ". Falling back to SOX4 row text parsing."
      )
      target_values <- apply(sox4_rows, 1, paste, collapse = " ")
    }
  } else {
    target_values <- sox4_rows[[target_col]]
  }

  tibble(
    target = extract_target_tokens(target_values),
    source = "regulons_csv"
  ) %>%
    filter(target != "SOX4")
}

read_sox4_adjacencies <- function(adj_tsv) {
  if (!file.exists(adj_tsv) || file.info(adj_tsv)$size == 0) {
    return(tibble(target = character(), importance = numeric(), source = character()))
  }

  adj <- suppressMessages(readr::read_tsv(adj_tsv, show_col_types = FALSE))
  tf_col <- intersect(c("TF", "tf", "regulator", "source"), colnames(adj))[1]
  target_col <- intersect(c("target", "Target", "gene", "Gene"), colnames(adj))[1]
  importance_col <- intersect(c("importance", "weight", "score"), colnames(adj))[1]
  if (is.na(tf_col) || is.na(target_col)) {
    message("[SKIP] Could not identify TF/target columns in ", adj_tsv)
    return(tibble(target = character(), importance = numeric(), source = character()))
  }
  if (is.na(importance_col)) {
    adj$importance <- 1
    importance_col <- "importance"
  }

  adj %>%
    filter(.data[[tf_col]] == "SOX4") %>%
    transmute(
      target = as.character(.data[[target_col]]),
      importance = as.numeric(.data[[importance_col]]),
      source = "adjacencies_tsv"
    ) %>%
    filter(!is.na(target), target != "SOX4")
}

get_expression_matrix <- function(seurat_obj, assay, layer, features, cells) {
  if (!assay %in% names(seurat_obj@assays)) {
    stop("Expression assay not found in Seurat object: ", assay, call. = FALSE)
  }
  expression_matrix <- tryCatch(
    GetAssayData(seurat_obj, assay = assay, layer = layer),
    error = function(e) GetAssayData(seurat_obj, assay = assay, slot = layer)
  )
  present_features <- intersect(features, rownames(expression_matrix))
  if (length(present_features) == 0L) {
    return(matrix(numeric(), nrow = 0L, ncol = length(cells), dimnames = list(NULL, cells)))
  }
  expression_matrix[present_features, cells, drop = FALSE]
}

safe_cor <- function(x, y) {
  if (length(x) < 3L || stats::sd(x) == 0 || stats::sd(y) == 0) {
    return(NA_real_)
  }
  suppressWarnings(stats::cor(x, y, method = "spearman", use = "complete.obs"))
}

dense_feature_vector <- function(expression_matrix, feature) {
  values <- as.numeric(as.matrix(expression_matrix[feature, , drop = FALSE]))
  names(values) <- colnames(expression_matrix)
  values
}

make_sox4_associations <- function(
  auc_matrix,
  seurat_obj,
  groups,
  regulons_csv,
  adj_tsv,
  assay,
  layer,
  top_targets
) {
  sox4_regulon <- choose_tf_regulon(auc_matrix, "SOX4")
  if (is.na(sox4_regulon)) {
    message("[SKIP] SOX4 regulon was not detected in AUCell matrix.")
    return(list(
      sox4_regulon = NA_character_,
      target_table = tibble(),
      expression_matrix = matrix(numeric(), nrow = 0L, ncol = ncol(auc_matrix))
    ))
  }

  regulon_targets <- read_sox4_regulon_targets(regulons_csv, sox4_regulon)
  adjacency_targets <- read_sox4_adjacencies(adj_tsv)
  target_table <- full_join(
    regulon_targets %>% mutate(in_regulon = TRUE),
    adjacency_targets %>% mutate(in_adjacencies = TRUE),
    by = "target",
    suffix = c("_regulon", "_adjacency")
  ) %>%
    mutate(
      in_regulon = ifelse(is.na(in_regulon), FALSE, in_regulon),
      in_adjacencies = ifelse(is.na(in_adjacencies), FALSE, in_adjacencies),
      importance = ifelse(is.na(importance), 0, importance)
    )

  if (nrow(target_table) == 0L) {
    message("[SKIP] No SOX4 target genes found in regulons/adjacencies outputs.")
    return(list(
      sox4_regulon = sox4_regulon,
      target_table = tibble(),
      expression_matrix = matrix(numeric(), nrow = 0L, ncol = ncol(auc_matrix))
    ))
  }

  expression_matrix <- get_expression_matrix(
    seurat_obj,
    assay,
    layer,
    target_table$target,
    colnames(auc_matrix)
  )
  present_targets <- rownames(expression_matrix)
  if (length(present_targets) == 0L) {
    message("[SKIP] SOX4 target genes were not present in the selected expression assay/layer.")
    return(list(
      sox4_regulon = sox4_regulon,
      target_table = target_table,
      expression_matrix = expression_matrix
    ))
  }

  sox4_auc <- as.numeric(auc_matrix[sox4_regulon, colnames(expression_matrix)])
  names(sox4_auc) <- colnames(expression_matrix)
  group_levels <- unique(groups[colnames(expression_matrix)])

  association_table <- bind_rows(lapply(present_targets, function(target) {
    expression <- dense_feature_vector(expression_matrix, target)
    bind_rows(lapply(group_levels, function(group_name) {
      cells <- names(groups)[groups == group_name]
      cells <- intersect(cells, colnames(expression_matrix))
      rho <- safe_cor(sox4_auc[cells], expression[cells])
      tibble(
        target = target,
        group = group_name,
        n_cells = length(cells),
        mean_expression = mean(expression[cells]),
        mean_sox4_auc = mean(sox4_auc[cells]),
        spearman_rho = rho,
        direction = case_when(
          is.na(rho) ~ "not_estimated",
          rho > 0 ~ "promotion_association",
          rho < 0 ~ "inhibition_association",
          TRUE ~ "neutral"
        )
      )
    }))
  })) %>%
    left_join(target_table, by = "target") %>%
    group_by(target) %>%
    mutate(
      max_abs_rho = max(abs(spearman_rho), na.rm = TRUE),
      max_mean_expression = max(mean_expression, na.rm = TRUE)
    ) %>%
    ungroup() %>%
    mutate(
      max_abs_rho = ifelse(is.infinite(max_abs_rho), NA_real_, max_abs_rho),
      plot_score = coalesce(max_abs_rho, 0) + log1p(importance)
    ) %>%
    arrange(desc(plot_score), desc(max_mean_expression), target)

  selected_targets <- association_table %>%
    distinct(target, plot_score, max_mean_expression) %>%
    slice_head(n = top_targets) %>%
    pull(target)

  list(
    sox4_regulon = sox4_regulon,
    target_table = association_table %>%
      mutate(selected_for_plot = target %in% selected_targets),
    expression_matrix = expression_matrix[selected_targets, , drop = FALSE]
  )
}

make_sox4_heatmap <- function(expression_matrix, groups, association_table) {
  if (nrow(expression_matrix) == 0L) {
    return(NULL)
  }
  group_levels <- unique(groups[colnames(expression_matrix)])
  mean_expression <- vapply(
    group_levels,
    function(group_name) {
      cells <- names(groups)[groups == group_name]
      cells <- intersect(cells, colnames(expression_matrix))
      Matrix::rowMeans(expression_matrix[, cells, drop = FALSE])
    },
    numeric(nrow(expression_matrix))
  )
  rownames(mean_expression) <- rownames(expression_matrix)
  colnames(mean_expression) <- group_levels

  row_sd <- apply(mean_expression, 1, sd)
  row_sd[is.na(row_sd) | row_sd == 0] <- 1
  z_matrix <- sweep(mean_expression, 1, rowMeans(mean_expression), "-")
  z_matrix <- sweep(z_matrix, 1, row_sd, "/")

  direction_by_target <- association_table %>%
    filter(target %in% rownames(z_matrix)) %>%
    group_by(target) %>%
    slice_max(order_by = abs(spearman_rho), n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    select(target, direction)
  row_direction <- direction_by_target$direction[match(rownames(z_matrix), direction_by_target$target)]
  row_direction[is.na(row_direction)] <- "not_estimated"

  row_ha <- rowAnnotation(
    association = row_direction,
    col = list(
      association = c(
        promotion_association = "#B2182B",
        inhibition_association = "#2166AC",
        neutral = "#BDBDBD",
        not_estimated = "#F0F0F0"
      )
    ),
    show_annotation_name = TRUE
  )

  Heatmap(
    z_matrix,
    name = "expr z-score",
    col = colorRamp2(c(-2, 0, 2), c("#2166AC", "#F7F7F7", "#B2182B")),
    cluster_rows = TRUE,
    cluster_columns = FALSE,
    show_row_names = TRUE,
    row_names_gp = gpar(fontsize = 7),
    column_names_gp = gpar(fontsize = 9),
    left_annotation = row_ha,
    row_title = "SOX4 target genes",
    column_title = "Mean target expression by T-cell state"
  )
}

# ---------------------------------------------------------------------------
# Main workflow
# ---------------------------------------------------------------------------

# 1. Resolve paths, validate inputs, and create output directories.
config <- parse_args(commandArgs(trailingOnly = TRUE))
setwd(config$project_dir)
source_plotting_shared()

check_file(config$auc_loom, "pySCENIC AUCell loom")
check_file(config$seurat_rds, "annotated Seurat RDS")
dir.create(config$figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(config$table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(config$result_rds), recursive = TRUE, showWarnings = FALSE)

# 2. Read AUCell activity and align it with Seurat metadata/UMAP cells.
message("Reading regulon AUC matrix: ", config$auc_loom)
auc_matrix <- read_regulon_auc(config$auc_loom)
seurat_obj <- readRDS(config$seurat_rds)
aligned <- align_auc_and_seurat(
  auc_matrix,
  seurat_obj,
  config$group_by,
  config$reduction
)
auc_matrix <- aligned$auc
seurat_obj <- aligned$seurat
groups <- aligned$groups
group_palette <- get_group_palette(config$group_by, groups)

message(
  "Aligned ", ncol(auc_matrix), " cells and ",
  nrow(auc_matrix), " regulons."
)

# 3. Rank regulons by between-group variation for grouped summaries.
summary_result <- summarise_regulons(auc_matrix, groups)
n_top <- min(config$n_top_regulons, nrow(auc_matrix))
top_regulons <- summary_result$top_regulons[seq_len(n_top)]

# 4. Write tidy summary statistics and the exact regulon order used in plots.
write_csv(
  summary_result$stats %>% arrange(desc(specificity), regulon, group),
  file.path(config$table_dir, "pyscenic_regulon_activity_by_group.csv")
)
write_lines(
  top_regulons,
  file.path(config$table_dir, "pyscenic_selected_regulons.txt")
)

# 5. Heatmap: relative group activity. Each row is standardized independently,
# so colors compare groups within a regulon, not absolute AUC across regulons.
heatmap_matrix <- summary_result$zscore[top_regulons, , drop = FALSE]
heatmap <- Heatmap(
  heatmap_matrix,
  name = "z-score",
  col = colorRamp2(c(-2, 0, 2), c("#2166AC", "#F7F7F7", "#B2182B")),
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  show_row_names = TRUE,
  row_names_gp = gpar(fontsize = 7),
  column_names_gp = gpar(fontsize = 9),
  border = TRUE
)
save_complex_heatmap(
  heatmap,
  file.path(config$figure_dir, paste0("pyscenic_regulon_activity_heatmap_by_", config$group_by)),
  width = 8,
  height = max(6, n_top * 0.22),
  overwrite = config$overwrite_plots
)

# 6. Dot plot: color shows absolute mean AUC; point size shows the percentage
# of cells above the regulon's global median AUC.
dotplot_data <- make_dotplot_data(auc_matrix, groups, top_regulons)
dotplot_data$regulon <- factor(dotplot_data$regulon, levels = rev(top_regulons))
dotplot_data$group <- factor(dotplot_data$group, levels = unique(groups))
dotplot <- ggplot(dotplot_data, aes(group, regulon)) +
  geom_point(aes(size = percent_active, color = mean_auc)) +
  scale_color_viridis_c(option = "plasma") +
  scale_size(range = c(1, 7)) +
  theme_bw() +
  theme(
    axis.title = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_line(linewidth = 0.2, color = "grey90")
  ) +
  labs(color = "Mean AUC", size = "% above median")
save_ggplot(
  dotplot,
  file.path(config$figure_dir, paste0("pyscenic_regulon_activity_dotplot_by_", config$group_by)),
  width = max(8, length(unique(groups)) * 0.8),
  height = max(7, n_top * 0.25),
  overwrite = config$overwrite_plots
)

# 7. Curated TF panel: preserve the supplied biological interpretation table,
# match TF symbols to pySCENIC regulon names, and draw focused activity plots.
curated_tf_table <- make_curated_tf_annotations() %>%
  match_curated_regulons(rownames(auc_matrix))
write_csv(
  curated_tf_table,
  file.path(config$table_dir, "pyscenic_curated_tf_interpretation.csv")
)

detected_curated <- curated_tf_table %>% filter(detected)
if (nrow(detected_curated) == 0L) {
  message("[SKIP] None of the curated TF regulons were detected in the AUCell loom.")
} else {
  curated_regulons <- detected_curated$matched_regulon
  curated_tf_labels <- ifelse(
    is.na(detected_curated$regulon_sign),
    detected_curated$tf,
    paste0(detected_curated$tf, " (", detected_curated$regulon_sign, ")")
  )

  curated_heatmap_matrix <- summary_result$zscore[curated_regulons, , drop = FALSE]
  rownames(curated_heatmap_matrix) <- curated_tf_labels
  if (config$group_by == "cell_type" && exists("cell_type_levels", inherits = TRUE)) {
    ordered_cell_types <- get("cell_type_levels", inherits = TRUE)
    ordered_columns <- ordered_cell_types[ordered_cell_types %in% colnames(curated_heatmap_matrix)]
    ordered_columns <- c(ordered_columns, setdiff(colnames(curated_heatmap_matrix), ordered_columns))
    curated_heatmap_matrix <- curated_heatmap_matrix[, ordered_columns, drop = FALSE]
  }
  if (!requireNamespace("SeuratExtend", quietly = TRUE)) {
    stop("SeuratExtend is required for the curated TF heatmap.", call. = FALSE)
  }
  curated_heatmap <- SeuratExtend::Heatmap(
    curated_heatmap_matrix,
    lab_fill = "zscore"
  )
  save_ggplot(
    curated_heatmap,
    file.path(config$figure_dir, "pyscenic_curated_tf_heatmap"),
    width = 9,
    height = max(5, nrow(detected_curated) * 0.34),
    overwrite = config$overwrite_plots
  )

}

curated_tf_expression_matrix <- get_expression_matrix(
  seurat_obj = seurat_obj,
  assay = config$assay,
  layer = config$expr_layer,
  features = curated_tf_table$tf,
  cells = colnames(auc_matrix)
)
curated_tf_expression_status <- curated_tf_table %>%
  transmute(
    tf,
    expression_detected = tf %in% rownames(curated_tf_expression_matrix)
  )
write_csv(
  curated_tf_expression_status,
  file.path(config$table_dir, "pyscenic_curated_tf_feature_status.csv")
)

# 8. Focused Vd2 effector-state comparison. This reuses the all-cell
# pySCENIC AUCell matrix already loaded above; it does not rerun pySCENIC on
# the three cell types. The subset is only for visualization/statistics.
vd2_effector_comparison <- make_vd2_effector_comparison(
  auc_matrix = auc_matrix,
  seurat_obj = seurat_obj,
  states = vd2_effector_compare_states,
  curated_tfs = vd2_effector_curated_tfs,
  top_n = config$n_top_regulons
)
if (is.null(vd2_effector_comparison)) {
  message("[SKIP] Vd2 effector pySCENIC comparison was not generated.")
} else {
  write_csv(
    vd2_effector_comparison$delta_table,
    file.path(config$table_dir, "pyscenic_vd2_effector_regulon_delta.csv")
  )

  vd2_effector_dotplot <- plot_vd2_effector_dotplot(vd2_effector_comparison$dotplot_data)
  save_ggplot(
    vd2_effector_dotplot,
    file.path(config$figure_dir, "pyscenic_vd2_effector_regulon_dotplot"),
    width = 8,
    height = max(6, length(unique(vd2_effector_comparison$dotplot_data$tf)) * 0.28),
    overwrite = config$overwrite_plots
  )

}

# 9. ZOL FOXP3+ Vd2 focused function-gene view. This combines the all-cell
# pySCENIC FOXP3 regulon AUC, if detected, with expression of regulatory,
# checkpoint, effector-context, and stress/fitness genes across all Vd2
# states. It is a visualization/statistics subset only.
zol_foxp3_function_result <- make_zol_foxp3_function_gene_analysis(
  auc_matrix = auc_matrix,
  seurat_obj = seurat_obj,
  assay = config$assay,
  layer = config$expr_layer,
  states = zol_foxp3_vd2_states,
  gene_sets = zol_foxp3_function_gene_sets
)
if (is.null(zol_foxp3_function_result)) {
  message("[SKIP] ZOL FOXP3+ Vd2 function-gene analysis was not generated.")
} else {
  if (nrow(zol_foxp3_function_result$expression_summary) > 0L) {
    write_csv(
      zol_foxp3_function_result$expression_summary,
      file.path(config$table_dir, "pyscenic_zol_foxp3_vd2_function_gene_summary.csv")
    )
  }
  if (nrow(zol_foxp3_function_result$delta_table) > 0L) {
    write_csv(
      zol_foxp3_function_result$delta_table,
      file.path(config$table_dir, "pyscenic_zol_foxp3_vd2_function_gene_delta.csv")
    )
  }
  if (is.null(zol_foxp3_function_result$heatmap)) {
    message("[SKIP] ZOL FOXP3+ function-gene heatmap has no expression data.")
  } else {
    save_complex_heatmap(
      zol_foxp3_function_result$heatmap,
      file.path(config$figure_dir, "pyscenic_zol_foxp3_vd2_function_gene_heatmap"),
      width = max(8, length(zol_foxp3_function_result$states) * 0.8),
      height = 8,
      overwrite = config$overwrite_plots
    )
  }
  if (is.null(zol_foxp3_function_result$dotplot)) {
    message("[SKIP] ZOL FOXP3+ function-gene dotplot has no expression data.")
  } else {
    save_ggplot(
      zol_foxp3_function_result$dotplot,
      file.path(config$figure_dir, "pyscenic_zol_foxp3_vd2_function_gene_dotplot"),
      width = max(8, length(zol_foxp3_function_result$states) * 0.8),
      height = 8,
      overwrite = config$overwrite_plots
    )
  }
}

# 10. SOX4 target relationship panel. Direction is inferred from Spearman
# correlation between SOX4 regulon AUC and target-gene expression within each
# T-cell state: positive rho = promotion association, negative rho =
# inhibition association. This is an association visualization, not proof of
# direct activation/repression.
sox4_result <- make_sox4_associations(
  auc_matrix = auc_matrix,
  seurat_obj = seurat_obj,
  groups = groups,
  regulons_csv = config$regulons_csv,
  adj_tsv = config$adj_tsv,
  assay = config$assay,
  layer = config$expr_layer,
  top_targets = config$sox4_top_targets
)

if (nrow(sox4_result$target_table) > 0L) {
  write_csv(
    sox4_result$target_table,
    file.path(config$table_dir, "pyscenic_sox4_target_associations.csv")
  )
}

if (!is.na(sox4_result$sox4_regulon) &&
  nrow(sox4_result$target_table) > 0L &&
  nrow(sox4_result$expression_matrix) > 0L) {
  sox4_heatmap <- make_sox4_heatmap(
    sox4_result$expression_matrix,
    groups,
    sox4_result$target_table
  )
  if (is.null(sox4_heatmap)) {
    message("[SKIP] SOX4 target expression heatmap has no available target genes.")
  } else {
    save_complex_heatmap(
      sox4_heatmap,
      file.path(config$figure_dir, "pyscenic_sox4_target_expression_heatmap"),
      width = 10,
      height = max(6, nrow(sox4_result$expression_matrix) * 0.24),
      overwrite = config$overwrite_plots
    )
  }

}

# 11. Save the aligned AUC matrix and derived summaries for custom follow-up
# plots without reopening the loom. This object can be large because it
# contains regulon activity for every aligned cell.
result <- list(
  auc_matrix = auc_matrix,
  group_by = config$group_by,
  groups = groups,
  mean_auc_by_group = summary_result$mean,
  zscore_by_group = summary_result$zscore,
  specificity = summary_result$specificity,
  selected_regulons = top_regulons,
  curated_tf_table = curated_tf_table,
  detected_curated_regulons = detected_curated,
  curated_tf_expression_status = curated_tf_expression_status,
  vd2_effector_comparison = vd2_effector_comparison,
  zol_foxp3_function_result = zol_foxp3_function_result,
  sox4_result = sox4_result,
  aligned_cells = colnames(auc_matrix)
)
if (config$overwrite || !file.exists(config$result_rds)) {
  saveRDS(result, config$result_rds)
} else {
  message("[SKIP] Existing result RDS: ", config$result_rds)
}

message("pySCENIC visualization complete.")
