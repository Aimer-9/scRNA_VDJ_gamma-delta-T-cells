# R version 4.5.2 (2025-10-31)
rm(list = ls())
# setwd("/path/to/project")
setwd("/path/to/project")
library(Seurat)
library(SeuratExtend)
library(genekitr) # gene annotation
library(gplots) # balloonplot
library(tidyverse)
library(patchwork) # for combining ggplot
options(
  tibble.width = Inf,
  print.width = Inf,
  max.print = 200
)
set.seed(1234)

metadata_file <- "config/samples.csv"
rds_dir <- "rds"
figure_dir <- file.path("figures", "3_cellAnnotation")
table_dir <- "table"

seurat_filtered_rds <- file.path(rds_dir, "all_seurat_2.rds")
vdj_stat_rds <- file.path(rds_dir, "VDJ_annotation_stat.rds")
seurat_celltype_rds <- file.path(rds_dir, "all_seurat_celltype.rds")
seurat_3_rds <- file.path(rds_dir, "all_seurat_3.rds")
seurat_4_rds <- file.path(rds_dir, "all_seurat_4.rds")
seurat_marker_file <- file.path(table_dir, "markers.csv")
seurat_3_marker_file <- file.path(table_dir, "all_seurat_3_seurat_markers.csv")
seurat_4_marker_file <- file.path(table_dir, "all_seurat_4_seurat_markers.csv")

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
# Blank-2 has less than 1K cells, exclude it.
excluded_sample_id <- c("Blank-2")
qc_features <- c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo")
cluster_dims <- 1:20
cluster_resolutions <- seq(0.1, 0.3, by = 0.1)
annotation_resolution <- "RNA_snn_res.0.3"
run_integration_check <- FALSE
integration_check_n_cells <- 2000
# Re-run filtering/table/RDS work from rds/all_seurat_2.rds. Leave FALSE to
# reuse the final annotated RDS when it already exists.
force_cell_annotation <- TRUE
# Re-draw figures even when their output files already exist. This is separate
# from force_cell_annotation so final plotting can be refreshed without
# rebuilding filtering outputs or marker tables.
force_cell_annotation_plot <- TRUE

dir.create(rds_dir, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, showWarnings = FALSE)

if (output_files_exist(seurat_celltype_rds) && !force_cell_annotation && !force_cell_annotation_plot) {
  message("[SKIP] Existing annotated Seurat RDS: ", seurat_celltype_rds)
  quit(save = "no")
} else if (output_files_exist(seurat_celltype_rds) && force_cell_annotation) {
  message("[FORCE] Rebuilding cell annotation outputs from current script settings: ", seurat_celltype_rds)
} else if (output_files_exist(seurat_celltype_rds) && force_cell_annotation_plot) {
  message("[PLOT] Re-drawing final cell annotation figures from existing annotated Seurat RDS: ", seurat_celltype_rds)
}

meta <- read.csv(metadata_file)
meta$group <- factor(meta$group, levels = group_levels)
meta$sample_name <- factor(
  paste0(meta$group, "_", meta$timepoint),
  levels = sample_name_levels
)
as_plain_character <- function(x) {
  if (is.null(x)) {
    return(character())
  }
  if (is.data.frame(x) || is.matrix(x)) {
    x <- x[, 1]
  }
  if (is.list(x)) {
    x <- unlist(x, use.names = FALSE)
  }
  if (isS4(x)) {
    s4_vector <- tryCatch(as.vector(x), error = function(e) NULL)
    if (!is.null(s4_vector) && !isS4(s4_vector)) {
      return(as_plain_character(s4_vector))
    }
    s4_names <- tryCatch(names(x), error = function(e) NULL)
    if (!is.null(s4_names)) {
      return(as_plain_character(s4_names))
    }
    stop("Cannot coerce S4 object of class ", paste(class(x), collapse = "/"), " to a character vector.", call. = FALSE)
  }
  as.character(x)
}

in_character_set <- function(x, table) {
  base::match(as_plain_character(x), as_plain_character(table), nomatch = 0L) > 0L
}

meta <- meta[!in_character_set(meta[["sample_id"]], excluded_sample_id), , drop = FALSE]

coerce_column_vector <- function(data, column, fallback = NULL) {
  if (!in_character_set(column, colnames(data))) {
    if (is.null(fallback)) {
      stop("Missing column: ", column, call. = FALSE)
    }
    return(as_plain_character(fallback))
  }
  value <- data[[column]]
  as_plain_character(value)
}

get_metadata_vector <- function(seurat_obj, column, fallback = NULL) {
  coerce_column_vector(seurat_obj@meta.data, column, fallback = fallback)
}

filter_included_samples <- function(seurat_obj, metadata) {
  orig_ident <- get_metadata_vector(seurat_obj, "orig.ident")
  sample_names <- get_metadata_vector(seurat_obj, "sample_name")
  included_sample_ids <- coerce_column_vector(metadata, "sample_id")
  included_sample_names <- coerce_column_vector(metadata, "sample_name")
  keep_cells <- colnames(seurat_obj)[
    in_character_set(orig_ident, included_sample_ids) &
      in_character_set(sample_names, included_sample_names)
  ]
  seurat_obj <- subset(
    seurat_obj,
    cells = keep_cells
  )
  seurat_obj@meta.data$group <- factor(get_metadata_vector(seurat_obj, "group"), levels = group_levels)
  seurat_obj@meta.data$sample_name <- factor(get_metadata_vector(seurat_obj, "sample_name"), levels = sample_name_levels)
  seurat_obj
}

add_v_gene_metadata <- function(seurat_obj, annotation) {
  v_gene_columns <- c("Vd1", "Vd2", "Vd3", "Vg2", "Vg3", "Vg4", "Vg5", "Vg8", "Vg9")
  missing_v_gene_columns <- setdiff(v_gene_columns, colnames(annotation))
  if (length(missing_v_gene_columns) > 0) {
    stop("Missing V gene annotation column(s): ", paste(missing_v_gene_columns, collapse = ", "), call. = FALSE)
  }
  seurat_barcodes <- get_metadata_vector(seurat_obj, "barcode", fallback = colnames(seurat_obj))
  annotation_barcodes <- coerce_column_vector(annotation, "barcode")
  vdj_match <- base::match(as_plain_character(seurat_barcodes), as_plain_character(annotation_barcodes))
  seurat_obj@meta.data[, v_gene_columns] <- annotation[vdj_match, v_gene_columns, drop = FALSE]
  seurat_obj
}

join_assay_layers_if_needed <- function(seurat_obj, assay = "RNA") {
  assay_names <- names(seurat_obj@assays)
  if (!in_character_set(assay, assay_names)) {
    return(seurat_obj)
  }
  assay_obj <- seurat_obj[[assay]]
  assay_layers <- tryCatch(Layers(assay_obj), error = function(e) character())
  layer_prefixes <- sub("\\..*$", "", assay_layers)
  has_split_layers <- length(unique(assay_layers)) > length(unique(layer_prefixes))
  has_multiple_same_type_layers <- any(table(layer_prefixes) > 1)

  if (inherits(assay_obj, "Assay5") && (has_split_layers || has_multiple_same_type_layers)) {
    message("[RUN] Joining Seurat v5 ", assay, " assay layers for layer-compatible downstream analysis.")
    seurat_obj <- JoinLayers(seurat_obj, assay = assay)
  }
  seurat_obj
}

run_resolution_clustering <- function(seurat_obj) {
  seurat_obj <- FindNeighbors(seurat_obj, dims = cluster_dims, reduction = "pca")
  for (resolution in cluster_resolutions) {
    seurat_obj <- FindClusters(seurat_obj, resolution = resolution)
  }
  seurat_obj
}

rerun_unintegrated_clustering <- function(seurat_obj) {
  seurat_obj <- join_assay_layers_if_needed(seurat_obj)
  seurat_obj <- FindVariableFeatures(seurat_obj)
  seurat_obj <- ScaleData(seurat_obj)
  seurat_obj <- RunPCA(seurat_obj, print = FALSE)
  ElbowPlot(seurat_obj)
  seurat_obj <- RunUMAP(seurat_obj,
    dims = cluster_dims,
    reduction = "pca",
    reduction.name = "umap.unintegrated"
  )
  run_resolution_clustering(seurat_obj)
}

find_markers_if_missing <- function(seurat_obj, ident_column, output_file, overwrite = FALSE) {
  if (output_files_exist(output_file) && !overwrite) {
    message("[SKIP] Existing marker table: ", output_file)
    return(read.csv(output_file))
  } else if (output_files_exist(output_file) && overwrite) {
    message("[FORCE] Rebuilding marker table: ", output_file)
  }

  seurat_obj <- join_assay_layers_if_needed(seurat_obj)
  Idents(seurat_obj) <- ident_column
  markers <- FindAllMarkers(seurat_obj)
  geneinfo <- genInfo(markers$gene) %>%
    unique()
  markers <- merge(markers, geneinfo,
    by.x = "gene",
    by.y = "symbol", all.x = TRUE
  )
  markers <- markers %>%
    mutate(fold.pct = pct.1 / pct.2, .before = 2) %>%
    filter(abs(avg_log2FC) > 1 & p_val_adj < 0.05)
  write.csv(markers, row.names = FALSE, file = output_file)
  markers
}

save_group_sample_umaps <- function(seurat_obj, prefix) {
  umap_sample <- make_umap_plot(seurat_obj, "sample_name", cols = color_sample_name)
  umap_group <- make_umap_plot(seurat_obj, "group", cols = color_group)
  umap_group_sample <- umap_group + umap_sample +
    plot_annotation(title = "UMAP of all cells by group and sample") &
    theme(plot.title = element_text(hjust = 0.5, size = 20, face = "bold"))
  save_plot(umap_group_sample, paste0(prefix, "_umap_group_sample"), 12, 6, overwrite = force_cell_annotation_plot)
}

save_cluster_umaps <- function(seurat_obj, prefix) {
  umap_cluster_res_0_1 <- make_umap_plot(seurat_obj, "RNA_snn_res.0.1")
  umap_cluster_res_0_2 <- make_umap_plot(seurat_obj, "RNA_snn_res.0.2")
  umap_cluster_res_0_3 <- make_umap_plot(seurat_obj, "RNA_snn_res.0.3")
  save_plot(umap_cluster_res_0_1, paste0(prefix, "_umap_cluster_res_0_1"), 8, 6, overwrite = force_cell_annotation_plot)
  save_plot(umap_cluster_res_0_2, paste0(prefix, "_umap_cluster_res_0_2"), 8, 6, overwrite = force_cell_annotation_plot)
  save_plot(umap_cluster_res_0_3, paste0(prefix, "_umap_cluster_res_0_3"), 8, 6, overwrite = force_cell_annotation_plot)
}

annotate_cell_types <- function(seurat_obj) {
  annotation_clusters <- get_metadata_vector(seurat_obj, annotation_resolution)
  seurat_obj@meta.data <- seurat_obj@meta.data %>%
    mutate(cell_type = case_when(
      in_character_set(annotation_clusters, c("0")) ~ "Effector Memory Vd2",
      in_character_set(annotation_clusters, c("1")) ~ "ZOL Effector Vd2",
      in_character_set(annotation_clusters, c("2", "6", "7")) ~ "Effector Vd1",
      in_character_set(annotation_clusters, c("3", "4", "9", "11")) ~ "PAN Effector Vd2",
      in_character_set(annotation_clusters, c("5")) ~ "Pre-activated Vd1",
      in_character_set(annotation_clusters, c("8")) ~ "Naive Vd1",
      in_character_set(annotation_clusters, c("10")) ~ "Pre-activated Vd2",
      in_character_set(annotation_clusters, c("12")) ~ "ZOL FOXP3+ Vd2"
    ))
  seurat_obj@meta.data$cell_type <- factor(get_metadata_vector(seurat_obj, "cell_type"), levels = cell_type_levels)
  Idents(seurat_obj) <- "cell_type"
  seurat_obj
}

grouped_features <- list(
  "Naive" = c("TCF7", "LEF1", "CCR7", "SELL", "IL7R"),
  "Vd2" = c("TNF", "CD40LG"),
  "Vd1" = c("MAL", "FCER1G"),
  "Proliferating" = c("MKI67", "TOP2A"),
  "Effector" = c("IL2RA", "GZMA", "GZMB", "GZMH", "IFNG"),
  "Exhaustion" = c("PDCD1", "CTLA4", "LAG3", "HAVCR2", "TIGIT"),
  "Regulatory" = c("FOXP3", "IL2RA", "IL10")
)

cd4_cd8_t_cell_markers <- c(
  "CD3D", "CD3E", "CD3G", "TRAC",
  "CD4", "IL7R", "CCR7", "LTB",
  "CD8A", "CD8B", "NKG7", "GNLY",
  "GZMA", "GZMB", "GZMH", "PRF1"
)

b_cell_contamination_markers <- c(
  "MS4A1", "CD79A", "CD79B", "CD74",
  "HLA-DRA", "BANK1", "CD37", "MZB1"
)

plasma_cell_contamination_markers <- c(
  "JCHAIN", "MZB1", "XBP1", "SDC1", "PRDM1",
  "DERL3", "HERPUD1", "FKBP11", "TNFRSF17"
)

myeloid_contamination_marker_sets <- list(
  "Monocyte markers" = c("LYZ", "LST1", "FCN1", "S100A8", "S100A9"),
  "DC markers" = c("FCER1A", "CLEC10A", "CST3", "LILRA4"),
  "Macrophage markers" = c("CD68", "CD163", "MSR1", "C1QA"),
  "Neutrophil markers" = c("FCGR3B", "CSF3R", "CXCR2", "MPO")
)

all_seurat_3_cluster10_contamination_marker_sets <- list(
  "Mast basophil eosinophil markers" = c(
    "TPSAB1", "TPSB2", "TPSD1", "CPA3", "MS4A2",
    "CTSG", "RNASE3", "RNASE2", "CLC", "PRG2", "GATA2"
  ),
  "Myeloid DC macrophage markers" = c(
    "APOE", "FCER1A", "LYZ", "CST3", "FCER1G", "IRF8",
    "MNDA", "TMEM176A", "TMEM176B", "MS4A3"
  ),
  "B plasma markers" = c(
    "JCHAIN", "BANK1", "MZB1", "XBP1",
    "CD79A", "CD79B", "MS4A1"
  )
)

save_marker_feature_plot <- function(seurat_obj, features, title, filename, width, height, ncol = 4) {
  # These panels document manual filtering decisions. Missing genes are common
  # across references, so filter the feature list before calling SeuratExtend's
  # DimPlot2 feature mode.
  seurat_obj <- join_assay_layers_if_needed(seurat_obj)
  present_features <- unique(features[in_character_set(features, rownames(seurat_obj))])
  missing_features <- setdiff(features, present_features)
  if (length(missing_features) > 0) {
    message("[SKIP] Missing marker feature(s) for ", title, ": ", paste(missing_features, collapse = ", "))
  }
  if (length(present_features) == 0) {
    message("[SKIP] No marker features available for ", title)
    return(invisible(FALSE))
  }

  marker_plot <- SeuratExtend::DimPlot2(
    seurat_obj,
    features = present_features,
    reduction = "umap.unintegrated",
    theme = NoAxes(),
    ncol = ncol
  )
  marker_plot <- add_fixed_umap_coordinates(marker_plot) +
    plot_annotation(title = title)
  save_plot(marker_plot, filename, width, height, overwrite = force_cell_annotation_plot)
  invisible(TRUE)
}

save_cluster_removal_plot <- function(seurat_obj, cluster_column, removed_clusters, title, filename) {
  # Highlight the exact cluster IDs being removed so the marker panels can be
  # interpreted together with cluster location on UMAP.
  cluster_values <- get_metadata_vector(seurat_obj, cluster_column)
  removed_label <- paste0("remove_", paste(removed_clusters, collapse = "_"))
  seurat_obj@meta.data$cluster_filter_status <- factor(
    ifelse(
      in_character_set(cluster_values, removed_clusters),
      removed_label,
      "kept"
    ),
    levels = c("kept", removed_label)
  )
  status_colors <- setNames(
    c("grey80", "#D62728"),
    c("kept", removed_label)
  )
  removal_plot <- DimPlot2(
    seurat_obj,
    group.by = "cluster_filter_status",
    reduction = "umap.unintegrated",
    theme = NoAxes(),
    label = FALSE,
    cols = status_colors
  )
  removal_plot <- add_fixed_umap_coordinates(removal_plot) +
    labs(title = title, color = "Filter status")
  save_plot(removal_plot, filename, 8, 6, overwrite = force_cell_annotation_plot)
  invisible(TRUE)
}

ensure_unintegrated_plot_inputs <- function(seurat_obj, required_cluster_column = NULL, label = "Seurat object") {
  reduction_names <- names(seurat_obj@reductions)
  if (is.null(reduction_names)) {
    reduction_names <- character()
  }
  required_cluster_columns <- unique(c(
    paste0("RNA_snn_res.", format(cluster_resolutions, trim = TRUE, scientific = FALSE)),
    required_cluster_column
  ))
  required_cluster_columns <- required_cluster_columns[!is.na(required_cluster_columns)]
  missing_umap <- !in_character_set("umap.unintegrated", reduction_names)
  missing_cluster <- !all(in_character_set(required_cluster_columns, colnames(seurat_obj@meta.data)))
  if (missing_umap || missing_cluster) {
    message(
      "[PLOT] Recomputing unintegrated PCA/UMAP/clusters for ", label,
      " because required plot inputs are missing."
    )
    seurat_obj <- rerun_unintegrated_clustering(seurat_obj)
  }
  seurat_obj
}

save_all_seurat_2_filter_diagnostics <- function(seurat_obj) {
  seurat_obj <- ensure_unintegrated_plot_inputs(
    seurat_obj,
    required_cluster_column = "RNA_snn_res.0.3",
    label = "all_seurat_2"
  )
  save_group_sample_umaps(seurat_obj, "all_seurat_2")
  save_cluster_umaps(seurat_obj, "all_seurat_2")

  save_cluster_removal_plot(
    seurat_obj,
    "RNA_snn_res.0.3",
    c("9", "11"),
    "all_seurat_2 clusters 9/11 selected for B-cell mixed-feature removal",
    "all_seurat_2_removed_cluster_highlight"
  )

  save_marker_feature_plot(
    seurat_obj,
    b_cell_contamination_markers,
    "B-cell markers before removing all_seurat_2 clusters 9/11",
    "all_seurat_2_b_cell_marker_diagnostic",
    14,
    8,
    ncol = 4
  )
  save_marker_feature_plot(
    seurat_obj,
    plasma_cell_contamination_markers,
    "Plasma-cell markers before removing all_seurat_2 clusters 9/11",
    "all_seurat_2_plasma_marker_diagnostic",
    12,
    8,
    ncol = 3
  )
  invisible(seurat_obj)
}

save_all_seurat_3_filter_diagnostics <- function(seurat_obj) {
  seurat_obj <- ensure_unintegrated_plot_inputs(
    seurat_obj,
    required_cluster_column = "RNA_snn_res.0.2",
    label = "all_seurat_3"
  )
  save_group_sample_umaps(seurat_obj, "all_seurat_3")
  save_cluster_umaps(seurat_obj, "all_seurat_3")

  save_cluster_removal_plot(
    seurat_obj,
    "RNA_snn_res.0.2",
    c("10"),
    "all_seurat_3 cluster 10 selected for myeloid mixed-feature removal",
    "all_seurat_3_removed_cluster_highlight"
  )

  purrr::iwalk(myeloid_contamination_marker_sets, function(features, marker_group) {
    save_marker_feature_plot(
      seurat_obj,
      features,
      paste(marker_group, "before removing all_seurat_3 cluster 10"),
      paste0("all_seurat_3_", str_replace_all(str_to_lower(marker_group), "[^a-z0-9]+", "_"), "_diagnostic"),
      10,
      6,
      ncol = 3
    )
  })
  purrr::iwalk(all_seurat_3_cluster10_contamination_marker_sets, function(features, marker_group) {
    save_marker_feature_plot(
      seurat_obj,
      features,
      paste(marker_group, "before removing all_seurat_3 cluster 10"),
      paste0("all_seurat_3_cluster10_", str_replace_all(str_to_lower(marker_group), "[^a-z0-9]+", "_"), "_diagnostic"),
      14,
      8,
      ncol = 4
    )
  })
  invisible(seurat_obj)
}

save_final_cell_annotation_outputs <- function(seurat_obj, annotation = NULL, save_final_rds = TRUE) {
  # Final plots can be refreshed from the annotated RDS without repeating the
  # filtering steps. Only add metadata that is missing from the loaded object.
  v_gene_plot_columns <- c("Vd1", "Vd2", "Vg9", "Vg4")
  if (!all(in_character_set(v_gene_plot_columns, colnames(seurat_obj@meta.data)))) {
    if (is.null(annotation)) {
      message("[SKIP] V gene UMAP needs V gene metadata but no annotation table was supplied.")
    } else {
      seurat_obj <- add_v_gene_metadata(seurat_obj, annotation)
    }
  }
  if (all(in_character_set(v_gene_plot_columns, colnames(seurat_obj@meta.data)))) {
    umap_vd12_vg49 <- make_umap_plot(seurat_obj, v_gene_plot_columns, ncol = 2)
    save_plot(umap_vd12_vg49, "all_seurat_4_umap_vd12_vg49", 12, 6, overwrite = force_cell_annotation_plot)
  }

  if (!in_character_set("cell_type", colnames(seurat_obj@meta.data))) {
    seurat_obj <- annotate_cell_types(seurat_obj)
  } else {
    seurat_obj@meta.data$cell_type <- factor(get_metadata_vector(seurat_obj, "cell_type"), levels = cell_type_levels)
  }

  if (!in_character_set("Phase", colnames(seurat_obj@meta.data))) {
    seurat_obj <- CellCycleScoring(seurat_obj,
      s.features = cc.genes$s.genes,
      g2m.features = cc.genes$g2m.genes,
      set.ident = TRUE
    )
  }
  seurat_obj@meta.data$Phase <- factor(get_metadata_vector(seurat_obj, "Phase"), levels = names(color_cellcycle))
  umap_cellcycle <- make_umap_plot(seurat_obj, "Phase", cols = color_cellcycle)
  save_plot(umap_cellcycle, "all_seurat_4_umap_cellcycle", 10, 10, overwrite = force_cell_annotation_plot)

  umap_celltype <- make_umap_plot(seurat_obj, "cell_type", cols = color_celltype)
  save_plot(umap_celltype, "all_seurat_4_umap_celltype", 10, 10, overwrite = force_cell_annotation_plot)

  save_marker_feature_plot(
    seurat_obj,
    cd4_cd8_t_cell_markers,
    "CD4/CD8 T-cell marker expression on annotated UMAP",
    "all_seurat_4_cd4_cd8_t_cell_marker_umap",
    16,
    16,
    ncol = 4
  )

  celltype_distr <- ClusterDistrBar(
    origin = get_metadata_vector(seurat_obj, "sample_name"),
    cluster = get_metadata_vector(seurat_obj, "cell_type"),
    cols = color_celltype
  )
  save_plot(celltype_distr, "celltype_distr", 10, 5, overwrite = force_cell_annotation_plot)

  qc_violin_plot <- VlnPlot2(seurat_obj,
    group.by = "cell_type",
    features = qc_features,
    nrow = 2,
    ncol = 2,
    pt = FALSE,
    cols = color_celltype
  ) +
    labs(title = "QC Metrics Violin Plot by Cell Type")
  save_plot(qc_violin_plot, "qc_violin_plot_celltype", 12, 8, overwrite = force_cell_annotation_plot)


  dotplot_celltype <- DotPlot2(seurat_obj,
    group.by = "cell_type",
    features = grouped_features,
    show_grid = FALSE, flip = TRUE
  )
  save_plot(dotplot_celltype, "dotplot_celltype", 12, 4, overwrite = force_cell_annotation_plot)

  if (save_final_rds) {
    save_rds_if_missing(
      seurat_obj,
      seurat_celltype_rds,
      "annotated Seurat RDS",
      overwrite = force_cell_annotation
    )
  }
  invisible(seurat_obj)
}

if (output_files_exist(seurat_celltype_rds) && !force_cell_annotation && force_cell_annotation_plot) {
  annotation_stat <- if (output_files_exist(vdj_stat_rds)) {
    readRDS(vdj_stat_rds)
  } else {
    message("[SKIP] Missing VDJ annotation table for optional V gene plot metadata: ", vdj_stat_rds)
    NULL
  }
  all_seurat_celltype <- readRDS(seurat_celltype_rds)
  save_final_cell_annotation_outputs(
    all_seurat_celltype,
    annotation = annotation_stat,
    save_final_rds = FALSE
  )
  quit(save = "no")
}

all_seurat_2 <- readRDS(seurat_filtered_rds)
annotation_stat <- readRDS(vdj_stat_rds)
all_seurat_2 <- filter_included_samples(all_seurat_2, meta)

if (run_integration_check) {
  set.seed(1)
  integration_check_barcodes <- get_metadata_vector(all_seurat_2, "barcode", fallback = colnames(all_seurat_2))
  selected_integration_barcodes <- sample(
    integration_check_barcodes,
    size = min(integration_check_n_cells, length(integration_check_barcodes))
  )
  all_seurat_2_sub <- subset(
    all_seurat_2,
    cells = colnames(all_seurat_2)[in_character_set(integration_check_barcodes, selected_integration_barcodes)]
  )
  all_seurat_2_sub <- SCTransform(all_seurat_2_sub, vars.to.regress = "percent.mt")
  all_seurat_2_sub <- RunPCA(all_seurat_2_sub, assay = "SCT")
  all_seurat_2_sub <- IntegrateLayers(
    all_seurat_2_sub,
    method = HarmonyIntegration,
    new.reduction = "harmony",
    orig.reduction = "pca",
    assay = "SCT"
  )
  all_seurat_2_sub <- IntegrateLayers(
    all_seurat_2_sub,
    method = CCAIntegration,
    new.reduction = "CCA",
    orig.reduction = "pca",
    assay = "SCT"
  )
}

# cluster
all_seurat_2 <- join_assay_layers_if_needed(all_seurat_2)
all_seurat_2 <- run_resolution_clustering(all_seurat_2)

# ...Process find Marker...
markers <- find_markers_if_missing(
  all_seurat_2,
  "RNA_snn_res.0.3",
  seurat_marker_file,
  overwrite = force_cell_annotation
)

# find mixed cells likely from B cells, remove them
# cluster 9 and 11 in 0.3 resolution
all_seurat_2 <- save_all_seurat_2_filter_diagnostics(all_seurat_2)
cluster_0_3 <- get_metadata_vector(all_seurat_2, "RNA_snn_res.0.3")
cells_to_remove <- colnames(all_seurat_2)[in_character_set(cluster_0_3, c("9", "11"))]
all_seurat_3 <- subset(all_seurat_2, cells = cells_to_remove, invert = TRUE)

# rerun all downstream manifold steps
all_seurat_3 <- rerun_unintegrated_clustering(all_seurat_3)

# ...Process find Marker...
markers <- find_markers_if_missing(
  all_seurat_3,
  "RNA_snn_res.0.2",
  seurat_3_marker_file,
  overwrite = force_cell_annotation
)

# cluster 10 and 11 separated into 2 parts
all_seurat_3 <- save_all_seurat_3_filter_diagnostics(all_seurat_3)

# filter
cluster_0_2 <- get_metadata_vector(all_seurat_3, "RNA_snn_res.0.2")
cells_to_remove_2 <- colnames(all_seurat_3)[in_character_set(cluster_0_2, c("10"))]
all_seurat_4 <- subset(all_seurat_3, cells = cells_to_remove_2, invert = TRUE)

# do downstream steps again
all_seurat_4 <- rerun_unintegrated_clustering(all_seurat_4)

save_group_sample_umaps(all_seurat_4, "all_seurat_4")
# samples in each group concentrate well, seems no need for integration

# visualize clustering results
save_cluster_umaps(all_seurat_4, "all_seurat_4")

markers <- find_markers_if_missing(
  all_seurat_4,
  annotation_resolution,
  seurat_4_marker_file,
  overwrite = force_cell_annotation
)

# Visualize V gene distribution, annotate cell types, draw final plots, and save
# the final annotated object. Figure overwrites use force_cell_annotation_plot;
# the final RDS still uses force_cell_annotation.
all_seurat_celltype <- save_final_cell_annotation_outputs(
  all_seurat_4,
  annotation = annotation_stat,
  save_final_rds = TRUE
)
