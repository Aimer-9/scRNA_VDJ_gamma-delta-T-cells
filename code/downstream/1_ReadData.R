# R version 4.5.2 (2025-10-31)
rm(list = ls())
setwd("/path/to/project")
library(Seurat)
library(SeuratExtend)
library(readxl)
library(tidyverse)
library(patchwork)
library(ggpubr)
library(gplots)
library(SoupX)
library(Matrix)
options(
  tibble.width = Inf,
  print.width = Inf,
  max.print = 200
)

cellranger_dir <- "/path/to/cellranger/output"
metadata_file <- "config/samples.csv"
rds_dir <- "rds"
figure_dir <- file.path("figures", "1_ReadData")
soupx_dir <- "soupx"

vdj_annotation_filtered_rds <- file.path(rds_dir, "All_VDJ_annotation_filtered.rds")
vdj_annotation_stat_rds <- file.path(rds_dir, "VDJ_annotation_stat.rds")
seurat_before_filtering_rds <- file.path(rds_dir, "all_seurat_before_filtering.rds")

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
# Blank-2 has less than 1K cells, exclude it
excluded_sample_id <- c("Blank-2")
v_gene_columns <- c(
  "Vd1" = "TRDV1",
  "Vd2" = "TRDV2",
  "Vd3" = "TRDV3",
  "Vg2" = "TRGV2",
  "Vg3" = "TRGV3",
  "Vg4" = "TRGV4",
  "Vg5" = "TRGV5",
  "Vg8" = "TRGV8",
  "Vg9" = "TRGV9"
)
qc_features <- c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo")
min_cells_per_gene <- 3
soupx_cluster_dims <- 1:20
soupx_cluster_resolution <- 0.5
force_read_data <- FALSE
force_read_data_plot <- TRUE

dir(cellranger_dir)
dir.create(rds_dir, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(soupx_dir, showWarnings = FALSE)

if (output_files_exist(c(vdj_annotation_filtered_rds, vdj_annotation_stat_rds, seurat_before_filtering_rds)) && !force_read_data && !force_read_data_plot) {
  message("[SKIP] Existing 1_ReadData checkpoint RDS files: ", paste(
    c(vdj_annotation_filtered_rds, vdj_annotation_stat_rds, seurat_before_filtering_rds),
    collapse = ", "
  ))
  quit(save = "no")
} else if (output_files_exist(c(vdj_annotation_filtered_rds, vdj_annotation_stat_rds, seurat_before_filtering_rds)) && force_read_data) {
  message("[FORCE] Rebuilding 1_ReadData outputs from current script settings.")
} else if (output_files_exist(c(vdj_annotation_filtered_rds, vdj_annotation_stat_rds, seurat_before_filtering_rds)) && force_read_data_plot) {
  message("[PLOT] Re-drawing 1_ReadData figures without forcing checkpoint RDS overwrites.")
}

read_vdj_annotation <- function(filepath) {
  annotation <- read.csv(filepath)
  annotation$sample <- str_split(annotation$sample, "_", simplify = TRUE)[, 1]
  annotation$barcode <- paste0(annotation$sample, "_", annotation$barcode)
  annotation
}

read_10x_counts <- function(matrix_dir) {
  counts <- Read10X(matrix_dir)
  if (is.list(counts)) {
    if ("Gene Expression" %in% names(counts)) {
      counts <- counts[["Gene Expression"]]
    } else {
      counts <- counts[[1]]
    }
  }
  counts
}

write_sparse_matrix_parts <- function(counts, prefix, out_dir, overwrite = FALSE) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  output_paths <- file.path(out_dir, paste0(prefix, c(".mtx", "_genes.csv", "_barcodes.csv")))
  if (output_files_exist(output_paths) && !overwrite) {
    message("[SKIP] Existing sparse matrix parts: ", paste(output_paths, collapse = ", "))
    return(invisible(NULL))
  } else if (output_files_exist(output_paths) && overwrite) {
    message("[FORCE] Rebuilding sparse matrix parts: ", paste(output_paths, collapse = ", "))
  }
  writeMM(counts, output_paths[1])
  write.csv(
    data.frame(gene = rownames(counts)),
    output_paths[2],
    row.names = FALSE
  )
  write.csv(
    data.frame(barcode = colnames(counts)),
    output_paths[3],
    row.names = FALSE
  )
}

make_soupx_plots <- function(soup_channel, sample_id, out_dir) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  contamination <- soup_channel$metaData
  contamination$barcode <- rownames(contamination)

  p <- ggplot(contamination, aes(x = nUMIs, y = rho)) +
    geom_point(size = 0.35, alpha = 0.35) +
    scale_x_log10() +
    labs(
      title = paste(sample_id, "SoupX contamination estimate"),
      x = "UMIs per cell",
      y = "Estimated contamination fraction"
    ) +
    theme_classic()

  output_paths <- file.path(out_dir, paste0(sample_id, "_soupx_contamination", c(".png", ".pdf")))
  if (output_files_exist(output_paths) && !force_read_data_plot) {
    message("[SKIP] Existing SoupX plot: ", paste(output_paths, collapse = ", "))
    return(invisible(NULL))
  } else if (output_files_exist(output_paths) && force_read_data_plot) {
    message("[FORCE] Rebuilding SoupX plot: ", paste(output_paths, collapse = ", "))
  }
  ggsave(output_paths[2],
    p,
    width = 6,
    height = 4
  )
  ggsave(output_paths[1],
    p,
    width = 6,
    height = 4
  )
}

run_soupx <- function(sample_id, raw_dir, filtered_dir, out_dir) {
  raw_counts <- read_10x_counts(raw_dir)
  filtered_counts <- read_10x_counts(filtered_dir)

  soup_channel <- SoupChannel(tod = raw_counts, toc = filtered_counts)

  clustering_obj <- CreateSeuratObject(
    counts = filtered_counts,
    project = sample_id,
    min.cells = min_cells_per_gene
  )
  clustering_obj[["percent.mt"]] <- PercentageFeatureSet(clustering_obj, pattern = "^MT-|^mt-")
  clustering_obj <- NormalizeData(clustering_obj, verbose = FALSE)
  clustering_obj <- FindVariableFeatures(clustering_obj, verbose = FALSE)
  clustering_obj <- ScaleData(clustering_obj, verbose = FALSE)
  clustering_obj <- RunPCA(clustering_obj, verbose = FALSE)
  clustering_obj <- FindNeighbors(clustering_obj, dims = soupx_cluster_dims, verbose = FALSE)
  clustering_obj <- FindClusters(clustering_obj,
    resolution = soupx_cluster_resolution,
    verbose = FALSE
  )
  clustering_obj <- RunUMAP(clustering_obj, dims = soupx_cluster_dims, verbose = FALSE)

  soup_channel <- setClusters(
    soup_channel,
    setNames(clustering_obj$seurat_clusters, colnames(clustering_obj))
  )
  soup_channel <- setDR(soup_channel, Embeddings(clustering_obj, "umap"))
  soup_channel <- autoEstCont(soup_channel, doPlot = FALSE)

  adjusted_counts <- adjustCounts(soup_channel, roundToInt = TRUE)

  write_sparse_matrix_parts(
    adjusted_counts,
    paste0(sample_id, "_soupx_adjusted_counts"),
    out_dir,
    overwrite = force_read_data
  )
  soupx_metadata_csv <- file.path(out_dir, paste0(sample_id, "_soupx_metadata.csv"))
  if (output_files_exist(soupx_metadata_csv) && !force_read_data) {
    message("[SKIP] Existing SoupX metadata CSV: ", soupx_metadata_csv)
  } else if (output_files_exist(soupx_metadata_csv) && force_read_data) {
    message("[FORCE] Rebuilding SoupX metadata CSV: ", soupx_metadata_csv)
    write.csv(soup_channel$metaData, soupx_metadata_csv)
  } else {
    write.csv(soup_channel$metaData, soupx_metadata_csv)
  }
  make_soupx_plots(soup_channel, sample_id, out_dir)

  list(
    adjusted_counts = adjusted_counts,
    soup_channel = soup_channel
  )
}

add_v_gene_flags <- function(annotation) {
  annotation %>%
    group_by(barcode) %>%
    mutate(
      Vd1 = ifelse(v_gene_columns["Vd1"] %in% v_gene, "true", "false"),
      Vd2 = ifelse(v_gene_columns["Vd2"] %in% v_gene, "true", "false"),
      Vd3 = ifelse(v_gene_columns["Vd3"] %in% v_gene, "true", "false"),
      Vg2 = ifelse(v_gene_columns["Vg2"] %in% v_gene, "true", "false"),
      Vg3 = ifelse(v_gene_columns["Vg3"] %in% v_gene, "true", "false"),
      Vg4 = ifelse(v_gene_columns["Vg4"] %in% v_gene, "true", "false"),
      Vg5 = ifelse(v_gene_columns["Vg5"] %in% v_gene, "true", "false"),
      Vg8 = ifelse(v_gene_columns["Vg8"] %in% v_gene, "true", "false"),
      Vg9 = ifelse(v_gene_columns["Vg9"] %in% v_gene, "true", "false")
    ) %>%
    ungroup()
}

merge_seurat_samples <- function(seurat_list, sample_ids) {
  if (length(seurat_list) == 1) {
    return(RenameCells(seurat_list[[1]], add.cell.id = sample_ids[1]))
  }

  merge(seurat_list[[1]],
    seurat_list[2:length(seurat_list)],
    add.cell.ids = sample_ids
  )
}

# rawdata and metadata PATH
sample_list <- dir(cellranger_dir) %>%
  str_extract("(.*)(_add_enrichment_primers)", group = 1) %>%
  na.omit()
# exp file paths
raw_matrix_filepath <-
  paste0(
    cellranger_dir, "/",
    sample_list,
    "_add_enrichment_primers/outs/multi/count/raw_feature_bc_matrix/"
  )
filtered_matrix_filepath <-
  paste0(
    cellranger_dir, "/",
    sample_list,
    "_add_enrichment_primers/outs/per_sample_outs/",
    sample_list, "_add_enrichment_primers/count/sample_filtered_feature_bc_matrix/"
  )
# vdj annotation file paths
data_annotation_filepath <-
  paste0(
    cellranger_dir, "/",
    sample_list,
    "_add_enrichment_primers/outs/multi/vdj_t_gd/all_contig_annotations.csv"
  )
# metadata file path
meta <- read.csv(metadata_file)
meta$group <-
  factor(meta$group,
    levels = group_levels
  )
meta$sample_name <-
  factor(paste0(meta$group, "_", meta$timepoint),
    levels = sample_name_levels
  )
meta <- meta %>%
  filter(!sample_id %in% excluded_sample_id)
keep_samples <- sample_list %in% meta$sample_id
sample_list <- sample_list[keep_samples]
raw_matrix_filepath <- raw_matrix_filepath[keep_samples]
filtered_matrix_filepath <- filtered_matrix_filepath[keep_samples]
data_annotation_filepath <- data_annotation_filepath[keep_samples]
# ...Process vdj annotation data...
# combine all annotation data into one dataframe
all_annotation <-
  lapply(seq_along(sample_list), function(i) {
    read_vdj_annotation(data_annotation_filepath[i])
  }) %>%
  Reduce(bind_rows, .)
# merge with metadata
# sample is the name in raw data when sequencing
# sample_name is the proper sample name, with group info and batch info
# e.g. AB3_2 means group AB3, batch 2
all_annotation <- merge(all_annotation, meta, by.x = "sample", by.y = "sample_id")
# stat Vδ and Vγ gene usage, create new columns for each V gene of interest
# because one cell can have multiple contigs, we need to group by barcode
# and check if the V gene of interest is present in any of the contigs
all_annotation <- add_v_gene_flags(all_annotation)
# plot productive vs non-productive contig number per group
# group level productive vs non-productive contig barplot
all_annotation %>%
  group_by(group, productive) %>%
  tally() %>%
  mutate(
    per = paste0(round(n / sum(n) * 100, digits = 1), "%"),
    productive = factor(productive, levels = c("true", "false")),
    n = n / 1e4,
    group = group
  ) %>%
  ggbarplot(
    data = .,
    x = "group",
    y = "n",
    fill = "productive",
    merge = FALSE
  ) +
  geom_label(aes(label = per, fill = productive), size = 8) +
  theme_test() +
  ylab("Cell Counts (10^4)") +
  xlab("") +
  labs(fill = "Productive (V(D)J Assembly Successfully)") +
  scale_fill_manual(values = c("#1fddff", "#ff4b1f")) +
  theme(
    title = element_text(size = 20),
    axis.text.x = element_text(size = 20),
    axis.text.y = element_text(size = 20),
    legend.text = element_text(size = 20),
    legend.position = "bottom"
  ) +
  coord_flip()

# sample level productive vs non-productive contig number
plot_sample_productive <- all_annotation %>%
  group_by(sample_name, productive) %>%
  tally() %>%
  mutate(
    per = paste0(round(n / sum(n) * 100, digits = 1), "%"),
    productive = factor(productive, levels = c("true", "false")),
    n = n / 1e4
  ) %>%
  ggbarplot(
    data = .,
    x = "sample_name",
    y = "n",
    fill = "productive",
    merge = FALSE
  ) +
  geom_label(aes(label = per, fill = productive), size = 8) +
  theme_test() +
  ylab("Cell Counts (10^4)") +
  xlab("") +
  labs(fill = "Productive (V(D)J Assembly Successfully)") +
  scale_fill_manual(values = c("#1fddff", "#ff4b1f")) +
  theme(
    title = element_text(size = 20),
    axis.text.x = element_text(size = 20),
    axis.text.y = element_text(size = 20),
    legend.text = element_text(size = 20),
    legend.position = "bottom"
  ) +
  coord_flip()
save_plot(plot_sample_productive, "sample_productive", 12, 8, overwrite = force_read_data_plot)

# filter according to productive contigs only
all_annotation_filter <- all_annotation %>%
  filter(productive == "true")
# stat TRD and TRG contig number per cell
annotation_stat <-
  all_annotation_filter %>%
  group_by(barcode) %>%
  summarise(trd = sum(chain == "TRD"), trg = sum(chain == "TRG"))
annotation_stat <-
  merge(annotation_stat,
    all_annotation_filter[, c(
      "barcode", "Vd1", "Vd2", "Vd3",
      "Vg2", "Vg3", "Vg4", "Vg5", "Vg8", "Vg9"
    )],
    by = "barcode"
  ) %>%
  unique()
# save filtered annotation data and stat data
save_rds_if_missing(all_annotation_filter, vdj_annotation_filtered_rds, "filtered VDJ annotation RDS", overwrite = force_read_data)
save_rds_if_missing(annotation_stat, vdj_annotation_stat_rds, "VDJ annotation stat RDS", overwrite = force_read_data)


# ...Process expression matrix...
# combine all expression data into one Seurat object
# use merge in Reduce will cause layer names terribly wrong
# like counts.ZTN-Vd1.SeuratProject counts.pan-3.SeuratProject.SeuratProject
# DONT KNOW WHY seperating CreateSeuratObject and then merging can work
all_seurat_list <-
  lapply(seq_along(sample_list), function(i) {
    message("Running SoupX and reading 10x data for sample: ", sample_list[i])
    soupx_result <- run_soupx(
      sample_id = sample_list[i],
      raw_dir = raw_matrix_filepath[i],
      filtered_dir = filtered_matrix_filepath[i],
      out_dir = soupx_dir
    )
    tmp_count <- soupx_result$adjusted_counts
    tmp_se <- CreateSeuratObject(
      tmp_count,
      assay = "RNA",
      names.field = 1,
      names.delim = "_",
      meta.data = NULL,
      project = sample_list[i]
    )
    tmp_se$soupx_contamination <- soupx_result$soup_channel$metaData[colnames(tmp_se), "rho"]
    tmp_se$soupx_umis <- soupx_result$soup_channel$metaData[colnames(tmp_se), "nUMIs"]
    tmp_se
  })
all_seurat <- merge_seurat_samples(all_seurat_list, sample_list)
rm(all_seurat_list)
gc()
# add metadata to seurat object
sample_meta <- meta[match(all_seurat$orig.ident, meta$sample_id), ]
if (any(is.na(sample_meta$sample_id))) {
  missing_samples <- unique(all_seurat$orig.ident[is.na(sample_meta$sample_id)])
  stop("Missing sample metadata for: ", paste(missing_samples, collapse = ", "))
}
all_seurat$sample_name <- sample_meta$sample_name
all_seurat$group <- sample_meta$group
all_seurat$group <-
  factor(all_seurat$group,
    levels = group_levels
  )
all_seurat$sample_name <-
  factor(all_seurat$sample_name,
    levels = sample_name_levels
  )
all_seurat$barcode <- colnames(all_seurat)
vdj_match <- match(all_seurat$barcode, annotation_stat$barcode)
all_seurat$trd <- annotation_stat$trd[vdj_match]
all_seurat$trg <- annotation_stat$trg[vdj_match]
all_seurat$trd[is.na(all_seurat$trd)] <- 0
all_seurat$trg[is.na(all_seurat$trg)] <- 0
# check the number of cells with TRD and TRG per group and sample
annotation_summary <- all_seurat@meta.data %>%
  group_by(group, sample_name) %>%
  summarise(
    n = n(),
    have_either = sum(trd > 0 | trg > 0),
    have_trd = sum(trd > 0),
    have_trg = sum(trg > 0),
    have_both = sum(trd > 0 & trg > 0),
    have_either_percent = have_either / n * 100,
    have_trd_percent = have_trd / n * 100,
    have_trg_percent = have_trg / n * 100,
    have_both_percent = have_both / n * 100
  )
# plot barplot for percentage of cells with TRD and TRG per sample
color_palette_contig <- c(
  "Either" = "#a5a5a5",
  "TRD" = "#1f77b4",
  "TRG" = "#ff7f0e",
  "Both" = "#2ca02c"
)
annotation_summary_barplot <- annotation_summary %>%
  pivot_longer(
    cols = c(
      "have_either_percent",
      "have_trd_percent",
      "have_trg_percent",
      "have_both_percent"
    ),
    names_to = "type",
    values_to = "percent"
  ) %>%
  mutate(type = case_when(
    type == "have_either_percent" ~ "Either",
    type == "have_trd_percent" ~ "TRD",
    type == "have_trg_percent" ~ "TRG",
    type == "have_both_percent" ~ "Both"
  )) %>%
  mutate(type = factor(type, levels = c("Either", "TRD", "TRG", "Both"))) %>%
  ggplot(aes(x = sample_name, y = percent, fill = type)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_bw() +
  coord_flip() +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(size = 10)
  ) +
  facet_wrap(~group, scales = "free_y", nrow = 1) +
  scale_fill_manual(values = color_palette_contig) +
  scale_y_continuous(breaks = seq(0, 100, by = 20)) +
  labs(
    title = "Percentage of cells with TRD and TRG contigs per sample",
    x = "",
    y = "Percentage (%)",
    fill = "Contig Type"
  )
save_plot(annotation_summary_barplot, "VDJ_contig_percentage_per_sample_barplot", 15, 3, overwrite = force_read_data_plot)
# Vd1 cells have much less cells with TRG contigs
# some reasons may be:
# 1. Vd1 T cells have lower expression of TRG genes
# 2. Vd1 T cells have abundant TRG repertoire
# Vd1 can pair with multiple TRGV genes while most Vd2 pair with TRGV9 only
# proved by qt-pcr experiment

# ...Process QC and filter...
# calculate percent.mt and percent.ribo
all_seurat[["percent.mt"]] <-
  PercentageFeatureSet(all_seurat, pattern = "^MT-")
all_seurat[["percent.ribo"]] <-
  PercentageFeatureSet(all_seurat, pattern = "^RP[LS]")
# violin plot for QC metrics
qc_violin_plot <-
  VlnPlot2(all_seurat,
    features = qc_features,
    group.by = "sample_name",
    nrow = 2,
    ncol = 2,
    pt = FALSE,
    cols = color_sample_name
  ) +
  labs(title = "QC Metrics Violin Plot")
save_plot(qc_violin_plot, "QC_metrics_violin_plot_before_filtering", 8, 6, overwrite = force_read_data_plot)
save_rds_if_missing(all_seurat, seurat_before_filtering_rds, "pre-filter Seurat RDS", overwrite = force_read_data)
