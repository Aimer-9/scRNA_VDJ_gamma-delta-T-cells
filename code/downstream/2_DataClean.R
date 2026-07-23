# R version 4.5.2 (2025-10-31)
rm(list = ls())
setwd("/data/huotong/project_tcr/2026May")
library(Seurat)
library(SeuratExtend)
library(scDblFinder)
library(SeuratDisk)
library(readxl)
library(tidyverse)
library(patchwork)
library(ggpubr)
library(gplots)
options(
  tibble.width = Inf,
  print.width = Inf,
  max.print = 200
)

metadata_file <- "config/samples.csv"
rds_dir <- "rds"
figure_dir <- file.path("figures", "2_DataClean")

vdj_annotation_rds <- file.path(rds_dir, "All_VDJ_annotation_filtered.rds")
vdj_stat_rds <- file.path(rds_dir, "VDJ_annotation_stat.rds")
seurat_before_filtering_rds <- file.path(rds_dir, "all_seurat_before_filtering.rds")
seurat_filtered_rds <- file.path(rds_dir, "all_seurat_2.rds")

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
# Blank-2 has less than 1K cells, exclude it.
excluded_sample_id <- c("Blank-2")
qc_features <- c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo")
max_trd_chains <- 2
max_trg_chains <- 2
top_cdr3_n <- 100
umap_dims <- 1:20
noise_colors <- c("FALSE" = "blue", "TRUE" = "red")
chain_count_colors <- c(
  "0" = "lightblue",
  "1" = "lightgreen",
  "2" = "lightyellow",
  "2+" = "lightcoral"
)
force_data_clean <- FALSE
force_data_clean_plot <- TRUE

dir.create(rds_dir, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

if (output_files_exist(seurat_filtered_rds) && !force_data_clean && !force_data_clean_plot) {
  message("[SKIP] Existing filtered Seurat RDS: ", seurat_filtered_rds)
  quit(save = "no")
} else if (output_files_exist(seurat_filtered_rds) && force_data_clean) {
  message("[FORCE] Rebuilding filtered Seurat RDS from current script settings: ", seurat_filtered_rds)
} else if (output_files_exist(seurat_filtered_rds) && force_data_clean_plot) {
  message("[PLOT] Re-drawing 2_DataClean figures without forcing filtered Seurat RDS overwrite.")
}

meta <- read.csv(metadata_file)
meta$group <- factor(meta$group, levels = group_levels)
meta$sample_name <- factor(
  paste0(meta$group, "_", meta$timepoint),
  levels = sample_name_levels
)
meta <- meta %>%
  filter(!sample_id %in% excluded_sample_id)

plot_qc_violin <- function(seurat_obj, title) {
  VlnPlot2(seurat_obj,
    features = qc_features,
    group.by = "sample_name",
    nrow = 2,
    ncol = 2,
    pt = FALSE,
    cols = color_sample_name
  ) +
    labs(title = title)
}

filter_included_samples <- function(seurat_obj, annotation, metadata) {
  included_sample_ids <- metadata$sample_id
  included_sample_names <- metadata$sample_name

  seurat_obj <- subset(
    seurat_obj,
    subset = orig.ident %in% included_sample_ids & sample_name %in% included_sample_names
  )
  seurat_obj$group <- factor(seurat_obj$group, levels = group_levels)
  seurat_obj$sample_name <- factor(seurat_obj$sample_name, levels = sample_name_levels)

  annotation <- annotation %>%
    filter(sample %in% included_sample_ids & sample_name %in% included_sample_names)
  annotation$group <- factor(annotation$group, levels = group_levels)
  annotation$sample_name <- factor(annotation$sample_name, levels = sample_name_levels)

  list(seurat = seurat_obj, annotation = annotation)
}

summarise_chain_counts <- function(annotation) {
  annotation %>%
    group_by(barcode, sample_name, group) %>%
    summarise(
      trd = sum(chain == "TRD"),
      trg = sum(chain == "TRG"),
      .groups = "drop"
    )
}

summarise_cdr3_by_sample <- function(annotation, count_col = "n") {
  annotation %>%
    group_by(sample_name, chain, cdr3, group) %>%
    summarise(total_umis = sum(umis), n = n(), .groups = "drop_last") %>%
    group_by(sample_name, chain) %>%
    mutate(percent = n / sum(n) * 100) %>%
    arrange(sample_name, chain, desc(.data[[count_col]])) %>%
    mutate(rank = row_number()) %>%
    ungroup()
}

get_noise_barcodes <- function(annotation_stat) {
  list(
    trd = annotation_stat %>%
      filter(trd > max_trd_chains) %>%
      pull(barcode),
    trg = annotation_stat %>%
      filter(trg > max_trg_chains) %>%
      pull(barcode)
  )
}

get_noise_cdr3_pool <- function(annotation, barcodes, chain_name) {
  annotation %>%
    filter(barcode %in% barcodes & chain == chain_name) %>%
    pull(cdr3) %>%
    unique()
}

plot_noise_cdr3 <- function(annotation, chain_name, noise_pool, legend_position) {
  summarise_cdr3_by_sample(annotation) %>%
    mutate(is_noise = ifelse(cdr3 %in% noise_pool, TRUE, FALSE)) %>%
    filter(chain == chain_name) %>%
    ggplot(aes(x = rank, y = percent, color = is_noise)) +
    geom_point(size = 0.5, alpha = 0.5) +
    facet_wrap(~sample_name, scales = "free_y") +
    theme_test() +
    scale_color_manual(values = noise_colors) +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      legend.position = legend_position
    ) +
    labs(title = paste(chain_name, "CDR3 Noise Marking"))
}

add_doublet_metadata <- function(seurat_obj, doublet_obj) {
  doublet_columns <- c(
    "scDblFinder.class", "scDblFinder.score",
    "scDblFinder.weighted", "scDblFinder.cxds_score"
  )
  seurat_obj@meta.data[, doublet_columns] <- doublet_obj@meta.data[, doublet_columns]
  seurat_obj
}

# read in
all_annotation_filter <- readRDS(vdj_annotation_rds)
annotation_stat <- readRDS(vdj_stat_rds)
all_seurat <- readRDS(seurat_before_filtering_rds)
included_data <- filter_included_samples(all_seurat, all_annotation_filter, meta)
all_seurat <- included_data$seurat
all_annotation_filter <- included_data$annotation
rm(included_data)

# 1. filter cells with no TRD/G
# cells without VDJ annotation are worthless in later analysis.
# what's more, cells without TRD/G annotation usually have low quality.
all_seurat_1 <- subset(all_seurat, trd > 0 | trg > 0)
# plot violn plot to show cells filtered with standard 1
qc_violin_plot_filter1 <- plot_qc_violin(
  all_seurat_1,
  "QC Metrics Violin Plot (Filtered: TRD>0 | TRG>0)"
)
save_plot(qc_violin_plot_filter1, "QC_metrics_violin_plot_after_filter1", 8, 6, overwrite = force_data_clean_plot)

# 2. filter cells with too many TRD/G
# impossible for diploid single γδT cell to have > 2 chains
# due to allelic exclusion, most cells should only 1 TRD/G
# cells wtih more than 2 TRD or TRG are likely doublets or multiplets
# 2.1 try to figure out where contaminant is from

# 2.1.1 stat TRD/TRG chain count in single cell
all_annotation_included1 <- all_annotation_filter %>%
  filter(barcode %in% all_seurat_1@meta.data$barcode)
annotation_stat <-
  summarise_chain_counts(all_annotation_included1)
# quick check stat data
chain_count_balloon_trd <- table(annotation_stat$trd, annotation_stat$sample_name) %>%
  balloonplot(ylab = "Sample Name", xlab = "TRD Count")
chain_count_balloon_trg <- table(annotation_stat$trg, annotation_stat$sample_name) %>%
  balloonplot(ylab = "Sample Name", xlab = "TRG Count")
# most cells have only 1 TRD
# except many AB3 cells present 0 TRG, most cells also have only 1 TRG
# TRG has less limited allelic exclusion than TRD found in researches
# cutoff setting 1 or 2 need further thoughts
# if possible, removing contaminant cdr3 could increase data quality
# let cells contatminated recovered.

# 2.1.2 assume cdr3s in cells with > 2 TRD/G marked with possible noise source
# figure out where the noise is from
noise_barcodes <- get_noise_barcodes(annotation_stat)
obviously_noise_barcode_trd <- noise_barcodes$trd
obviously_noise_barcode_trg <- noise_barcodes$trg

noise_trd_pool <- get_noise_cdr3_pool(
  all_annotation_included1,
  obviously_noise_barcode_trd,
  "TRD"
)
noise_trd_pool %>% length()
all_annotation_included1 %>%
  filter(chain == "TRD") %>%
  pull(cdr3) %>%
  unique() %>%
  length()
noise_markered_plot_trd <- plot_noise_cdr3(
  all_annotation_included1,
  "TRD",
  noise_trd_pool,
  "none"
)

noise_trg_pool <- get_noise_cdr3_pool(
  all_annotation_included1,
  obviously_noise_barcode_trg,
  "TRG"
)
noise_trg_pool %>% length()
all_annotation_included1 %>%
  filter(chain == "TRG") %>%
  pull(cdr3) %>%
  unique() %>%
  length()
noise_markered_plot_trg <- plot_noise_cdr3(
  all_annotation_included1,
  "TRG",
  noise_trg_pool,
  "bottom"
)
noise_marked_plot <- noise_markered_plot_trd + noise_markered_plot_trg +
  plot_layout(guides = "collect")
save_plot(noise_marked_plot, "noise_cdr3_marked_plot", 12, 8, overwrite = force_data_clean_plot)
# it shows noise cdr3s are usually highly frequent cdr3 in each sample.
# it makes sense but also makes it difficult to distinguish from normal cdr3.
# it is also doubtful that cells with only 1 high-frequent cdr3 are true or not

# 2.1.3 check top100 frequency cdr3s' umis in cells with 2 more TRD/G
# if noise cdr3s have lower umis than true cdr3, it makes easier to denoise
top_cdr3_samples <- summarise_cdr3_by_sample(
  all_annotation_included1,
  count_col = "total_umis"
) %>%
  group_by(sample_name, chain) %>%
  slice_head(n = top_cdr3_n) %>%
  ungroup() %>%
  mutate(pair = paste0(sample_name, "_", cdr3))
all_annotation_included_unfiltered <- all_annotation_included1 %>%
  filter(barcode %in% c(
    obviously_noise_barcode_trd,
    obviously_noise_barcode_trg
  )) %>%
  group_by(sample_name, chain) %>%
  mutate(
    barcode_rank = dense_rank(barcode),
    cdr3_top = ifelse(paste0(sample_name, "_", cdr3) %in%
      top_cdr3_samples$pair,
      TRUE, FALSE
    )
  )

noise_cdr3_umi_plot <- all_annotation_included_unfiltered %>%
  ggplot(aes(x = barcode_rank, y = umis, color = cdr3_top)) +
  geom_point(size = 0.3, alpha = 0.5) +
  theme(
    axis.text.x = element_blank()
  ) +
  facet_grid(chain ~ sample_name, scales = "free") +
  theme_test() +
  scale_color_manual(values = noise_colors) +
  labs(color = "Is Top100 CDR3")
save_plot(noise_cdr3_umi_plot, "noise_cdr3_umi_plot", 12, 8, overwrite = force_data_clean_plot)
# boxplot and compare with wilcox test
boxplot_umi_cdr3_top <- all_annotation_included_unfiltered %>%
  ggplot(aes(x = cdr3_top, y = umis, color = cdr3_top)) +
  geom_boxplot() +
  facet_grid(chain ~ sample_name, scales = "free") +
  theme_test() +
  stat_compare_means(
    label = "p.signif",
    comparisons = list(c("FALSE", "TRUE"))
  ) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    strip.background = element_blank()
  ) +
  scale_color_manual(values = noise_colors) +
  labs(color = "Is Top100 CDR3") +
  xlab("")
save_plot(boxplot_umi_cdr3_top, "boxplot_umi_cdr3_top", 12, 8, overwrite = force_data_clean_plot)

# sad
# 1. contaminant cdr3s are usually top frequency cdr3s in each samples.
# 2. contaminant cdr3s are not consistently with lower umis in noise cells.
# it seems impossible to denoise.
# Or just remove TRD > 2 cells or TRG > 2 cells
# well.

# 2.1.4 check proportion of cells to be removed
annotation_stat_2 <- annotation_stat %>%
  mutate(
    trd = ifelse(trd > max_trd_chains, "2+", as.character(trd)),
    trg = ifelse(trg > max_trg_chains, "2+", as.character(trg))
  ) %>%
  group_by(sample_name, trd, trg) %>%
  tally()
annotation_stat_2_trd <- annotation_stat_2 %>%
  ggplot(aes(x = sample_name, y = n, fill = trd)) +
  geom_col(position = "fill") +
  theme_test() +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  ) +
  labs(fill = "TRD count", title = "TRD counts across samples") +
  ylab("Proportion of cells") +
  scale_fill_manual(values = chain_count_colors)
annotation_stat_2_trg <- annotation_stat_2 %>%
  ggplot(aes(x = sample_name, y = n, fill = trg)) +
  geom_col(position = "fill") +
  theme_test() +
  theme(axis.title.x = element_blank()) +
  labs(fill = "TRG count", title = "TRG counts across samples") +
  ylab("Proportion of cells") +
  scale_fill_manual(values = chain_count_colors)
annotation_plot_nTRDG <- annotation_stat_2_trd / annotation_stat_2_trg +
  plot_annotation(title = "TRD and TRG counts across samples") +
  plot_layout(guides = "collect")
save_plot(annotation_plot_nTRDG, "annotation_plot_nTRDG", 12, 8, overwrite = force_data_clean_plot)

# 2.2 remove cells with 2 more TRD/G
all_seurat_2 <- subset(all_seurat_1, trd <= max_trd_chains & trg <= max_trg_chains)

# cell number
all_seurat %>% dim()
all_seurat_1 %>% dim()
all_seurat_2 %>% dim()

qc_violin_plot_filtered <- plot_qc_violin(
  all_seurat_2,
  "QC Metrics Violin Plot Filtered"
)
save_plot(qc_violin_plot_filtered, "qc_violin_plot_filtered", 12, 8, overwrite = force_data_clean_plot)

# 3. doublet finder
all_seurat_2 <- NormalizeData(all_seurat_2)
all_seurat_2 <- FindVariableFeatures(all_seurat_2)
all_seurat_2 <- ScaleData(all_seurat_2)

all_seurat_2_join <-
  JoinLayers(all_seurat_2)
all_singlecellexperiment <-
  as.SingleCellExperiment(all_seurat_2_join)
all_singlecellexperiment <-
  scDblFinder(all_singlecellexperiment, samples = "orig.ident")
all_seurat_doublet <- as.Seurat(all_singlecellexperiment)
all_seurat_2 <- add_doublet_metadata(all_seurat_2, all_seurat_doublet)
# run pca and umap, check doublet status
all_seurat_2 <- RunPCA(all_seurat_2)
ElbowPlot(all_seurat_2, ndims = 30)
all_seurat_2 <- RunUMAP(all_seurat_2,
  dims = umap_dims,
  reduction = "pca",
  reduction.name = "umap.unintegrated"
)
umap_doublet_class <- DimPlot2(all_seurat_2,
  group.by = c("scDblFinder.class"),
  reduction = "umap.unintegrated",
  theme = NoAxes(),
  label = TRUE, box = TRUE, label.color = "black", repel = TRUE
)
umap_doublet_class <- add_fixed_umap_coordinates(umap_doublet_class)
umap_doublet_score <- DimPlot2(all_seurat_2,
  group.by = c("scDblFinder.score"),
  reduction = "umap.unintegrated",
  theme = NoAxes(),
  label = TRUE, box = TRUE, label.color = "black", repel = TRUE
)
umap_doublet_score <- add_fixed_umap_coordinates(umap_doublet_score)
umap_doublet <- umap_doublet_class + umap_doublet_score
save_plot(umap_doublet, "umap_doublet", 12, 8, overwrite = force_data_clean_plot)

save_rds_if_missing(all_seurat_2, seurat_filtered_rds, "filtered Seurat RDS", overwrite = force_data_clean)
