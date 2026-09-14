# R version 4.5.2 (2025-10-31)
rm(list = ls())
setwd("/path/to/project")
library(Seurat)
library(SeuratExtend)
library(gplots)
library(tidyverse)
library(patchwork)
library(ggpubr)
options(
  tibble.width = Inf,
  print.width = Inf,
  max.print = 200
)

rds_dir <- "rds"
figure_dir <- file.path("figures", "4_VDJgene")

vdj_annotation_rds <- file.path(rds_dir, "All_VDJ_annotation_filtered.rds")
vdj_stat_rds <- file.path(rds_dir, "VDJ_annotation_stat.rds")
seurat_celltype_rds <- file.path(rds_dir, "all_seurat_celltype.rds")
included_annotation_rds <- file.path(rds_dir, "all_annotation_included.rds")

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
trgc_expression_groups <- c("AB3", "PAN")
force_vdj_gene <- FALSE
force_vdj_gene_plot <- TRUE

dir.create(rds_dir, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

add_celltype_to_annotation <- function(annotation, seurat_obj) {
  annotation %>%
    filter(barcode %in% colnames(seurat_obj)) %>%
    merge(
      seurat_obj@meta.data[, c("barcode", "cell_type")],
      by = "barcode"
    ) %>%
    mutate(
      group = factor(group, levels = group_levels),
      cell_type = factor(cell_type, levels = cell_type_levels)
    )
}

summarise_trdg_per_cell <- function(annotation) {
  annotation %>%
    group_by(barcode, group) %>%
    summarise(
      n_trd = sum(chain == "TRD"),
      n_trg = sum(chain == "TRG"),
      .groups = "drop"
    ) %>%
    group_by(n_trd, n_trg, group) %>%
    tally() %>%
    group_by(group) %>%
    mutate(per = round(n / sum(n) * 100, digits = 2)) %>%
    ungroup()
}

summarise_vj_gene_distribution <- function(annotation) {
  annotation %>%
    group_by(chain, group, v_gene, j_gene) %>%
    tally() %>%
    group_by(chain, group) %>%
    mutate(per = round(n / sum(n) * 100, digits = 2)) %>%
    arrange(chain, group, desc(per)) %>%
    ungroup()
}

summarise_jc_gene_distribution <- function(annotation) {
  annotation %>%
    group_by(chain, group, j_gene, c_gene) %>%
    tally() %>%
    group_by(chain, group) %>%
    mutate(per = round(n / sum(n) * 100, digits = 2)) %>%
    arrange(chain, group, desc(per)) %>%
    ungroup()
}

plot_trgc_expression <- function(seurat_obj) {
  VlnPlot2(
    subset(seurat_obj, group %in% trgc_expression_groups),
    group.by = "group",
    features = c("TRGC1", "TRGC2"),
    pt = FALSE,
    nrow = 2,
    cols = color_group,
    stat.method = "wilcox.test"
  ) +
    ylab("Scaled Expression")
}

plot_trdg_per_cell <- function(trdg_summary) {
  ggplot(trdg_summary, aes(x = n_trd, y = n_trg)) +
    geom_point(
      aes(color = group, size = per),
      alpha = 0.7,
      position = position_dodge(width = 0.5)
    ) +
    scale_size(range = c(1, 20)) +
    scale_color_manual(values = color_group) +
    scale_x_continuous(breaks = 0:3) +
    scale_y_continuous(breaks = 0:3) +
    xlim(-0.5, 2.5) +
    ylim(-0.5, 2.5) +
    theme_bw() +
    labs(
      title = "TRD and TRG number per cell",
      x = "Number of TRD",
      y = "Number of TRG",
      size = "Percentage (%)",
      color = ""
    )
}

plot_vj_gene_distribution <- function(vj_summary, chain_name, height_title) {
  vj_summary %>%
    filter(chain == chain_name) %>%
    ggplot(aes(
      x = v_gene,
      y = j_gene,
      size = per,
      group = group
    )) +
    geom_point(aes(color = group),
      alpha = 0.7,
      position = position_dodge(0.7)
    ) +
    scale_size(range = c(1, 20)) +
    scale_color_manual(values = color_group) +
    theme_bw() +
    labs(
      title = height_title,
      x = "V Gene",
      y = "J Gene",
      size = "Percentage (%)",
      color = ""
    )
}

plot_trg_jc_gene_distribution <- function(jc_summary) {
  jc_summary %>%
    filter(chain == "TRG") %>%
    mutate(c_gene = ifelse(c_gene == "", "NA", c_gene)) %>%
    ggplot(aes(
      x = j_gene,
      y = c_gene,
      size = per,
      group = group
    )) +
    geom_point(aes(color = group),
      alpha = 0.7,
      position = position_dodge(0.7)
    ) +
    scale_size(range = c(1, 20)) +
    scale_color_manual(values = color_group) +
    theme_bw() +
    labs(
      title = "TRG JC Gene Distribution",
      x = "J Gene",
      y = "C Gene",
      size = "Percentage (%)",
      color = ""
    )
}

# 1. Read VDJ annotations and final annotated Seurat object.
all_annotation_filter <- readRDS(vdj_annotation_rds)
annotation_stat <- readRDS(vdj_stat_rds)
all_seurat_celltype <- readRDS(seurat_celltype_rds)
all_seurat_celltype$group <- factor(all_seurat_celltype$group, levels = group_levels)
all_seurat_celltype$cell_type <- factor(all_seurat_celltype$cell_type, levels = cell_type_levels)

# 2. Compare TRGC expression between selected treatment groups.
plot_vln_gc <- plot_trgc_expression(all_seurat_celltype)
save_plot(plot_vln_gc, "VDJ_vln_GC", 8, 5, overwrite = force_vdj_gene_plot)

# 3. Keep only productive VDJ annotations from cells present in the final object.
all_annotation_included <- add_celltype_to_annotation(
  all_annotation_filter,
  all_seurat_celltype
)
save_rds_if_missing(all_annotation_included, included_annotation_rds, "included annotation RDS", overwrite = force_vdj_gene)

# 4. Summarise TRD/TRG chain counts per cell.
stat_ntrdg_percell <- summarise_trdg_per_cell(all_annotation_included)
plot_ntrdg_percell <- plot_trdg_per_cell(stat_ntrdg_percell)
save_plot(plot_ntrdg_percell, "VDJ_nTRDG_percell", 8, 3, overwrite = force_vdj_gene_plot)

# 5. Summarise V/J and J/C gene usage by group.
stat_vj_gene_distribution <- summarise_vj_gene_distribution(all_annotation_included)
plot_trd_vj_gene_distribution <- plot_vj_gene_distribution(
  stat_vj_gene_distribution,
  "TRD",
  "TRD VJ Gene Distribution"
)
plot_trg_vj_gene_distribution <- plot_vj_gene_distribution(
  stat_vj_gene_distribution,
  "TRG",
  "TRG VJ Gene Distribution"
)
save_plot(plot_trd_vj_gene_distribution, "VDJ_TRD_VDJ_gene_distribution", 8, 4, overwrite = force_vdj_gene_plot)
save_plot(plot_trg_vj_gene_distribution, "VDJ_TRG_VDJ_gene_distribution", 8, 3, overwrite = force_vdj_gene_plot)

# TRGJ1 has the same sequence as TRGJ2, so alignment cannot distinguish them.
# TRGJ2-TRGC1 likely represents TRGJ1-TRGC1 because TRGJ2 should not pair with TRGC1.
stat_jc_gene_distribution <- summarise_jc_gene_distribution(all_annotation_included)
plot_trg_jc_gene_distribution <- plot_trg_jc_gene_distribution(stat_jc_gene_distribution)
save_plot(plot_trg_jc_gene_distribution, "VDJ_TRG_JC_gene_distribution", 10, 8, overwrite = force_vdj_gene_plot)

# 6. Save combined VDJ summary figure.
combined_plots <- (plot_ntrdg_percell / plot_trd_vj_gene_distribution) |
  (plot_trg_vj_gene_distribution / plot_trg_jc_gene_distribution) +
    plot_layout(guides = "collect")
save_plot(combined_plots, "VDJ_combined_plots", 24, 8, overwrite = force_vdj_gene_plot)
