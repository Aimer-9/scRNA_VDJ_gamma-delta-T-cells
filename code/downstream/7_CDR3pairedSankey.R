# R version 4.5.2 (2025-10-31)
rm(list = ls())
setwd("/data/huotong/project_tcr/2026May")
library(tidyverse)
library(ggpubr)
library(RColorBrewer)
library(ggalluvial)
options(
  max.print = 100,
  tibble.width = Inf,
  spe = "human"
)

rds_dir <- "rds"
figure_dir <- file.path("figures", "7_CDR3pairedSankey")

included_annotation_rds <- file.path(rds_dir, "all_annotation_included.rds")
paired_cdr3_rds <- file.path(rds_dir, "barcode_trgd_paired.rds")

top_pair_n_per_sample <- 30
force_cdr3_paired_sankey <- FALSE
force_cdr3_paired_sankey_plot <- TRUE

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

dir.create(rds_dir, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

normalise_annotation_levels <- function(annotation) {
  annotation %>%
    mutate(
      group = factor(group, levels = group_levels),
      sample_name = factor(sample_name, levels = sample_name_levels),
      cell_type = factor(cell_type, levels = cell_type_levels)
    )
}

make_paired_cdr3 <- function(annotation) {
  annotation %>%
    select(barcode, chain, cdr3) %>%
    distinct() %>%
    pivot_wider(names_from = chain, values_from = cdr3, values_fn = list) %>%
    split(.$barcode) %>%
    lapply(function(record) {
      trd <- unlist(record$TRD)
      trg <- unlist(record$TRG)
      if (length(trd) == 0 || length(trg) == 0) {
        return(NULL)
      }
      cbind(
        barcode = unique(record$barcode),
        expand.grid(TRD = trd, TRG = trg)
      )
    }) %>%
    bind_rows()
}

add_pair_metadata <- function(paired_cdr3, annotation) {
  cell_metadata <- annotation %>%
    select(barcode, group, sample_name, cell_type) %>%
    distinct()
  paired_cdr3 %>%
    merge(cell_metadata, by = "barcode") %>%
    distinct()
}

summarise_top_pairs <- function(paired_cdr3) {
  paired_cdr3 %>%
    group_by(sample_name, group, TRD, TRG) %>%
    tally() %>%
    group_by(sample_name) %>%
    mutate(percent = n / sum(n) * 100) %>%
    slice_max(order_by = percent, n = top_pair_n_per_sample) %>%
    ungroup()
}

get_pair_factor_levels <- function(top_pairs, pair_column) {
  top_pairs %>%
    group_by(.data[[pair_column]]) %>%
    summarise(total_percent = sum(percent), .groups = "drop") %>%
    arrange(desc(total_percent)) %>%
    pull(.data[[pair_column]])
}

summarise_pair_sample_counts <- function(top_pairs, group_columns) {
  top_pairs %>%
    group_by(across(all_of(group_columns))) %>%
    summarise(sample_n = length(unique(sample_name)), .groups = "drop") %>%
    arrange(desc(sample_n))
}

plot_pair_sankey_within_group <- function(top_pairs) {
  trd_factor <- get_pair_factor_levels(top_pairs, "TRD")
  trg_factor <- get_pair_factor_levels(top_pairs, "TRG")
  sample_counts <- summarise_pair_sample_counts(top_pairs, c("TRD", "TRG", "group"))

  top_pairs %>%
    merge(sample_counts, by = c("TRD", "TRG", "group")) %>%
    ggplot(aes(
      y = percent,
      axis1 = factor(TRD, levels = trd_factor),
      axis2 = sample_name,
      axis3 = factor(TRG, levels = trg_factor)
    )) +
    geom_stratum(width = 1 / 6) +
    geom_alluvium(aes(
      fill = sample_name,
      alpha = ifelse(sample_n > 1, 0.5, 0)
    ), width = 1 / 12) +
    geom_text(
      stat = "stratum",
      aes(label = after_stat(ifelse(
        stratum %in% as.character(top_pairs$sample_name),
        as.character(stratum),
        ""
      ))),
      size = 3
    ) +
    theme_test() +
    scale_x_continuous(
      breaks = 1:3,
      labels = c("TRD", "Sample Name", "TRG")
    ) +
    scale_fill_manual(values = color_sample_name) +
    facet_wrap(~group, ncol = 1, scales = "free_y") +
    theme(
      legend.position = "none",
      strip.background = element_blank(),
      axis.title.y = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank()
    )
}

plot_pair_sankey_among_groups <- function(top_pairs) {
  trd_factor <- get_pair_factor_levels(top_pairs, "TRD")
  trg_factor <- get_pair_factor_levels(top_pairs, "TRG")
  sample_counts <- summarise_pair_sample_counts(top_pairs, c("TRD", "TRG"))

  top_pairs %>%
    merge(sample_counts, by = c("TRD", "TRG")) %>%
    ggplot(aes(
      y = n,
      axis1 = factor(TRD, levels = trd_factor),
      axis2 = sample_name,
      axis3 = factor(TRG, levels = trg_factor)
    )) +
    geom_stratum(width = 1 / 6) +
    geom_alluvium(aes(
      fill = group,
      alpha = ifelse(sample_n > 1, 0.5, 0)
    ), width = 1 / 12) +
    geom_text(
      stat = "stratum",
      aes(label = after_stat(ifelse(
        stratum %in% as.character(top_pairs$sample_name),
        as.character(stratum),
        ""
      ))),
      size = 5
    ) +
    theme_test() +
    scale_x_continuous(
      breaks = 1:3,
      labels = c("TRD", "Sample Name", "TRG")
    ) +
    scale_fill_manual(values = color_group) +
    theme(
      legend.position = "none",
      axis.title.y = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank()
    )
}

# 1. Read productive VDJ annotations and build paired TRD/TRG CDR3 records.
all_annotation_included <- readRDS(included_annotation_rds)
all_annotation_included <- normalise_annotation_levels(all_annotation_included)

barcode_trgd_paired <- make_paired_cdr3(all_annotation_included)
paired_cell_fraction <- length(unique(barcode_trgd_paired$barcode)) /
  length(unique(all_annotation_included$barcode))

barcode_trgd_paired <- add_pair_metadata(barcode_trgd_paired, all_annotation_included)
save_rds_if_missing(
  barcode_trgd_paired,
  paired_cdr3_rds,
  "paired TRD/TRG CDR3 RDS",
  overwrite = force_cdr3_paired_sankey
)

# 2. Plot top TRD/TRG pair sharing within and across groups.
top_trgd_pair <- summarise_top_pairs(barcode_trgd_paired)
plot_sankey_trdg_within_group <- plot_pair_sankey_within_group(top_trgd_pair)
save_plot(plot_sankey_trdg_within_group, "sankey_TRDG_within_group", 5, 15, overwrite = force_cdr3_paired_sankey_plot)

plot_sankey_trdg_among_group <- plot_pair_sankey_among_groups(top_trgd_pair)
save_plot(plot_sankey_trdg_among_group, "sankey_TRDG_among_group", 10, 6, overwrite = force_cdr3_paired_sankey_plot)
