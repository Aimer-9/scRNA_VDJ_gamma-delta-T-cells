legacy_steps <- c(
  "1_ReadData.R", "2_DataClean.R", "3_cellAnnotation.R", "4_VDJgene.R", "5_CDR3stat.R",
  "6_CDR3betweenSamples.R", "7_CDR3pairedSankey.R", "8_Vd1vs2.R", "9_CDR3paired.R",
  "10_MSH2.R"
)

legacy_step_outputs <- function(config, step) {
  rds <- config$artifacts$rds_dir
  table <- config$artifacts$table_dir
  switch(step,
    "1_ReadData.R" = file.path(rds, c("All_VDJ_annotation_filtered.rds", "VDJ_annotation_stat.rds", "all_seurat_before_filtering.rds")),
    "2_DataClean.R" = file.path(rds, "all_seurat_2.rds"),
    "3_cellAnnotation.R" = file.path(rds, "all_seurat_celltype.rds"),
    "4_VDJgene.R" = file.path(rds, "all_annotation_included.rds"),
    "7_CDR3pairedSankey.R" = file.path(rds, "barcode_trgd_paired.rds"),
    "9_CDR3paired.R" = file.path(rds, c("trdg_pair_rank_metadata.rds", "all_seurat_celltype_toprank_cells.rds", "pr_graph_test_res_toprank.rds")),
    character()
  )
}

assert_nonempty_files <- function(paths, label) {
  missing <- paths[!file.exists(paths) | file.info(paths)$size == 0]
  if (length(missing)) pipeline_abort(label, " did not create expected artifact(s): ", paste(missing, collapse = ", "))
  paths
}

run_legacy_step <- function(config, step, prerequisites = character()) {
  assert_nonempty_files(prerequisites, paste0("Prerequisite for ", step))
  script <- file.path(config$root, "code", "downstream", step)
  if (!file.exists(script)) pipeline_abort("Unknown legacy downstream step: ", step)
  output_paths <- legacy_step_outputs(config, step)
  if (!length(output_paths)) {
    marker_dir <- file.path(config$root, "_targets", "legacy-markers")
    dir.create(marker_dir, recursive = TRUE, showWarnings = FALSE)
    output_paths <- file.path(marker_dir, paste0(tools::file_path_sans_ext(step), ".done"))
  }
  targets_dir <- file.path(config$root, "_targets")
  dir.create(targets_dir, recursive = TRUE, showWarnings = FALSE)
  run_dir <- tempfile(paste0("tcr-", tools::file_path_sans_ext(step), "-"), tmpdir = targets_dir)
  dir.create(run_dir, recursive = TRUE, showWarnings = FALSE)
  prepared_script <- file.path(run_dir, basename(script))
  if (!file.copy(script, prepared_script, overwrite = TRUE)) pipeline_abort("Unable to prepare legacy script: ", script)
  contents <- paste(readLines(prepared_script, warn = FALSE), collapse = "\n")
  escaped_root <- gsub("\\\\", "\\\\\\\\", config$root)
  escaped_root <- gsub('"', '\\\\"', escaped_root)
  escaped_cellranger <- gsub("\\\\", "\\\\\\\\", config$cellranger_output)
  escaped_rds <- gsub("\\\\", "\\\\\\\\", config$artifacts$rds_dir)
  escaped_table <- gsub("\\\\", "\\\\\\\\", config$artifacts$table_dir)
  escaped_figure <- gsub("\\\\", "\\\\\\\\", config$artifacts$figure_dir)
  escaped_soupx <- gsub("\\\\", "\\\\\\\\", config$artifacts$soupx_dir)
  excluded <- paste(shQuote(config$excluded_sample_ids), collapse = ", ")
  contents <- gsub('setwd\\("/path/to/project"\\)', paste0('setwd("', escaped_root, '")'), contents)
  contents <- gsub('project_dir <- "/path/to/project"', paste0('project_dir <- "', escaped_root, '"'), contents)
  contents <- gsub('cellranger_dir <- "/path/to/cellranger/output"', paste0('cellranger_dir <- "', escaped_cellranger, '"'), contents)
  contents <- gsub('rds_dir <- "rds"', paste0('rds_dir <- "', escaped_rds, '"'), contents, fixed = TRUE)
  contents <- gsub('table_dir <- "table"', paste0('table_dir <- "', escaped_table, '"'), contents, fixed = TRUE)
  contents <- gsub('soupx_dir <- "soupx"', paste0('soupx_dir <- "', escaped_soupx, '"'), contents, fixed = TRUE)
  contents <- gsub('figure_dir <- file.path\\("figures",', paste0('figure_dir <- file.path("', escaped_figure, '",'), contents)
  contents <- gsub('excluded_sample_id <- c\\("Blank-2"\\)', paste0('excluded_sample_id <- c(', excluded, ')'), contents)
  contents <- gsub('set.seed\\(1234\\)', paste0('set.seed(', config$seed, ')'), contents)
  writeLines(contents, prepared_script)
  status <- system2("Rscript", c("--vanilla", prepared_script), stdout = TRUE, stderr = TRUE)
  if (!identical(attr(status, "status") %||% 0L, 0L)) pipeline_abort("Legacy step failed: ", step, "\n", paste(status, collapse = "\n"))
  if (grepl("_targets/legacy-markers", output_paths[1], fixed = TRUE)) writeLines(format(Sys.time(), tz = "UTC"), output_paths)
  assert_nonempty_files(output_paths, paste0("Legacy step ", step))
}
