pipeline_abort <- function(...) {
  stop(..., call. = FALSE)
}

`%||%` <- function(x, y) if (is.null(x)) y else x

pipeline_root <- function() {
  configured <- Sys.getenv("TCR_PROJECT_DIR", unset = "")
  if (nzchar(configured)) return(normalizePath(configured, mustWork = TRUE))
  normalizePath(getwd(), mustWork = TRUE)
}

resolve_pipeline_path <- function(path, root) {
  if (!nzchar(path)) return("")
  if (grepl("^/", path)) return(normalizePath(path, mustWork = FALSE))
  normalizePath(file.path(root, path), mustWork = FALSE)
}

read_pipeline_config <- function(root = pipeline_root()) {
  if (!requireNamespace("yaml", quietly = TRUE)) {
    pipeline_abort("Package `yaml` is required. Restore the project environment first.")
  }
  analysis_path <- file.path(root, "code", "config", "analysis.yaml")
  params_path <- file.path(root, "code", "config", "params.yaml")
  samples_path <- file.path(root, "code", "config", "samples.csv")
  for (path in c(analysis_path, params_path, samples_path)) {
    if (!file.exists(path)) pipeline_abort("Required pipeline configuration is missing: ", path)
  }
  analysis <- yaml::read_yaml(analysis_path) %||% list()
  params <- yaml::read_yaml(params_path) %||% list()
  root <- resolve_pipeline_path(analysis$project_dir %||% ".", root)
  output_root <- resolve_pipeline_path(params$outdir %||% "output", root)
  cellranger_output <- analysis$cellranger_output %||% ""
  if (!nzchar(cellranger_output)) cellranger_output <- file.path(output_root, "cellranger_output")
  artifacts <- analysis$artifacts %||% list()
  defaults <- list(rds_dir = "rds", table_dir = "table", figure_dir = "figures", soupx_dir = "soupx", pyscenic_dir = "pyscenic_10k")
  list(
    root = root,
    config_files = c(analysis_path, params_path, samples_path),
    samples_path = samples_path,
    cellranger_output = resolve_pipeline_path(cellranger_output, root),
    excluded_sample_ids = unlist(analysis$excluded_sample_ids %||% character(), use.names = FALSE),
    seed = as.integer(analysis$seed %||% 1234L),
    artifacts = lapply(modifyList(defaults, artifacts), resolve_pipeline_path, root = root)
  )
}

read_and_validate_samples <- function(config) {
  samples <- utils::read.csv(config$samples_path, stringsAsFactors = FALSE, check.names = FALSE)
  required <- c("sample_id", "group", "timepoint", "gex_fastq_path", "vdj_t_gd_fastq_path", "gex_fastq_prefix", "vdj_t_gd_fastq_prefix")
  missing <- setdiff(required, names(samples))
  if (length(missing)) pipeline_abort("samples.csv is missing required column(s): ", paste(missing, collapse = ", "))
  if (anyNA(samples$sample_id) || any(!nzchar(samples$sample_id)) || anyDuplicated(samples$sample_id)) {
    pipeline_abort("samples.csv must contain unique, non-empty sample_id values.")
  }
  samples
}
