source(testthat::test_path("..", "..", "code", "pipeline", "R", "config.R"))
source(testthat::test_path("..", "..", "code", "pipeline", "R", "legacy.R"))

testthat::test_that("sample metadata rejects duplicate sample identifiers", {
  root <- tempfile("tcr-config-")
  dir.create(file.path(root, "code", "config"), recursive = TRUE)
  writeLines("outdir: output", file.path(root, "code/config/params.yaml"))
  writeLines(paste0("project_dir: ", root), file.path(root, "code/config/analysis.yaml"))
  writeLines("sample_id,group,timepoint,gex_fastq_path,vdj_t_gd_fastq_path,gex_fastq_prefix,vdj_t_gd_fastq_prefix\na,A,1,x,x,x,x\na,A,2,x,x,x,x", file.path(root, "code/config/samples.csv"))
  config <- read_pipeline_config(root)
  testthat::expect_error(read_and_validate_samples(config), "unique")
})

testthat::test_that("legacy artifact contract preserves established RDS paths", {
  config <- list(artifacts = list(rds_dir = "rds", table_dir = "table"))
  testthat::expect_equal(legacy_step_outputs(config, "2_DataClean.R"), "rds/all_seurat_2.rds")
  testthat::expect_true("rds/trdg_pair_rank_metadata.rds" %in% legacy_step_outputs(config, "9_CDR3paired.R"))
})
