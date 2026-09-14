library(targets)

source("code/pipeline/R/config.R")
source("code/pipeline/R/legacy.R")

tar_option_set(packages = c("yaml"), format = "rds", seed = 1234)

list(
  tar_target(config_files, c("code/config/analysis.yaml", "code/config/params.yaml", "code/config/samples.csv"), format = "file"),
  tar_target(legacy_sources, file.path("code", "downstream", legacy_steps), format = "file"),
  tar_target(config, { config_files; read_pipeline_config() }),
  tar_target(samples, read_and_validate_samples(config)),
  tar_target(step_01_read_data, { legacy_sources; samples; run_legacy_step(config, "1_ReadData.R") }, format = "file"),
  tar_target(step_02_data_clean, { legacy_sources; run_legacy_step(config, "2_DataClean.R", step_01_read_data) }, format = "file"),
  tar_target(step_03_cell_annotation, { legacy_sources; run_legacy_step(config, "3_cellAnnotation.R", step_02_data_clean) }, format = "file"),
  tar_target(step_04_vdj_gene, { legacy_sources; run_legacy_step(config, "4_VDJgene.R", step_03_cell_annotation) }, format = "file"),
  tar_target(step_05_cdr3_stat, { legacy_sources; run_legacy_step(config, "5_CDR3stat.R", step_04_vdj_gene) }, format = "file"),
  tar_target(step_06_cdr3_between_samples, { legacy_sources; run_legacy_step(config, "6_CDR3betweenSamples.R", step_04_vdj_gene) }, format = "file"),
  tar_target(step_07_cdr3_paired_sankey, { legacy_sources; run_legacy_step(config, "7_CDR3pairedSankey.R", step_04_vdj_gene) }, format = "file"),
  tar_target(step_08_vd1_vs_vd2, { legacy_sources; run_legacy_step(config, "8_Vd1vs2.R", step_03_cell_annotation) }, format = "file"),
  tar_target(step_09_cdr3_paired, { legacy_sources; run_legacy_step(config, "9_CDR3paired.R", step_07_cdr3_paired_sankey) }, format = "file"),
  tar_target(step_10_msh2, { legacy_sources; run_legacy_step(config, "10_MSH2.R", step_09_cdr3_paired) }, format = "file"),
  tar_target(step_13_vd2_pseudotime, { legacy_sources; run_legacy_step(config, "13_Vd2_pseudotime.R", step_09_cdr3_paired) }, format = "file"),
  tar_target(step_14_vd1_vd2_pairwise, { legacy_sources; run_legacy_step(config, "14_Vd1Vd2_pairwise.R", step_03_cell_annotation) }, format = "file"),
  tar_target(step_15_cd80_cd86, { legacy_sources; run_legacy_step(config, "15_CD80_CD86_expression.R", step_03_cell_annotation) }, format = "file"),
  tar_target(step_16_vd1_vd2_extra, { legacy_sources; run_legacy_step(config, "16_Vd1Vd2_extra_visualization.R", c(step_04_vdj_gene, step_07_cdr3_paired_sankey)) }, format = "file"),
  tar_target(step_17_zol_pan_effector_vd2, { legacy_sources; run_legacy_step(config, "17_ZOL_PAN_effector_Vd2_comparison.R", c(step_04_vdj_gene, step_07_cdr3_paired_sankey)) }, format = "file")
)
