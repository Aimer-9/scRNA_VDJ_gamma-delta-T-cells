# `code/config/params.yaml` Parameter Reference

This document explains every key currently used in `code/config/params.yaml`.

## Top-level

| Key | Type | Meaning |
|---|---|---|
| `outdir` | string | Root output directory for pipeline outputs. |

## `cellranger`

| Key | Type | Meaning |
|---|---|---|
| `cellranger.reference_gex` | string | Path to Cell Ranger gene-expression reference (GRCh38 reference package). |
| `cellranger.reference_vdj` | string | Path to Cell Ranger V(D)J reference package. |
| `cellranger.cellranger_path` | string | Path to the `cellranger` executable. |
| `cellranger.fastqc_path` | string | Path to the `fastqc` executable used in FASTQ QC. |
| `cellranger.fastqc_threads` | integer | CPU threads passed to FastQC. |
| `cellranger.inner_enrichment_primers` | string | Path to primer list used for relevant V(D)J enrichment workflows. |
| `cellranger.localcores` | integer | CPU cores given to Cell Ranger runs. |
| `cellranger.localmem` | integer | Memory (GB) given to Cell Ranger runs. |

## `setup`

Optional local setup controls used by `setup.sh`.

| Key | Type | Meaning |
|---|---|---|
| `setup.parts` | list | Default setup parts run by `bash scripts/setup.sh all`; valid values are `env`, `software`, `ref`, and `r`. |
| `setup.project_dir` | string | Project working directory to create/check. |
| `setup.software_dir` | string | Directory where local software archives are extracted. |
| `setup.reference_dir` | string | Directory where local Cell Ranger reference archives are extracted. |
| `setup.bin_dir` | string | Local executable symlink directory, relative paths are resolved from the repository root. |
| `setup.conda_env_prefix` | string | Optional conda environment prefix to create. Leave empty to skip conda environment creation. |
| `setup.conda_env_yaml` | string | Optional conda environment YAML used with `setup.conda_env_prefix`. |
| `setup.r_lib` | string | R library path for local package installation; relative paths are resolved from the repository root. |
| `setup.r_package_source_dir` | string | Directory containing local R package source archives such as `pkg_1.0.0.tar.gz`. Leave empty to only report missing packages. |
| `setup.cellranger_archive` | string | Optional local Cell Ranger archive to extract if `cellranger.cellranger_path` is missing. |
| `setup.fastqc_archive` | string | Optional local FastQC archive to extract if `cellranger.fastqc_path` is missing. |
| `setup.gex_reference_archive` | string | Optional local GEX reference archive to extract if `cellranger.reference_gex` is missing. |
| `setup.vdj_reference_archive` | string | Optional local VDJ reference archive to extract if `cellranger.reference_vdj` is missing. |
