# Local Setup

`scripts/setup.sh` prepares the local project environment from `code/config/params.yaml`.
It is designed for offline/local setup: software, references, and R packages are
installed only from paths or archives you provide in params.

## Commands

```bash
bash scripts/setup.sh check
bash scripts/setup.sh env
bash scripts/setup.sh software
bash scripts/setup.sh ref
bash scripts/setup.sh r
bash scripts/setup.sh all
```

Run selected parts:

```bash
bash scripts/setup.sh --part software --part ref
```

Preview without changing files:

```bash
bash scripts/setup.sh all --dry-run
```

Use a different params file:

```bash
bash scripts/setup.sh all --params code/config/params.yaml
```

## Params

The setup section in `code/config/params.yaml` controls optional local setup:

```yaml
setup:
  parts: ["env", "software", "ref", "r"]
  project_dir: "/path/to/project"
  software_dir: "/path/to/software"
  reference_dir: "/path/to/reference/cellranger"
  bin_dir: "tools/bin"
  conda_env_prefix: ""
  conda_env_yaml: ""
  r_lib: "renv/library"
  r_package_source_dir: ""
  cellranger_archive: ""
  fastqc_archive: ""
  gex_reference_archive: ""
  vdj_reference_archive: ""
```

Existing Cell Ranger keys are still used:

```yaml
cellranger:
  cellranger_path: "/path/to/cellranger"
  fastqc_path: "/path/to/fastqc"
  reference_gex: "/path/to/refdata-gex"
  reference_vdj: "/path/to/refdata-vdj"
```

## What Each Part Does

`env`

- Creates project directories: `rds`, `figures`, `table`, `log`, `tmp`
- Creates `outdir`, `setup.bin_dir`, and `setup.r_lib`
- Optionally creates a conda environment when both `setup.conda_env_prefix` and
  `setup.conda_env_yaml` are set

`software`

- Checks `cellranger.cellranger_path` and `cellranger.fastqc_path`
- If missing, extracts `setup.cellranger_archive` or `setup.fastqc_archive`
- Creates local symlinks in `setup.bin_dir`

`ref`

- Checks `cellranger.reference_gex` and `cellranger.reference_vdj`
- If missing, extracts `setup.gex_reference_archive` or
  `setup.vdj_reference_archive` into `setup.reference_dir`

`r`

- Uses `renv.lock` as the required package list
- Creates `setup.r_lib`
- Installs missing packages only from local archives in
  `setup.r_package_source_dir`

`check`

- Reports whether configured software, reference, and primer paths exist

## Notes

- `scripts/setup.sh` requires `python3` and Python package `PyYAML` to read params.
- R setup requires `Rscript` and package `jsonlite`.
- The script does not use CRAN/Bioconductor/network package installation.
- By default existing software/reference paths are not re-extracted. Use
  `--no-skip-existing` to force extraction from local archives.
