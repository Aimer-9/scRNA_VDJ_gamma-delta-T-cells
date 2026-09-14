# TCR single-cell analysis pipeline

This repository contains the code used to process 10x Genomics GEX plus gamma-delta TCR V(D)J data, build annotated Seurat objects, and run downstream repertoire, clonotype, trajectory, and pySCENIC analyses.

The code is organized around two stages:

1. Cell Ranger and FASTQ preparation under `code/cellranger/`.
2. R and shell downstream analysis under `code/downstream/`.

Top-level helper files are grouped by purpose:

```text
scripts/        Shell entry points for setup, FASTQ upload, WeCom messages, and downstream runs
scripts/env/    Local secret/config files ignored by git
docs/           Setup notes, review notes, and repository structure notes
metadata/       Submission/helper metadata such as file lists and MCP tool schemas
code/           Analysis source code and pipeline configuration
renv/           Local R package library
```

They call `setwd()` to that path and expect runtime folders such as `config/`, `rds/`, `figures/`, `table/`, `cache/`, and `pyscenic/` to exist there.

## Quick Start

Prepare or update the sample/config files:

```text
code/config/samples.csv
code/config/params.yaml
```

Prepare local runtime folders, software links, references, and R package library
from the same params file when needed:

```bash
bash scripts/setup.sh check
bash scripts/setup.sh all --dry-run
bash scripts/setup.sh all
```

`setup.sh` is parameter-driven. The default setup parts are defined by
`setup.parts` in `code/config/params.yaml`, and individual parts can be run
without running everything:

```bash
bash scripts/setup.sh env
bash scripts/setup.sh software
bash scripts/setup.sh ref
bash scripts/setup.sh r
bash scripts/setup.sh --part software --part ref
```

The setup workflow is local/offline by design: Cell Ranger, FastQC, references,
and R packages are installed only from paths or local archives configured in
`params.yaml`. See `docs/setup.md` for all setup keys and examples.

Link FASTQs into the project raw-data directory:

```bash
bash code/cellranger/src/00_find_fq.sh \
  --dir /path/to/raw_fastq_root \
  --outdir /path/to/project/rawdata
```

Run the Cell Ranger wrapper:

```bash
bash scripts/cellranger.sh dry-run \
  --params code/config/params.yaml \
  --samples code/config/samples.csv

bash scripts/cellranger.sh qc \
  --params code/config/params.yaml \
  --samples code/config/samples.csv

bash scripts/cellranger.sh multi \
  --params code/config/params.yaml \
  --samples code/config/samples.csv
```

Prepare FASTQs for submission with `ftp.sh` in three steps. The default output root is the current directory, so integrated FASTQs are kept under `./integrated_fastq/` and upload files are written under `./upload_tar/`:

```bash
bash scripts/ftp.sh integrate --rawdata rawdata --result . --jobs 20
bash scripts/ftp.sh tar --rawdata rawdata --result . --jobs 20
bash scripts/ftp.sh check --result .
bash scripts/ftp.sh ftp --result .
```

If integration is interrupted, rerun the integrate step. `ftp.sh` removes incomplete `*.fastq.gz.tmp.*` files, checks expected integrated FASTQs in parallel, skips existing non-empty files, and rebuilds missing sample/read outputs. The tar step also checks that all expected integrated FASTQs exist before packaging.

Example directory layout after integration and tar:

```text
.
+-- rawdata/
|   +-- sample-5LIB_S7_L004_I1_001.fastq.gz
|   +-- sample-5LIB_S7_L004_I2_001.fastq.gz
|   +-- sample-5LIB_S7_L004_R1_001.fastq.gz
|   +-- sample-5LIB_S7_L004_R2_001.fastq.gz
+-- integrated_fastq/
|   +-- sample-5LIB_I1_001.fastq.gz
|   +-- sample-5LIB_I2_001.fastq.gz
|   +-- sample-5LIB_R1_001.fastq.gz
|   +-- sample-5LIB_R2_001.fastq.gz
+-- upload_tar/
    +-- sample-5LIB.tar.gz
    +-- manifest.tsv
```

For samples that only have `R1/R2`, `scripts/ftp.sh tar` copies the integrated files to `./upload_tar/plain_pe/` by default. Add `--tar-pe` when paired-end-only samples must also be packed into `.tar.gz`.

To send automatic emails after each step succeeds or fails, configure `msmtp` first. Keep account values in `scripts/env/ftp.env` or pass them through environment variables:

```bash
bash scripts/ftp.sh test-email
bash scripts/ftp.sh ftp --result .
```

The full workflow sends one email after each completed or failed step:

```bash
bash scripts/ftp.sh all --rawdata rawdata --result . --jobs 20
```

FTP and mail credentials are loaded from `scripts/env/ftp.env` by default. Use `FTP_CONFIG=/path/to/ftp.env` or `--config /path/to/ftp.env` to use another config file. Override the recipient with `NOTIFY_EMAIL` or `--email`, and override the sender with `MAIL_FROM` or `--mail-from` when needed.

Run downstream scripts in numeric order after Cell Ranger outputs are available. The downstream scripts are intended to be run from the project runtime directory or with matching deployed paths:

```bash
Rscript code/downstream/1_ReadData.R
Rscript code/downstream/2_DataClean.R
Rscript code/downstream/3_cellAnnotation.R
Rscript code/downstream/4_VDJgene.R
Rscript code/downstream/5_CDR3stat.R
Rscript code/downstream/6_CDR3betweenSamples.R
Rscript code/downstream/7_CDR3pairedSankey.R
Rscript code/downstream/8_Vd1vs2.R
Rscript code/downstream/9_CDR3paired.R
Rscript code/downstream/10_MSH2.R
Rscript code/downstream/13_Vd2_pseudotime.R
Rscript code/downstream/14_Vd1Vd2_pairwise.R
Rscript code/downstream/15_CD80_CD86_expression.R
Rscript code/downstream/16_Vd1Vd2_extra_visualization.R
Rscript code/downstream/17_ZOL_PAN_effector_Vd2_comparison.R
```

To run the downstream R scripts in order and send one final email on success or failure, use:

```bash
bash scripts/run_downstream.sh
```

## Dependency-Tracked Downstream Builds

The compatibility pipeline under `code/pipeline/` uses `targets` to record
dependencies and rebuild only affected stages. It preserves the current
`rds/`, `table/`, and `figures/` artifact paths while the numbered scripts are
being migrated into modules.

Configure project-specific downstream paths and settings in
`code/config/analysis.yaml`, then run:

```bash
bash scripts/run_pipeline.sh
bash scripts/run_pipeline.sh --target step_10_msh2
```

Existing automation can opt in without changing its step syntax:

```bash
bash scripts/run_downstream.sh --targets --step 10
```

`targets` must be available in the local R package source before running the
new entry point; `scripts/setup.sh r` reads `renv.lock` and installs it through
the existing offline setup process.

Each run writes to a new directory named with time and PID under `downstream_runs/`, for example `downstream_runs/20260825_142301_pid12345/`. `run_downstream.sh` copies the downstream scripts into that run directory and patches their working directory, so R outputs such as `rds/`, `figures/`, `table/`, `soupx/`, and optional `pyscenic_10k/` are created under the run directory. Logs are written under `logs/` in the same run directory, with one log file per script plus `run_summary.tsv`, `run_info.tsv`, and `run_message.txt`. On failure, the email includes the failed script, failed log path, and the last 60 log lines. Set `NOTIFY_EMAIL` and `MAIL_FROM` to enable notifications.

To run one step only, use `--step`. To write into an existing run directory, pass either the run name under `downstream_runs/` or the full path:

```bash
bash scripts/run_downstream.sh --step 10
bash scripts/run_downstream.sh --step 10_MSH2.R
bash scripts/run_downstream.sh --existing-run 20260825_213811_pid1045941 --step 10
bash scripts/run_downstream.sh --resume-run downstream_runs/20260825_213811_pid1045941 --step 10_MSH2.R
bash scripts/run_downstream.sh --existing-run 20260825_213811_pid1045941 --step 11 \
  --tf-list /path/to/allTFs_hg38.txt \
  --ranking-db /path/to/hg38_10kb.feather \
  --ranking-db /path/to/hg38_500bp.feather \
  --motif-annotations /path/to/motifs.tbl \
  --pyscenic-python /path/to/conda/env/bin/python \
  --force-selection \
  --force-export
```

Optional pySCENIC scripts can be included after the main R workflow:

```bash
bash scripts/run_downstream.sh --include-pyscenic
```

Most downstream R scripts now skip existing non-empty data outputs by default. Plot and heatmap overwrite flags default to `TRUE`, so figures are refreshed when the script body runs; set the relevant `*_plot` flag to `FALSE` when you want existing figures to be kept. Major scripts such as `1_ReadData.R`, `2_DataClean.R`, and `3_cellAnnotation.R` still exit early when their final RDS outputs already exist unless their plot flag or main force flag is enabled.

To rebuild a script's outputs without manually deleting files, set that script's internal force value to `TRUE` before running it:

- `1_ReadData.R`: `force_read_data`; `force_read_data_plot` redraws figures when the script is run
- `2_DataClean.R`: `force_data_clean`; `force_data_clean_plot` redraws figures when the script is run
- `3_cellAnnotation.R`: `force_cell_annotation`; use `force_cell_annotation_plot` to redraw figures without rebuilding filtering outputs
- `4_VDJgene.R`: `force_vdj_gene`; `force_vdj_gene_plot`
- `5_CDR3stat.R`: `force_cdr3_stat`; `force_cdr3_stat_plot`
- `6_CDR3betweenSamples.R`: `force_cdr3_between_samples`; `force_cdr3_between_samples_plot`
- `7_CDR3pairedSankey.R`: `force_cdr3_paired_sankey`; `force_cdr3_paired_sankey_plot`
- `8_Vd1vs2.R`: `force_vd1_vs_vd2`; `force_vd1_vs_vd2_plot`
- `9_CDR3paired.R`: `force_ranked_rds`; `force_ranked_plot`
- `10_MSH2.R`: `force_msh2`; `force_msh2_plot`
- `13_Vd2_pseudotime.R`: `force_vd2_pseudotime`; `force_vd2_pseudotime_plot`
- `14_Vd1Vd2_pairwise.R`: `force_vd1_vd2_pairwise`; `force_vd1_vd2_pairwise_plot`
- `15_CD80_CD86_expression.R`: `force_cd80_cd86_plot`
- `16_Vd1Vd2_extra_visualization.R`: `force_vd1_vd2_extra`; `force_vd1_vd2_extra_plot`
- `17_ZOL_PAN_effector_Vd2_comparison.R`: `force_zol_pan_effector_vd2`; `force_zol_pan_effector_vd2_plot`
- `12_pySCENIC_visualization.R`: `force_pyscenic_visualization` or CLI `--overwrite`; use `force_pyscenic_visualization_plot` or CLI `--overwrite-plots` for figures only

`code/downstream/11_pySCENIC.sh` is the pySCENIC workflow. It selects a balanced 10,000-cell subset by `cell_type`, exports expression with sparse-aware chunks, and runs the Python pySCENIC steps through `code/downstream/lib/pyscenic_run.py`; use `--pyscenic-python` to run the Python part from a specific conda environment.

## Runtime Inputs

Core inputs expected by the pipeline:

- `config/samples.csv`: sample metadata and FASTQ prefixes.
- `config/params.yaml`: Cell Ranger references, executable paths, output root, thread/memory settings, and optional local setup settings.
- `scripts/setup.sh` and `docs/setup.md`: parameter-driven local setup for runtime folders, software links, references, and R package library preparation.
- Cell Ranger multi outputs under the path configured in `1_ReadData.R` and `params.yaml`.
- `code/downstream/lib/plotting_shared.R`: shared group/sample/cell-type levels, palettes, and plot-saving helpers.
- Optional pySCENIC resources: TF list, cisTarget ranking databases, and motif annotation table.

The repository currently stores the shared plotting file at:

```text
code/downstream/lib/plotting_shared.R
```

The downstream scripts first look for this repository path, then `lib/plotting_shared.R` inside copied run scripts, then legacy `_cache_`/`cache` fallback paths.

## Main Outputs

Important RDS outputs:

- `rds/All_VDJ_annotation_filtered.rds`
- `rds/VDJ_annotation_stat.rds`
- `rds/all_seurat_before_filtering.rds`
- `rds/all_seurat_2.rds`
- `rds/all_seurat_3.rds`
- `rds/all_seurat_4.rds`
- `rds/all_seurat_celltype.rds`
- `rds/all_annotation_included.rds`
- `rds/barcode_trgd_paired.rds`
- `rds/trdg_pair_rank_metadata.rds`
- `rds/all_seurat_celltype_toprank_cells.rds`
- `rds/pr_graph_test_res_toprank.rds`
- `table/toprank_vd2_trdg_pairs_for_monocle3.csv`

`3_cellAnnotation.R` now reruns from `rds/all_seurat_2.rds` when `force_cell_annotation <- TRUE` and saves full intermediate Seurat objects as `rds/all_seurat_3.rds` and `rds/all_seurat_4.rds`. Set `force_cell_annotation_plot <- TRUE` to overwrite final figures from an existing `rds/all_seurat_celltype.rds` without rebuilding filtering outputs or marker tables.

For scripts that already start from saved RDS/table inputs, the `*_plot` flags redraw figures without forcing cached RDS/CSV outputs to be rebuilt. `1_ReadData.R` and `2_DataClean.R` include raw/SoupX/scDblFinder diagnostic plots whose intermediate objects are not cached separately, so their plot flags prevent data-output overwrites but still require rerunning the script body.

Figures are written to step-specific subdirectories under `figures/`, for example `figures/5_CDR3stat/`; marker and summary tables are written to `table/`; pySCENIC outputs are written to `pyscenic/` by default.

`9_CDR3paired.R` treats Monocle3 `graph_test()` as optional. If R package `sf` cannot load because the system PROJ library is missing, for example `libproj.so.15`, the script skips `graph_test`, the pseudotime heatmap, and graph-test Hallmark enrichment while still saving the top-pair cell metadata, hallmark heatmaps, pseudotime UMAP, and cell-cycle UMAP. To produce `rds/pr_graph_test_res_toprank.rds`, fix the R environment by installing the matching PROJ runtime library or reinstalling `sf` against the available PROJ version.

`9_CDR3paired.R` also writes `table/trdg_clone_dispersion_metrics.csv` to describe how dispersed each exact paired TRD+TRG clone is. For each clone, the script first finds all cells with the same `TRD||TRG` pair. UMAP dispersion is calculated as the mean Euclidean distance from each clone cell to that clone's UMAP centroid, using `umap.unintegrated` when available and otherwise `umap`. Hallmark pathway dispersion is calculated from Hallmark AUCell scores: for each Hallmark pathway, the script calculates the standard deviation across cells in the clone, then averages those pathway-level standard deviations. The two raw metrics are min-max scaled across clones as `umap_mean_distance_scaled` and `hallmark_mean_pathway_sd_scaled`; `dispersion_score` is the mean of those two scaled values. Group, sample, and cell-type entropy/evenness values are saved as descriptive annotations, but they do not contribute to `dispersion_score`.

`13_Vd2_pseudotime.R` runs a Monocle3 trajectory on about 5,000 Vd2 cells from exact paired TRD+TRG clones, using the same clone-ranking rules as `9_CDR3paired.R`: each selected clone must be a single-cell TRD/TRG pair and the same pair must appear in Naive, ZOL, and PAN groups. It roots pseudotime in `Effector Memory Vd2` and saves UMAP, violin, density, sample-summary, gene-expression, gene-pseudotime, and pseudotime-expression heatmap visualizations.

`14_Vd1Vd2_pairwise.R` performs non-pseudotime pairwise comparisons for `Naive Vd1` versus `Effector Memory Vd2`, and `Effector Vd1` versus `PAN Effector Vd2`. It saves differential-expression tables, marker-expression summaries, UMAP highlights, volcano plots, marker dotplots, expression heatmaps, sample-level expression plots, and Hallmark AUCell delta plots.

`15_CD80_CD86_expression.R` plots CD80 and CD86 expression on fixed-coordinate UMAPs and condition-group dotplots from the final annotated Seurat object.

`16_Vd1Vd2_extra_visualization.R` adds figure-first Vd1/Vd2 comparison panels: UMAP state highlights, sample-level cell fractions, curated marker dotplots and heatmaps, module-score summaries, expanded DE volcano/overlap/top-gene heatmaps, and repertoire-aware clone-size, top-CDR3, paired-clone alluvial, and paired-clone sharing heatmap outputs.

`17_ZOL_PAN_effector_Vd2_comparison.R` focuses on `ZOL Effector Vd2` versus `PAN Effector Vd2`: UMAP context/density, sample-level abundance and ratio plots, curated marker and module-score summaries, DE volcano/lollipop/top-gene heatmap, Hallmark pathway delta/selected heatmap, and clone-size/top-CDR3/paired-clone repertoire views. Its module-score violin additionally includes `ZOL FOXP3+ Vd2` and shows all pairwise state comparisons; the remaining Step 17 analyses retain the ZOL-versus-PAN scope.

## Downstream Order

The downstream dependency chain is:

```text
1_ReadData.R
  -> 2_DataClean.R
    -> 3_cellAnnotation.R
      -> 4_VDJgene.R
        -> 5_CDR3stat.R
        -> 6_CDR3betweenSamples.R
        -> 7_CDR3pairedSankey.R
        -> 8_Vd1vs2.R
          -> 9_CDR3paired.R
          -> 10_MSH2.R
          -> 13_Vd2_pseudotime.R
          -> 14_Vd1Vd2_pairwise.R
          -> 15_CD80_CD86_expression.R
          -> 16_Vd1Vd2_extra_visualization.R
          -> 17_ZOL_PAN_effector_Vd2_comparison.R
      -> 11_pySCENIC.sh
        -> 12_pySCENIC_visualization.R
```

## Cell Type Convention

The current cell type vocabulary is:

- `Effector Memory Vd2`
- `Pre-activated Vd2`
- `ZOL Effector Vd2`
- `ZOL FOXP3+ Vd2`
- `PAN Effector Vd2`
- `Naive Vd1`
- `Pre-activated Vd1`
- `Effector Vd1`

The palette and factor levels are centralized in `code/downstream/lib/plotting_shared.R`.

## pySCENIC

Run `code/downstream/11_pySCENIC.sh` after `rds/all_seurat_celltype.rds` exists and an environment with `pyscenic` is active.
The default run is intentionally conservative: it selects a reproducible balanced 10,000-cell subset, runs GRN with one worker, and uses `custom_multiprocessing` for ctx to avoid the common local Dask failures seen on full-size matrices.

Create the pinned Python environment:

```bash
conda env create -f code/downstream/envs/pyscenic_environment.yml
conda activate pyscenic
```

Download default human pySCENIC references when needed:

```bash
bash code/downstream/11_pySCENIC.sh download-ref \
  --out-dir /path/to/pyscenic_ref \
  --species human

source /path/to/pyscenic_ref/pyscenic_reference.env
```

```bash
bash code/downstream/11_pySCENIC.sh \
  --tf-list /path/to/allTFs_hg38.txt \
  --ranking-db /path/to/hg38_10kb.feather \
  --ranking-db /path/to/hg38_500bp.feather \
  --motif-annotations /path/to/motifs.tbl \
  --pyscenic-python /path/to/conda/env/bin/python \
  --target-cells 10000 \
  --group-by cell_type \
  --min-cells-per-group 200 \
  --grn-num-workers 1 \
  --num-workers 16
```

The script writes `pyscenic_10k/selected_cells_10000_by_cell_type.tsv`, exports exactly those cells, and then runs GRN, ctx, and AUCell.

The same values can still be supplied through environment variables (`TF_LIST`, `RANKING_DBS`, `MOTIF_ANNOTATIONS`, `NUM_WORKERS`, and the other options shown by `--help`).
GRN inference uses Arboreto/Dask internally, so the wrapper defaults to `--grn-num-workers 1` and `--seed 1` for reproducibility and robustness; `--num-workers` still controls ctx/AUCell.
Use `--force-selection --force-export` when changing `--target-cells`, `--group-by`, or filtering thresholds, otherwise existing selected cells or exports may be reused.
For a larger run after the 10k-cell run succeeds, increase `--target-cells`, but expect much longer runtime and a higher chance of Dask worker failure.
The default `--ctx-mode custom_multiprocessing` avoids Dask failures such as `TypeError: Must supply at least one delayed object`; use `--ctx-mode dask_multiprocessing` only after confirming the pinned Dask environment works.
The Python runner also checks that TF symbols overlap expression-matrix gene columns before launching pySCENIC, since a species or gene-symbol mismatch can surface as the same Dask error.
If Dask raises `AttributeError: module 'pandas.core.strings' has no attribute 'StringMethods'`, the Python part is running with pandas 2.x or another incompatible environment; rebuild or update from `code/downstream/envs/pyscenic_environment.yml` and pass that environment with `--pyscenic-python`.
If `pyscenic_grn.log` shows Dask warnings such as `Could not find data` or `Worker ... failed to acquire keys`, reduce `--grn-num-workers` to 2 or 1 and consider lowering `--target-cells`; these warnings usually mean local Dask workers died or restarted during GRN.
If `pyscenic_10k/adjacencies.tsv` is empty or header-only, GRN produced zero TF-target edges. Check `pyscenic_10k/logs/pyscenic_grn.log`, confirm the reported TF overlap is not near zero, then rerun with `--grn-num-workers 1 --seed 1`; if it still has zero edges, try `--grn-method genie3` or relax expression filtering with lower `--min-counts-per-gene`/`--min-cells-per-gene`.
This pySCENIC version requires AUCell expression input in loom format. The runner automatically converts `expression_matrix.csv` to `expression_matrix.loom` in bounded cell chunks, while GRN and ctx continue using CSV. If visualization reports that `auc_mtx.loom` is not an HDF5 file, inspect `pyscenic_10k/logs/pyscenic_aucell.log`, remove the invalid `pyscenic_10k/auc_mtx.loom`, and rerun `11_pySCENIC.sh`. The runner removes stale AUCell output and validates the HDF5 signature and required loom datasets before reporting success.
The runner suppresses the known Python multiprocessing shutdown warning about leaked semaphore objects; it is a cleanup warning from worker shutdown, not a pySCENIC output validation failure.
The Seurat expression matrix is exported in sparse-aware cell chunks to avoid a full sparse-to-dense conversion; lower `--export-chunk-cells` if memory is still tight.

Default outputs are:

- `pyscenic_10k/selected_cells_10000_by_cell_type.tsv`
- `pyscenic_10k/expression_matrix.csv`
- `pyscenic_10k/expression_matrix.loom`
- `pyscenic_10k/cell_metadata.tsv`
- `pyscenic_10k/adjacencies.tsv`
- `pyscenic_10k/regulons.csv`
- `pyscenic_10k/auc_mtx.loom`
- `pyscenic_10k/logs/pyscenic_select_cells.log`
- `pyscenic_10k/logs/pyscenic_export.log`
- `pyscenic_10k/logs/pyscenic_grn.log`
- `pyscenic_10k/logs/pyscenic_ctx.log`
- `pyscenic_10k/logs/pyscenic_aucell.log`

Visualize pySCENIC regulon activity after `auc_mtx.loom` is available:

```bash
Rscript code/downstream/12_pySCENIC_visualization.R \
  --auc-loom pyscenic_10k/auc_mtx.loom \
  --seurat-rds rds/all_seurat_celltype.rds \
  --group-by cell_type \
  --n-top-regulons 30 \
  --n-umap-regulons 12
```

Visualization outputs are written to `figures/12_pySCENIC_visualization/`, `table/pyscenic/`, and `rds/pyscenic_visualization_result.rds`.
The script generates a grouped regulon z-score heatmap, mean-AUC/activity dot plot, top-regulon UMAP panels, grouped activity tables, selected-regulon list, and a reusable result RDS.
It also generates a curated panel for `BACH2`, `SOX4`, `NFATC1`, `MTF1`, `TFEB`, `NR3C1`, `GABPB1`, `NFATC3`, `IRF1`, `FOXP3`, `STAT3`, `XBP1`, `CEBPB`, `TFDP1`, and `HMGA1`. This panel includes a focused heatmap, dot plot, per-cell regulon-AUC UMAP grid, TF gene-expression Feature UMAP grid, and publication-style interpretation table containing the associated cell types, established T-cell functions, and predicted Vd1/Vd2 interpretations. Regulon names such as `BACH2(+)` and `BACH2_extended(+)` are matched to the TF symbol; TFs absent from the loom remain in the interpretation CSV/table and are marked as not detected.
For SOX4 specifically, the script also generates three relationship views: a faceted SOX4-target network across T-cell states, a SOX4 target-gene expression heatmap, and scatter plots of SOX4 regulon AUC versus target-gene expression. Positive and negative labels are inferred from within-state Spearman correlation and should be interpreted as association, not proof of direct activation or repression.
It requires the R packages `hdf5r`, `ComplexHeatmap`, `circlize`, `Seurat`, `dplyr`, `tidyr`, `readr`, `ggplot2`, and `patchwork`.

## Review Notes

- `scripts/cellranger.sh` is the Cell Ranger entry point and uses implementation files under `code/cellranger/src`.
- Downstream R scripts use hardcoded `setwd("/path/to/project")`; change this in each script or reproduce that runtime path.
- Downstream R scripts source `code/downstream/lib/plotting_shared.R` with fallbacks for `_cache_/plotting_shared.R` and `cache/plotting_shared.R`.
