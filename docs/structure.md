# Repository Structure

This document reviews every script under `code/` and describes how the pipeline pieces connect.

## Top-Level Layout

```text
.
|-- README.md
|-- code/
|   |-- run.md
|   |-- config/
|   |-- cellranger/
|   `-- downstream/
|-- docs/
|   |-- setup.md
|   |-- structure.md
|   `-- review.txt
|-- metadata/
|   |-- files.txt
|   `-- tools.json
|-- scripts/
|   |-- cellranger.sh
|   |-- setup.sh
|   |-- ftp.sh
|   |-- run_downstream.sh
|   |-- wecom_mcp.sh
|   `-- env/
|-- figures/
|-- table/
`-- manuscript/
```

## Dependency-Tracked Pipeline

`code/pipeline/_targets.R` is the new downstream dependency graph. It first
validates the merged configuration and sample manifest, then records one
compatibility target per numbered R stage. The target layer invokes a prepared
copy of the legacy script with configured project and Cell Ranger paths, and
checks its established output artifacts. This keeps current output contracts
stable while analysis functions are extracted into `code/pipeline/R/`.

`code/config/analysis.yaml` owns downstream-only settings such as the project
root, Cell Ranger output, exclusions, seed, and artifact directories.

`scripts/env/` contains local account and credential files and is ignored by git. `metadata/` contains helper file lists and MCP tool schema exports. `figures/` and `table/` contain generated outputs. Downstream R scripts write figures under step-specific subdirectories such as `figures/5_CDR3stat/`. The main runtime object directory, `rds/`, is expected by the scripts but is not present in this repository snapshot.

## `code/config`

| Path | Role |
|---|---|
| `code/config/params.yaml` | Runtime parameters for Cell Ranger and FASTQ QC: output root, references, executable paths, cores, memory, FastQC settings. |
| `code/config/samples.csv` | Sample metadata and FASTQ path/prefix columns for GEX and VDJ libraries. |
| `code/config/inner-enrichment-primers.txt` | Primer list for VDJ-T-GD Cell Ranger multi configs. |
| `code/config/params_reference.md` | Broader parameter dictionary. It includes keys not used by the current scripts. |

## `scripts/cellranger.sh` And `code/cellranger`

### `scripts/cellranger.sh`

Wrapper for the Cell Ranger preparation workflow.

Subcommands:

- `dry-run`: validate params/samples and print planned commands.
- `qc`: run FASTQ validation, FastQC, and MultiQC.
- `multi`: generate per-sample Cell Ranger multi configs, then run `cellranger multi`.
- `all`: run `qc` then `multi`.

Inputs:

- `--params code/config/params.yaml`
- `--samples code/config/samples.csv`
- `--template code/cellranger/src/multi_template.csv`

Outputs:

- `<outdir>/fastq_qc/`
- `<outdir>/cellranger_configs/`
- `<outdir>/cellranger_output/`

Review note: the wrapper was corrected to use `code/cellranger/src`, matching the repository layout.

### `src/00_find_fq.sh`

Recursively finds `.fastq.gz` and `.fq.gz` files and creates symlinks in a target directory.

Important behavior:

- Optionally prompts to rename source `.fq.gz` files to `.fastq.gz` in interactive mode.
- Detects symlinks containing `Clean/clean` and prompts before removing them.
- Does not overwrite existing files or conflicting symlinks.

### `src/01_fastq_qc.sh`

Validates FASTQ paths and prefixes from `samples.csv`, writes manifests, and optionally runs FastQC plus MultiQC.

Inputs:

- `samples.csv`
- `params.yaml`
- FastQC executable from `params.yaml` or `--fastqc-path`
- Conda environment named `multiqc`

Outputs:

- `<outdir>/fastq_manifest.tsv`
- `<outdir>/fastq_validation_summary.tsv`
- `<outdir>/fastqc/`
- `<outdir>/multiqc/`

### `src/02_generate_cellranger_configs.py`

Generates one Cell Ranger multi config CSV per sample using `samples.csv`, `params.yaml`, and `multi_template.csv`.

Inputs:

- `params.yaml`
- `samples.csv`
- `multi_template.csv`

Outputs:

- `<outdir>/cellranger_configs/<sample_id>_multi_config.csv`

Review note: the project-root detector was corrected to look for `cellranger/src`.

### `src/03_cellranger_multi.sh`

Runs `cellranger multi` sequentially for generated config CSVs.

Important behavior:

- Uses a lock file, `<outdir>/cellranger_output/cellranger_done.lock`, to skip completed samples.
- Skips samples with existing `_SUCCESS` or `outs/` directories.
- Writes per-sample logs under `<outdir>/cellranger_output/logs/`.

### `src/multi_template.csv`

Template for Cell Ranger multi config generation. Contains placeholders for GEX, VDJ-T, VDJ-B, and VDJ-T-GD libraries.

## Shared Downstream Helpers

### `downstream/_cache_/plotting_shared.R`

Defines shared:

- `group_levels`
- `sample_name_levels`
- `cell_type_levels`
- `color_group`
- `color_sample_name`
- `color_celltype`
- `save_plot()`
- `save_heatmap()`
- `make_umap_plot()`

Runtime note: downstream scripts search for `code/downstream/_cache_/plotting_shared.R`, `_cache_/plotting_shared.R`, then the legacy `cache/plotting_shared.R` path.

## Downstream R Pipeline

Most downstream R scripts skip existing non-empty data outputs by default. Plot and heatmap overwrite flags default to `TRUE`, so figures are refreshed when the script body runs; set the relevant `*_plot` flag to `FALSE` when you want existing figures to be kept. To regenerate one step's data outputs, set that script's internal main force value to `TRUE` before running:

- `1_ReadData.R`: `force_read_data`; `force_read_data_plot`
- `2_DataClean.R`: `force_data_clean`; `force_data_clean_plot`
- `3_cellAnnotation.R`: `force_cell_annotation`; `force_cell_annotation_plot` redraws figures only
- `4_VDJgene.R`: `force_vdj_gene`; `force_vdj_gene_plot`
- `5_CDR3stat.R`: `force_cdr3_stat`; `force_cdr3_stat_plot`
- `6_CDR3betweenSamples.R`: `force_cdr3_between_samples`; `force_cdr3_between_samples_plot`
- `7_CDR3pairedSankey.R`: `force_cdr3_paired_sankey`; `force_cdr3_paired_sankey_plot`
- `8_Vd1vs2.R`: `force_vd1_vs_vd2`; `force_vd1_vs_vd2_plot`
- `9_CDR3paired.R`: `force_ranked_rds`; `force_ranked_plot`
- `10_MSH2.R`: `force_msh2`; `force_msh2_plot`
- `12_pySCENIC_visualization.R`: `force_pyscenic_visualization` or CLI `--overwrite`; `force_pyscenic_visualization_plot` or CLI `--overwrite-plots`

The `*_plot` flags overwrite figure files without forcing RDS/CSV data outputs. These plot flags now default to `TRUE`. For `1_ReadData.R` and `2_DataClean.R`, some diagnostic plots depend on raw intermediate objects that are not cached separately, so the script body still has to run; data outputs are not overwritten unless the main force flag is also `TRUE`.

### `downstream/1_ReadData.R`

Builds initial VDJ annotations and a merged Seurat object from Cell Ranger outputs.

Inputs:

- `config/samples.csv`
- Cell Ranger outputs under `cellranger_dir`
- Raw and filtered GEX matrices
- `all_contig_annotations.csv` files from VDJ-T-GD outputs

Main processing:

- Reads and merges productive/nonproductive VDJ annotations.
- Adds V gene flags such as `Vd1`, `Vd2`, `Vg4`, `Vg9`.
- Runs SoupX correction per sample.
- Builds and merges Seurat objects.
- Adds sample metadata, TRD/TRG chain counts, mitochondrial/ribosomal QC metrics.

Outputs:

- `rds/All_VDJ_annotation_filtered.rds`
- `rds/VDJ_annotation_stat.rds`
- `rds/all_seurat_before_filtering.rds`
- QC and contig summary figures.

### `downstream/2_DataClean.R`

Filters cells and adds doublet metadata.

Inputs:

- `rds/All_VDJ_annotation_filtered.rds`
- `rds/VDJ_annotation_stat.rds`
- `rds/all_seurat_before_filtering.rds`
- `config/samples.csv`

Main processing:

- Excludes `Blank-2`.
- Keeps cells with TRD or TRG evidence.
- Investigates high-frequency CDR3 contamination in cells with more than two TRD/TRG chains.
- Removes cells with `trd > 2` or `trg > 2`.
- Runs `scDblFinder`.
- Runs PCA and UMAP for doublet visualization.

Outputs:

- `rds/all_seurat_2.rds`
- Noise, chain-count, QC, and doublet figures.

### `downstream/3_cellAnnotation.R`

Clusters filtered cells, removes suspect clusters, annotates cell types, and saves the final annotated Seurat object.

Inputs:

- `rds/all_seurat_2.rds`
- `rds/VDJ_annotation_stat.rds`
- `config/samples.csv`

Main processing:

- `force_cell_annotation <- TRUE` rebuilds step-3 marker tables, full intermediate RDS files, and the final annotated Seurat RDS from `rds/all_seurat_2.rds`.
- `force_cell_annotation_plot <- TRUE` overwrites final annotation figure files from an existing final annotated RDS without rebuilding filtering outputs or marker tables.
- If `force_cell_annotation <- FALSE` and the final annotated RDS exists, the script exits early. If the final annotated RDS is missing, the workflow starts from `rds/all_seurat_2.rds`.
- Runs Seurat clustering over resolutions `0.1` to `0.3`.
- Removes likely B-cell/mixed clusters and another suspect cluster.
- Writes Seurat marker tables with `FindAllMarkers()`.
- Saves diagnostic UMAPs for removed clusters: cluster-highlight plots plus B-cell markers and a compact plasma-cell marker panel including `JCHAIN` before removing `all_seurat_2` clusters 9/11; monocyte/DC/macrophage/neutrophil marker panels and cluster-10 mast/basophil/eosinophil, myeloid/DC/macrophage, and B/plasma marker panels before removing `all_seurat_3` cluster 10.
- Adds V gene metadata.
- Annotates cells into the eight current cell types.
- Saves a CD4/CD8 T-cell marker UMAP panel on the annotated object.
- Runs cell-cycle scoring and selected hallmark/gene waterfall plots.

Outputs:

- `rds/all_seurat_3.rds`
- `rds/all_seurat_4.rds`
- `rds/all_seurat_celltype.rds`
- `table/markers.csv`
- `table/all_seurat_3_seurat_markers.csv`
- `table/all_seurat_4_seurat_markers.csv`
- `figures/3_cellAnnotation/all_seurat_2_removed_cluster_highlight.png/pdf`
- `figures/3_cellAnnotation/all_seurat_2_b_cell_marker_diagnostic.png/pdf`
- `figures/3_cellAnnotation/all_seurat_2_plasma_marker_diagnostic.png/pdf`
- `figures/3_cellAnnotation/all_seurat_3_removed_cluster_highlight.png/pdf`
- `figures/3_cellAnnotation/all_seurat_3_monocyte_markers_diagnostic.png/pdf`
- `figures/3_cellAnnotation/all_seurat_3_dc_markers_diagnostic.png/pdf`
- `figures/3_cellAnnotation/all_seurat_3_macrophage_markers_diagnostic.png/pdf`
- `figures/3_cellAnnotation/all_seurat_3_neutrophil_markers_diagnostic.png/pdf`
- `figures/3_cellAnnotation/all_seurat_3_cluster10_mast_basophil_eosinophil_markers_diagnostic.png/pdf`
- `figures/3_cellAnnotation/all_seurat_3_cluster10_myeloid_dc_macrophage_markers_diagnostic.png/pdf`
- `figures/3_cellAnnotation/all_seurat_3_cluster10_b_plasma_markers_diagnostic.png/pdf`
- `figures/3_cellAnnotation/all_seurat_4_cd4_cd8_t_cell_marker_umap.png/pdf`
- UMAP, marker, cell-type, QC, and waterfall figures.

### `downstream/4_VDJgene.R`

Adds final cell-type labels to productive VDJ annotations and summarizes V/J/C gene usage.

Inputs:

- `rds/All_VDJ_annotation_filtered.rds`
- `rds/VDJ_annotation_stat.rds`
- `rds/all_seurat_celltype.rds`

Outputs:

- `rds/all_annotation_included.rds`
- `VDJ_vln_GC`
- `VDJ_nTRDG_percell`
- `VDJ_TRD_VDJ_gene_distribution`
- `VDJ_TRG_VDJ_gene_distribution`
- `VDJ_TRG_JC_gene_distribution`
- `VDJ_combined_plots`

### `downstream/5_CDR3stat.R`

Summarizes CDR3 clone frequencies, TCR diversity/clonality metrics, CDR3 lengths, and sequence logos.

Inputs:

- `rds/all_seurat_celltype.rds`
- `rds/all_annotation_included.rds`

Main processing:

- Inspects duplicated TRD CDR3 calls per cell.
- Marks and excludes the expanded CDR3 `CALDTTFPIGDRGYTDKLIF`.
- Plots top CDR3 frequencies by sample and cell type.
- Calculates Shannon diversity and Gini-coefficient clonality by sample, condition group, and chain.
- Plots TCR diversity and clonality boxplots across condition groups with and without Naive-versus-other Wilcoxon significance annotations.
- Classifies clones as rare `<0.1%`, large `0.1-1%`, or hyperexpanded `>1%` and plots grouped stacked bars.
- Plots TRD/TRG CDR3 length distributions.
- Builds selected sequence logos with `ggseqlogo`.

Outputs:

- `top200_cdr3_samples`
- `top200_cdr3_samples_excluded`
- `top200_cdr3_cell_type_excluded`
- `table/tcr_diversity_clonality_metrics.csv`
- `table/tcr_clone_size_classes.csv`
- `TCR_diversity_clonality_boxplot`
- `TCR_diversity_clonality_boxplot_with_stats`
- `TCR_clone_size_class_stacked_bar`
- `VDJ_CDR3_length_TRD_lineplot`
- `VDJ_CDR3_length_TRG_lineplot`
- sequence-logo figures.

### `downstream/6_CDR3betweenSamples.R`

Computes overlap of top 100 CDR3s between cell types and samples.

Inputs:

- `rds/all_seurat_celltype.rds`
- `rds/all_annotation_included.rds`

Outputs:

- `shared_cdr3_between_cell_type`
- `shared_cdr3_between_samples`

### `downstream/7_CDR3pairedSankey.R`

Builds paired TRD/TRG CDR3 records and visualizes top pair sharing with Sankey/alluvial plots.

Inputs:

- `rds/all_annotation_included.rds`

Outputs:

- `rds/barcode_trgd_paired.rds`
- `sankey_TRDG_within_group`
- `sankey_TRDG_among_group`

### `downstream/8_Vd1vs2.R`

Scores MSigDB Hallmark pathways and compares current Vd1/Vd2 cell types.

Inputs:

- `rds/all_seurat_celltype.rds`

Main processing:

- Loads shared group, sample, and cell-type levels from `downstream/_cache_/plotting_shared.R`.
- Runs `GeneSetAnalysis()` with `hall50$human`.
- Saves cell-type-level Hallmark statistics.
- Saves a Hallmark heatmap and selected Vd1 versus Vd2 waterfall plots with `SeuratExtend::WaterfallPlot()`.

Outputs:

- `table/vd1_vd2_hallmark_celltype_stats.csv`
- `figures/8_Vd1vs2/vd1_vd2_hallmark_celltype_heatmap.png`
- `figures/8_Vd1vs2/vd1_vd2_hallmark_celltype_heatmap.pdf`
- `figures/8_Vd1vs2/waterfall_*`

### `downstream/9_CDR3paired.R`

Ranks shared paired TRD/TRG clones, subsets top Vd2 clones, and performs hallmark/trajectory analyses.

Inputs:

- `rds/all_seurat_celltype.rds`
- `rds/barcode_trgd_paired.rds`

Internal settings:

- `force_ranked_rds`: when set to `TRUE`, rebuilds paired-clone dependent outputs from current exact TRD+TRG ranks, including `rds/trdg_pair_rank_metadata.rds`, `table/trdg_pair_rank_metadata.csv`, `rds/all_seurat_celltype_toprank_cells.rds`, graph-test cache when available, and all figures generated by this script.

Main processing:

- Adds exact paired TRD+TRG clone rank metadata to the Seurat object; listed top clones must contain the same TRD and TRG pair and must be present in `Naive`, `ZOL`, and `PAN`, while `Blank` and `MSH2` are optional.
- Writes a TRD+TRG clone dispersion dictionary where the score is calculated from UMAP position spread and Hallmark AUCell pathway-score spread, then compares ranked clones with an overall dispersion bar plot and component heatmap.
- Draws top paired-clone stacked bar plots by condition group and by cell-type composition within each ranked clone.
- Builds Hallmark heatmaps with Hallmark pathways as rows, cell types as column splits/titles, and compact top paired TRD+TRG clone labels as columns inside each cell-type block.
- Selects top-frequency Vd2 TRD/TRG pairs until the cumulative selected cells reach about 5,000, excluding `AB3`.
- Subsets those top-pair Vd2 cells for Monocle3.
- Runs cell-cycle scoring, hallmark heatmaps, Monocle3 pseudotime, graph test, pseudotime heatmap, and hallmark enrichment.
- Skips `graph_test`-dependent outputs if `sf` cannot load because the linked system PROJ library is missing.

Outputs:

- `rds/trdg_pair_rank_metadata.rds`
- `rds/all_seurat_celltype_toprank_cells.rds`
- `rds/pr_graph_test_res_toprank.rds`
- `table/trdg_pair_rank_metadata.csv`
- `table/toprank_vd2_trdg_pairs_for_monocle3.csv`
- `table/trdg_clone_dispersion_metrics.csv`
- `figures/9_CDR3paired/trdg_clone_dispersion_score.png/pdf`
- `figures/9_CDR3paired/trdg_clone_dispersion_components.png/pdf`
- `figures/9_CDR3paired/plot_rank_top10_celltype_proportion.png/pdf`
- `figures/9_CDR3paired/toprank_umap_TRDG_pair_clone.png/pdf`
- `figures/9_CDR3paired/toprank_subset_umap_TRDG_pair_clone.png/pdf`
- `figures/9_CDR3paired/hallmark_heatmap_toprank_TRDG_pair_clones.png/pdf`
- top-pair proportion, top paired-clone UMAP, hallmark, pseudotime, and enrichment figures.

### `downstream/10_MSH2.R`

Analyzes cells carrying the MSH2-expanded TRD CDR3 sequence `CALDTTFPIGDRGYTDKLIF`.

Inputs:

- `rds/all_seurat_celltype.rds`
- `rds/all_annotation_included.rds`

Main processing:

- Marks cells with the target CDR3 sequence.
- Subsets Vd1 cells and groups them as `Naive`, `Non-Marker`, `Marker`, or `Effector`.
- Quantifies whether dominant CDR3-bearing Vd1 cells resemble effector Vd1 cells using cell-type composition, Vd1 marker-group UMAP, effector activation module score, UMAP distance to Vd1 state centroids, and marker-group differential expression.
- Runs Hallmark AUCell scoring for the Vd1 subset.
- Saves Hallmark heatmaps, marker-versus-reference delta plots, selected-gene dotplot, and effector-gene heatmap.
- Uses skip-aware plot/table writers from `downstream/_cache_/plotting_shared.R`.

Outputs:

- `table/msh2_marker_cdr3_cells.csv`
- `table/msh2_marker_group_summary.csv`
- `table/msh2_marker_effector_score_summary.csv`
- `table/msh2_marker_umap_centroid_distance.csv`
- `table/msh2_marker_umap_centroid_distance_summary.csv`
- `table/msh2_marker_de_markers.csv`
- `table/msh2_marker_hallmark_stats.csv`
- `table/msh2_marker_hallmark_delta.csv`
- `figures/10_MSH2/msh2_marker_*`


### `downstream/11_pySCENIC.sh`

Selects a balanced proper-size cell subset, exports a Seurat expression matrix in sparse-aware chunks, and runs pySCENIC GRN, cisTarget pruning, and AUCell scoring.
The default run uses a reproducible 10,000-cell subset balanced by `cell_type`, one GRN worker, and `custom_multiprocessing` ctx mode for local robustness.
The cell selector lives in `downstream/_cache_/pyscenic_select_cells.R`, the R export code lives in `downstream/_cache_/pyscenic_export.R`, and the Python pySCENIC runner lives in `downstream/_cache_/pyscenic_run.py` so it can be run with `--pyscenic-python` from a dedicated conda environment.
The `download-ref` subcommand downloads default human or mouse pySCENIC reference resources and writes a `pyscenic_reference.env` file.

Main behavior:

- Selects up to 10,000 cells by default.
- Balances sampling by `cell_type` by default.
- Keeps at least 200 cells per present group when possible.
- Writes the selected cells to `pyscenic_10k/selected_cells_10000_by_cell_type.tsv`.
- Runs pySCENIC on exactly those selected cells.

Helper:

- `downstream/_cache_/pyscenic_select_cells.R`

Inputs:

- `rds/all_seurat_celltype.rds`
- `code/downstream/envs/pyscenic_environment.yml` for the pinned pySCENIC Python environment
- `--tf-list` or `TF_LIST`
- `--ranking-db`/`--ranking-dbs` or `RANKING_DBS`
- `--motif-annotations` or `MOTIF_ANNOTATIONS`
- `--pyscenic-python` or `PYSCENIC_PYTHON` pointing to an environment containing `pyscenic`
- `--target-cells 10000` by default for the robust first run
- `--grn-num-workers 1` and `--seed 1` by default for stable GRN inference
- `--ctx-mode custom_multiprocessing` by default to avoid incompatible Dask scheduler behavior

Outputs:

- `pyscenic/expression_matrix.csv`
- `pyscenic/cell_metadata.tsv`
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

### `downstream/12_pySCENIC_visualization.R`

Reads pySCENIC AUCell loom output, aligns regulon activity to the annotated Seurat object, and visualizes regulon activity on UMAPs and grouped summaries.

Inputs:

- `rds/all_seurat_celltype.rds`
- `pyscenic_10k/auc_mtx.loom`

Outputs:

- `figures/12_pySCENIC_visualization/pyscenic_regulon_activity_umap_top.png`
- `figures/12_pySCENIC_visualization/pyscenic_regulon_activity_heatmap_by_<group>.png`
- `figures/12_pySCENIC_visualization/pyscenic_regulon_activity_dotplot_by_<group>.png`
- `figures/12_pySCENIC_visualization/pyscenic_curated_tf_heatmap.png`
- `figures/12_pySCENIC_visualization/pyscenic_curated_tf_dotplot.png`
- `figures/12_pySCENIC_visualization/pyscenic_curated_tf_umap.png`
- `figures/12_pySCENIC_visualization/pyscenic_curated_tf_feature_umap.png`
- `figures/12_pySCENIC_visualization/pyscenic_curated_tf_interpretation_table.png`
- `figures/12_pySCENIC_visualization/pyscenic_vd2_effector_regulon_heatmap.png`
- `figures/12_pySCENIC_visualization/pyscenic_vd2_effector_regulon_dotplot.png`
- `figures/12_pySCENIC_visualization/pyscenic_vd2_effector_regulon_delta_plot.png`
- `figures/12_pySCENIC_visualization/pyscenic_zol_foxp3_vd2_function_gene_heatmap.png`
- `figures/12_pySCENIC_visualization/pyscenic_zol_foxp3_vd2_function_gene_dotplot.png`
- `figures/12_pySCENIC_visualization/pyscenic_zol_foxp3_vd2_function_gene_umap.png`
- `figures/12_pySCENIC_visualization/pyscenic_zol_foxp3_vd2_foxp3_regulon_boxplot.png`
- `figures/12_pySCENIC_visualization/pyscenic_sox4_target_network_by_state.png`
- `figures/12_pySCENIC_visualization/pyscenic_sox4_target_expression_heatmap.png`
- `figures/12_pySCENIC_visualization/pyscenic_sox4_target_scatter.png`
- `table/pyscenic/pyscenic_regulon_activity_by_group.csv`
- `table/pyscenic/pyscenic_selected_regulons.txt`
- `table/pyscenic/pyscenic_curated_tf_interpretation.csv`
- `table/pyscenic/pyscenic_curated_tf_feature_status.csv`
- `table/pyscenic/pyscenic_vd2_effector_regulon_delta.csv`
- `table/pyscenic/pyscenic_zol_foxp3_vd2_function_gene_summary.csv`
- `table/pyscenic/pyscenic_zol_foxp3_vd2_function_gene_delta.csv`
- `table/pyscenic/pyscenic_sox4_target_associations.csv`
- `rds/pyscenic_visualization_result.rds`

## Review Findings

1. Runtime path assumptions are hardcoded in downstream R scripts through `setwd("/path/to/project")`. This is fine for the deployed machine but makes local reuse brittle.
2. Downstream R scripts now source `code/downstream/_cache_/plotting_shared.R` with fallbacks for `_cache_/plotting_shared.R` and `cache/plotting_shared.R`.
3. `code/config/params_reference.md` is broader than the current `params.yaml` and mentions analysis sections not consumed by the visible scripts.
4. The scripts are order-dependent and communicate mainly through RDS files in `rds/`; rerunning a middle step requires checking that upstream RDS files match the current annotation and palette conventions.
5. `code/downstream/3_cellAnnotation.R` has removed COSG and now relies on Seurat `FindAllMarkers()` for marker tables.

## Verification Performed

- Parsed all R scripts under `code/downstream/*.R` with `Rscript`.
- Ran shell syntax checks with `bash -n` on Cell Ranger helper scripts and `11_pySCENIC.sh`.
- Ran `python3 -m py_compile` on `02_generate_cellranger_configs.py`.
