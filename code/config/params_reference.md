# `config/params.yaml` Parameter Reference

This document explains every key currently used in `config/params.yaml`.

## Top-level

| Key | Type | Meaning |
|---|---|---|
| `outdir` | string | Root output directory for the full pipeline. Subfolders like `cellranger_configs`, `cellranger_output`, and `r_output` are created under this path. |

## `cellranger`

| Key | Type | Meaning |
|---|---|---|
| `cellranger.reference_gex` | string | Path to Cell Ranger gene-expression reference (GRCh38 reference package). |
| `cellranger.reference_vdj` | string | Path to Cell Ranger V(D)J reference package. |
| `cellranger.cellranger_path` | string | Path to the `cellranger` executable. |
| `cellranger.fastqc_path` | string | Path to the `fastqc` executable used in FASTQ QC. |
| `cellranger.fastqc_threads` | integer | CPU threads passed to FastQC when `--threads` is not provided. |
| `cellranger.inner_enrichment_primers` | string | Path to primer list used for relevant V(D)J enrichment workflows. |
| `cellranger.localcores` | integer | CPU cores given to Cell Ranger runs. |
| `cellranger.localmem` | integer | Memory (GB) given to Cell Ranger runs. |

## `qc`

| Key | Type | Meaning |
|---|---|---|
| `qc.soupx` | boolean | Enable ambient RNA correction with SoupX in QC step. |
| `qc.min_genes` | integer | Minimum detected genes per cell to keep. |
| `qc.max_genes` | integer | Maximum detected genes per cell to keep (high values can indicate doublets). |
| `qc.max_mito_pct` | number | Maximum mitochondrial percentage allowed per cell. |
| `qc.min_counts` | integer | Minimum UMI counts per cell to keep. |
| `qc.max_counts` | integer | Maximum UMI counts per cell to keep. |
| `qc.doublet_method` | string | Doublet detection backend (for example `scDblFinder`). |
| `qc.doublet_filter` | boolean | Whether to remove predicted doublets. |
| `qc.filter_method` | string | QC filter policy. Typical values: `standard`, `vdj_detected_if_exists`. |

## `normalization`

| Key | Type | Meaning |
|---|---|---|
| `normalization.method` | string | Final-mode normalization method. Supported by current pipeline logic: `LogNormalization`, `SCTransform`. |
| `normalization.parallel_workers` | integer | Worker count for Step-02 group-level parallelism (non-Windows fork path). |

### `normalization.log_normalization`

| Key | Type | Meaning |
|---|---|---|
| `normalization.log_normalization.method` | string | Log-normalization method label, typically `LogNormalize`. |

### `normalization.sctransform`

| Key | Type | Meaning |
|---|---|---|
| `normalization.sctransform.vst_flavor` | string | SCTransform flavor (commonly `v2`). |
| `normalization.sctransform.variable_features` | integer | Number of variable genes for SCTransform in final run. |
| `normalization.sctransform.vars_to_regress` | string list | Metadata covariates to regress during SCTransform (for example `percent.mt`). |

### `normalization.scvi`

| Key | Type | Meaning |
|---|---|---|
| `normalization.scvi.conda_env` | string | Conda environment name containing scvi-tools dependencies. |
| `normalization.scvi.python_bin` | string | Optional explicit Python path. Empty string means auto-resolve from environment. |
| `normalization.scvi.n_latent` | integer | scVI latent dimensionality. |
| `normalization.scvi.max_epochs` | integer | Maximum training epochs. |
| `normalization.scvi.batch_size` | integer | Mini-batch size for training. |
| `normalization.scvi.train_size` | number | Fraction of cells used for training split. |
| `normalization.scvi.validation_size` | number | Fraction of cells used for validation split. |
| `normalization.scvi.early_stopping` | boolean | Enable early stopping during training. |
| `normalization.scvi.early_stopping_monitor` | string | Metric used for early stopping (for example `elbo_validation`). |
| `normalization.scvi.early_stopping_patience` | integer | Number of evaluation windows to wait before stopping. |
| `normalization.scvi.check_val_every_n_epoch` | integer | Validation interval in epochs. |
| `normalization.scvi.target_sum` | integer | Target library size used in preprocessing/normalization steps. |
| `normalization.scvi.span` | number | Smoothing span used in relevant preprocessing stages. |
| `normalization.scvi.n_top_genes` | integer | Number of top genes used in scVI preprocessing feature selection. |

## `integration`

| Key | Type | Meaning |
|---|---|---|
| `integration.method` | string | Final-mode integration method. Supported by current logic: `Harmony`, `CCA`, `NoIntegration`. |
| `integration.harmony_theta` | number | Harmony diversity penalty parameter (`theta`). |
| `integration.dims` | integer | Number of dimensions used for neighbors/UMAP/clustering in final run. |
| `integration.batch_column` | string | Metadata column treated as batch variable (for example `sample_id`). |

## `clustering`

| Key | Type | Meaning |
|---|---|---|
| `clustering.resolution` | number (scalar) | Final clustering resolution. Must be a single positive numeric value. |
| `clustering.algorithm` | integer | Seurat `FindClusters` algorithm code (for example `4` for Leiden in current setup). |
| `clustering.min_cells_per_cluster_for_marker` | integer | Minimum cluster size for marker DE inclusion. Small clusters under this threshold are skipped in marker testing. |

## `annotation`

| Key | Type | Meaning |
|---|---|---|
| `annotation.method` | string | Legacy/compatibility annotation method key. Current Step-03 is CSV-driven; Step-02 marker sweep still uses marker-related keys below. |
| `annotation.azimuth_reference` | string | Legacy Azimuth reference key kept for compatibility. |
| `annotation.singler_reference` | string | Legacy SingleR reference key kept for compatibility. |
| `annotation.annotation_level` | string | Preferred annotation level label (for example `celltype.l2`) for compatibility/reporting conventions. |
| `annotation.marker_top_n` | integer | Top marker genes retained per cluster in marker outputs. |
| `annotation.geneset_top_n` | integer | Number of top-scoring gene sets shown per cluster in geneset dotplots. |
| `annotation.marker_assay` | string | Marker DE assay selection: `auto`, `RNA`, or `SCT`. |
| `annotation.enable_cosg` | boolean | Enable COSG marker ranking path. |
| `annotation.enable_findmarkers` | boolean | Enable Seurat `FindAllMarkers(wilcox)` path (Presto-accelerated if available). |
| `annotation.resolution_eval_max_cells` | integer | Max cells sampled for resolution/silhouette evaluation to control runtime. |

## `vdj`

| Key | Type | Meaning |
|---|---|---|
| `vdj.clone_call_tcr` | string | Clonotype definition mode for TCR (for example `strict`, `aa`, `nt`, `gene`). |
| `vdj.clone_call_bcr` | string | Clonotype definition mode for BCR. |
| `vdj.filter_nonproductive` | boolean | Remove nonproductive chains from downstream analyses. |
| `vdj.diversity_metrics` | string list | Diversity metrics to compute (for example `shannon`, `inv.simpson`, `chao1`). |
| `vdj.overlap_method` | string | Overlap similarity method (for example `morisita`). |

## `differential`

| Key | Type | Meaning |
|---|---|---|
| `differential.de_method` | string | Differential expression strategy (for example `pseudobulk_deseq2`). |
| `differential.min_cells_per_group` | integer | Minimum cells required per group for DE/composition testing. |
| `differential.min_samples_per_group` | integer | Minimum biological samples per group for pseudobulk-style tests. |
| `differential.lfc_threshold` | number | Log fold-change threshold used in differential filtering. |
| `differential.padj_threshold` | number | Adjusted p-value cutoff used in differential filtering. |
| `differential.comparison_column` | string | Metadata column used to define groups to compare (commonly `group`). |
| `differential.composition_method` | string | Cell composition testing method (for example `propeller`). |

## Notes

- Keep exploration sweeps in `config/analysis.yaml`.
- Keep `config/params.yaml` for final/default scalar settings.
- Current contract requires `clustering.resolution` as a single scalar value (not a list).
