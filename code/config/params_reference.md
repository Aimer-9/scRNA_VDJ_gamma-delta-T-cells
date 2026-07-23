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
