# Downstream figure inventory

This is an inventory of the figures *defined by* `code/downstream/`. It is a
review aid for pruning the workflow: figure files are written as both `.png`
and `.pdf` under `figures/<step>/` (unless a plot is skipped because required
data are unavailable). No rendered `figures/` directory is currently tracked
in this checkout.

## How to use this list

- **Yes**: retain in the default, compact figure set.
- **No**: safe first candidate to remove from the default set.
- **Review**: retain only when it supports the specific biological claim.
- Names containing `<comparison>`, `<set>`, or `<group>` represent one output
  for each configured item. Removing a family means removing/guarding the
  corresponding `save_plot()` or `save_heatmap()` call in the named script.

The natural pruning order is: diagnostic filtering panels, duplicate views of
the same result, then specialised/optional analyses (steps 9--17). Do not
delete intermediate RDS/CSV outputs merely because a plot is removed.

## Core processing, QC, and annotation

| Step and output directory | Figure basename(s) | What it shows | Keep? | Suggested role |
| --- | --- | --- | --- | --- |
| 1 `1_ReadData` | `sample_productive`; `VDJ_contig_percentage_per_sample_barplot`; `QC_metrics_violin_plot_before_filtering` | Productive VDJ counts, contig composition, pre-filter RNA QC | Yes | Keep the contig and QC panels; `sample_productive` optional |
| 2 `2_DataClean` | `QC_metrics_violin_plot_after_filter1`; `noise_cdr3_marked_plot`; `noise_cdr3_umi_plot`; `boxplot_umi_cdr3_top`; `annotation_plot_nTRDG`; `qc_violin_plot_filtered`; `umap_doublet` | Sequential QC, suspected contaminant CDR3, chain-count, and doublet checks | Review | Keep one pre/post QC panel and `umap_doublet`; the three noise-CDR3 plots and `annotation_plot_nTRDG` are diagnostic candidates to remove |
| 3 `3_cellAnnotation` | `<all_seurat_2/3>_umap_group_sample`; `<all_seurat_2/3>_umap_cluster_res_0_[1,2,3]`; removal-highlights; `*_diagnostic` marker UMAPs | Intermediate clustering and B/plasma/myeloid contamination-removal evidence | No | Diagnostic: remove after the filtering decision is documented |
| 3 `3_cellAnnotation` | `all_seurat_4_umap_vd12_vg49`; `all_seurat_4_umap_cellcycle`; `all_seurat_4_umap_celltype`; `all_seurat_4_cd4_cd8_t_cell_marker_umap`; `celltype_distr`; `qc_violin_plot_celltype`; `dotplot_celltype` | Final annotation, marker support, abundance, and per-type QC | Yes | Keep cell-type UMAP, marker dotplot, and distribution; cell-cycle, V-gene, CD4/CD8, and QC panels are optional |
| 4 `4_VDJgene` | `VDJ_vln_GC`; `VDJ_nTRDG_percell`; `VDJ_TRD_VDJ_gene_distribution`; `VDJ_TRG_VDJ_gene_distribution`; `VDJ_TRG_JC_gene_distribution`; `VDJ_combined_plots` | VDJ gene and chain usage | Review | Prefer `VDJ_combined_plots`; the four panels embedded in it are duplicate candidates |

## CDR3 repertoire and clone analyses

| Step and output directory | Figure basename(s) | What it shows | Keep? | Suggested role |
| --- | --- | --- | --- | --- |
| 5 `5_CDR3stat` | `top200_cdr3_samples`; `umap_mark`; `cluster_distr_mark`; `top200_cdr3_samples_excluded`; `top200_cdr3_cell_type_excluded` | Dominant clones and the excluded extreme MSH2 clone | Review | Keep one post-exclusion frequency view; the before/marking panels are diagnostic |
| 5 `5_CDR3stat` | `TCR_diversity_clonality_boxplot`; `TCR_diversity_clonality_boxplot_with_stats`; `TCR_clone_size_class_stacked_bar` | Diversity/clonality and clone-size classes by condition | Yes | Keep the stats version *or* the clean version, not both; clone-size bar optional |
| 5 `5_CDR3stat` | `VDJ_CDR3_length_TRG_lineplot`; `VDJ_CDR3_length_TRD_lineplot`; `<cdr3_logo_set prefix>` | CDR3 lengths and sequence logos at selected peaks | No | Optional/supplementary |
| 6 `6_CDR3betweenSamples` | `shared_cdr3_between_cell_type`; `shared_cdr3_between_samples` | Pair of triangular heatmaps: shared top-100 TRD and TRG CDR3s | Review | Keep only the comparison level relevant to the story |
| 7 `7_CDR3pairedSankey` | `sankey_TRDG_within_group`; `sankey_TRDG_among_group` | TRD--TRG pairing within/across conditions | No | Optional; Sankey panels can be visually dense |
| 9 `9_CDR3paired` | `trdg_clone_dispersion_score`; `trdg_clone_dispersion_components`; `plot_rank_top10_proportion`; `plot_rank_top10_celltype_proportion`; `toprank_umap_TRDG_pair_clone` | Exact paired-clone dispersion, burden, cell-type composition, and location | Review | Keep a burden/composition panel and UMAP if paired clones are a key result; the two dispersion plots are optional |
| 9 `9_CDR3paired` | `toprank_umap_celltype`; `toprank_subset_umap_TRDG_pair_clone`; `hallmark_heatmap_top_vd2_pairs`; `hallmark_heatmap_toprank_TRDG_pair_clones`; `hallmark_heatmap_all_vd2`; `hallmark_heatmap_all_vd1` | Selected top-pair Vd2 subset and Hallmark activity | Review | Several overlapping heatmaps: retain the clone-split heatmap plus at most one broad Vd1/Vd2 view |
| 9 `9_CDR3paired` | `plot_umap_pseudotime`; `plot_umap_cellcycle`; `pseudotime_heatmap`; `hallmark_enrichment_cluster_3`; `all_seurat_celltype_toprank_heatmap_hallmark` | Top-pair Vd2 trajectory and functional programs | No | Optional specialised analysis; keep only if pseudotime is a result |

## State and condition comparisons

| Step and output directory | Figure basename(s) | What it shows | Keep? | Suggested role |
| --- | --- | --- | --- | --- |
| 8 `8_Vd1vs2` | `waterfall_<comparison>`; `waterfall_vd1_vd2_hallmark_comparisons`; `vd1_vd2_hallmark_celltype_heatmap`; `vd1_vd2_pairwise_marker_dotplot`; `vd1_vd2_curated_marker_dotplot`; `vd1_vd2_TRDV1_TRDV2_expression_violin`; `vd1_vd2_module_score_violin`; `zol_pan_effector_vd2_module_score_violin` | Hallmark, marker, TRDV expression, and module-score views of Vd1/Vd2 states | Yes | Consolidated Vd1/Vd2 step |
| 10 `10_MSH2` | `msh2_marker_cdr3g_to_caldttfpigdrgytdklif_top10_sankey`; `msh2_marker_top_cdr3g_to_cdr3d_top10_sankey`; `msh2_marker_cdr3g_cdr3d_pairing_combined_sankey` | MSH2-associated CDR3/TRD--TRG pairing | No | Optional specialised repertoire evidence; combined Sankey supersedes its two components |
| 10 `10_MSH2` | `msh2_marker_cdr3_umap`; `msh2_marker_celltype_composition`; `msh2_marker_group_vd1_umap`; `msh2_marker_effector_activation_score`; `msh2_marker_umap_centroid_distance` | Location, composition, activation, and spatial separation of marker Vd1 cells | Review | Keep group UMAP and activation score if MSH2 is central; centroid distance is optional |
| 10 `10_MSH2` | `msh2_marker_hallmark_heatmap`; `msh2_marker_hallmark_heatmap_selected`; `msh2_marker_vs_nonmarker_hallmark_delta`; `msh2_marker_vs_effector_hallmark_delta`; `msh2_marker_vd1_gene_dotplot`; `msh2_marker_effector_gene_heatmap` | MSH2 marker-group pathways and genes | Review | Retain one pathway view and one gene view; full and selected heatmaps are redundant |
| 12 `12_pySCENIC_visualization` | `pyscenic_regulon_activity_heatmap_by_<group>`; `pyscenic_regulon_activity_dotplot_by_<group>`; `pyscenic_curated_tf_heatmap` | Global and curated transcription-factor regulon activity | No | Optional regulatory analysis; choose heatmap or dotplot for the global view |
| 12 `12_pySCENIC_visualization` | `pyscenic_vd2_effector_regulon_dotplot`; `pyscenic_zol_foxp3_vd2_function_gene_heatmap`; `pyscenic_zol_foxp3_vd2_function_gene_dotplot`; `pyscenic_sox4_target_expression_heatmap` | Effector-Vd2 and ZOL FOXP3+/SOX4 regulatory interpretation | No | Optional specialised analysis; function-gene heatmap/dotplot are alternate views |

## Fast, low-risk pruning set

If the goal is a tighter main figure set without changing any analysis, first
hide or remove the plotting calls for these groups:

### Exact removal map (one figure per row)

Each row below is an individual figure currently marked **No**. To stop making
it, comment out or delete the matching `save_plot(...)` / `save_heatmap(...)`
call in the listed script. This leaves the analysis data and tables intact.

| Figure basename | Source code to remove/comment | Keep? |
| --- | --- | --- |
| `noise_cdr3_marked_plot` | `code/downstream/2_DataClean.R`: `save_plot(noise_marked_plot, "noise_cdr3_marked_plot", ...)` | No |
| `noise_cdr3_umi_plot` | `code/downstream/2_DataClean.R`: `save_plot(noise_cdr3_umi_plot, "noise_cdr3_umi_plot", ...)` | No |
| `boxplot_umi_cdr3_top` | `code/downstream/2_DataClean.R`: `save_plot(boxplot_umi_cdr3_top, "boxplot_umi_cdr3_top", ...)` | No |
| `annotation_plot_nTRDG` | `code/downstream/2_DataClean.R`: `save_plot(annotation_plot_nTRDG, "annotation_plot_nTRDG", ...)` | No |
| `<all_seurat_2/3>_umap_group_sample` | `code/downstream/3_cellAnnotation.R`: call to `save_group_sample_umaps(...)` | No |
| `<all_seurat_2/3>_umap_cluster_res_0_[1,2,3]` | `code/downstream/3_cellAnnotation.R`: call to `save_cluster_umaps(...)` | No |
| `<all_seurat_2/3>_removed_cluster_highlight` | `code/downstream/3_cellAnnotation.R`: `save_cluster_removal_plot(...)` calls | No |
| `all_seurat_2_b_cell_marker_diagnostic` | `code/downstream/3_cellAnnotation.R`: `save_marker_feature_plot(..., "all_seurat_2_b_cell_marker_diagnostic", ...)` | No |
| `all_seurat_2_plasma_marker_diagnostic` | `code/downstream/3_cellAnnotation.R`: `save_marker_feature_plot(..., "all_seurat_2_plasma_marker_diagnostic", ...)` | No |
| `all_seurat_3_*_diagnostic` | `code/downstream/3_cellAnnotation.R`: both `purrr::iwalk(...)` diagnostic blocks | No |
| `VDJ_CDR3_length_TRG_lineplot` | `code/downstream/5_CDR3stat.R`: `save_plot(..., "VDJ_CDR3_length_TRG_lineplot", ...)` | No |
| `VDJ_CDR3_length_TRD_lineplot` | `code/downstream/5_CDR3stat.R`: `save_plot(..., "VDJ_CDR3_length_TRD_lineplot", ...)` | No |
| each `<cdr3_logo_set prefix>` | `code/downstream/5_CDR3stat.R`: the `for (logo_set in cdr3_logo_sets)` block | No |
| `sankey_TRDG_within_group` | `code/downstream/7_CDR3pairedSankey.R`: `save_plot(..., "sankey_TRDG_within_group", ...)` | No |
| `sankey_TRDG_among_group` | `code/downstream/7_CDR3pairedSankey.R`: `save_plot(..., "sankey_TRDG_among_group", ...)` | No |
| `plot_umap_pseudotime` | `code/downstream/9_CDR3paired.R`: `save_plot(..., "plot_umap_pseudotime", ...)` | No |
| `plot_umap_cellcycle` | `code/downstream/9_CDR3paired.R`: `save_plot(..., "plot_umap_cellcycle", ...)` | No |
| `pseudotime_heatmap` | `code/downstream/9_CDR3paired.R`: `save_heatmap(..., "pseudotime_heatmap", ...)` | No |
| `hallmark_enrichment_cluster_3` | `code/downstream/9_CDR3paired.R`: `save_plot(..., "hallmark_enrichment_cluster_3", ...)` | No |
| `all_seurat_celltype_toprank_heatmap_hallmark` | `code/downstream/9_CDR3paired.R`: `save_plot(..., "all_seurat_celltype_toprank_heatmap_hallmark", ...)` | No |
| `msh2_marker_cdr3g_to_caldttfpigdrgytdklif_top10_sankey` | `code/downstream/10_MSH2.R`: matching `save_plot(...)` call | No |
| `msh2_marker_top_cdr3g_to_cdr3d_top10_sankey` | `code/downstream/10_MSH2.R`: matching `save_plot(...)` call | No |
| `msh2_marker_cdr3g_cdr3d_pairing_combined_sankey` | `code/downstream/10_MSH2.R`: matching `save_plot(...)` call | No |
| `pyscenic_regulon_activity_heatmap_by_<group>` | `code/downstream/12_pySCENIC_visualization.R`: `save_heatmap(...)` near the global regulon section | No |
| `pyscenic_regulon_activity_dotplot_by_<group>` | `code/downstream/12_pySCENIC_visualization.R`: `save_ggplot(...)` near the global regulon section | No |
| `pyscenic_curated_tf_heatmap` | `code/downstream/12_pySCENIC_visualization.R`: matching `save_heatmap(...)` call | No |
| `pyscenic_vd2_effector_regulon_dotplot` | `code/downstream/12_pySCENIC_visualization.R`: matching `save_ggplot(...)` call | No |
| `pyscenic_zol_foxp3_vd2_function_gene_heatmap` | `code/downstream/12_pySCENIC_visualization.R`: matching `save_heatmap(...)` call | No |
| `pyscenic_zol_foxp3_vd2_function_gene_dotplot` | `code/downstream/12_pySCENIC_visualization.R`: matching `save_ggplot(...)` call | No |
| `pyscenic_sox4_target_expression_heatmap` | `code/downstream/12_pySCENIC_visualization.R`: matching `save_heatmap(...)` call | No |

1. Step 2 noise-CDR3 diagnostics and step 3 intermediate contamination-removal
   diagnostics.
2. Duplicate component panels already contained in `VDJ_combined_plots`, the
   combined Vd1/Vd2 waterfall, or a combined Sankey.
3. One of each alternate visual encoding: full vs selected Hallmark heatmap,
   heatmap vs dotplot, and clean vs statistics-labelled diversity boxplot.
4. Entire optional modules: pySCENIC (12), and repertoire panels if they are not part of the biological claim.

Steps 0 and 11 do not define figures; step 11 runs pySCENIC and step 12 renders
its figures.
