# 07 — Orthogonal molecular support

`multiomics_de_support/` holds the numbered chain that assembles differential-expression
support across sources: model overlap, per-contrast single-cell support, gene-level
direction, cross-cell-type support, up/down summaries, the linkage diagram, the CPTAC
limma re-analysis and proteome-wide differential expression, the CosMx highlight panels
and the CosMx supplementary table. `pseudobulk_de/` carries the pseudobulk aggregation and
the differential-expression contrast drivers. `spatial_geomx_cosmx/` holds the GeoMx
extraction and spatial differential expression, and the CosMx per-cell, composition,
gene-level paired-test and panel code.

## Position in the pipeline

```
04_association_testing
        |
        v
   >>> 07_orthogonal_support <<<
        |
        v
09_functional_screens  12_supplementary_tables
```

The arrows above are the logical order of the published workflow, taken from the Methods (`docs/workflow.md`). The file-level links in *From other sections* and *Consumed by* below are the ones actually observed in the deposited code.

## Scripts, in run order (35 files)

Ordered by sub-directory, then by the numeric prefix the script carries in its original working directory. Descriptions are the script's own header line; blank where the script has no header.

| # | File | Lang | What it does | Declared display item |
|---|---|---|---|---|
| 1 | `multiomics_de_support/00_sctwas_model_overlap_triangle_13models.R` | R | 09_sctwas_model_overlap_triangle_13models.R Gene-set overlap between 13 scTWAS cell-type models: 7 normal-colon (NL) models + 6 Polyps models. |  |
| 2 | `multiomics_de_support/01_sctwas_scRNAseq_DEsupport_C_E_L.R` | R | 07_sctwas_scRNAseq_DEsupport_C_E_L.R --------------------------------------- scRNA-seq differential-expression support for the scTWAS gene x cell-type. |  |
| 3 | `multiomics_de_support/01a_sctwas_gene_level_direction_and_support.R` | R | 10_sctwas_gene_level_direction_and_support.R -------------------------------- Collapse the 546 scTWAS gene x cell-type pairs to the 347 unique genes, . |  |
| 4 | `multiomics_de_support/02_sctwas_ASC_SSC_crossceltype_support.R` | R | Cross-cell-type scRNA-seq DE support for the scTWAS genes nominated in the lesion-specific ADSER models, ASC and SSC. |  |
| 5 | `multiomics_de_support/03_sctwas_DEsupport_summary_UpDown.R` | R | Summary of how many scTWAS genes each cell-type model has supported by scRNA-seq differential expression, split Up / Down, for the eight C / E / L con. |  |
| 6 | `multiomics_de_support/04_sctwas_ASC_SSC_linkage_diagram.R` | R | Linkage diagram of how the lesion-specific scTWAS genes connect the normal epithelial cell types to the polyp compartment and to carcinoma. |  |
| 7 | `multiomics_de_support/06_sctwas_multisource_DEsupport_diverging.R` | R | 05_sctwas_multisource_DEsupport_diverging.R -------------------------------- Panel C: how many of the 347 scTWAS genes are supported, and in which dir. |  |
| 8 | `multiomics_de_support/06a_cptac_limma_reanalysis.R` | R | 06_cptac_limma_reanalysis.R ------------------------------------------------- Alternative differential-abundance analysis for the CPTAC-COAD protein a. |  |
| 9 | `multiomics_de_support/06b_cptac_proteome_wide_DE.py` | Python | 07_cptac_proteome_wide_DE.py |  |
| 10 | `multiomics_de_support/08_prep_cosmx_highlight_panels.py` | Python |  | [internal report numbering: Figure 16; Figure 18] |
| 11 | `multiomics_de_support/09_fig_NOR_CORE_highlight_4panel.R` | R | One row, four columns: the CosMx slide's cell composition, then three scTWAS risk genes up in the normal -> carcinoma transition, each annotated with . | [internal report numbering: Figure 16; Figure 18] |
| 12 | `multiomics_de_support/09_sctwas_model_overlap_triangle.R` | R | Gene-set overlap between the 16 scTWAS cell-type models (7 normal colon + 9 polyp). |  |
| 13 | `multiomics_de_support/11_sctwas_DEsupport_UpDown_bars.R` | R | Redesign of Fig_scTWAS_DEsupport_summary_UpDown as a grouped bar chart. |  |
| 14 | `multiomics_de_support/12_sctwas_multisource_bubble_heatmap.R` | R | Redesign of Fig_scTWAS_multisource_DEsupport_diverging as two bubble heatmaps, split by the direction of the genetic effect. |  |
| 15 | `multiomics_de_support/13_figure5_legend_and_review.py` | Python |  | Figure 5 |
| 16 | `multiomics_de_support/14_build_supptable_S20_cosmx.py` | Python |  | Supplementary Table 20; Figure 5d [internal report numbering: Figure 30] |
| 17 | `pseudobulk_de/aggregate_pseudobulk.py` | Python | Rebuild specimen x cell-type pseudobulk count matrices from the four VUMC/HTAN h5ad files. |  |
| 18 | `pseudobulk_de/refit_L_expanded.py` | Python | Refit L1-L2 with the ASC / SSC arms expanded to all four files (ASC n=75, SSC n=46). |  |
| 19 | `pseudobulk_de/run_de.py` | Python | Specimen-level pseudobulk DE: early stage E1-E3, late stage L1-L6. |  |
| 20 | `pseudobulk_de/run_new_contrasts.py` | Python | New contrast sets requested on top of res/new_DE: |  |
| 21 | `spatial_geomx_cosmx/18_spatial_discover.py` | Python |  |  |
| 22 | `spatial_geomx_cosmx/19_geomx_extract.py` | Python |  |  |
| 23 | `spatial_geomx_cosmx/20_spatial_de.py` | Python |  |  |
| 24 | `spatial_geomx_cosmx/21_cosmx_composition.py` | Python |  |  |
| 25 | `spatial_geomx_cosmx/22_spatial_figures.py` | Python | 22 - Spatial figures (12-14). | [internal report numbering: Figure 12] |
| 26 | `spatial_geomx_cosmx/23_validated_set.py` | Python |  |  |
| 27 | `spatial_geomx_cosmx/24_rds_reader.py` | Python |  |  |
| 28 | `spatial_geomx_cosmx/25_cosmx_percell.py` | Python |  |  |
| 29 | `spatial_geomx_cosmx/26_cosmx_gene_tests.py` | Python |  |  |
| 30 | `spatial_geomx_cosmx/27_validated_priority.py` | Python |  |  |
| 31 | `spatial_geomx_cosmx/28_validated_figures.py` | Python | 28 - Figures 15-17: the validated gene set and the CosMx single-cell evidence. | [internal report numbering: Figure 15] |
| 32 | `spatial_geomx_cosmx/29_cosmx_plot_extract.py` | Python |  |  |
| 33 | `spatial_geomx_cosmx/30_cosmx_gene_panels.py` | Python |  | [internal report numbering: Figure 18; Figure 16b-e] |
| 34 | `spatial_geomx_cosmx/31_summary_statistics.py` | Python |  |  |
| 35 | `spatial_geomx_cosmx/47_cosmx_requested_genes.py` | Python |  | [internal report numbering: Figure 28; Figure 18] |

## Inputs

**From other sections of this deposit**

| File | Produced by |
|---|---|
| `risk_genes_harmonized.csv` | `08_druggability_and_connectivity/gene_sets_and_resources/01_parse_risk_genes.py` |
| `cosmx_selected_genes.csv` | `09_functional_screens/48_screened_and_selected.py` |
| `cosmx_selected_genes_stats.csv` | `09_functional_screens/48_screened_and_selected.py` |

**Not produced anywhere in this deposit** (26) — cohort data, public resources, or output of a step that is not deposited:

- `CPTAC_COAD_protein_matrix.csv`
- `CPTAC_COAD_psm_counts.csv`
- `CPTAC_COAD_samples.csv`
- `TCGA_CPTAC_support.csv`
- `cosmx_cells.tsv.gz`
- `cosmx_composition.tsv`
- `cosmx_composition_epi.tsv`
- `cosmx_composition_niche.tsv`
- `cosmx_percell_summary.json`
- `de_C_E_L_contrasts_ASC75_SSC46.csv`
- `de_C_E_L_contrasts_ASC75_SSC46.parquet`
- `de_all_contrasts.parquet`
- `de_exp1.tsv.gz`
- `de_exp2_epi.tsv.gz`
- `de_exp2_str.tsv.gz`
- `exp2_annot_Metadata.tsv`
- `exp2_annot_SegmentProperties.tsv`
- `gene_druggability.csv`
- `geomx_exp2_cms.tsv`
- `geomx_expr_exp2.tsv.gz`
- `integrated_gene_priority.csv`
- `meta_sel.tsv.gz`
- `pseudobulk_units.csv`
- `risk_gene_driver_proximity.csv`
- `scTWAS_ADSER_NL_collected_results.csv`
- `validated_summary.json`

## Outputs

| File | Written by | Consumed by |
|---|---|---|
| `ASC_SSC_contrast_legend.csv` | `multiomics_de_support/02_sctwas_ASC_SSC_crossceltype_support.R` | — (terminal) |
| `ASC_SSC_genes_crossceltype_support_counts.csv` | `multiomics_de_support/02_sctwas_ASC_SSC_crossceltype_support.R` | — (terminal) |
| `ASC_SSC_genes_crossceltype_support_long.csv` | `multiomics_de_support/02_sctwas_ASC_SSC_crossceltype_support.R` | `13_figure5_legend_and_review.py` |
| `ASC_SSC_genes_crossceltype_support_matrix.csv` | `multiomics_de_support/02_sctwas_ASC_SSC_crossceltype_support.R` | — (terminal) |
| `ASC_SSC_genes_crossceltype_support_matrix_flat.csv` | `multiomics_de_support/02_sctwas_ASC_SSC_crossceltype_support.R` | — (terminal) |
| `ASC_SSC_genes_crossceltype_support_summary.csv` | `multiomics_de_support/02_sctwas_ASC_SSC_crossceltype_support.R` | — (terminal) |
| `CPTAC_COAD_protein_DE_limma.csv` | `multiomics_de_support/06a_cptac_limma_reanalysis.R` | — (terminal) |
| `CPTAC_COAD_proteome_wide_DE_tumor_vs_normal.csv` | `multiomics_de_support/06b_cptac_proteome_wide_DE.py` | — (terminal) |
| `Fig_NOR_CORE_highlight_4panel_sourcedata.csv` | `multiomics_de_support/09_fig_NOR_CORE_highlight_4panel.R` | — (terminal) |
| `Fig_scTWAS_DEsupport_UpDown_bars_sourcedata.csv` | `multiomics_de_support/11_sctwas_DEsupport_UpDown_bars.R` | — (terminal) |
| `Fig_scTWAS_multisource_bubble_sourcedata.csv` | `multiomics_de_support/12_sctwas_multisource_bubble_heatmap.R` | — (terminal) |
| `Figure5A_panelA_section_counts.csv` | `multiomics_de_support/13_figure5_legend_and_review.py` | — (terminal) |
| `Figure5A_panelA_support_counts.csv` | `multiomics_de_support/13_figure5_legend_and_review.py` | — (terminal) |
| `SupplementaryTable_S20_CosMx_singlecell.csv` | `multiomics_de_support/14_build_supptable_S20_cosmx.py` | — (terminal) |
| `_sourcedata.csv` | `multiomics_de_support/00_sctwas_model_overlap_triangle_13models.R` | — (terminal) |
| `cosmx_celltype_proportions.csv` | `spatial_geomx_cosmx/21_cosmx_composition.py` | `22_spatial_figures.py` |
| `cosmx_celltype_specificity.csv` | `spatial_geomx_cosmx/26_cosmx_gene_tests.py` | — (terminal) |
| `cosmx_celltype_trend.csv` | `spatial_geomx_cosmx/21_cosmx_composition.py` | `22_spatial_figures.py`, `10_supplementary_tables.py`, `11_report_data.py` |
| `cosmx_crc_vs_normal_epi.csv` | `spatial_geomx_cosmx/26_cosmx_gene_tests.py` | `13_figure5_legend_and_review.py`, `14_build_supptable_S20_cosmx.py`, `27_validated_priority.py` (+7) |
| `cosmx_epithelium_trend.csv` | `spatial_geomx_cosmx/21_cosmx_composition.py` | — (terminal) |
| `cosmx_gene_stage_trend.csv` | `spatial_geomx_cosmx/26_cosmx_gene_tests.py` | `10_supplementary_tables.py` |
| `cosmx_geomx_replication.csv` | `spatial_geomx_cosmx/26_cosmx_gene_tests.py` | `08_prep_cosmx_highlight_panels.py`, `28_validated_figures.py`, `48_screened_and_selected.py` (+1) |
| `cosmx_niche_trend.csv` | `spatial_geomx_cosmx/21_cosmx_composition.py` | `10_supplementary_tables.py` |
| `cosmx_panelA_composition.csv.gz` | `multiomics_de_support/08_prep_cosmx_highlight_panels.py` | — (terminal) |
| `cosmx_panel_expr.csv.gz` | `multiomics_de_support/08_prep_cosmx_highlight_panels.py` | — (terminal) |
| `cosmx_requested_genes.csv` | `spatial_geomx_cosmx/47_cosmx_requested_genes.py` | `10_supplementary_tables.py`, `11_report_data.py` |
| `cx_celltype_cells.csv` | `spatial_geomx_cosmx/25_cosmx_percell.py` | — (terminal) |
| `cx_pseudobulk.csv.gz` | `spatial_geomx_cosmx/25_cosmx_percell.py` | `14_build_supptable_S20_cosmx.py`, `26_cosmx_gene_tests.py`, `28_validated_figures.py` (+2) |
| `de_senescence_low_grade_vs_normal.csv` | `multiomics_de_support/01a_sctwas_gene_level_direction_and_support.R` | — (terminal) |
| `de_stemness_high_vs_low.csv` | `multiomics_de_support/01a_sctwas_gene_level_direction_and_support.R` | — (terminal) |
| `exp2_cms.tsv` | `spatial_geomx_cosmx/19_geomx_extract.py` | — (terminal) |
| `found.json` | `spatial_geomx_cosmx/18_spatial_discover.py` | — (terminal) |
| `geomx_DE_compact.tsv` | `spatial_geomx_cosmx/19_geomx_extract.py` | — (terminal) |
| `panel_stats.csv` | `multiomics_de_support/08_prep_cosmx_highlight_panels.py` | — (terminal) |
| `risk_gene_all_comparisons.csv` | `spatial_geomx_cosmx/31_summary_statistics.py` | `10_supplementary_tables.py` |
| `risk_gene_summary_statistics.csv` | `spatial_geomx_cosmx/31_summary_statistics.py` | `10_supplementary_tables.py`, `11_report_data.py` |
| `scTWAS_347genes_direction_summary.csv` | `multiomics_de_support/01a_sctwas_gene_level_direction_and_support.R` | — (terminal) |
| `scTWAS_347genes_model_Zscore_matrix.csv` | `multiomics_de_support/01a_sctwas_gene_level_direction_and_support.R` | — (terminal) |
| `scTWAS_347genes_multisource_DE_support.csv` | `multiomics_de_support/01a_sctwas_gene_level_direction_and_support.R` | `08_prep_cosmx_highlight_panels.py`, `12_sctwas_multisource_bubble_heatmap.R`, `13_figure5_legend_and_review.py` (+1) |
| `scTWAS_347genes_support_counts.csv` | `multiomics_de_support/01a_sctwas_gene_level_direction_and_support.R` | — (terminal) |
| `scTWAS_ADSER_NL_collected_results_scRNAseq_DEsupport.csv` | `multiomics_de_support/01_sctwas_scRNAseq_DEsupport_C_E_L.R` | `13_figure5_legend_and_review.py` |
| `scTWAS_ASC_SSC_linkage_edges.csv` | `multiomics_de_support/04_sctwas_ASC_SSC_linkage_diagram.R` | — (terminal) |
| `scTWAS_DEsupport_summary_UpDown.csv` | `multiomics_de_support/03_sctwas_DEsupport_summary_UpDown.R` | — (terminal) |
| `scTWAS_DEsupport_summary_UpDown_long.csv` | `multiomics_de_support/03_sctwas_DEsupport_summary_UpDown.R` | — (terminal) |
| `scTWAS_bubble_summary.csv` | `multiomics_de_support/12_sctwas_multisource_bubble_heatmap.R` | — (terminal) |
| `scTWAS_bubble_summary_byclass.csv` | `multiomics_de_support/12_sctwas_multisource_bubble_heatmap.R` | — (terminal) |
| `scTWAS_model_gene_overlap_matrix.csv` | `multiomics_de_support/00_sctwas_model_overlap_triangle_13models.R` | — (terminal) |
| `scTWAS_model_gene_overlap_pairs.csv` | `multiomics_de_support/00_sctwas_model_overlap_triangle_13models.R` | — (terminal) |
| `scTWAS_multisource_DEsupport_by_model.csv` | `multiomics_de_support/06_sctwas_multisource_DEsupport_diverging.R` | — (terminal) |
| `spatial_compartment_specificity.csv` | `spatial_geomx_cosmx/20_spatial_de.py` | — (terminal) |
| `spatial_de_summary.json` | `spatial_geomx_cosmx/20_spatial_de.py` | — (terminal) |
| `spatial_direction_concordance.csv` | `spatial_geomx_cosmx/20_spatial_de.py` | `22_spatial_figures.py`, `10_supplementary_tables.py`, `11_report_data.py` |
| `spatial_replication.csv` | `spatial_geomx_cosmx/20_spatial_de.py` | `22_spatial_figures.py`, `10_supplementary_tables.py`, `11_report_data.py` |
| `spatial_risk_gene_DE.csv` | `spatial_geomx_cosmx/20_spatial_de.py` | `22_spatial_figures.py`, `23_validated_set.py`, `31_summary_statistics.py` (+2) |
| `spatial_risk_vs_background.csv` | `spatial_geomx_cosmx/20_spatial_de.py` | `10_supplementary_tables.py`, `11_report_data.py` |
| `validated_by_celltype.csv` | `spatial_geomx_cosmx/23_validated_set.py` | `28_validated_figures.py`, `10_supplementary_tables.py`, `11_report_data.py` |
| `validated_priority.csv` | `spatial_geomx_cosmx/27_validated_priority.py` | `28_validated_figures.py`, `33_functional_stats.py`, `35_crispr_by_cellline.py` (+3) |
| `validated_risk_genes.csv` | `spatial_geomx_cosmx/23_validated_set.py` | `25_cosmx_percell.py`, `26_cosmx_gene_tests.py` |
| `validated_risk_genes_all.csv` | `spatial_geomx_cosmx/23_validated_set.py` | `26_cosmx_gene_tests.py`, `27_validated_priority.py`, `28_validated_figures.py` (+12) |
| `validated_vs_failed.json` | `spatial_geomx_cosmx/27_validated_priority.py` | `28_validated_figures.py`, `11_report_data.py` |

## How this file was built

Automatically generated 2026-10-19 from MANIFEST.csv, the Methods order in `docs/workflow.md`, and the file reads and writes parsed out of the deposited scripts themselves. The extraction reads only string literals and path variables that appear at a read or write call site, so three classes of file are invisible to it: paths assembled from shell variables (most `.sh` and `.sbatch` wrappers), paths built by string formatting at run time, and files handled by external tools that the scripts invoke rather than open themselves (tensorQTL, PLINK, MACS, SMR, PEER). The lists above are therefore a lower bound, not a complete inventory. `MANIFEST.csv` at the repository root is authoritative for the file list itself.

## Added in this release

These files were merged in from the project repository after the first deposit was assembled, so they are not described in the run order above.

The CPTAC and TCGA annotation scripts, the cell-type-matched spatial support check, the general functional-query script, and the 8-panel CosMx highlight figures (the 4-panel versions are also here).

- `FunctionalQuery_updated_v3.R`
- `multiomics_de_support/08_prep_cosmx_highlight_8panel_linux.py`
- `multiomics_de_support/08_run_cosmx_8panel_linux.sh`
- `multiomics_de_support/09_fig_NOR_CORE_highlight_8panel_linux.R`
- `multiomics_de_support/CPTAC_COAD_TWAS_annotation_streamfix_updated.py`
- `multiomics_de_support/CPTAC_env.sh`
- `multiomics_de_support/TCGA_COAD_TWAS_annotation_updated.R`
- `spatial_geomx_cosmx/check_scTWAS_spatial_support_4modules_celltypeMatched.R`
