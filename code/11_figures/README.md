# 11 — Figure assembly

`manhattan/` renders the stacked per-cell-type Manhattan panels for the normal-colon and
polyp model sets, in European-ancestry and multi-ancestry-only variants, with the panel A
source-data extractor. `workflow/` draws the two analytic workflow schematics.
`figure7_depmap_screens/` is the Figure 7 panel layout; note that it is a rendering
function only, and its analysis caller is missing (see `known_gaps.md`).

## Position in the pipeline

```
04_association_testing  02_accessible_chromatin  01_single_cell_and_xqtl
        |
        v
   >>> 11_figures <<<
        |
        v
   (terminal: display items / supplementary tables)
```

The arrows above are the logical order of the published workflow, taken from the Methods (`docs/workflow.md`). The file-level links in *From other sections* and *Consumed by* below are the ones actually observed in the deposited code.

## Scripts, in run order (12 files)

Ordered by sub-directory, then by the numeric prefix the script carries in its original working directory. Descriptions are the script's own header line; blank where the script has no header.

| # | File | Lang | What it does | Declared display item |
|---|---|---|---|---|
| 1 | `figure7_depmap_screens/build_figure7.py` | Python |  |  |
| 2 | `manhattan/Manhattan_ADSER_EUR.R` | R | Manhattan plots for ADSER cell-type-specific TWAS (EUR only) Plotting rules: 1. Plot ALL tested EUR TWAS genes. |  |
| 3 | `manhattan/Manhattan_ADSER_Multi_EURnotSig_MultiSig_only.R` | R | Manhattan plots for ADSER cell-type-specific TWAS (Multi ancestry) Plotting rules: 1. Plot ALL tested Multi-ancestry TWAS genes. |  |
| 4 | `manhattan/Manhattan_NL_EUR.R` | R | Manhattan plots for HTAN_NL + WellcomeSanger_NL Cauchy meta-TWAS (EUR) Input files: HTAN_NL_WellcomeSanger_NL_{ct}_peakSNP.EUR.CauchyMeta.csv Cell typ. |  |
| 5 | `manhattan/Manhattan_NL_Multi_EURnotSig_MultiSig_only.R` | R | Manhattan plots for HTAN_NL + WellcomeSanger_NL Cauchy meta-TWAS (Multi ancestry) Input files: HTAN_NL_WellcomeSanger_NL_{ct}_peakSNP.Multi.CauchyMeta. |  |
| 6 | `manhattan/extract_panelA_sourcedata.R` | R | WHAT IT DOES Reads a Seurat .rds and writes the small cell-level table that Final_Viz_4panel.R needs in order to DRAW panel A itself: cell_id, UMAPHAR. |  |
| 7 | `manhattan/extract_panelA_sourcedata.sh` | Shell |  |  |
| 8 | `manhattan/run_umap_WS_allcells.sh` | Shell |  |  |
| 9 | `model_performance/Fig_peaks_R2_violin.R` | R | PURPOSE TWO figures, one per model group, each one row x two columns: Figure 1 Normal colon a) scATAC-seq peaks per cell type (bar) b) global R^2 of t. | Figure 1; Figure 2; Supplementary Table 5 |
| 10 | `umap_annotation/VUMC_NL/05_annotation_UMAP.R` | R | 05_annotation.R Annotate clusters using the dominant ORIGINAL cell type in each cluster. |  |
| 11 | `umap_annotation/WellcomeSanger_NL/04_umap_epithelial_LauLabel_04c_filtered.R` | R | Wellcome Sanger healthy epithelium: BPCells -> Seurat -> three labelled UMAPs. |  |
| 12 | `workflow/Fig_scTWAS_workflow_NL_ADSER.py` | Python |  | Supplementary Figure 2 |

## Inputs

**Not produced anywhere in this deposit** (2) — cohort data, public resources, or output of a step that is not deposited:

- `05_Epithelial_cell_metadata_LauLabel.csv.gz`
- `SupplementaryTables.xlsx`

## Outputs

| File | Written by | Consumed by |
|---|---|---|
| `04_umap_R_sessionInfo.txt` | `ell_and_xqtl/normal_colon_WellcomeSanger/04_umap_epithelial_LauLabel_revised.R` | — (terminal) |
| `04c_fine_annotation_filter_counts.csv` | `umap_annotation/WellcomeSanger_NL/04_umap_epithelial_LauLabel_04c_filtered.R` | — (terminal) |
| `ADSER_EUR_TWAS_R2gt0.01_Manhattan_Bonferroni_summary.csv` | `manhattan/Manhattan_ADSER_EUR.R` | — (terminal) |
| `ADSER_EUR_TWAS_R2gt0.01_Manhattan_plotted_genes.csv` | `manhattan/Manhattan_ADSER_EUR.R` | — (terminal) |
| `ADSER_Multi_TWAS_R2gt0.01_Manhattan_highlight_summary.csv` | `manhattan/Manhattan_ADSER_Multi_EURnotSig_MultiSig_only.R` | — (terminal) |
| `ADSER_Multi_TWAS_R2gt0.01_Manhattan_plotted_genes.csv` | `manhattan/Manhattan_ADSER_Multi_EURnotSig_MultiSig_only.R` | — (terminal) |
| `Fig_peaks_R2_violin_sessionInfo.txt` | `model_performance/Fig_peaks_R2_violin.R` | — (terminal) |
| `HTAN_NL_WellcomeSanger_NL_EUR_CauchyMeta_7CT_Manhattan_Bonferroni_summary.csv` | `manhattan/Manhattan_NL_EUR.R` | — (terminal) |
| `HTAN_NL_WellcomeSanger_NL_EUR_CauchyMeta_7CT_Manhattan_plotted_genes.csv` | `manhattan/Manhattan_NL_EUR.R` | — (terminal) |
| `HTAN_NL_WellcomeSanger_NL_Multi_CauchyMeta_7CT_Manhattan_highlight_summary.csv` | `manhattan/Manhattan_NL_Multi_EURnotSig_MultiSig_only.R` | — (terminal) |
| `HTAN_NL_WellcomeSanger_NL_Multi_CauchyMeta_7CT_Manhattan_plotted_genes.csv` | `manhattan/Manhattan_NL_Multi_EURnotSig_MultiSig_only.R` | — (terminal) |
| `_R2_sourcedata.csv` | `model_performance/Fig_peaks_R2_violin.R` | — (terminal) |
| `_peaks_sourcedata.csv` | `model_performance/Fig_peaks_R2_violin.R` | — (terminal) |
| `_percelltype.csv` | `model_performance/Fig_peaks_R2_violin.R` | — (terminal) |
| `_sampling.csv` | `umap_annotation/VUMC_NL/05_annotation_UMAP.R` | — (terminal) |
| `cell_type_colors.csv` | `ell_and_xqtl/normal_colon_WellcomeSanger/04_umap_epithelial_LauLabel_revised.R` | — (terminal) |
| `cluster_LauLabel_composition.csv` | `ell_and_xqtl/normal_colon_WellcomeSanger/04_umap_epithelial_LauLabel_revised.R` | — (terminal) |
| `cluster_by_original_CellType_cb_annot_pct_res0.8.rm23_sub.csv` | `ell_and_xqtl/normal_colon_GSE266546/05_annotation_UMAP.R` | — (terminal) |
| `cluster_by_original_CellType_cb_res0.8.rm23_sub.csv` | `ell_and_xqtl/normal_colon_GSE266546/05_annotation_UMAP.R` | — (terminal) |
| `cluster_celltype_composition_res0.8.rm23_sub.csv` | `ell_and_xqtl/normal_colon_GSE266546/05_annotation_UMAP.R` | — (terminal) |
| `cluster_dominant_celltype_res0.8.rm23_sub.csv` | `ell_and_xqtl/normal_colon_GSE266546/05_annotation_UMAP.R` | — (terminal) |
| `cluster_majority_LauLabel.csv` | `ell_and_xqtl/normal_colon_WellcomeSanger/04_umap_epithelial_LauLabel_revised.R` | — (terminal) |
| `cluster_reviewed_annotation_res0.8.rm23_sub.csv` | `ell_and_xqtl/normal_colon_GSE266546/05_annotation_UMAP.R` | — (terminal) |
| `fine_annotation_mapping.csv` | `ell_and_xqtl/normal_colon_WellcomeSanger/04_umap_epithelial_LauLabel_revised.R` | — (terminal) |
| `integrated_seurat_final2.res0.8.rm23_sub.rds` | `ell_and_xqtl/normal_colon_GSE266546/05_annotation_UMAP.R` | — (terminal) |
| `integrated_seurat_final2_reviewed.res0.8.rm23_sub.rds` | `ell_and_xqtl/normal_colon_GSE266546/05_annotation_UMAP.R` | `explore_specimen_GCA.R` |
| `integrated_seurat_final2_reviewed.res0.8.rm23_sub_metadata.csv` | `ell_and_xqtl/normal_colon_GSE266546/05_annotation_UMAP.R` | `explore_specimen_GCA.R` |
| `umap_WellcomeSanger_epithelial_clusters_LauLabel_sourcedata.csv` | `ell_and_xqtl/normal_colon_WellcomeSanger/04_umap_epithelial_LauLabel_revised.R` | — (terminal) |

## How this file was built

Automatically generated 2026-10-19 from MANIFEST.csv, the Methods order in `docs/workflow.md`, and the file reads and writes parsed out of the deposited scripts themselves. The extraction reads only string literals and path variables that appear at a read or write call site, so three classes of file are invisible to it: paths assembled from shell variables (most `.sh` and `.sbatch` wrappers), paths built by string formatting at run time, and files handled by external tools that the scripts invoke rather than open themselves (tensorQTL, PLINK, MACS, SMR, PEER). The lists above are therefore a lower bound, not a complete inventory. `MANIFEST.csv` at the repository root is authoritative for the file list itself.

## Added in this release

These files were merged in from the project repository after the first deposit was assembled, so they are not described in the run order above.

`Final_Viz_4panel.R` and the model R-squared violin with peak barplots.

- `Final_Viz_4panel.R`
- `model_performance/Violin_ADSER_NL_EUR_R2gt0.01_with_peak_barplots.R`
