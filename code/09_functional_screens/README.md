# 09 — Arrayed functional screens

Parsing and quality control of the experimenters' screen workbook, per-cell-line siRNA and
CRISPR-Cas9 statistics with paired t-tests across replicate gene-control pairs, replicate
quality control, concordance, risk-direction and worked-example figures, and the
screened-versus-selected accounting.

## Position in the pipeline

```
07_orthogonal_support  08_druggability_and_connectivity
        |
        v
   >>> 09_functional_screens <<<
        |
        v
10_pathway_enrichment  12_supplementary_tables
```

The arrows above are the logical order of the published workflow, taken from the Methods (`docs/workflow.md`). The file-level links in *From other sections* and *Consumed by* below are the ones actually observed in the deposited code.

## Scripts, in run order (18 files)

Ordered by sub-directory, then by the numeric prefix the script carries in its original working directory. Descriptions are the script's own header line; blank where the script has no header.

| # | File | Lang | What it does | Declared display item |
|---|---|---|---|---|
| 1 | `32_functional_parse_qc.py` | Python |  |  |
| 2 | `33_functional_stats.py` | Python |  |  |
| 3 | `34_functional_figures.py` | Python | 34 - Figure 19: the siRNA and CRISPR viability screens. | [internal report numbering: Figure 19] |
| 4 | `35_crispr_by_cellline.py` | Python |  | [internal report numbering: Figure 23] |
| 5 | `36_crispr_cellline_figure.py` | Python |  | [internal report numbering: Figure 20] |
| 6 | `37_sirna_by_cellline.py` | Python |  |  |
| 7 | `38_sirna_cellline_figure.py` | Python |  | [internal report numbering: Figure 21; Figure 20] |
| 8 | `39_nontwas_figure.py` | Python |  | [internal report numbering: Figure 22; Figure 20; Figure 20 - a] |
| 9 | `40_ttest_diagnostics.py` | Python |  | [internal report numbering: Figure 20] |
| 10 | `41_rps5_worked_example.py` | Python |  |  |
| 11 | `42_crispr_replicate_qc.py` | Python |  |  |
| 12 | `43_concordance_heatmap.py` | Python |  | [internal report numbering: Figure 23; Figure 24] |
| 13 | `43_crispr_concordance_heatmap.py` | Python |  | [internal report numbering: Figure 23] |
| 14 | `44_shared_genes_figure.py` | Python |  | [internal report numbering: Figure 25; Figure 20a; Figure 20] |
| 15 | `45_risk_increasing_figure.py` | Python |  | [internal report numbering: Figure 26; Figure 20] |
| 16 | `46_risk_increasing_crispr.py` | Python |  | [internal report numbering: Figure 27; Figure 26; Figure 20] |
| 17 | `48_screened_and_selected.py` | Python |  | [internal report numbering: Figure 29; Figure 16f; Figure 30; Figure 28] |
| 18 | `_ttest_support.py` | Python |  | [internal report numbering: Figure 20] |

## Inputs

**From other sections of this deposit**

| File | Produced by |
|---|---|
| `cosmx_crc_vs_normal_epi.csv` | `07_orthogonal_support/spatial_geomx_cosmx/26_cosmx_gene_tests.py` |
| `cosmx_geomx_replication.csv` | `07_orthogonal_support/spatial_geomx_cosmx/26_cosmx_gene_tests.py` |
| `cx_pseudobulk.csv.gz` | `07_orthogonal_support/spatial_geomx_cosmx/25_cosmx_percell.py` |
| `validated_priority.csv` | `07_orthogonal_support/spatial_geomx_cosmx/27_validated_priority.py` |
| `validated_risk_genes_all.csv` | `07_orthogonal_support/spatial_geomx_cosmx/23_validated_set.py` |
| `risk_genes_harmonized.csv` | `08_druggability_and_connectivity/gene_sets_and_resources/01_parse_risk_genes.py` |

**Not produced anywhere in this deposit** (9) — cohort data, public resources, or output of a step that is not deposited:

- `crispr_by_cellline.json`
- `crispr_replicate_qc.json`
- `func_qc.json`
- `gene_druggability.csv`
- `high-confidence-genes.xlsx`
- `integrated_gene_priority.csv`
- `risk_gene_driver_proximity.csv`
- `siRNA_CRISPR-116.xlsx`
- `sirna_by_cellline.json`

## Outputs

| File | Written by | Consumed by |
|---|---|---|
| `cosmx_selected_genes.csv` | `48_screened_and_selected.py` | `14_build_supptable_S20_cosmx.py`, `10_supplementary_tables.py` |
| `cosmx_selected_genes_stats.csv` | `48_screened_and_selected.py` | `14_build_supptable_S20_cosmx.py` |
| `crispr_by_cellline.csv` | `35_crispr_by_cellline.py` | `36_crispr_cellline_figure.py`, `39_nontwas_figure.py`, `45_risk_increasing_figure.py` (+5) |
| `crispr_concordance_all_plates.json` | `42_crispr_replicate_qc.py` | `11_report_data.py` |
| `crispr_plate_correlations.csv` | `43_crispr_concordance_heatmap.py` | `10_supplementary_tables.py` |
| `crispr_replicate_qc.csv` | `42_crispr_replicate_qc.py` | `10_supplementary_tables.py`, `11_report_data.py` |
| `crispr_risk36_by_cellline.csv` | `35_crispr_by_cellline.py` | `36_crispr_cellline_figure.py`, `37_sirna_by_cellline.py`, `38_sirna_cellline_figure.py` (+5) |
| `func_crispr_gene_stats.csv` | `33_functional_stats.py` | `34_functional_figures.py`, `10_supplementary_tables.py`, `11_report_data.py` |
| `func_crispr_long.csv` | `32_functional_parse_qc.py` | `33_functional_stats.py`, `35_crispr_by_cellline.py`, `40_ttest_diagnostics.py` (+3) |
| `func_qc.csv` | `32_functional_parse_qc.py` | `34_functional_figures.py`, `37_sirna_by_cellline.py`, `42_crispr_replicate_qc.py` (+1) |
| `func_qc_crossline.csv` | `32_functional_parse_qc.py` | — (terminal) |
| `func_qc_replicates.csv` | `32_functional_parse_qc.py` | — (terminal) |
| `func_risk_gene_results.csv` | `33_functional_stats.py` | `34_functional_figures.py`, `10_supplementary_tables.py`, `11_report_data.py` |
| `func_sirna_gene_stats.csv` | `33_functional_stats.py` | `34_functional_figures.py`, `10_supplementary_tables.py`, `11_report_data.py` |
| `func_sirna_long.csv` | `32_functional_parse_qc.py` | `33_functional_stats.py`, `37_sirna_by_cellline.py` |
| `func_sirna_vs_crispr.csv` | `33_functional_stats.py` | `34_functional_figures.py`, `10_supplementary_tables.py`, `11_report_data.py` |
| `func_summary.json` | `33_functional_stats.py` | `34_functional_figures.py`, `11_report_data.py` |
| `highconf_replication.csv` | `48_screened_and_selected.py` | `10_supplementary_tables.py`, `11_report_data.py` |
| `highconf_replication.json` | `48_screened_and_selected.py` | `11_report_data.py` |
| `nontwas_by_cellline.csv` | `39_nontwas_figure.py` | `10_supplementary_tables.py`, `11_report_data.py` |
| `nontwas_comparison.json` | `39_nontwas_figure.py` | `11_report_data.py` |
| `risk_increasing_both_assays.csv` | `45_risk_increasing_figure.py` | `10_supplementary_tables.py` |
| `risk_increasing_crispr.csv` | `46_risk_increasing_crispr.py` | `10_supplementary_tables.py` |
| `risk_increasing_crispr.json` | `46_risk_increasing_crispr.py` | `11_report_data.py` |
| `risk_increasing_summary.json` | `45_risk_increasing_figure.py` | `11_report_data.py` |
| `rps5_input_data.csv` | `41_rps5_worked_example.py` | — (terminal) |
| `shared_genes_both_assays.csv` | `44_shared_genes_figure.py` | `10_supplementary_tables.py`, `11_report_data.py` |
| `sirna_by_cellline.csv` | `37_sirna_by_cellline.py` | `48_screened_and_selected.py`, `49_pathway_by_celltype.py`, `53_all_genes_etiology.py` (+1) |
| `sirna_risk_by_cellline.csv` | `37_sirna_by_cellline.py` | `38_sirna_cellline_figure.py`, `44_shared_genes_figure.py`, `45_risk_increasing_figure.py` (+2) |
| `ttest_diagnostics.csv` | `40_ttest_diagnostics.py` | — (terminal) |

## How this file was built

Automatically generated 2026-10-19 from MANIFEST.csv, the Methods order in `docs/workflow.md`, and the file reads and writes parsed out of the deposited scripts themselves. The extraction reads only string literals and path variables that appear at a read or write call site, so three classes of file are invisible to it: paths assembled from shell variables (most `.sh` and `.sbatch` wrappers), paths built by string formatting at run time, and files handled by external tools that the scripts invoke rather than open themselves (tensorQTL, PLINK, MACS, SMR, PEER). The lists above are therefore a lower bound, not a complete inventory. `MANIFEST.csv` at the repository root is authoritative for the file list itself.

## Added in this release

These files were merged in from the project repository after the first deposit was assembled, so they are not described in the run order above.

`experimental_validation/` -- the DESeq2 differential-expression analysis of the knockout validation RNA-seq, with its EnrichR and CRC-pathway mapping steps.

- `experimental_validation/15174-XG.l2_bioinfo1.sbatch.sh`
- `experimental_validation/1_DESeq2_DE_Genes_corrected_final.R`
- `experimental_validation/2_DESeq2_EnrichR_corrected_final.R`
- `experimental_validation/3a_DESeq2_genes_EnrichR_mapto_CRC_relevance.R`
- `experimental_validation/3b_DESeq2_genes_map_to_CRC_KEGG_pathways.R`
