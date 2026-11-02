# 12 — Supplementary table assembly

Assembly of the supplementary workbooks and the consolidated results record spanning the
orthogonal, druggability and screening layers.

## Position in the pipeline

```
07_orthogonal_support  08_druggability_and_connectivity  09_functional_screens  10_pathway_enrichment
        |
        v
   >>> 12_supplementary_tables <<<
        |
        v
   (terminal: display items / supplementary tables)
```

The arrows above are the logical order of the published workflow, taken from the Methods (`docs/workflow.md`). The file-level links in *From other sections* and *Consumed by* below are the ones actually observed in the deposited code.

## Scripts, in run order (3 files)

Ordered by sub-directory, then by the numeric prefix the script carries in its original working directory. Descriptions are the script's own header line; blank where the script has no header.

| # | File | Lang | What it does | Declared display item |
|---|---|---|---|---|
| 1 | `10_supplementary_tables.py` | Python | 10 - Assemble the supplementary workbook. |  |
| 2 | `11_report_data.py` | Python | 11 - Emit the numbers the report needs, so the document never restates a value by hand. |  |
| 3 | `_patch_validated.py` | Python | Splice the validated-gene section (2.14) and its methods into the report builder. |  |

## Inputs

**From other sections of this deposit**

| File | Produced by |
|---|---|
| `cosmx_celltype_trend.csv` | `07_orthogonal_support/spatial_geomx_cosmx/21_cosmx_composition.py` |
| `cosmx_crc_vs_normal_epi.csv` | `07_orthogonal_support/spatial_geomx_cosmx/26_cosmx_gene_tests.py` |
| `cosmx_gene_stage_trend.csv` | `07_orthogonal_support/spatial_geomx_cosmx/26_cosmx_gene_tests.py` |
| `cosmx_geomx_replication.csv` | `07_orthogonal_support/spatial_geomx_cosmx/26_cosmx_gene_tests.py` |
| `cosmx_niche_trend.csv` | `07_orthogonal_support/spatial_geomx_cosmx/21_cosmx_composition.py` |
| `cosmx_requested_genes.csv` | `07_orthogonal_support/spatial_geomx_cosmx/47_cosmx_requested_genes.py` |
| `risk_gene_all_comparisons.csv` | `07_orthogonal_support/spatial_geomx_cosmx/31_summary_statistics.py` |
| `risk_gene_summary_statistics.csv` | `07_orthogonal_support/spatial_geomx_cosmx/31_summary_statistics.py` |
| `spatial_direction_concordance.csv` | `07_orthogonal_support/spatial_geomx_cosmx/20_spatial_de.py` |
| `spatial_replication.csv` | `07_orthogonal_support/spatial_geomx_cosmx/20_spatial_de.py` |
| `spatial_risk_gene_DE.csv` | `07_orthogonal_support/spatial_geomx_cosmx/20_spatial_de.py` |
| `spatial_risk_vs_background.csv` | `07_orthogonal_support/spatial_geomx_cosmx/20_spatial_de.py` |
| `validated_by_celltype.csv` | `07_orthogonal_support/spatial_geomx_cosmx/23_validated_set.py` |
| `validated_priority.csv` | `07_orthogonal_support/spatial_geomx_cosmx/27_validated_priority.py` |
| `validated_risk_genes_all.csv` | `07_orthogonal_support/spatial_geomx_cosmx/23_validated_set.py` |
| `validated_vs_failed.json` | `07_orthogonal_support/spatial_geomx_cosmx/27_validated_priority.py` |
| `risk_genes_harmonized.csv` | `08_druggability_and_connectivity/gene_sets_and_resources/01_parse_risk_genes.py` |
| `cosmx_selected_genes.csv` | `09_functional_screens/48_screened_and_selected.py` |
| `crispr_by_cellline.csv` | `09_functional_screens/35_crispr_by_cellline.py` |
| `crispr_concordance_all_plates.json` | `09_functional_screens/42_crispr_replicate_qc.py` |
| `crispr_plate_correlations.csv` | `09_functional_screens/43_crispr_concordance_heatmap.py` |
| `crispr_replicate_qc.csv` | `09_functional_screens/42_crispr_replicate_qc.py` |
| `crispr_risk36_by_cellline.csv` | `09_functional_screens/35_crispr_by_cellline.py` |
| `func_crispr_gene_stats.csv` | `09_functional_screens/33_functional_stats.py` |
| `func_qc.csv` | `09_functional_screens/32_functional_parse_qc.py` |
| `func_risk_gene_results.csv` | `09_functional_screens/33_functional_stats.py` |
| `func_sirna_gene_stats.csv` | `09_functional_screens/33_functional_stats.py` |
| `func_sirna_vs_crispr.csv` | `09_functional_screens/33_functional_stats.py` |
| `func_summary.json` | `09_functional_screens/33_functional_stats.py` |
| `highconf_replication.csv` | `09_functional_screens/48_screened_and_selected.py` |
| `highconf_replication.json` | `09_functional_screens/48_screened_and_selected.py` |
| `nontwas_by_cellline.csv` | `09_functional_screens/39_nontwas_figure.py` |
| `nontwas_comparison.json` | `09_functional_screens/39_nontwas_figure.py` |
| `risk_increasing_both_assays.csv` | `09_functional_screens/45_risk_increasing_figure.py` |
| `risk_increasing_crispr.csv` | `09_functional_screens/46_risk_increasing_crispr.py` |
| `risk_increasing_crispr.json` | `09_functional_screens/46_risk_increasing_crispr.py` |
| `risk_increasing_summary.json` | `09_functional_screens/45_risk_increasing_figure.py` |
| `shared_genes_both_assays.csv` | `09_functional_screens/44_shared_genes_figure.py` |
| `sirna_by_cellline.csv` | `09_functional_screens/37_sirna_by_cellline.py` |
| `sirna_risk_by_cellline.csv` | `09_functional_screens/37_sirna_by_cellline.py` |
| `all347_by_direction.csv` | `10_pathway_enrichment/53_all_genes_etiology.py` |
| `all347_enrichment.csv` | `10_pathway_enrichment/53_all_genes_etiology.py` |
| `all347_etiology.json` | `10_pathway_enrichment/53_all_genes_etiology.py` |
| `canonical_membership.csv` | `10_pathway_enrichment/53_all_genes_etiology.py` |
| `canonical_membership_summary.csv` | `10_pathway_enrichment/53_all_genes_etiology.py` |
| `etiology_gene_assignment.csv` | `10_pathway_enrichment/53_all_genes_etiology.py` |
| `etiology_themes.csv` | `10_pathway_enrichment/53_all_genes_etiology.py` |
| `gene_pathway_proximity.csv` | `10_pathway_enrichment/53_all_genes_etiology.py` |
| `pathway_celltype.json` | `10_pathway_enrichment/49_pathway_by_celltype.py` |
| `pathway_celltype_enrichment.csv` | `10_pathway_enrichment/49_pathway_by_celltype.py` |
| `pathway_celltype_null.csv` | `10_pathway_enrichment/49_pathway_by_celltype.py` |
| `pathway_celltype_significant.csv` | `10_pathway_enrichment/49_pathway_by_celltype.py` |
| `pathway_leading_genes.csv` | `10_pathway_enrichment/49_pathway_by_celltype.py` |
| `pathway_merged.json` | `10_pathway_enrichment/51_pathway_merged.py` |
| `pathway_merged_enrichment.csv` | `10_pathway_enrichment/51_pathway_merged.py` |
| `pathway_merged_null.csv` | `10_pathway_enrichment/51_pathway_merged.py` |
| `pathway_merged_significant.csv` | `10_pathway_enrichment/51_pathway_merged.py` |
| `pathway_merged_vs_split.csv` | `10_pathway_enrichment/51_pathway_merged.py` |
| `pathway_novel_recurrent.csv` | `10_pathway_enrichment/49_pathway_by_celltype.py` |
| `pathway_pooled_enrichment.csv` | `10_pathway_enrichment/49_pathway_by_celltype.py` |

**Not produced anywhere in this deposit** (49) — cohort data, public resources, or output of a step that is not deposited:

- `cosmx_percell_summary.json`
- `cosmx_selected_genes.json`
- `cosmx_summary.json`
- `crispr_by_cellline.json`
- `crispr_replicate_qc.json`
- `driver_overlap.csv`
- `driver_overlap.json`
- `driver_perturbation_by_ko.csv`
- `dti_predictions_validated_regime.csv`
- `epistasis_control.json`
- `epistasis_controlled.csv`
- `func_qc.json`
- `gdsc_proximity_bydrug.csv`
- `gene_drug_candidates.csv`
- `gene_druggability.csv`
- `gene_ko_impact.csv`
- `gene_tractability_score.csv`
- `integrated_gene_priority.csv`
- `ko_convergence.json`
- `ko_direction_vs_crc_state.csv`
- `ko_downstream_control.json`
- `ko_downstream_pathways.csv`
- `ko_downstream_pathways_recurrent.csv`
- `ko_downstream_vs_background.csv`
- `ko_driver_analysis.json`
- `ko_driver_spectrum.json`
- `ko_impact_spectrum_full.csv`
- `ko_module_convergence.csv`
- `lineage_network.json`
- `lineage_network_proximity.csv`
- `lineage_own_vs_other.csv`
- `module_knockout_impact.csv`
- `module_pathway_proximity.csv`
- `network_proximity.csv.gz`
- `network_proximity_significant.csv`
- `pathway_enrichment.csv`
- `reversal_positive_control.csv`
- `reversal_positive_control_knowndrugs.csv`
- `risk_gene_driver_proximity.csv`
- `signature_reversal_all.csv.gz`
- `signature_reversal_recurrent.csv`
- `signature_reversal_summary.csv`
- `signature_reversal_top20.csv`
- `sirna_by_cellline.json`
- `sirna_plate_correlations.csv`
- `stage_contrast.csv`
- `stage_epistasis.json`
- `summary_statistics_overview.json`
- `validated_summary.json`

## Outputs

| File | Written by | Consumed by |
|---|---|---|
| `report_data.json` | `11_report_data.py` | — (terminal) |

## How this file was built

Automatically generated 2026-10-19 from MANIFEST.csv, the Methods order in `docs/workflow.md`, and the file reads and writes parsed out of the deposited scripts themselves. The extraction reads only string literals and path variables that appear at a read or write call site, so three classes of file are invisible to it: paths assembled from shell variables (most `.sh` and `.sbatch` wrappers), paths built by string formatting at run time, and files handled by external tools that the scripts invoke rather than open themselves (tensorQTL, PLINK, MACS, SMR, PEER). The lists above are therefore a lower bound, not a complete inventory. `MANIFEST.csv` at the repository root is authoritative for the file list itself.
