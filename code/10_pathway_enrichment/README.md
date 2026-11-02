# 10 — Pathway and gene-set enrichment

Enrichment by cell type and by risk direction, the merged analyses and their figures.

## Position in the pipeline

```
08_druggability_and_connectivity  09_functional_screens
        |
        v
   >>> 10_pathway_enrichment <<<
        |
        v
12_supplementary_tables
```

The arrows above are the logical order of the published workflow, taken from the Methods (`docs/workflow.md`). The file-level links in *From other sections* and *Consumed by* below are the ones actually observed in the deposited code.

## Scripts, in run order (6 files)

Ordered by sub-directory, then by the numeric prefix the script carries in its original working directory. Descriptions are the script's own header line; blank where the script has no header.

| # | File | Lang | What it does | Declared display item |
|---|---|---|---|---|
| 1 | `49_pathway_by_celltype.py` | Python |  |  |
| 2 | `50_pathway_figure.py` | Python |  | [internal report numbering: Figure 31] |
| 3 | `51_pathway_merged.py` | Python |  |  |
| 4 | `52_pathway_merged_figures.py` | Python |  | [internal report numbering: Figure 33; Figure 31; Figure 34; Figure 33 u] |
| 5 | `53_all_genes_etiology.py` | Python |  |  |
| 6 | `54_etiology_figures.py` | Python |  | [internal report numbering: Figure 35; Figure 36] |

## Inputs

**From other sections of this deposit**

| File | Produced by |
|---|---|
| `validated_risk_genes_all.csv` | `07_orthogonal_support/spatial_geomx_cosmx/23_validated_set.py` |
| `gene_id_map.json` | `08_druggability_and_connectivity/gene_sets_and_resources/02_build_resources.py` |
| `risk_genes_harmonized.csv` | `08_druggability_and_connectivity/gene_sets_and_resources/01_parse_risk_genes.py` |
| `crispr_by_cellline.csv` | `09_functional_screens/35_crispr_by_cellline.py` |
| `sirna_by_cellline.csv` | `09_functional_screens/37_sirna_by_cellline.py` |

**Not produced anywhere in this deposit** (1) — cohort data, public resources, or output of a step that is not deposited:

- `gene_druggability.csv`

## Outputs

| File | Written by | Consumed by |
|---|---|---|
| `all347_by_direction.csv` | `53_all_genes_etiology.py` | `10_supplementary_tables.py` |
| `all347_enrichment.csv` | `53_all_genes_etiology.py` | `54_etiology_figures.py`, `10_supplementary_tables.py` |
| `all347_etiology.json` | `53_all_genes_etiology.py` | `54_etiology_figures.py`, `11_report_data.py` |
| `all347_nulls.json` | `53_all_genes_etiology.py` | `54_etiology_figures.py` |
| `canonical_membership.csv` | `53_all_genes_etiology.py` | `10_supplementary_tables.py`, `11_report_data.py` |
| `canonical_membership_summary.csv` | `53_all_genes_etiology.py` | `54_etiology_figures.py`, `11_report_data.py` |
| `etiology_gene_assignment.csv` | `53_all_genes_etiology.py` | `54_etiology_figures.py`, `10_supplementary_tables.py`, `11_report_data.py` |
| `etiology_themes.csv` | `53_all_genes_etiology.py` | `54_etiology_figures.py`, `10_supplementary_tables.py`, `11_report_data.py` |
| `gene_pathway_proximity.csv` | `53_all_genes_etiology.py` | `10_supplementary_tables.py` |
| `pathway_celltype.json` | `49_pathway_by_celltype.py` | `50_pathway_figure.py`, `11_report_data.py` |
| `pathway_celltype_enrichment.csv` | `49_pathway_by_celltype.py` | `10_supplementary_tables.py` |
| `pathway_celltype_null.csv` | `49_pathway_by_celltype.py` | `50_pathway_figure.py`, `51_pathway_merged.py`, `10_supplementary_tables.py` (+1) |
| `pathway_celltype_significant.csv` | `49_pathway_by_celltype.py` | `50_pathway_figure.py`, `51_pathway_merged.py`, `52_pathway_merged_figures.py` (+1) |
| `pathway_families.csv` | `49_pathway_by_celltype.py` | — (terminal) |
| `pathway_gene_core.csv` | `49_pathway_by_celltype.py` | — (terminal) |
| `pathway_leading_genes.csv` | `49_pathway_by_celltype.py` | `50_pathway_figure.py`, `10_supplementary_tables.py`, `11_report_data.py` |
| `pathway_merged.json` | `51_pathway_merged.py` | `52_pathway_merged_figures.py`, `11_report_data.py` |
| `pathway_merged_enrichment.csv` | `51_pathway_merged.py` | `10_supplementary_tables.py` |
| `pathway_merged_families.csv` | `51_pathway_merged.py` | `52_pathway_merged_figures.py` |
| `pathway_merged_null.csv` | `51_pathway_merged.py` | `52_pathway_merged_figures.py`, `10_supplementary_tables.py`, `11_report_data.py` |
| `pathway_merged_significant.csv` | `51_pathway_merged.py` | `52_pathway_merged_figures.py`, `11_report_data.py` |
| `pathway_merged_vs_split.csv` | `51_pathway_merged.py` | `52_pathway_merged_figures.py`, `10_supplementary_tables.py`, `11_report_data.py` |
| `pathway_novel_recurrent.csv` | `49_pathway_by_celltype.py` | `10_supplementary_tables.py` |
| `pathway_pooled_enrichment.csv` | `49_pathway_by_celltype.py` | `10_supplementary_tables.py`, `11_report_data.py` |

## How this file was built

Automatically generated 2026-10-19 from MANIFEST.csv, the Methods order in `docs/workflow.md`, and the file reads and writes parsed out of the deposited scripts themselves. The extraction reads only string literals and path variables that appear at a read or write call site, so three classes of file are invisible to it: paths assembled from shell variables (most `.sh` and `.sbatch` wrappers), paths built by string formatting at run time, and files handled by external tools that the scripts invoke rather than open themselves (tensorQTL, PLINK, MACS, SMR, PEER). The lists above are therefore a lower bound, not a complete inventory. `MANIFEST.csv` at the repository root is authoritative for the file list itself.

## Added in this release

These files were merged in from the project repository after the first deposit was assembled, so they are not described in the run order above.

`EnrichR.R`, the general EnrichR driver.

- `EnrichR.R`
