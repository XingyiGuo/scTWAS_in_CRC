# Data flow

Two views of the same pipeline. The **logical flow** is the order of the published workflow, taken from the Methods. The **observed flow** is what the deposited code actually reads and writes, parsed from the scripts; it is a lower bound, because paths built from shell variables or string formatting at run time, and files handled by external tools (tensorQTL, PLINK, MACS, SMR, PEER), are invisible to the parse.

## Logical flow

```
  cohort scRNA-seq + genotypes            scATAC-seq (VUMC, GSE201336)
            |                                        |
            v                                        v
  01_single_cell_and_xqtl                 02_accessible_chromatin
  per-donor pseudobulk, PEER,             ArchR / MACS peak calling,
  cis-eQTL mapping (tensorQTL)            cross-cohort summit merge
            |                                        |
            +--------------------+-------------------+
                                 v
                   03_grex_model_training   [NO CODE]
                   elastic net, alpha = 0.5, peak-SNP predictors
                   -> 23 *_peakSNP.db model databases
                                 |
                                 v
                   04_association_testing   [NO CODE]
                   S-PrediXcan-style Z, ACAT Cauchy meta
                   -> significant gene x cell-type tables
                                 |
     +----------------+----------+-----------+------------------+
     v                v                      v                  v
 05_coloc_        06_tf_network      08_druggability_and_    11_figures
 smr_heidi        TF-target          connectivity            Manhattan,
 coloc, SMR,      correlation,       01_parse_risk_genes.py  peaks/R2,
 HEIDI            lead-variant       -> risk_genes_          UMAPs,
                  occupancy,         harmonized.csv          workflow
                  networks                   |
                                             v
                                   07_orthogonal_support
                                   scRNA/pseudobulk DE, CPTAC,
                                   GeoMx, CosMx
                                             |
                                             v
                                   09_functional_screens
                                   siRNA + CRISPR-Cas9 viability
                                             |
                                             v
                                   10_pathway_enrichment
                                             |
                                             v
                                   12_supplementary_tables
```

The break at steps 3 and 4 is real: those two directories hold no code (`known_gaps.md`). In the observed flow below, sections 01-02 and 05-12 are therefore two disconnected components — nothing in the deposit carries data across the gap.

## Observed flow between sections

75 file-level links across 12 section pairs.

| Producer | Consumer | Files |
|---|---|---|
| `07_orthogonal_support` | `08_druggability_and_connectivity` | 1 |
| `07_orthogonal_support` | `09_functional_screens` | 5 |
| `07_orthogonal_support` | `10_pathway_enrichment` | 1 |
| `07_orthogonal_support` | `12_supplementary_tables` | 16 |
| `08_druggability_and_connectivity` | `07_orthogonal_support` | 1 |
| `08_druggability_and_connectivity` | `09_functional_screens` | 1 |
| `08_druggability_and_connectivity` | `10_pathway_enrichment` | 2 |
| `08_druggability_and_connectivity` | `12_supplementary_tables` | 1 |
| `09_functional_screens` | `07_orthogonal_support` | 2 |
| `09_functional_screens` | `10_pathway_enrichment` | 2 |
| `09_functional_screens` | `12_supplementary_tables` | 23 |
| `10_pathway_enrichment` | `12_supplementary_tables` | 20 |

Sections `01`, `02`, `05`, `06`, `11` and `manuscript_text` have **no observed link to any other section**. For `01`, `02` and `06` that is the 03/04 gap. For `05` it is because the colocalization and SMR drivers build every path from shell and config variables, so not one resolves to a filename. For `11` and `manuscript_text` it is because they read exported result tables and supplementary workbooks that live outside the deposit.

## A note on sections 07-10

These four sections were cut from one working directory (`L:\AI-gene-drug\code`) whose scripts share a single numeric sequence and a single `results/` folder. The true run order runs *across* the four directories by numeric prefix, not within them:

```
   1  [08]  01_parse_risk_genes.py
   2  [08]  02_build_resources.py
  18  [07]  18_spatial_discover.py
  19  [07]  19_geomx_extract.py
  20  [07]  20_spatial_de.py
  21  [07]  21_cosmx_composition.py
  22  [07]  22_spatial_figures.py
  23  [07]  23_validated_set.py
  24  [07]  24_rds_reader.py
  25  [07]  25_cosmx_percell.py
  26  [07]  26_cosmx_gene_tests.py
  27  [07]  27_validated_priority.py
  28  [07]  28_validated_figures.py
  29  [07]  29_cosmx_plot_extract.py
  30  [07]  30_cosmx_gene_panels.py
  31  [07]  31_summary_statistics.py
  32  [09]  32_functional_parse_qc.py
  33  [09]  33_functional_stats.py
  34  [09]  34_functional_figures.py
  35  [09]  35_crispr_by_cellline.py
  36  [09]  36_crispr_cellline_figure.py
  37  [09]  37_sirna_by_cellline.py
  38  [09]  38_sirna_cellline_figure.py
  39  [09]  39_nontwas_figure.py
  40  [09]  40_ttest_diagnostics.py
  41  [09]  41_rps5_worked_example.py
  42  [09]  42_crispr_replicate_qc.py
  43  [09]  43_concordance_heatmap.py
  43  [09]  43_crispr_concordance_heatmap.py
  44  [09]  44_shared_genes_figure.py
  45  [09]  45_risk_increasing_figure.py
  46  [09]  46_risk_increasing_crispr.py
  47  [07]  47_cosmx_requested_genes.py
  48  [09]  48_screened_and_selected.py
  49  [10]  49_pathway_by_celltype.py
  50  [10]  50_pathway_figure.py
  51  [10]  51_pathway_merged.py
  52  [10]  52_pathway_merged_figures.py
  53  [10]  53_all_genes_etiology.py
  54  [10]  54_etiology_figures.py
```

Two gene lists are the spine of this block: `risk_genes_harmonized.csv` (written by `08_druggability_and_connectivity/gene_sets_and_resources/01_parse_risk_genes.py`, read by 15 scripts) and `validated_risk_genes_all.csv` (written by `07_orthogonal_support/spatial_geomx_cosmx/23_validated_set.py`, read by 15).

## Per-section detail

Each numbered directory carries a `README.md` with its scripts in run order, its resolved inputs (separated into files from other sections and files from outside the deposit), and its outputs with the downstream scripts that consume each one.

| Section | Scripts | Inputs resolved | Outputs resolved |
|---|---|---|---|
| [`01_single_cell_and_xqtl/`](../01_single_cell_and_xqtl/README.md) | 51 | 21 | 85 |
| [`02_accessible_chromatin/`](../02_accessible_chromatin/README.md) | 19 | 2 | 6 |
| [`03_grex_model_training/`](../03_grex_model_training/README.md) | 0 | 0 | 0 |
| [`04_association_testing/`](../04_association_testing/README.md) | 0 | 0 | 0 |
| [`05_colocalization_smr_heidi/`](../05_colocalization_smr_heidi/README.md) | 7 | 0 | 0 |
| [`06_tf_network/`](../06_tf_network/README.md) | 23 | 11 | 32 |
| [`07_orthogonal_support/`](../07_orthogonal_support/README.md) | 35 | 45 | 60 |
| [`08_druggability_and_connectivity/`](../08_druggability_and_connectivity/README.md) | 19 | 41 | 41 |
| [`09_functional_screens/`](../09_functional_screens/README.md) | 18 | 27 | 30 |
| [`10_pathway_enrichment/`](../10_pathway_enrichment/README.md) | 6 | 21 | 24 |
| [`11_figures/`](../11_figures/README.md) | 12 | 2 | 28 |
| [`12_supplementary_tables/`](../12_supplementary_tables/README.md) | 3 | 109 | 1 |
| [`manuscript_text/`](../manuscript_text/README.md) | 1 | 0 | 0 |

Totals: 307 distinct output files and 218 distinct input files resolved by name across the 194 deposited scripts.

## How this file was built

Automatically generated 2026-10-19 from MANIFEST.csv, the Methods order in `workflow.md`, and the file reads and writes parsed out of the deposited scripts themselves. The parse reads only string literals and path variables that appear at a read or write call site, and matches files between sections on basename. Paths assembled from shell variables, paths built by string formatting at run time, and files handled by external tools the scripts invoke rather than open themselves (tensorQTL, PLINK, MACS, SMR, PEER) are invisible to it, so every count here is a lower bound.
