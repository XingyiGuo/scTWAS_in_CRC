# 06 — Transcription-factor-target networks

`correlation_and_occupancy/` extracts the model SNPs for significant TWAS genes from the
per-cell-type model databases and performs the one-pass GWAS rsID lookup (`1a`), computes
TF-target pseudobulk correlations for the polyp (9 cell types) and normal (7 cell types)
model sets (`2a`), and tests lead-variant TF occupancy against target-gene correlation
(`2b`). `network_figures/` builds the matrices, orderings and network panels and
reconstructs per-panel source data. `pathway_integration/` joins the network with the
pathway enrichment results and tests cell-type specificity.

## Position in the pipeline

```
03_grex_model_training  04_association_testing  01_single_cell_and_xqtl
        |
        v
   >>> 06_tf_network <<<
        |
        v
   (terminal: display items / supplementary tables)
```

The arrows above are the logical order of the published workflow, taken from the Methods (`docs/workflow.md`). The file-level links in *From other sections* and *Consumed by* below are the ones actually observed in the deposited code.

## Scripts, in run order (23 files)

Ordered by sub-directory, then by the numeric prefix the script carries in its original working directory. Descriptions are the script's own header line; blank where the script has no header.

| # | File | Lang | What it does | Declared display item |
|---|---|---|---|---|
| 1 | `correlation_and_occupancy/1a_Extract_sigTWASGenes_SNPs_gtf.R` | R | Fast extraction of model SNPs for significant TWAS genes from either ADSER models or normal-colon (NL) models, followed by a one-pass GWAS RSID lookup. |  |
| 2 | `correlation_and_occupancy/2a_HTAN_ADSER_TWASgene_TFs_pseudobulk_correlation_ADSER.R` | R | HTAN ADSER: one cell-type-specific TF <-> TWAS-risk-gene pseudobulk correlation Biological question Within each epithelial cell type, are the 51 CRC-r. |  |
| 3 | `correlation_and_occupancy/2a_WellcomeSanger_NL_TWASgene_TFs_pseudobulk_correlation_NL.R` | R | WellcomeSanger healthy/NL: cell-type-specific CRC TF <-> NL TWAS-risk-gene pseudobulk correlations Biological question Within each healthy epithelial . |  |
| 4 | `correlation_and_occupancy/2b_HTAN_ADSER_TWASgene_leadVariant_TFs_occupancy_vs_geneCorrelation_9CT.R` | R | ADSER, 9 cell types. 1) Every model-build SNP of every significant TWAS gene is read from the cell-type specific prediction weights (one weight… |  |
| 5 | `correlation_and_occupancy/2b_WellcomeSanger_NL_TWASgene_leadVariant_TFs_occupancy_vs_geneCorrelation_7CT.R` | R | Normal colon (NL), 7 cell types. 1) Every model-build SNP of every significant TWAS gene is read from the cell-type specific prediction weights… |  |
| 6 | `network_figures/LQ_review/build_sourcedata.py` | Python | Reconstruct per-panel source data for Figures 1-3 and Supplementary Figure 1 | Figure 1; Supplementary Figure 1; Figure 2; Figure 3 |
| 7 | `network_figures/code/barfull.py` | Python |  |  |
| 8 | `network_figures/code/bars2.py` | Python |  |  |
| 9 | `network_figures/code/build_mats.py` | Python |  |  |
| 10 | `network_figures/code/check_order.py` | Python |  |  |
| 11 | `network_figures/code/extract2.py` | Python |  |  |
| 12 | `network_figures/code/makefig1b.py` | Python |  |  |
| 13 | `network_figures/code/makefig2.py` | Python |  |  |
| 14 | `network_figures/code/makefig3.py` | Python |  |  |
| 15 | `network_figures/code/makefig4.py` | Python |  |  |
| 16 | `network_figures/code/netfig2.py` | Python |  |  |
| 17 | `network_figures/code/orders.py` | Python |  |  |
| 18 | `network_figures/code/slots.py` | Python |  |  |
| 19 | `network_figures/code/spec.py` | Python |  |  |
| 20 | `network_figures/code/stats.py` | Python |  |  |
| 21 | `network_figures/code/tidy.py` | Python |  |  |
| 22 | `pathway_integration/celltype_specificity_test.py` | Python | Cell-type specificity of the 546 scTWAS associations (Results paragraph 1, Figure 4A). | Figure 4A |
| 23 | `pathway_integration/tf_pathway_integration.py` | Python | Integration of the 51-TF risk-gene network with the scTWAS EnrichR pathway results. | Figure 4e,f |

## Inputs

**Not produced anywhere in this deposit** (7) — cohort data, public resources, or output of a step that is not deposited:

- `ADSER_lead_variant_per_celltype_gene_with_occupying_TFs.csv`
- `NL_ADSER_cell_genes_EnrichR.xlsx`
- `NL_lead_variant_per_celltype_gene_with_occupying_TFs.csv`
- `S1_TF_gene_correlations_all.csv`
- `S2_TF_leadvariant_occupancy.csv`
- `scTWAS_ADSER_NL_collected_results.csv`
- `tf_leadvariant_occupancy.csv`

## Outputs

| File | Written by | Consumed by |
|---|---|---|
| `ALL7CT_WellcomeSanger_NL_CRC51TF_x_TWASriskGene_FDR005.csv` | `correlation_and_occupancy/2a_WellcomeSanger_NL_TWASgene_TFs_pseudobulk_correlation_NL.R` | — (terminal) |
| `ALL7CT_WellcomeSanger_NL_CRC51TF_x_TWASriskGene_correlation_long.csv.gz` | `correlation_and_occupancy/2a_WellcomeSanger_NL_TWASgene_TFs_pseudobulk_correlation_NL.R` | — (terminal) |
| `ALL7CT_WellcomeSanger_NL_TF_correlation_summary.csv` | `correlation_and_occupancy/2a_WellcomeSanger_NL_TWASgene_TFs_pseudobulk_correlation_NL.R` | — (terminal) |
| `ALL7CT_WellcomeSanger_NL_correlation_QC_summary.csv` | `correlation_and_occupancy/2a_WellcomeSanger_NL_TWASgene_TFs_pseudobulk_correlation_NL.R` | — (terminal) |
| `ALL9CT_oneCorrelation_CRC51TF_x_TWASriskGene_FDR005.csv` | `correlation_and_occupancy/2a_HTAN_ADSER_TWASgene_TFs_pseudobulk_correlation_ADSER.R` | — (terminal) |
| `ALL9CT_oneCorrelation_CRC51TF_x_TWASriskGene_long.csv.gz` | `correlation_and_occupancy/2a_HTAN_ADSER_TWASgene_TFs_pseudobulk_correlation_ADSER.R` | — (terminal) |
| `ALL9CT_oneCorrelation_QC_summary.csv` | `correlation_and_occupancy/2a_HTAN_ADSER_TWASgene_TFs_pseudobulk_correlation_ADSER.R` | — (terminal) |
| `CRC51TFs_used.csv` | `correlation_and_occupancy/2a_WellcomeSanger_NL_TWASgene_TFs_pseudobulk_correlation_NL.R` | — (terminal) |
| `CRC51TFs_used_alphabetical.csv` | `correlation_and_occupancy/2a_HTAN_ADSER_TWASgene_TFs_pseudobulk_correlation_ADSER.R` | — (terminal) |
| `Figure4FG_sourcedata.csv` | `pathway_integration/tf_pathway_integration.py` | — (terminal) |
| `T2_pathway_gene_TF_evidence.csv` | `pathway_integration/tf_pathway_integration.py` | — (terminal) |
| `T3_pathway_modules_collapsed.csv` | `pathway_integration/tf_pathway_integration.py` | — (terminal) |
| `T4_TF_pathway_pairs_hypergeom.csv` | `pathway_integration/tf_pathway_integration.py` | — (terminal) |
| `T7_associations_by_model.csv` | `pathway_integration/celltype_specificity_test.py` | — (terminal) |
| `TF_counts_genes_with_occupied_lead_variant.csv` | `correlation_and_occupancy/2b_WellcomeSanger_NL_TWASgene_leadVariant_TFs_occupancy_vs_geneCorrelation_7CT.R` | — (terminal) |
| `TF_leadVariantOccupation_vs_geneCorrelation_plot_long.csv` | `correlation_and_occupancy/2b_WellcomeSanger_NL_TWASgene_leadVariant_TFs_occupancy_vs_geneCorrelation_7CT.R` | — (terminal) |
| `TF_leadVariantOccupation_vs_geneCorrelation_wide.csv` | `correlation_and_occupancy/2b_WellcomeSanger_NL_TWASgene_leadVariant_TFs_occupancy_vs_geneCorrelation_7CT.R` | — (terminal) |
| `_CRC51TF_feature_presence_alphabetical.csv` | `correlation_and_occupancy/2a_HTAN_ADSER_TWASgene_TFs_pseudobulk_correlation_ADSER.R` | — (terminal) |
| `_TF_counts_genes_with_occupied_lead_variant.csv` | `correlation_and_occupancy/2b_HTAN_ADSER_TWASgene_leadVariant_TFs_occupancy_vs_geneCorrelation_9CT.R` | — (terminal) |
| `_TF_leadVariantOccupation_vs_geneCorrelation_plot_long.csv` | `correlation_and_occupancy/2b_HTAN_ADSER_TWASgene_leadVariant_TFs_occupancy_vs_geneCorrelation_9CT.R` | — (terminal) |
| `_TF_leadVariantOccupation_vs_geneCorrelation_wide.csv` | `correlation_and_occupancy/2b_HTAN_ADSER_TWASgene_leadVariant_TFs_occupancy_vs_geneCorrelation_9CT.R` | — (terminal) |
| `_TWAS_risk_genes_used_alphabetical.csv` | `correlation_and_occupancy/2a_HTAN_ADSER_TWASgene_TFs_pseudobulk_correlation_ADSER.R` | — (terminal) |
| `_TWASriskGene_feature_presence_alphabetical.csv` | `correlation_and_occupancy/2a_HTAN_ADSER_TWASgene_TFs_pseudobulk_correlation_ADSER.R` | — (terminal) |
| `_lead_variant_per_celltype_gene_with_occupying_TFs.csv` | `correlation_and_occupancy/2b_HTAN_ADSER_TWASgene_leadVariant_TFs_occupancy_vs_geneCorrelation_9CT.R` | — (terminal) |
| `_modelSNPs_allGWASinfo_smallestP_flagged.csv` | `correlation_and_occupancy/2b_HTAN_ADSER_TWASgene_leadVariant_TFs_occupancy_vs_geneCorrelation_9CT.R` | — (terminal) |
| `bar_values.json` | `network_figures/code/barfull.py` | — (terminal) |
| `enrichment_stats.csv` | `network_figures/code/stats.py` | — (terminal) |
| `file_meta.csv` | `network_figures/code/tidy.py` | `stats.py` |
| `key_tfs.csv` | `network_figures/code/makefig1b.py` | — (terminal) |
| `slot_values.json` | `network_figures/code/slots.py` | `check_order.py` |
| `tf_gene_correlations.csv` | `network_figures/code/tidy.py` | `barfull.py`, `check_order.py`, `makefig1b.py` (+2) |
| `tf_gene_correlations_annot.csv` | `network_figures/code/stats.py` | `makefig4.py`, `spec.py` |

## How this file was built

Automatically generated 2026-10-19 from MANIFEST.csv, the Methods order in `docs/workflow.md`, and the file reads and writes parsed out of the deposited scripts themselves. The extraction reads only string literals and path variables that appear at a read or write call site, so three classes of file are invisible to it: paths assembled from shell variables (most `.sh` and `.sbatch` wrappers), paths built by string formatting at run time, and files handled by external tools that the scripts invoke rather than open themselves (tensorQTL, PLINK, MACS, SMR, PEER). The lists above are therefore a lower bound, not a complete inventory. `MANIFEST.csv` at the repository root is authoritative for the file list itself.

## Added in this release

These files were merged in from the project repository after the first deposit was assembled, so they are not described in the run order above.

Two TF figure scripts: the JASPAR2024 motif logos for the 51 CRC TFs, and the gene-SNP-TF Cistrome panel.

- `network_figures/CRC51TF_JASPAR2024_motif_logos.R`
- `network_figures/Plot_Gene_SNP_TF_Cistrome_ggplot2_updated_36x14_TSS_TFcategory_topRSlabels.R`
