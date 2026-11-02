# 08 — Druggability and compound connectivity

`gene_sets_and_resources/` harmonises the prioritized gene list from the association
results and builds the public resource tables. `druggability/` integrates DrugBank,
ChEMBL, TTD and Open Targets into the druggable-gene and compound tables, the KEGG
enrichment for druggable genes, the clinical-phase and indication annotation, and the
figure export and source data. `ilincs_connectivity/` scores the prespecified compound
panel against the iLINCS LIB5 large-intestine library, including the size-matched
empirical null built from 5,000 randomly sampled signature sets, the consensus check and
the panel figure.

## Position in the pipeline

```
04_association_testing
        |
        v
   >>> 08_druggability_and_connectivity <<<
        |
        v
09_functional_screens  10_pathway_enrichment  12_supplementary_tables
```

The arrows above are the logical order of the published workflow, taken from the Methods (`docs/workflow.md`). The file-level links in *From other sections* and *Consumed by* below are the ones actually observed in the deposited code.

## Scripts, in run order (19 files)

Ordered by sub-directory, then by the numeric prefix the script carries in its original working directory. Descriptions are the script's own header line; blank where the script has no header.

| # | File | Lang | What it does | Declared display item |
|---|---|---|---|---|
| 1 | `druggability/01_druggability_ChEMBL_DB_TTD.py` | Python | Druggability of the 347 scTWAS colorectal risk genes. |  |
| 2 | `druggability/02_druggable_support_table.py` | Python | Support table for the druggable scTWAS risk genes, in 43_druggable_genes.xlsx format. |  |
| 3 | `druggability/03_druggability_paragraph.py` | Python | Manuscript paragraph: drug-target relationships and clinically advanced compounds. |  |
| 4 | `druggability/04_advanced_pairs_indication_json.py` | Python | Radial-tree JSON for the clinically advanced drug-target pairs (Figure 6d). | Figure 6d |
| 5 | `druggability/06_export_figure_png_pdf.py` | Python | Export the radial drug tree from SVG to a 300-dpi PNG and a vector PDF. |  |
| 6 | `druggability/07_figure_sourcedata.py` | Python | Source data for CRC_scTWAS_CellType_Gene_Drug_Indication.{svg,png,pdf}. |  |
| 7 | `druggability/protein_indication_json.py` | Python |  |  |
| 8 | `gene_sets_and_resources/01_parse_risk_genes.py` | Python |  |  |
| 9 | `gene_sets_and_resources/02_build_resources.py` | Python |  |  |
| 10 | `ilincs_connectivity/04_ilincs_reversal.py` | Python | Do iLINCS perturbagen signatures reverse the cell-type TWAS Z profile? |  |
| 11 | `ilincs_connectivity/05_ilincs_reversal_maxstat.py` | Python | Corrected inference for the iLINCS reversal scan, plus perturbagen names. |  |
| 12 | `ilincs_connectivity/06_ilincs_cmap_connectivity.py` | Python | CMap weighted connectivity (WTCS/NCS) between cell-type TWAS profiles and iLINCS. |  |
| 13 | `ilincs_connectivity/07_ilincs_drug_panel.py` | Python | Do the CRC risk-gene drugs themselves reverse the TWAS profile? |  |
| 14 | `ilincs_connectivity/08_sigcom_consensus.py` | Python | Consensus reversal query against the full LINCS collection via SigCom LINCS. |  |
| 15 | `ilincs_connectivity/09_ms_ilincs_panel_text.py` | Python | Methods and Results text for Fig_ilincs_drug_panel panel a. |  |
| 16 | `ilincs_connectivity/10_fig6e_panel.py` | Python | Figure 6e - iLINCS connectivity of the pre-specified risk-gene drug panel. | Figure 6e |
| 17 | `ilincs_connectivity/11_ms_ilincs_NG_text.py` | Python | Nature Genetics-format Results, figure legend and Methods for Figure 6e. | Figure 6e |
| 18 | `ilincs_connectivity/12_edfig7_caption.py` | Python | Extended Data Fig. 7 caption and reader's explanation. | Extended Data Fig. 7; Figure 6e; Figure 6; Supplementary Table 22 |
| 19 | `ilincs_connectivity/13_ms_ilincs_methods_revised.py` | Python | Revised iLINCS connectivity Methods for the manuscript, with Q19/Q20 resolved. |  |

## Inputs

**From other sections of this deposit**

| File | Produced by |
|---|---|
| `scTWAS_347genes_multisource_DE_support.csv` | `07_orthogonal_support/multiomics_de_support/01a_sctwas_gene_level_direction_and_support.R` |

**Not produced anywhere in this deposit** (23) — cohort data, public resources, or output of a step that is not deposited:

- `43_druggable_genes.xlsx`
- `63_druggable_genes.xlsx`
- `63_druggable_genes_KEGG_enrichment.csv`
- `CEGv2.txt`
- `Fig_ilincs_drug_panel_sourcedata.csv`
- `HTAN_NL_WellcomeSanger_NL_7celltypes_peakSNP.EUR_Multi.CauchyMeta.sig.csv`
- `LIB5_signature_expr_logDiffExpr_large_intestine.csv`
- `NEGv1.txt`
- `Network_STRING_v10.txt`
- `Proteins_Drugs_20250321.xlsx`
- `VUMC_HTAN_ADSER_9CTs_peakSNP.EUR_Multi.sig.csv`
- `bindings-drugbank-gene.tsv`
- `drugbank.tsv`
- `edges.sif.gz`
- `geneinfo.tsv.gz`
- `genes.tsv`
- `handoff_panel_cids.json`
- `hetionet-v1.0-nodes.tsv`
- `meta_sigSearch_LIB_5_143374_20231227.csv`
- `pathways.tsv`
- `proteins.tsv`
- `scTWAS_347genes_COSMIC_CGC_v103.csv`
- `scTWAS_347genes_DepMap_largeintestine.csv`

## Outputs

| File | Written by | Consumed by |
|---|---|---|
| `63_druggable_genes_advanced_records.csv` | `druggability/04_advanced_pairs_indication_json.py` | `07_figure_sourcedata.py` |
| `63_druggable_genes_drug_records.csv` | `druggability/02_druggable_support_table.py` | `03_druggability_paragraph.py`, `07_ilincs_drug_panel.py`, `11_ms_ilincs_NG_text.py` |
| `63_druggable_genes_summary.csv` | `druggability/02_druggable_support_table.py` | `03_druggability_paragraph.py` |
| `Fig6e_ilincs_drug_panel_sourcedata.csv` | `ilincs_connectivity/10_fig6e_panel.py` | — (terminal) |
| `Fig_scTWAS_druggability_sourcedata.csv` | `druggability/01_druggability_ChEMBL_DB_TTD.py` | — (terminal) |
| `bindingdb.tsv` | `gene_sets_and_resources/02_build_resources.py` | — (terminal) |
| `compound_treats_disease.tsv` | `gene_sets_and_resources/02_build_resources.py` | — (terminal) |
| `crc_DaG.tsv` | `gene_sets_and_resources/02_build_resources.py` | — (terminal) |
| `crc_disease_signature.tsv` | `gene_sets_and_resources/02_build_resources.py` | — (terminal) |
| `drugbank_drugs.tsv` | `gene_sets_and_resources/02_build_resources.py` | — (terminal) |
| `drugbank_targets.tsv` | `gene_sets_and_resources/02_build_resources.py` | — (terminal) |
| `essentiality_refsets.json` | `gene_sets_and_resources/02_build_resources.py` | — (terminal) |
| `gene_celltype_summary.csv` | `gene_sets_and_resources/01_parse_risk_genes.py` | — (terminal) |
| `gene_id_map.json` | `gene_sets_and_resources/02_build_resources.py` | `49_pathway_by_celltype.py`, `51_pathway_merged.py`, `53_all_genes_etiology.py` |
| `hetionet_CbG.tsv` | `gene_sets_and_resources/02_build_resources.py` | — (terminal) |
| `ilincs_all_compounds_mean_NCS.csv` | `ilincs_connectivity/07_ilincs_drug_panel.py` | `09_ms_ilincs_panel_text.py`, `11_ms_ilincs_NG_text.py` |
| `ilincs_cmap_all_signatures.csv` | `ilincs_connectivity/06_ilincs_cmap_connectivity.py` | `07_ilincs_drug_panel.py`, `09_ms_ilincs_panel_text.py`, `11_ms_ilincs_NG_text.py` |
| `ilincs_cmap_model_summary.csv` | `ilincs_connectivity/06_ilincs_cmap_connectivity.py` | `09_ms_ilincs_panel_text.py`, `11_ms_ilincs_NG_text.py` |
| `ilincs_cmap_top20_per_model.csv` | `ilincs_connectivity/06_ilincs_cmap_connectivity.py` | — (terminal) |
| `ilincs_cmap_vs_pearson.csv` | `ilincs_connectivity/06_ilincs_cmap_connectivity.py` | — (terminal) |
| `ilincs_drug_panel_matching.csv` | `ilincs_connectivity/07_ilincs_drug_panel.py` | `09_ms_ilincs_panel_text.py` |
| `ilincs_drug_panel_per_compound.csv` | `ilincs_connectivity/07_ilincs_drug_panel.py` | `09_ms_ilincs_panel_text.py`, `11_ms_ilincs_NG_text.py` |
| `ilincs_drug_panel_per_model.csv` | `ilincs_connectivity/07_ilincs_drug_panel.py` | — (terminal) |
| `ilincs_drug_panel_summary.csv` | `ilincs_connectivity/07_ilincs_drug_panel.py` | `09_ms_ilincs_panel_text.py` |
| `ilincs_query_genes.csv` | `ilincs_connectivity/04_ilincs_reversal.py` | `05_ilincs_reversal_maxstat.py`, `06_ilincs_cmap_connectivity.py`, `08_sigcom_consensus.py` (+1) |
| `ilincs_query_model_inventory.csv` | `ilincs_connectivity/11_ms_ilincs_NG_text.py` | `13_ms_ilincs_methods_revised.py` |
| `ilincs_reversal_all_signatures.csv` | `ilincs_connectivity/04_ilincs_reversal.py` | `05_ilincs_reversal_maxstat.py`, `06_ilincs_cmap_connectivity.py` |
| `ilincs_reversal_all_signatures_annotated.csv` | `ilincs_connectivity/05_ilincs_reversal_maxstat.py` | — (terminal) |
| `ilincs_reversal_maxstat_null.csv` | `ilincs_connectivity/05_ilincs_reversal_maxstat.py` | — (terminal) |
| `ilincs_reversal_model_summary.csv` | `ilincs_connectivity/04_ilincs_reversal.py` | `09_ms_ilincs_panel_text.py` |
| `ilincs_reversal_top20_per_model.csv` | `ilincs_connectivity/05_ilincs_reversal_maxstat.py` | — (terminal) |
| `ilincs_reversal_top_signatures.csv` | `ilincs_connectivity/04_ilincs_reversal.py` | — (terminal) |
| `ms_ilincs_panel_NG_numbers.csv` | `ilincs_connectivity/11_ms_ilincs_NG_text.py` | `12_edfig7_caption.py`, `13_ms_ilincs_methods_revised.py` |
| `ms_ilincs_panel_numbers.csv` | `ilincs_connectivity/09_ms_ilincs_panel_text.py` | `11_ms_ilincs_NG_text.py` |
| `risk_genes_harmonized.csv` | `gene_sets_and_resources/01_parse_risk_genes.py` | `19_geomx_extract.py`, `20_spatial_de.py`, `21_cosmx_composition.py` (+12) |
| `scTWAS_347genes_druggability_ChEMBL_DB_TTD.csv` | `druggability/01_druggability_ChEMBL_DB_TTD.py` | `02_druggable_support_table.py` |
| `scTWAS_druggability_tier_summary.csv` | `druggability/01_druggability_ChEMBL_DB_TTD.py` | — (terminal) |
| `scTWAS_druggable_gene_drug_pairs.csv` | `druggability/01_druggability_ChEMBL_DB_TTD.py` | — (terminal) |
| `sigcom_model_summary.csv` | `ilincs_connectivity/08_sigcom_consensus.py` | `09_ms_ilincs_panel_text.py` |
| `sigcom_recurrent_reversers.csv` | `ilincs_connectivity/08_sigcom_consensus.py` | — (terminal) |
| `sigcom_reversers_top.csv` | `ilincs_connectivity/08_sigcom_consensus.py` | — (terminal) |

## How this file was built

Automatically generated 2026-10-19 from MANIFEST.csv, the Methods order in `docs/workflow.md`, and the file reads and writes parsed out of the deposited scripts themselves. The extraction reads only string literals and path variables that appear at a read or write call site, so three classes of file are invisible to it: paths assembled from shell variables (most `.sh` and `.sbatch` wrappers), paths built by string formatting at run time, and files handled by external tools that the scripts invoke rather than open themselves (tensorQTL, PLINK, MACS, SMR, PEER). The lists above are therefore a lower bound, not a complete inventory. `MANIFEST.csv` at the repository root is authoritative for the file list itself.
