# Druggability and compound connectivity

Annotate scTWAS risk genes with drug–target evidence and test whether
clinically advanced compounds reverse their expression profiles in iLINCS.
Script numbering is local to each folder. Check the input and output paths at
the top of each script. Use the [Python dependencies](../../env/python_requirements.txt)
with pandas 2.2 or newer.

## Drug–target annotation

1. Run [01_druggability_ChEMBL_DB_TTD.py](druggability/01_druggability_ChEMBL_DB_TTD.py)
   with the `ChEMBL_DB_TTD` sheet of `Proteins_Drugs_20250321.xlsx`,
   [section 07's](../07_orthogonal_support/README.md) gene-level
   `scTWAS_347genes_multisource_DE_support.csv`, and the staged
   `scTWAS_347genes_DepMap_largeintestine.csv`. The script falls back to the
   DepMap export in `Experiments/Archive`. It assigns development tiers and
   writes `scTWAS_347genes_druggability_ChEMBL_DB_TTD.csv`, gene–drug pairs,
   and a tier summary under `Drugs/table/`.
2. Run [02_druggable_support_table.py](druggability/02_druggable_support_table.py)
   with the gene-level export from step 01 and both workbook sheets
   (`ChEMBL_DB_TTD` and `OpenTargets`). It reads the original
   `43_druggable_genes.xlsx` or the archived membership audit
   `audit_43_druggable_genes_vs_ChEMBL_DB_TTD.csv` for the legacy gene flag.
   It writes `Drugs/63_druggable_genes.xlsx` and CSVs under `Drugs/table/`.
   The drug records supply Table 22 and contain 303 distinct drug names
   after trimming and case normalization. The workbook feeds Figure 6,
   and `63_druggable_genes_drug_records.csv` feeds the iLINCS panel.

## Enrichr pathway analysis

[03_enrich_druggable_genes_enrichr.R](druggability/03_enrich_druggable_genes_enrichr.R)
submits the unique gene symbols from `63_druggable_genes_summary.csv` to
Enrichr. Install the R package `enrichR`. The default library is
`KEGG_2021_Human`, and an optional third argument selects another library.

```powershell
Rscript code/08_druggability_and_connectivity/druggability/03_enrich_druggable_genes_enrichr.R `
  'W:/VUMC_projects/scTWAS/Drugs/table/63_druggable_genes_summary.csv' `
  results/63_druggable_genes_Enrichr.tsv
```

The example writes all Enrichr results to `63_druggable_genes_Enrichr.tsv`
and pathways with Enrichr adjusted P < 0.05 to
`63_druggable_genes_Enrichr_significant.tsv`. Both files retain the statistics
returned by Enrichr, including its adjusted P values and overlapping genes.

The current Supplementary Table 23 matches the archived KEGG REST export
`W:/VUMC_projects/scTWAS/Drugs/table/63_druggable_genes_KEGG_enrichment.csv`.
The original script remains in the project archive at
`W:/VUMC_projects/scTWAS/ms-guo/res/Drugs/src/02c_kegg_enrichment.py`.
Using the Enrichr results in the manuscript requires updating Table 23 and
the Methods to match the Enrichr run.

## iLINCS compound connectivity

Stage the cohort-level `*_peakSNP.EUR.TWAS.csv` files,
`LIB5_signature_expr_logDiffExpr_large_intestine.csv`, and
`meta_sigSearch_LIB_5_143374_20231227.csv`. Run the following three scripts
in order. Outputs go under `iLINCS_reversal/table/`.

1. [01_scan_ilincs_reversal.py](ilincs_connectivity/01_scan_ilincs_reversal.py)
   selects the 50 highest and 50 lowest TWAS Z scores among shared landmark
   genes for each model with at least 200 shared genes. It writes
   `ilincs_query_genes.csv` and `ilincs_reversal_all_signatures.csv`.
   The reported run retained 18 of 23 cohort-level models.
2. [03_score_cmap_connectivity.py](ilincs_connectivity/03_score_cmap_connectivity.py)
   uses both exports to calculate weighted connectivity (WTCS) and
   cell-line-normalized connectivity (NCS). It writes
   `ilincs_cmap_all_signatures.csv` for the compound panel.
3. [04_test_ilincs_drug_panel.py](ilincs_connectivity/04_test_ilincs_drug_panel.py)
   selects approved or phase II/III compounds from the drug records. It uses
   PubChem CIDs to match compounds and collapse aliases. Pass `--cids` with
   `handoff_panel_cids.json` or a staged copy of the archived
   `ilincs_drug_panel_matching.csv`. For each matched compound, the script
   averages NCS across models and then signatures, compares it with 5,000
   signature-count-matched draws with replacement, and applies BH correction
   across the matched compounds. The archived matching export contains 184
   panel entries after alias collapse, of which 43 are in the library. The
   `ilincs_drug_panel_per_compound.csv` supplies Table 24 and Extended Data Figure 9.

[02_test_reversal_maxstat.py](ilincs_connectivity/02_test_reversal_maxstat.py)
adds a Pearson max-statistic test and a Spearman sensitivity check.
[05_check_sigcom_consensus.py](ilincs_connectivity/05_check_sigcom_consensus.py)
queries the live SigCom LINCS libraries. Both use exports from step 01 and
run separately from the compound panel.

See [section 11](../11_figures/README.md)
for figure scripts and the
[supplementary workbook](../manuscript_text/Li_etal_Supplementary_Tables_11062026.xlsx)
for the final tables.
