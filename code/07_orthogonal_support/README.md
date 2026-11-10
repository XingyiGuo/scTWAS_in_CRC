# Orthogonal molecular support

Evaluate scTWAS risk genes using independent colorectal cancer expression and
protein data. Script numbering is local to each folder. Check each script's
configured paths before running it. Run scripts with relative paths from the
analysis workspace that contains their inputs.

## Shared risk-gene list

Run [01_parse_risk_genes.py](01_parse_risk_genes.py) with the significant
normal-colon and polyp scTWAS tables in `data/`. It writes
`results/risk_genes_harmonized.csv` for the spatial analyses and
`results/gene_celltype_summary.csv` for the gene-by-cell-type summary.

## scRNA-seq differential expression

1. In [pseudobulk_de](pseudobulk_de), set `ROOT` and `UNITS_CSV` in script 01
   to the four COLON MAP H5AD files' directory and the original
   `pseudobulk_units.csv`. Run script 01 to create `pseudobulk_raw.npz`.
2. Run [04_refit_expanded_late_contrasts.py](pseudobulk_de/04_refit_expanded_late_contrasts.py)
   with `--units` pointing to the same CSV and `--counts` to the NPZ.
   It uses the fitting and contrast definitions from scripts 02–03 to export
   all 32 C/E/L contrasts to `data/de_C_E_L_contrasts_ASC75_SSC46.csv`.
   The manuscript specifies PyDESeq2 0.5.4.
3. In [multiomics_de_support](multiomics_de_support), scripts 01 and 02 use
   that DE export and `data/scTWAS_ADSER_NL_collected_results.csv` to produce
   the per-association and ASC/SSC cross-cell-type support for Supplementary
   Tables 18 and 19. They prefer a matching Parquet file when Arrow is
   available. Script 03 counts direction agreement from the staged support
   CSV or the output of script 01. Set `SC_TWAS_SUPPORT_ROOT` for script 03
   if using a different workspace.

## TCGA, CPTAC and gene-level support

In [multiomics_de_support](multiomics_de_support), scripts 04 and 05
analyze TCGA-COAD RNA and CPTAC-COAD protein abundance, respectively. Scripts
06 and 07 run separate CPTAC analyses using a proteome-wide paired test and
limma. The Python CPTAC scripts used Python 3.11, `cptac==1.5.14`, and
`setuptools==80.8.0`.

[08_summarize_gene_level_support.R](multiomics_de_support/08_summarize_gene_level_support.R)
reads the collected scTWAS CSV, `data/risk_gene_all_comparisons.csv` from
the GeoMx workflow below, and a separately staged `data/TCGA_CPTAC_support.csv`
with its column names on row 2. The extended export also requires the two
published Visium DE CSVs at the paths configured in the script. It writes
`table/scTWAS_347genes_multisource_DE_support.csv` and accompanying summaries.

## GeoMx and CosMx

Run the following scripts in [spatial_geomx_cosmx](spatial_geomx_cosmx).

- **GeoMx (01, 02, 08).** Set `SCTWAS_GEOMX_ROOT` to the extracted directory
  containing `exp1/` and `exp2/`, `SCTWAS_GEOMX_OUT` to the workspace's
  `spatial/` directory, and `SCTWAS_RISK_GENES` to the shared risk-gene CSV.
  Script 01 stages `de_exp1.tsv.gz` (CTA) and `de_exp2_epi.tsv.gz` and
  `de_exp2_str.tsv.gz` (WTA). Script 02 tests support using the published
  contrasts and writes `results/spatial_risk_gene_DE.csv` and the validation
  tables. Use `--validation-only` with the existing DE and shared risk-gene
  CSVs to refresh validation. Run script 08 after 02 and stage its
  `results/risk_gene_all_comparisons.csv` under `data/` for multiomics script 08.
- **CosMx expression (04–06).** Run script 04 with `map` and then `genes` on
  the Seurat RDS. Use a list with one risk gene per line and output prefix
  `cx` to create `cx_counts.npz` and `cx_cells.txt.gz`. Script 05 reads these
  files plus barcode-matched `meta_sel.tsv.gz` from its configured `CX`
  directory. It uses the GeoMx `validated_risk_genes.csv` list to choose
  plotting genes. Stage its `cx_pseudobulk.csv.gz` under `spatial/`. Script 06
  reads it and `results/validated_risk_genes_all.csv` and writes the paired
  Wilcoxon results to `results/cosmx_crc_vs_normal_epi.csv`. It requires at
  least 50 cells per patient and compartment and six eligible pairs per gene.
- **Composition and follow-up (03, 07, 09).** Script 03 uses staged CosMx
  composition tables and cell metadata. Script 07 joins the earlier GeoMx
  validation set to integrated gene priorities and the paired CosMx results.
  Script 09 uses the paired results and staged per-cell expression archives
  to produce follow-up plots for requested genes.
