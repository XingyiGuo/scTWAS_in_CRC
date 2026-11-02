# 04 - Association testing

Each cell-type GReX model is tested against the European-ancestry and multi-ancestry CRC
GWAS meta-analyses with the summary-statistic framework of S-PrediXcan as implemented for
sTF-TWAS: the gene-level Z score is a weighted sum over model variants of the GWAS effect
estimates, weighted by the model weights and scaled by the ratio of the variant and
predicted-expression standard deviations. For the two normal-colon training cohorts,
gene-level P values are combined across cohorts with the Cauchy combination test (ACAT).
Significance is set within each cell type by Bonferroni correction over the number of
retained gene models.

## Scripts, in run order

| File | What it does |
|---|---|
| `asso_selectCT_addR2global.sh` | The S-PrediXcan runner. `sbatch asso_selectCT_addR2global.sh <NL\|ADSER> <CT>` -- runs the association for one cell type and appends `R2_global` and `r_squared_global` from the model DB `extra` table to each `*.TWAS.csv` by gene. |
| `Genes_CauchyMeta_NL_cohort_fixed.R` | Cauchy combination (ACAT) across the two normal-colon cohorts. `Rscript Genes_CauchyMeta_NL_cohort_fixed.R all <EUR\|Multi>`; missing P values are dropped before ACAT is called. |
| `Annotate_CauchyMeta_NL.py` | Annotates the Cauchy-meta summary across the 7 normal-colon cell types, applying the `R2_global > 0.01` model-performance filter. |
| `Select_representitive_genes_for_each_CellType.R` | Picks the representative gene per cell type from the significant EUR and multi-ancestry sets, for the normal-colon and polyp arms. |

`config/asso.config.txt` supplies the model, GWAS and output paths to the runner.

## Inputs and outputs

**Inputs.** Model weights from `03_grex_model_training/`; the CRC GWAS meta-analysis
summary statistics (GWAS Catalog -- see `data/GWAS_SS/README.md`).

**Outputs.** Per-cell-type association tables (`*.TWAS.csv`) and the Cauchy-meta
significant-gene tables.

**Consumed by.** `05_colocalization_smr_heidi/`, `06_tf_network/`, the Manhattan panels in
`11_figures/`, and
`08_druggability_and_connectivity/gene_sets_and_resources/01_parse_risk_genes.py`, which
builds `risk_genes_harmonized.csv` -- the gene list every downstream section depends on.

## Note on cluster paths

`Annotate_CauchyMeta_NL.py` and `Select_representitive_genes_for_each_CellType.R` carry
absolute ACCRE input paths at the top of the file; repoint `BASE_DIR` / `NL_FILE` /
`ADSER_FILE` before running elsewhere.
