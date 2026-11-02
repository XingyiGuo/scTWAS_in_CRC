# 03 - GReX model training

Per-gene elastic-net models (alpha = 0.5) of genetically regulated expression, trained for
each cell type over variants that lie inside the cell-type accessible chromatin and within
1 Mb of the transcription start site, under the sTF-TWAS framework and following PrediXcan.
Predictive performance is the coefficient of determination between predicted and observed
expression on held-out samples across three repeated tenfold cross-validation runs; models
meeting the retention threshold are carried forward.

## Scripts, in run order

| File | What it does |
|---|---|
| `TWAS_modelbuilding.sh` | Slurm wrapper. `sbatch TWAS_modelbuilding.sh <CT_index> <chrom> <config>` -- one job per cell type and chromosome. |
| `TWAS_modelbuilding_scTWAS.R` | The trainer. Reads residualised pseudobulk expression and peak-restricted genotypes, fits the per-gene elastic net with `glmnet`, and writes per-gene weights and cross-validation summaries. Uses balanced non-empty CV folds, which matters for the small ASC and SSC cohorts. |
| `WeightsSummaryToDB.R` | Collects the per-gene weights and summaries into the PredictDB-format SQLite `.db` used downstream. Retains `R2_global` and `r_squared_global` in the DB `extra` table. |

## Configuration

`config/` holds one file per training cohort; it is the third argument to the trainer and
the second to `WeightsSummaryToDB.R`.

| File | Cohort |
|---|---|
| `config/HTAN_GSE266546_peakSNP.ADSER.config.txt` | VUMC HTAN polyps (adenoma / serrated), 9 cell types |
| `config/HTAN_GSE266546_peakSNP.NL.config.txt` | VUMC HTAN + GSE266546 normal mucosa, 7 cell types |
| `config/WellcomeSanger_peakSNP.NL.config.txt` | Wellcome Sanger Gut Cell Atlas normal mucosa, 7 cell types |

## Inputs and outputs

**Inputs.** Residualised per-donor pseudobulk expression per cell type from
`01_single_cell_and_xqtl/`; the merged cell-type peak sets from `02_accessible_chromatin/`,
which are also shipped in `data/scATAC-seq/used_peaks/`; cohort genotypes.

**Outputs.** 23 per-cell-type `*_peakSNP.db` model databases, released through Zenodo and
documented in `scTWAS_models/README.md`.

**Consumed by.** `04_association_testing/` (model weights);
`06_tf_network/correlation_and_occupancy/1a_Extract_sigTWASGenes_SNPs_gtf.R`, which opens
the `.db` files directly via RSQLite; and
`11_figures/model_performance/Fig_peaks_R2_violin.R`, which reads the model summaries.

## Note on cluster paths

The wrapper and configs carry the absolute ACCRE paths the analysis ran under. They record
what was executed and need repointing to a local layout.
