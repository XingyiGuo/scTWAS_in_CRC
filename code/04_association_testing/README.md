# Association testing

Test cell-type GReX models against the European-ancestry (EUR) and
multi-ancestry (Multi) colorectal cancer GWAS.

## Inputs

Edit [`config/asso.config.txt`](config/asso.config.txt) to point to MetaXcan,
the [model databases and covariances](../03_grex_model_training), and the
[chromosome-level GWAS files](../../data/GWAS_SS). Set `ASSO_OUTPUT_DIR` to
share a results directory across the scripts. R and Python dependencies are
listed in [`env/`](../../env).

## Run order

1. **S-PrediXcan.** From this directory, submit one job per cell type and model
   set. Each job tests both ancestries and writes `*.TWAS.csv` with model
   performance fields.

   ```bash
   sbatch asso_selectCT_addR2global.sh NL ABS
   MODEL_DATASET_PREFIX=WellcomeSanger sbatch asso_selectCT_addR2global.sh NL ABS
   sbatch asso_selectCT_addR2global.sh ADSER ASC
   ```

2. **Normal-mucosa meta-analysis.** Run the following for each ancestry.
   The R script keeps models with `R2_global > 0.01` and combines the two
   normal-mucosa cohorts using ACAT. The Python script applies per-cell-type
   Bonferroni thresholds and writes annotated tables and counts.

   ```bash
   for ancestry in EUR Multi; do
       Rscript Genes_CauchyMeta_NL_cohort_fixed.R all "$ancestry"
       python Annotate_CauchyMeta_NL.py "$ancestry"
   done
   ```

3. **Collect significant results.** Supply the preassembled EUR/Multi
   `*.sig.csv` tables through `NL_FILE` and `ADSER_FILE`, then run
   `Rscript Select_representitive_genes_for_each_CellType.R`. It exports
   one Z score and P value per input gene–cell-type pair.
