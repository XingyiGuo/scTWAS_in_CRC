# GReX model training

This section trains per-gene elastic-net models of genetically regulated expression
for each cell type. The trainer uses SNPs in the cell type's peak rsID list and
within 1 Mb of the gene's strand-specific transcription start site. The exporter
combines chromosome outputs into the PredictDB model and covariance file used by
[S-PrediXcan](../04_association_testing).

## Choose a cohort

| Config | Cohort | Cell types | Minimum matched samples |
|---|---|---:|---:|
| [HTAN_GSE266546_peakSNP.ADSER.config.txt](config/HTAN_GSE266546_peakSNP.ADSER.config.txt) | VUMC HTAN COLONMAP polyps | 9 | 30 |
| [HTAN_GSE266546_peakSNP.NL.config.txt](config/HTAN_GSE266546_peakSNP.NL.config.txt) | VUMC HTAN and GSE266546 normal mucosa | 7 | 80 |
| [WellcomeSanger_peakSNP.NL.config.txt](config/WellcomeSanger_peakSNP.NL.config.txt) | Wellcome Sanger healthy tissue | 7 | 250 |

The first argument to both R scripts is the **1-based index** in the config's
`cell_types` vector. Index 1 selects ABS in all three supplied configs. Each
config covers autosomes 1–22 and uses `EUR` peak rsID lists, a minimum minor
allele frequency of 0.05, at least two cis SNPs, and elastic-net alpha of 0.5.

## Stage the inputs

Use the R packages in [env](../../env) and update the cohort config's cluster
paths for your data layout. The trainer uses the bundled
[`GenoPhenoBasic.R`](../lib/GenoPhenoBasic.R) by default. Set
`TWAS_CONFIG$paths$helper_functions` to load a specific copy.

- **Expression and sample mapping.** Use gene-by-sample PEER residuals and the
  corresponding expression-to-genotype ID mapping from
  [section 01](../01_single_cell_and_xqtl). The trainer aligns donors in
  mapping-file order. COLONMAP uses a separate mapping file per cell type,
  and Wellcome Sanger uses `individual_id` on both sides.
- **Genotypes and annotations.** Supply one PLINK `.traw` dosage table and
  matching master variant annotation per chromosome, plus the GENCODE GTF.
  The supplied configs align dosages to the annotation's effect allele and
  remove variants with incompatible alleles or missing dosages before MAF filtering.
- **Peak rsIDs.** Supply the precomputed cell-type lists named by the config,
  such as `EUR.NL_ABS.rsid.txt.gz`. The trainer uses the first column as
  rsIDs. [Section 02](../02_accessible_chromatin) provides the peak BEDs,
  and the script that turns those intervals into the configured rsID lists
  is external to this repository.

## Train and export

From this directory, run one job per cell type and chromosome. For example,
these commands train chromosome 1 for ABS in the normal HTAN/GSE266546 cohort
and export the ABS model after all configured chromosome jobs have finished:

```bash
cfg=config/HTAN_GSE266546_peakSNP.NL.config.txt
sbatch TWAS_modelbuilding.sh 1 1 "$cfg"
# After all ABS chromosome jobs complete:
Rscript WeightsSummaryToDB.R 1 "$cfg"
```

For a direct run, use `Rscript TWAS_modelbuilding_scTWAS.R <CT_index> <chrom> <config>`.
To submit the Slurm wrapper from elsewhere, set `TWAS_SCRIPT_DIR` to the
training script directory. Review the Slurm account and resource settings.

The trainer selects lambda by the mean cross-validation MSE across three
balanced tenfold repeats. It calculates the reported held-out prediction
metrics from the first repeat at the selected lambda. Each chromosome job
writes model summaries, weights, SNP covariances and sample-matching audits
under the configured output directory.

The exporter reads those chromosome summaries, weights and covariances.
Check any missing-file warnings before using the assembled outputs.

| File | Contents |
|---|---|
| `<prefix>.db` | PredictDB `weights` and `extra` tables. The exporter retains summary rows with a finite `r_squared_global`. That legacy field stores Pearson's r. `R2_global` stores held-out 1 − RSS/TSS. |
| `<prefix>_cov.txt.gz` | Genotype covariances for the model SNPs, used alongside the database by S-PrediXcan. |
| `<prefix>_model_summaries.csv` | The combined per-gene performance summaries. |

See [the released model inventory](../../scTWAS_models/README.md) for the 23
cohort and cell-type databases.
