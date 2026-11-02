# lib - shared helpers sourced by the analysis scripts

These four files are not analysis steps; they are the shared helper modules that several
scripts load at run time. They lived outside the project tree on the VUMC cluster, which is
why they were missing from the first deposit. They are included here so the scripts they
back can actually run.

| File | Loaded by | What it provides |
|---|---|---|
| `GenoPhenoBasic.R` | `06_tf_network/correlation_and_occupancy/1a_Extract_sigTWASGenes_SNPs_gtf.R`, `09_functional_screens/experimental_validation/1_DESeq2_DE_Genes_corrected_final.R` | Shared genotype/phenotype helpers. |
| `FunctionalBasic.R` | `07_orthogonal_support/FunctionalQuery_updated_v3.R` | Shared gene functional-annotation helpers. |
| `QTLS/annotation.py` | `01_single_cell_and_xqtl/normal_colon_WellcomeSanger/Prepare_pseudobulk_expression.py`, `01_single_cell_and_xqtl/polyps_COLONMAP/02_prepare_pseudobulk_expression_ASC_SSC_bashconfig.py` | Gene annotation used in pseudobulk preparation. |
| `QTLS/phenotype.py` | the same two scripts | Quantile normalisation and the inverse-normal transform described in Methods. |

## Repointing

The scripts load these by absolute cluster path, exactly as they were run:

```r
source("/data/l2_bioinfo1/liq17/LQ_frameworks/GenoPhenoBasic.R")
source("/data/l2_bioinfo1/liq17/LQ_references/GeneFunctioalSupports/FunctionalBasic.R")
```

```python
sys.path.insert(0, "/data/l2_bioinfo1/liq17/LQ_frameworks/QTLS")   # annotation, phenotype
```

Those lines are left verbatim rather than rewritten, so the deposited code still matches
what produced the results. Repoint them at this directory before running elsewhere.
