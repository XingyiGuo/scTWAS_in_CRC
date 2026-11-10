# lib - shared helpers sourced by the analysis scripts

These four files are not analysis steps; they are the shared helper modules that several
scripts load at run time. They lived outside the project tree on the VUMC cluster, which is
why they were missing from the first deposit. They are included here so the scripts they
back can actually run.

| File | Loaded by | What it provides |
|---|---|---|
| `GenoPhenoBasic.R` | `03_grex_model_training/TWAS_modelbuilding_scTWAS.R` | Shared genotype/phenotype helpers. |
| `FunctionalBasic.R` | Standalone helper | Shared gene functional-annotation helpers. |
| `QTLS/annotation.py` | `01_single_cell_and_xqtl/normal_colon_WellcomeSanger/Prepare_pseudobulk_expression.py`, `01_single_cell_and_xqtl/polyps_COLONMAP/02_prepare_pseudobulk_expression_ASC_SSC_bashconfig.py` | Gene annotation used in pseudobulk preparation. |
| `QTLS/phenotype.py` | the same two scripts | Quantile normalisation and the inverse-normal transform described in Methods. |

## Repointing

The model trainer uses this directory by default. The other R scripts load
their helpers by absolute cluster path. Repoint those paths when running elsewhere.

```r
source("/data/l2_bioinfo1/liq17/LQ_frameworks/GenoPhenoBasic.R")
source("/data/l2_bioinfo1/liq17/LQ_references/GeneFunctioalSupports/FunctionalBasic.R")
```

The Python scripts try the original cluster paths and fall back to the bundled
`QTLS/` directory. The Wellcome Sanger script also accepts `--qtls-dir`.
