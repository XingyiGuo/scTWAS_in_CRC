# Trained GReX models

The 23 cell-type-specific models of genetically regulated expression (`.db`, PredictDB
format, 75 MB total) are **released through Zenodo**, not git -- see the DOI in the
top-level `README.md`.

Three model sets are provided, one per training cohort:

| Prefix | Cohort | Tissue | Cell types |
|---|---|---|---|
| `VUMC_HTAN_ADSER_<CT>_peakSNP.db` | VUMC HTAN COLON MAP | colorectal polyps (adenoma / serrated) | 9 -- ABS, ASC, CT, EE, GOB, SSC, STM, TAC, TUF |
| `VUMC_HTAN_GSE266546_NL_<CT>_peakSNP.db` | VUMC HTAN + GSE266546 | normal colorectal mucosa | 7 -- ABS, CT, EE, GOB, STM, TAC, TUF |
| `WellcomeSanger_NL_<CT>_peakSNP.db` | Wellcome Sanger Gut Cell Atlas | normal colorectal mucosa | 7 -- ABS, CT, EE, GOB, STM, TAC, TUF |

Each database holds the elastic-net predictor weights and the per-gene cross-validation
performance written by `code/03_grex_model_training/WeightsSummaryToDB.R`. Predictor
variants are restricted to the cell-type accessible-chromatin peaks in
`data/scATAC-seq/used_peaks/`.

To use them with S-PrediXcan, download the Zenodo archive and unpack it into this
directory; the scripts in `code/04_association_testing/` resolve model paths relative to
here.
