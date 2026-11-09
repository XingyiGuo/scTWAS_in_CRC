# Known gaps - what cannot be re-run from this repository alone

Recorded here rather than silently omitted.

## 1. DepMap dependency analysis (Figure 7a,b; Supplementary Table 23)

**Status: partially present.** `11_figures/figure7_depmap_screens/build_figure7.py` is a
pure rendering function: it receives the gene-effect distributions, the pairwise test table
and the per-gene values as arguments and draws the panels. The code that retrieved Chronos
gene-effect scores from the Open Targets GraphQL API (release 26.09,
`depMapEssentiality`), restricted them to the 63 colorectal adenocarcinoma lines in the
DepMap intestine category, assembled the four comparison gene sets and ran the
Kruskal-Wallis and Holm-corrected Mann-Whitney tests was not found.

The retrieved data and the test results survive as exported tables next to the figure
script in the working directory (`DepMap_essentiality_by_geneset_sourcedata.csv`,
`depmap_intestine_*_percellline.csv`,
`SupplementaryTable_DepMap_essentiality_by_geneset.xlsx`), so the figure can be redrawn,
but the retrieval and statistics cannot be re-run from code.

<<<<<<< HEAD
This is the only remaining gap in the analysis chain.
=======
## 2. The shared `QTLS` helper framework

**Status: not present.** Two pseudobulk preparation scripts --
`01_single_cell_and_xqtl/normal_colon_WellcomeSanger/Prepare_pseudobulk_expression.py` and
`01_single_cell_and_xqtl/polyps_COLONMAP/02_prepare_pseudobulk_expression_ASC_SSC_bashconfig.py`
-- import `annotation.py` and `phenotype.py` from a shared framework directory, resolved at
run time from the cluster path `/data/l2_bioinfo1/liq17/LQ_frameworks/QTLS` (with a
`/v5000/...` alternative). These modules implement the gene annotation, quantile
normalisation and inverse-normal transform that the Methods describe, so the normalisation
is defined in code that is not deposited. The two modules are small and should be added to
`01_single_cell_and_xqtl/` before release.
>>>>>>> 2de678c2bfef89f7c1f9743c89a392869bce8c37

## Also worth noting

- **GWAS meta-analysis.** The European-ancestry and multi-ancestry CRC meta-analyses were
  run with external meta-analysis software; only the script that assembles the Methods
  subsection and the combined cohort table is deposited, under `manuscript_text/`.
<<<<<<< HEAD
- **Cluster paths.** The `.sh` and `.sbatch` wrappers, and a handful of R and Python
  scripts, carry absolute ACCRE paths (`/data/l2_bioinfo1/liq17/...`,
  `/v5000/data/l2_bioinfo1/liq17/...`). They record what was run; they will not resolve
  elsewhere without editing. This includes the three `source()` / `sys.path` lines that
  load the shared helpers now shipped in `code/lib/` - see `code/lib/README.md`.
=======
- **Cluster paths.** The `.sh` and `.sbatch` wrappers carry absolute ACCRE paths
  (`/v5000/data/l2_bioinfo1/liq17/...`). They record what was run; they will not resolve
  elsewhere without editing.
>>>>>>> 2de678c2bfef89f7c1f9743c89a392869bce8c37
- **One hardcoded input path.** `09_functional_screens/48_screened_and_selected.py` points
  at an upload path from a different host and will not resolve as written.
- **A config filename typo.** `05_colocalization_smr_heidi/run_all_celltypes.sh` passes
  `coloc_peakSNP_ADSER.config.config.txt`; the file in `config/` is
  `coloc_peakSNP_ADSER.config.txt`. Left verbatim as it was run.
<<<<<<< HEAD
- **Stale filenames in usage text (25 places).** Many scripts print or comment a usage
  line naming an earlier version of themselves, e.g. `03_grex_model_training/WeightsSummaryToDB.R`
  prints `Usage: Rscript WeightsSummaryToDB_peakSNP.R ...`, and several scripts in
  `01_single_cell_and_xqtl/` carry header names from before the deposit was renumbered
  (`03_cluster_annot.R` says `03_cluster_plots.R`, `01_match_HTAN_ASC_SSC_...py` says
  `08_match_...py`). These are help strings and comments only: as of this audit no script
  executes a file that is absent. Left verbatim so the code matches what was run.
- **Shell syntax unverified.** The Python and R files were all machine-checked for
  parseability (106 Python, 76 R, 0 failures). The 44 shell and Slurm scripts were not:
  no shell was available to run `bash -n` against them.

## Closed since the first deposit

- **GReX model training** (`03_grex_model_training/`) and **association testing**
  (`04_association_testing/`) were listed as missing. The elastic-net trainer, the
  PredictDB builder, the S-PrediXcan runner and the ACAT Cauchy combination were recovered
  from the project repository and are now present.
- **The shared `QTLS` helper framework** was listed as missing. `annotation.py` and
  `phenotype.py` - which define the gene annotation, quantile normalisation and
  inverse-normal transform the Methods describe - were recovered and are now in
  `code/lib/QTLS/`, together with the two R helper files (`GenoPhenoBasic.R`,
  `FunctionalBasic.R`) that three scripts `source()` and that no earlier audit had noticed.
=======

## Closed since the first deposit

The GReX model training code (`03_grex_model_training/`) and the association-testing code
(`04_association_testing/`) were listed as missing when the deposit was first assembled.
Both were recovered from the project repository and are now present, so the workflow has a
continuous path from raw data through to the downstream sections.
>>>>>>> 2de678c2bfef89f7c1f9743c89a392869bce8c37
