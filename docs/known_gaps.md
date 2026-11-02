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

## Also worth noting

- **GWAS meta-analysis.** The European-ancestry and multi-ancestry CRC meta-analyses were
  run with external meta-analysis software; only the script that assembles the Methods
  subsection and the combined cohort table is deposited, under `manuscript_text/`.
- **Cluster paths.** The `.sh` and `.sbatch` wrappers carry absolute ACCRE paths
  (`/v5000/data/l2_bioinfo1/liq17/...`). They record what was run; they will not resolve
  elsewhere without editing.
- **One hardcoded input path.** `09_functional_screens/48_screened_and_selected.py` points
  at an upload path from a different host and will not resolve as written.
- **A config filename typo.** `05_colocalization_smr_heidi/run_all_celltypes.sh` passes
  `coloc_peakSNP_ADSER.config.config.txt`; the file in `config/` is
  `coloc_peakSNP_ADSER.config.txt`. Left verbatim as it was run.

## Closed since the first deposit

The GReX model training code (`03_grex_model_training/`) and the association-testing code
(`04_association_testing/`) were listed as missing when the deposit was first assembled.
Both were recovered from the project repository and are now present, so the workflow has a
continuous path from raw data through to the downstream sections.
