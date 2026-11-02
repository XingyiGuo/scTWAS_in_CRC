# Manuscript display items and the code that produces them

Two kinds of evidence are distinguished below. *Declared* means the script's own header
names the display item; those entries also appear in the `declared_display_items_in_header`
column of `MANIFEST.csv`. *Methods* means the correspondence follows from the Methods
section rather than from a statement in the code.

| Display item | Section | Evidence |
|---|---|---|
| Figure 1 (study design, evidence layers) | `11_figures/workflow/` | Methods |
| Figure 2a–g (normal colon, European ancestry) | `11_figures/manhattan/Manhattan_NL_EUR.R` | Methods |
| Figure 3a–f (polyps, European ancestry) | `11_figures/manhattan/Manhattan_ADSER_EUR.R` | Methods |
| Figure 4 (TF–target networks) | `06_tf_network/` | Declared |
| Figure 4e,f (TF–pathway integration) | `06_tf_network/pathway_integration/tf_pathway_integration.py` | Declared |
| Figure 5c (GeoMx, TCGA, CPTAC support) | `07_orthogonal_support/multiomics_de_support/`, `07_orthogonal_support/spatial_geomx_cosmx/` | Declared |
| Figure 5d (CosMx, six per-cell genes) | `07_orthogonal_support/multiomics_de_support/08_prep_cosmx_highlight_panels.py`, `09_fig_NOR_CORE_highlight_4panel.R` | Declared |
| Figure 6 (druggability) | `08_druggability_and_connectivity/druggability/` | Declared |
| Figure 6e (iLINCS connectivity) | `08_druggability_and_connectivity/ilincs_connectivity/10_fig6e_panel.py` | Declared |
| Figure 7a,b (DepMap dependency) | `11_figures/figure7_depmap_screens/build_figure7.py` — rendering only, analysis code missing | Methods; see `known_gaps.md` |
| Figure 7c (arrayed CRISPR–Cas9) | `09_functional_screens/`, `11_figures/figure7_depmap_screens/` | Methods |
| Extended Data Fig. 3 (normal colon, multi-ancestry) | `11_figures/manhattan/Manhattan_NL_Multi_EURnotSig_MultiSig_only.R` | Methods |
| Extended Data Fig. 4 (polyps, multi-ancestry) | `11_figures/manhattan/Manhattan_ADSER_Multi_EURnotSig_MultiSig_only.R` | Methods |
| Extended Data Fig. 7 (iLINCS connectivity detail) | `08_druggability_and_connectivity/ilincs_connectivity/12_edfig7_caption.py` | Declared |
| Supplementary Figure 2 (sTF-TWAS workflow style reference) | `11_figures/workflow/Fig_scTWAS_workflow_NL_ADSER.py` | Declared |
| Accessible-chromatin peak counts and model R-squared violins (normal colon and polyps) | `11_figures/model_performance/Fig_peaks_R2_violin.R` | Methods; reads Supplementary Table 5 peak counts and the per-cell-type model summaries |
| Cell annotation UMAP panels (VUMC normal colon, Wellcome Sanger normal colon) | `11_figures/umap_annotation/` | Methods |
| Supplementary Table 2 (donors and cells) | `01_single_cell_and_xqtl/cohort_metadata/` | Declared |
| Supplementary Table 5 (scATAC-seq peaks) | `02_accessible_chromatin/` | Methods |
| Supplementary Table 4 (51 anchor TFs) | `06_tf_network/correlation_and_occupancy/` | Methods |
| Supplementary Table 6 (GWAS cohorts) | `manuscript_text/build_GWAS_meta_NG.py` | Methods |
| Supplementary Tables 7–10 (scTWAS associations, coloc, SMR–HEIDI) | `05_colocalization_smr_heidi/`; association step missing | Methods; see `known_gaps.md` |
| Supplementary Table 18 (orthogonal support) | `07_orthogonal_support/multiomics_de_support/` | Methods |
| Supplementary Table 20 (CosMx single-cell support) | `07_orthogonal_support/multiomics_de_support/14_build_supptable_S20_cosmx.py` | Declared |
| Supplementary Table 22 (iLINCS connectivity) | `08_druggability_and_connectivity/ilincs_connectivity/` | Declared |
| Supplementary Table 23 (DepMap dependency) | analysis code missing | see `known_gaps.md` |
| Supplementary Table 24 (arrayed screens) | `09_functional_screens/`, `12_supplementary_tables/` | Methods |

## A caution on figure numbering inside the code

Several scripts in `07_orthogonal_support/spatial_geomx_cosmx/`, `09_functional_screens/`,
`10_pathway_enrichment/` and `12_supplementary_tables/` were written as part of an internal
analysis report with its own figure and table numbering, which runs well past the
manuscript's. Where a header in those directories says "Figure 22" or "S36" it refers to
that internal report, not to a manuscript display item; `MANIFEST.csv` flags those entries
as `[internal report numbering]`. The same applies to the TF-network figure code in
`06_tf_network/network_figures/`, whose "Figures 1–3" are the panels of the TF analysis
report that became manuscript Figure 4.
