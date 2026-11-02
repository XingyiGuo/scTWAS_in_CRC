# Analysis workflow, in execution order

Each step names the directory that holds its code. Steps 3 and 4 are the two components
that are not in this repository (see `known_gaps.md`); they are listed in place so the
order of the published workflow stays readable.

## 1. Single-cell references and cis-eQTL mapping — `01_single_cell_and_xqtl/`

Three training cohorts, each in its own directory and each following the same arc.

- `normal_colon_GSE266546/` — VUMC Gut Cell Atlas normal colon. Cell calling and ambient
  correction, integration, clustering and annotation, sub-clustering, marker ranking and
  sample-level pseudobulk counts (`01`-`07`); expression-genotype matching with segment
  priority (`08`); PEER covariate residualisation; genotype merge and conversion
  (`g1`, `g2`); and tensorQTL input preparation and per-chromosome runs
  (`xQTL_tensorqtls/`).
- `normal_colon_WellcomeSanger/` — Wellcome Sanger normal colon. Conversion of the
  published h5ad to a BPCells-backed Seurat object, extraction of the healthy terminal
  ileum and rectum donors, epithelial refresh including transit-amplifying cells, UMAP
  relabelling, then pseudobulk preparation, PEER and genotype conversion.
- `polyps_COLONMAP/` — VUMC HTAN COLON MAP polyps. Cell-type assembly over the COLON MAP
  workflow, expression-genotype matching, pseudobulk preparation, PEER residualisation
  including specimen class, and tensorQTL preparation and runs.
- `polyps_COLONMAP_annotation/` — the adenoma/serrated annotation and UMAP work for the
  polyp models, including the merged-QC PCA and per-cell-type metadata assembly.
- `cohort_metadata/` — donor lists, cells per sample by cell type, specimen selection for
  the Gut Cell Atlas pairs, and the expression-identifier provenance table.

Pseudobulk normalisation follows the Methods: raw UMI counts aggregated per donor and cell
type, converted to transcripts per million and log2-transformed, quantile-normalised across
donors and inverse-normal transformed by gene, then residualised for age, sex, five
genotype principal components, PEER factors and (for polyps) specimen class. The
quantile-normalisation, inverse-normal-transform and gene-annotation helpers are imported
from a shared `QTLS` framework that is **not** in this repository — see `known_gaps.md`.

## 2. Accessible-chromatin peak sets — `02_accessible_chromatin/`

ArchR arrow-file construction, label transfer and MACS peak calling per cell type for the
VUMC scATAC-seq cohort (`vumc_scATAC/`, MACS3) and the Stanford HTAN GSE201336 cohort
(`GSE201336_scATAC/`, ArchR plus MACS2), then the summit-based merge across cohorts
(`merged_peaks/`). The merged cell-type peak sets define the candidate predictor variants
for model training.

## 3. GReX model training — not in this repository

Elastic-net (alpha = 0.5) models per gene and cell type over variants inside the cell-type
accessible chromatin and within 1 Mb of the transcription start site, trained under the
sTF-TWAS framework following PrediXcan, with performance assessed by three repeated
tenfold cross-validation runs. See `03_grex_model_training/README.md`.

## 4. Association testing — not in this repository

S-PrediXcan-style summary-statistic association of each model against the
European-ancestry and multi-ancestry CRC GWAS meta-analyses, with the two normal-colon
training cohorts combined by the Cauchy combination test (ACAT) and Bonferroni correction
applied within each cell type. See `04_association_testing/README.md`.

## 5. Colocalization and SMR-HEIDI — `05_colocalization_smr_heidi/`

`01_coloc_gwas_xqtl_sigSource_peakSNP.R` runs coloc with the approximate Bayes factor
approach for the significant gene-cell-type pairs against the source-matched cis-eQTL
data, driven per cell type by a config file; `01a_*_summarize.R` collapse the per-cell-type
output for the normal-colon (7 cell types) and polyp (9 populations) model sets.
`02_SMR_HEIDI_inputs_targetTopEQTL_500kb_directTWAS.R` builds the SMR inputs, selecting the
top cis-eQTL per gene as the instrument and writing the locus-specific PLINK reference;
`02a`/`02b` prepare and run SMR v1.4.2 with HEIDI.

## 6. TF-target networks — `06_tf_network/`

`correlation_and_occupancy/` extracts the model SNPs for significant TWAS genes from the
per-cell-type model databases and performs the one-pass GWAS rsID lookup (`1a`), computes
TF-target pseudobulk correlations for the polyp (9 cell types) and normal (7 cell types)
model sets (`2a`), and tests lead-variant TF occupancy against target-gene correlation
(`2b`). `network_figures/` builds the matrices, orderings and network panels and
reconstructs per-panel source data. `pathway_integration/` joins the network with the
pathway enrichment results and tests cell-type specificity.

## 7. Orthogonal molecular support — `07_orthogonal_support/`

`multiomics_de_support/` holds the numbered chain that assembles differential-expression
support across sources: model overlap, per-contrast single-cell support, gene-level
direction, cross-cell-type support, up/down summaries, the linkage diagram, the CPTAC
limma re-analysis and proteome-wide differential expression, the CosMx highlight panels
and the CosMx supplementary table. `pseudobulk_de/` carries the pseudobulk aggregation and
the differential-expression contrast drivers. `spatial_geomx_cosmx/` holds the GeoMx
extraction and spatial differential expression, and the CosMx per-cell, composition,
gene-level paired-test and panel code.

## 8. Druggability and connectivity — `08_druggability_and_connectivity/`

`gene_sets_and_resources/` harmonises the prioritized gene list from the association
results and builds the public resource tables. `druggability/` integrates DrugBank,
ChEMBL, TTD and Open Targets into the druggable-gene and compound tables, the KEGG
enrichment for druggable genes, the clinical-phase and indication annotation, and the
figure export and source data. `ilincs_connectivity/` scores the prespecified compound
panel against the iLINCS LIB5 large-intestine library, including the size-matched
empirical null built from 5,000 randomly sampled signature sets, the consensus check and
the panel figure.

## 9. Functional screens — `09_functional_screens/`

Parsing and quality control of the experimenters' screen workbook, per-cell-line siRNA and
CRISPR-Cas9 statistics with paired t-tests across replicate gene-control pairs, replicate
quality control, concordance, risk-direction and worked-example figures, and the
screened-versus-selected accounting.

## 10. Pathway enrichment — `10_pathway_enrichment/`

Enrichment by cell type and by risk direction, the merged analyses and their figures.

## 11. Figures — `11_figures/`

`manhattan/` renders the stacked per-cell-type Manhattan panels for the normal-colon and
polyp model sets, in European-ancestry and multi-ancestry-only variants, with the panel A
source-data extractor. `workflow/` draws the two analytic workflow schematics.
`figure7_depmap_screens/` is the Figure 7 panel layout; note that it is a rendering
function only, and its analysis caller is missing (see `known_gaps.md`).

## 12. Supplementary tables — `12_supplementary_tables/`

Assembly of the supplementary workbooks and the consolidated results record spanning the
orthogonal, druggability and screening layers.
