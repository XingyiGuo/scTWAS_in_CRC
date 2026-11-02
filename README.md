# scTWAS in colorectal cancer

Code and curated data for **"Single-cell transcriptome-wide association study in colorectal
cancer"** (Li, Li, Lyu, Wang, Yuan, Song, Yan, Choi, Xu, Wang, Sheng, Shu, Singh, Coburn,
Wilson, Ciombor, Eng, Peters, Yin, Zheng, Lau & Guo).

The study builds cell-type-specific models of genetically regulated expression (GReX) in
normal colorectal mucosa and in colorectal polyps, integrates them with colorectal cancer
GWAS, and follows the prioritized genes through regulatory, spatial, proteomic,
druggability and functional-screen analyses.

## Layout

```
<<<<<<< HEAD
code/     234 analysis files. 12 Methods-ordered sections (01_single_cell_and_xqtl
          through 12_supplementary_tables) plus manuscript_text; Viz/ for the figure
          scripts; and lib/ for the shared helpers several scripts load at run time
=======
code/     230 analysis files in 13 sections, ordered as in Methods
>>>>>>> 2de678c2bfef89f7c1f9743c89a392869bce8c37
data/     curated inputs that the code reads directly (bulk data via Zenodo)
docs/     data flow, workflow order, figure/table map, known gaps
env/      Python, R and external-tool environment specifications
scTWAS_models/   trained GReX models (released via Zenodo)
MANIFEST.csv     every file with its type, size and MD5
```

## Analysis sections

Each section directory carries its own `README.md` listing its scripts in run order, their
resolved inputs, and the downstream scripts that consume each output. Start there.

| Section | Files | What it does |
|---|---|---|
| [`01_single_cell_and_xqtl/`](code/01_single_cell_and_xqtl) | 51 | Single-cell reference construction and cis-xQTL mapping for the three training cohorts: cell calling and annotation, pseudobulk aggregation, PEER residualisation, genotype preparation and tensorQTL. |
| [`02_accessible_chromatin/`](code/02_accessible_chromatin) | 19 | Cell-type accessible-chromatin peak calling for the VUMC and GSE201336 scATAC-seq cohorts, and the merge that defines the candidate predictor variants. |
| [`03_grex_model_training/`](code/03_grex_model_training) | 6 | Elastic-net training of the cell-type GReX models restricted to peak variants, and assembly of the weights into the PredictDB `.db` files. |
| [`04_association_testing/`](code/04_association_testing) | 5 | S-PrediXcan association testing against the CRC GWAS meta-analyses, the Cauchy combination across cell types, and selection of representative genes. |
| [`05_colocalization_smr_heidi/`](code/05_colocalization_smr_heidi) | 14 | Bayesian colocalization, SMR with the HEIDI heterogeneity test, and the summaries that join TWAS, coloc and SMR into the reported gene sets. |
| [`06_tf_network/`](code/06_tf_network) | 25 | Cell-type TF-target networks anchored on the 51 CRC susceptibility-associated transcription factors: model-SNP extraction, pseudobulk correlation, lead-variant TF occupancy, motif logos and network figures. |
| [`07_orthogonal_support/`](code/07_orthogonal_support) | 43 | Orthogonal molecular support: single-cell and pseudobulk differential expression, CPTAC and TCGA re-analysis, and the GeoMx and CosMx spatial analyses. |
| [`08_druggability_and_connectivity/`](code/08_druggability_and_connectivity) | 19 | Drug-target annotation from DrugBank, ChEMBL, TTD and Open Targets, and the iLINCS LIB5 connectivity analysis with its size-matched empirical null. |
| [`09_functional_screens/`](code/09_functional_screens) | 23 | Arrayed siRNA and CRISPR-Cas9 viability screens, and the DESeq2 differential-expression analysis of the knockout validation experiment. |
| [`10_pathway_enrichment/`](code/10_pathway_enrichment) | 7 | Pathway and gene-set enrichment across cell types and risk directions. |
| [`11_figures/`](code/11_figures) | 14 | Figure assembly: stacked Manhattan panels, workflow schematics, model R-squared violins and peak counts, annotation UMAPs. |
| [`12_supplementary_tables/`](code/12_supplementary_tables) | 3 | Assembly of the supplementary tables and the consolidated results record. |
| [`manuscript_text/`](code/manuscript_text) | 1 | Scripts that only recompute reported numbers and assemble manuscript prose. No analysis is performed here. |

## Data

`data/` holds only what is small and directly read by the code -- the cell-type
accessible-chromatin peak sets used for model training, and the peak-count summaries.
See [`data/README.md`](data/README.md).

The trained GReX models, the raw per-group peak calls and the GWAS summary statistics are
**not** in git: models and raw peaks are released through Zenodo (see the DOI below), and
the GWAS meta-analyses are deposited in the GWAS Catalog
(see [`data/GWAS_SS/README.md`](data/GWAS_SS/README.md)). Single-cell, genotype and spatial
source data are controlled by their originating studies and are obtained from the
repositories named in the manuscript's Data availability statement.

## Reproducing

1. Install the environments in [`env/`](env) (`python_requirements.txt`, `r_packages.txt`,
   and the external tools in `external_tools.md` -- PLINK, tensorQTL, S-PrediXcan, SMR,
   MACS3, ArchR).
2. R dependencies are pinned in `env/renv.lock` (`renv::restore()`). Note it was built on
   Linux (R 4.5); on another platform renv resolves platform-appropriate builds.
3. Obtain the controlled source data from the accessions in the Data availability
   statement, and download the Zenodo archive for the trained models.
4. Run sections in numeric order. [`docs/workflow.md`](docs/workflow.md) gives the
   execution order and [`docs/data_flow.md`](docs/data_flow.md) the file-level links
   between sections.

Scripts keep the paths they ran under on the VUMC cluster; those are recorded in
`MANIFEST.csv` and in each section README, and need to be repointed to a local layout.
[`docs/known_gaps.md`](docs/known_gaps.md) records what is not reproducible from this
repository alone.

## Citation

See [`CITATION.cff`](CITATION.cff). Archived on Zenodo at DOI *(to be added on release)*.

Released under the terms in [`LICENSE`](LICENSE).

## Contact

Qing Li -- liqingbioinfo@gmail.com \
Xingyi Guo -- xingyi.guo@vumc.org
