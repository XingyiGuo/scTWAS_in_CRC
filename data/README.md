# Data

Small, curated data files that the code in `code/` reads directly. Bulk files are released
through Zenodo instead of git -- see **Released through Zenodo** below.

## In this repository

### `scATAC-seq/used_peaks/`
The cell-type accessible-chromatin peak sets actually used to define candidate predictor
variants for GReX model training. These are the merged, filtered peak sets referred to in
Methods; they are the inputs to model training, not the raw MACS output.

- `normal/` -- 7 peak sets for normal colorectal mucosa cell types
  (ABS, CT, EE, GOB, STM, TAC from the VUMC + GSE201336 merge; TUF from the VUMC robust set).
- `disease/` -- 9 peak sets for polyp cell types
  (the VUMC adenoma/serrated + GSE201336 merge, plus `Adenoma_ASC` and `Serrated_SSC`).

BED format, one interval per line, genome build GRCh38.

### `scATAC-seq/peak_calling_summaries/`
Per-cell-type peak and cell counts written by the peak-calling and merging steps in
`code/02_accessible_chromatin/`. Tab-separated; these are the numbers behind the
peak-count panels.

### `GWAS_SS/`
See `GWAS_SS/README.md`. The CRC GWAS meta-analysis summary statistics generated in this
study are deposited in the GWAS Catalog; that directory records the accessions rather than
duplicating the files, which are far too large for git.

## Released through Zenodo

Not in this repository, to keep clones small. Download from the Zenodo record cited in the
top-level `README.md`:

| Item | Files | Size |
|---|---|---|
| `scTWAS_models/` -- the 23 trained GReX models: each `.db` plus its `*_cov.txt.gz` covariance matrix and `*_model_summaries.csv` | 69 | 208 MB |
| `CRC_GWAS_SS/` -- CRC GWAS meta-analysis summary statistics, per chromosome, EUR and multi-ancestry | 45 | 812 MB |
| GSE201336 per-group narrowPeak calls | 12 | 110 MB |
| VUMC robust peak regions (BED) | 24 | 41 MB |
| VUMC adenoma/serrated chaining-merged peaks (BED) | 7 | 14 MB |

The GWAS summary statistics are also deposited in the NHGRI-EBI GWAS Catalog, which is the
citable source of record; see `GWAS_SS/README.md`.

The `.counts.tsv` and summary tables for those peak sets *are* kept here, under
`scATAC-seq/peak_calling_summaries/`, so the reported counts can be checked without the
download.

## Data not redistributed here

Single-cell, genotype and spatial source data are controlled by their originating studies
and are obtained from the repositories listed in the manuscript's Data availability
statement (HTAN Data Portal, GEO, EGA, and the spatial Zenodo record). Nothing in this
directory substitutes for those.
