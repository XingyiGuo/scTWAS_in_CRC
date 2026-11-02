# 02 — Cell-type accessible-chromatin peak sets

ArchR arrow-file construction, label transfer and MACS peak calling per cell type for the
VUMC scATAC-seq cohort (`vumc_scATAC/`, MACS3) and the Stanford HTAN GSE201336 cohort
(`GSE201336_scATAC/`, ArchR plus MACS2), then the summit-based merge across cohorts
(`merged_peaks/`). The merged cell-type peak sets define the candidate predictor variants
for model training.

## Position in the pipeline

```
(external data — see Inputs)
        |
        v
   >>> 02_accessible_chromatin <<<
        |
        v
03_grex_model_training  11_figures
```

The arrows above are the logical order of the published workflow, taken from the Methods (`docs/workflow.md`). The file-level links in *From other sections* and *Consumed by* below are the ones actually observed in the deposited code.

## Scripts, in run order (19 files)

Ordered by sub-directory, then by the numeric prefix the script carries in its original working directory. Descriptions are the script's own header line; blank where the script has no header.

| # | File | Lang | What it does | Declared display item |
|---|---|---|---|---|
| 1 | `GSE201336_scATAC/01_build_arrows.R` | R | 01 -- Build ArchR Arrow files from GSE201336 fragment files. Arrow names are set to the authors' ATAC_Sample IDs (e.g. |  |
| 2 | `GSE201336_scATAC/02_relabel.R` | R | 02 -- Re-label Becker et al. scATAC epithelial cells into the Chen/Lau (Cell 2021) nomenclature and attach sample metadata. |  |
| 3 | `GSE201336_scATAC/02_relabel_and_transfer.R` | R | 02 -- Re-label Becker et al. scATAC epithelial cells into the Chen/Lau (Cell 2021) nomenclature, and recover the adenoma-specific (ASC) population… |  |
| 4 | `GSE201336_scATAC/03_call_peaks.R` | R | 03 -- Reproducible peak set for every CELL TYPE, pooling cells across tissue conditions (normal + adenoma/polyp). |  |
| 5 | `GSE201336_scATAC/03_call_peaks.sbatch` | Slurm |  |  |
| 6 | `GSE201336_scATAC/04_export_peaks.R` | R | 04 -- Export the reproducible peaks to BED, the union peak set, and a pseudobulk accessibility matrix (broken out by Condition x CellType). |  |
| 7 | `GSE201336_scATAC/04_export_peaks.sbatch` | Slurm |  |  |
| 8 | `GSE201336_scATAC/05_GSE201336_WellcomeSanger_peak_calling_summary.sh` | Shell |  |  |
| 9 | `GSE201336_scATAC/build_manifest.py` | Python | Build GSE201336_atac_manifest.csv joining GEO GSM accessions to Becker et al. sample metadata. |  |
| 10 | `GSE201336_scATAC/check_annotation_join.R` | R | PRE-FLIGHT CHECK -- run this before 01/02. Costs seconds; saves hours. |  |
| 11 | `GSE201336_scATAC/macs2_wrapper.sh` | Shell |  |  |
| 12 | `GSE201336_scATAC/run_pipeline.sbatch` | Slurm |  |  |
| 13 | `merged_peaks/01_VUMC_GSE201336_merge_with_peak_counts.sh` | Shell |  |  |
| 14 | `vumc_scATAC/01_VUMC_ATAC_MACS3_callpeak_parallel8.sh` | Shell |  |  |
| 15 | `vumc_scATAC/01b_VUMC_ATAC_prepare_Normal_ready_peaks.sh` | Shell |  |  |
| 16 | `vumc_scATAC/01c_VUMC_add_cellcounts_to_peak_calling_summary.sh` | Shell |  |  |
| 17 | `vumc_scATAC/02_VUMC_ATAC_MACS3_calls_chaining_merge.sh` | Shell |  |  |
| 18 | `vumc_scATAC/02_VUMC_ATAC_MACS3_calls_summitmerge_with_summit500.sh` | Shell |  |  |
| 19 | `vumc_scATAC/03_VUMC_GSE201336_merge.sh` | Shell |  |  |

## Inputs

**Not produced anywhere in this deposit** (2) — cohort data, public resources, or output of a step that is not deposited:

- `HuBMAPCellAnnotations_hubmap_colon_epithelial.tsv`
- `all_samples_cell_types_projections.csv`

## Outputs

| File | Written by | Consumed by |
|---|---|---|
| `adenoma_label_transfer_Lau.csv` | `GSE201336_scATAC/02_relabel_and_transfer.R` | — (terminal) |
| `celltype_counts_lau.csv` | `GSE201336_scATAC/02_relabel.R` | — (terminal) |
| `peak_by_condition_celltype_cpm.csv.gz` | `GSE201336_scATAC/04_export_peaks.R` | — (terminal) |
| `peak_calling_group_sizes.csv` | `GSE201336_scATAC/03_call_peaks.R` | — (terminal) |
| `peak_summary_by_group.csv` | `GSE201336_scATAC/04_export_peaks.R` | — (terminal) |
| `union_peakset_epithelial.rds` | `GSE201336_scATAC/04_export_peaks.R` | — (terminal) |

## How this file was built

Automatically generated 2026-10-19 from MANIFEST.csv, the Methods order in `docs/workflow.md`, and the file reads and writes parsed out of the deposited scripts themselves. The extraction reads only string literals and path variables that appear at a read or write call site, so three classes of file are invisible to it: paths assembled from shell variables (most `.sh` and `.sbatch` wrappers), paths built by string formatting at run time, and files handled by external tools that the scripts invoke rather than open themselves (tensorQTL, PLINK, MACS, SMR, PEER). The lists above are therefore a lower bound, not a complete inventory. `MANIFEST.csv` at the repository root is authoritative for the file list itself.

## Added in this release

These files were merged in from the project repository after the first deposit was assembled, so they are not described in the run order above.

`GSE201336_scATAC/01_build_arrows.R` is the newer copy that carries the GEO download URL for GSE201336.

- `GSE201336_scATAC/01_build_arrows.R`
