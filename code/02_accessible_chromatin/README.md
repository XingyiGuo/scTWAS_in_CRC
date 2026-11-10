# Cell-type accessible chromatin

This section calls epithelial scATAC-seq peaks in VUMC and GSE201336, then
combines the cohorts by cell type and disease group. The released BEDs define
the accessible regions used to select candidate variants for
[GReX model training](../03_grex_model_training).

Use the [released peak BEDs](../../data/scATAC-seq/used_peaks) if you are starting
at model training. The [peak-calling summaries](../../data/scATAC-seq/peak_calling_summaries)
record the VUMC peak counts and overlap statistics.

## VUMC scATAC-seq

Run these scripts from [vumc_scATAC](vumc_scATAC). Set the project paths in each
script and provide hg38 BAMs split by condition and cell type, plus the hg38
blacklist BED. This branch uses SAMtools, MACS3 and BEDTools.

1. `01_VUMC_ATAC_MACS3_callpeak_parallel8.sh` filters the split BAMs, calls
   peaks at MACS3 q < 0.01 and removes blacklisted intervals. For groups with
   fewer than one million usable fragments, it also calls peaks in three 80%
   subsamples and retains intervals supported by at least two runs. It writes
   the per-group calls under `peaks/` and a `peak_calling_summary.tsv` that
   includes cell counts from BAM barcode tags.
2. `02_VUMC_ATAC_MACS3_calls_chaining_merge.sh` selects robust BEDs for the
   deposited cohort. It uses subsampling-supported peaks for Normal TAC,
   Adenoma ASC and TUF, and Serrated SSC and TUF, with blacklist-clean
   full-depth calls for the other groups. It writes `robust_peak_regions/`
   and combines adenoma and serrated intervals for the seven shared cell
   types under `Adenoma_Serrated_chaining_merged_peaks/`. ASC and SSC retain
   their separate robust BEDs.

## GSE201336 scATAC-seq

Use the ArchR environment in [the GSE201336 R lockfile](../../env/renv_GSE201336_scATAC.lock),
the hg38 genome package and MACS2. Stage the GEO fragment files in `WORK/fragments/`
with names beginning `<GSM>_` and ending `_fragments.tsv.gz`. Also stage these
annotation inputs in `WORK`:

- `epithelial_celltypes_atac.tsv` with `Cell` and `CellType` columns. Cell IDs
  use the `<ATAC_Sample>#<barcode>` format.
- `lau_label_map.csv` with `BeckerATAC` and `LauLabel` columns.

Run `build_manifest.py` from `WORK` to create
`GSE201336_atac_manifest.csv`. Use `check_annotation_join.R` to check the
annotation sample names against that manifest before building Arrow files.

| Order | Script | Handoff |
|---|---|---|
| 1 | `01_build_arrows.R` | Matches fragments by GSM accession, names the Arrow files with the annotation sample IDs, and saves `archr/proj_all_cells`. |
| 2 | `02_relabel.R` | Joins the epithelial annotations, applies the Lau label map, and saves `archr/proj_epithelial_lau`. The six retained lineages used for the merge are ABS, CT, EE, GOB, STM and TAC. |
| 3 | `03_call_peaks.R` | Calls reproducible peaks from ArchR pseudobulk replicates with MACS2. Pass `--byCondition --conditions Normal,Adenoma` to produce the condition-specific groups used in the merge. The script saves `archr/proj_epithelial_peaks` and `peak_calling_group_sizes.csv`. |
| 4 | `04_export_peaks.R` | Exports `peaks/per_group/{Condition}__{CellType}.narrowPeak.bed`, the union peak set and a pseudobulk accessibility matrix. |
| 5 | `05_GSE201336_WellcomeSanger_peak_calling_summary.sh` | Reports cell counts, blacklist-clean peak counts and interval widths. This reporting step writes separate filtered BEDs under `blacklist_clean/`. |

To run steps 1–4 together, submit from [GSE201336_scATAC](GSE201336_scATAC):

```bash
mkdir -p logs
WORK=/path/to/staged/work sbatch run_pipeline.sbatch
```

Review the Slurm account and module settings first. Set `MACS2` to the
absolute executable path if ArchR needs help locating it. The included
`macs2_wrapper.sh` loads the original cluster environment and needs executable
permissions if used as that path. The full wrapper uses a 200-cell cutoff per
group. The standalone `03_call_peaks.sbatch` is configured for Normal only
with a 100-cell cutoff, and `04_export_peaks.sbatch` exports the saved project.
Check the group-size report and make sure all six cell types have BEDs for both
conditions before merging.

## Cross-cohort merge

Run [the cross-cohort merge script](merged_peaks/01_VUMC_GSE201336_merge_with_peak_counts.sh)
after setting its VUMC, GSE201336 and output directories. It pairs the VUMC
robust Normal BEDs with the GSE201336 Normal exports, and the VUMC
adenoma–serrated BEDs with the GSE201336 Adenoma exports. It reads the exported
`per_group/` files directly.

For ABS, CT, EE, GOB, STM and TAC, the script unions the original interval
boundaries with `bedtools merge -d 0`, joining overlapping and book-ended
intervals. It writes separate normal and disease BEDs and
`VUMC_GSE201336_merge_summary.tsv` with input paths, peak counts and overlap
statistics. The released set also includes the VUMC-only Normal TUF, disease
TUF, Adenoma ASC and Serrated SSC BEDs from the VUMC branch.
