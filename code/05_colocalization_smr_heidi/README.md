# 05 — Colocalization and SMR-HEIDI

`01_coloc_gwas_xqtl_sigSource_peakSNP.R` runs coloc with the approximate Bayes factor
approach for the significant gene-cell-type pairs against the source-matched cis-eQTL
data, driven per cell type by a config file; `01a_*_summarize.R` collapse the per-cell-type
output for the normal-colon (7 cell types) and polyp (9 populations) model sets.
`02_SMR_HEIDI_inputs_targetTopEQTL_500kb_directTWAS.R` builds the SMR inputs, selecting the
top cis-eQTL per gene as the instrument and writing the locus-specific PLINK reference;
`02a`/`02b` prepare and run SMR v1.4.2 with HEIDI.

## Position in the pipeline

```
04_association_testing  01_single_cell_and_xqtl
        |
        v
   >>> 05_colocalization_smr_heidi <<<
        |
        v
   (terminal: display items / supplementary tables)
```

The arrows above are the logical order of the published workflow, taken from the Methods (`docs/workflow.md`). The file-level links in *From other sections* and *Consumed by* below are the ones actually observed in the deposited code.

## Scripts, in run order (7 files)

Ordered by sub-directory, then by the numeric prefix the script carries in its original working directory. Descriptions are the script's own header line; blank where the script has no header.

| # | File | Lang | What it does | Declared display item |
|---|---|---|---|---|
| 1 | `01_coloc_gwas_xqtl_sigSource_peakSNP.R` | R |  |  |
| 2 | `01_run_coloc.sh` | Shell |  |  |
| 3 | `01a_coloc_ADSER_9cts_summarize.R` | R |  |  |
| 4 | `01a_coloc_NL_7cts_summarize.R` | R |  |  |
| 5 | `02_SMR_HEIDI_inputs_targetTopEQTL_500kb_directTWAS.R` | R |  |  |
| 6 | `02a_prepare_SMR_HEIDI.sh` | Shell |  |  |
| 7 | `02b_run_SMR_HEIDI.sh` | Shell |  |  |

## Inputs

No input file could be resolved by name: this section's scripts build their paths from shell or configuration variables. See the caveat at the foot of this file.

## Outputs

No output file could be resolved by name; see the caveat below.

## How this file was built

Automatically generated 2026-10-19 from MANIFEST.csv, the Methods order in `docs/workflow.md`, and the file reads and writes parsed out of the deposited scripts themselves. The extraction reads only string literals and path variables that appear at a read or write call site, so three classes of file are invisible to it: paths assembled from shell variables (most `.sh` and `.sbatch` wrappers), paths built by string formatting at run time, and files handled by external tools that the scripts invoke rather than open themselves (tensorQTL, PLINK, MACS, SMR, PEER). The lists above are therefore a lower bound, not a complete inventory. `MANIFEST.csv` at the repository root is authoritative for the file list itself.

## Added in this release

These files were merged in from the project repository after the first deposit was assembled, so they are not described in the run order above.

The coloc and SMR/HEIDI configuration files, the all-cell-type driver, and the `summary/` scripts that join TWAS, coloc and SMR into the reported gene sets.

- `config/SMR_HEIDI_ADSER_targetTopEQTL_500kb.config.txt`
- `config/SMR_HEIDI_NL_targetTopEQTL_500kb.config.txt`
- `config/coloc_peakSNP_ADSER.config.txt`
- `config/coloc_peakSNP_NL.config.txt`
- `run_all_celltypes.sh`
- `summary/Summarize_ADSER_TWAS_coloc_SMR_EUR_then_Multi_with_reportedCRC.R`
- `summary/Summarize_NL_TWAS_coloc_SMR_EUR_then_Multi_with_reportedCRC_v2.R`
