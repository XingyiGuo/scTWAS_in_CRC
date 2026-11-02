# manuscript_text — Manuscript prose and reported-number recomputation

## Position in the pipeline

```
(external data — see Inputs)
        |
        v
   >>> manuscript_text <<<
        |
        v
   (terminal: display items / supplementary tables)
```

The arrows above are the logical order of the published workflow, taken from the Methods (`docs/workflow.md`). The file-level links in *From other sections* and *Consumed by* below are the ones actually observed in the deposited code.

## Scripts, in run order (1 files)

Ordered by sub-directory, then by the numeric prefix the script carries in its original working directory. Descriptions are the script's own header line; blank where the script has no header.

| # | File | Lang | What it does | Declared display item |
|---|---|---|---|---|
| 1 | `build_GWAS_meta_NG.py` | Python | -*- coding: utf-8 -*- |  |

## Inputs

No input file could be resolved by name: this section's scripts build their paths from shell or configuration variables. See the caveat at the foot of this file.

## Outputs

No output file could be resolved by name; see the caveat below.

## How this file was built

Automatically generated 2026-10-19 from MANIFEST.csv, the Methods order in `docs/workflow.md`, and the file reads and writes parsed out of the deposited scripts themselves. The extraction reads only string literals and path variables that appear at a read or write call site, so three classes of file are invisible to it: paths assembled from shell variables (most `.sh` and `.sbatch` wrappers), paths built by string formatting at run time, and files handled by external tools that the scripts invoke rather than open themselves (tensorQTL, PLINK, MACS, SMR, PEER). The lists above are therefore a lower bound, not a complete inventory. `MANIFEST.csv` at the repository root is authoritative for the file list itself.
