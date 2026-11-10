# Software setup

Install the Python packages from [`python_requirements.txt`](python_requirements.txt).
Use the R environment notes below and install the external tools for the sections you plan to run.

| Tool | Used for |
|---|---|
| Slurm (`sbatch`) | Job submission on Vanderbilt ACCRE |
| PEER and tensorQTL | Expression residualization and cis-eQTL mapping |
| PLINK, BCFtools, GCTA and tabix | Genotype processing, indexing and LD references |
| MACS2, MACS3, BEDTools and SAMtools | scATAC-seq peak calling and interval or alignment processing |
| MetaXcan / S-PrediXcan | TWAS association testing |
| SMR v1.4.2 | SMR–HEIDI tests |
| METAL | Upstream GWAS meta-analysis |
| liftOver | Coordinate conversion |
| Node.js and npm | Figure 6 radial tree, with `d3` and `jsdom` |

Versions: SMR v1.4.2 is named in the Methods. The remaining tools were used at the
versions available in the cluster module environment over the analysis period and were not
pinned; where a version could plausibly change a reported number it is named in the
Methods.

## R dependency locks

- [`renv.lock`](renv.lock) records the project environment from ACCRE (Linux, R 4.5.0).
  From the repository root, restore it with `renv::restore(lockfile = "env/renv.lock")`.
- [`renv_GSE201336_scATAC.lock`](renv_GSE201336_scATAC.lock) records the separate
  GSE201336 scATAC-seq environment. It pins different versions for 15 shared packages.
- [`r_packages.txt`](r_packages.txt) lists R packages used by the scripts, including
  packages outside the project lockfile.
