# External command-line tools invoked by the deposited scripts

Counts are the number of deposited scripts that invoke each tool.

| Tool | Scripts | Used for |
|---|---|---|
| Slurm (`sbatch`) | 19 | job submission on the Vanderbilt ACCRE cluster |
| PEER | 17 | probabilistic estimation of expression residuals (covariate factors) |
| tensorQTL | 13 | cis-eQTL mapping per chromosome |
| PLINK | 9 | genotype conversion, LD reference construction |
| BEDTools | 6 | interval operations on peak sets |
| SMR (v1.4.2) | 5 | summary-data-based Mendelian randomisation with HEIDI |
| MACS3 | 5 | peak calling, VUMC scATAC-seq |
| MACS2 | 4 | peak calling, GSE201336 scATAC-seq |
| BCFtools | 4 | VCF handling |
| METAL | 4 | referenced in the GWAS meta-analysis context |
| SAMtools | 2 | alignment file handling |
| tabix | 1 | indexing |
| GCTA | 1 | genotype utilities |
| liftOver | 1 | coordinate conversion |

Versions: SMR v1.4.2 is named in the Methods. The remaining tools were used at the
versions available in the cluster module environment over the analysis period and were not
pinned; where a version could plausibly change a reported number it is named in the
Methods.

## R dependency locks

- `renv.lock` -- the project-level lockfile from the ACCRE analysis environment. This is the
  authoritative R package set; restore with `renv::restore()`. It was generated on
  linux-rocky-9.8 / R-4.5 / x86_64-pc-linux-gnu, so the recorded package *records* are what
  transfers; binaries are resolved per platform on restore.
- `renv_GSE201336_scATAC.lock` -- the narrower lockfile that shipped with the GSE201336
  scATAC-seq peak-calling pipeline.
- `r_packages.txt` -- a human-readable summary of the same package set, grouped by purpose.
