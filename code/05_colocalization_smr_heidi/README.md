# Colocalization and SMR-HEIDI

Both analyses follow up significant TWAS gene–cell-type pairs. Coloc estimates
support for a shared GWAS and cis-eQTL signal. SMR-HEIDI tests the association
using a selected cis-eQTL instrument and checks for heterogeneity. SMR-HEIDI
starts directly from the TWAS hits, regardless of the coloc result.

## Inputs

Edit the [NL and ADSER configs](config) to point to the TWAS significance
tables, ancestry-matched CRC GWAS statistics, tensorQTL cis-eQTL Parquet files
and a GENCODE GTF. Use per-cell-type TWAS tables for coloc and ADSER SMR-HEIDI,
and the combined seven-cell-type table for NL SMR-HEIDI. SMR-HEIDI also needs
source-matched PLINK 2 genotypes and SMR v1.4.2. NL uses the eQTL source that
drove the TWAS result. ADSER uses HTAN polyp eQTLs.

## Run order

1. **Coloc.** Run `01_run_coloc.sh <CT|all> <coloc_config>` with Bash or Slurm.
   Run the corresponding `01a_coloc_*_summarize.R` with the same config to
   combine cell-type results.
2. **SMR-HEIDI.** Run `02a_prepare_SMR_HEIDI.sh <CT|all> <smr_config>`.
   It selects the lowest-P eQTL among harmonized SNPs within ±500 kb of the
   transcription start site. Run `bash 02b_run_SMR_HEIDI.sh <manifest.tsv>`
   on each cell-type manifest to produce the result and status tables.
3. **Report.** Run the matching R script in [`summary/`](summary) with the
   TWAS results root and the previously reported CRC gene CSV. It joins
   coloc and SMR-HEIDI evidence to the significant TWAS pairs.
