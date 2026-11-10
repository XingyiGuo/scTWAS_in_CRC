# TF–target networks

This section follows up significant TWAS genes with two kinds of evidence:
co-expression with the 51 CRC-associated transcription factors (TFs), and
Cistrome peak overlap at each gene's lead model variant.

## Run order

1. **Model SNPs.** In [`correlation_and_occupancy/`](correlation_and_occupancy),
   run `1a_Extract_sigTWASGenes_SNPs_gtf.R` for ADSER and again with
   `TF_NETWORK_GROUP=NL`. It reads the significant TWAS tables, PredictDB
   weights, GENCODE annotation and EUR GWAS statistics.
2. **Correlations and occupancy.** Supply the 51-TF list and cohort expression
   inputs to the two `2a_*` scripts. They calculate TF–gene correlations across
   donor pseudobulk samples for seven normal and nine polyp cell types. Then run
   the matching `2b_*` scripts to select the model SNP with the smallest GWAS
   P value for each gene–cell-type pair and compare its TF occupancy with
   co-expression. Supply the supplementary TWAS workbook, chromosome weight
   files and precomputed Cistrome overlap tables. Set the input and output paths
   at the top of the R scripts.
3. **Pathway analysis.** The [pathway integration scripts](pathway_integration)
   use staged network tables under `TF_NETWORK_DIR`. Set `TF_ENRICHR_XLSX`
   and `TF_PATHWAY_OUT` for pathway integration, or `TF_TWAS_RESULTS` and
   `TF_SPECIFICITY_OUT` for the cell-type specificity permutation test.
   [EnrichR.R](pathway_integration/EnrichR.R) generates the upstream
   cell-type Enrichr TSVs. Integration reads the staged combined workbook.

The [network panels and source data](../11_figures/README.md)
are in section 11. The motif logos and gene–SNP–TF panels are in [`code/Viz/`](../Viz).
