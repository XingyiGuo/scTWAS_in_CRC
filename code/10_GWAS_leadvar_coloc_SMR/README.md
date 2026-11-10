# GWAS lead-variant follow-up

This workflow tests 275 CRC GWAS lead variants for associations with cell-type-specific gene expression in Wellcome Sanger normal colon, VUMC normal colon, and VUMC polyps. Supplementary Table 12 lists the variants. Supplementary Table 13 reports the cis results. [Section 05](../05_colocalization_smr_heidi/README.md) follows up significant TWAS gene–cell-type pairs.

## Run order

1. Run the cohort's `LeadVariants_Profile_lm_complete_*.R` for each cell type using PEER-residualized pseudobulk expression, matched genotypes, and the lead-variant list.
2. Run [`LeadVariants_BF_significant.R`](LeadVariants_BF_significant.R) with `WellcomeSanger`, `VUMC_HTAN_GSE266546`, or `VUMC_Polyps` to select Bonferroni-significant pairs within each cell type. Cis pairs involve genes overlapping the ±500 kb window around the lead variant.
3. Run [`GWAS_leadvar_coloc_protein.sh`](GWAS_leadvar_coloc_protein.sh) with the cohort (`WellcomeSanger`, `VUMC`, or `Polyps`) and `cisprotein`. In parallel, run the cohort's `GWAS_leadvar_SMR_HEIDI_inputs*.R` with `cisprotein`, then pass its manifest to [`GWAS_leadvar_run_SMR_HEIDI.sh`](GWAS_leadvar_run_SMR_HEIDI.sh). Both analyses use regional GWAS and tensorQTL statistics. SMR–HEIDI also uses cohort genotypes for LD.
4. Merge the results with `GWAS_leadvar_coloc_SMR_merge.R`, archived at `W:/VUMC_projects/scTWAS/GWAS_leadvar_coloc_SMR/`. Table 13 includes 17 of the 18 significant VUMC normal-colon pairs. The omitted pair had neither follow-up result.

The scripts use ACCRE cluster paths. Set the input and output paths for your environment before running.
