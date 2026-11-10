# Single-cell references and cis-eQTL mapping

This section builds cell-type-specific expression and genotype inputs for scTWAS in
normal colorectal tissue and colorectal polyps. The GSE266546 and COLONMAP workflows
also include per-chromosome tensorQTL mapping for downstream colocalization and
SMR-HEIDI.

## Before running

Use a Linux environment with the R and Python packages listed in [env](../../env).
The scripts use Seurat, Harmony and cohort-specific single-cell packages, including
BPCells for Wellcome Sanger. Install PEER (`peertool`), bcftools, PLINK 2 and
tensorQTL for the steps that invoke them. The Slurm wrappers use the original cluster
modules, account and paths. Review those settings and the input and output directories
in each script before submitting jobs.

Stage the cohort count matrices or H5AD files, sample and donor metadata, genotype
files, and the matching GENCODE annotation and dbSNP reference. The GSE266546 and
COLONMAP tensorQTL wrappers take a cohort-specific Bash config. Supply the cell types,
genotype and phenotype path templates, mapping and annotation files, and tensorQTL
settings using the variables read by the scripts. The original cohort configs are
external inputs. The two pseudobulk preparation scripts use the original QTLS helper
directory when available and fall back to the copies in [code/lib/QTLS](../lib/QTLS).

## Gut Cell Atlas and HTAN normal tissue

Run the following from [normal_colon_GSE266546](normal_colon_GSE266546).

| Step | Scripts and handoff |
|---|---|
| Prepare and integrate cells | `01_pref_AC.R` prepares the GSE266546 and HTAN Seurat objects. `02_merge_integration.R` merges the cells and runs Harmony. |
| Cluster and annotate | `03_cluster_annot.R` exports the initial cluster markers and saves `integrated_seurat_final.res0.8.rm23.rds`. Stage that checkpoint and `merged.rds` in the results directory expected by `04_sub_cluster.R` (or align the paths). Step 04 subclusters the retained cells, and `05_annotation_UMAP.R` applies the reviewed labels. |
| Export expression and marker summaries | `07_sample_level_pseudobulk_counts.R` sums raw RNA counts by specimen and reviewed cell type. `06_cell_type_marker_genes_rank.R` exports marker rankings and annotation checks from steps 04 and 05. |
| Match specimens and prepare genotypes | `08_match_N157_expression_genotype_segment_priority_NL.py` selects one eligible expression specimen per donor by segment priority and writes genotype keep lists and match audits. Supply the GEO lookup cache to run without fetching GEO records. Run `g1_HTAN_GSE266546_merge_chr1_22.sh` and then `g2_HTAN_GSE266546_vcf_to_traw.sh` to merge the cohort VCFs and generate genotype QC, PLINK files, dosages and PCA outputs. |
| Adjust expression | Run `Process_peer_HTAN_GSE266546_fixed.py <cell_type>` for each cell type. It expects preprocessed, inverse-normalized expression and a matched age, sex, segment and genotype-PC table. Stage these inputs separately from the raw counts and match lists. |
| Map cis-eQTLs | In `xQTL_tensorqtls/`, run `prepare_tensorqtl_genotype.sh <config>` and `prepare_tensorqtl_beds.sh <config>`, then run `xQTL_per_chrom.sh <cell_type> <chromosome> <config>` for each cell type and autosome. The BED prep uses `PEER_residual_to_tensorqtl_bed.py` to map expression sample IDs to genotype IDs and place genes at their transcription start sites. |

`04_sub_cluster.sh`, `sbatch.sh`, and `process_peer_sbatch.sh` submit subclustering,
VCF merging, and a single cell type's PEER job, respectively.

## Wellcome Sanger healthy tissue

Run the following from [normal_colon_WellcomeSanger](normal_colon_WellcomeSanger).

1. Set up the R environment with `00_install_BPCell.sh` if needed. Use
   `01_h5ad_to_seurat_bpcells.R` or its `run_pipeline.sh` wrapper to convert
   E-MTAB-16986 H5AD to a Seurat object backed by an on-disk BPCells matrix. Keep
   the matrix at the path referenced by the saved RDS.
2. Run `02_extract_healthy_N278_ti_re.R` to select one healthy specimen per donor,
   preferring terminal ileum over rectum. Run
   `03_Healthy_TI_then_rectum_278samples_updated.R` with the LauLabel mapping to
   write epithelial cell metadata and per-donor, per-cell-type pseudobulk counts.
3. For epithelial UMAPs, use `03b_refresh_Epithelial_BPCELLS_with_TAC.R` to refresh
   the BPCells matrix when it needs to include cycling transit-amplifying cells,
   then run `04_umap_epithelial_LauLabel_revised.R`.
4. Run `Genotype_vcf_to_traw_202608.sh` with the selected donors' genotype keep
   list to produce QC-filtered PLINK files, dosage tables and genotype PCs. Run
   `Prepare_pseudobulk_expression.py --input <counts> --gtf <annotation>` for each
   cell type to generate TPM, log-transformed, quantile-normalized and
   inverse-normalized expression. Supply the matched individual covariate table
   with age, sex and five genotype PCs to `Process_peer_WellcomeSanger.py` (or
   submit it with `process_peer_sbatch.sh`).

The deposited Wellcome Sanger scripts end at PEER residuals and prepared genotype
files. The cohort-specific tensorQTL jobs used to generate its nominal eQTL results
are external to this section.

## HTAN COLONMAP polyps

The [polyp annotation](polyps_COLONMAP_annotation) workflow combines DIS, VAL and
PCA83 cells, applies adenoma/serrated and lineage filters, and exports cell metadata
and Seurat objects. Run
`HTAN_ADSER_CellOracleRsave_mergeQC_PCA83_HPtoSER_ASCSSCfilter.R` (or submit its
`HTAN_ADSER_CellOracleRsave_ADSET_9CT_perCT.sh` wrapper), then use
`combine_PCA83_generate_metadata_add3.sh` to add the three late PCA83 specimens
to the cell metadata. The separate `08_umap_ADSER_model_cells.R` and
`09_umap05e_ADSER_polyp_celltype.R` steps use the expression-matched cell metadata
to generate UMAP coordinates and figures. Each has a Slurm wrapper.

In [polyps_COLONMAP](polyps_COLONMAP), run these steps for ASC, SSC, ABS, EE, GOB,
CT, TUF, TAC and STM as appropriate for the analysis:

1. Run `01_match_HTAN_ASC_SSC_expression_genotype_HTANdirect.py <cell_type>` with
   the existing per-cell-type count matrix, HTAN metadata, and the intended
   COLONMAP genotype PSAM and PCA files. This selects eligible AD and SER specimens
   by donor and segment priority, then writes matched counts and covariates.
2. Run `02_prepare_pseudobulk_expression_ASC_SSC_bashconfig.py --config <config>`
   to annotate genes, calculate TPM, apply expression QC, and produce log,
   quantile-normalized and inverse-normalized expression. Point the config at
   the matched counts from step 1.
3. Run `03_Process_peer_HTAN_ASC_SSC_HTANdirect_geneSymbol.py <cell_type>` or its
   `03a_step3_sbatch.sh` wrapper. The supplied PEER script uses age, sex, segment
   indicators and five genotype PCs. AD and SER are specimen eligibility filters
   in the matching script.
4. Run `prepare_tensorqtl_genotype.sh <config>` for the shared genotypes and
   `prepare_tensorqtl_beds_ASC_SSC.sh <config>` for the cell-type-specific
   phenotypes. Set `PEER_TO_BED_SCRIPT` to `05_PEER_residual_to_tensorqtl_bed.py`
   using its full path. Run `06_xQTL_per_chrom_ASC_SSC.sh <cell_type> <chromosome> <config>`
   for each configured cell type and autosome.

Stage the original polyp pseudobulk count matrices and the
`9CT_expression_id_matched_*` metadata exports separately. The annotation workflow
exports cell-level data, and the deposited matching script starts from existing
sample-level count matrices.

## Cohort metadata and downstream outputs

In [cohort_metadata](cohort_metadata), `cells_per_sample_by_celltype.R` builds the
donor-by-cell-type count workbook and CSV using the Supplementary Table 2 donor
list and cohort metadata. `provenance_9CT_expression_id_matched.R` verifies the
polyp cell counts and reconstructs the join from the available metadata and
expression sample order.

- Matched expression, genotype dosages and sample maps support
  [GReX model training](../03_grex_model_training).
- Nominal SNP–gene pair statistics are written as chromosome-specific Parquet files
  for [colocalization and SMR-HEIDI](../05_colocalization_smr_heidi).
- Cell-type annotations, pseudobulk expression, cell counts and UMAP coordinates
  support downstream network analyses and study figures.
