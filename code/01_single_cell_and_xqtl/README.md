# 01 — Single-cell references and cis-eQTL mapping

Three training cohorts, each in its own directory and each following the same arc.

## Position in the pipeline

```
(external data — see Inputs)
        |
        v
   >>> 01_single_cell_and_xqtl <<<
        |
        v
03_grex_model_training  05_colocalization_smr_heidi  06_tf_network  11_figures
```

The arrows above are the logical order of the published workflow, taken from the Methods (`docs/workflow.md`). The file-level links in *From other sections* and *Consumed by* below are the ones actually observed in the deposited code.

## Scripts, in run order (51 files)

Ordered by sub-directory, then by the numeric prefix the script carries in its original working directory. Descriptions are the script's own header line; blank where the script has no header.

| # | File | Lang | What it does | Declared display item |
|---|---|---|---|---|
| 1 | `cohort_metadata/cells_per_sample_by_celltype.R` | R | Build a samples x cell-types matrix of CELL COUNTS for every donor in the scTWAS donor list (Supplementary Table S2), covering both tissues: rows = Ex. | Supplementary Table 2 |
| 2 | `cohort_metadata/check_sample_cells_ADSER.R` | R | Count the cells one ADSER specimen contributes, per cell type, and cross-check that number against the bookkeeping tables that are supposed to agree w. |  |
| 3 | `cohort_metadata/explore_specimen_GCA.R` | R | Interrogate one normal-colon specimen in the GSE266546 + HTAN integrated object: how many cells it contributes per cell type, how deep those cells are. |  |
| 4 | `cohort_metadata/provenance_9CT_expression_id_matched.R` | R | Results/TWAS_scRNA_meta/9CT_expression_id_matched_sample_metadata.csv Results/TWAS_scRNA_meta/9CT_expression_id_matched_cell_metadata.csv Results/TWAS. |  |
| 5 | `normal_colon_GSE266546/01_pref_AC.R` | R |  |  |
| 6 | `normal_colon_GSE266546/02_merge_integration.R` | R |  |  |
| 7 | `normal_colon_GSE266546/03_cluster_annot.R` | R | 03_cluster_plots.R Loads the Harmony-integrated Seurat object, runs clustering, makes UMAP panels (study / sample / condition / cluster), computes per. |  |
| 8 | `normal_colon_GSE266546/04_sub_cluster.R` | R |  |  |
| 9 | `normal_colon_GSE266546/04_sub_cluster.sh` | Shell |  |  |
| 10 | `normal_colon_GSE266546/05_annotation_UMAP.R` | R | 05_annotation.R Annotate clusters using the dominant ORIGINAL cell type in each cluster. |  |
| 11 | `normal_colon_GSE266546/06_cell_type_marker_genes_rank.R` | R | Rank cluster marker genes from: cluster_markers_res0.8.rm23_sub.csv Ranking rules: 1) smaller P value first 2) if tied, larger avg_log2FC first Two ra. |  |
| 12 | `normal_colon_GSE266546/07_sample_level_pseudobulk_counts.R` | R | Create sample-level pseudobulk raw gene-count matrices for each reviewed cell type. |  |
| 13 | `normal_colon_GSE266546/08_match_N157_expression_genotype_segment_priority_NL.py` | Python |  |  |
| 14 | `normal_colon_GSE266546/Process_peer_HTAN_GSE266546_fixed.py` | Python |  |  |
| 15 | `normal_colon_GSE266546/g1_HTAN_GSE266546_merge_chr1_22.sh` | Shell |  |  |
| 16 | `normal_colon_GSE266546/g2_HTAN_GSE266546_vcf_to_traw.sh` | Shell |  |  |
| 17 | `normal_colon_GSE266546/process_peer_sbatch.sh` | Shell |  |  |
| 18 | `normal_colon_GSE266546/sbatch.sh` | Shell |  |  |
| 19 | `normal_colon_GSE266546/xQTL_tensorqtls/PEER_residual_to_tensorqtl_bed.py` | Python |  |  |
| 20 | `normal_colon_GSE266546/xQTL_tensorqtls/prepare_tensorqtl_beds.sh` | Shell |  |  |
| 21 | `normal_colon_GSE266546/xQTL_tensorqtls/prepare_tensorqtl_genotype.sh` | Shell |  |  |
| 22 | `normal_colon_GSE266546/xQTL_tensorqtls/xQTL_per_chrom.sh` | Shell |  |  |
| 23 | `normal_colon_WellcomeSanger/00_install_BPCell.sh` | Shell |  |  |
| 24 | `normal_colon_WellcomeSanger/01_h5ad_to_seurat_bpcells.R` | R | Wellcome Sanger E-MTAB-16986 H5AD -> BPCells on-disk matrix -> Seurat v5 The dataset contains >2 million cells and >7.3 billion non-zero entries. |  |
| 25 | `normal_colon_WellcomeSanger/02_extract_healthy_N278_ti_re.R` | R |  |  |
| 26 | `normal_colon_WellcomeSanger/03_Healthy_TI_then_rectum_278samples_updated.R` | R | Script: 03_Healthy_TI_then_rectum_278samples.R Purpose 1. Load the selected 278-sample Wellcome Sanger Seurat object. 2. |  |
| 27 | `normal_colon_WellcomeSanger/03b_refresh_Epithelial_BPCELLS_with_TAC.R` | R | Script: 03b_refresh_Epithelial_BPCELLS_with_TAC.R >>> RUN ON THE CLUSTER (needs /v5000 and the renv library) <<< Purpose Refresh ONLY the on-disk epit. |  |
| 28 | `normal_colon_WellcomeSanger/04_umap_epithelial_LauLabel_revised.R` | R | Wellcome Sanger healthy epithelium: BPCells -> Seurat -> three labelled UMAPs. |  |
| 29 | `normal_colon_WellcomeSanger/Genotype_vcf_to_traw_202608.sh` | Shell |  |  |
| 30 | `normal_colon_WellcomeSanger/Prepare_pseudobulk_expression.py` | Python |  |  |
| 31 | `normal_colon_WellcomeSanger/Process_peer_WellcomeSanger.py` | Python |  |  |
| 32 | `normal_colon_WellcomeSanger/bpcells_reader.py` | Python | Minimal reader for a BPCells on-disk 'packed-float-matrix-v2' directory. |  |
| 33 | `normal_colon_WellcomeSanger/process_peer_sbatch.sh` | Shell |  |  |
| 34 | `normal_colon_WellcomeSanger/run_pipeline.sh` | Shell |  |  |
| 35 | `polyps_COLONMAP/00_HTAN_ADSER_7cts_COLONMAP95.sh` | Shell | STEP 1: match only the remaining 7 CTs to COLONMAP_95 genotype/PCA. |  |
| 36 | `polyps_COLONMAP/01_match_HTAN_ASC_SSC_expression_genotype_HTANdirect.py` | Python |  |  |
| 37 | `polyps_COLONMAP/02_prepare_pseudobulk_expression_ASC_SSC_bashconfig.py` | Python |  |  |
| 38 | `polyps_COLONMAP/03_Process_peer_HTAN_ASC_SSC_HTANdirect_geneSymbol.py` | Python |  |  |
| 39 | `polyps_COLONMAP/03a_step3_sbatch.sh` | Shell |  |  |
| 40 | `polyps_COLONMAP/05_PEER_residual_to_tensorqtl_bed.py` | Python |  |  |
| 41 | `polyps_COLONMAP/06_xQTL_per_chrom_ASC_SSC.sh` | Shell |  |  |
| 42 | `polyps_COLONMAP/VUMC_HTAN_3h5ad_metadata_clean.R` | R | Explore metadata from three VUMC HTAN h5ad files Files: 1. VUMC_HTAN_DIS_EPI_V2.h5ad 2. VUMC_HTAN_VAL_EPI_V2.h5ad 3. |  |
| 43 | `polyps_COLONMAP/prepare_tensorqtl_beds_ASC_SSC.sh` | Shell |  |  |
| 44 | `polyps_COLONMAP/prepare_tensorqtl_genotype.sh` | Shell |  |  |
| 45 | `polyps_COLONMAP_annotation/08_umap_ADSER_model_cells.R` | R |  |  |
| 46 | `polyps_COLONMAP_annotation/08a_umap_ADSER_model_cells_sbatch.sh` | Shell |  |  |
| 47 | `polyps_COLONMAP_annotation/09_umap05e_ADSER_polyp_celltype.R` | R | 05e-style UMAP for the VUMC HTAN AD/SER polyp panel (DIS + VAL + PCA83), drawn over the cells that built the 9 scTWAS cell-type models. |  |
| 48 | `polyps_COLONMAP_annotation/09a_umap05e_ADSER_polyp_celltype_sbatch.sh` | Shell |  |  |
| 49 | `polyps_COLONMAP_annotation/HTAN_ADSER_CellOracleRsave_ADSET_9CT_perCT.sh` | Shell |  |  |
| 50 | `polyps_COLONMAP_annotation/HTAN_ADSER_CellOracleRsave_mergeQC_PCA83_HPtoSER_ASCSSCfilter.R` | R |  |  |
| 51 | `polyps_COLONMAP_annotation/combine_PCA83_generate_metadata_add3.sh` | Shell |  |  |

## Inputs

**Not produced anywhere in this deposit** (10) — cohort data, public resources, or output of a step that is not deposited:

- `05_Epithelial_cell_metadata_LauLabel.csv.gz`
- `9CT_expression_id_match_summary.csv`
- `9CT_expression_id_matched_sample_metadata.csv`
- `HTAN_DISVAL_PCA83_ADSER_HPtoSER_AllCT_ASCSSCfiltered_CellOracleR_metadata.csv`
- `PCA83samples.h5ad`
- `PCA83samples.metadata.csv`
- `integrated_seurat_final.res0.8.rm23.rds`
- `matrix_source.txt`
- `objH_DISVAL_PCA83_EPI_afterCell200_Gene3QC_beforeConditionFilter.rds`
- `residuals.csv`

## Outputs

| File | Written by | Consumed by |
|---|---|---|
| `01_LauLabel_mapping_file_used.csv` | `normal_colon_WellcomeSanger/03_Healthy_TI_then_rectum_278samples_updated.R` | — (terminal) |
| `01_h5ad_file_summary.tsv` | `polyps_COLONMAP/VUMC_HTAN_3h5ad_metadata_clean.R` | — (terminal) |
| `02_all_metadata_column_summary.tsv` | `polyps_COLONMAP/VUMC_HTAN_3h5ad_metadata_clean.R` | — (terminal) |
| `02_finelabel_to_LauLabel_mapping_audit.csv` | `normal_colon_WellcomeSanger/03_Healthy_TI_then_rectum_278samples_updated.R` | — (terminal) |
| `03_original_ad.obs_column_presence.tsv` | `polyps_COLONMAP/VUMC_HTAN_3h5ad_metadata_clean.R` | — (terminal) |
| `03_unmapped_finelabels.csv` | `normal_colon_WellcomeSanger/03_Healthy_TI_then_rectum_278samples_updated.R` | — (terminal) |
| `03b_BPCELLS_refresh_log.csv` | `normal_colon_WellcomeSanger/03b_refresh_Epithelial_BPCELLS_with_TAC.R` | — (terminal) |
| `03b_sessionInfo.txt` | `normal_colon_WellcomeSanger/03b_refresh_Epithelial_BPCELLS_with_TAC.R` | — (terminal) |
| `04_DROPPED_finelabels.csv` | `normal_colon_WellcomeSanger/03_Healthy_TI_then_rectum_278samples_updated.R` | — (terminal) |
| `04_standardized_metadata_column_presence.tsv` | `polyps_COLONMAP/VUMC_HTAN_3h5ad_metadata_clean.R` | — (terminal) |
| `04_umap_R_sessionInfo.txt` | `normal_colon_WellcomeSanger/04_umap_epithelial_LauLabel_revised.R` | — (terminal) |
| `05_categorical_metadata_counts.tsv` | `polyps_COLONMAP/VUMC_HTAN_3h5ad_metadata_clean.R` | — (terminal) |
| `06_important_metadata_counts.tsv` | `polyps_COLONMAP/VUMC_HTAN_3h5ad_metadata_clean.R` | — (terminal) |
| `06_subject_x_LauLabel_ncells.csv` | `normal_colon_WellcomeSanger/03_Healthy_TI_then_rectum_278samples_updated.R` | — (terminal) |
| `07_LauLabel_cells_subjects_summary.csv` | `normal_colon_WellcomeSanger/03_Healthy_TI_then_rectum_278samples_updated.R` | — (terminal) |
| `07_sample_level_overview.tsv` | `polyps_COLONMAP/VUMC_HTAN_3h5ad_metadata_clean.R` | — (terminal) |
| `08_Epithelial_subject_metadata.csv` | `normal_colon_WellcomeSanger/03_Healthy_TI_then_rectum_278samples_updated.R` | — (terminal) |
| `10_LauLabel_pseudobulk_manifest.csv` | `normal_colon_WellcomeSanger/03_Healthy_TI_then_rectum_278samples_updated.R` | — (terminal) |
| `12_sessionInfo.txt` | `normal_colon_WellcomeSanger/03_Healthy_TI_then_rectum_278samples_updated.R` | — (terminal) |
| `HTAN_DISVAL_PCA83_ADSER_HPtoSER_ASCSSC_lineage_drop_AUDIT.csv` | `polyps_COLONMAP_annotation/HTAN_ADSER_CellOracleRsave_mergeQC_PCA83_HPtoSER_ASCSSCfilter.R` | — (terminal) |
| `HTAN_DISVAL_PCA83_ADSER_HPtoSER_AllCT_ASCSSCfiltered_CellOracleR_metadata.addHTA11_6407_8890_8930.csv` | `polyps_COLONMAP_annotation/combine_PCA83_generate_metadata_add3.sh` | — (terminal) |
| `HTAN_DISVAL_PCA83_ADSER_HPtoSER_AllCT_ASCSSCfiltered_cell_counts_AUDIT.csv` | `polyps_COLONMAP_annotation/HTAN_ADSER_CellOracleRsave_mergeQC_PCA83_HPtoSER_ASCSSCfilter.R` | — (terminal) |
| `PCA83_Polyp_ADTYPE_HPTYPE_to_ADSER_HPtoSER_AUDIT.csv` | `polyps_COLONMAP_annotation/HTAN_ADSER_CellOracleRsave_mergeQC_PCA83_HPtoSER_ASCSSCfilter.R` | — (terminal) |
| `QC_model_cells_absent_from_object.csv` | `polyps_COLONMAP_annotation/08_umap_ADSER_model_cells.R` | — (terminal) |
| `UMAP_ADSER_model_cells_coords.csv` | `polyps_COLONMAP_annotation/08_umap_ADSER_model_cells.R` | `09_umap05e_ADSER_polyp_celltype.R` |
| `UMAP_ADSER_polyp_coords.csv` | `polyps_COLONMAP_annotation/09_umap05e_ADSER_polyp_celltype.R` | — (terminal) |
| `_celltype_composition.csv` | `polyps_COLONMAP_annotation/09_umap05e_ADSER_polyp_celltype.R` | — (terminal) |
| `barcodes.tsv` | `normal_colon_GSE266546/01_pref_AC.R` | `VUMC_HTAN_3h5ad_metadata_clean.R` |
| `cell_type_colors.csv` | `normal_colon_WellcomeSanger/04_umap_epithelial_LauLabel_revised.R` | — (terminal) |
| `cells_by_celltype_cohort.csv` | `polyps_COLONMAP_annotation/08_umap_ADSER_model_cells.R` | — (terminal) |
| `cells_per_cluster_by_condition.res0.8.rm23.csv` | `normal_colon_GSE266546/03_cluster_annot.R` | — (terminal) |
| `cells_per_cluster_by_condition.res0.8.rm23_sub.csv` | `normal_colon_GSE266546/04_sub_cluster.R` | — (terminal) |
| `cells_per_cluster_by_study.res0.8.rm23.csv` | `normal_colon_GSE266546/03_cluster_annot.R` | — (terminal) |
| `cells_per_cluster_by_study.res0.8.rm23_sub.csv` | `normal_colon_GSE266546/04_sub_cluster.R` | — (terminal) |
| `cells_per_final_label.csv` | `normal_colon_WellcomeSanger/01_h5ad_to_seurat_bpcells.R` | — (terminal) |
| `cells_per_sample.csv` | `normal_colon_WellcomeSanger/01_h5ad_to_seurat_bpcells.R` | — (terminal) |
| `cluster_LauLabel_composition.csv` | `normal_colon_WellcomeSanger/04_umap_epithelial_LauLabel_revised.R` | — (terminal) |
| `cluster_by_original_CellType_cb_annot_pct_res0.8.rm23_sub.csv` | `normal_colon_GSE266546/05_annotation_UMAP.R` | — (terminal) |
| `cluster_by_original_CellType_cb_res0.8.rm23_sub.csv` | `normal_colon_GSE266546/05_annotation_UMAP.R` | — (terminal) |
| `cluster_celltype_composition_res0.8.rm23_sub.csv` | `normal_colon_GSE266546/05_annotation_UMAP.R` | — (terminal) |
| `cluster_dominant_celltype_res0.8.rm23_sub.csv` | `normal_colon_GSE266546/05_annotation_UMAP.R` | — (terminal) |
| `cluster_majority_LauLabel.csv` | `normal_colon_WellcomeSanger/04_umap_epithelial_LauLabel_revised.R` | — (terminal) |
| `cluster_markers_res0.8.rm23.csv` | `normal_colon_GSE266546/03_cluster_annot.R` | — (terminal) |
| `cluster_markers_res0.8.rm23_sub.csv` | `normal_colon_GSE266546/04_sub_cluster.R` | — (terminal) |
| `cluster_reviewed_annotation_res0.8.rm23_sub.csv` | `normal_colon_GSE266546/05_annotation_UMAP.R` | — (terminal) |
| `counts.qc.csv` | `polyps_COLONMAP/02_prepare_pseudobulk_expression_ASC_SSC_bashconfig.py` | — (terminal) |
| `counts_matrix_res0.8.rm23_sub.rds` | `normal_colon_GSE266546/04_sub_cluster.R` | — (terminal) |
| `features.tsv` | `normal_colon_GSE266546/01_pref_AC.R` | `VUMC_HTAN_3h5ad_metadata_clean.R`, `HTAN_ADSER_CellOracleRsave_mergeQC_PCA83_HPtoSER_ASCSSCfilter.R` |
| `fine_annotation_mapping.csv` | `normal_colon_WellcomeSanger/04_umap_epithelial_LauLabel_revised.R` | — (terminal) |
| `integrated_seurat_final.res0.8.rm23_metadata_sub.csv` | `normal_colon_GSE266546/04_sub_cluster.R` | — (terminal) |
| `integrated_seurat_final.res0.8.rm23_sub.rds` | `normal_colon_GSE266546/04_sub_cluster.R` | — (terminal) |
| `integrated_seurat_final1.res0.8.rm23_sub.rds` | `normal_colon_GSE266546/04_sub_cluster.R` | — (terminal) |
| `integrated_seurat_final1.res0.8.rm23_sub_metadata.csv` | `normal_colon_GSE266546/04_sub_cluster.R` | — (terminal) |
| `integrated_seurat_final2.res0.8.rm23_sub.rds` | `normal_colon_GSE266546/05_annotation_UMAP.R` | — (terminal) |
| `integrated_seurat_final2_reviewed.res0.8.rm23_sub.rds` | `normal_colon_GSE266546/05_annotation_UMAP.R` | `explore_specimen_GCA.R` |
| `integrated_seurat_final2_reviewed.res0.8.rm23_sub_metadata.csv` | `normal_colon_GSE266546/05_annotation_UMAP.R` | `explore_specimen_GCA.R` |
| `integrated_seurat_postharmony.rds` | `normal_colon_GSE266546/02_merge_integration.R` | `03_cluster_annot.R` |
| `integrated_seurat_postharmony_sub.rds` | `normal_colon_GSE266546/04_sub_cluster.R` | — (terminal) |
| `merged.rds` | `normal_colon_GSE266546/02_merge_integration.R` | `04_sub_cluster.R` |
| `merged_metadata.csv` | `normal_colon_GSE266546/02_merge_integration.R` | — (terminal) |
| `merged_postSCT_PCA_sub.rds` | `normal_colon_GSE266546/04_sub_cluster.R` | — (terminal) |
| `merged_sub.rds` | `normal_colon_GSE266546/04_sub_cluster.R` | — (terminal) |
| `merged_sub_metadata.csv` | `normal_colon_GSE266546/04_sub_cluster.R` | — (terminal) |
| `metadata_columns.csv` | `normal_colon_WellcomeSanger/01_h5ad_to_seurat_bpcells.R` | — (terminal) |
| `obj266_metadata_pref_AC_immu.csv` | `normal_colon_GSE266546/02_merge_integration.R` | — (terminal) |
| `obj266_metadata_preqc.csv` | `normal_colon_GSE266546/01_pref_AC.R` | — (terminal) |
| `obj266_metadata_preqc_pref_AC.csv` | `normal_colon_GSE266546/01_pref_AC.R` | — (terminal) |
| `obj266_pref_AC_immu.rds` | `normal_colon_GSE266546/02_merge_integration.R` | — (terminal) |
| `obj266_preqc.rds` | `normal_colon_GSE266546/01_pref_AC.R` | — (terminal) |
| `obj266_preqc_pref_AC.rds` | `normal_colon_GSE266546/01_pref_AC.R` | `02_merge_integration.R` |
| `objH_ADSER_model_cells_umap.rds` | `polyps_COLONMAP_annotation/08_umap_ADSER_model_cells.R` | — (terminal) |
| `objH_immu.rds` | `normal_colon_GSE266546/02_merge_integration.R` | — (terminal) |
| `objH_metadata_immu.csv` | `normal_colon_GSE266546/02_merge_integration.R` | — (terminal) |
| `objH_metadata_preqc.csv` | `normal_colon_GSE266546/01_pref_AC.R` | — (terminal) |
| `objH_preqc.rds` | `normal_colon_GSE266546/01_pref_AC.R` | `02_merge_integration.R` |
| `obs.tsv` | `normal_colon_GSE266546/01_pref_AC.R` | `VUMC_HTAN_3h5ad_metadata_clean.R` |
| `obs_columns.tsv` | `polyps_COLONMAP/VUMC_HTAN_3h5ad_metadata_clean.R` | — (terminal) |
| `sampling_by_celltype.csv` | `polyps_COLONMAP_annotation/08_umap_ADSER_model_cells.R` | — (terminal) |
| `scTWAS_cells_per_sample_by_celltype.csv` | `cohort_metadata/cells_per_sample_by_celltype.R` | — (terminal) |
| `sessionInfo.txt` | `polyps_COLONMAP_annotation/08_umap_ADSER_model_cells.R` | — (terminal) |
| `sessionInfo_09.txt` | `polyps_COLONMAP_annotation/09_umap05e_ADSER_polyp_celltype.R` | — (terminal) |
| `sessionInfo_h5ad_to_seurat.txt` | `normal_colon_WellcomeSanger/01_h5ad_to_seurat_bpcells.R` | — (terminal) |
| `tpmlog2.qc.csv` | `polyps_COLONMAP/02_prepare_pseudobulk_expression_ASC_SSC_bashconfig.py` | — (terminal) |
| `tpmlog2.qc.qn.csv` | `polyps_COLONMAP/02_prepare_pseudobulk_expression_ASC_SSC_bashconfig.py` | — (terminal) |
| `umap_WellcomeSanger_epithelial_clusters_LauLabel_sourcedata.csv` | `normal_colon_WellcomeSanger/04_umap_epithelial_LauLabel_revised.R` | — (terminal) |

## How this file was built

Automatically generated 2026-10-19 from MANIFEST.csv, the Methods order in `docs/workflow.md`, and the file reads and writes parsed out of the deposited scripts themselves. The extraction reads only string literals and path variables that appear at a read or write call site, so three classes of file are invisible to it: paths assembled from shell variables (most `.sh` and `.sbatch` wrappers), paths built by string formatting at run time, and files handled by external tools that the scripts invoke rather than open themselves (tensorQTL, PLINK, MACS, SMR, PEER). The lists above are therefore a lower bound, not a complete inventory. `MANIFEST.csv` at the repository root is authoritative for the file list itself.
