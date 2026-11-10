library(ACAT)
library(data.table)
library(dplyr)

args <- commandArgs(trailingOnly = TRUE)

# Cohort can be supplied as the second command-line argument.
# Examples:
#   Rscript Genes_CauchyMeta_NL_cohort_fixed.R all EUR
#   Rscript Genes_CauchyMeta_NL_cohort_fixed.R all Multi
# Defaults to EUR for backward compatibility.
cohort <- if (length(args) >= 2) args[2] else "EUR"

allowed_cohorts <- c("EUR", "Multi")
if (!cohort %in% allowed_cohorts) {
  stop(
    "Unknown cohort: ", cohort,
    ". Allowed cohorts: ", paste(allowed_cohorts, collapse = ", ")
  )
}
combine_acat <- function(p1, p2) {
  p <- c(p1, p2)
  p <- p[!is.na(p)]

  if (length(p) == 0) return(NA_real_)
  if (length(p) == 1) return(as.numeric(p[1]))

  ACAT(p)
}

# -----------------------------------------------------------------------------
# Cell types
# Usage:
#   Rscript Genes_CauchyMeta_NL_cohort_fixed.R all EUR
#   Rscript Genes_CauchyMeta_NL_cohort_fixed.R 1 EUR
#   Rscript Genes_CauchyMeta_NL_cohort_fixed.R ABS Multi
# where 1-7 correspond to ABS, CT, EE, GOB, STM, TAC, TUF.
# -----------------------------------------------------------------------------
CTlist <- c("ABS", "CT", "EE", "GOB", "STM", "TAC", "TUF")

if (length(args) < 1) {
  stop("Please provide 'all', a cell type index (1-7), or a cell type name: ",
       paste(CTlist, collapse = ", "))
}

arg1 <- toupper(args[1])

if (arg1 == "ALL") {
  CTs_to_run <- CTlist
} else if (grepl("^[0-9]+$", args[1])) {
  idx <- suppressWarnings(as.integer(args[1]))
  if (is.na(idx) || idx < 1 || idx > length(CTlist)) {
    stop("Cell type index must be between 1 and ", length(CTlist))
  }
  CTs_to_run <- CTlist[idx]
} else {
  if (!arg1 %in% CTlist) {
    stop("Unknown cell type: ", arg1,
         ". Allowed: all, ", paste(CTlist, collapse = ", "))
  }
  CTs_to_run <- arg1
}

# -----------------------------------------------------------------------------
# Input/output paths
# -----------------------------------------------------------------------------
BASE_DIR <- Sys.getenv(
  "ASSO_OUTPUT_DIR",
  "/data/l2_bioinfo1/liq17/scTWAS/TWAS_HTAN_GSE266546_WellcomeSanger_res/SPrediXcan_peakSNP_association"
)
OUT_DIR <- file.path(BASE_DIR, paste0("Cauchy_Meta_NL_WellcomeSanger_NL_", cohort))
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

run_meta <- function(CT) {
  HTAN_NL_file <- file.path(
    BASE_DIR,
    paste0("VUMC_HTAN_GSE266546_NL_", CT, "_peakSNP.",cohort,".TWAS.csv")
  )

  Wellcome_NL_file <- file.path(
    BASE_DIR,
    paste0("WellcomeSanger_NL_", CT, "_peakSNP.",cohort,".TWAS.csv")
  )

  out_file <- file.path(
    OUT_DIR,
    paste0("HTAN_NL_WellcomeSanger_NL_", CT, "_peakSNP.",cohort,".CauchyMeta.csv")
  )

  cat("============================================================\n")
  cat("Cauchy meta-analysis: VUMC HTAN NL + WellcomeSanger NL [", cohort, "]\n", sep = "")
  cat("============================================================\n")
  cat("Cohort           :", cohort, "\n")
  cat("Cell type        :", CT, "\n")
  cat("HTAN NL file     :", HTAN_NL_file, "\n")
  cat("Wellcome NL file :", Wellcome_NL_file, "\n")
  cat("Output           :", out_file, "\n\n")

  if (!file.exists(HTAN_NL_file)) {
    stop("Missing HTAN NL file: ", HTAN_NL_file)
  }
  if (!file.exists(Wellcome_NL_file)) {
    stop("Missing WellcomeSanger NL file: ", Wellcome_NL_file)
  }

  # Read gene-level SPrediXcan TWAS results.
  HTAN_NL_df <- as.data.frame(fread(HTAN_NL_file))
  Wellcome_NL_df <- as.data.frame(fread(Wellcome_NL_file))

  # Required gene-level SPrediXcan columns.
  # Only models with R2_global > 0.01 are retained for meta-analysis.
  required_cols <- c("gene", "pvalue", "zscore", "R2_global")

  missing_HTAN <- setdiff(required_cols, colnames(HTAN_NL_df))
  missing_Wellcome <- setdiff(required_cols, colnames(Wellcome_NL_df))

  if (length(missing_HTAN) > 0) {
    stop(
      "HTAN NL file is missing required column(s): ",
      paste(missing_HTAN, collapse = ", "),
      "\nAvailable columns: ",
      paste(colnames(HTAN_NL_df), collapse = ", ")
    )
  }

  if (length(missing_Wellcome) > 0) {
    stop(
      "WellcomeSanger NL file is missing required column(s): ",
      paste(missing_Wellcome, collapse = ", "),
      "\nAvailable columns: ",
      paste(colnames(Wellcome_NL_df), collapse = ", ")
    )
  }

  # Filter prediction models by global model performance before merging.
  n_HTAN_raw <- nrow(HTAN_NL_df)
  n_Wellcome_raw <- nrow(Wellcome_NL_df)

  HTAN_NL_df <- HTAN_NL_df %>%
    filter(!is.na(R2_global), R2_global > 0.01)

  Wellcome_NL_df <- Wellcome_NL_df %>%
    filter(!is.na(R2_global), R2_global > 0.01)

  cat("R2_global filter  : > 0.01\n")
  cat("HTAN NL kept      :", nrow(HTAN_NL_df), "/", n_HTAN_raw, "genes\n")
  cat("Wellcome NL kept  :", nrow(Wellcome_NL_df), "/", n_Wellcome_raw, "genes\n\n")

  # SPrediXcan should already contain one row per gene. If duplicates occur,
  # retain the smallest-p row to avoid an accidental many-to-many merge.
  if (anyDuplicated(HTAN_NL_df$gene)) {
    warning(
      "Duplicate genes detected in HTAN NL input for ", CT,
      "; keeping the row with the smallest pvalue per gene."
    )
    HTAN_NL_df <- HTAN_NL_df %>%
      arrange(pvalue) %>%
      distinct(gene, .keep_all = TRUE)
  }

  if (anyDuplicated(Wellcome_NL_df$gene)) {
    warning(
      "Duplicate genes detected in WellcomeSanger NL input for ", CT,
      "; keeping the row with the smallest pvalue per gene."
    )
    Wellcome_NL_df <- Wellcome_NL_df %>%
      arrange(pvalue) %>%
      distinct(gene, .keep_all = TRUE)
  }

  # Add dataset suffixes while retaining 'gene' as the merge key.
  colnames(HTAN_NL_df) <- ifelse(
    colnames(HTAN_NL_df) == "gene",
    "gene",
    paste0(colnames(HTAN_NL_df), "_HTAN_NL")
  )

  colnames(Wellcome_NL_df) <- ifelse(
    colnames(Wellcome_NL_df) == "gene",
    "gene",
    paste0(colnames(Wellcome_NL_df), "_WellcomeSanger_NL")
  )

  # Merge by gene and Cauchy-combine the two cohort-specific TWAS p-values.
  # If a gene occurs in only one dataset, retain that dataset's p-value.
  meta_df <- merge(
    HTAN_NL_df,
    Wellcome_NL_df,
    by = "gene",
    all = TRUE
  )

  meta_df <- meta_df %>%
    rowwise() %>%
    mutate(
      META_P = combine_acat(
        pvalue_HTAN_NL,
        pvalue_WellcomeSanger_NL
      )
    ) %>%
    ungroup() %>%
    mutate(
      META_SOURCE = case_when(
        !is.na(pvalue_HTAN_NL) & !is.na(pvalue_WellcomeSanger_NL) ~
          "VUMC_HTAN_GSE266546_NL+WellcomeSanger_NL",
        !is.na(pvalue_HTAN_NL) ~ "VUMC_HTAN_GSE266546_NL_only",
        !is.na(pvalue_WellcomeSanger_NL) ~ "WellcomeSanger_NL_only",
        TRUE ~ NA_character_
      ),
      `P driven` = case_when(
        !is.na(pvalue_HTAN_NL) &
          (is.na(pvalue_WellcomeSanger_NL) |
             pvalue_HTAN_NL <= pvalue_WellcomeSanger_NL) ~
          "VUMC_HTAN_GSE266546",
        !is.na(pvalue_WellcomeSanger_NL) ~ "WellcomeSanger",
        TRUE ~ NA_character_
      ),
      TWAS_zscore_direction = case_when(
        !is.na(zscore_HTAN_NL) & !is.na(zscore_WellcomeSanger_NL) &
          sign(zscore_HTAN_NL) == sign(zscore_WellcomeSanger_NL) ~
          "Same_direction",
        !is.na(zscore_HTAN_NL) & !is.na(zscore_WellcomeSanger_NL) &
          sign(zscore_HTAN_NL) != sign(zscore_WellcomeSanger_NL) ~
          "Different_direction",
        TRUE ~ NA_character_
      )
    ) %>%
    arrange(META_P)

  fwrite(meta_df, out_file)

  n_both <- sum(
    !is.na(meta_df$pvalue_HTAN_NL) &
      !is.na(meta_df$pvalue_WellcomeSanger_NL)
  )

  n_both_zscores <- sum(
    !is.na(meta_df$zscore_HTAN_NL) &
      !is.na(meta_df$zscore_WellcomeSanger_NL)
  )

  n_same_direction <- sum(
    meta_df$TWAS_zscore_direction == "Same_direction",
    na.rm = TRUE
  )

  n_different_direction <- sum(
    meta_df$TWAS_zscore_direction == "Different_direction",
    na.rm = TRUE
  )

  cat("Genes in HTAN NL           :", nrow(HTAN_NL_df), "\n")
  cat("Genes in WellcomeSanger NL :", nrow(Wellcome_NL_df), "\n")
  cat("Genes in merged output     :", nrow(meta_df), "\n")
  cat("Genes present in both      :", n_both, "\n")
  cat("Genes with both TWAS zscores :", n_both_zscores, "\n")
  cat("Same TWAS zscore direction   :", n_same_direction, "\n")
  cat("Different TWAS zscore direction :", n_different_direction, "\n")

  cat("Finished: ", out_file, "\n\n", sep = "")

  # Return one summary row for the cross-cell-type summary table.
  data.frame(
    Cohort = cohort,
    Cell_type = CT,
    R2_global_filter = "> 0.01",
    HTAN_NL_kept = nrow(HTAN_NL_df),
    HTAN_NL_total = n_HTAN_raw,
    WellcomeSanger_NL_kept = nrow(Wellcome_NL_df),
    WellcomeSanger_NL_total = n_Wellcome_raw,
    Genes_in_HTAN_NL = nrow(HTAN_NL_df),
    Genes_in_WellcomeSanger_NL = nrow(Wellcome_NL_df),
    Genes_in_merged_output = nrow(meta_df),
    Genes_present_in_both = n_both,
    Genes_with_both_TWAS_zscores = n_both_zscores,
    Genes_same_TWAS_zscore_direction = n_same_direction,
    Genes_different_TWAS_zscore_direction = n_different_direction,
    stringsAsFactors = FALSE
  )
}

summary_list <- lapply(CTs_to_run, run_meta)
summary_df <- rbindlist(summary_list, fill = TRUE)

if (length(CTs_to_run) == length(CTlist) && setequal(CTs_to_run, CTlist)) {
  summary_name <- paste0(
    "HTAN_NL_WellcomeSanger_NL_",
    cohort,
    "_CauchyMeta_summary_7celltypes_R2global_gt0.01.csv"
  )
} else {
  summary_name <- paste0(
    "HTAN_NL_WellcomeSanger_NL_",
    cohort,
    "_CauchyMeta_summary_",
    paste(CTs_to_run, collapse = "_"),
    "_R2global_gt0.01.csv"
  )
}

summary_file <- file.path(OUT_DIR, summary_name)

fwrite(summary_df, summary_file)

cat("============================================================\n")
cat("Summary table\n")
cat("============================================================\n")
print(summary_df)
cat("Cohort: ", cohort, "\n", sep = "")
cat("Summary written to: ", summary_file, "\n", sep = "")
