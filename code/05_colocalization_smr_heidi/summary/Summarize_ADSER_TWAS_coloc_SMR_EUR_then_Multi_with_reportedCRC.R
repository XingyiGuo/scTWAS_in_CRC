#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
})

# ============================================================================
# ADSER cell-type TWAS summary + coloc/SMR integration
#
# Goal
#   1) Summarize EUR first, then Multi ancestry.
#   2) For each cell type and ancestry:
#        - keep TWAS models with R2_global > 0.01
#        - define n_TWAS_tests as the number of unique genes passing that filter
#        - define TWAS significance as pvalue < 0.05 / n_TWAS_tests
#        - collect all significant TWAS genes in a dataframe
#   3) For those TWAS-significant genes only, attach ancestry-matched coloc and
#      SMR/HEIDI results.
#   4) Match coloc_gene_name against previously reported CRC genes in
#      a_Sum_Reported_CRC_genes_202601.csv. If present, add Resources_10Res
#      to the output row. Non-matches remain NA. No output rows are removed.
#
# Rules inherited from the supplied scripts
#   - If duplicate TWAS rows exist for a gene after R2 filtering, retain the row
#     with the smallest pvalue (same behavior as the supplied Python script).
#   - SMR number of tests is the number of SMR result rows for that cell type
#     and ancestry, independently.
#   - P_SMR_BF = min(p_SMR * n_SMR_tests, 1).
#   - SMR support: P_SMR_BF < 0.05 AND p_HEIDI >= 0.05.
#   - coloc support: PP.H4.abf > 0.5 or > 0.8.
#   - coloc/SMR are left-joined onto the TWAS-significant gene backbone;
#     missing evidence remains NA.
#
# Usage
#   Rscript Summarize_ADSER_TWAS_coloc_SMR_EUR_then_Multi.R
#
# Optional custom root directory
#   Rscript Summarize_ADSER_TWAS_coloc_SMR_EUR_then_Multi.R \
#     /v5000/data/l2_bioinfo1/liq17/scTWAS/TWAS_HTAN_GSE266546_WellcomeSanger_res
# ============================================================================

# -----------------------------
# User settings
# -----------------------------
CELL_TYPES <- c("ABS", "CT", "EE", "GOB", "STM", "TAC", "TUF", "ASC", "SSC")
COHORTS <- c("EUR", "Multi")  # IMPORTANT: EUR is processed and reported first.
R2_CUTOFF <- 0.01

args <- commandArgs(trailingOnly = TRUE)

resolve_twas_dir <- function(args) {
  if (length(args) >= 1L && nzchar(args[1])) {
    return(normalizePath(args[1], winslash = "/", mustWork = FALSE))
  }

  env_dir <- Sys.getenv("TWAS_DIR", unset = "")
  if (nzchar(env_dir)) {
    return(normalizePath(env_dir, winslash = "/", mustWork = FALSE))
  }

  candidates <- c(
    "/v5000/data/l2_bioinfo1/liq17/scTWAS/TWAS_HTAN_GSE266546_WellcomeSanger_res",
    "/data/l2_bioinfo1/liq17/scTWAS/TWAS_HTAN_GSE266546_WellcomeSanger_res"
  )

  has_twas_inputs <- vapply(
    candidates,
    function(x) dir.exists(file.path(x, "SPrediXcan_peakSNP_association")),
    logical(1)
  )

  if (any(has_twas_inputs)) {
    return(candidates[which(has_twas_inputs)[1]])
  }

  # Fall back to the path used by the supplied R integration script.
  candidates[1]
}

TWAS_DIR <- resolve_twas_dir(args)

# Previously reported CRC genes used for annotation.
# The second command-line argument or SUM_REPORT_FILE environment variable
# can override this default cluster path.
resolve_sum_report_file <- function(args) {
  if (length(args) >= 2L && nzchar(args[2])) {
    return(normalizePath(args[2], winslash = "/", mustWork = FALSE))
  }

  env_file <- Sys.getenv("SUM_REPORT_FILE", unset = "")
  if (nzchar(env_file)) {
    return(normalizePath(env_file, winslash = "/", mustWork = FALSE))
  }

  "/data/l2_bioinfo1/liq17/LQ_references/pmid_reported_genes/crc_genes_raw/a_Sum_Reported_CRC_genes_202601.csv"
}

SUM_REPORT_FILE <- resolve_sum_report_file(args)

TWAS_INPUT_DIR <- file.path(TWAS_DIR, "SPrediXcan_peakSNP_association")
SMR_INPUT_DIR <- file.path(
  TWAS_DIR,
  "SMR_HEIDI_GWAS_eQTL_sigSource/HTAN_ADSER_targetTopEQTL_500kb"
)
COLOC_FILE <- file.path(
  TWAS_DIR,
  "coloc_GWAS_eQTL_sigSource/VUMC_HTAN_ADSER",
  "coloc.ADSER.9cts.GWAS_eQTL_sigSource.summary.csv"
)

OUT_DIR <- file.path(
  TWAS_DIR,
  "TWAS_summary_EUR_then_Multi_coloc_SMR"
)
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

dir.create(file.path(OUT_DIR, "EUR"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(OUT_DIR, "Multi"), recursive = TRUE, showWarnings = FALSE)

# -----------------------------
# Helpers
# -----------------------------
check_columns <- function(dt, required, label) {
  missing_cols <- setdiff(required, names(dt))
  if (length(missing_cols) > 0L) {
    stop(
      label, " is missing required column(s): ",
      paste(missing_cols, collapse = ", ")
    )
  }
}

count_present <- function(x) sum(!is.na(x))
count_true <- function(x) sum(x %in% TRUE)

prefix_except_gene <- function(dt, prefix) {
  x <- copy(dt)
  old <- setdiff(names(x), "gene")
  if (length(old) > 0L) {
    setnames(x, old, paste0(prefix, old))
  }
  x
}

# Annotate an existing TWAS+coloc+SMR dataframe using coloc_gene_name.
# Matching is case-insensitive and ignores leading/trailing whitespace.
# Existing rows are preserved; unmatched genes remain NA.
add_previous_crc_resources <- function(dt, sum_report) {
  x <- copy(dt)

  # Always create the requested column.
  x[, Resources_10Res := NA_character_]

  if (!"coloc_gene_name" %in% names(x) || nrow(x) == 0L) {
    return(x)
  }

  lookup <- copy(sum_report[, .(Gene, Resources_10Res)])
  lookup[, .gene_key := toupper(trimws(as.character(Gene)))]
  lookup <- lookup[!is.na(.gene_key) & nzchar(.gene_key)]

  if (anyDuplicated(lookup$.gene_key)) {
    stop(
      "SUM_REPORT_FILE contains duplicated Gene values after case-insensitive ",
      "normalization. Resolve duplicates before annotation."
    )
  }

  resource_map <- setNames(
    as.character(lookup$Resources_10Res),
    lookup$.gene_key
  )

  coloc_key <- toupper(trimws(as.character(x$coloc_gene_name)))
  x[, Resources_10Res := unname(resource_map[coloc_key])]

  x
}

# Python-equivalent TWAS filtering:
#   - R2_global > 0.01
#   - if gene duplicates exist, keep the smallest pvalue row
read_twas <- function(ct, cohort) {
  f <- file.path(
    TWAS_INPUT_DIR,
    paste0("VUMC_HTAN_ADSER_", ct, "_peakSNP.", cohort, ".TWAS.csv")
  )

  if (!file.exists(f)) {
    stop("Missing TWAS file: ", f)
  }

  x <- fread(f)
  check_columns(x, c("gene", "pvalue", "R2_global"), paste0("TWAS [", ct, "/", cohort, "]"))

  n_raw_rows <- nrow(x)

  x <- x[!is.na(R2_global) & R2_global > R2_CUTOFF]
  n_R2_rows_before_dedup <- nrow(x)

  n_duplicate_rows_removed <- 0L
  if (anyDuplicated(x$gene)) {
    # Same intent as Python sort_values("pvalue").drop_duplicates("gene", keep="first").
    x <- x[order(pvalue, na.last = TRUE)]
    n_before <- nrow(x)
    x <- x[!duplicated(gene)]
    n_duplicate_rows_removed <- n_before - nrow(x)
  }

  n_tests <- nrow(x)
  bf_cutoff <- if (n_tests > 0L) 0.05 / n_tests else NA_real_

  x[, TWAS_sig := !is.na(pvalue) & pvalue < bf_cutoff]

  sig <- copy(x[TWAS_sig == TRUE])
  if (nrow(sig) > 0L) {
    setorder(sig, pvalue)
  }

  list(
    all_tested = x,
    significant = sig,
    n_raw_rows = n_raw_rows,
    n_R2_rows_before_dedup = n_R2_rows_before_dedup,
    n_duplicate_rows_removed = n_duplicate_rows_removed,
    n_tests = n_tests,
    bf_cutoff = bf_cutoff,
    n_sig = nrow(sig),
    file = f
  )
}

prepare_coloc <- function(coloc_all, ct, cohort) {
  # Use local scalar names to avoid data.table column/argument name ambiguity.
  ct_value <- ct
  cohort_value <- cohort
  coloc <- copy(coloc_all[cell_type == ct_value & cohort == cohort_value])

  if (nrow(coloc) == 0L) {
    return(data.table(gene = character()))
  }

  if (anyDuplicated(coloc$gene)) {
    stop(
      "coloc results contain duplicated gene rows for ", ct, "/", cohort,
      ". Resolve duplicates before joining."
    )
  }

  coloc[, support_PPH4_gt_0p5 := !is.na(PP.H4.abf) & PP.H4.abf > 0.5]
  coloc[, support_PPH4_gt_0p8 := !is.na(PP.H4.abf) & PP.H4.abf > 0.8]

  keep <- c(
    "gene", "gene_name", "chromosome", "TWAS_p", "qtl_N",
    "n_eqtl_variants", "n_overlap_harmonized", "nsnps",
    "PP.H0.abf", "PP.H1.abf", "PP.H2.abf", "PP.H3.abf", "PP.H4.abf",
    "status", "support_PPH4_gt_0p5", "support_PPH4_gt_0p8"
  )
  keep <- intersect(keep, names(coloc))

  prefix_except_gene(coloc[, ..keep], "coloc_")
}

prepare_smr <- function(ct, cohort) {
  f <- file.path(
    SMR_INPUT_DIR,
    paste0("SMR_HEIDI.ADSER.", ct, ".SMR_HEIDI.results.tsv")
  )

  if (!file.exists(f)) {
    stop("Missing SMR/HEIDI file: ", f)
  }

  smr <- fread(f)
  check_columns(
    smr,
    c("cell_type", "gene", "cohort", "p_SMR", "p_HEIDI"),
    paste0("SMR [", ct, "/", cohort, "]")
  )

  ct_value <- ct
  cohort_value <- cohort
  smr <- smr[cell_type == ct_value & cohort == cohort_value]

  if (anyDuplicated(smr$gene)) {
    stop(
      "SMR results contain duplicated gene rows for ", ct, "/", cohort,
      ". Resolve duplicates before joining."
    )
  }

  n_SMR_tests <- nrow(smr)
  SMR_BF_cutoff <- if (n_SMR_tests > 0L) 0.05 / n_SMR_tests else NA_real_

  if (n_SMR_tests > 0L) {
    n_smr_tests_value <- n_SMR_tests
    smr_bf_cutoff_value <- SMR_BF_cutoff
    smr[, n_SMR_tests := n_smr_tests_value]
    smr[, P_SMR_BF := pmin(p_SMR * n_smr_tests_value, 1)]
    smr[, SMR_BF_cutoff := smr_bf_cutoff_value]
    smr[, SMR_support :=
          !is.na(P_SMR_BF) & P_SMR_BF < 0.05 &
          !is.na(p_HEIDI) & p_HEIDI >= 0.05]
  }

  keep <- c(
    "gene", "gene_name", "target_SNP", "run_status", "n_ld_snps",
    "TWAS_p", "target_pos",
    "target_eQTL_p", "target_eQTL_beta", "target_eQTL_se", "target_eQTL_freq",
    "target_GWAS_p", "target_GWAS_beta", "target_GWAS_se", "target_GWAS_freq",
    "b_GWAS", "se_GWAS", "p_GWAS",
    "b_eQTL", "se_eQTL", "p_eQTL",
    "b_SMR", "se_SMR", "p_SMR",
    "P_SMR_BF", "SMR_BF_cutoff", "p_HEIDI", "nsnp_HEIDI",
    "n_SMR_tests", "SMR_support"
  )
  keep <- intersect(keep, names(smr))

  smr_out <- if (nrow(smr) > 0L) {
    prefix_except_gene(smr[, ..keep], "SMR_")
  } else {
    data.table(gene = character())
  }

  list(
    data = smr_out,
    n_SMR_tests = n_SMR_tests,
    SMR_BF_cutoff = SMR_BF_cutoff,
    file = f
  )
}

attach_evidence <- function(sig_twas, coloc_dt, smr_dt) {
  x <- copy(sig_twas)
  x[, .TWAS_row_order := .I]

  x <- merge(x, coloc_dt, by = "gene", all.x = TRUE, sort = FALSE)
  x <- merge(x, smr_dt, by = "gene", all.x = TRUE, sort = FALSE)

  setorder(x, .TWAS_row_order)
  x[, .TWAS_row_order := NULL]

  coloc05 <- if ("coloc_support_PPH4_gt_0p5" %in% names(x)) {
    x[["coloc_support_PPH4_gt_0p5"]]
  } else {
    rep(NA, nrow(x))
  }

  coloc08 <- if ("coloc_support_PPH4_gt_0p8" %in% names(x)) {
    x[["coloc_support_PPH4_gt_0p8"]]
  } else {
    rep(NA, nrow(x))
  }

  smr_support <- if ("SMR_SMR_support" %in% names(x)) {
    x[["SMR_SMR_support"]]
  } else {
    rep(NA, nrow(x))
  }

  # Preserve three-state logic for missing evidence, as in the supplied R script.
  x[, support_coloc0p5_or_SMR := coloc05 | smr_support]
  x[, support_coloc0p8_or_SMR := coloc08 | smr_support]
  x[, support_coloc0p5_and_SMR := coloc05 & smr_support]
  x[, support_coloc0p8_and_SMR := coloc08 & smr_support]

  x
}

summarize_evidence <- function(integrated, ct, cohort, twas_info, smr_info) {
  get_or_na <- function(column, mode = "logical") {
    if (column %in% names(integrated)) return(integrated[[column]])
    if (mode == "numeric") return(rep(NA_real_, nrow(integrated)))
    rep(NA, nrow(integrated))
  }

  pph4 <- get_or_na("coloc_PP.H4.abf", "numeric")
  coloc05 <- get_or_na("coloc_support_PPH4_gt_0p5")
  coloc08 <- get_or_na("coloc_support_PPH4_gt_0p8")
  smr_p <- get_or_na("SMR_p_SMR", "numeric")
  smr_support <- get_or_na("SMR_SMR_support")
  either05 <- get_or_na("support_coloc0p5_or_SMR")
  either08 <- get_or_na("support_coloc0p8_or_SMR")
  both05 <- get_or_na("support_coloc0p5_and_SMR")
  both08 <- get_or_na("support_coloc0p8_and_SMR")

  data.table(
    ancestry = cohort,
    cell_type = ct,
    R2_global_filter = paste0("> ", R2_CUTOFF),
    n_TWAS_raw_rows = twas_info$n_raw_rows,
    n_TWAS_R2_rows_before_dedup = twas_info$n_R2_rows_before_dedup,
    n_TWAS_duplicate_rows_removed = twas_info$n_duplicate_rows_removed,
    n_TWAS_tests = twas_info$n_tests,
    TWAS_Bonferroni_cutoff = twas_info$bf_cutoff,
    n_TWAS_sig = twas_info$n_sig,
    n_SMR_tests = smr_info$n_SMR_tests,
    SMR_Bonferroni_cutoff = smr_info$SMR_BF_cutoff,
    n_sig_with_coloc_result = count_present(pph4),
    n_sig_coloc_PPH4_gt_0p5 = count_true(coloc05),
    n_sig_coloc_PPH4_gt_0p8 = count_true(coloc08),
    n_sig_with_SMR_result = count_present(smr_p),
    n_sig_SMR_support = count_true(smr_support),
    n_sig_either_coloc0p5_or_SMR = count_true(either05),
    n_sig_either_coloc0p8_or_SMR = count_true(either08),
    n_sig_both_coloc0p5_and_SMR = count_true(both05),
    n_sig_both_coloc0p8_and_SMR = count_true(both08)
  )
}

# -----------------------------
# Read shared coloc once
# -----------------------------
if (!file.exists(COLOC_FILE)) {
  stop("Missing coloc input file: ", COLOC_FILE)
}

coloc_all <- fread(COLOC_FILE)
check_columns(
  coloc_all,
  c("cell_type", "gene", "cohort", "PP.H4.abf"),
  "coloc"
)

if (!file.exists(SUM_REPORT_FILE)) {
  stop(
    "Missing previous reported CRC-gene summary file: ", SUM_REPORT_FILE, "\n",
    "Provide it as the second command-line argument or set SUM_REPORT_FILE."
  )
}

Sum_report <- fread(SUM_REPORT_FILE)
check_columns(
  Sum_report,
  c("Gene", "Resources_10Res"),
  "Previous reported CRC-gene summary"
)

cat("Previous reported CRC-gene file: ", SUM_REPORT_FILE, "\n", sep = "")
cat("Previously reported genes loaded: ", nrow(Sum_report), "\n", sep = "")

# -----------------------------
# Run: EUR first, then Multi
# -----------------------------
all_twas_summary <- list()
all_sig_twas <- list()
all_integrated <- list()
all_evidence_summary <- list()

for (cohort in COHORTS) {
  cat("\n", strrep("=", 88), "\n", sep = "")
  cat("ANCESTRY: ", cohort, "\n", sep = "")
  cat(strrep("=", 88), "\n", sep = "")

  for (ct in CELL_TYPES) {
    cat("\n--- ", cohort, " | ", ct, " ---\n", sep = "")

    twas_info <- read_twas(ct, cohort)

    # Add analysis metadata to the significant-gene dataframe.
    sig_twas <- copy(twas_info$significant)
    sig_twas[, ancestry := cohort]
    sig_twas[, cell_type := ct]
    sig_twas[, n_TWAS_tests := twas_info$n_tests]
    sig_twas[, TWAS_Bonferroni_cutoff := twas_info$bf_cutoff]

    # For Multi-ancestry TWAS-significant genes, indicate whether the same gene
    # is also TWAS-significant in EUR for the SAME cell type.
    #
    # Because COHORTS = c("EUR", "Multi"), EUR results have already been
    # processed before the corresponding Multi results are reached here.
    if (cohort == "Multi") {
      eur_key <- paste("EUR", ct, sep = "__")

      if (!eur_key %in% names(all_sig_twas)) {
        stop(
          "EUR significant-gene results were not available before Multi for cell type: ",
          ct
        )
      }

      eur_sig_genes <- unique(all_sig_twas[[eur_key]]$gene)

      sig_twas[, Sig_resource := fifelse(
        gene %chin% eur_sig_genes,
        "EUR_sig",
        "EUR_notSig_Multi_Sig"
      )]
    }

    # Put key identifiers first.
    front <- c(
      "ancestry", "cell_type", "gene", "Sig_resource",
      "pvalue", "R2_global",
      "n_TWAS_tests", "TWAS_Bonferroni_cutoff", "TWAS_sig"
    )
    front <- intersect(front, names(sig_twas))
    setcolorder(sig_twas, c(front, setdiff(names(sig_twas), front)))

    coloc_dt <- prepare_coloc(coloc_all, ct, cohort)
    smr_info <- prepare_smr(ct, cohort)

    integrated <- attach_evidence(sig_twas, coloc_dt, smr_info$data)

    # Annotate previously reported CRC genes using coloc_gene_name.
    integrated <- add_previous_crc_resources(integrated, Sum_report)

    # Put Resources_10Res immediately after coloc_gene_name for readability.
    if ("coloc_gene_name" %in% names(integrated) &&
        "Resources_10Res" %in% names(integrated)) {
      cols <- names(integrated)
      cols <- cols[cols != "Resources_10Res"]
      pos <- match("coloc_gene_name", cols)
      cols <- append(cols, "Resources_10Res", after = pos)
      setcolorder(integrated, cols)
    }

    twas_summary <- data.table(
      ancestry = cohort,
      cell_type = ct,
      R2_global_filter = paste0("> ", R2_CUTOFF),
      n_raw_rows = twas_info$n_raw_rows,
      n_R2_rows_before_dedup = twas_info$n_R2_rows_before_dedup,
      n_duplicate_rows_removed = twas_info$n_duplicate_rows_removed,
      n_TWAS_tests = twas_info$n_tests,
      TWAS_Bonferroni_cutoff = twas_info$bf_cutoff,
      n_significant_genes = twas_info$n_sig
    )

    evidence_summary <- summarize_evidence(
      integrated = integrated,
      ct = ct,
      cohort = cohort,
      twas_info = twas_info,
      smr_info = smr_info
    )

    key <- paste(cohort, ct, sep = "__")
    all_twas_summary[[key]] <- twas_summary
    all_sig_twas[[key]] <- sig_twas
    all_integrated[[key]] <- integrated
    all_evidence_summary[[key]] <- evidence_summary

    # Per-cell-type outputs, separated by ancestry.
    out_subdir <- file.path(OUT_DIR, cohort)

    fwrite(
      sig_twas,
      file.path(
        out_subdir,
        paste0("VUMC_HTAN_ADSER_", ct, ".", cohort, ".TWAS_significant.csv")
      ),
      na = ""
    )

    fwrite(
      integrated,
      file.path(
        out_subdir,
        paste0(
          "VUMC_HTAN_ADSER_", ct, ".", cohort,
          ".TWAS_significant.coloc_SMR.csv"
        )
      ),
      na = ""
    )

    cat("TWAS tests (R2_global > 0.01): ", twas_info$n_tests, "\n", sep = "")
    cat("TWAS Bonferroni cutoff       : ", signif(twas_info$bf_cutoff, 6), "\n", sep = "")
    cat("TWAS significant genes       : ", twas_info$n_sig, "\n", sep = "")
    cat("SMR tests                    : ", smr_info$n_SMR_tests, "\n", sep = "")
    cat("coloc PP.H4 > 0.5 among sig  : ", evidence_summary$n_sig_coloc_PPH4_gt_0p5, "\n", sep = "")
    cat("SMR-supported among sig      : ", evidence_summary$n_sig_SMR_support, "\n", sep = "")
  }
}

# -----------------------------
# Combined dataframes / outputs
# -----------------------------
TWAS_summary_df <- rbindlist(all_twas_summary, use.names = TRUE, fill = TRUE)
Sig_TWAS_df <- rbindlist(all_sig_twas, use.names = TRUE, fill = TRUE)
Sig_TWAS_coloc_SMR_df <- rbindlist(all_integrated, use.names = TRUE, fill = TRUE)
Evidence_summary_df <- rbindlist(all_evidence_summary, use.names = TRUE, fill = TRUE)

# -------------------------------------------------------------------------
# FAIL-SAFE re-annotation of the FINAL combined TWAS+coloc+SMR dataframe.
#
# This guarantees Resources_10Res is present in:
#   ADSER_EUR_TWAS_significant_genes_coloc_SMR.csv
#   ADSER_Multi_TWAS_significant_genes_coloc_SMR.csv
#   ADSER_TWAS_significant_genes_coloc_SMR_all.csv
#
# Matching: coloc_gene_name -> Sum_report$Gene
# Existing rows are preserved; unmatched genes remain NA.
# -------------------------------------------------------------------------
Sig_TWAS_coloc_SMR_df <- add_previous_crc_resources(
  Sig_TWAS_coloc_SMR_df,
  Sum_report
)

if (!"Resources_10Res" %in% names(Sig_TWAS_coloc_SMR_df)) {
  stop("Internal error: Resources_10Res was not added to final ADSER integrated dataframe.")
}

if ("coloc_gene_name" %in% names(Sig_TWAS_coloc_SMR_df)) {
  cols <- names(Sig_TWAS_coloc_SMR_df)
  cols <- cols[cols != "Resources_10Res"]
  pos <- match("coloc_gene_name", cols)
  cols <- append(cols, "Resources_10Res", after = pos)
  setcolorder(Sig_TWAS_coloc_SMR_df, cols)
}

n_reported_matches <- Sig_TWAS_coloc_SMR_df[
  !is.na(Resources_10Res) & nzchar(trimws(Resources_10Res)),
  .N
]

cat(
  "\nPreviously reported CRC-gene matches in final ADSER TWAS+coloc+SMR table: ",
  n_reported_matches, " / ", nrow(Sig_TWAS_coloc_SMR_df), "\n",
  sep = ""
)

# Force requested report order: EUR first, then Multi; cell types in CELL_TYPES order.
TWAS_summary_df[, .ancestry_order := match(ancestry, COHORTS)]
TWAS_summary_df[, .cell_order := match(cell_type, CELL_TYPES)]
setorder(TWAS_summary_df, .ancestry_order, .cell_order)
TWAS_summary_df[, c(".ancestry_order", ".cell_order") := NULL]

Evidence_summary_df[, .ancestry_order := match(ancestry, COHORTS)]
Evidence_summary_df[, .cell_order := match(cell_type, CELL_TYPES)]
setorder(Evidence_summary_df, .ancestry_order, .cell_order)
Evidence_summary_df[, c(".ancestry_order", ".cell_order") := NULL]

# The significant-gene tables were already accumulated EUR then Multi.
if (nrow(Sig_TWAS_df) > 0L) {
  Sig_TWAS_df[, .ancestry_order := match(ancestry, COHORTS)]
  Sig_TWAS_df[, .cell_order := match(cell_type, CELL_TYPES)]
  setorder(Sig_TWAS_df, .ancestry_order, .cell_order, pvalue)
  Sig_TWAS_df[, c(".ancestry_order", ".cell_order") := NULL]
}

if (nrow(Sig_TWAS_coloc_SMR_df) > 0L) {
  Sig_TWAS_coloc_SMR_df[, .ancestry_order := match(ancestry, COHORTS)]
  Sig_TWAS_coloc_SMR_df[, .cell_order := match(cell_type, CELL_TYPES)]
  setorder(Sig_TWAS_coloc_SMR_df, .ancestry_order, .cell_order, pvalue)
  Sig_TWAS_coloc_SMR_df[, c(".ancestry_order", ".cell_order") := NULL]
}

# Final check before writing combined manuscript/reporting CSVs.
stopifnot("Resources_10Res" %in% names(Sig_TWAS_coloc_SMR_df))

cat(
  "Resources_10Res column confirmed for ADSER EUR, Multi, and all-ancestry ",
  "TWAS+coloc+SMR CSV outputs.\n"
)

# Main combined outputs.
fwrite(
  TWAS_summary_df,
  file.path(OUT_DIR, "ADSER_TWAS_summary_by_ancestry_celltype.csv"),
  na = ""
)

fwrite(
  Sig_TWAS_df,
  file.path(OUT_DIR, "ADSER_TWAS_significant_genes_all.csv"),
  na = ""
)

fwrite(
  Sig_TWAS_coloc_SMR_df,
  file.path(OUT_DIR, "ADSER_TWAS_significant_genes_coloc_SMR_all.csv"),
  na = ""
)

fwrite(
  Evidence_summary_df,
  file.path(OUT_DIR, "ADSER_TWAS_significant_genes_coloc_SMR_summary.csv"),
  na = ""
)

# Separate ancestry-level combined outputs for convenient manuscript summaries.
for (cohort in COHORTS) {
  fwrite(
    TWAS_summary_df[ancestry == cohort],
    file.path(OUT_DIR, paste0("ADSER_", cohort, "_TWAS_summary_by_celltype.csv")),
    na = ""
  )

  fwrite(
    Sig_TWAS_df[ancestry == cohort],
    file.path(OUT_DIR, paste0("ADSER_", cohort, "_TWAS_significant_genes.csv")),
    na = ""
  )

  fwrite(
    Sig_TWAS_coloc_SMR_df[ancestry == cohort],
    file.path(
      OUT_DIR,
      paste0("ADSER_", cohort, "_TWAS_significant_genes_coloc_SMR.csv")
    ),
    na = ""
  )

  fwrite(
    Evidence_summary_df[ancestry == cohort],
    file.path(
      OUT_DIR,
      paste0("ADSER_", cohort, "_TWAS_significant_genes_coloc_SMR_summary.csv")
    ),
    na = ""
  )
}

cat("\n", strrep("=", 88), "\n", sep = "")
cat("FINISHED\n")
cat(strrep("=", 88), "\n", sep = "")
cat("Output directory: ", OUT_DIR, "\n", sep = "")

cat("\nTWAS summary (EUR first, then Multi):\n")
print(TWAS_summary_df)

cat("\nColoc/SMR summary among TWAS-significant genes:\n")
print(Evidence_summary_df)
