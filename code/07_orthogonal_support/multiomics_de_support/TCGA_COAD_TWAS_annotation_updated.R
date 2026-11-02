#!/usr/bin/env Rscript

# =============================================================================
# TWAS + TCGA-COAD RNA-seq DE annotation
# Input : TWAS_ADSER_NL_collected_results.csv
# Output: TWAS_ADSER_NL_TCGA_COAD_DE_support.csv
# =============================================================================

suppressPackageStartupMessages({
  library(TCGAbiolinks)
  library(SummarizedExperiment)
  library(DESeq2)
  library(dplyr)
  library(readr)
  library(tibble)
})

# =============================================================================
# 1. Paths
# =============================================================================

BASE_DIR <- paste0(
  "/data/l2_bioinfo1/liq17/scTWAS/",
  "TWAS_HTAN_GSE266546_WellcomeSanger_res/",
  "SPrediXcan_peakSNP_association"
)

TWAS_FILE <- file.path(
  BASE_DIR,
  "TWAS_ADSER_NL_collected_results.csv"
)

OUT_DIR <- file.path(
  "/data/l2_bioinfo1/liq17/scTWAS/Functional_res/TCGA_CPTAC_COAD_DE"
)

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# If your GDC files were downloaded somewhere else, change this path.
TCGA_GDC_DIR <- file.path( "/data/l2_bioinfo1/liq17/LQ_references/TCGA_CPTAC_DepMap/GDC_TCGA_COAD")

dir.create(TCGA_GDC_DIR, recursive = TRUE, showWarnings = FALSE)

TCGA_SE_CACHE <- file.path(
  OUT_DIR,
  "TCGA_COAD_STARcounts_SE.rds"
)

TCGA_DE_CACHE <- file.path(
  OUT_DIR,
  "TCGA_COAD_DE_allgenes.rds"
)

TCGA_DE_CSV <- file.path(
  OUT_DIR,
  "TCGA_COAD_DE_allgenes.csv"
)

FINAL_OUT <- file.path(
  OUT_DIR,
  "TWAS_ADSER_NL_TCGA_COAD_DE_support.csv"
)

LOG_FILE <- file.path(
  OUT_DIR,
  "TCGA_COAD_TWAS_annotation.log"
)

# =============================================================================
# Logging and software versions
# =============================================================================

log_con <- file(LOG_FILE, open = "wt")
sink(log_con, type = "output", split = TRUE)
sink(log_con, type = "message")

cat("============================================================\n")
cat("TCGA-COAD TWAS DE SUPPORT ANALYSIS\n")
cat("============================================================\n")
cat("Run started:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), "\n")
cat("R version:", R.version.string, "\n")
cat("Platform:", R.version$platform, "\n")
cat("\nPackage versions:\n")

packages_to_log <- c(
  "TCGAbiolinks",
  "SummarizedExperiment",
  "DESeq2",
  "dplyr",
  "readr",
  "tibble"
)

for (pkg in packages_to_log) {
  cat(sprintf("  %-22s %s\n", pkg, as.character(packageVersion(pkg))))
}

cat("\nDE support thresholds:\n")
cat("  FDR < 0.05\n")
cat("  |log2FC| > 0.5\n")
cat("  sign(log2FC) == sign(scTWAS Z score)\n")
cat("\nSupport labels:\n")
cat("  support        : passes FDR/effect-size thresholds and direction agrees\n")
cat("  DE_discordant  : passes FDR/effect-size thresholds but direction disagrees\n")
cat("  NS             : tested but does not pass FDR/effect-size thresholds\n")
cat("  not_tested     : DE result unavailable\n")
cat("============================================================\n")

# =============================================================================
# 2. Read TWAS results
# =============================================================================

twas <- read_csv(TWAS_FILE, show_col_types = FALSE)

required_cols <- c(
  "Dataset", "CellType", "gene", "gene_name",
  "Sig_resource", "Z_score", "P_value"
)

missing_cols <- setdiff(required_cols, colnames(twas))
if (length(missing_cols) > 0) {
  stop(
    paste(
      "Missing required TWAS columns:",
      paste(missing_cols, collapse = ", ")
    )
  )
}

twas <- twas %>%
  mutate(
    gene = sub("\\..*$", "", as.character(gene)),
    gene_name = as.character(gene_name),
    Z_score = as.numeric(Z_score),
    P_value = as.numeric(P_value)
  )

cat("\n============================================================\n")
cat("TWAS INPUT\n")
cat("============================================================\n")
cat("Rows:", nrow(twas), "\n")
cat("Unique Ensembl genes:", n_distinct(twas$gene), "\n")
cat("Unique gene symbols:", n_distinct(twas$gene_name), "\n")

# =============================================================================
# 3. DE support
# =============================================================================

get_DE_support <- function(
  Z_score,
  log2FC,
  FDR,
  fdr_cutoff = 0.05,
  abs_log2fc_cutoff = 0.5
) {
  case_when(
    # No usable DE result for this gene.
    is.na(log2FC) | is.na(FDR) | is.na(Z_score) ~ "not_tested",

    # A gene is considered DE only when BOTH thresholds are passed.
    FDR >= fdr_cutoff | abs(log2FC) <= abs_log2fc_cutoff ~ "NS",

    # scTWAS Z scores are expected to be directional; guard against exact zero.
    Z_score == 0 ~ "NS",

    # Significant DE with sufficient effect size and concordant direction.
    sign(Z_score) == sign(log2FC) ~ "support",

    # Significant DE with sufficient effect size but discordant direction.
    TRUE ~ "DE_discordant"
  )
}

# =============================================================================
# 4. Load/query TCGA-COAD
# =============================================================================

if (file.exists(TCGA_SE_CACHE)) {

  cat("\nLoading cached TCGA-COAD SummarizedExperiment...\n")
  se_coad <- readRDS(TCGA_SE_CACHE)

} else {

  cat("\nQuerying GDC for TCGA-COAD metadata...\n")

  query_coad <- GDCquery(
    project = "TCGA-COAD",
    data.category = "Transcriptome Profiling",
    data.type = "Gene Expression Quantification",
    workflow.type = "STAR - Counts",
    sample.type = c("Primary Tumor", "Solid Tissue Normal")
  )

  # First try to prepare from files already present locally.
  se_coad <- tryCatch(
    GDCprepare(
      query = query_coad,
      directory = TCGA_GDC_DIR,
      summarizedExperiment = TRUE
    ),
    error = function(e) {
      message("Existing local files could not be prepared: ", conditionMessage(e))
      NULL
    }
  )

  # Only download if local files are missing/incomplete.
  if (is.null(se_coad)) {

    cat("\nTCGA local files missing/incomplete. Downloading...\n")

    GDCdownload(
      query = query_coad,
      directory = TCGA_GDC_DIR,
      method = "api",
      files.per.chunk = 100
    )

    se_coad <- GDCprepare(
      query = query_coad,
      directory = TCGA_GDC_DIR,
      summarizedExperiment = TRUE
    )
  }

  saveRDS(se_coad, TCGA_SE_CACHE)
}

# =============================================================================
# 5. TCGA differential expression
# =============================================================================

if (file.exists(TCGA_DE_CACHE)) {

  cat("\nLoading cached TCGA DE results...\n")
  tcga_de <- readRDS(TCGA_DE_CACHE)

} else {

  cat("\nTCGA assay names:\n")
  print(assayNames(se_coad))

  if ("unstranded" %in% assayNames(se_coad)) {
    count_assay <- "unstranded"
  } else {
    count_candidates <- grep(
      "unstranded",
      assayNames(se_coad),
      value = TRUE,
      ignore.case = TRUE
    )
    if (length(count_candidates) == 0) {
      stop(
        paste(
          "Cannot identify raw STAR count assay. Available assays:",
          paste(assayNames(se_coad), collapse = ", ")
        )
      )
    }
    count_assay <- count_candidates[1]
  }

  cat("\nUsing count assay:", count_assay, "\n")

  tcga_counts <- assay(se_coad, count_assay)

  # Clean Ensembl versions and collapse duplicates if needed.
  tcga_gene_ids <- sub("\\..*$", "", rownames(tcga_counts))

  if (anyDuplicated(tcga_gene_ids) > 0) {
    cat("\nCollapsing duplicated Ensembl IDs by summing counts.\n")
    tcga_counts <- rowsum(
      tcga_counts,
      group = tcga_gene_ids,
      reorder = FALSE
    )
  } else {
    rownames(tcga_counts) <- tcga_gene_ids
  }

  # Sample type from TCGA barcode:
  # 01 = Primary Tumor, 11 = Solid Tissue Normal
  sample_ids <- colnames(tcga_counts)

  tcga_meta <- tibble(
    sample = sample_ids,
    patient = substr(sample_ids, 1, 12),
    sample_code = substr(sample_ids, 14, 15)
  ) %>%
    mutate(
      condition = case_when(
        sample_code == "01" ~ "Tumor",
        sample_code == "11" ~ "Normal",
        TRUE ~ NA_character_
      )
    ) %>%
    filter(!is.na(condition))

  cat("\nTCGA samples before deduplication:\n")
  print(table(tcga_meta$condition))

  # Keep one sample per patient per condition.
  tcga_meta <- tcga_meta %>%
    group_by(patient, condition) %>%
    slice(1) %>%
    ungroup() %>%
    mutate(
      condition = factor(
        condition,
        levels = c("Normal", "Tumor")
      )
    )

  tcga_counts <- tcga_counts[
    , tcga_meta$sample, drop = FALSE
  ]

  rownames(tcga_meta) <- tcga_meta$sample

  cat("\nTCGA samples used in DESeq2:\n")
  print(table(tcga_meta$condition))

  dds <- DESeqDataSetFromMatrix(
    countData = round(tcga_counts),
    colData = as.data.frame(tcga_meta),
    design = ~ condition
  )

  # Keep genes with >=10 counts in at least 5 samples.
  keep <- rowSums(counts(dds) >= 10) >= 5
  dds <- dds[keep, ]

  cat("\nGenes entering DESeq2:", nrow(dds), "\n")

  dds <- DESeq(dds, quiet = FALSE)

  # Tumor vs Normal. No apeglm is needed here.
  tcga_res <- results(
    dds,
    contrast = c("condition", "Tumor", "Normal"),
    alpha = 0.05
  )

  tcga_de <- as.data.frame(tcga_res) %>%
    rownames_to_column("gene") %>%
    transmute(
      gene = sub("\\..*$", "", gene),
      TCGA_log2FC = as.numeric(log2FoldChange),
      TCGA_P = as.numeric(pvalue),
      TCGA_FDR = as.numeric(padj)
    )

  saveRDS(tcga_de, TCGA_DE_CACHE)
}

write_csv(tcga_de, TCGA_DE_CSV)
cat("\nTCGA DE genes:", nrow(tcga_de), "\n")

# =============================================================================
# 6. Annotate TWAS
# =============================================================================

final <- twas %>%
  left_join(
    tcga_de %>%
      select(gene, TCGA_log2FC, TCGA_P, TCGA_FDR),
    by = "gene"
  ) %>%
  mutate(
    TCGA_DE_support = get_DE_support(
      Z_score = Z_score,
      log2FC = TCGA_log2FC,
      FDR = TCGA_FDR
    )
  ) %>%
  select(
    Dataset,
    CellType,
    gene,
    gene_name,
    Sig_resource,
    Z_score,
    P_value,
    TCGA_log2FC,
    TCGA_P,
    TCGA_FDR,
    TCGA_DE_support
  ) %>%
  arrange(Dataset, CellType, P_value)

# =============================================================================
# 7. QC + save
# =============================================================================

expected_cols <- c(
  "Dataset", "CellType", "gene", "gene_name", "Sig_resource",
  "Z_score", "P_value", "TCGA_log2FC", "TCGA_P", "TCGA_FDR",
  "TCGA_DE_support"
)

if (!identical(colnames(final), expected_cols)) {
  stop("Final TCGA column structure does not match expected output.")
}

cat("\n============================================================\n")
cat("TCGA DE SUPPORT SUMMARY\n")
cat("============================================================\n")
print(table(final$TCGA_DE_support, useNA = "ifany"))

cat(
  "\nTWAS rows with TCGA DE data:",
  sum(!is.na(final$TCGA_log2FC)),
  "/", nrow(final), "\n"
)

cat(
  "TWAS rows supported by TCGA:",
  sum(final$TCGA_DE_support == "support", na.rm = TRUE),
  "\n"
)

write_csv(final, FINAL_OUT)

cat("\n============================================================\n")
cat("TCGA ANNOTATED FILE\n")
cat("============================================================\n")
cat(FINAL_OUT, "\n")

cat("\nTCGA support labels:\n")
print(table(final$TCGA_DE_support, useNA = "ifany"))

cat("\nRun completed:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), "\n")
cat("Log file:", LOG_FILE, "\n")

sink(type = "message")
sink(type = "output")
close(log_con)
