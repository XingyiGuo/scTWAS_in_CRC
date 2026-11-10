#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
})

# ==============================================================================
# Summarize coloc results across 7 NL cell types
#
# Usage:
#   Rscript 01a_coloc_NL_7cts_summarize.R config/coloc_peakSNP_NL.config.txt
# ==============================================================================

args <- commandArgs(trailingOnly = TRUE)

CONFIG <- if (length(args) >= 1) args[1] else "coloc_peakSNP_NL.config.txt"
if (!length(args) && !file.exists(CONFIG)) {
  script_arg <- grep("^--file=", commandArgs(), value = TRUE)
  script_dir <- dirname(sub("^--file=", "", script_arg[1]))
  CONFIG <- file.path(script_dir, "config", CONFIG)
}

if (!file.exists(CONFIG)) {
  stop("Config file not found: ", CONFIG)
}

source(CONFIG)

if (!exists("COLOC_CONFIG")) {
  stop("COLOC_CONFIG was not defined after sourcing: ", CONFIG)
}

CTs <- COLOC_CONFIG$cell_types
ROOT <- COLOC_CONFIG$output_dir
GROUP <- COLOC_CONFIG$analysis_group

PPH4_CUTOFF <- if (!is.null(COLOC_CONFIG$summary_pph4_threshold)) {
  COLOC_CONFIG$summary_pph4_threshold
} else {
  0.8
}

message("Config: ", CONFIG)
message("Analysis group: ", GROUP)
message("Cell types: ", paste(CTs, collapse = ", "))
message("Output root: ", ROOT)
message("PP.H4 cutoff: ", PPH4_CUTOFF)

# ==============================================================================
# Read and combine coloc summary files
# ==============================================================================

all_coloc <- rbindlist(
  lapply(CTs, function(CT) {

    infile <- file.path(
      ROOT,
      CT,
      sprintf(
        "coloc.%s.%s.GWAS_eQTL_sigSource.summary.csv",
        GROUP,
        CT
      )
    )

    if (!file.exists(infile)) {
      warning("Missing file: ", infile)
      return(NULL)
    }

    message("Reading: ", infile)

    x <- fread(infile)

    required_cols <- c("gene", "cohort", "PP.H4.abf")
    missing_cols <- setdiff(required_cols, names(x))

    if (length(missing_cols) > 0) {
      stop(
        "Missing required columns in ", infile, ": ",
        paste(missing_cols, collapse = ", ")
      )
    }

    x[, cell_type := CT]
    x
  }),
  use.names = TRUE,
  fill = TRUE
)

if (nrow(all_coloc) == 0) {
  stop("No coloc summary files were loaded.")
}

# Keep EUR and Multi only
all_coloc <- all_coloc[cohort %in% c("EUR", "Multi")]

if (nrow(all_coloc) == 0) {
  stop("No EUR or Multi rows were found in the loaded coloc summaries.")
}

# ==============================================================================
# Summary by cell type and cohort
# ==============================================================================

summary_by_ct <- all_coloc[
  ,
  .(
    n_genes_total = uniqueN(gene),
    n_genes_PPH4_gt_cutoff = uniqueN(
      gene[!is.na(PP.H4.abf) & PP.H4.abf > PPH4_CUTOFF]
    )
  ),
  by = .(cell_type, cohort)
]

summary_by_ct[
  ,
  pct_PPH4_gt_cutoff := round(
    100 * n_genes_PPH4_gt_cutoff / n_genes_total,
    1
  )
]

summary_by_ct[
  ,
  cell_type := factor(cell_type, levels = CTs)
]

setorder(summary_by_ct, cell_type, cohort)
summary_by_ct[, cell_type := as.character(cell_type)]

# ==============================================================================
# Wide EUR vs Multi table
# ==============================================================================

summary_wide <- dcast(
  summary_by_ct,
  cell_type ~ cohort,
  value.var = c(
    "n_genes_total",
    "n_genes_PPH4_gt_cutoff",
    "pct_PPH4_gt_cutoff"
  )
)

# ==============================================================================
# Overall summary across all 7 cell types
#
# n_unique_genes_*:
#   each gene is counted once per cohort across all cell types
#
# n_celltype_gene_pairs_*:
#   the same gene in two cell types is counted twice
# ==============================================================================

summary_overall <- all_coloc[
  ,
  .(
    n_unique_genes_total = uniqueN(gene),

    n_unique_genes_PPH4_gt_cutoff = uniqueN(
      gene[!is.na(PP.H4.abf) & PP.H4.abf > PPH4_CUTOFF]
    ),

    n_celltype_gene_pairs = uniqueN(
      paste(cell_type, gene, sep = "|")
    ),

    n_celltype_gene_pairs_PPH4_gt_cutoff = uniqueN(
      paste(
        cell_type[!is.na(PP.H4.abf) & PP.H4.abf > PPH4_CUTOFF],
        gene[!is.na(PP.H4.abf) & PP.H4.abf > PPH4_CUTOFF],
        sep = "|"
      )
    )
  ),
  by = cohort
]

summary_overall[
  ,
  pct_unique_genes_PPH4_gt_cutoff := round(
    100 * n_unique_genes_PPH4_gt_cutoff / n_unique_genes_total,
    1
  )
]

# ==============================================================================
# Save outputs
# ==============================================================================

combined_file <- file.path(
  ROOT,
  sprintf("coloc.%s.7cts.GWAS_eQTL_sigSource.summary.csv", GROUP)
)

by_ct_file <- file.path(
  ROOT,
  sprintf("coloc.%s.7cts.EUR_Multi.summary_by_celltype.csv", GROUP)
)

wide_file <- file.path(
  ROOT,
  sprintf("coloc.%s.7cts.EUR_Multi.summary_by_celltype_WIDE.csv", GROUP)
)

overall_file <- file.path(
  ROOT,
  sprintf("coloc.%s.7cts.EUR_Multi.summary_overall.csv", GROUP)
)

fwrite(all_coloc, combined_file)
fwrite(summary_by_ct, by_ct_file)
fwrite(summary_wide, wide_file)
fwrite(summary_overall, overall_file)

# ==============================================================================
# Print results
# ==============================================================================

cat("\n============================================================\n")
cat("Coloc summary by cell type\n")
cat("PP.H4 >", PPH4_CUTOFF, "\n")
cat("============================================================\n\n")
print(summary_by_ct)

cat("\n============================================================\n")
cat("EUR vs Multi side-by-side\n")
cat("============================================================\n\n")
print(summary_wide)

cat("\n============================================================\n")
cat("Overall across all cell types\n")
cat("============================================================\n\n")
print(summary_overall)

cat("\n============================================================\n")
cat("Saved files\n")
cat("============================================================\n")
cat(combined_file, "\n")
cat(by_ct_file, "\n")
cat(wide_file, "\n")
cat(overall_file, "\n")
