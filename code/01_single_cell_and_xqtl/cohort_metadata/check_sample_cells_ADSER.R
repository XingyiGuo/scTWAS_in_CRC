#!/usr/bin/env Rscript

# check_sample_cells_ADSER.R
# ------------------------------------------------------------------------------
# Count the cells one ADSER specimen contributes, per cell type, and cross-check
# that number against the bookkeeping tables that are supposed to agree with it.
#
#   keys     SampleID x CellType_cb
#
# WHICH SOURCE IS AUTHORITATIVE
# -----------------------------
# Not the .rds. The script that actually summed cells into
# pseudobulk_counts_<CT>.csv for ADSER is not on any mounted drive - the local
# pipeline step (src/02_prepare_pseudobulk_expression_ASC_SSC_bashconfig.py)
# CONSUMES those files and does not create them - so no file here proves which
# object the summation ran on.
#
# What IS established numerically: grouping
#   ..._CellOracleR_metadata.addHTA11_6407_8890_8930.csv  (100,992 cells,
#                                                          107 specimens)
# by CellType_cb x SampleID reproduces all 628 n_cells values in
# 9CT_expression_id_matched_sample_metadata.csv and the match summary exactly
# (81,929 cells). That makes this CSV the correct cell-count reference whatever
# object the pseudobulk was computed from. The companion .rds shares its name
# stem but not its cell set - see READ_METADATA_ONLY below.
#
# Cross-checks (each optional - skipped with a message if the file is absent):
#   1. ..._CellOracleR_metadata.addHTA11_6407_8890_8930.csv
#        the cell-level metadata exported from this same object; must agree
#        exactly.
#   2. 9CT_expression_id_matched_sample_metadata.csv
#        the TWAS panel. Agrees only for cell types this specimen has an
#        expression column for; a cell type present in the object but absent
#        here was not modelled.
#   3. scTWAS_cells_per_sample_by_celltype.csv
#        the delivered table. NA = not modelled, so object counts can exceed it.
#
# Usage:
#   Rscript check_sample_cells_ADSER.R                      # default specimen
#   Rscript check_sample_cells_ADSER.R HTA11_8755_2000001011
#
# Memory: readRDS on a 445 MB Seurat object needs roughly 4-6 GB. On a small
# machine use READ_METADATA_ONLY <- TRUE, which reads the exported cell-level
# CSV instead and gives the same counts without materialising the object.
# ------------------------------------------------------------------------------

# Seurat is needed only to deserialise the object; the metadata-only path below
# has no package dependencies at all. Load it lazily so this script still runs
# on a machine without Seurat installed.

SAMPLE_ID <- {
  a <- commandArgs(trailingOnly = TRUE)
  if (length(a)) a[1] else "HTA11_10167_2000001011"
}

# TRUE (default) = count from the exported cell-level metadata CSV.
# FALSE           = count from the .rds object.
#
# The CSV is the DEFAULT on purpose. The two are not the same cell set:
#   ..._CellOracleR_input.rds                       2026-09-07 17:25
#   by_celltype_.../*_CellOracleR_metadata.csv      2026-09-07 17:25-17:26
#                                                     -> 97,277 cells, and the
#                                                        three PCA83 add-on
#                                                        specimens are ABSENT
#   ..._CellOracleR_metadata.csv                    2026-09-07 20:00
#   ..._CellOracleR_metadata.addHTA11_...csv        2026-09-07 20:33
#                                                     -> 100,992 cells / 107
#                                                        specimens, add-ons
#                                                        INCLUDED (3,715 cells)
# This is not just a timestamp inference. The two producing scripts in
# VUMC_HTAN_DIS_VAL_ADSER/src say it outright:
#
#   HTAN_ADSER_CellOracleRsave_mergeQC_PCA83_HPtoSER_ASCSSCfilter.R
#     "Purpose: Prepare Vanderbilt HTAN epithelial DIS + VAL + PCA83 data for
#      downstream CellOracleR analyses."  It writes the AllCT .rds and its
#      metadata CSV, plus the nine per-cell-type .rds/.csv pairs. It contains
#      no aggregation step of any kind - no AggregateExpression, no pseudobulk.
#
#   combine_PCA83_generate_metadata_add3.sh
#     Appends the three PCA83 add-on specimens. Its only writes are
#     obs.to_csv(OUT) and combined.to_csv(OUT) - the file contains no saveRDS
#     and no .rds path at all. The add-on merge therefore updated the METADATA
#     CSV only; the .rds was never regenerated and is 3 specimens short.
#
# The RDS branch below still tests this at run time.
READ_METADATA_ONLY <- TRUE

# ---- paths -------------------------------------------------------------------
# ACCRE and the mounted L: drive hold the same tree; take whichever exists.

ADSER_ROOT <- Filter(dir.exists, c(
  "/data/l2_bioinfo1/liq17/scTWAS/data/VUMC_HTAN/VUMC_HTAN_DIS_VAL_ADSER/Results",
  file.path("L:", "scRNA-seq", "VUMC_HTAN", "VUMC_HTAN_DIS_VAL_ADSER", "Results")
))[1]
if (is.na(ADSER_ROOT)) stop("ADSER Results directory not found")

OBJ_RDS <- file.path(
  ADSER_ROOT,
  "HTAN_DISVAL_PCA83_ADSER_HPtoSER_AllCT_ASCSSCfiltered_CellOracleR_input.rds"
)
OBJ_META_CSV <- file.path(
  ADSER_ROOT,
  paste0("HTAN_DISVAL_PCA83_ADSER_HPtoSER_AllCT_ASCSSCfiltered_",
         "CellOracleR_metadata.addHTA11_6407_8890_8930.csv")
)

META_DIR <- Filter(dir.exists, c(
  file.path("L:", "scRNA-seq", "TWAS_scRNA_samples_cells_meta"),
  file.path(ADSER_ROOT, "TWAS_scRNA_meta")
))[1]

PANEL_CSV <- if (!is.na(META_DIR))
  file.path(META_DIR, "9CT_expression_id_matched_sample_metadata.csv") else NA
TABLE_CSV <- if (!is.na(META_DIR))
  file.path(META_DIR, "scTWAS_cells_per_sample_by_celltype.csv") else NA

CT9 <- c("ABS", "GOB", "CT", "EE", "STM", "TUF", "TAC", "ASC", "SSC")

# ==============================================================================
# 1. Counts from the object
# ==============================================================================

if (READ_METADATA_ONLY) {

  stopifnot(file.exists(OBJ_META_CSV))
  md <- read.csv(OBJ_META_CSV, colClasses = "character", check.names = FALSE)
  md <- md[, c("SampleID", "CellType_cb")]
  src <- paste("metadata CSV:", basename(OBJ_META_CSV))

} else {

  stopifnot(file.exists(OBJ_RDS))

  if (!requireNamespace("SeuratObject", quietly = TRUE)) {
    stop("SeuratObject is required to deserialise ", basename(OBJ_RDS),
         ".\nEither install it, or set READ_METADATA_ONLY <- TRUE to read the ",
         "exported cell-level metadata CSV instead (same counts, no object).")
  }
  suppressPackageStartupMessages(library(SeuratObject))

  message("Reading ", OBJ_RDS, " (needs roughly 4-6 GB) ...")
  obj <- readRDS(OBJ_RDS)

  md <- obj@meta.data[, c("SampleID", "CellType_cb")]
  md$SampleID    <- as.character(md$SampleID)
  md$CellType_cb <- as.character(md$CellType_cb)

  message("Object: ", nrow(md), " cells x ",
          length(unique(md$SampleID)), " specimens")
  src <- paste("Seurat object:", basename(OBJ_RDS))

  # Is this object the same cell set as the exported metadata, or the earlier
  # pre-add-on build? Report, do not guess.
  if (file.exists(OBJ_META_CSV)) {
    ref <- read.csv(OBJ_META_CSV, colClasses = "character",
                    check.names = FALSE)[, c("SampleID", "CellType_cb")]
    missing_spec <- setdiff(unique(ref$SampleID), unique(md$SampleID))
    if (length(missing_spec) || nrow(ref) != nrow(md)) {
      warning(
        "Object is NOT the same cell set as ", basename(OBJ_META_CSV), ": ",
        nrow(md), " vs ", nrow(ref), " cells",
        if (length(missing_spec))
          paste0("; specimens absent from the object: ",
                 paste(missing_spec, collapse = ", ")),
        "\nCounts below are from the object and will understate any specimen ",
        "it is missing. Set READ_METADATA_ONLY <- TRUE for the complete set.",
        call. = FALSE, immediate. = TRUE
      )
    } else {
      message("Object matches ", basename(OBJ_META_CSV), " exactly.")
    }
  }
}

hit <- md[md$SampleID == SAMPLE_ID, , drop = FALSE]

if (!nrow(hit)) {
  stop("SampleID '", SAMPLE_ID, "' not found. Present ids look like: ",
       paste(head(sort(unique(md$SampleID)), 3), collapse = ", "))
}

obj_counts <- table(factor(hit$CellType_cb, levels = sort(unique(md$CellType_cb))))
obj_counts <- obj_counts[obj_counts > 0]

cat("\n=== ", SAMPLE_ID, " ===\n", sep = "")
cat("source: ", src, "\n", sep = "")
cat("total cells in object: ", nrow(hit), "\n\n", sep = "")
print(as.data.frame(obj_counts, stringsAsFactors = FALSE,
                    responseName = "cells_in_object"))

# ==============================================================================
# 2. Cross-checks
# ==============================================================================

cmp <- data.frame(
  CellType        = CT9,
  cells_in_object = as.integer(table(factor(hit$CellType_cb, levels = CT9))),
  stringsAsFactors = FALSE
)

# 2a. exported cell-level metadata (only meaningful when we read the .rds)
if (!READ_METADATA_ONLY && file.exists(OBJ_META_CSV)) {
  em <- read.csv(OBJ_META_CSV, colClasses = "character", check.names = FALSE)
  em <- em[em$SampleID == SAMPLE_ID, "CellType_cb"]
  cmp$exported_metadata <- as.integer(table(factor(em, levels = CT9)))
}

# 2b. TWAS expression panel
if (!is.na(PANEL_CSV) && file.exists(PANEL_CSV)) {
  pm <- read.csv(PANEL_CSV, check.names = FALSE)
  pm <- pm[as.character(pm$SampleID) == SAMPLE_ID, ]
  cmp$panel_n_cells <- pm$n_cells[match(CT9, pm$expression_celltype)]
} else {
  message("Panel table not found - skipping that check")
}

# 2c. delivered table
if (!is.na(TABLE_CSV) && file.exists(TABLE_CSV)) {
  dt <- read.csv(TABLE_CSV, check.names = FALSE)
  dt <- dt[as.character(dt$Expression_id) == SAMPLE_ID, ]
  if (nrow(dt) == 1) {
    cmp$delivered_table <- as.integer(unlist(dt[, CT9]))
  } else {
    message("Specimen is not a row of the delivered table ",
            "(second specimen under ROW_SET = 'sheet_only')")
  }
} else {
  message("Delivered table not found - skipping that check")
}

cat("\n")
print(cmp, row.names = FALSE)

# A cell type the object has but the panel does not = present in the tissue,
# never modelled. The reverse would be a real inconsistency.
if ("panel_n_cells" %in% names(cmp)) {
  not_modelled <- cmp$CellType[cmp$cells_in_object > 0 & is.na(cmp$panel_n_cells)]
  disagree     <- cmp$CellType[!is.na(cmp$panel_n_cells) &
                                 cmp$panel_n_cells != cmp$cells_in_object]
  cat("\nnot modelled (cells present, no expression column): ",
      if (length(not_modelled)) paste(not_modelled, collapse = ", ") else "none",
      "\ndisagreements between object and panel: ",
      if (length(disagree)) paste(disagree, collapse = ", ") else "none", "\n", sep = "")
}
