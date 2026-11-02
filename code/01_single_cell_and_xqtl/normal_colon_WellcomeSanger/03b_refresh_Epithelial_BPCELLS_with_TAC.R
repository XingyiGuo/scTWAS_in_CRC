#!/usr/bin/env Rscript

# ==============================================================================
# Script: 03b_refresh_Epithelial_BPCELLS_with_TAC.R
#
#   >>> RUN ON THE CLUSTER (needs /v5000 and the renv library) <<<
#
# Purpose
# -------
# Refresh ONLY the on-disk epithelial expression matrix
#
#     Healthy_TI_then_rectum_278samples_outputs_LauLabel/
#         Epithelial_expression_BPCELLS
#
# so that it covers every retained epithelial cell of the current label map -
# in particular the 13,297 TAC cells and the 10,844 ABS cells that the existing
# copy lacks - and thereby allows a cell-level UMAP of all SEVEN LauLabels.
#
# This script is a companion to 03_Healthy_TI_then_rectum_278samples_updated.R,
# which is deliberately left untouched so the published pipeline stays
# reproducible. It re-uses that script's cell-selection logic verbatim
# (Lineage == "Epithelial", Final_labels -> LauLabel via lau7_label_map_29_1.csv,
# drop unmapped and "DROP") and then stops. It does NOT write pseudobulk.
#
# Why the existing matrix is incomplete
# -------------------------------------
# Step 10 of script 03 writes the matrix under `if (!dir.exists(RAW_BP_DIR))`,
# i.e. only when the directory is absent. After the label map changed, the
# directory already existed and was therefore never rewritten: it holds 334,066
# of the 358,207 currently retained cells. The pseudobulk matrices were not
# affected, because script 03 builds them from the RDS counts layer rather than
# from this directory.
#
# What it does NOT touch (guaranteed, and verified at the end)
# -----------------------------------------------------------
#   pseudobulk_Epithelial_LauLabel_x_individual/*      (the 7 {ct}.rds, the
#                                                       .csv.gz, subject_ncells)
#   01..12 audit files, including 05_Epithelial_cell_metadata_LauLabel.csv.gz
# The pseudobulk directory is fingerprinted (size + mtime of every file) before
# and after the rewrite, and the script stops with an error if anything moved.
#
# Outputs
# -------
#   Epithelial_expression_BPCELLS                 refreshed, all 7 LauLabels
#   Epithelial_expression_BPCELLS_stale_<stamp>   previous copy, kept unless
#                                                 KEEP_STALE_COPY <- FALSE
#   03b_BPCELLS_refresh_log.csv                   cells/subjects per LauLabel
#   03b_sessionInfo.txt                           session record (own filename,
#                                                 does not clobber 12_*)
# ==============================================================================


# ==============================================================================
# 0. Paths and settings (identical to script 03 where they overlap)
# ==============================================================================

RENV_ACTIVATE <- paste0(
    "/data/l2_bioinfo1/liq17/scTWAS/src/renv/activate.R"
)

# Expression root. The E-MTAB-16986 tree was copied to the scTWAS data area on
# 2026-09-23, so that copy is tried first; the original /v5000 location and an
# override in $WS_EXPRESSION_DIR remain as alternatives. First hit wins.
DATA_DIR_CANDIDATES <- c(
    Sys.getenv("WS_EXPRESSION_DIR"),
    "/data/l2_bioinfo1/liq17/scTWAS/data/WellcomeSanger/expression",
    "/v5000/data/l2_bioinfo1/liq17/scTWAS/data/WellcomeSanger/expression",
    "/v5000/data/l2_bioinfo1/scRNA-seq/WellcomeSanger/processed/E-MTAB-16986"
)

DATA_DIR <- DATA_DIR_CANDIDATES[
    nzchar(DATA_DIR_CANDIDATES) & dir.exists(DATA_DIR_CANDIDATES)
][1]

if (is.na(DATA_DIR)) {
    stop(
        "Could not find the expression directory. Tried: ",
        paste(DATA_DIR_CANDIDATES[nzchar(DATA_DIR_CANDIDATES)], collapse = "; "),
        ". Set WS_EXPRESSION_DIR and re-run."
    )
}

RDS_FILE <- file.path(
    DATA_DIR,
    "Healthy_TI_then_rectum_278samples.seurat.rds"
)

# Label map. These copies are identical in content; the 2026-09-23 run
# re-derived every LauLabel from it and matched
# 05_Epithelial_cell_metadata_LauLabel.csv.gz exactly.
LABEL_MAP_CANDIDATES <- c(
    Sys.getenv("WS_LAU_LABEL_MAP"),
    "/data/l2_bioinfo1/liq17/scTWAS/data/WellcomeSanger/src/lau7_label_map_29_1.csv",
    "/v5000/data/l2_bioinfo1/liq17/scTWAS/data/WellcomeSanger/src/lau7_label_map_29_1.csv",
    "/v5000/data/l2_bioinfo1/liq17/scTWAS/src/WellcomeSanger/lau7_label_map_29_1.csv"
)

LABEL_MAP_FILE <- LABEL_MAP_CANDIDATES[
    nzchar(LABEL_MAP_CANDIDATES) & file.exists(LABEL_MAP_CANDIDATES)
][1]

if (is.na(LABEL_MAP_FILE)) {
    stop(
        "Could not find lau7_label_map_29_1.csv. Tried: ",
        paste(LABEL_MAP_CANDIDATES[nzchar(LABEL_MAP_CANDIDATES)], collapse = "; "),
        ". Set WS_LAU_LABEL_MAP and re-run."
    )
}

OUT_DIR <- file.path(
    DATA_DIR,
    "Healthy_TI_then_rectum_278samples_outputs_LauLabel"
)

RAW_BP_DIR <- file.path(OUT_DIR, "Epithelial_expression_BPCELLS")

PB_DIR <- file.path(OUT_DIR, "pseudobulk_Epithelial_LauLabel_x_individual")

LINEAGE_COL <- "Lineage"
EPITHELIAL_VALUE <- "Epithelial"

SOURCE_LABEL_CANDIDATES <- c("finelabel", "fine_label", "Final_labels")

MAP_FINE_COL <- "fine_label"
MAP_LAU_COL <- "LauLabel"

SUBJECT_COL <- "individual_id"

# Keep the current (incomplete) matrix as Epithelial_expression_BPCELLS_stale_*
# instead of deleting it. The refreshed matrix needs roughly the same space
# again, so set this to FALSE if the filesystem is tight.
KEEP_STALE_COPY <- TRUE

# Cross-check the cells selected here against the cell list that produced the
# pseudobulk files (05_Epithelial_cell_metadata_LauLabel.csv.gz). Leave TRUE:
# it is what guarantees the refreshed matrix describes the same cells as the
# pseudobulk already on disk.
CHECK_AGAINST_05_METADATA <- TRUE

# Number of random cells whose column sums are compared against nCount_RNA
# after writing, as a numeric check that the write is faithful.
N_VERIFY_CELLS <- 200L


# ==============================================================================
# 1. Environment
# ==============================================================================

if (!file.exists(RENV_ACTIVATE)) {
    stop("Cannot find renv activation file: ", RENV_ACTIVATE)
}

source(RENV_ACTIVATE)

suppressPackageStartupMessages({
    library(Seurat)
    library(SeuratObject)
    library(BPCells)
    library(dplyr)
})

normalize_label <- function(x) {
    x <- as.character(x)
    x <- trimws(x)
    x <- gsub("\\s+", " ", x)
    tolower(x)
}

fingerprint <- function(dir) {
    if (!dir.exists(dir)) {
        return(data.frame())
    }
    info <- file.info(list.files(dir, full.names = TRUE, recursive = TRUE))
    data.frame(
        path = rownames(info),
        size = info$size,
        mtime = as.character(info$mtime),
        stringsAsFactors = FALSE
    )
}

cat("============================================================\n")
cat("Refresh epithelial BPCells matrix (adds TAC); pseudobulk untouched\n")
cat("============================================================\n")
cat("Input RDS   :", RDS_FILE, "\n")
cat("Label map   :", LABEL_MAP_FILE, "\n")
cat("Target dir  :", RAW_BP_DIR, "\n")
cat("Pseudobulk  :", PB_DIR, " (read-only here)\n\n")

pb_before <- fingerprint(PB_DIR)
cat("pseudobulk files fingerprinted:", nrow(pb_before), "\n\n")


# ==============================================================================
# 2. Load Seurat object and map Final_labels -> LauLabel (as in script 03)
# ==============================================================================

if (!file.exists(RDS_FILE)) {
    stop("Input RDS does not exist: ", RDS_FILE)
}

obj <- readRDS(RDS_FILE)
DefaultAssay(obj) <- "RNA"

meta <- obj@meta.data

SOURCE_LABEL_COL <- SOURCE_LABEL_CANDIDATES[
    SOURCE_LABEL_CANDIDATES %in% colnames(meta)
][1]

if (is.na(SOURCE_LABEL_COL)) {
    stop(
        "Cannot find a cell-level fine annotation column. Tried: ",
        paste(SOURCE_LABEL_CANDIDATES, collapse = ", ")
    )
}

cat("Fine-label column:", SOURCE_LABEL_COL, "\n")
cat("Cells in object  :", nrow(meta), "\n")

label_map <- read.csv(
    LABEL_MAP_FILE,
    stringsAsFactors = FALSE,
    check.names = FALSE
)

if (!all(c(MAP_FINE_COL, MAP_LAU_COL) %in% colnames(label_map))) {
    stop("Label map must contain columns ", MAP_FINE_COL, " and ", MAP_LAU_COL)
}

label_map$match_key <- normalize_label(label_map[[MAP_FINE_COL]])

if (anyDuplicated(label_map$match_key)) {
    stop("Duplicate fine-label keys in the mapping file after normalization.")
}

epi_cells <- rownames(meta)[
    !is.na(meta[[LINEAGE_COL]]) &
    meta[[LINEAGE_COL]] == EPITHELIAL_VALUE
]

epi_meta <- meta[epi_cells, , drop = FALSE]

epi_meta$finelabel <- as.character(epi_meta[[SOURCE_LABEL_COL]])

epi_meta$LauLabel <- label_map[[MAP_LAU_COL]][
    match(normalize_label(epi_meta$finelabel), label_map$match_key)
]

keep_lau <- (
    !is.na(epi_meta$LauLabel) &
    trimws(epi_meta$LauLabel) != "" &
    epi_meta$LauLabel != "DROP"
)

cat("Epithelial cells :", nrow(epi_meta), "\n")
cat("  unmapped       :", sum(is.na(epi_meta$LauLabel) |
                              trimws(epi_meta$LauLabel) == ""), "\n")
cat("  DROP           :", sum(!is.na(epi_meta$LauLabel) &
                              epi_meta$LauLabel == "DROP"), "\n")

epi_meta <- epi_meta[keep_lau, , drop = FALSE]
epi_cells <- rownames(epi_meta)

cat("  retained       :", nrow(epi_meta), "\n")
print(table(epi_meta$LauLabel))


# ==============================================================================
# 3. Cross-check the retained set against the pseudobulk cell list
# ==============================================================================

if (isTRUE(CHECK_AGAINST_05_METADATA)) {

    meta05_file <- file.path(OUT_DIR, "05_Epithelial_cell_metadata_LauLabel.csv.gz")

    if (file.exists(meta05_file)) {

        meta05 <- read.csv(meta05_file, row.names = 1, stringsAsFactors = FALSE)

        only_now <- setdiff(epi_cells, rownames(meta05))
        only_05 <- setdiff(rownames(meta05), epi_cells)

        cat("\nAgainst 05_Epithelial_cell_metadata_LauLabel.csv.gz:\n")
        cat("  cells in that file      :", nrow(meta05), "\n")
        cat("  selected here but not there:", length(only_now), "\n")
        cat("  there but not selected here:", length(only_05), "\n")

        if (length(only_now) > 0 || length(only_05) > 0) {
            stop(
                "The retained cell set differs from the one behind the existing ",
                "pseudobulk files. The label map has changed since the pseudobulk ",
                "was built, so refreshing only this matrix would make the two ",
                "inconsistent. Re-run script 03 in full instead."
            )
        }

        stopifnot(
            identical(
                as.character(meta05[epi_cells, "LauLabel"]),
                as.character(epi_meta$LauLabel)
            )
        )

        cat("  cell sets and LauLabels agree exactly\n")

    } else {
        warning("05_Epithelial_cell_metadata_LauLabel.csv.gz not found; ",
                "cross-check skipped.")
    }
}


# ==============================================================================
# 4. Report what the existing matrix is missing
# ==============================================================================

if (dir.exists(RAW_BP_DIR)) {

    stored_cells <- colnames(BPCells::open_matrix_dir(RAW_BP_DIR))

    absent <- setdiff(epi_cells, stored_cells)
    extra <- setdiff(stored_cells, epi_cells)

    cat("\nExisting matrix:", length(stored_cells), "cells\n")
    cat("  retained cells absent from it:", length(absent), "\n")

    if (length(absent) > 0) {
        print(table(epi_meta[absent, "LauLabel"]))
    }

    cat("  cells in it that are no longer retained:", length(extra), "\n")

    if (length(absent) == 0 && length(extra) == 0) {
        cat("\nMatrix is already complete; nothing to do.\n")
        quit(save = "no", status = 0)
    }

} else {
    cat("\nNo existing matrix at the target path; it will be created.\n")
}


# ==============================================================================
# 5. Slice the RNA/counts layer for the retained cells
# ==============================================================================

counts_all <- SeuratObject::LayerData(obj, assay = "RNA", layer = "counts")

epi_idx <- match(epi_cells, colnames(counts_all))

if (anyNA(epi_idx)) {
    stop(sum(is.na(epi_idx)), " retained cells are missing from RNA/counts.")
}

epi_counts <- counts_all[, epi_idx, drop = FALSE]

stopifnot(identical(colnames(epi_counts), rownames(epi_meta)))

cat("\nMatrix to write:", nrow(epi_counts), "genes x", ncol(epi_counts), "cells\n")


# ==============================================================================
# 6. Move the stale copy aside and write the refreshed matrix
# ==============================================================================

if (dir.exists(RAW_BP_DIR)) {

    stale_dir <- paste0(
        RAW_BP_DIR, "_stale_", format(Sys.time(), "%Y%m%d%H%M%S")
    )

    if (isTRUE(KEEP_STALE_COPY)) {

        cat("Moving the existing matrix to:", basename(stale_dir), "\n")

        if (!file.rename(RAW_BP_DIR, stale_dir)) {
            stop("Could not move the existing directory: ", RAW_BP_DIR)
        }

    } else {

        cat("Deleting the existing matrix (KEEP_STALE_COPY = FALSE)\n")
        unlink(RAW_BP_DIR, recursive = TRUE)

        if (dir.exists(RAW_BP_DIR)) {
            stop("Could not remove the existing directory: ", RAW_BP_DIR)
        }
    }
}

cat("Writing:", RAW_BP_DIR, "\n")

BPCells::write_matrix_dir(mat = epi_counts, dir = RAW_BP_DIR)


# ==============================================================================
# 7. Verify the refreshed matrix
# ==============================================================================

new_mat <- BPCells::open_matrix_dir(RAW_BP_DIR)

stopifnot(identical(colnames(new_mat), epi_cells))
stopifnot(identical(rownames(new_mat), rownames(epi_counts)))
stopifnot(setequal(unique(epi_meta$LauLabel), unique(epi_meta[colnames(new_mat), "LauLabel"])))

cat("\nRefreshed matrix:", nrow(new_mat), "genes x", ncol(new_mat), "cells\n")
cat("LauLabels present:",
    paste(sort(unique(epi_meta[colnames(new_mat), "LauLabel"])), collapse = ", "),
    "\n")
print(table(epi_meta[colnames(new_mat), "LauLabel"]))

# numeric check: column sums must equal nCount_RNA for a random set of cells
set.seed(20260923)

verify_idx <- sample(
    seq_len(ncol(new_mat)),
    min(N_VERIFY_CELLS, ncol(new_mat))
)

verify_sums <- Matrix::colSums(
    as(new_mat[, verify_idx, drop = FALSE], "dgCMatrix")
)

expected_sums <- epi_meta[colnames(new_mat)[verify_idx], "nCount_RNA"]

max_diff <- max(abs(verify_sums - expected_sums))

cat("column-sum check on", length(verify_idx),
    "cells: max |colSum - nCount_RNA| =", max_diff, "\n")

if (!isTRUE(max_diff < 1e-6)) {
    stop("Column sums do not reproduce nCount_RNA; the write is suspect.")
}


# ==============================================================================
# 8. Confirm the pseudobulk directory was not touched
# ==============================================================================

pb_after <- fingerprint(PB_DIR)

pb_cmp <- merge(
    pb_before,
    pb_after,
    by = "path",
    all = TRUE,
    suffixes = c("_before", "_after")
)

pb_changed <- pb_cmp$path[
    is.na(pb_cmp$size_before) |
    is.na(pb_cmp$size_after) |
    pb_cmp$size_before != pb_cmp$size_after |
    pb_cmp$mtime_before != pb_cmp$mtime_after
]

if (length(pb_changed) > 0) {
    stop(
        "Pseudobulk directory changed during this run, which must not happen: ",
        paste(head(basename(pb_changed), 10), collapse = ", ")
    )
}

cat("\npseudobulk directory unchanged:", nrow(pb_after),
    "files, identical size and mtime\n")


# ==============================================================================
# 9. Log
# ==============================================================================

log_df <- epi_meta %>%
    count(LauLabel, name = "n_cells") %>%
    left_join(
        epi_meta %>%
            group_by(LauLabel) %>%
            summarise(n_subjects = n_distinct(.data[[SUBJECT_COL]]),
                      .groups = "drop"),
        by = "LauLabel"
    ) %>%
    mutate(
        bpcells_dir = RAW_BP_DIR,
        refreshed_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
        n_genes = nrow(new_mat)
    ) %>%
    arrange(desc(n_cells))

write.csv(
    log_df,
    file.path(OUT_DIR, "03b_BPCELLS_refresh_log.csv"),
    row.names = FALSE
)

print(as.data.frame(log_df))

writeLines(
    capture.output(sessionInfo()),
    file.path(OUT_DIR, "03b_sessionInfo.txt")
)

cat("\n============================================================\n")
cat("Done. Copy this directory to the Windows side for the UMAP:\n")
cat("  ", RAW_BP_DIR, "\n", sep = "")
cat("Pseudobulk files were not modified.\n")
cat("============================================================\n")
