#!/usr/bin/env Rscript

# ==============================================================================
# Script: 03_Healthy_TI_then_rectum_278samples.R
#
# Purpose
# -------
# 1. Load the selected 278-sample Wellcome Sanger Seurat object.
# 2. Keep epithelial cells (Lineage == "Epithelial").
# 3. Read the external Lau-label map:
#       fine_label -> LauLabel
# 4. Map epithelial cell annotations to LauLabel.
# 5. Exclude LauLabel == "DROP" and unmapped cells from pseudobulk.
# 6. Summarize cells and subjects for each LauLabel.
# 7. Generate one genes x subjects pseudobulk expression matrix per LauLabel,
#    aggregated by individual_id.
#
# Important
# ---------
# The source H5AD contained AnnData X but no separate raw-count layer.
# The Seurat RNA/counts layer therefore stores the original X values.
# This script preserves and sums those values exactly; it also checks whether
# sampled non-zero values are integer-like.
# ==============================================================================


# ==============================================================================
# 0. Paths and settings
# ==============================================================================

# --- path resolution ----------------------------------------------------------
# The expression tree was copied to the scTWAS data area on 2026-09-23 and the
# label map now lives under .../scTWAS/data/WellcomeSanger/src/. Each setting
# below takes the first candidate that exists, so the script runs unchanged
# against the new layout or the original /v5000 one. Overrides:
# $WS_RENV_ACTIVATE, $WS_EXPRESSION_DIR, $WS_LAU_LABEL_MAP.

resolve_path <- function(candidates, what, dir = TRUE) {
    candidates <- candidates[nzchar(candidates)]
    hit <- candidates[if (dir) dir.exists(candidates) else file.exists(candidates)][1]
    if (is.na(hit)) {
        stop("Could not find ", what, ". Tried: ",
             paste(candidates, collapse = "; "))
    }
    hit
}

RENV_ACTIVATE <- resolve_path(
    c(
        Sys.getenv("WS_RENV_ACTIVATE"),
        "/data/l2_bioinfo1/liq17/scTWAS/src/renv/activate.R",
        "/data/l2_bioinfo1/liq17/scTWAS/data/VUMC_HTAN_GSE266546/src/renv/activate.R"
    ),
    "renv activate.R",
    dir = FALSE
)

DATA_DIR <- resolve_path(
    c(
        Sys.getenv("WS_EXPRESSION_DIR"),
        "/data/l2_bioinfo1/liq17/scTWAS/data/WellcomeSanger/expression",
        "/v5000/data/l2_bioinfo1/liq17/scTWAS/data/WellcomeSanger/expression",
        "/v5000/data/l2_bioinfo1/scRNA-seq/WellcomeSanger/processed/E-MTAB-16986"
    ),
    "the expression directory"
)

RDS_FILE <- file.path(
    DATA_DIR,
    "Healthy_TI_then_rectum_278samples.seurat.rds"
)

LABEL_MAP_FILE <- resolve_path(
    c(
        Sys.getenv("WS_LAU_LABEL_MAP"),
        "/data/l2_bioinfo1/liq17/scTWAS/data/WellcomeSanger/src/lau7_label_map_29_1.csv",
        "/v5000/data/l2_bioinfo1/liq17/scTWAS/data/WellcomeSanger/src/lau7_label_map_29_1.csv",
        "/v5000/data/l2_bioinfo1/liq17/scTWAS/src/WellcomeSanger/lau7_label_map_29_1.csv"
    ),
    "lau7_label_map_29_1.csv",
    dir = FALSE
)

OUT_DIR <- file.path(
    DATA_DIR,
    "Healthy_TI_then_rectum_278samples_outputs_LauLabel"
)

RAW_BP_DIR <- file.path(
    OUT_DIR,
    "Epithelial_expression_BPCELLS"
)

PB_DIR <- file.path(
    OUT_DIR,
    "pseudobulk_Epithelial_LauLabel_x_individual"
)

LINEAGE_COL <- "Lineage"
EPITHELIAL_VALUE <- "Epithelial"

# The Seurat object previously used Final_labels as the fine cell annotation.
# If a "finelabel" column is already present, it will be preferred.
SOURCE_LABEL_CANDIDATES <- c(
    "finelabel",
    "fine_label",
    "Final_labels"
)

MAP_FINE_COL <- "fine_label"
MAP_LAU_COL <- "LauLabel"

SUBJECT_COL <- "individual_id"
SAMPLE_COL <- "sanger_sample_id"

MIN_CELLS_PER_SUBJECT_CELLTYPE <- 1L
WRITE_PSEUDOBULK_CSV <- TRUE

THREADS <- suppressWarnings(
    as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset = "8"))
)

if (is.na(THREADS) || THREADS < 1L) {
    THREADS <- 8L
}


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

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(PB_DIR, recursive = TRUE, showWarnings = FALSE)


# ==============================================================================
# Helper functions
# ==============================================================================

safe_filename <- function(x) {
    x <- gsub("[^A-Za-z0-9._-]+", "_", x)
    x <- gsub("_+", "_", x)
    gsub("^_|_$", "", x)
}


# Normalize only capitalization and whitespace for matching.
# We intentionally do NOT perform fuzzy matching or silently alter gene names.
normalize_label <- function(x) {
    x <- as.character(x)
    x <- trimws(x)
    x <- gsub("\\s+", " ", x)
    tolower(x)
}


write_csv_gz <- function(x, filename, row_names = FALSE) {
    con <- gzfile(filename, "wt")
    on.exit(close(con), add = TRUE)

    write.csv(
        x,
        con,
        row.names = row_names,
        quote = TRUE
    )
}


# ==============================================================================
# 2. Load Seurat object
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

required_cols <- c(
    LINEAGE_COL,
    SOURCE_LABEL_COL,
    SUBJECT_COL,
    SAMPLE_COL
)

missing_cols <- setdiff(required_cols, colnames(meta))

if (length(missing_cols) > 0) {
    stop(
        "Missing required metadata columns: ",
        paste(missing_cols, collapse = ", ")
    )
}

cat("============================================================\n")
cat("Healthy TI/Rectum LauLabel pseudobulk workflow\n")
cat("============================================================\n")
cat("Input RDS          :", RDS_FILE, "\n")
cat("Label map          :", LABEL_MAP_FILE, "\n")
cat("Fine-label column  :", SOURCE_LABEL_COL, "\n")
cat("Output directory   :", OUT_DIR, "\n")
cat("Threads            :", THREADS, "\n\n")

cat("Whole object:\n")
cat("  Genes       :", nrow(obj), "\n")
cat("  Cells       :", ncol(obj), "\n")
cat("  Samples     :", length(unique(meta[[SAMPLE_COL]])), "\n")
cat("  Individuals :", length(unique(meta[[SUBJECT_COL]])), "\n")


# ==============================================================================
# 3. Read fine_label -> LauLabel map
# ==============================================================================

if (!file.exists(LABEL_MAP_FILE)) {
    stop("Label map does not exist: ", LABEL_MAP_FILE)
}

label_map <- read.csv(
    LABEL_MAP_FILE,
    stringsAsFactors = FALSE,
    check.names = FALSE
)

required_map_cols <- c(
    MAP_FINE_COL,
    MAP_LAU_COL
)

missing_map_cols <- setdiff(
    required_map_cols,
    colnames(label_map)
)

if (length(missing_map_cols) > 0) {
    stop(
        "Label map is missing columns: ",
        paste(missing_map_cols, collapse = ", ")
    )
}

label_map$match_key <- normalize_label(
    label_map[[MAP_FINE_COL]]
)

if (anyDuplicated(label_map$match_key)) {
    dup <- unique(
        label_map[[MAP_FINE_COL]][
            duplicated(label_map$match_key) |
            duplicated(label_map$match_key, fromLast = TRUE)
        ]
    )

    stop(
        "Duplicate fine-label keys in mapping file after whitespace/case ",
        "normalization: ",
        paste(dup, collapse = "; ")
    )
}

write.csv(
    label_map,
    file.path(
        OUT_DIR,
        "01_LauLabel_mapping_file_used.csv"
    ),
    row.names = FALSE
)

cat("\nLau mapping file:\n")
cat("  Mapping rows :", nrow(label_map), "\n")
cat(
    "  LauLabels    :",
    paste(sort(unique(label_map[[MAP_LAU_COL]])), collapse = ", "),
    "\n"
)


# ==============================================================================
# 4. Keep epithelial cells and map fine labels to LauLabel
# ==============================================================================

epi_cells <- rownames(meta)[
    !is.na(meta[[LINEAGE_COL]]) &
    meta[[LINEAGE_COL]] == EPITHELIAL_VALUE
]

if (length(epi_cells) == 0) {
    stop("No epithelial cells found.")
}

epi_meta <- meta[
    epi_cells,
    ,
    drop = FALSE
]

# Keep the original annotation in an explicit output column called finelabel.
epi_meta$finelabel <- as.character(
    epi_meta[[SOURCE_LABEL_COL]]
)

epi_meta$label_match_key <- normalize_label(
    epi_meta$finelabel
)

map_idx <- match(
    epi_meta$label_match_key,
    label_map$match_key
)

epi_meta$LauLabel <- label_map[[MAP_LAU_COL]][
    map_idx
]

# Include mapping note if present.
if ("Note" %in% colnames(label_map)) {
    epi_meta$LauLabel_Note <- label_map$Note[
        map_idx
    ]
}

cat("\nEpithelial cells before LauLabel filtering:\n")
cat("  Cells       :", nrow(epi_meta), "\n")
cat("  Samples     :", length(unique(epi_meta[[SAMPLE_COL]])), "\n")
cat("  Individuals :", length(unique(epi_meta[[SUBJECT_COL]])), "\n")


# ==============================================================================
# 5. Mapping audit
# ==============================================================================

mapping_audit <- epi_meta %>%
    count(
        finelabel,
        LauLabel,
        name = "n_cells"
    ) %>%
    arrange(
        is.na(LauLabel),
        finelabel,
        LauLabel
    )

write.csv(
    mapping_audit,
    file.path(
        OUT_DIR,
        "02_finelabel_to_LauLabel_mapping_audit.csv"
    ),
    row.names = FALSE
)


# Fine labels that were not found in the map.
unmapped <- epi_meta %>%
    filter(
        is.na(LauLabel) |
        trimws(LauLabel) == ""
    ) %>%
    count(
        finelabel,
        name = "n_cells"
    ) %>%
    arrange(
        desc(n_cells)
    )

write.csv(
    unmapped,
    file.path(
        OUT_DIR,
        "03_unmapped_finelabels.csv"
    ),
    row.names = FALSE
)

if (nrow(unmapped) > 0) {
    warning(
        sum(unmapped$n_cells),
        " epithelial cells have fine labels not found in ",
        "lau7_label_map_29_1.csv. See 03_unmapped_finelabels.csv."
    )
}


# Summarize DROP cells separately.
drop_summary <- epi_meta %>%
    filter(
        !is.na(LauLabel),
        LauLabel == "DROP"
    ) %>%
    count(
        finelabel,
        name = "n_cells"
    ) %>%
    arrange(
        desc(n_cells)
    )

write.csv(
    drop_summary,
    file.path(
        OUT_DIR,
        "04_DROPPED_finelabels.csv"
    ),
    row.names = FALSE
)


# ==============================================================================
# 6. Keep analysis LauLabels
#
# Valid Lau labels are mapped labels other than "DROP".
# Unmapped cells are also excluded.
# ==============================================================================

keep_lau <- (
    !is.na(epi_meta$LauLabel) &
    trimws(epi_meta$LauLabel) != "" &
    epi_meta$LauLabel != "DROP"
)

n_epi_before <- nrow(epi_meta)
n_unmapped_cells <- sum(
    is.na(epi_meta$LauLabel) |
    trimws(epi_meta$LauLabel) == ""
)
n_drop_cells <- sum(
    !is.na(epi_meta$LauLabel) &
    epi_meta$LauLabel == "DROP"
)

epi_meta <- epi_meta[
    keep_lau,
    ,
    drop = FALSE
]

epi_cells <- rownames(epi_meta)

cat("\nAfter LauLabel filtering:\n")
cat("  Epithelial cells before :", n_epi_before, "\n")
cat("  Unmapped cells excluded :", n_unmapped_cells, "\n")
cat("  DROP cells excluded     :", n_drop_cells, "\n")
cat("  Cells retained          :", nrow(epi_meta), "\n")
cat("  Individuals retained    :", length(unique(epi_meta[[SUBJECT_COL]])), "\n")
cat(
    "  LauLabels retained     :",
    paste(sort(unique(epi_meta$LauLabel)), collapse = ", "),
    "\n"
)

write_csv_gz(
    epi_meta,
    file.path(
        OUT_DIR,
        "05_Epithelial_cell_metadata_LauLabel.csv.gz"
    ),
    row_names = TRUE
)


# ==============================================================================
# 7. Summarize cells and subjects for each LauLabel
# ==============================================================================

subject_lau_ncells <- epi_meta %>%
    count(
        .data[[SUBJECT_COL]],
        LauLabel,
        name = "n_cells"
    ) %>%
    arrange(
        LauLabel,
        .data[[SUBJECT_COL]]
    )

write.csv(
    subject_lau_ncells,
    file.path(
        OUT_DIR,
        "06_subject_x_LauLabel_ncells.csv"
    ),
    row.names = FALSE
)


lau_summary <- subject_lau_ncells %>%
    group_by(
        LauLabel
    ) %>%
    summarise(
        n_cells = sum(n_cells),
        n_subjects = n_distinct(.data[[SUBJECT_COL]]),
        min_cells_per_subject = min(n_cells),
        median_cells_per_subject = median(n_cells),
        max_cells_per_subject = max(n_cells),
        .groups = "drop"
    ) %>%
    mutate(
        percent_retained_epithelial_cells =
            100 * n_cells / sum(n_cells)
    ) %>%
    arrange(
        desc(n_cells)
    )

write.csv(
    lau_summary,
    file.path(
        OUT_DIR,
        "07_LauLabel_cells_subjects_summary.csv"
    ),
    row.names = FALSE
)

cat("\n============================================================\n")
cat("LauLabel cell / subject summary\n")
cat("============================================================\n")
print(lau_summary)


# Subject metadata
subject_metadata_cols <- intersect(
    c(
        SUBJECT_COL,
        SAMPLE_COL,
        "biopsy_type",
        "tissue",
        "disease_status",
        "age_binned_10",
        "sex",
        "vcf_sample_order"
    ),
    colnames(epi_meta)
)

subject_metadata <- epi_meta[
    !duplicated(epi_meta[[SUBJECT_COL]]),
    subject_metadata_cols,
    drop = FALSE
]

subject_metadata <- subject_metadata[
    order(subject_metadata[[SUBJECT_COL]]),
    ,
    drop = FALSE
]

write.csv(
    subject_metadata,
    file.path(
        OUT_DIR,
        "08_Epithelial_subject_metadata.csv"
    ),
    row.names = FALSE
)


# ==============================================================================
# 8. Extract epithelial RNA/counts layer
# ==============================================================================

counts_all <- SeuratObject::LayerData(
    object = obj,
    assay = "RNA",
    layer = "counts"
)

epi_idx <- match(
    epi_cells,
    colnames(counts_all)
)

if (anyNA(epi_idx)) {
    stop(
        sum(is.na(epi_idx)),
        " retained epithelial cells are missing from RNA/counts."
    )
}

epi_counts <- counts_all[
    ,
    epi_idx,
    drop = FALSE
]

stopifnot(
    identical(
        colnames(epi_counts),
        rownames(epi_meta)
    )
)

cat(
    "\nRetained epithelial matrix:",
    nrow(epi_counts),
    "genes x",
    ncol(epi_counts),
    "cells\n"
)


# ==============================================================================
# 9. Check whether stored values are integer-like
# ==============================================================================

set.seed(12345)

gene_check <- sample(
    seq_len(nrow(epi_counts)),
    min(100L, nrow(epi_counts))
)

cell_check <- sample(
    seq_len(ncol(epi_counts)),
    min(1000L, ncol(epi_counts))
)

check_mat <- as.matrix(
    epi_counts[
        gene_check,
        cell_check,
        drop = FALSE
    ]
)

nz <- check_mat[
    is.finite(check_mat) &
    check_mat != 0
]

integer_like <- if (length(nz) == 0) {
    NA
} else {
    all(abs(nz - round(nz)) < 1e-7)
}

write.csv(
    data.frame(
        sampled_nonzero_values = length(nz),
        integer_like = integer_like,
        min_value = ifelse(length(nz) > 0, min(nz), NA),
        max_value = ifelse(length(nz) > 0, max(nz), NA)
    ),
    file.path(
        OUT_DIR,
        "09_expression_integer_check.csv"
    ),
    row.names = FALSE
)

cat("Integer-like sampled values:", integer_like, "\n")

if (identical(integer_like, FALSE)) {
    warning(
        "Non-integer values detected in RNA/counts. ",
        "Pseudobulk files will be sums of stored AnnData X values, ",
        "not verified raw UMI counts."
    )
}


# ==============================================================================
# 10. Save retained epithelial expression matrix as BPCells
# ==============================================================================

if (!dir.exists(RAW_BP_DIR)) {
    BPCells::write_matrix_dir(
        mat = epi_counts,
        dir = RAW_BP_DIR
    )
}


# ==============================================================================
# 11. Pseudobulk each LauLabel by individual_id
#
# Output:
#   genes x subjects for each LauLabel
# ==============================================================================

lau_labels <- sort(
    unique(
        epi_meta$LauLabel
    )
)

pb_manifest <- vector(
    "list",
    length(lau_labels)
)


for (i in seq_along(lau_labels)) {

    lau <- lau_labels[[i]]
    safe_lau <- safe_filename(lau)

    cat(
        "\n[",
        i,
        "/",
        length(lau_labels),
        "] ",
        lau,
        "\n",
        sep = ""
    )

    ct_idx <- which(
        epi_meta$LauLabel == lau
    )

    ct_meta <- epi_meta[
        ct_idx,
        ,
        drop = FALSE
    ]

    ct_subject_ncells <- ct_meta %>%
        count(
            .data[[SUBJECT_COL]],
            name = "n_cells"
        )

    keep_subjects <- ct_subject_ncells %>%
        filter(
            n_cells >= MIN_CELLS_PER_SUBJECT_CELLTYPE
        ) %>%
        pull(
            .data[[SUBJECT_COL]]
        )

    keep_ct <- ct_meta[[SUBJECT_COL]] %in% keep_subjects

    ct_meta <- ct_meta[
        keep_ct,
        ,
        drop = FALSE
    ]

    ct_idx <- ct_idx[
        keep_ct
    ]

    if (nrow(ct_meta) == 0) {
        next
    }

    ct_counts <- epi_counts[
        ,
        ct_idx,
        drop = FALSE
    ]

    stopifnot(
        identical(
            colnames(ct_counts),
            rownames(ct_meta)
        )
    )

    # Sum all cells for the same individual within this LauLabel.
    pb <- BPCells::pseudobulk_matrix(
        mat = ct_counts,
        cell_groups = as.character(
            ct_meta[[SUBJECT_COL]]
        ),
        method = "sum",
        threads = THREADS
    )

    pb <- pb[
        ,
        order(colnames(pb)),
        drop = FALSE
    ]

    # RDS: efficient R-native genes x subjects matrix.
    rds_file <- file.path(
        PB_DIR,
        paste0(
            "pseudobulk_expression_",
            safe_lau,
            ".rds"
        )
    )

    saveRDS(
        pb,
        rds_file,
        compress = FALSE
    )

    # CSV.GZ: convenient portable genes x subjects matrix.
    if (WRITE_PSEUDOBULK_CSV) {

        csv_file <- file.path(
            PB_DIR,
            paste0(
                "pseudobulk_expression_",
                safe_lau,
                ".csv.gz"
            )
        )

        pb_df <- data.frame(
            gene = rownames(pb),
            as.data.frame(
                pb,
                check.names = FALSE
            ),
            check.names = FALSE
        )

        write_csv_gz(
            pb_df,
            csv_file,
            row_names = FALSE
        )

        rm(pb_df)
    }

    ct_subject_ncells <- ct_subject_ncells[
        ct_subject_ncells[[SUBJECT_COL]] %in% colnames(pb),
        ,
        drop = FALSE
    ]

    write.csv(
        ct_subject_ncells,
        file.path(
            PB_DIR,
            paste0(
                "subject_ncells_",
                safe_lau,
                ".csv"
            )
        ),
        row.names = FALSE
    )

    pb_manifest[[i]] <- data.frame(
        LauLabel = lau,
        n_cells = sum(ct_subject_ncells$n_cells),
        n_subjects = ncol(pb),
        min_cells_per_subject = min(ct_subject_ncells$n_cells),
        median_cells_per_subject = median(ct_subject_ncells$n_cells),
        max_cells_per_subject = max(ct_subject_ncells$n_cells),
        pseudobulk_rds = basename(rds_file),
        stringsAsFactors = FALSE
    )

    cat(
        "  cells   : ",
        sum(ct_subject_ncells$n_cells),
        "\n",
        "  subjects: ",
        ncol(pb),
        "\n",
        sep = ""
    )

    rm(ct_counts, ct_meta, pb)
    gc(verbose = FALSE)
}


# ==============================================================================
# 12. Pseudobulk manifest + workflow summary
# ==============================================================================

pb_manifest <- Filter(
    Negate(is.null),
    pb_manifest
)

if (length(pb_manifest) > 0) {

    pb_manifest <- bind_rows(
        pb_manifest
    ) %>%
        arrange(
            desc(n_cells)
        )

    write.csv(
        pb_manifest,
        file.path(
            OUT_DIR,
            "10_LauLabel_pseudobulk_manifest.csv"
        ),
        row.names = FALSE
    )

    cat("\n============================================================\n")
    cat("Pseudobulk summary\n")
    cat("============================================================\n")
    print(pb_manifest)
}


writeLines(
    c(
        paste0("Input RDS: ", RDS_FILE),
        paste0("Label map: ", LABEL_MAP_FILE),
        paste0("Source fine-label column: ", SOURCE_LABEL_COL),
        paste0("Total cells: ", ncol(obj)),
        paste0("Total subjects: ", length(unique(meta[[SUBJECT_COL]]))),
        paste0("Epithelial cells before mapping filter: ", n_epi_before),
        paste0("Unmapped epithelial cells excluded: ", n_unmapped_cells),
        paste0("DROP epithelial cells excluded: ", n_drop_cells),
        paste0("Retained epithelial cells: ", nrow(epi_meta)),
        paste0(
            "Retained epithelial subjects: ",
            length(unique(epi_meta[[SUBJECT_COL]]))
        ),
        paste0(
            "Retained LauLabels: ",
            paste(sort(unique(epi_meta$LauLabel)), collapse = ", ")
        ),
        paste0("Integer-like stored expression: ", integer_like),
        paste0("Pseudobulk directory: ", PB_DIR)
    ),
    file.path(
        OUT_DIR,
        "11_workflow_summary.txt"
    )
)

writeLines(
    capture.output(sessionInfo()),
    file.path(
        OUT_DIR,
        "12_sessionInfo.txt"
    )
)


cat("\n============================================================\n")
cat("Completed successfully\n")
cat("============================================================\n")
cat(
    "LauLabel summary  : ",
    file.path(
        OUT_DIR,
        "07_LauLabel_cells_subjects_summary.csv"
    ),
    "\n",
    sep = ""
)
cat("Pseudobulk output : ", PB_DIR, "\n", sep = "")
