#!/usr/bin/env Rscript

# ==============================================================================
# Wellcome Sanger E-MTAB-16986
# H5AD -> BPCells on-disk matrix -> Seurat v5
#
# The dataset contains >2 million cells and >7.3 billion non-zero entries.
# BPCells is used so the expression matrix remains on disk.
#
# IMPORTANT:
# The saved Seurat RDS references the BPCells directory. Do not delete, rename,
# or move that directory after creating the RDS.
# ==============================================================================

# ------------------------------------------------------------------------------
# Paths
# ------------------------------------------------------------------------------

# --- path resolution ----------------------------------------------------------
# The expression tree was copied to the scTWAS data area on 2026-09-23; that
# copy is preferred and the original /v5000 location is kept as a fallback, so
# this script runs unchanged against either. Environment overrides:
# $WS_RENV_ACTIVATE, $WS_EXPRESSION_DIR.

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

h5ad_file <- file.path(
    DATA_DIR,
    "multi_tissue-predicted.h5ad"
)

bp_dir <- file.path(
    DATA_DIR,
    "multi_tissue_predicted_BPCELLS"
)

out_rds <- file.path(
    DATA_DIR,
    "multi_tissue-predicted.seurat.rds"
)

summary_dir <- file.path(
    DATA_DIR,
    "multi_tissue_predicted_summary"
)

# ------------------------------------------------------------------------------
# Activate renv and load packages
# ------------------------------------------------------------------------------

if (!file.exists(RENV_ACTIVATE)) {
    stop("renv activation file not found: ", RENV_ACTIVATE)
}
source(RENV_ACTIVATE)

suppressPackageStartupMessages({
    library(Seurat)
    library(BPCells)
    library(dplyr)
})

if (!requireNamespace("Azimuth", quietly = TRUE)) {
    stop("Package 'Azimuth' is required. Run 00_install_BPCell.sh first.")
}

options(Seurat.object.assay.version = "v5")

dir.create(
    summary_dir,
    recursive = TRUE,
    showWarnings = FALSE
)

if (!file.exists(h5ad_file)) {
    stop("Input H5AD file does not exist: ", h5ad_file)
}

cat("============================================================\n")
cat("Wellcome Sanger H5AD -> BPCells -> Seurat\n")
cat("============================================================\n")
cat("Input H5AD : ", h5ad_file, "\n", sep = "")
cat("BPCells dir: ", bp_dir, "\n", sep = "")
cat("Output RDS : ", out_rds, "\n", sep = "")
cat("Summary dir: ", summary_dir, "\n\n", sep = "")

# ------------------------------------------------------------------------------
# 1. Open AnnData X directly with BPCells
# ------------------------------------------------------------------------------

cat("============================================================\n")
cat("Step 1: Open H5AD expression matrix\n")
cat("============================================================\n")

mat <- BPCells::open_matrix_anndata_hdf5(
    h5ad_file
)

print(mat)

cat(
    "Dimensions: ",
    nrow(mat), " features x ",
    ncol(mat), " cells\n\n",
    sep = ""
)

# ------------------------------------------------------------------------------
# 2. Convert to BPCells on-disk format
#
# The conversion is performed only if bp_dir does not already exist.
# If a previous run was interrupted, remove or rename the incomplete directory
# before rerunning this script.
# ------------------------------------------------------------------------------

cat("============================================================\n")
cat("Step 2: Prepare BPCells on-disk matrix\n")
cat("============================================================\n")

if (!dir.exists(bp_dir)) {

    cat("Writing BPCells matrix to:\n", bp_dir, "\n\n", sep = "")

    BPCells::write_matrix_dir(
        mat = mat,
        dir = bp_dir
    )

} else {

    cat(
        "BPCells directory already exists; skipping conversion:\n",
        bp_dir,
        "\n\n",
        sep = ""
    )
}

counts <- BPCells::open_matrix_dir(
    bp_dir
)

print(counts)

cat(
    "Dimensions: ",
    nrow(counts), " features x ",
    ncol(counts), " cells\n\n",
    sep = ""
)

if (nrow(counts) != nrow(mat) || ncol(counts) != ncol(mat)) {
    stop(
        "Existing BPCells matrix dimensions do not match H5AD. ",
        "The BPCells directory may be incomplete or from another dataset."
    )
}

# ------------------------------------------------------------------------------
# 3. Read adata.obs metadata only
#
# This avoids loading the gigantic X matrix through a conventional R sparse
# matrix reader.
# ------------------------------------------------------------------------------

cat("============================================================\n")
cat("Step 3: Read H5AD cell metadata\n")
cat("============================================================\n")

metadata <- Azimuth::LoadH5ADobs(
    path = h5ad_file
)

cat(
    "Metadata dimensions: ",
    nrow(metadata), " cells x ",
    ncol(metadata), " columns\n",
    sep = ""
)

cat("\nMetadata columns:\n")
print(colnames(metadata))

# ------------------------------------------------------------------------------
# 4. Validate cell IDs and reorder metadata
# ------------------------------------------------------------------------------

cat("\n============================================================\n")
cat("Step 4: Validate cell IDs\n")
cat("============================================================\n")

matrix_cells <- colnames(counts)
metadata_cells <- rownames(metadata)

if (is.null(matrix_cells)) {
    stop("The BPCells matrix does not contain cell names.")
}

if (is.null(metadata_cells)) {
    stop("The H5AD metadata does not contain cell row names.")
}

cat("Matrix cells  : ", length(matrix_cells), "\n", sep = "")
cat("Metadata cells: ", length(metadata_cells), "\n", sep = "")

missing_metadata <- setdiff(
    matrix_cells,
    metadata_cells
)

if (length(missing_metadata) > 0) {
    stop(
        length(missing_metadata),
        " matrix cells are absent from metadata. Example: ",
        paste(head(missing_metadata, 5), collapse = ", ")
    )
}

metadata <- metadata[
    matrix_cells,
    ,
    drop = FALSE
]

stopifnot(
    identical(
        rownames(metadata),
        matrix_cells
    )
)

cat("Cell IDs are aligned correctly.\n")

# ------------------------------------------------------------------------------
# 5. Create Seurat object
# ------------------------------------------------------------------------------

cat("\n============================================================\n")
cat("Step 5: Create Seurat object\n")
cat("============================================================\n")

obj <- Seurat::CreateSeuratObject(
    counts = counts,
    meta.data = metadata,
    project = "WellcomeSanger"
)

print(obj)

cat(
    "\nNumber of cells   : ",
    format(ncol(obj), big.mark = ","),
    "\n",
    sep = ""
)

cat(
    "Number of features: ",
    format(nrow(obj), big.mark = ","),
    "\n",
    sep = ""
)

# ------------------------------------------------------------------------------
# 6. Generate metadata summaries
# ------------------------------------------------------------------------------

cat("\n============================================================\n")
cat("Step 6: Generate metadata summaries\n")
cat("============================================================\n")

metadata_columns <- data.frame(
    variable = colnames(obj@meta.data),
    n_unique = vapply(
        obj@meta.data,
        function(x) length(unique(x)),
        integer(1)
    ),
    stringsAsFactors = FALSE
)

metadata_columns <- metadata_columns[
    order(metadata_columns$n_unique),
    ,
    drop = FALSE
]

write.csv(
    metadata_columns,
    file.path(summary_dir, "metadata_columns.csv"),
    row.names = FALSE
)

if ("sanger_sample_id" %in% colnames(obj@meta.data)) {

    cells_per_sample <- obj@meta.data %>%
        count(
            sanger_sample_id,
            name = "n_cells"
        ) %>%
        arrange(desc(n_cells))

    write.csv(
        cells_per_sample,
        file.path(summary_dir, "cells_per_sample.csv"),
        row.names = FALSE
    )

    cat(
        "Number of Sanger samples: ",
        length(unique(obj$sanger_sample_id)),
        "\n",
        sep = ""
    )
}

if ("Final_labels" %in% colnames(obj@meta.data)) {

    cells_per_final_label <- obj@meta.data %>%
        count(
            Final_labels,
            name = "n_cells"
        ) %>%
        mutate(
            percent_cells = 100 * n_cells / sum(n_cells)
        ) %>%
        arrange(desc(n_cells))

    write.csv(
        cells_per_final_label,
        file.path(summary_dir, "cells_per_final_label.csv"),
        row.names = FALSE
    )

    cat("\nLargest Final_labels categories:\n")
    print(head(cells_per_final_label, 20))
}

if (
    all(
        c("sanger_sample_id", "Final_labels") %in%
        colnames(obj@meta.data)
    )
) {

    sample_celltype <- obj@meta.data %>%
        count(
            sanger_sample_id,
            Final_labels,
            name = "n_cells"
        ) %>%
        arrange(
            sanger_sample_id,
            desc(n_cells)
        )

    write.csv(
        sample_celltype,
        file.path(summary_dir, "sample_celltype_counts.csv"),
        row.names = FALSE
    )
}

if (
    all(
        c("tissue", "Final_labels") %in%
        colnames(obj@meta.data)
    )
) {

    tissue_celltype <- obj@meta.data %>%
        count(
            tissue,
            Final_labels,
            name = "n_cells"
        ) %>%
        arrange(
            tissue,
            desc(n_cells)
        )

    write.csv(
        tissue_celltype,
        file.path(summary_dir, "tissue_celltype_counts.csv"),
        row.names = FALSE
    )
}

for (field in c("Final_labels", "Lineage", "predicted_labels", "tissue")) {

    if (field %in% colnames(obj@meta.data)) {

        cat("\n---------------- ", field, " ----------------\n", sep = "")

        print(
            sort(
                table(
                    obj@meta.data[[field]],
                    useNA = "ifany"
                ),
                decreasing = TRUE
            )
        )
    }
}

# ------------------------------------------------------------------------------
# 7. Save Seurat object
# ------------------------------------------------------------------------------

cat("\n============================================================\n")
cat("Step 7: Save Seurat object\n")
cat("============================================================\n")

saveRDS(
    obj,
    out_rds,
    compress = FALSE
)

cat("Saved Seurat RDS:\n", out_rds, "\n\n", sep = "")
cat("BPCells matrix directory:\n", bp_dir, "\n\n", sep = "")
cat("Summary directory:\n", summary_dir, "\n\n", sep = "")

# ------------------------------------------------------------------------------
# 8. Save R session information for reproducibility
# ------------------------------------------------------------------------------

writeLines(
    capture.output(sessionInfo()),
    file.path(
        summary_dir,
        "sessionInfo_h5ad_to_seurat.txt"
    )
)

cat("Pipeline completed successfully.\n")
