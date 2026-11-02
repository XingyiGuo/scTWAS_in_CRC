library(Seurat)
library(BPCells)
library(Azimuth)

# --- path resolution ----------------------------------------------------------
# Expression tree copied to the scTWAS data area on 2026-09-23; the original
# /v5000 location is kept as a fallback, so this script runs unchanged against
# either. Overrides: $WS_RENV_ACTIVATE, $WS_EXPRESSION_DIR.

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
source(RENV_ACTIVATE)

wk <- resolve_path(
    c(
        Sys.getenv("WS_EXPRESSION_DIR"),
        "/data/l2_bioinfo1/liq17/scTWAS/data/WellcomeSanger/expression",
        "/v5000/data/l2_bioinfo1/liq17/scTWAS/data/WellcomeSanger/expression",
        "/v5000/data/l2_bioinfo1/scRNA-seq/WellcomeSanger/processed/E-MTAB-16986"
    ),
    "the expression directory"
)

# ------------------------------------------------
# Files
# ------------------------------------------------

sample_file <- file.path(
    wk,
    "multi_tissue_predicted_summary",
    "one_sample_per_individual_Healthy_TI_then_rectum_individual_metadata.csv"
)

h5ad_file <- file.path(
    wk,
    "multi_tissue-predicted.h5ad"
)

bp_dir <- file.path(
    wk,
    "multi_tissue_predicted_BPCELLS"
)

subset_bp_dir <- file.path(
    wk,
    "Healthy_TI_then_rectum_BPCELLS"
)

out_rds <- file.path(
    wk,
    "Healthy_TI_then_rectum_278samples.seurat.rds"
)


# ------------------------------------------------
# 1. Read 278 selected sample IDs
# ------------------------------------------------

sample_info <- read.csv(
    sample_file,
    stringsAsFactors = FALSE
)

keep_samples <- unique(
    sample_info$sanger_sample_id
)

length(keep_samples)
# 278


# ------------------------------------------------
# 2. Read H5AD metadata
# ------------------------------------------------

metadata <- Azimuth::LoadH5ADobs(
    path = h5ad_file
)

keep_cells <- rownames(metadata)[
    metadata$sanger_sample_id %in% keep_samples
]

length(keep_cells)
# expected: 695199


# Check all 278 samples are represented
found_samples <- unique(
    metadata[
        keep_cells,
        "sanger_sample_id"
    ]
)

setdiff(
    keep_samples,
    found_samples
)


# ------------------------------------------------
# 3. Open full BPCells matrix
# ------------------------------------------------

counts <- BPCells::open_matrix_dir(
    bp_dir
)

# Make sure cell IDs match
stopifnot(
    all(keep_cells %in% colnames(counts))
)


# ------------------------------------------------
# 4. Subset expression matrix
# ------------------------------------------------

cell_idx <- which(
    colnames(counts) %in% keep_cells
)

counts_sub <- counts[
    ,
    cell_idx
]

dim(counts_sub)


# ------------------------------------------------
# 5. Write selected cells to new BPCells directory
# ------------------------------------------------

if (!dir.exists(subset_bp_dir)) {

    BPCells::write_matrix_dir(
        mat = counts_sub,
        dir = subset_bp_dir
    )
}

counts_sub <- BPCells::open_matrix_dir(
    subset_bp_dir
)


# ------------------------------------------------
# 6. Subset and align metadata
# ------------------------------------------------

metadata_sub <- metadata[
    colnames(counts_sub),
    ,
    drop = FALSE
]

stopifnot(
    identical(
        rownames(metadata_sub),
        colnames(counts_sub)
    )
)


# ------------------------------------------------
# 7. Add individual-level selected-sample metadata
# ------------------------------------------------

sample_info2 <- sample_info[
    match(
        metadata_sub$sanger_sample_id,
        sample_info$sanger_sample_id
    ),
    ,
    drop = FALSE
]

# Add selected individual/sample metadata with "selected_" prefix
add_cols <- c(
    "individual_id",
    "biopsy_type",
    "disease_status",
    "age_binned_10",
    "sex",
    "vcf_sample_order"
)

for (x in add_cols) {
    metadata_sub[[x]] <- sample_info2[[x]]
}


# ------------------------------------------------
# 8. Create Seurat object
# ------------------------------------------------

obj_Healthy <- CreateSeuratObject(
    counts = counts_sub,
    meta.data = metadata_sub,
    project = "WellcomeSanger_Healthy"
)

obj_Healthy

ncol(obj_Healthy)
length(unique(obj_Healthy$sanger_sample_id))
length(unique(obj_Healthy$individual_id))


# ------------------------------------------------
# 9. Save
# ------------------------------------------------

saveRDS(
    obj_Healthy,
    out_rds,
    compress = FALSE
)

cat("Saved:", out_rds, "\n")