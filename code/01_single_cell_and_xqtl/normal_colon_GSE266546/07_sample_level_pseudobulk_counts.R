#!/usr/bin/env Rscript

# 07_sample_level_pseudobulk_counts.R
# ------------------------------------------------------------------------------
# Goal:
#   Create sample-level pseudobulk raw gene-count matrices for each reviewed
#   cell type.
#
# For each cell type:
#   - rows    = genes
#   - columns = SampleID
#   - values  = SUM of raw RNA counts across all cells from that sample
#
# Example outputs:
#   pseudobulk_counts_ABS.csv
#   pseudobulk_counts_GOB.csv
#   pseudobulk_counts_STM.csv
#   pseudobulk_counts_EE.csv
#   pseudobulk_counts_CT.csv
#   pseudobulk_counts_TAC.csv
#   pseudobulk_counts_TUF.csv
#
# Notes:
#   - Uses RAW RNA counts, not SCT values.
#   - No normalization is performed here.
#   - "Review" and NA annotations are excluded by default.
# ------------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(dplyr)
})

set.seed(123)

# ---- paths -------------------------------------------------------------------
RESULTS_DIR <- "/data/l2_bioinfo1/liq17/scTWAS/data/VUMC_HTAN_GSE266546/raw/rds2"

INPUT_RDS <- file.path(
  RESULTS_DIR,
  "integrated_seurat_final2_reviewed.res0.8.rm23_sub.rds"
)

OUT_DIR <- file.path(
  RESULTS_DIR,
  "pseudobulk_counts_by_celltype"
)

dir.create(
  OUT_DIR,
  showWarnings = FALSE,
  recursive = TRUE
)

stopifnot(file.exists(INPUT_RDS))

# ---- load object -------------------------------------------------------------
message("Loading: ", INPUT_RDS)
obj <- readRDS(INPUT_RDS)

message(
  "Object: ",
  nrow(obj), " genes x ",
  ncol(obj), " cells"
)

# ---- required metadata -------------------------------------------------------
required_meta <- c(
  "SampleID",
  "celltype_annot_reviewed"
)

missing_meta <- setdiff(
  required_meta,
  colnames(obj@meta.data)
)

if (length(missing_meta) > 0) {
  stop(
    "Missing required metadata columns: ",
    paste(missing_meta, collapse = ", ")
  )
}

# ---- use raw RNA counts ------------------------------------------------------
if (!"RNA" %in% Assays(obj)) {
  stop("RNA assay not found in Seurat object.")
}

DefaultAssay(obj) <- "RNA"

# Seurat v5 objects may contain multiple RNA count layers.
# JoinLayers() ensures GetAssayData(layer = "counts") returns one count matrix.
if (inherits(obj[["RNA"]], "Assay5")) {
  message("Joining RNA layers...")
  obj <- JoinLayers(
    obj,
    assay = "RNA"
  )
}

counts <- GetAssayData(
  obj,
  assay = "RNA",
  layer = "counts"
)

message(
  "RNA count matrix: ",
  nrow(counts), " genes x ",
  ncol(counts), " cells"
)

# ---- metadata aligned to count matrix ----------------------------------------
meta <- obj@meta.data[
  colnames(counts),
  ,
  drop = FALSE
]

meta$SampleID <- as.character(meta$SampleID)
meta$celltype_annot_reviewed <- as.character(
  meta$celltype_annot_reviewed
)

# Exclude cells without a final usable annotation.
keep <- !is.na(meta$SampleID) &
        meta$SampleID != "" &
        !is.na(meta$celltype_annot_reviewed) &
        meta$celltype_annot_reviewed != "" &
        meta$celltype_annot_reviewed != "Review"

counts <- counts[, keep, drop = FALSE]
meta   <- meta[keep, , drop = FALSE]

cell_types <- sort(
  unique(meta$celltype_annot_reviewed)
)

message(
  "Cell types to pseudobulk: ",
  paste(cell_types, collapse = ", ")
)

# ==============================================================================
# Aggregate raw counts: gene x cell -> gene x sample, separately by cell type
# ==============================================================================
for (ct in cell_types) {

  message("\nProcessing cell type: ", ct)

  ct_cells <- rownames(meta)[
    meta$celltype_annot_reviewed == ct
  ]

  if (length(ct_cells) == 0) {
    warning("No cells found for cell type: ", ct)
    next
  }

  ct_meta <- meta[
    ct_cells,
    ,
    drop = FALSE
  ]

  ct_counts <- counts[
    ,
    ct_cells,
    drop = FALSE
  ]

  samples <- sort(
    unique(ct_meta$SampleID)
  )

  message(
    "  ",
    length(ct_cells), " cells across ",
    length(samples), " samples"
  )

  # Initialize genes x samples pseudobulk matrix.
  pb <- matrix(
    0,
    nrow = nrow(ct_counts),
    ncol = length(samples),
    dimnames = list(
      rownames(ct_counts),
      samples
    )
  )

  # Sum raw counts across cells for each sample.
  for (sid in samples) {

    sid_cells <- rownames(ct_meta)[
      ct_meta$SampleID == sid
    ]

    pb[, sid] <- Matrix::rowSums(
      ct_counts[
        ,
        sid_cells,
        drop = FALSE
      ]
    )
  }

  # Store gene names explicitly as first column for CSV.
  pb_df <- data.frame(
    gene = rownames(pb),
    pb,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )

  # Make filename safe.
  safe_ct <- gsub(
    "[^A-Za-z0-9_-]",
    "_",
    ct
  )

  out_file <- file.path(
    OUT_DIR,
    paste0(
      "pseudobulk_counts_",
      safe_ct,
      ".csv"
    )
  )

  write.csv(
    pb_df,
    out_file,
    row.names = FALSE,
    quote = FALSE
  )

  message(
    "  Wrote: ",
    out_file
  )

  message(
    "  Matrix: ",
    nrow(pb), " genes x ",
    ncol(pb), " samples"
  )

  rm(
    ct_counts,
    ct_meta,
    pb,
    pb_df
  )

  gc()
}

message(
  "\nDone. Pseudobulk count matrices are in:\n",
  OUT_DIR
)
