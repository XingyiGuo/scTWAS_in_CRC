#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(dplyr)
  library(data.table)
})

set.seed(123)

# =============================================================================
# HTAN_ADSER_CellOracleRsave_mergeQC_PCA83_HPtoSER_ASCSSCfilter.R
# =============================================================================
# Purpose
#   Prepare Vanderbilt HTAN epithelial DIS + VAL + PCA83 data for downstream
#   CellOracleR analyses.
#
# Input h5ad files
#   1) VUMC_HTAN_DIS_EPI_V2.h5ad
#   2) VUMC_HTAN_VAL_EPI_V2.h5ad
#   3) PCA83samples.h5ad
#
# Filtering strategy
#   CELL QC (per h5ad, before merge):
#     - keep cells with >= 200 detected genes
#
#   GENE QC (after DIS + VAL + PCA83 are merged):
#     - keep genes detected in >= 3 cells across the COMBINED dataset
#
#   CONDITION + CELL-TYPE FILTER (after merge/QC):
#     - DIS/VAL: keep sample_condition in c("AD", "SER") and the requested 9 epithelial cell types
#     - PCA83: keep TISSUE == "Polyp" and classify pathology hierarchically:
#         AD  = ADTYPE Tubular, Tubulovillous
#         SER = ADTYPE Sessile Serrated, Traditional Serrated
#         SER = HPTYPE Microvesicular, Goblet-cell rich ONLY when ADTYPE is
#               missing or Not stated (broad serrated-polyp definition)
#       HPTYPE Unknown or missing without an informative ADTYPE remains
#       unclassified and is excluded. PCA83 classified cells are merged into
#       the same canonical AD/SER groups as DIS/VAL. All non-missing PCA83
#       CellType_cb labels are otherwise allowed.
#     - Lineage/tissue consistency rule applied across ALL cohorts:
#         ASC cells are retained only from AD tissue
#         SSC cells are retained only from SER tissue
#
#   CELL-TYPE OUTPUT:
#     - save one Seurat RDS per retained CellType_cb across all retained cohorts
#     - also save one combined Seurat RDS for audit/reuse
#
# Important
#   - No upper nFeature_RNA/nCount_RNA cutoff is applied here.
#   - No mitochondrial-percent or doublet filtering is added here.
#   - If ad.raw exists, ad.raw.X is used; otherwise ad.X is used.
#   - Python is auto-detected; no single /mnt_home interpreter path is required.
#   - The final Seurat object is intentionally saved as a reusable input for
#     later CellOracleR normalization/subsetting/network analyses.
# =============================================================================

# -----------------------------------------------------------------------------
# Paths
# -----------------------------------------------------------------------------
RENV_ACTIVATE <- "/data/l2_bioinfo1/liq17/scTWAS/src/renv/activate.R"
source(RENV_ACTIVATE)

options(stringsAsFactors = FALSE)
set.seed(202609)

ROOT <- "/data/l2_bioinfo1/liq17/scTWAS/data/VUMC_HTAN/VUMC_HTAN_DIS_VAL_ADSER"
RESULTS_DIR <- file.path(ROOT, "Results")
CELLORACLE_INPUT_DIR <- RESULTS_DIR
HTAN_DIR <- "/data/l2_bioinfo1/scRNA-seq/VUMC_HTAN/data/"

dir.create(RESULTS_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(CELLORACLE_INPUT_DIR, showWarnings = FALSE, recursive = TRUE)

message("ROOT                 = ", ROOT)
message("HTAN_DIR             = ", HTAN_DIR)
message("RESULTS_DIR          = ", RESULTS_DIR)
message("CELLORACLE_INPUT_DIR = ", CELLORACLE_INPUT_DIR)

# -----------------------------------------------------------------------------
# Requested disease conditions and epithelial cell types
# -----------------------------------------------------------------------------
TARGET_CONDITIONS <- c("AD", "SER")

TARGET_CELLTYPES <- c(
  "ABS", "GOB", "CT", "EE", "TUF", "TAC", "STM", "ASC", "SSC"
)

message("TARGET_CONDITIONS = ", paste(TARGET_CONDITIONS, collapse = ", "))
message("TARGET_CELLTYPES (DIS/VAL) = ", paste(TARGET_CELLTYPES, collapse = ", "))
message("PCA83 pathology mapping: Tubular/Tubulovillous -> AD; Sessile/Traditional Serrated -> SER; Microvesicular/Goblet-cell rich HP -> SER when ADTYPE is missing/Not stated")
message("PCA83 rule: broad SER includes classified hyperplastic-polyp subtypes; Unknown/unclassified polyps are excluded")
message("Global lineage rule: ASC -> AD only; SSC -> SER only")

# Explicit whitelist: do not recursively load all h5ad files in HTAN_DIR.
h5ad_files <- file.path(
  HTAN_DIR,
  c("VUMC_HTAN_DIS_EPI_V2.h5ad", # DIS_VAL for batch1
    "VUMC_HTAN_VAL_EPI_V2.h5ad",
	"PCA83samples.h5ad"  # For batch2
  )
)

missing_h5ad <- h5ad_files[!file.exists(h5ad_files)]
if (length(missing_h5ad) > 0) {
  stop(
    "Missing required HTAN h5ad file(s):\n  ",
    paste(missing_h5ad, collapse = "\n  ")
  )
}

message("\nUsing HTAN h5ad files:")
message(paste0("  ", basename(h5ad_files), collapse = "\n"))

# -----------------------------------------------------------------------------
# Helpers
# -----------------------------------------------------------------------------
first_meta_col <- function(md, candidates) {
  hit <- intersect(candidates, colnames(md))[1]
  if (is.na(hit)) return(NULL)
  hit
}

# Canonicalize HTAN disease labels before disease-condition filtering.
canonicalize_sample_condition <- function(x) {
  raw <- trimws(as.character(x))
  key <- toupper(raw)
  key[!nzchar(key)] <- NA_character_

  out <- rep(NA_character_, length(key))

  out[
    !is.na(key) &
      (
        key %in% c("AD", "ADENOMA") |
          grepl("ADENOMA", key, fixed = TRUE)
      )
  ] <- "AD"

  out[
    !is.na(key) &
      (
        key %in% c("SER", "SET", "SERRATED") |
          grepl("SERRAT", key)
      )
  ] <- "SER"

  out
}

# Map PCA83 polyp pathology to canonical AD/SER groups.
# Priority is ADTYPE. HPTYPE is used only as a fallback when ADTYPE is missing
# or "Not stated". This implements the requested broad serrated-polyp grouping:
#   AD:  Tubular, Tubulovillous
#   SER: Sessile Serrated, Traditional Serrated
#   SER: Microvesicular HP, Goblet-cell-rich HP (fallback only)
# Unknown/missing pathology remains unclassified and is excluded downstream.
canonicalize_pca83_pathology <- function(adtype, hptype) {
  ad <- toupper(trimws(as.character(adtype)))
  hp <- toupper(trimws(as.character(hptype)))

  ad[is.na(ad) | !nzchar(ad) | ad %in% c("NA", "NAN")] <- NA_character_
  hp[is.na(hp) | !nzchar(hp) | hp %in% c("NA", "NAN")] <- NA_character_

  out <- rep(NA_character_, length(ad))

  explicit_ad <- !is.na(ad) & ad %in% c("TUBULAR", "TUBULOVILLOUS")
  explicit_ser <- !is.na(ad) & ad %in% c("SESSILE SERRATED", "TRADITIONAL SERRATED")
  hp_fallback_ok <- is.na(ad) | ad == "NOT STATED"
  hp_ser <- !is.na(hp) & hp %in% c("MICROVESICULAR", "GOBLET-CELL RICH")

  out[explicit_ad] <- "AD"
  out[explicit_ser] <- "SER"
  out[is.na(out) & hp_fallback_ok & hp_ser] <- "SER"

  out
}

# Standardize HTAN metadata before creating the final merged Seurat object.
standardize_htan_meta_df <- function(md, nm) {
  md <- as.data.frame(md, stringsAsFactors = FALSE)

  md$study <- "HTAN_islam2024"
  md$htan_file <- nm

  if (grepl("_DIS_", nm, fixed = TRUE)) {
    md$cohort <- "DIS"
    md$batch <- "HTAN_DIS"
  } else if (grepl("_VAL_", nm, fixed = TRUE)) {
    md$cohort <- "VAL"
    md$batch <- "HTAN_VAL"
  } else if (identical(nm, "PCA83samples")) {
    md$cohort <- "PCA83"
    md$batch <- "HTAN_PCA83"
  } else {
    md$cohort <- nm
    md$batch <- "HTAN_VUMC"
  }

  if (!"sample" %in% colnames(md)) {
    md$sample <- nm
  }
  # h5ad files may encode sample IDs as integer in one file and character in
  # another. Always treat identifiers as character so metadata can be merged.
  md$sample <- as.character(md$sample)

  # Patient ID
  if ("STUDYID.htan" %in% colnames(md)) {
    md$PatientID <- as.character(md[["STUDYID.htan"]])
  } else if ("HTAN Parent Data File ID" %in% colnames(md)) {
    md$PatientID <- sub(
      "^(([^_]+_[^_]+)).*$",
      "\\1",
      as.character(md[["HTAN Parent Data File ID"]])
    )
  } else if ("HTAN Specimen ID" %in% colnames(md)) {
    md$PatientID <- sub(
      "^(([^_]+_[^_]+)).*$",
      "\\1",
      as.character(md[["HTAN Specimen ID"]])
    )
  } else {
    md$PatientID <- nm
  }

  # Sample ID
  if ("SPECID.htan" %in% colnames(md)) {
    md$SampleID <- gsub(
      "^(HTA[0-9]+_[0-9]+_[0-9]+)_([0-9]+)$",
      "\\1\\2",
      as.character(md[["SPECID.htan"]])
    )
  } else if ("HTAN Specimen ID" %in% colnames(md)) {
    md$SampleID <- as.character(md[["HTAN Specimen ID"]])
  } else {
    md$SampleID <- NA_character_
  }

  # Cell type
  md$CellType_cb <- NA_character_

  if ("Cell_Type" %in% colnames(md)) {
    idx <- is.na(md$CellType_cb) & !is.na(md[["Cell_Type"]])
    md$CellType_cb[idx] <- as.character(md[["Cell_Type"]][idx])
  }

  if ("Cell Type" %in% colnames(md)) {
    idx <- is.na(md$CellType_cb) & !is.na(md[["Cell Type"]])
    md$CellType_cb[idx] <- as.character(md[["Cell Type"]][idx])
  }

  # Preserve PCA83 pathology/tissue fields explicitly for disease mapping and audit.
  if ("ADTYPE" %in% colnames(md)) {
    md$PCA83_ADTYPE <- as.character(md[["ADTYPE"]])
  } else {
    md$PCA83_ADTYPE <- NA_character_
  }

  if ("HPTYPE" %in% colnames(md)) {
    md$PCA83_HPTYPE <- as.character(md[["HPTYPE"]])
  } else {
    md$PCA83_HPTYPE <- NA_character_
  }

  if ("TISSUE" %in% colnames(md)) {
    md$TISSUE_original <- as.character(md[["TISSUE"]])
  } else {
    md$TISSUE_original <- NA_character_
  }

  # Polyp/tumor type
  type_col <- first_meta_col(
    md,
    c("HPTYPE", "Tumor_Type", "Polyp_Type", "Type")
  )

  if (!is.null(type_col)) {
    md$Type <- as.character(md[[type_col]])
  } else {
    md$Type <- NA_character_
  }

  # Disease/sample condition.
  # Prefer Sample_Classification because this is the desired HTAN disease
  # classification field when present. TISSUE is used only as a later fallback.
  cond_col <- first_meta_col(
    md,
    c(
      "Sample_Classification",
      "sample_condition",
      "condition",
      "Status",
      "Histology",
      "TISSUE"
    )
  )

  if (!is.null(cond_col)) {
    md$sample_condition <- as.character(md[[cond_col]])
    md$sample_condition_source <- cond_col
  } else {
    md$sample_condition <- NA_character_
    md$sample_condition_source <- NA_character_
  }

  # Standardized identifiers must remain character across datasets.
  for (id_col in intersect(
    c("sample", "PatientID", "SampleID", "CellType_cb", "Type",
      "PCA83_ADTYPE", "PCA83_HPTYPE", "TISSUE_original",
      "sample_condition", "sample_condition_source", "study", "htan_file",
      "cohort", "batch"),
    colnames(md)
  )) {
    md[[id_col]] <- as.character(md[[id_col]])
  }

  md
}

# Harmonize source metadata columns before bind_rows(). h5ad/TSV import may infer
# the same source field as numeric in one file and character in another. Only
# incompatible shared fields are coerced to character; otherwise native types
# are preserved.
harmonize_meta_for_bind <- function(meta_list) {
  all_cols <- unique(unlist(lapply(meta_list, colnames), use.names = FALSE))

  for (col in all_cols) {
    present <- which(vapply(meta_list, function(md) col %in% colnames(md), logical(1)))
    if (length(present) < 2L) next

    types <- unique(vapply(
      meta_list[present],
      function(md) typeof(md[[col]]),
      character(1)
    ))

    # integer + double are compatible. Any character-vs-numeric/logical
    # mismatch is safest as character for metadata preservation.
    numeric_only <- all(types %in% c("integer", "double"))
    if (length(types) > 1L && !numeric_only) {
      message(
        "  Harmonizing metadata column '", col, "' to character: ",
        paste(types, collapse = " / ")
      )
      for (i in present) meta_list[[i]][[col]] <- as.character(meta_list[[i]][[col]])
    }
  }

  meta_list
}

# Add zero rows so all matrices have the same gene set/order before cbind().
align_to_genes <- function(mat, genes) {
  missing_genes <- setdiff(genes, rownames(mat))

  if (length(missing_genes) > 0) {
    zero_mat <- Matrix::Matrix(
      0,
      nrow = length(missing_genes),
      ncol = ncol(mat),
      sparse = TRUE
    )
    rownames(zero_mat) <- missing_genes
    colnames(zero_mat) <- colnames(mat)
    mat <- rbind(mat, zero_mat)
  }

  mat[genes, , drop = FALSE]
}

# -----------------------------------------------------------------------------
# Detect a usable Python for h5ad reading
# -----------------------------------------------------------------------------
# Priority:
#   1) HTAN_H5AD_PYTHON environment variable, if supplied
#   2) known Vanderbilt/cluster virtualenv locations
#   3) python3 / python available on PATH
#
# A candidate is accepted only if it can import:
#   anndata, pandas, scipy, scipy.sparse, scipy.io
detect_h5ad_python <- function() {
  env_py <- Sys.getenv("HTAN_H5AD_PYTHON", unset = "")

  candidates <- unique(c(
    env_py,
    "/home/lyul1/venvs/r-h5ad/bin/python",
    "/mnt_home/venvs/r-h5ad/bin/python",
    "/home/liq17/venvs/r-h5ad/bin/python",
    Sys.which("python3"),
    Sys.which("python")
  ))

  candidates <- candidates[nzchar(candidates)]

  import_test <- paste(
    "import anndata",
    "import pandas",
    "import scipy",
    "import scipy.sparse",
    "import scipy.io",
    sep = "; "
  )

  message("\nSearching for a Python environment that can read h5ad...")

  checked <- character(0)

  for (py in candidates) {
    # Sys.which() may return a path; explicit candidates may not exist.
    if (!file.exists(py)) next

    checked <- c(checked, py)

    status <- suppressWarnings(
      system2(
        py,
        c("-c", shQuote(import_test)),
        stdout = FALSE,
        stderr = FALSE
      )
    )

    if (identical(status, 0L)) {
      message("  Using Python: ", py)
      return(py)
    } else {
      message("  Rejected Python (missing required package[s]): ", py)
    }
  }

  stop(
    paste0(
      "No usable Python interpreter was found for h5ad reading.\n",
      "A usable interpreter must import: anndata, pandas, scipy.\n\n",
      if (length(checked) > 0) {
        paste0("Python executables checked:\n  ",
               paste(checked, collapse = "\n  "), "\n\n")
      } else {
        "No candidate Python executable existed.\n\n"
      },
      "If you already have a working environment, rerun with:\n",
      "  export HTAN_H5AD_PYTHON=/full/path/to/python\n",
      "  Rscript 01a_HTAN_ADSER_CellOracleRsave_mergeQC_pythonfix.R\n\n",
      "Useful shell checks:\n",
      "  which python3\n",
      "  python3 -c \"import anndata, pandas, scipy; print('OK')\"\n"
    )
  )
}

H5AD_PYTHON <- detect_h5ad_python()

# -----------------------------------------------------------------------------
# Read one h5ad
# -----------------------------------------------------------------------------
# The Python helper exports:
#   - expression matrix
#   - feature names
#   - cell barcodes
#   - ad.obs metadata
#
# Cell QC is applied here: >= 200 detected genes per cell.
# Gene >=3-cells QC is intentionally NOT applied until after DIS+VAL are merged.
load_one_h5ad <- function(f) {
  nm <- tools::file_path_sans_ext(basename(f))

  py <- H5AD_PYTHON

  outdir <- tempfile(paste0("h5ad_", nm, "_"))
  dir.create(outdir, recursive = TRUE)

  script <- tempfile(fileext = ".py")

  writeLines(c(
    "import sys",
    "import anndata",
    "import pandas as pd",
    "import scipy.sparse as sp",
    "import scipy.io as sio",
    "",
    "f = sys.argv[1]",
    "outdir = sys.argv[2]",
    "ad = anndata.read_h5ad(f)",
    "",
    "if ad.raw is not None:",
    "    X = ad.raw.X",
    "    genes = ad.raw.var_names",
    "    matrix_source = 'ad.raw.X'",
    "else:",
    "    X = ad.X",
    "    genes = ad.var_names",
    "    matrix_source = 'ad.X'",
    "",
    "# AnnData is cells x genes; Seurat uses genes x cells",
    "X = X.T",
    "if not sp.issparse(X):",
    "    X = sp.csr_matrix(X)",
    "",
    "sio.mmwrite(outdir + '/matrix.mtx', X)",
    "pd.Series(genes).to_csv(outdir + '/features.tsv', sep='\\t', index=False, header=False)",
    "pd.Series(ad.obs_names).to_csv(outdir + '/barcodes.tsv', sep='\\t', index=False, header=False)",
    "ad.obs.to_csv(outdir + '/obs.tsv', sep='\\t')",
    "with open(outdir + '/matrix_source.txt', 'w') as fh:",
    "    fh.write(matrix_source + '\\n')"
  ), script)

  status <- system2(
    py,
    c(script, f, outdir)
  )

  if (status != 0) {
    stop(
      "Python h5ad conversion failed for ", basename(f),
      "\nPython used: ", py,
      "\nTest the environment directly with:\n  ",
      py, " -c \"import anndata, pandas, scipy; print('OK')\""
    )
  }

  mat <- Matrix::readMM(file.path(outdir, "matrix.mtx"))

  genes <- read.delim(
    file.path(outdir, "features.tsv"),
    header = FALSE,
    stringsAsFactors = FALSE
  )[[1]]

  barcodes <- read.delim(
    file.path(outdir, "barcodes.tsv"),
    header = FALSE,
    stringsAsFactors = FALSE
  )[[1]]

  obs <- read.delim(
    file.path(outdir, "obs.tsv"),
    row.names = 1,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )

  matrix_source <- readLines(
    file.path(outdir, "matrix_source.txt"),
    n = 1,
    warn = FALSE
  )

  rownames(mat) <- make.unique(as.character(genes))
  colnames(mat) <- as.character(barcodes)
  mat <- as(mat, "CsparseMatrix")

  # Ensure obs and matrix contain exactly the same cells in the same order.
  if (!all(colnames(mat) %in% rownames(obs))) {
    missing_obs <- setdiff(colnames(mat), rownames(obs))
    stop(
      "Metadata missing for ", length(missing_obs),
      " cell(s) in ", basename(f)
    )
  }
  obs <- obs[colnames(mat), , drop = FALSE]

  message("\nLoaded: ", basename(f))
  message("  matrix source: ", matrix_source)
  message("  before cell QC: ", nrow(mat), " genes x ", ncol(mat), " cells")

  # ---------------------------------------------------------------------------
  # Cell QC: equivalent to CreateSeuratObject(..., min.features = 200),
  # but performed manually so gene filtering can be delayed until after merge.
  # ---------------------------------------------------------------------------
  nFeature <- Matrix::colSums(mat > 0)
  keep_cell <- nFeature >= 200

  message("  cell QC >=200 detected genes: retained ",
          sum(keep_cell), " / ", length(keep_cell), " cells")

  mat <- mat[, keep_cell, drop = FALSE]
  obs <- obs[keep_cell, , drop = FALSE]

  # Standardize metadata after cell QC.
  obs <- standardize_htan_meta_df(obs, nm)

  # Keep original h5ad barcode and create globally unique merged cell IDs.
  obs$h5ad_barcode <- rownames(obs)
  merged_cell_id <- paste0(nm, "_", obs$h5ad_barcode)

  colnames(mat) <- merged_cell_id
  rownames(obs) <- merged_cell_id

  stopifnot(
    ncol(mat) == nrow(obs),
    all(colnames(mat) == rownames(obs))
  )

  message("  after cell QC : ", nrow(mat), " genes x ", ncol(mat), " cells")
  message("  sample_condition source: ",
          paste(unique(obs$sample_condition_source), collapse = ", "))
  message("  sample_condition:")
  print(table(obs$sample_condition, useNA = "ifany"))
  message("  CellType_cb:")
  print(table(obs$CellType_cb, useNA = "ifany"))

  unlink(outdir, recursive = TRUE, force = TRUE)
  unlink(script, force = TRUE)

  list(
    counts = mat,
    meta = obs,
    name = nm,
    matrix_source = matrix_source
  )
}

# -----------------------------------------------------------------------------
# 1) Load DIS, VAL, and PCA83; apply cell >=200 genes QC separately
# -----------------------------------------------------------------------------
message("\n============================================================")
message("1) LOAD HTAN DIS + VAL + PCA83 EPITHELIAL H5AD")
message("============================================================")

htan_list <- lapply(h5ad_files, load_one_h5ad)
names(htan_list) <- vapply(htan_list, `[[`, character(1), "name")

if (length(htan_list) == 0) {
  stop("No HTAN h5ad objects were successfully loaded.")
}

# -----------------------------------------------------------------------------
# 2) Merge matrices first using the union of genes across all three h5ad files
# -----------------------------------------------------------------------------
message("\n============================================================")
message("2) MERGE DIS + VAL + PCA83")
message("============================================================")

all_genes <- Reduce(
  union,
  lapply(htan_list, function(x) rownames(x$counts))
)

message("Union of genes before combined gene QC: ", length(all_genes))

aligned_mats <- lapply(
  htan_list,
  function(x) align_to_genes(x$counts, all_genes)
)

combined_counts <- do.call(cbind, aligned_mats)

meta_list <- lapply(htan_list, function(x) x$meta)
meta_list <- harmonize_meta_for_bind(meta_list)
combined_meta <- bind_rows(meta_list)
rownames(combined_meta) <- unlist(
  lapply(htan_list, function(x) rownames(x$meta)),
  use.names = FALSE
)

# Reorder metadata to exactly match the combined count matrix.
combined_meta <- combined_meta[colnames(combined_counts), , drop = FALSE]

stopifnot(
  ncol(combined_counts) == nrow(combined_meta),
  all(colnames(combined_counts) == rownames(combined_meta))
)

message("Combined before gene QC: ",
        nrow(combined_counts), " genes x ",
        ncol(combined_counts), " cells")

# -----------------------------------------------------------------------------
# 3) Combined gene QC: gene detected in >=3 cells across DIS + VAL + PCA83
# -----------------------------------------------------------------------------
message("\n============================================================")
message("3) COMBINED GENE QC")
message("============================================================")

gene_ncells <- Matrix::rowSums(combined_counts > 0)
keep_gene <- gene_ncells >= 3

message("Gene QC >=3 cells across combined DIS+VAL+PCA83: retained ",
        sum(keep_gene), " / ", length(keep_gene), " genes")

combined_counts <- combined_counts[keep_gene, , drop = FALSE]

# -----------------------------------------------------------------------------
# 4) Create a single merged Seurat object AFTER the requested QC
# -----------------------------------------------------------------------------
# min.cells/min.features are set to 0 here because both requested filters were
# already applied explicitly above:
#   - cell >=200 genes: before merge
#   - gene detected in >=3 cells: after merge
objH <- CreateSeuratObject(
  counts = combined_counts,
  meta.data = combined_meta,
  project = "HTAN_ADSER_DISVAL_EPI",
  min.cells = 0,
  min.features = 0
)

message("\nMerged Seurat object after cell + gene QC:")
message("  ", nrow(objH), " genes x ", ncol(objH), " cells")
message("  samples : ", length(unique(na.omit(objH$SampleID))))
message("  patients: ", length(unique(na.omit(objH$PatientID))))
message("  cohorts:")
print(table(objH$cohort, useNA = "ifany"))
message("  sample_condition:")
print(table(objH$sample_condition, useNA = "ifany"))
message("  CellType_cb:")
print(table(objH$CellType_cb, useNA = "ifany"))

# Save the merged object BEFORE disease-condition filtering for audit/reuse.
MERGED_QC_RDS <- file.path(
  RESULTS_DIR,
  "objH_DISVAL_PCA83_EPI_afterCell200_Gene3QC_beforeConditionFilter.rds"
)
saveRDS(objH, MERGED_QC_RDS, compress = TRUE)
message("Saved merged QC object before condition filter: ", MERGED_QC_RDS)

# -----------------------------------------------------------------------------
# 5) Cohort-specific disease mapping + lineage/tissue filtering
# -----------------------------------------------------------------------------
message("\n============================================================")
message("4) MAP PCA83 ADTYPE/HPTYPE TO AD/SER + APPLY ASC/SSC TISSUE RULES")
message("============================================================")

# Keep source values for auditing.
objH$sample_condition_original <- as.character(objH$sample_condition)
objH$CellType_cb <- trimws(as.character(objH$CellType_cb))
objH$CellType_cb[!nzchar(objH$CellType_cb)] <- NA_character_

is_pca83 <- as.character(objH$cohort) == "PCA83" |
  as.character(objH$htan_file) == "PCA83samples"
is_disval <- as.character(objH$cohort) %in% c("DIS", "VAL")

# DIS/VAL: canonicalize Sample_Classification to AD/SER.
disval_condition <- canonicalize_sample_condition(objH$sample_condition_original)

# PCA83: require TISSUE == Polyp, then classify hierarchically.
#   ADTYPE Tubular / Tubulovillous                  -> AD
#   ADTYPE Sessile Serrated / Traditional Serrated -> SER
#   If ADTYPE is missing or Not stated:
#       HPTYPE Microvesicular / Goblet-cell rich   -> SER
#   HPTYPE Unknown or otherwise unclassified        -> drop
pca83_tissue <- toupper(trimws(as.character(objH$TISSUE_original)))
pca83_adtype_raw <- trimws(as.character(objH$PCA83_ADTYPE))
pca83_hptype_raw <- trimws(as.character(objH$PCA83_HPTYPE))
pca83_condition <- canonicalize_pca83_pathology(
  pca83_adtype_raw,
  pca83_hptype_raw
)
is_pca83_polyp <- is_pca83 & !is.na(pca83_tissue) & pca83_tissue == "POLYP"

objH$sample_condition <- NA_character_
objH$sample_condition[is_disval] <- disval_condition[is_disval]
objH$sample_condition[is_pca83_polyp] <- pca83_condition[is_pca83_polyp]

# Keep a dedicated source-level field showing how the final AD/SER label arose.
ad_key <- toupper(trimws(pca83_adtype_raw))
hp_key <- toupper(trimws(pca83_hptype_raw))
ad_key[is.na(ad_key) | !nzchar(ad_key) | ad_key %in% c("NA", "NAN")] <- NA_character_
hp_key[is.na(hp_key) | !nzchar(hp_key) | hp_key %in% c("NA", "NAN")] <- NA_character_

pca83_explicit_ad <- is_pca83_polyp & !is.na(ad_key) &
  ad_key %in% c("TUBULAR", "TUBULOVILLOUS")
pca83_explicit_ser <- is_pca83_polyp & !is.na(ad_key) &
  ad_key %in% c("SESSILE SERRATED", "TRADITIONAL SERRATED")
pca83_hp_ser <- is_pca83_polyp & (is.na(ad_key) | ad_key == "NOT STATED") &
  !is.na(hp_key) & hp_key %in% c("MICROVESICULAR", "GOBLET-CELL RICH")

objH$condition_mapping_source <- NA_character_
objH$condition_mapping_source[is_disval & !is.na(objH$sample_condition)] <-
  "DISVAL_Sample_Classification"
objH$condition_mapping_source[pca83_explicit_ad] <-
  "PCA83_ADTYPE_Tubular_or_Tubulovillous"
objH$condition_mapping_source[pca83_explicit_ser] <-
  "PCA83_ADTYPE_Sessile_or_Traditional_Serrated"
objH$condition_mapping_source[pca83_hp_ser] <-
  "PCA83_HPTYPE_Microvesicular_or_GobletCellRich_to_SER"

celltype_upper <- toupper(objH$CellType_cb)
target_celltypes_upper <- toupper(TARGET_CELLTYPES)

# DIS/VAL: canonical AD/SER + requested 9 epithelial cell types.
keep_disval_base <- is_disval &
  !is.na(objH$sample_condition) &
  objH$sample_condition %in% TARGET_CONDITIONS &
  !is.na(celltype_upper) &
  celltype_upper %in% target_celltypes_upper

# PCA83: only pathologically classified Polyp cells; retain every non-missing
# CellType_cb label (including any additional labels outside the standard 9).
keep_pca83_base <- is_pca83_polyp &
  !is.na(objH$sample_condition) &
  objH$sample_condition %in% TARGET_CONDITIONS &
  !is.na(objH$CellType_cb)

keep_before_lineage_rule <- keep_disval_base | keep_pca83_base

# ---------------------------------------------------------------------------
# Global lineage/tissue consistency rule requested by the user:
#   ASC must come from AD tissue only.
#   SSC must come from SER tissue only.
# Applied AFTER disease mapping so it affects DIS, VAL, and PCA83 consistently.
# ---------------------------------------------------------------------------
is_asc <- !is.na(celltype_upper) & celltype_upper == "ASC"
is_ssc <- !is.na(celltype_upper) & celltype_upper == "SSC"

drop_asc_nonad <- keep_before_lineage_rule & is_asc &
  objH$sample_condition != "AD"
drop_ssc_nonser <- keep_before_lineage_rule & is_ssc &
  objH$sample_condition != "SER"

keep_cells <- keep_before_lineage_rule & !drop_asc_nonad & !drop_ssc_nonser

message("\nRaw source condition / tissue values:")
print(table(objH$cohort, objH$sample_condition_original, useNA = "ifany"))

message("\nPCA83 ADTYPE values among Polyp cells:")
print(table(
  objH$PCA83_ADTYPE[is_pca83_polyp],
  useNA = "ifany"
))

message("\nPCA83 HPTYPE values among Polyp cells:")
print(table(
  objH$PCA83_HPTYPE[is_pca83_polyp],
  useNA = "ifany"
))

message("\nPCA83 ADTYPE x HPTYPE -> canonical condition:")
print(table(
  ADTYPE = objH$PCA83_ADTYPE[is_pca83_polyp],
  HPTYPE = objH$PCA83_HPTYPE[is_pca83_polyp],
  assigned_condition = objH$sample_condition[is_pca83_polyp],
  useNA = "ifany"
))

message("\nCanonical conditions by cohort BEFORE ASC/SSC tissue rule:")
print(table(
  objH$cohort[keep_before_lineage_rule],
  objH$sample_condition[keep_before_lineage_rule],
  useNA = "ifany"
))

message("\nFilter counts:")
message("  DIS/VAL AD/SER + 9CT before ASC/SSC rule : ", sum(keep_disval_base))
message("  PCA83 classified Polyp + ALL CT before rule: ", sum(keep_pca83_base))
message("    PCA83 -> AD                              : ",
        sum(keep_pca83_base & objH$sample_condition == "AD"))
message("    PCA83 -> SER                             : ",
        sum(keep_pca83_base & objH$sample_condition == "SER"))
message("  PCA83 Polyp unclassified after ADTYPE + HPTYPE mapping: ",
        sum(is_pca83_polyp & is.na(pca83_condition)))
message("    PCA83 SER assigned by HPTYPE fallback        : ",
        sum(pca83_hp_ser))
message("  Dropped ASC cells from non-AD tissue      : ", sum(drop_asc_nonad))
message("  Dropped SSC cells from non-SER tissue     : ", sum(drop_ssc_nonser))
message("  Final retained cells                      : ", sum(keep_cells))
message("  PCA83 Normal excluded                     : ",
        sum(is_pca83 & !is.na(pca83_tissue) & pca83_tissue == "NORMAL"))
message("  Cells with missing CellType_cb            : ", sum(is.na(objH$CellType_cb)))

# Extra audit by cohort for the lineage-specific exclusions.
lineage_drop_audit <- data.frame(
  cohort = c("DIS", "VAL", "PCA83"),
  dropped_ASC_nonAD = vapply(
    c("DIS", "VAL", "PCA83"),
    function(cc) sum(drop_asc_nonad & as.character(objH$cohort) == cc),
    integer(1)
  ),
  dropped_SSC_nonSER = vapply(
    c("DIS", "VAL", "PCA83"),
    function(cc) sum(drop_ssc_nonser & as.character(objH$cohort) == cc),
    integer(1)
  ),
  stringsAsFactors = FALSE
)
message("\nASC/SSC lineage-rule drop audit by cohort:")
print(lineage_drop_audit)

objH_keep <- subset(
  objH,
  cells = colnames(objH)[keep_cells]
)

if (ncol(objH_keep) == 0) {
  stop(
    "No cells remain after disease mapping and lineage/tissue filtering. ",
    "Check DIS/VAL disease labels and PCA83 TISSUE/ADTYPE/HPTYPE values."
  )
}

# Safety assertions for the requested tissue-specific cell types.
stopifnot(
  all(objH_keep$sample_condition[objH_keep$CellType_cb == "ASC"] == "AD"),
  all(objH_keep$sample_condition[objH_keep$CellType_cb == "SSC"] == "SER")
)

# Standardize capitalization for the known 9 epithelial labels while preserving
# any additional PCA83 labels exactly as supplied (e.g. current label '20').
celltype_lookup <- setNames(TARGET_CELLTYPES, toupper(TARGET_CELLTYPES))
ct_key <- toupper(as.character(objH_keep$CellType_cb))
known_ct <- ct_key %in% names(celltype_lookup)
objH_keep$CellType_cb[known_ct] <- unname(celltype_lookup[ct_key[known_ct]])

# Dynamic output cell-type list: requested 9 first, then any additional PCA83
# labels retained after PCA83 pathology mapping and ASC/SSC disease filtering.
observed_ct <- unique(as.character(objH_keep$CellType_cb))
observed_ct <- observed_ct[!is.na(observed_ct) & nzchar(observed_ct)]
OUTPUT_CELLTYPES <- c(
  TARGET_CELLTYPES[TARGET_CELLTYPES %in% observed_ct],
  sort(setdiff(observed_ct, TARGET_CELLTYPES))
)

message("\nFiltered combined AD/SER epithelial object:")
message("  ", nrow(objH_keep), " genes x ", ncol(objH_keep), " cells")
message("  samples : ", length(unique(na.omit(objH_keep$SampleID))))
message("  patients: ", length(unique(na.omit(objH_keep$PatientID))))
message("  cohorts:")
print(table(objH_keep$cohort, useNA = "ifany"))
message("  sample_condition:")
print(table(objH_keep$sample_condition, useNA = "ifany"))
message("  CellType_cb:")
print(table(objH_keep$CellType_cb, useNA = "ifany"))
message("  cohort x sample_condition:")
print(table(objH_keep$cohort, objH_keep$sample_condition, useNA = "ifany"))
message("  sample_condition x CellType_cb:")
print(table(
  objH_keep$sample_condition,
  objH_keep$CellType_cb,
  useNA = "ifany"
))

# Explicit final assertions after cell-type capitalization.
if (any(objH_keep$CellType_cb == "ASC" & objH_keep$sample_condition != "AD")) {
  stop("ASC tissue rule violated: retained ASC cell(s) outside AD.")
}
if (any(objH_keep$CellType_cb == "SSC" & objH_keep$sample_condition != "SER")) {
  stop("SSC tissue rule violated: retained SSC cell(s) outside SER.")
}

condition_celltype_audit <- as.data.frame(
  table(
    sample_condition = factor(
      objH_keep$sample_condition,
      levels = c("AD", "SER")
    ),
    CellType_cb = factor(
      objH_keep$CellType_cb,
      levels = OUTPUT_CELLTYPES
    )
  ),
  stringsAsFactors = FALSE
)

write.csv(
  condition_celltype_audit,
  file.path(
    CELLORACLE_INPUT_DIR,
    "HTAN_DISVAL_PCA83_ADSER_HPtoSER_AllCT_ASCSSCfiltered_cell_counts_AUDIT.csv"
  ),
  row.names = FALSE
)

write.csv(
  lineage_drop_audit,
  file.path(
    CELLORACLE_INPUT_DIR,
    "HTAN_DISVAL_PCA83_ADSER_HPtoSER_ASCSSC_lineage_drop_AUDIT.csv"
  ),
  row.names = FALSE
)

pca83_pathology_audit <- objH@meta.data[
  is_pca83_polyp,
  c("PCA83_ADTYPE", "PCA83_HPTYPE", "sample_condition", "condition_mapping_source"),
  drop = FALSE
] %>%
  mutate(
    ADTYPE = ifelse(is.na(PCA83_ADTYPE) | !nzchar(trimws(PCA83_ADTYPE)), NA_character_, PCA83_ADTYPE),
    HPTYPE = ifelse(is.na(PCA83_HPTYPE) | !nzchar(trimws(PCA83_HPTYPE)), NA_character_, PCA83_HPTYPE),
    assigned_condition = sample_condition,
    mapping_source = condition_mapping_source
  ) %>%
  count(ADTYPE, HPTYPE, assigned_condition, mapping_source, name = "cells", .drop = FALSE) %>%
  arrange(assigned_condition, ADTYPE, HPTYPE)

write.csv(
  pca83_pathology_audit,
  file.path(
    CELLORACLE_INPUT_DIR,
    "PCA83_Polyp_ADTYPE_HPTYPE_to_ADSER_HPtoSER_AUDIT.csv"
  ),
  row.names = FALSE
)


# -----------------------------------------------------------------------------
# 6) Summary tables
# -----------------------------------------------------------------------------
count_by_group <- function(obj, group_col) {
  obj@meta.data %>%
    group_by(.data[[group_col]]) %>%
    summarise(
      genes = nrow(obj),
      cells = n(),
      patients = n_distinct(PatientID, na.rm = TRUE),
      samples = n_distinct(SampleID, na.rm = TRUE),
      .groups = "drop"
    )
}

summary_by_file <- count_by_group(objH_keep, "htan_file")
summary_by_cohort <- count_by_group(objH_keep, "cohort")
summary_by_condition <- count_by_group(objH_keep, "sample_condition")
summary_by_celltype <- count_by_group(objH_keep, "CellType_cb")

summary_condition_celltype <- objH_keep@meta.data %>%
  group_by(sample_condition, CellType_cb) %>%
  summarise(
    genes = nrow(objH_keep),
    cells = n(),
    patients = n_distinct(PatientID, na.rm = TRUE),
    samples = n_distinct(SampleID, na.rm = TRUE),
    .groups = "drop"
  )

# Dynamic retained cell-type audit, including requested 9 plus any additional classified PCA83 labels.
celltype_audit <- data.frame(
  CellType_cb = OUTPUT_CELLTYPES,
  present = OUTPUT_CELLTYPES %in% unique(as.character(objH_keep$CellType_cb)),
  cells = vapply(
    OUTPUT_CELLTYPES,
    function(ct) {
      sum(as.character(objH_keep$CellType_cb) == ct, na.rm = TRUE)
    },
    integer(1)
  ),
  patients = vapply(
    OUTPUT_CELLTYPES,
    function(ct) {
      md <- objH_keep@meta.data[
        as.character(objH_keep$CellType_cb) == ct,
        ,
        drop = FALSE
      ]
      if (!nrow(md)) 0L else dplyr::n_distinct(md$PatientID, na.rm = TRUE)
    },
    integer(1)
  ),
  samples = vapply(
    OUTPUT_CELLTYPES,
    function(ct) {
      md <- objH_keep@meta.data[
        as.character(objH_keep$CellType_cb) == ct,
        ,
        drop = FALSE
      ]
      if (!nrow(md)) 0L else dplyr::n_distinct(md$SampleID, na.rm = TRUE)
    },
    integer(1)
  ),
  stringsAsFactors = FALSE
)

message("\nSummary by h5ad file:")
print(summary_by_file)

message("\nSummary by cohort:")
print(summary_by_cohort)

message("\nSummary by sample_condition:")
print(summary_by_condition)

message("\nSummary by CellType_cb:")
print(summary_by_celltype)

message("\nSummary by sample_condition x CellType_cb:")
print(summary_condition_celltype)

message("\nRetained cell-type audit:")
print(celltype_audit)


# -----------------------------------------------------------------------------
# 7) Save combined retained object
# -----------------------------------------------------------------------------
message("\n============================================================")
message("5) SAVE COMBINED CELLORACLER AD/SER INPUT")
message("============================================================")

COMBINED_RDS <- file.path(
  CELLORACLE_INPUT_DIR,
  "HTAN_DISVAL_PCA83_ADSER_HPtoSER_AllCT_ASCSSCfiltered_CellOracleR_input.rds"
)

COMBINED_META <- file.path(
  CELLORACLE_INPUT_DIR,
  "HTAN_DISVAL_PCA83_ADSER_HPtoSER_AllCT_ASCSSCfiltered_CellOracleR_metadata.csv"
)

saveRDS(
  objH_keep,
  COMBINED_RDS,
  compress = TRUE
)

write.csv(
  objH_keep@meta.data,
  COMBINED_META,
  row.names = FALSE
)

message("Saved combined CellOracleR object: ", COMBINED_RDS)
message("Saved combined metadata          : ", COMBINED_META)


# -----------------------------------------------------------------------------
# 8) Save one CellOracleR Seurat RDS per cell type
# -----------------------------------------------------------------------------
message("\n============================================================")
message("6) SAVE ONE RDS PER CELL TYPE")
message("============================================================")

PER_CT_DIR <- file.path(
  CELLORACLE_INPUT_DIR,
  "by_celltype_DISVAL_PCA83_ADSER_HPtoSER_ASCSSCfiltered"
)
dir.create(
  PER_CT_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)

per_ct_manifest <- lapply(OUTPUT_CELLTYPES, function(CT) {

  ct_cells <- colnames(objH_keep)[
    as.character(objH_keep$CellType_cb) == CT
  ]

  if (!length(ct_cells)) {
    warning(
      "No retained cells for CellType_cb = ", CT,
      "; no RDS will be written for this cell type."
    )

    return(data.frame(
      CellType_cb = CT,
      present = FALSE,
      cells = 0L,
      AD_cells = 0L,
      SER_cells = 0L,
      AD_patients = 0L,
      SER_patients = 0L,
      patients = 0L,
      samples = 0L,
      rds_file = NA_character_,
      metadata_file = NA_character_,
      stringsAsFactors = FALSE
    ))
  }

  obj_ct <- subset(
    objH_keep,
    cells = ct_cells
  )

  # Safety check: each output must contain exactly one CellType_cb.
  stopifnot(
    ncol(obj_ct) == length(ct_cells),
    length(unique(as.character(obj_ct$CellType_cb))) == 1L,
    unique(as.character(obj_ct$CellType_cb)) == CT
  )

  ct_rds <- file.path(
    PER_CT_DIR,
    paste0(
      "HTAN_DISVAL_PCA83_ADSER_HPtoSER_",
      CT,
      "_CellOracleR_input.rds"
    )
  )

  ct_meta <- file.path(
    PER_CT_DIR,
    paste0(
      "HTAN_DISVAL_PCA83_ADSER_HPtoSER_",
      CT,
      "_CellOracleR_metadata.csv"
    )
  )

  saveRDS(
    obj_ct,
    ct_rds,
    compress = TRUE
  )

  write.csv(
    obj_ct@meta.data,
    ct_meta,
    row.names = FALSE
  )

  message(
    "Saved ", CT, ": ",
    nrow(obj_ct), " genes x ",
    ncol(obj_ct), " cells -> ",
    ct_rds
  )

  ct_cond <- table(
    factor(
      obj_ct$sample_condition,
      levels = c("AD", "SER")
    )
  )

  out <- data.frame(
    CellType_cb = CT,
    present = TRUE,
    cells = ncol(obj_ct),
    AD_cells = as.integer(ct_cond["AD"]),
    SER_cells = as.integer(ct_cond["SER"]),
    AD_patients = dplyr::n_distinct(
      obj_ct$PatientID[obj_ct$sample_condition == "AD"],
      na.rm = TRUE
    ),
    SER_patients = dplyr::n_distinct(
      obj_ct$PatientID[obj_ct$sample_condition == "SER"],
      na.rm = TRUE
    ),
    patients = dplyr::n_distinct(obj_ct$PatientID, na.rm = TRUE),
    samples = dplyr::n_distinct(obj_ct$SampleID, na.rm = TRUE),
    rds_file = ct_rds,
    metadata_file = ct_meta,
    stringsAsFactors = FALSE
  )

  rm(obj_ct)
  gc(verbose = FALSE)

  out
})

per_ct_manifest <- bind_rows(per_ct_manifest)


# -----------------------------------------------------------------------------
# 9) Save summaries and manifest
# -----------------------------------------------------------------------------
SUMMARY_FILE <- file.path(
  CELLORACLE_INPUT_DIR,
  "HTAN_DISVAL_PCA83_ADSER_HPtoSER_AllCT_ASCSSCfiltered_summary_by_file.csv"
)

SUMMARY_COHORT <- file.path(
  CELLORACLE_INPUT_DIR,
  "HTAN_DISVAL_PCA83_ADSER_HPtoSER_AllCT_ASCSSCfiltered_summary_by_cohort.csv"
)

SUMMARY_CONDITION <- file.path(
  CELLORACLE_INPUT_DIR,
  "HTAN_DISVAL_PCA83_ADSER_HPtoSER_AllCT_ASCSSCfiltered_summary_by_condition.csv"
)

SUMMARY_CELLTYPE <- file.path(
  CELLORACLE_INPUT_DIR,
  "HTAN_DISVAL_PCA83_ADSER_HPtoSER_AllCT_ASCSSCfiltered_summary_by_celltype.csv"
)

SUMMARY_CONDITION_CELLTYPE <- file.path(
  CELLORACLE_INPUT_DIR,
  "HTAN_DISVAL_PCA83_ADSER_HPtoSER_AllCT_ASCSSCfiltered_summary_by_condition_celltype.csv"
)

CELLTYPE_AUDIT_FILE <- file.path(
  CELLORACLE_INPUT_DIR,
  "HTAN_DISVAL_PCA83_ADSER_HPtoSER_AllCT_ASCSSCfiltered_presence_audit.csv"
)

PER_CT_MANIFEST <- file.path(
  CELLORACLE_INPUT_DIR,
  "HTAN_DISVAL_PCA83_ADSER_HPtoSER_AllCT_ASCSSCfiltered_perCellType_RDS_manifest.csv"
)

write.csv(summary_by_file, SUMMARY_FILE, row.names = FALSE)
write.csv(summary_by_cohort, SUMMARY_COHORT, row.names = FALSE)
write.csv(summary_by_condition, SUMMARY_CONDITION, row.names = FALSE)
write.csv(summary_by_celltype, SUMMARY_CELLTYPE, row.names = FALSE)
write.csv(
  summary_condition_celltype,
  SUMMARY_CONDITION_CELLTYPE,
  row.names = FALSE
)
write.csv(celltype_audit, CELLTYPE_AUDIT_FILE, row.names = FALSE)
write.csv(per_ct_manifest, PER_CT_MANIFEST, row.names = FALSE)

message("\nSaved summary by file                 : ", SUMMARY_FILE)
message("Saved summary by cohort              : ", SUMMARY_COHORT)
message("Saved summary by condition           : ", SUMMARY_CONDITION)
message("Saved summary by cell type           : ", SUMMARY_CELLTYPE)
message("Saved summary by condition x celltype: ", SUMMARY_CONDITION_CELLTYPE)
message("Saved retained cell-type audit      : ", CELLTYPE_AUDIT_FILE)
message("Saved per-cell-type RDS manifest     : ", PER_CT_MANIFEST)

message("\nFinal CellOracleR preparation:")
message("  combined object: ", nrow(objH_keep),
        " genes x ", ncol(objH_keep), " cells")
message("  cell QC: >=200 detected genes per input h5ad")
message("  gene QC: detected in >=3 cells across merged DIS+VAL+PCA83")
message("  final conditions: AD, SER")
message("  PCA83 AD mapping : Tubular + Tubulovillous -> AD")
message("  PCA83 SER mapping: Sessile Serrated + Traditional Serrated -> SER")
message("  PCA83 HP fallback : Microvesicular + Goblet-cell rich -> SER when ADTYPE is missing/Not stated")
message("  PCA83 HPTYPE Unknown or otherwise unclassified polyps: excluded")
message("  global lineage rule: ASC retained only in AD; SSC retained only in SER")
message("  DIS/VAL CellType_cb kept: ", paste(TARGET_CELLTYPES, collapse = ", "))
message("  PCA83 classified CellType_cb: ALL non-missing; final labels = ",
        paste(OUTPUT_CELLTYPES, collapse = ", "))
message("  per-cell-type RDS directory: ", PER_CT_DIR)

disval_present_ct <- unique(as.character(
  objH_keep$CellType_cb[as.character(objH_keep$cohort) %in% c("DIS", "VAL")]
))
missing_target_ct <- setdiff(TARGET_CELLTYPES, disval_present_ct)

if (length(missing_target_ct)) {
  warning(
    "The following requested DIS/VAL cell types had zero retained cells and therefore no RDS: ",
    paste(missing_target_ct, collapse = ", ")
  )
} else {
  message("All retained cell types have individual RDS files.")
}

rm(
  htan_list,
  aligned_mats,
  combined_counts,
  combined_meta,
  meta_list,
  objH
)
gc()
