#!/usr/bin/env Rscript

# =============================================================================
# Explore metadata from three VUMC HTAN h5ad files
#
# Files:
#   1. VUMC_HTAN_DIS_EPI_V2.h5ad
#   2. VUMC_HTAN_VAL_EPI_V2.h5ad
#   3. PCA83samples.h5ad
#
# Uses the same load_one_h5ad() strategy as the main CellOracleR pipeline:
#   - ad.raw.X if available, otherwise ad.X
#   - exports ad.obs metadata
#   - cell QC: >= 200 detected genes
#   - standardize_htan_meta_df()
#
# Outputs:
#   - full metadata for each h5ad
#   - metadata-column summary
#   - raw/original metadata-column comparison
#   - standardized metadata-column comparison
#   - categorical-value counts
#   - important metadata distributions
#   - key cross-tabulations
# =============================================================================


# -----------------------------------------------------------------------------
# 0. Packages / environment
# -----------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(Matrix)
  library(dplyr)
  library(data.table)
})

RENV_ACTIVATE <- "/data/l2_bioinfo1/liq17/scTWAS/src/renv/activate.R"

if (file.exists(RENV_ACTIVATE)) {
  source(RENV_ACTIVATE)
}

options(stringsAsFactors = FALSE)

set.seed(202609)


# -----------------------------------------------------------------------------
# 1. Paths
# -----------------------------------------------------------------------------

ROOT <- "/data/l2_bioinfo1/liq17/scTWAS/data/VUMC_HTAN/VUMC_HTAN_DIS_VAL_ADSER"

RESULTS_DIR <- file.path(
  ROOT,
  "Results"
)

HTAN_DIR <- "/data/l2_bioinfo1/scRNA-seq/VUMC_HTAN/data/"

META_OUT_DIR <- file.path(
  RESULTS_DIR,
  "h5ad_metadata_exploration"
)

dir.create(
  RESULTS_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  META_OUT_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)


message("ROOT         = ", ROOT)
message("HTAN_DIR     = ", HTAN_DIR)
message("RESULTS_DIR  = ", RESULTS_DIR)
message("META_OUT_DIR = ", META_OUT_DIR)


# -----------------------------------------------------------------------------
# 2. H5AD files
# -----------------------------------------------------------------------------

h5ad_files <- file.path(
  HTAN_DIR,
  c(
    "VUMC_HTAN_DIS_EPI_V2.h5ad",   # DIS batch1
    "VUMC_HTAN_VAL_EPI_V2.h5ad",   # VAL batch1
    "PCA83samples.h5ad"             # batch2
  )
)


missing_h5ad <- h5ad_files[
  !file.exists(h5ad_files)
]

if (length(missing_h5ad) > 0) {

  stop(
    "Missing required h5ad file(s):\n  ",
    paste(
      missing_h5ad,
      collapse = "\n  "
    )
  )
}


message("\nUsing h5ad files:")

message(
  paste0(
    "  ",
    basename(h5ad_files),
    collapse = "\n"
  )
)


# =============================================================================
# 3. Metadata helper functions
# =============================================================================


# -----------------------------------------------------------------------------
# Find the first available metadata column among candidate names
# -----------------------------------------------------------------------------

first_meta_col <- function(md, candidates) {

  hit <- intersect(
    candidates,
    colnames(md)
  )[1]

  if (is.na(hit)) {
    return(NULL)
  }

  hit
}


# -----------------------------------------------------------------------------
# Standardize HTAN metadata
#
# This follows the logic in your main CellOracleR script.
# Original metadata columns are retained.
# New standardized columns are added.
# -----------------------------------------------------------------------------

standardize_htan_meta_df <- function(md, nm) {

  md <- as.data.frame(
    md,
    stringsAsFactors = FALSE
  )


  # ---------------------------------------------------------------------------
  # Dataset identifiers
  # ---------------------------------------------------------------------------

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


  # ---------------------------------------------------------------------------
  # sample
  # ---------------------------------------------------------------------------

  if (!"sample" %in% colnames(md)) {
    md$sample <- nm
  }

  md$sample <- as.character(md$sample)


  # ---------------------------------------------------------------------------
  # Patient ID
  # ---------------------------------------------------------------------------

  if ("STUDYID.htan" %in% colnames(md)) {

    md$PatientID <- as.character(
      md[["STUDYID.htan"]]
    )

  } else if ("HTAN Parent Data File ID" %in% colnames(md)) {

    md$PatientID <- sub(
      "^(([^_]+_[^_]+)).*$",
      "\\1",
      as.character(
        md[["HTAN Parent Data File ID"]]
      )
    )

  } else if ("HTAN Specimen ID" %in% colnames(md)) {

    md$PatientID <- sub(
      "^(([^_]+_[^_]+)).*$",
      "\\1",
      as.character(
        md[["HTAN Specimen ID"]]
      )
    )

  } else {

    md$PatientID <- nm
  }


  # ---------------------------------------------------------------------------
  # Sample ID
  # ---------------------------------------------------------------------------

  if ("SPECID.htan" %in% colnames(md)) {

    md$SampleID <- gsub(
      "^(HTA[0-9]+_[0-9]+_[0-9]+)_([0-9]+)$",
      "\\1\\2",
      as.character(
        md[["SPECID.htan"]]
      )
    )

  } else if ("HTAN Specimen ID" %in% colnames(md)) {

    md$SampleID <- as.character(
      md[["HTAN Specimen ID"]]
    )

  } else {

    md$SampleID <- NA_character_
  }


  # ---------------------------------------------------------------------------
  # Cell type
  # ---------------------------------------------------------------------------

  md$CellType_cb <- NA_character_


  if ("Cell_Type" %in% colnames(md)) {

    idx <- (
      is.na(md$CellType_cb) &
        !is.na(md[["Cell_Type"]])
    )

    md$CellType_cb[idx] <-
      as.character(
        md[["Cell_Type"]][idx]
      )
  }


  if ("Cell Type" %in% colnames(md)) {

    idx <- (
      is.na(md$CellType_cb) &
        !is.na(md[["Cell Type"]])
    )

    md$CellType_cb[idx] <-
      as.character(
        md[["Cell Type"]][idx]
      )
  }


  # ---------------------------------------------------------------------------
  # PCA83 pathology fields
  # ---------------------------------------------------------------------------

  if ("ADTYPE" %in% colnames(md)) {

    md$PCA83_ADTYPE <- as.character(
      md[["ADTYPE"]]
    )

  } else {

    md$PCA83_ADTYPE <- NA_character_
  }


  if ("HPTYPE" %in% colnames(md)) {

    md$PCA83_HPTYPE <- as.character(
      md[["HPTYPE"]]
    )

  } else {

    md$PCA83_HPTYPE <- NA_character_
  }


  if ("TISSUE" %in% colnames(md)) {

    md$TISSUE_original <- as.character(
      md[["TISSUE"]]
    )

  } else {

    md$TISSUE_original <- NA_character_
  }


  # ---------------------------------------------------------------------------
  # Polyp / tumor type
  # ---------------------------------------------------------------------------

  type_col <- first_meta_col(
    md,
    c(
      "HPTYPE",
      "Tumor_Type",
      "Polyp_Type",
      "Type"
    )
  )


  if (!is.null(type_col)) {

    md$Type <- as.character(
      md[[type_col]]
    )

  } else {

    md$Type <- NA_character_
  }


  # ---------------------------------------------------------------------------
  # Disease / sample condition
  # ---------------------------------------------------------------------------

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

    md$sample_condition <- as.character(
      md[[cond_col]]
    )

    md$sample_condition_source <- cond_col

  } else {

    md$sample_condition <- NA_character_
    md$sample_condition_source <- NA_character_
  }


  # ---------------------------------------------------------------------------
  # Force standardized identifiers to character
  # ---------------------------------------------------------------------------

  for (
    id_col in intersect(
      c(
        "sample",
        "PatientID",
        "SampleID",
        "CellType_cb",
        "Type",
        "PCA83_ADTYPE",
        "PCA83_HPTYPE",
        "TISSUE_original",
        "sample_condition",
        "sample_condition_source",
        "study",
        "htan_file",
        "cohort",
        "batch"
      ),
      colnames(md)
    )
  ) {

    md[[id_col]] <-
      as.character(
        md[[id_col]]
      )
  }


  md
}


# =============================================================================
# 4. Detect Python environment
# =============================================================================

detect_h5ad_python <- function() {

  env_py <- Sys.getenv(
    "HTAN_H5AD_PYTHON",
    unset = ""
  )


  candidates <- unique(
    c(
      env_py,

      "/home/lyul1/venvs/r-h5ad/bin/python",

      "/mnt_home/venvs/r-h5ad/bin/python",

      "/home/liq17/venvs/r-h5ad/bin/python",

      Sys.which("python3"),

      Sys.which("python")
    )
  )


  candidates <- candidates[
    nzchar(candidates)
  ]


  import_test <- paste(
    "import anndata",
    "import pandas",
    "import scipy",
    "import scipy.sparse",
    "import scipy.io",
    sep = "; "
  )


  message(
    "\nSearching for a Python environment that can read h5ad..."
  )


  checked <- character(0)


  for (py in candidates) {

    if (!file.exists(py)) {
      next
    }


    checked <- c(
      checked,
      py
    )


    status <- suppressWarnings(
      system2(
        py,
        c(
          "-c",
          shQuote(import_test)
        ),
        stdout = FALSE,
        stderr = FALSE
      )
    )


    if (identical(status, 0L)) {

      message(
        "  Using Python: ",
        py
      )

      return(py)

    } else {

      message(
        "  Rejected Python: ",
        py
      )
    }
  }


  stop(
    paste0(

      "No usable Python interpreter was found.\n",

      "Required Python packages:\n",

      "  anndata\n",
      "  pandas\n",
      "  scipy\n",
      "  scipy.sparse\n",
      "  scipy.io\n\n",

      if (length(checked) > 0) {

        paste0(
          "Python executables checked:\n  ",
          paste(
            checked,
            collapse = "\n  "
          ),
          "\n\n"
        )

      } else {

        "No candidate Python executable existed.\n\n"
      },

      "You can specify a Python environment with:\n",

      "  export HTAN_H5AD_PYTHON=/full/path/to/python\n"
    )
  )
}


H5AD_PYTHON <- detect_h5ad_python()


# =============================================================================
# 5. load_one_h5ad()
# =============================================================================

load_one_h5ad <- function(f) {

  nm <- tools::file_path_sans_ext(
    basename(f)
  )


  py <- H5AD_PYTHON


  outdir <- tempfile(
    paste0(
      "h5ad_",
      nm,
      "_"
    )
  )

  dir.create(
    outdir,
    recursive = TRUE
  )


  script <- tempfile(
    fileext = ".py"
  )


  # ---------------------------------------------------------------------------
  # Python: extract matrix + metadata
  # ---------------------------------------------------------------------------

  writeLines(
    c(

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

      "# AnnData = cells x genes; Seurat/R matrix = genes x cells",

      "X = X.T",

      "",

      "if not sp.issparse(X):",

      "    X = sp.csr_matrix(X)",

      "",

      "sio.mmwrite(outdir + '/matrix.mtx', X)",

      "",

      "pd.Series(genes).to_csv(",

      "    outdir + '/features.tsv',",

      "    sep='\\t',",

      "    index=False,",

      "    header=False",

      ")",

      "",

      "pd.Series(ad.obs_names).to_csv(",

      "    outdir + '/barcodes.tsv',",

      "    sep='\\t',",

      "    index=False,",

      "    header=False",

      ")",

      "",

      "# Original ad.obs metadata",

      "ad.obs.to_csv(",

      "    outdir + '/obs.tsv',",

      "    sep='\\t'",

      ")",

      "",

      "# Save original ad.obs column names separately",

      "pd.Series(ad.obs.columns).to_csv(",

      "    outdir + '/obs_columns.tsv',",

      "    sep='\\t',",

      "    index=False,",

      "    header=False",

      ")",

      "",

      "with open(outdir + '/matrix_source.txt', 'w') as fh:",

      "    fh.write(matrix_source + '\\n')"

    ),
    script
  )


  status <- system2(
    py,
    c(
      script,
      f,
      outdir
    )
  )


  if (status != 0) {

    stop(
      "Python h5ad conversion failed for ",
      basename(f),
      "\nPython used: ",
      py
    )
  }


  # ---------------------------------------------------------------------------
  # Read expression matrix
  # ---------------------------------------------------------------------------

  mat <- Matrix::readMM(
    file.path(
      outdir,
      "matrix.mtx"
    )
  )


  genes <- read.delim(
    file.path(
      outdir,
      "features.tsv"
    ),
    header = FALSE,
    stringsAsFactors = FALSE
  )[[1]]


  barcodes <- read.delim(
    file.path(
      outdir,
      "barcodes.tsv"
    ),
    header = FALSE,
    stringsAsFactors = FALSE
  )[[1]]


  # ---------------------------------------------------------------------------
  # Read original ad.obs metadata
  # ---------------------------------------------------------------------------

  obs <- read.delim(
    file.path(
      outdir,
      "obs.tsv"
    ),
    row.names = 1,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )


  original_metadata_columns <- read.delim(
    file.path(
      outdir,
      "obs_columns.tsv"
    ),
    header = FALSE,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )[[1]]


  matrix_source <- readLines(
    file.path(
      outdir,
      "matrix_source.txt"
    ),
    n = 1,
    warn = FALSE
  )


  # ---------------------------------------------------------------------------
  # Matrix names
  # ---------------------------------------------------------------------------

  rownames(mat) <-
    make.unique(
      as.character(genes)
    )

  colnames(mat) <-
    as.character(barcodes)


  mat <- as(
    mat,
    "CsparseMatrix"
  )


  # ---------------------------------------------------------------------------
  # Make sure metadata and matrix align
  # ---------------------------------------------------------------------------

  if (!all(
    colnames(mat) %in% rownames(obs)
  )) {

    missing_obs <- setdiff(
      colnames(mat),
      rownames(obs)
    )

    stop(
      "Metadata missing for ",
      length(missing_obs),
      " cell(s) in ",
      basename(f)
    )
  }


  obs <- obs[
    colnames(mat),
    ,
    drop = FALSE
  ]


  n_cells_before_qc <- ncol(mat)
  n_genes_before_qc <- nrow(mat)


  message("\n============================================================")

  message(
    "Loaded: ",
    basename(f)
  )

  message(
    "  matrix source: ",
    matrix_source
  )

  message(
    "  before cell QC: ",
    nrow(mat),
    " genes x ",
    ncol(mat),
    " cells"
  )

  message(
    "  original ad.obs columns: ",
    length(original_metadata_columns)
  )


  # ---------------------------------------------------------------------------
  # Cell QC: >=200 detected genes
  # ---------------------------------------------------------------------------

  nFeature <- Matrix::colSums(
    mat > 0
  )


  keep_cell <- (
    nFeature >= 200
  )


  message(
    "  cell QC >=200 detected genes: retained ",
    sum(keep_cell),
    " / ",
    length(keep_cell),
    " cells"
  )


  mat <- mat[
    ,
    keep_cell,
    drop = FALSE
  ]


  obs <- obs[
    keep_cell,
    ,
    drop = FALSE
  ]


  # Keep QC value in metadata for exploration
  obs$nFeature_h5ad <- nFeature[
    keep_cell
  ]


  # ---------------------------------------------------------------------------
  # Standardize metadata
  # ---------------------------------------------------------------------------

  obs <- standardize_htan_meta_df(
    obs,
    nm
  )


  # ---------------------------------------------------------------------------
  # Cell IDs
  # ---------------------------------------------------------------------------

  obs$h5ad_barcode <- rownames(obs)


  merged_cell_id <- paste0(
    nm,
    "_",
    obs$h5ad_barcode
  )


  colnames(mat) <- merged_cell_id
  rownames(obs) <- merged_cell_id


  stopifnot(

    ncol(mat) == nrow(obs),

    all(
      colnames(mat) ==
        rownames(obs)
    )
  )


  message(
    "  after cell QC: ",
    nrow(mat),
    " genes x ",
    ncol(mat),
    " cells"
  )


  message(
    "  standardized metadata columns: ",
    ncol(obs)
  )


  message(
    "============================================================"
  )


  # ---------------------------------------------------------------------------
  # Clean temporary files
  # ---------------------------------------------------------------------------

  unlink(
    outdir,
    recursive = TRUE,
    force = TRUE
  )


  unlink(
    script,
    force = TRUE
  )


  # ---------------------------------------------------------------------------
  # Return
  # ---------------------------------------------------------------------------

  list(

    counts = mat,

    meta = obs,

    name = nm,

    matrix_source = matrix_source,

    original_metadata_columns =
      original_metadata_columns,

    n_cells_before_qc =
      n_cells_before_qc,

    n_cells_after_qc =
      ncol(mat),

    n_genes =
      n_genes_before_qc
  )
}


# =============================================================================
# 6. Load all three h5ad files
# =============================================================================

message("\n")
message("############################################################")
message("LOAD THREE H5AD FILES")
message("############################################################")


h5ad_objects <- lapply(
  h5ad_files,
  load_one_h5ad
)


names(h5ad_objects) <- vapply(
  h5ad_objects,
  `[[`,
  character(1),
  "name"
)


message(
  "\nLoaded objects: ",
  paste(
    names(h5ad_objects),
    collapse = ", "
  )
)


# =============================================================================
# 7. Basic file summary
# =============================================================================

file_summary <- bind_rows(
  lapply(
    h5ad_objects,
    function(x) {

      data.frame(

        h5ad_file = x$name,

        matrix_source =
          x$matrix_source,

        genes =
          x$n_genes,

        cells_before_QC =
          x$n_cells_before_qc,

        cells_after_QC =
          x$n_cells_after_qc,

        original_metadata_columns =
          length(
            x$original_metadata_columns
          ),

        standardized_metadata_columns =
          ncol(x$meta),

        stringsAsFactors = FALSE
      )
    }
  )
)


message("\n============================================================")
message("FILE SUMMARY")
message("============================================================")

print(
  file_summary,
  row.names = FALSE
)


data.table::fwrite(
  file_summary,
  file.path(
    META_OUT_DIR,
    "01_h5ad_file_summary.tsv"
  ),
  sep = "\t"
)


# =============================================================================
# 8. Extract metadata
# =============================================================================

metadata_list <- lapply(
  h5ad_objects,
  function(x) {
    x$meta
  }
)


# =============================================================================
# 9. Print metadata columns and first rows
# =============================================================================

for (nm in names(metadata_list)) {

  md <- metadata_list[[nm]]

  original_cols <-
    h5ad_objects[[nm]]$original_metadata_columns


  message("\n\n")
  message("############################################################")
  message("FILE: ", nm)
  message("############################################################")


  message(
    "\nCells after QC: ",
    nrow(md)
  )


  message(
    "Metadata columns after standardization: ",
    ncol(md)
  )


  message(
    "\n----------------------------------------"
  )

  message(
    "ORIGINAL ad.obs COLUMNS"
  )

  message(
    "----------------------------------------"
  )


  print(
    original_cols
  )


  added_cols <- setdiff(
    colnames(md),
    original_cols
  )


  message(
    "\n----------------------------------------"
  )

  message(
    "COLUMNS ADDED BY R PIPELINE"
  )

  message(
    "----------------------------------------"
  )


  print(
    added_cols
  )


  message(
    "\n----------------------------------------"
  )

  message(
    "FIRST 6 METADATA ROWS"
  )

  message(
    "----------------------------------------"
  )


  print(
    head(md)
  )
}


# =============================================================================
# 10. Metadata summary function
# =============================================================================

summarize_metadata <- function(
  md,
  file_name,
  original_cols
) {

  out <- lapply(
    colnames(md),
    function(col) {

      x <- md[[col]]

      x_char <- as.character(x)


      missing_idx <- (
        is.na(x_char) |
          !nzchar(
            trimws(x_char)
          )
      )


      valid_values <- x_char[
        !missing_idx
      ]


      unique_values <- unique(
        valid_values
      )


      example_values <- head(
        unique_values,
        10
      )


      data.frame(

        h5ad_file =
          file_name,

        metadata_column =
          col,

        column_origin =
          ifelse(
            col %in% original_cols,
            "original_ad.obs",
            "added_by_R_pipeline"
          ),

        class =
          paste(
            class(x),
            collapse = ";"
          ),

        n_cells =
          length(x),

        n_nonmissing =
          sum(!missing_idx),

        n_missing =
          sum(missing_idx),

        percent_missing =
          round(
            100 *
              sum(missing_idx) /
              length(x),
            3
          ),

        n_unique =
          length(unique_values),

        example_values =
          paste(
            example_values,
            collapse = " | "
          ),

        stringsAsFactors = FALSE
      )
    }
  )


  bind_rows(out)
}


# =============================================================================
# 11. Generate metadata summary for all three files
# =============================================================================

metadata_summary <- bind_rows(
  lapply(
    names(metadata_list),
    function(nm) {

      summarize_metadata(

        md =
          metadata_list[[nm]],

        file_name =
          nm,

        original_cols =
          h5ad_objects[[nm]]$
            original_metadata_columns
      )
    }
  )
)


message("\n============================================================")
message("METADATA SUMMARY")
message("============================================================")


message(
  "Metadata summary contains ",
  nrow(metadata_summary),
  " rows. Showing the first 30 rows in the console."
)

print(
  head(as.data.frame(metadata_summary), 30),
  row.names = FALSE
)


data.table::fwrite(
  metadata_summary,
  file.path(
    META_OUT_DIR,
    "02_all_metadata_column_summary.tsv"
  ),
  sep = "\t"
)


# =============================================================================
# 12. Save per-file metadata summaries
# =============================================================================

for (nm in names(metadata_list)) {

  temp <- metadata_summary[
    metadata_summary$h5ad_file == nm,
    ,
    drop = FALSE
  ]


  data.table::fwrite(
    temp,
    file.path(
      META_OUT_DIR,
      paste0(
        nm,
        "_metadata_column_summary.tsv"
      )
    ),
    sep = "\t"
  )
}


# =============================================================================
# 13. Compare ORIGINAL ad.obs columns across three files
# =============================================================================

all_original_cols <- sort(
  unique(
    unlist(
      lapply(
        h5ad_objects,
        function(x) {
          x$original_metadata_columns
        }
      ),
      use.names = FALSE
    )
  )
)


original_column_presence <- data.frame(
  metadata_column =
    all_original_cols,
  stringsAsFactors = FALSE
)


for (nm in names(h5ad_objects)) {

  original_column_presence[[nm]] <-
    ifelse(
      all_original_cols %in%
        h5ad_objects[[nm]]$
          original_metadata_columns,
      "Present",
      "Absent"
    )
}


message("\n============================================================")
message("ORIGINAL ad.obs COLUMN COMPARISON")
message("============================================================")


print(
  original_column_presence,
  row.names = FALSE
)


data.table::fwrite(
  original_column_presence,
  file.path(
    META_OUT_DIR,
    "03_original_ad.obs_column_presence.tsv"
  ),
  sep = "\t"
)


# =============================================================================
# 14. Compare ALL standardized metadata columns
# =============================================================================

all_standardized_cols <- sort(
  unique(
    unlist(
      lapply(
        metadata_list,
        colnames
      ),
      use.names = FALSE
    )
  )
)


standardized_column_presence <- data.frame(

  metadata_column =
    all_standardized_cols,

  stringsAsFactors = FALSE
)


for (nm in names(metadata_list)) {

  standardized_column_presence[[nm]] <-
    ifelse(
      all_standardized_cols %in%
        colnames(
          metadata_list[[nm]]
        ),
      "Present",
      "Absent"
    )
}


data.table::fwrite(
  standardized_column_presence,
  file.path(
    META_OUT_DIR,
    "04_standardized_metadata_column_presence.tsv"
  ),
  sep = "\t"
)


# =============================================================================
# 15. Categorical metadata counts
#
# Only metadata columns with <= MAX_CATEGORICAL_LEVELS distinct values
# are included.
# =============================================================================

MAX_CATEGORICAL_LEVELS <- 50


categorical_counts_list <- list()

counter <- 1L


for (nm in names(metadata_list)) {

  md <- metadata_list[[nm]]


  for (col in colnames(md)) {

    x <- as.character(
      md[[col]]
    )


    x[
      is.na(x) |
        !nzchar(
          trimws(x)
        )
    ] <- "<NA>"


    n_unique <- length(
      unique(x)
    )


    if (
      n_unique >= 1 &&
        n_unique <= MAX_CATEGORICAL_LEVELS
    ) {

      tab <- as.data.frame(
        table(
          x,
          useNA = "ifany"
        ),
        stringsAsFactors = FALSE
      )


      colnames(tab) <- c(
        "value",
        "n_cells"
      )


      tab$h5ad_file <-
        nm


      tab$metadata_column <-
        col


      tab <- tab[
        order(
          -tab$n_cells
        ),
        ,
        drop = FALSE
      ]


      categorical_counts_list[[counter]] <-
        tab[
          ,
          c(
            "h5ad_file",
            "metadata_column",
            "value",
            "n_cells"
          )
        ]


      counter <- counter + 1L
    }
  }
}


categorical_counts <- bind_rows(
  categorical_counts_list
)


data.table::fwrite(
  categorical_counts,
  file.path(
    META_OUT_DIR,
    "05_categorical_metadata_counts.tsv"
  ),
  sep = "\t"
)


# =============================================================================
# 16. Important metadata fields
# =============================================================================

important_cols <- c(

  "study",

  "htan_file",

  "cohort",

  "batch",

  "sample",

  "PatientID",

  "SampleID",

  "STUDYID.htan",

  "SPECID.htan",

  "HTAN Parent Data File ID",

  "HTAN Specimen ID",

  "Cell_Type",

  "Cell Type",

  "CellType_cb",

  "Sample_Classification",

  "sample_condition",

  "sample_condition_source",

  "TISSUE",

  "TISSUE_original",

  "ADTYPE",

  "HPTYPE",

  "PCA83_ADTYPE",

  "PCA83_HPTYPE",

  "Tumor_Type",

  "Polyp_Type",

  "Type"
)


important_counts <- categorical_counts[
  categorical_counts$metadata_column %in%
    important_cols,
  ,
  drop = FALSE
]


data.table::fwrite(
  important_counts,
  file.path(
    META_OUT_DIR,
    "06_important_metadata_counts.tsv"
  ),
  sep = "\t"
)


# -----------------------------------------------------------------------------
# Print important metadata distributions
# -----------------------------------------------------------------------------

for (nm in names(metadata_list)) {

  md <- metadata_list[[nm]]


  message("\n\n")
  message("############################################################")
  message("IMPORTANT METADATA: ", nm)
  message("############################################################")


  available <- intersect(
    important_cols,
    colnames(md)
  )


  for (col in available) {

    message(
      "\n--- ",
      col,
      " ---"
    )


    print(
      sort(
        table(
          md[[col]],
          useNA = "ifany"
        ),
        decreasing = TRUE
      )
    )
  }
}


# =============================================================================
# 17. Key cross-tabulations
# =============================================================================

cross_tab_dir <- file.path(
  META_OUT_DIR,
  "cross_tabs"
)


dir.create(
  cross_tab_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


for (nm in names(metadata_list)) {

  md <- metadata_list[[nm]]


  # ---------------------------------------------------------------------------
  # sample_condition x CellType_cb
  # ---------------------------------------------------------------------------

  if (
    all(
      c(
        "sample_condition",
        "CellType_cb"
      ) %in% colnames(md)
    )
  ) {

    tab <- as.data.frame.matrix(
      table(
        md$sample_condition,
        md$CellType_cb,
        useNA = "ifany"
      )
    )


    tab <-
      cbind(
        sample_condition =
          rownames(tab),
        tab
      )


    rownames(tab) <- NULL


    data.table::fwrite(
      tab,
      file.path(
        cross_tab_dir,
        paste0(
          nm,
          "_sampleCondition_x_CellType.tsv"
        )
      ),
      sep = "\t"
    )
  }


  # ---------------------------------------------------------------------------
  # Patient x Sample
  # ---------------------------------------------------------------------------

  if (
    all(
      c(
        "PatientID",
        "SampleID"
      ) %in% colnames(md)
    )
  ) {

    patient_sample <- md %>%
      mutate(
        PatientID =
          as.character(PatientID),

        SampleID =
          as.character(SampleID)
      ) %>%
      group_by(
        PatientID,
        SampleID
      ) %>%
      summarise(
        cells = n(),
        .groups = "drop"
      ) %>%
      arrange(
        PatientID,
        SampleID
      )


    data.table::fwrite(
      patient_sample,
      file.path(
        cross_tab_dir,
        paste0(
          nm,
          "_PatientID_x_SampleID.tsv"
        )
      ),
      sep = "\t"
    )
  }


  # ---------------------------------------------------------------------------
  # PCA83: ADTYPE x HPTYPE
  # ---------------------------------------------------------------------------

  if (
    all(
      c(
        "PCA83_ADTYPE",
        "PCA83_HPTYPE"
      ) %in% colnames(md)
    ) &&
      any(
        !is.na(
          md$PCA83_ADTYPE
        ) |
          !is.na(
            md$PCA83_HPTYPE
          )
      )
  ) {

    pca83_pathology <- as.data.frame(
      table(
        ADTYPE =
          md$PCA83_ADTYPE,

        HPTYPE =
          md$PCA83_HPTYPE,

        useNA =
          "ifany"
      ),
      stringsAsFactors = FALSE
    )


    data.table::fwrite(
      pca83_pathology,
      file.path(
        cross_tab_dir,
        paste0(
          nm,
          "_ADTYPE_x_HPTYPE.tsv"
        )
      ),
      sep = "\t"
    )
  }


  # ---------------------------------------------------------------------------
  # PCA83: TISSUE x CellType
  # ---------------------------------------------------------------------------

  if (
    all(
      c(
        "TISSUE_original",
        "CellType_cb"
      ) %in% colnames(md)
    ) &&
      any(
        !is.na(
          md$TISSUE_original
        )
      )
  ) {

    tissue_celltype <- as.data.frame(
      table(
        TISSUE =
          md$TISSUE_original,

        CellType =
          md$CellType_cb,

        useNA =
          "ifany"
      ),
      stringsAsFactors = FALSE
    )


    data.table::fwrite(
      tissue_celltype,
      file.path(
        cross_tab_dir,
        paste0(
          nm,
          "_TISSUE_x_CellType.tsv"
        )
      ),
      sep = "\t"
    )
  }
}


# =============================================================================
# 18. Sample-level overview
# =============================================================================

sample_overview_list <- list()


for (nm in names(metadata_list)) {

  md <- metadata_list[[nm]]


  # We use standardized SampleID here.
  # If SampleID is unavailable, use "sample".
  sample_col <- if (
    "SampleID" %in% colnames(md) &&
      any(
        !is.na(md$SampleID)
      )
  ) {

    "SampleID"

  } else {

    "sample"
  }


  md$.__sample_id__ <-
    as.character(
      md[[sample_col]]
    )


  sample_overview <- md %>%
    group_by(
      .__sample_id__
    ) %>%
    summarise(

      h5ad_file =
        nm,

      sample_column_used =
        sample_col,

      cells =
        n(),

      patients =
        if (
          "PatientID" %in%
            colnames(md)
        ) {
          n_distinct(
            PatientID,
            na.rm = TRUE
          )
        } else {
          NA_integer_
        },

      cell_types =
        if (
          "CellType_cb" %in%
            colnames(md)
        ) {
          paste(
            sort(
              unique(
                na.omit(
                  as.character(
                    CellType_cb
                  )
                )
              )
            ),
            collapse = " | "
          )
        } else {
          NA_character_
        },

      sample_conditions =
        if (
          "sample_condition" %in%
            colnames(md)
        ) {
          paste(
            sort(
              unique(
                na.omit(
                  as.character(
                    sample_condition
                  )
                )
              )
            ),
            collapse = " | "
          )
        } else {
          NA_character_
        },

      tissue =
        if (
          "TISSUE_original" %in%
            colnames(md)
        ) {
          paste(
            sort(
              unique(
                na.omit(
                  as.character(
                    TISSUE_original
                  )
                )
              )
            ),
            collapse = " | "
          )
        } else {
          NA_character_
        },

      ADTYPE =
        if (
          "PCA83_ADTYPE" %in%
            colnames(md)
        ) {
          paste(
            sort(
              unique(
                na.omit(
                  as.character(
                    PCA83_ADTYPE
                  )
                )
              )
            ),
            collapse = " | "
          )
        } else {
          NA_character_
        },

      HPTYPE =
        if (
          "PCA83_HPTYPE" %in%
            colnames(md)
        ) {
          paste(
            sort(
              unique(
                na.omit(
                  as.character(
                    PCA83_HPTYPE
                  )
                )
              )
            ),
            collapse = " | "
          )
        } else {
          NA_character_
        },

      .groups = "drop"
    )


  colnames(sample_overview)[
    colnames(sample_overview) ==
      ".__sample_id__"
  ] <- "sample_id"


  sample_overview_list[[nm]] <-
    sample_overview
}


sample_overview_all <- bind_rows(
  sample_overview_list
)


data.table::fwrite(
  sample_overview_all,
  file.path(
    META_OUT_DIR,
    "07_sample_level_overview.tsv"
  ),
  sep = "\t"
)


# =============================================================================
# 19. Save full metadata from each h5ad
#
# Use .tsv.gz because cell-level metadata may be large.
# =============================================================================

for (nm in names(metadata_list)) {

  md <- metadata_list[[nm]]


  md_export <- data.frame(

    cell_id =
      rownames(md),

    md,

    check.names =
      FALSE,

    stringsAsFactors =
      FALSE
  )


  outfile <- file.path(
    META_OUT_DIR,
    paste0(
      nm,
      "_metadata_afterCell200_standardized.tsv.gz"
    )
  )


  data.table::fwrite(
    md_export,
    outfile,
    sep = "\t"
  )


  message(
    "Saved full metadata: ",
    outfile
  )
}


# =============================================================================
# 20. Original metadata columns only
# =============================================================================

for (nm in names(metadata_list)) {

  md <- metadata_list[[nm]]


  original_cols <-
    intersect(
      h5ad_objects[[nm]]$
        original_metadata_columns,
      colnames(md)
    )


  md_original <- data.frame(

    cell_id =
      rownames(md),

    md[
      ,
      original_cols,
      drop = FALSE
    ],

    check.names = FALSE,

    stringsAsFactors = FALSE
  )


  data.table::fwrite(
    md_original,
    file.path(
      META_OUT_DIR,
      paste0(
        nm,
        "_ORIGINAL_ad.obs_afterCell200.tsv.gz"
      )
    ),
    sep = "\t"
  )
}


# =============================================================================
# 21. Console summary of biologically important fields
# =============================================================================

message("\n\n")
message("############################################################")
message("FINAL METADATA OVERVIEW")
message("############################################################")


for (nm in names(metadata_list)) {

  md <- metadata_list[[nm]]


  message("\n============================================================")
  message(nm)
  message("============================================================")


  message(
    "Cells: ",
    format(
      nrow(md),
      big.mark = ","
    )
  )


  if ("PatientID" %in% colnames(md)) {

    message(
      "Patients: ",
      dplyr::n_distinct(
        md$PatientID,
        na.rm = TRUE
      )
    )
  }


  if ("SampleID" %in% colnames(md)) {

    message(
      "Samples: ",
      dplyr::n_distinct(
        md$SampleID,
        na.rm = TRUE
      )
    )
  }


  if ("CellType_cb" %in% colnames(md)) {

    message("\nCellType_cb:")

    print(
      sort(
        table(
          md$CellType_cb,
          useNA = "ifany"
        ),
        decreasing = TRUE
      )
    )
  }


  if ("sample_condition" %in% colnames(md)) {

    message("\nsample_condition:")

    print(
      sort(
        table(
          md$sample_condition,
          useNA = "ifany"
        ),
        decreasing = TRUE
      )
    )
  }


  if ("TISSUE_original" %in% colnames(md)) {

    if (
      any(
        !is.na(
          md$TISSUE_original
        )
      )
    ) {

      message("\nTISSUE_original:")

      print(
        sort(
          table(
            md$TISSUE_original,
            useNA = "ifany"
          ),
          decreasing = TRUE
        )
      )
    }
  }


  if ("PCA83_ADTYPE" %in% colnames(md)) {

    if (
      any(
        !is.na(
          md$PCA83_ADTYPE
        )
      )
    ) {

      message("\nPCA83_ADTYPE:")

      print(
        sort(
          table(
            md$PCA83_ADTYPE,
            useNA = "ifany"
          ),
          decreasing = TRUE
        )
      )
    }
  }


  if ("PCA83_HPTYPE" %in% colnames(md)) {

    if (
      any(
        !is.na(
          md$PCA83_HPTYPE
        )
      )
    ) {

      message("\nPCA83_HPTYPE:")

      print(
        sort(
          table(
            md$PCA83_HPTYPE,
            useNA = "ifany"
          ),
          decreasing = TRUE
        )
      )
    }
  }
}


# =============================================================================
# 22. Free count matrices
#
# Metadata exploration no longer needs the large expression matrices.
# =============================================================================

for (i in seq_along(h5ad_objects)) {

  h5ad_objects[[i]]$counts <- NULL
}


gc(
  verbose = FALSE
)


# =============================================================================
# 23. Final output information
# =============================================================================

message("\n\n")
message("############################################################")
message("METADATA EXPLORATION COMPLETE")
message("############################################################")


message(
  "\nOutput directory:\n  ",
  META_OUT_DIR
)


message(
  "\nMain output files:"
)


message(
  "  01_h5ad_file_summary.tsv"
)


message(
  "  02_all_metadata_column_summary.tsv"
)


message(
  "  03_original_ad.obs_column_presence.tsv"
)


message(
  "  04_standardized_metadata_column_presence.tsv"
)


message(
  "  05_categorical_metadata_counts.tsv"
)


message(
  "  06_important_metadata_counts.tsv"
)


message(
  "  07_sample_level_overview.tsv"
)


message(
  "\nPer-file outputs:"
)


for (nm in names(metadata_list)) {

  message(
    "  ",
    nm,
    "_metadata_afterCell200_standardized.tsv.gz"
  )

  message(
    "  ",
    nm,
    "_ORIGINAL_ad.obs_afterCell200.tsv.gz"
  )

  message(
    "  ",
    nm,
    "_metadata_column_summary.tsv"
  )
}


message(
  "\nCross-tabs directory:\n  ",
  cross_tab_dir
)


message("\nDone.")