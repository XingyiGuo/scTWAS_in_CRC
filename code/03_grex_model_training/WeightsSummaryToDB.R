library("data.table")
library("dplyr")
library("DBI")
library("RSQLite")

# Revision: robust named-column DB construction.
# - data.table outputs from read_strict_fread() are converted to data.frame
#   before character-vector column selection
# - model summaries are validated/selected by column name, not numeric index
# - R2_global and r_squared_global are retained in the DB extra table
# - weights/covariances are normalized using explicit column names
# - DB schema is checked after writing

args <- commandArgs(trailingOnly = TRUE)

# =============================================================================
# I. Load configuration
# Usage:
#   Rscript WeightsSummaryToDB.R 1 config/HTAN_GSE266546_peakSNP.ADSER.config.txt
# =============================================================================
get_script_dir <- function() {
  full_args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", full_args, value = TRUE)
  if (length(file_arg) == 0) return(getwd())
  dirname(normalizePath(sub("^--file=", "", file_arg[1]), mustWork = FALSE))
}

fill_template <- function(template, values) {
  out <- template
  for (nm in names(values)) {
    out <- gsub(paste0("{", nm, "}"), as.character(values[[nm]]), out, fixed = TRUE)
  }
  out
}

read_strict_fread <- function(file_path, expected_cols) {
  if (!file.exists(file_path)) {
    cat("Warning: File does not exist ->", file_path, "\n")
    return(NULL)
  }
  lines <- readLines(file_path)
  split_lines <- strsplit(gsub("\\s+", "\t", lines), "\t")
  valid_lines <- split_lines[sapply(split_lines, length) == expected_cols]
  if (length(valid_lines) == 0) {
    cat("Warning: No valid", expected_cols, "column rows in", file_path, "\n")
    return(NULL)
  }
  if (length(valid_lines) == 1L) return(NULL)

  texts <- vapply(valid_lines, function(x) paste(x, collapse = "\t"), character(1))
  fread(text = paste(texts, collapse = "\n"), header = TRUE)
}

if (length(args) < 1) {
  stop("Usage: Rscript WeightsSummaryToDB.R <CT_index> [*.config.txt]")
}

if (length(args) >= 2) {
  config_file <- args[2]
} else {
  script_dir <- get_script_dir()
  config_candidates <- unique(c(
    file.path(script_dir, c("TWAS.config.txt", "TWAS.config")),
    list.files(script_dir, pattern = "\\.config\\.txt$", full.names = TRUE)
  ))
  config_candidates <- config_candidates[file.exists(config_candidates)]

  if (length(config_candidates) == 0) {
    stop("No *.config.txt/TWAS.config file found beside the script; pass the config as argument 2.")
  }
  if (length(config_candidates) > 1) {
    stop(
      "Multiple config files found beside the script; pass the intended config as argument 2: ",
      paste(basename(config_candidates), collapse = ", ")
    )
  }
  config_file <- config_candidates[1]
}
if (!file.exists(config_file)) stop("TWAS config file not found: ", config_file)
source(config_file)
if (!exists("TWAS_CONFIG")) stop("TWAS_CONFIG was not defined in: ", config_file)

tissue <- TWAS_CONFIG$tissue
if (!tissue %in% names(TWAS_CONFIG$tissues)) {
  stop("Configured tissue is unavailable: ", tissue,
       ". Available: ", paste(names(TWAS_CONFIG$tissues), collapse = ", "))
}

tissue_cfg <- TWAS_CONFIG$tissues[[tissue]]
CTlist <- tissue_cfg$cell_types
ct_index <- suppressWarnings(as.integer(args[1]))
if (is.na(ct_index) || ct_index < 1 || ct_index > length(CTlist)) {
  stop("CT_index must be between 1 and ", length(CTlist), " for tissue ", tissue)
}

CT <- CTlist[ct_index]
values <- list(
  work_dir = tissue_cfg$work_dir,
  dataset_prefix = TWAS_CONFIG$dataset_prefix,
  tissue = tissue,
  CT = CT,
  chrom = ""
)
model_prefix <- fill_template(TWAS_CONFIG$model_prefix_template, values)
output_dir <- fill_template(TWAS_CONFIG$paths$output_dir_template, values)
chromosomes <- as.character(TWAS_CONFIG$chromosomes)

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

cat("Config       :", normalizePath(config_file, mustWork = FALSE), "\n")
cat("Tissue       :", tissue, "\n")
cat("Cell type    :", CT, "\n")
cat("Model prefix :", model_prefix, "\n")
cat("Input/output :", output_dir, "\n")

# =============================================================================
# II. Collect model summary files
# =============================================================================
model_summary <- NULL
model_weights <- NULL
model_covars <- NULL

# Required model-summary fields are selected by NAME rather than column position.
# This protects the DB build if TWAS_modelbuilding changes column order.
required_summary_cols <- c(
  "gene",
  "gene_name",
  "pred_perf_R2",
  "pred_perf_pval",
  "n_snps_in_model",
  "R2_global",
  "r_squared_global"
)

for (chrom in chromosomes) {
  f <- file.path(output_dir, paste0(model_prefix, "_chr", chrom, "_model_summaries.txt"))

  if (file.exists(f)) {
    model_summary_perchr <- as.data.frame(
      fread(f, check.names = FALSE),
      stringsAsFactors = FALSE
    )

    missing_summary_cols <- setdiff(
      required_summary_cols,
      names(model_summary_perchr)
    )

    if (length(missing_summary_cols) > 0) {
      stop(
        "Model summary for chr", chrom,
        " is missing required column(s): ",
        paste(missing_summary_cols, collapse = ", "),
        "\nFile: ", f,
        "\nObserved columns: ",
        paste(names(model_summary_perchr), collapse = ", ")
      )
    }

    # Select by explicit names only. Extra columns are allowed and ignored.
    model_summary_perchr_remain <- model_summary_perchr[
      ,
      required_summary_cols,
      drop = FALSE
    ]

    model_summary <- rbind(
      model_summary,
      model_summary_perchr_remain
    )
  } else {
    cat("Warning: missing model summary ->", f, "\n")
  }
}

if (is.null(model_summary) || nrow(model_summary) == 0) {
  stop("No model summary rows were collected from: ", output_dir)
}

# Normalize types before filtering/deduplication.
numeric_summary_cols <- c(
  "pred_perf_R2",
  "pred_perf_pval",
  "n_snps_in_model",
  "R2_global",
  "r_squared_global"
)

for (nm in numeric_summary_cols) {
  model_summary[[nm]] <- suppressWarnings(
    as.numeric(model_summary[[nm]])
  )
}

model_summary$gene <- as.character(model_summary$gene)
model_summary$gene_name <- as.character(model_summary$gene_name)

# Retain models with a finite r_squared_global, matching the previous
# requirement that the model have a valid global correlation statistic.
model_summary <- model_summary[
  is.finite(model_summary$r_squared_global),
  ,
  drop = FALSE
]

if (nrow(model_summary) == 0) {
  stop(
    "No model-summary rows remain after requiring finite r_squared_global."
  )
}

# Preserve the previous duplicate-resolution behavior, but reference the
# statistic by name rather than through the former r.global positional alias.
# Duplicates are not expected in the normal one-gene/one-chromosome workflow.
model_summary_unique <- model_summary %>%
  arrange(r_squared_global) %>%
  distinct(gene, .keep_all = TRUE)

# Build the SPrediXcan-compatible "extra" table.
#
# Keep the legacy field names used by the current MetaXcan workflow, and add
# the two exact model-performance fields requested here.
#
# NOTE: the source chromosome summary does not contain a pred_perf_qval field.
# For backward compatibility with the previous DB builder, pred.perf.qval is
# retained as a copy of pred_perf_pval. This preserves existing behavior rather
# than silently changing downstream SPrediXcan output.
model_summary_unique <- model_summary_unique %>%
  transmute(
    gene = gene,
    genename = gene_name,
    pred.perf.r2 = pred_perf_R2,
    pred.perf.pval = pred_perf_pval,
    pred.perf.qval = pred_perf_pval,
    n.snps.in.model = n_snps_in_model,

    # Exact fields requested from *_chrX_model_summaries.txt
    R2_global = R2_global,
    r_squared_global = r_squared_global,

    # Backward-compatible alias retained from the previous DB schema.
    r.global = r_squared_global
  )

summary_csv <- file.path(
  output_dir,
  paste0(model_prefix, "_model_summaries.csv")
)

# =============================================================================
# III. Collect model weights
# =============================================================================
for (chrom in chromosomes) {
  f <- file.path(output_dir, paste0(model_prefix, "_chr", chrom, "_weights.txt"))
  if (file.exists(f)) {
    model_weights_perchr <- read_strict_fread(f, 6)
    if (!is.null(model_weights_perchr)) {
      # read_strict_fread() returns a data.table. Convert to an ordinary
      # data.frame so character-vector column selection behaves consistently
      # with base R and does not require data.table's .. prefix / with=FALSE.
      model_weights_perchr <- as.data.frame(
        model_weights_perchr,
        stringsAsFactors = FALSE
      )

      required_weight_cols <- c(
        "gene_id", "rsid", "varID", "ref", "effect", "beta"
      )
      missing_weight_cols <- setdiff(
        required_weight_cols,
        names(model_weights_perchr)
      )

      if (length(missing_weight_cols) > 0) {
        stop(
          "Model weights for chr", chrom,
          " are missing required column(s): ",
          paste(missing_weight_cols, collapse = ", "),
          "\nFile: ", f,
          "\nObserved columns: ",
          paste(names(model_weights_perchr), collapse = ", ")
        )
      }

      # Select/reorder by name. Input column order may change safely.
      model_weights_perchr_remain <- model_weights_perchr[
        ,
        c("rsid", "gene_id", "beta", "ref", "effect"),
        drop = FALSE
      ]
      model_weights <- rbind(model_weights, model_weights_perchr_remain)
    }
  } else {
    cat("Warning: missing model weights ->", f, "\n")
  }
}

if (is.null(model_weights) || nrow(model_weights) == 0) {
  stop("No model weight rows were collected from: ", output_dir)
}
model_weights <- model_weights %>%
  transmute(
    rsid = as.character(rsid),
    gene = as.character(gene_id),
    weight = suppressWarnings(as.numeric(beta)),
    ref_allele = toupper(as.character(ref)),
    eff_allele = toupper(as.character(effect))
  )

model_weights <- model_weights[!is.na(model_weights$weight), ]

# SPrediXcan interprets weight relative to eff_allele.  In this pipeline,
# eff_allele is exactly the "effect" allele whose dosage entered elastic net.
bad_allele <- is.na(model_weights$ref_allele) | is.na(model_weights$eff_allele) |
              model_weights$ref_allele == "" | model_weights$eff_allele == "" |
              toupper(model_weights$ref_allele) == toupper(model_weights$eff_allele)
if (any(bad_allele)) {
  stop("Invalid ref/effect allele values in ", sum(bad_allele), " model-weight row(s).")
}

model_weights_unique <- model_weights %>% distinct(gene, rsid, .keep_all = TRUE)

# =============================================================================
# IV. Collect model covariances
# =============================================================================
for (chrom in chromosomes) {
  f <- file.path(output_dir, paste0(model_prefix, "_chr", chrom, "_covariances.txt"))
  if (file.exists(f)) {
    model_covars_perchr <- read_strict_fread(f, 4)
    if (!is.null(model_covars_perchr)) {
      # Same reason as for weights above: normalize data.table -> data.frame
      # before selecting columns through a character vector.
      model_covars_perchr <- as.data.frame(
        model_covars_perchr,
        stringsAsFactors = FALSE
      )

      required_cov_cols <- c(
        "gene_id", "rsid1", "rsid2", "corvarianceValues"
      )
      missing_cov_cols <- setdiff(
        required_cov_cols,
        names(model_covars_perchr)
      )

      if (length(missing_cov_cols) > 0) {
        stop(
          "Model covariances for chr", chrom,
          " are missing required column(s): ",
          paste(missing_cov_cols, collapse = ", "),
          "\nFile: ", f,
          "\nObserved columns: ",
          paste(names(model_covars_perchr), collapse = ", ")
        )
      }

      # Select/reorder by name. Input column order may change safely.
      model_covars_perchr <- model_covars_perchr[
        ,
        required_cov_cols,
        drop = FALSE
      ]
      model_covars <- rbind(model_covars, model_covars_perchr)
    }
  } else {
    cat("Warning: missing model covariances ->", f, "\n")
  }
}

if (is.null(model_covars) || nrow(model_covars) == 0) {
  stop("No model covariance rows were collected from: ", output_dir)
}
model_covars <- model_covars %>%
  transmute(
    GENE = as.character(gene_id),
    RSID1 = as.character(rsid1),
    RSID2 = as.character(rsid2),
    VALUE = suppressWarnings(as.numeric(corvarianceValues))
  )

model_covars <- model_covars[!is.na(model_covars$VALUE), ]
model_covars_unique <- model_covars %>% distinct(GENE, RSID1, RSID2, .keep_all = TRUE)

cov_file <- file.path(output_dir, paste0(model_prefix, "_cov.txt.gz"))
fwrite(model_covars_unique, cov_file, sep = " ", compress = "gzip")
fwrite(model_summary_unique, summary_csv)

# =============================================================================
# V. Generate SPrediXcan SQLite DB
# =============================================================================
dbfile <- file.path(output_dir, paste0(model_prefix, ".db"))

# Replace both tables together so a failed rerun preserves the existing model DB.
db <- dbConnect(SQLite(), dbname = dbfile)
extra_fields <- tryCatch(
  dbWithTransaction(db, {
    dbWriteTable(db, "extra", model_summary_unique, overwrite = TRUE, row.names = FALSE)
    dbWriteTable(db, "weights", model_weights_unique, overwrite = TRUE, row.names = FALSE)

    fields <- dbListFields(db, "extra")
    missing_fields <- setdiff(c("R2_global", "r_squared_global"), fields)
    if (length(missing_fields) > 0) {
      stop("DB extra table is missing field(s): ", paste(missing_fields, collapse = ", "))
    }
    fields
  }),
  finally = dbDisconnect(db)
)

cat("DB extra fields:", paste(extra_fields, collapse = ", "), "\n")
cat("Wrote summary :", summary_csv, "\n")
cat("Wrote cov     :", cov_file, "\n")
cat("Wrote DB      :", dbfile, "\n")
