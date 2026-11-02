#!/usr/bin/env Rscript

suppressMessages({
  library(data.table)
  library(dplyr)
  library(arrow)
})

# ==============================================================================
# Prepare SMR + HEIDI inputs DIRECTLY from TWAS-significant genes.
#
# Usage:
#   Rscript 02_SMR_HEIDI_inputs_targetTopEQTL_500kb.R ABS \
#     SMR_HEIDI_NL_targetTopEQTL_500kb.config.txt
#
#   Rscript 02_SMR_HEIDI_inputs_targetTopEQTL_500kb.R all \
#     SMR_HEIDI_ADSER_targetTopEQTL_500kb.config.txt
#
# For each TWAS-derived task (gene x cell type x GWAS cohort x eQTL source):
#   1) define TSS +/- flank_length (default 500 kb)
#   2) extract source-matched TensorQTL cis-eQTLs for the gene
#   3) intersect with the corresponding CRC GWAS by rsID
#   4) harmonize GWAS A1/A2 with TensorQTL ALT/REF
#   5) retain usable SNPs in the +/-500-kb region
#   6) choose target SNP = smallest eQTL p-value among harmonized overlaps
#   7) write gene-specific GWAS .ma, eQTL .esd/.flist, SNP list, target QC
#   8) write a manifest consumed by run_SMR_HEIDI_targetTopEQTL_500kb.sh
#
# IMPORTANT:
#   TensorQTL slope and af are assumed ALT-oriented, exactly as in the coloc
#   workflow. For SMR files, eQTL beta/frequency are aligned to GWAS A1.
# ==============================================================================

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2L) {
  stop(
    "Usage: Rscript 02_SMR_HEIDI_inputs_targetTopEQTL_500kb.R ",
    "<CT|all> <SMR_config.R>"
  )
}

CT_ARG <- toupper(trimws(args[1]))
CONFIG_FILE <- args[2]

if (!file.exists(CONFIG_FILE)) stop("Config file not found: ", CONFIG_FILE)
source(CONFIG_FILE)
if (!exists("SMR_CONFIG")) stop("SMR_CONFIG was not defined in: ", CONFIG_FILE)
CFG <- SMR_CONFIG

# ------------------------------------------------------------------------------
# Generic helpers
# ------------------------------------------------------------------------------
fill_template <- function(template, values) {
  out <- template
  for (nm in names(values)) {
    out <- gsub(
      paste0("{", nm, "}"),
      as.character(values[[nm]]),
      out,
      fixed = TRUE
    )
  }
  out
}

canonical_gene_id <- function(x) {
  x <- trimws(as.character(x))
  sub("\\.[0-9]+$", "", x)
}

normalize_chr <- function(x) {
  sub("^chr", "", as.character(x), ignore.case = TRUE)
}

safe_numeric <- function(x) suppressWarnings(as.numeric(x))

ensure_dir <- function(x) {
  dir.create(x, recursive = TRUE, showWarnings = FALSE)
  x
}

sanitize_component <- function(x) {
  x <- as.character(x)
  x[is.na(x) | !nzchar(x)] <- "NA"
  gsub("[^A-Za-z0-9._-]+", "_", x)
}

abs_path <- function(x) {
  d <- normalizePath(dirname(x), mustWork = TRUE)
  file.path(d, basename(x))
}

first_nonempty <- function(x) {
  x <- as.character(x)
  x <- x[!is.na(x) & nzchar(trimws(x))]
  if (length(x)) x[1] else NA_character_
}

# Keep manifest rows valid TSV records even when error messages contain
# embedded newlines or tabs.
clean_message <- function(x) {
  if (length(x) == 0L || is.na(x[1])) return(NA_character_)
  x <- as.character(x[1])
  x <- gsub("[\\r\\n\\t]+", " ", x)
  trimws(x)
}

# ------------------------------------------------------------------------------
# GENCODE gene map with strand-aware TSS
# ------------------------------------------------------------------------------
read_gene_map <- function(gtf_file) {
  if (!file.exists(gtf_file)) stop("Gene annotation not found: ", gtf_file)

  cmd <- paste0(
    "awk 'BEGIN{FS=\"\\t\"; OFS=\"\\t\"} ",
    "$3==\"gene\" {print $1,$4,$5,$7,$9}' ",
    shQuote(gtf_file)
  )

  g <- fread(
    cmd = cmd,
    header = FALSE,
    sep = "\t",
    data.table = FALSE,
    quote = ""
  )

  if (ncol(g) < 5L) stop("Could not parse gene records from GTF: ", gtf_file)
  names(g)[1:5] <- c("chr", "start", "end", "strand", "attribute")

  g$gene <- sub('.*gene_id "([^"]+)".*', "\\1", g$attribute)
  g$gene_name <- ifelse(
    grepl('gene_name "', g$attribute, fixed = TRUE),
    sub('.*gene_name "([^"]+)".*', "\\1", g$attribute),
    NA_character_
  )

  g$gene <- canonical_gene_id(g$gene)
  g$chr <- normalize_chr(g$chr)
  chr_num <- suppressWarnings(as.integer(g$chr))

  keep <- (
    !is.na(chr_num) & chr_num >= 1L & chr_num <= 22L &
      grepl("^ENSG[0-9]+$", g$gene) &
      g$strand %in% c("+", "-")
  )

  g <- g[keep, c("gene", "gene_name", "chr", "start", "end", "strand")]
  g <- g[!duplicated(g$gene), , drop = FALSE]
  g$tss <- ifelse(g$strand == "+", g$start, g$end)

  message("GENCODE gene map: ", nrow(g), " chr1-22 genes")
  g
}

GENE_MAP <- read_gene_map(CFG$gene_annotation)

# ------------------------------------------------------------------------------
# Build SMR tasks DIRECTLY from the TWAS-significant input tables.
#
# This intentionally does NOT depend on coloc-generated task files.
#
# NL:
#   one combined 7-cell-type Cauchy-meta significant table:
#   .../Cauchy_Meta_NL_WellcomeSanger_NL_EUR_Multi/
#      HTAN_NL_WellcomeSanger_NL_7celltypes_peakSNP.EUR_Multi.CauchyMeta.sig.csv
#
# ADSER:
#   the user-specified 9-cell summary is a COUNT/QC table and does not contain
#   gene-level rows. Therefore gene-level SMR tasks are read from the same
#   per-cell-type EUR/Multi significant-gene files used by coloc.
# ------------------------------------------------------------------------------

NL_SIG_FILE_DEFAULT <- paste0(
  "/data/l2_bioinfo1/liq17/scTWAS/",
  "TWAS_HTAN_GSE266546_WellcomeSanger_res/",
  "SPrediXcan_peakSNP_association/",
  "Cauchy_Meta_NL_WellcomeSanger_NL_EUR_Multi/",
  "HTAN_NL_WellcomeSanger_NL_7celltypes_peakSNP.",
  "EUR_Multi.CauchyMeta.sig.csv"
)

ADSER_SIG_DIR_DEFAULT <- paste0(
  "/data/l2_bioinfo1/liq17/scTWAS/",
  "TWAS_HTAN_GSE266546_WellcomeSanger_res/",
  "SPrediXcan_peakSNP_association/",
  "VUMC_HTAN_ADSER_EUR_Multi_sig"
)

ADSER_SIG_SUMMARY_QC <- file.path(
  ADSER_SIG_DIR_DEFAULT,
  "VUMC_HTAN_ADSER_EUR_Multi_sig_summary_9celltypes_R2global_gt0.01.csv"
)

is_yes <- function(x) {
  tolower(trimws(as.character(x))) %in%
    c("yes", "y", "true", "t", "1", "significant", "sig")
}

first_existing <- function(paths) {
  paths <- unique(as.character(paths))
  paths <- paths[!is.na(paths) & nzchar(paths)]
  hit <- paths[file.exists(paths)]
  if (length(hit)) hit[1] else NA_character_
}

find_column <- function(x, candidates, required = TRUE, label = "column") {
  hit <- candidates[candidates %in% names(x)]
  if (length(hit)) return(hit[1])
  if (required) {
    stop(
      "Could not identify ", label, ". Tried: ",
      paste(candidates, collapse = ", "),
      "\nAvailable columns: ", paste(names(x), collapse = ", ")
    )
  }
  NA_character_
}

get_gene_name_from_sig_row <- function(row_df) {
  candidates <- c(
    grep("^gene_name", names(row_df), value = TRUE),
    "Gene", "gene_symbol", "symbol", "SYMBOL"
  )
  candidates <- unique(candidates[candidates %in% names(row_df)])
  if (!length(candidates)) return(NA_character_)

  vals <- unlist(
    row_df[1, candidates, drop = FALSE],
    use.names = FALSE
  )
  first_nonempty(vals)
}

resolve_configured_source <- function(kind) {
  available <- names(CFG$eqtl_sources)
  if (!length(available)) {
    stop("CFG$eqtl_sources has no named eQTL sources.")
  }

  if (identical(kind, "wellcome")) {
    hit <- available[grepl("wellcome|sanger", available, ignore.case = TRUE)]
  } else if (identical(kind, "htan")) {
    hit <- available[grepl("htan|vumc", available, ignore.case = TRUE)]
  } else {
    hit <- character()
  }

  if (length(hit) == 1L) return(hit)

  stop(
    "Could not uniquely map eQTL source type '", kind,
    "' to CFG$eqtl_sources. Available sources: ",
    paste(available, collapse = ", ")
  )
}

normalize_source_to_config <- function(x) {
  x <- trimws(as.character(x)[1])
  if (is.na(x) || !nzchar(x)) return(NA_character_)

  available <- names(CFG$eqtl_sources)

  # Exact config name first.
  exact <- available[tolower(available) == tolower(x)]
  if (length(exact) == 1L) return(exact)

  if (grepl("wellcome|sanger", x, ignore.case = TRUE)) {
    return(resolve_configured_source("wellcome"))
  }

  if (grepl("htan|vumc", x, ignore.case = TRUE)) {
    return(resolve_configured_source("htan"))
  }

  NA_character_
}

resolve_nl_source <- function(row_df, cohort) {

  cfg_col <- NULL
  if (!is.null(CFG$source_columns) &&
      !is.null(CFG$source_columns[[cohort]])) {
    cfg_col <- CFG$source_columns[[cohort]]
  }

  source_candidates <- unique(c(
    cfg_col,
    paste0("P driven_", cohort),
    paste0("P_driven_", cohort),
    paste0("P.driven_", cohort)
  ))

  for (nm in source_candidates) {
    if (!is.null(nm) && !is.na(nm) && nm %in% names(row_df)) {
      src <- normalize_source_to_config(row_df[[nm]][1])
      if (!is.na(src) && nzchar(src)) return(src)
    }
  }

  # Same fallback logic as coloc: source text.
  meta_col <- paste0("META_SOURCE_", cohort)
  if (meta_col %in% names(row_df)) {
    src <- normalize_source_to_config(row_df[[meta_col]][1])
    if (!is.na(src) && nzchar(src)) return(src)
  }

  # Final fallback: choose smaller source-specific TWAS P.
  p_htan <- paste0("pvalue_HTAN_NL_", cohort)
  p_ws   <- paste0("pvalue_WellcomeSanger_NL_", cohort)

  if (all(c(p_htan, p_ws) %in% names(row_df))) {
    ph <- safe_numeric(row_df[[p_htan]][1])
    pw <- safe_numeric(row_df[[p_ws]][1])

    if (is.finite(ph) && is.finite(pw)) {
      if (ph <= pw) return(resolve_configured_source("htan"))
      return(resolve_configured_source("wellcome"))
    }
    if (is.finite(ph)) return(resolve_configured_source("htan"))
    if (is.finite(pw)) return(resolve_configured_source("wellcome"))
  }

  NA_character_
}

resolve_adser_sig_file <- function(CT) {

  # Optional config override.
  cfg_candidates <- character()

  if (!is.null(CFG$twas_sig_file_template)) {
    cfg_candidates <- c(
      cfg_candidates,
      fill_template(CFG$twas_sig_file_template, list(CT = CT))
    )
  }

  # Support both basename conventions seen in this project.
  candidates <- c(
    cfg_candidates,
    file.path(
      ADSER_SIG_DIR_DEFAULT,
      paste0(
        "VUMC_HTAN_ADSER_", CT,
        "_peakSNP.EUR_Multi.sig.csv"
      )
    ),
    file.path(
      ADSER_SIG_DIR_DEFAULT,
      paste0(
        "VUMC_HTAN_GSE266546_ADSER_", CT,
        "_peakSNP.EUR_Multi.sig.csv"
      )
    ),
    file.path(
      "/v5000/data/l2_bioinfo1/liq17/scTWAS/",
      "TWAS_HTAN_GSE266546_WellcomeSanger_res/",
      "SPrediXcan_peakSNP_association/",
      "VUMC_HTAN_ADSER_EUR_Multi_sig/",
      paste0(
        "VUMC_HTAN_ADSER_", CT,
        "_peakSNP.EUR_Multi.sig.csv"
      )
    ),
    file.path(
      "/v5000/data/l2_bioinfo1/liq17/scTWAS/",
      "TWAS_HTAN_GSE266546_WellcomeSanger_res/",
      "SPrediXcan_peakSNP_association/",
      "VUMC_HTAN_ADSER_EUR_Multi_sig/",
      paste0(
        "VUMC_HTAN_GSE266546_ADSER_", CT,
        "_peakSNP.EUR_Multi.sig.csv"
      )
    )
  )

  path <- first_existing(candidates)

  if (is.na(path)) {
    stop(
      "ADSER gene-level significant-gene file not found for ", CT, ".\n",
      "The 9-cell-type *_sig_summary_*.csv is count-only and cannot be used ",
      "to generate gene-level SMR tasks.\nTried:\n",
      paste("  -", unique(candidates), collapse = "\n")
    )
  }

  path
}

read_sig_for_celltype <- function(CT) {

  GROUP <- toupper(trimws(as.character(CFG$analysis_group)))

  if (GROUP == "NL") {

    path <- if (!is.null(CFG$twas_sig_file) &&
                length(CFG$twas_sig_file) == 1L &&
                !is.na(CFG$twas_sig_file) &&
                nzchar(CFG$twas_sig_file)) {
      CFG$twas_sig_file
    } else {
      NL_SIG_FILE_DEFAULT
    }

    if (!file.exists(path)) {
      stop("NL combined TWAS significant-gene file not found:\n", path)
    }

    message("TWAS task source: ", path)

    sig <- fread(
      path,
      data.table = FALSE,
      check.names = FALSE
    )

    ct_col <- find_column(
      sig,
      c(
        "CellType", "cell_type", "Cell_Type",
        "Cell Type", "risk_celltype", "CT"
      ),
      TRUE,
      "NL cell-type column"
    )

    sig <- sig[
      toupper(trimws(as.character(sig[[ct_col]]))) == CT,
      ,
      drop = FALSE
    ]

  } else if (GROUP == "ADSER") {

    path <- resolve_adser_sig_file(CT)

    message("TWAS task source: ", path)

    if (file.exists(ADSER_SIG_SUMMARY_QC)) {
      message(
        "ADSER 9-cell summary retained as QC only (not task input): ",
        ADSER_SIG_SUMMARY_QC
      )
    }

    sig <- fread(
      path,
      data.table = FALSE,
      check.names = FALSE
    )

  } else {
    stop(
      "Unsupported analysis_group for direct TWAS task generation: ",
      CFG$analysis_group,
      ". Expected NL or ADSER."
    )
  }

  attr(sig, "sig_file") <- path
  sig
}

build_tasks_for_celltype <- function(CT) {

  sig <- read_sig_for_celltype(CT)
  sig_file <- attr(sig, "sig_file")

  gene_col <- if (!is.null(CFG$gene_column) &&
                  CFG$gene_column %in% names(sig)) {
    CFG$gene_column
  } else {
    find_column(
      sig,
      c("gene", "ENSG", "ENSG_no_version", "phenotype_id"),
      TRUE,
      "gene column"
    )
  }

  eur_flag_col <- if (!is.null(CFG$sig_flag_columns) &&
                      !is.null(CFG$sig_flag_columns[["EUR"]]) &&
                      CFG$sig_flag_columns[["EUR"]] %in% names(sig)) {
    CFG$sig_flag_columns[["EUR"]]
  } else {
    find_column(sig, c("EUR_sig", "EUR.significant"), TRUE, "EUR significance flag")
  }

  multi_flag_col <- if (!is.null(CFG$sig_flag_columns) &&
                        !is.null(CFG$sig_flag_columns[["Multi"]]) &&
                        CFG$sig_flag_columns[["Multi"]] %in% names(sig)) {
    CFG$sig_flag_columns[["Multi"]]
  } else {
    find_column(sig, c("Multi_sig", "Multi.significant"), TRUE, "Multi significance flag")
  }

  sig$gene_key <- canonical_gene_id(sig[[gene_col]])
  sig <- sig[
    !is.na(sig$gene_key) &
      nzchar(sig$gene_key),
    ,
    drop = FALSE
  ]

  # One biological row per gene for this cell type.
  sig <- sig[!duplicated(sig$gene_key), , drop = FALSE]

  map_idx <- match(sig$gene_key, GENE_MAP$gene)
  sig$chromosome_gencode <- GENE_MAP$chr[map_idx]
  sig$gene_name_gencode <- GENE_MAP$gene_name[map_idx]

  tasks <- list()
  k <- 0L

  for (i in seq_len(nrow(sig))) {

    row_df <- sig[i, , drop = FALSE]

    gene <- row_df$gene_key[1]
    chrom <- normalize_chr(row_df$chromosome_gencode[1])

    if (is.na(chrom) || !chrom %in% as.character(1:22)) {
      message("Skipping gene without chr1-22 GENCODE mapping: ", gene)
      next
    }

    gene_name <- get_gene_name_from_sig_row(row_df)
    if (is.na(gene_name) || !nzchar(gene_name)) {
      gene_name <- row_df$gene_name_gencode[1]
    }

    for (cohort in c("EUR", "Multi")) {

      flag_col <- if (cohort == "EUR") eur_flag_col else multi_flag_col

      if (!is_yes(row_df[[flag_col]][1])) next

      source <- if (toupper(CFG$analysis_group) == "NL") {

        resolve_nl_source(row_df, cohort)

      } else {

        # ADSER is HTAN/VUMC eQTL only.
        if (!is.null(CFG$fixed_eqtl_source) &&
            CFG$fixed_eqtl_source %in% names(CFG$eqtl_sources)) {
          CFG$fixed_eqtl_source
        } else {
          resolve_configured_source("htan")
        }
      }

      if (is.na(source) || !nzchar(source)) {
        message(
          "Skipping unresolved eQTL source: ",
          CT, " ", gene, " ", cohort
        )
        next
      }

      if (!source %in% names(CFG$eqtl_sources)) {
        stop(
          "Resolved eQTL source '", source,
          "' is not configured. Available: ",
          paste(names(CFG$eqtl_sources), collapse = ", ")
        )
      }

      p_candidates <- if (toupper(CFG$analysis_group) == "NL") {
        c(
          paste0("META_P_", cohort),
          paste0("pvalue_", cohort)
        )
      } else {
        c(
          paste0("pvalue_", cohort),
          paste0("P_", cohort),
          paste0("META_P_", cohort)
        )
      }

      p_col <- p_candidates[p_candidates %in% names(row_df)]
      TWAS_p <- if (length(p_col)) {
        safe_numeric(row_df[[p_col[1]]][1])
      } else {
        NA_real_
      }

      k <- k + 1L

      tasks[[k]] <- data.frame(
        cell_type = CT,
        gene = gene,
        gene_name = gene_name,
        chromosome = chrom,
        cohort = cohort,
        eqtl_source = source,
        TWAS_p = TWAS_p,
        sig_file = sig_file,
        stringsAsFactors = FALSE,
        check.names = FALSE
      )
    }
  }

  if (!length(tasks)) {
    return(data.frame())
  }

  x <- bind_rows(tasks)

  x <- x[
    !duplicated(
      paste(
        x$cell_type,
        x$gene,
        x$cohort,
        x$eqtl_source,
        sep = "|"
      )
    ),
    ,
    drop = FALSE
  ]

  x
}

# ------------------------------------------------------------------------------
# GWAS reader
# ------------------------------------------------------------------------------
GWAS_CACHE <- new.env(parent = emptyenv())

read_gwas_chr <- function(cohort, chrom) {
  key <- paste(cohort, chrom, sep = "|")
  if (exists(key, envir = GWAS_CACHE, inherits = FALSE)) {
    return(get(key, envir = GWAS_CACHE, inherits = FALSE))
  }

  gcfg <- CFG$gwas[[cohort]]
  if (is.null(gcfg)) stop("GWAS cohort not configured: ", cohort)

  path <- fill_template(gcfg$file_template, list(chrom = chrom))
  if (!file.exists(path)) stop("GWAS file not found: ", path)

  g <- fread(path, data.table = FALSE, check.names = FALSE)
  req <- c(
    gcfg$rsid_col, gcfg$a1_col, gcfg$a2_col, gcfg$freq_col,
    gcfg$beta_col, gcfg$se_col, gcfg$p_col
  )
  miss <- setdiff(req, names(g))
  if (length(miss)) {
    stop(
      "GWAS file missing required column(s): ", paste(miss, collapse = ", "),
      "\nFile: ", path
    )
  }

  out <- data.frame(
    SNP = tolower(trimws(as.character(g[[gcfg$rsid_col]]))),
    A1 = toupper(trimws(as.character(g[[gcfg$a1_col]]))),
    A2 = toupper(trimws(as.character(g[[gcfg$a2_col]]))),
    FRQ = safe_numeric(g[[gcfg$freq_col]]),
    BETA_g = safe_numeric(g[[gcfg$beta_col]]),
    SE_g = safe_numeric(g[[gcfg$se_col]]),
    P_g = safe_numeric(g[[gcfg$p_col]]),
    stringsAsFactors = FALSE
  )

  keep <- (
    grepl("^rs[0-9]+$", out$SNP, ignore.case = TRUE) &
      out$A1 %in% c("A", "C", "G", "T") &
      out$A2 %in% c("A", "C", "G", "T") &
      out$A1 != out$A2 &
      is.finite(out$BETA_g) &
      is.finite(out$SE_g) & out$SE_g > 0 &
      is.finite(out$P_g) & out$P_g > 0 & out$P_g <= 1
  )

  out <- out[keep, , drop = FALSE]
  out <- out[order(out$P_g), , drop = FALSE]
  out <- out[!duplicated(out$SNP), , drop = FALSE]

  assign(key, out, envir = GWAS_CACHE)
  message(cohort, " chr", chrom, " GWAS variants loaded: ", nrow(out))
  out
}

# ------------------------------------------------------------------------------
# TensorQTL eQTL reader
# ------------------------------------------------------------------------------
EQTL_CACHE <- new.env(parent = emptyenv())

parse_tensorqtl_variant <- function(variant_id) {
  x <- as.character(variant_id)
  rsid <- ifelse(
    grepl(";", x, fixed = TRUE),
    sub(";.*$", "", x),
    NA_character_
  )
  key <- ifelse(
    grepl(";", x, fixed = TRUE),
    sub("^[^;]+;", "", x),
    x
  )
  parts <- strsplit(key, "_", fixed = TRUE)

  data.frame(
    SNP = tolower(trimws(rsid)),
    variant_chr = normalize_chr(vapply(
      parts, function(z) if (length(z) >= 4L) z[1] else NA_character_, character(1)
    )),
    variant_pos = safe_numeric(vapply(
      parts, function(z) if (length(z) >= 4L) z[2] else NA_character_, character(1)
    )),
    REF = toupper(trimws(vapply(
      parts, function(z) if (length(z) >= 4L) z[3] else NA_character_, character(1)
    ))),
    ALT = toupper(trimws(vapply(
      parts, function(z) if (length(z) >= 4L) z[4] else NA_character_, character(1)
    ))),
    stringsAsFactors = FALSE
  )
}

read_eqtl_chr <- function(source, CT, chrom) {
  key <- paste(source, CT, chrom, sep = "|")
  if (exists(key, envir = EQTL_CACHE, inherits = FALSE)) {
    return(get(key, envir = EQTL_CACHE, inherits = FALSE))
  }

  scfg <- CFG$eqtl_sources[[source]]
  if (is.null(scfg)) stop("eQTL source not configured: ", source)

  path <- fill_template(
    scfg$file_template,
    list(CT = CT, chrom = chrom)
  )
  if (!file.exists(path)) {
    stop("eQTL parquet not found for source=", source, ", CT=", CT, ", chr=", chrom, ":\n", path)
  }

  ds <- arrow::open_dataset(path, format = "parquet")
  schema_names <- as.character(ds$schema$names)
  rm(ds)

  required <- c("phenotype_id", "variant_id", "slope", "slope_se", "pval_nominal")
  miss <- setdiff(required, schema_names)
  if (length(miss)) {
    stop(
      "eQTL parquet missing required column(s): ", paste(miss, collapse = ", "),
      "\nFile: ", path
    )
  }

  optional <- intersect(c("tss_distance", "af", "ma_samples", "ma_count"), schema_names)
  eq <- as.data.frame(
    arrow::read_parquet(path, col_select = c(required, optional))
  )

  eq$phenotype_key <- canonical_gene_id(eq$phenotype_id)
  parsed <- parse_tensorqtl_variant(eq$variant_id)
  eq <- cbind(eq, parsed)

  eq$BETA_x <- safe_numeric(eq$slope)
  eq$SE_x <- safe_numeric(eq$slope_se)
  eq$P_x <- safe_numeric(eq$pval_nominal)
  eq$AF <- if ("af" %in% names(eq)) safe_numeric(eq$af) else NA_real_

  keep <- (
    grepl("^rs[0-9]+$", eq$SNP, ignore.case = TRUE) &
      eq$variant_chr == as.character(chrom) &
      is.finite(eq$variant_pos) &
      eq$REF %in% c("A", "C", "G", "T") &
      eq$ALT %in% c("A", "C", "G", "T") &
      eq$REF != eq$ALT &
      is.finite(eq$BETA_x) &
      is.finite(eq$SE_x) & eq$SE_x > 0 &
      is.finite(eq$P_x) & eq$P_x > 0 & eq$P_x <= 1
  )

  eq <- eq[keep, , drop = FALSE]
  eq <- eq[order(eq$P_x), , drop = FALSE]
  eq <- eq[!duplicated(paste(eq$phenotype_key, eq$SNP, sep = "|")), , drop = FALSE]

  assign(key, eq, envir = EQTL_CACHE)
  message(source, " ", CT, " chr", chrom, " eQTL rows loaded: ", nrow(eq))
  eq
}

# ------------------------------------------------------------------------------
# Allele harmonization
# ------------------------------------------------------------------------------
complement_allele <- function(x) chartr("ACGT", "TGCA", toupper(x))

is_ambiguous_pair <- function(a1, a2) {
  (a1 == "A" & a2 == "T") |
    (a1 == "T" & a2 == "A") |
    (a1 == "C" & a2 == "G") |
    (a1 == "G" & a2 == "C")
}

harmonize_gwas_eqtl <- function(gwas, eqtl) {
  x <- merge(gwas, eqtl, by = "SNP")
  if (!nrow(x)) return(x)

  a1g <- toupper(x$A1)
  a2g <- toupper(x$A2)
  a1x <- toupper(x$ALT)  # TensorQTL effect allele
  a2x <- toupper(x$REF)

  ambiguous <- is_ambiguous_pair(a1g, a2g)
  same <- a1x == a1g & a2x == a2g
  flip <- a1x == a2g & a2x == a1g
  comp_same <- a1x == complement_allele(a1g) & a2x == complement_allele(a2g)
  comp_flip <- a1x == complement_allele(a2g) & a2x == complement_allele(a1g)

  keep <- same | flip | comp_same | comp_flip
  if (isTRUE(CFG$drop_strand_ambiguous)) keep <- keep & !ambiguous

  x <- x[keep, , drop = FALSE]
  if (!nrow(x)) return(x)

  a1g <- toupper(x$A1)
  a2g <- toupper(x$A2)
  a1x <- toupper(x$ALT)
  a2x <- toupper(x$REF)

  same <- a1x == a1g & a2x == a2g
  flip <- a1x == a2g & a2x == a1g
  comp_same <- a1x == complement_allele(a1g) & a2x == complement_allele(a2g)
  comp_flip <- a1x == complement_allele(a2g) & a2x == complement_allele(a1g)
  beta_flip <- flip | comp_flip

  x$alignment <- ifelse(
    same, "same",
    ifelse(flip, "flip", ifelse(comp_same, "comp_same", "comp_flip"))
  )

  # Align eQTL effect to GWAS A1.
  x$BETA_x_aligned <- x$BETA_x
  x$BETA_x_aligned[beta_flip] <- -x$BETA_x_aligned[beta_flip]

  # TensorQTL AF is ALT frequency. After beta flipping, the aligned eQTL
  # effect allele is REF (= GWAS A1), hence frequency must become 1-AF.
  x$AF_x_aligned <- x$AF
  idx <- beta_flip & is.finite(x$AF_x_aligned)
  x$AF_x_aligned[idx] <- 1 - x$AF_x_aligned[idx]

  # Prefer each study's own frequency; use the other study as fallback only if
  # the required frequency is absent/invalid.
  x$GWAS_A1_freq <- x$FRQ
  bad_g <- !is.finite(x$GWAS_A1_freq) | x$GWAS_A1_freq <= 0 | x$GWAS_A1_freq >= 1
  x$GWAS_A1_freq[bad_g] <- x$AF_x_aligned[bad_g]

  x$EQTL_A1_freq <- x$AF_x_aligned
  bad_x <- !is.finite(x$EQTL_A1_freq) | x$EQTL_A1_freq <= 0 | x$EQTL_A1_freq >= 1
  x$EQTL_A1_freq[bad_x] <- x$GWAS_A1_freq[bad_x]

  x <- x[
    is.finite(x$GWAS_A1_freq) & x$GWAS_A1_freq > 0 & x$GWAS_A1_freq < 1 &
      is.finite(x$EQTL_A1_freq) & x$EQTL_A1_freq > 0 & x$EQTL_A1_freq < 1 &
      is.finite(x$BETA_x_aligned) &
      is.finite(x$SE_x) & x$SE_x > 0 &
      is.finite(x$P_x) & x$P_x > 0 & x$P_x <= 1,
    ,
    drop = FALSE
  ]

  x <- x[order(x$P_x, x$P_g), , drop = FALSE]
  x <- x[!duplicated(x$SNP), , drop = FALSE]
  x
}

# ------------------------------------------------------------------------------
# Resolve chromosome PGEN prefix used as HEIDI LD reference.
# ------------------------------------------------------------------------------
resolve_pfile_prefix <- function(source, chrom) {
  scfg <- CFG$eqtl_sources[[source]]

  if (!is.null(scfg$pfile_template) &&
      length(scfg$pfile_template) == 1L &&
      !is.na(scfg$pfile_template) &&
      nzchar(scfg$pfile_template)) {
    p <- fill_template(scfg$pfile_template, list(chrom = chrom))
    if (file.exists(paste0(p, ".pgen")) &&
        file.exists(paste0(p, ".pvar")) &&
        file.exists(paste0(p, ".psam"))) {
      return(p)
    }
    warning("Configured pfile_template does not resolve to .pgen/.pvar/.psam: ", p)
  }

  dirs <- scfg$pgen_dirs
  if (is.null(dirs)) return(NA_character_)
  dirs <- as.character(dirs)
  dirs <- dirs[dir.exists(dirs)]
  if (!length(dirs)) return(NA_character_)

  candidates <- unlist(lapply(
    dirs,
    function(d) list.files(
      d,
      pattern = paste0("^chr", chrom, "(\\.|_|$).*\\.pgen$"),
      full.names = TRUE,
      ignore.case = TRUE
    )
  ), use.names = FALSE)

  if (!length(candidates)) return(NA_character_)

  # Prefer rsID-ready files when more than one chromosome PGEN is present.
  preferred <- candidates[grepl("rsid", basename(candidates), ignore.case = TRUE)]
  if (length(preferred) == 1L) candidates <- preferred
  else if (length(preferred) > 1L) candidates <- preferred

  # Prefer names containing qc.rsid if available.
  preferred2 <- candidates[grepl("qc.*rsid|rsid.*qc", basename(candidates), ignore.case = TRUE)]
  if (length(preferred2) == 1L) candidates <- preferred2

  # Only accept a unique prefix with companion .pvar/.psam.
  prefixes <- sub("\\.pgen$", "", candidates, ignore.case = TRUE)
  ok <- file.exists(paste0(prefixes, ".pvar")) & file.exists(paste0(prefixes, ".psam"))
  prefixes <- unique(prefixes[ok])

  if (length(prefixes) == 1L) return(prefixes)

  if (length(prefixes) > 1L) {
    warning(
      "Multiple PGEN prefixes found for source=", source, " chr", chrom,
      ". Set eqtl_sources[[source]]$pfile_template explicitly. Candidates: ",
      paste(prefixes, collapse = "; ")
    )
  }
  NA_character_
}

# ------------------------------------------------------------------------------
# Prepare one task
# ------------------------------------------------------------------------------
empty_manifest_row <- function(task, status, message_text = NA_character_) {
  data.frame(
    analysis_group = CFG$analysis_group,
    cell_type = task$cell_type,
    gene = task$gene,
    gene_name = task$gene_name,
    chromosome = task$chromosome,
    cohort = task$cohort,
    eqtl_source = task$eqtl_source,
    TWAS_p = task$TWAS_p,
    tss = NA_real_,
    region_start = NA_real_,
    region_end = NA_real_,
    n_eqtl_region = NA_integer_,
    n_overlap_harmonized = NA_integer_,
    target_SNP = NA_character_,
    target_pos = NA_real_,
    target_eQTL_p = NA_real_,
    target_eQTL_beta = NA_real_,
    target_eQTL_se = NA_real_,
    target_eQTL_freq = NA_real_,
    target_GWAS_p = NA_real_,
    target_GWAS_beta = NA_real_,
    target_GWAS_se = NA_real_,
    target_GWAS_freq = NA_real_,
    gwas_ma = NA_character_,
    eqtl_esd = NA_character_,
    eqtl_flist = NA_character_,
    beqtl_prefix = NA_character_,
    snplist = NA_character_,
    target_qc = NA_character_,
    pfile_prefix = NA_character_,
    task_out_prefix = NA_character_,
    peqtl_smr = CFG$smr$peqtl_smr,
    peqtl_heidi = CFG$smr$peqtl_heidi,
    heidi_mtd = CFG$smr$heidi_mtd,
    heidi_min_m = CFG$smr$heidi_min_m,
    heidi_max_m = CFG$smr$heidi_max_m,
    ld_lower_limit = CFG$smr$ld_lower_limit,
    ld_upper_limit = CFG$smr$ld_upper_limit,
    maf = CFG$smr$maf,
    diff_freq = CFG$smr$diff_freq,
    diff_freq_prop = CFG$smr$diff_freq_prop,
    thread_num = CFG$smr$thread_num,
    status = status,
    message = clean_message(message_text),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}

prepare_one_task <- function(task) {
  gene <- canonical_gene_id(task$gene)
  CT <- toupper(task$cell_type)
  cohort <- as.character(task$cohort)
  source <- as.character(task$eqtl_source)
  chrom <- normalize_chr(task$chromosome)

  message(
    "\n------------------------------------------------------------\n",
    "CT=", CT, " gene=", gene, " cohort=", cohort,
    " eQTL=", source, " chr=", chrom
  )

  gm <- GENE_MAP[GENE_MAP$gene == gene, , drop = FALSE]
  if (!nrow(gm)) return(empty_manifest_row(task, "gene_not_in_gencode"))
  if (nrow(gm) > 1L) gm <- gm[1, , drop = FALSE]

  tss <- as.numeric(gm$tss[1])
  region_start <- max(1, tss - as.numeric(CFG$flank_length))
  region_end <- tss + as.numeric(CFG$flank_length)

  gwas <- read_gwas_chr(cohort, chrom)
  eqchr <- read_eqtl_chr(source, CT, chrom)

  eqgene <- eqchr[eqchr$phenotype_key == gene, , drop = FALSE]
  if (!nrow(eqgene)) {
    out <- empty_manifest_row(task, "gene_not_found_in_eqtl")
    out$tss <- tss; out$region_start <- region_start; out$region_end <- region_end
    out$n_eqtl_region <- 0L
    return(out)
  }

  # Use both explicit genomic position and tensorQTL tss_distance when present.
  keep_region <- (
    eqgene$variant_pos >= region_start & eqgene$variant_pos <= region_end
  )
  if ("tss_distance" %in% names(eqgene)) {
    td <- safe_numeric(eqgene$tss_distance)
    keep_region <- keep_region & is.finite(td) & abs(td) <= CFG$flank_length
  }
  eqgene <- eqgene[keep_region, , drop = FALSE]

  if (!nrow(eqgene)) {
    out <- empty_manifest_row(task, "no_eqtl_in_500kb_region")
    out$tss <- tss; out$region_start <- region_start; out$region_end <- region_end
    out$n_eqtl_region <- 0L
    return(out)
  }

  matched <- harmonize_gwas_eqtl(gwas, eqgene)
  n_overlap <- nrow(matched)

  if (n_overlap < CFG$min_overlap_snps) {
    out <- empty_manifest_row(
      task,
      "too_few_harmonized_overlap_snps",
      paste0(n_overlap, " < ", CFG$min_overlap_snps)
    )
    out$tss <- tss; out$region_start <- region_start; out$region_end <- region_end
    out$n_eqtl_region <- nrow(eqgene)
    out$n_overlap_harmonized <- n_overlap
    return(out)
  }

  # The user's requested target definition: smallest eQTL p-value among the
  # harmonized GWAS x eQTL overlap within gene TSS +/-500 kb.
  matched <- matched[order(matched$P_x, matched$P_g, matched$SNP), , drop = FALSE]
  target <- matched[1, , drop = FALSE]

  task_dir <- ensure_dir(file.path(
    CFG$output_dir,
    CT,
    cohort,
    sanitize_component(source),
    sanitize_component(gene)
  ))

  base <- paste(
    CFG$analysis_group, CT, cohort, sanitize_component(source), gene,
    sep = "."
  )

  gwas_ma <- file.path(task_dir, paste0(base, ".gwas_500kb.ma"))
  eqtl_esd <- file.path(task_dir, paste0(base, ".eqtl_500kb.esd"))
  eqtl_flist <- file.path(task_dir, paste0(base, ".eqtl_500kb.flist"))
  beqtl_prefix <- file.path(task_dir, paste0(base, ".eqtl_500kb"))
  snplist <- file.path(task_dir, paste0(base, ".harmonized_snps.txt"))
  target_qc <- file.path(task_dir, paste0(base, ".target_snp_qc.tsv"))
  task_out_prefix <- file.path(task_dir, paste0(base, ".SMR_HEIDI"))

  # GWAS .ma format required by SMR/GCTA-COJO:
  # SNP A1 A2 freq b se p n
  gwas_out <- data.frame(
    SNP = matched$SNP,
    A1 = matched$A1,
    A2 = matched$A2,
    freq = matched$GWAS_A1_freq,
    b = matched$BETA_g,
    se = matched$SE_g,
    p = matched$P_g,
    n = CFG$gwas[[cohort]]$N,
    stringsAsFactors = FALSE
  )
  fwrite(gwas_out, gwas_ma, sep = "\t", quote = FALSE, na = "NA")

  # ESD format required by SMR. We intentionally align eQTL effect to the same
  # A1/A2 used in GWAS, so Beta/Freq are also A1-oriented.
  esd_out <- data.frame(
    Chr = as.integer(chrom),
    SNP = matched$SNP,
    Bp = as.integer(round(matched$variant_pos)),
    A1 = matched$A1,
    A2 = matched$A2,
    Freq = matched$EQTL_A1_freq,
    Beta = matched$BETA_x_aligned,
    se = matched$SE_x,
    p = matched$P_x,
    stringsAsFactors = FALSE
  )
  esd_out <- esd_out[order(esd_out$Bp, esd_out$SNP), , drop = FALSE]
  fwrite(esd_out, eqtl_esd, sep = "\t", quote = FALSE, na = "NA")

  orientation <- as.character(gm$strand[1])
  gene_name <- first_nonempty(c(task$gene_name, gm$gene_name[1], gene))
  flist_out <- data.frame(
    Chr = as.integer(chrom),
    ProbeID = gene,
    GeneticDistance = 0,
    ProbeBp = as.integer(round(tss)),
    Gene = gene_name,
    Orientation = orientation,
    PathOfEsd = abs_path(eqtl_esd),
    stringsAsFactors = FALSE
  )
  fwrite(flist_out, eqtl_flist, sep = "\t", quote = FALSE, na = "NA")

  fwrite(data.frame(SNP = matched$SNP), snplist, sep = "\t", col.names = FALSE)

  target_out <- data.frame(
    analysis_group = CFG$analysis_group,
    cell_type = CT,
    gene = gene,
    gene_name = gene_name,
    chromosome = chrom,
    cohort = cohort,
    eqtl_source = source,
    tss = tss,
    region_start = region_start,
    region_end = region_end,
    target_SNP = target$SNP,
    target_pos = target$variant_pos,
    GWAS_A1 = target$A1,
    GWAS_A2 = target$A2,
    original_eQTL_REF = target$REF,
    original_eQTL_ALT = target$ALT,
    allele_alignment = target$alignment,
    target_eQTL_beta_aligned_A1 = target$BETA_x_aligned,
    target_eQTL_se = target$SE_x,
    target_eQTL_p = target$P_x,
    target_eQTL_A1_freq = target$EQTL_A1_freq,
    target_GWAS_beta = target$BETA_g,
    target_GWAS_se = target$SE_g,
    target_GWAS_p = target$P_g,
    target_GWAS_A1_freq = target$GWAS_A1_freq,
    n_eqtl_region = nrow(eqgene),
    n_overlap_harmonized = n_overlap,
    stringsAsFactors = FALSE
  )
  fwrite(target_out, target_qc, sep = "\t", quote = FALSE, na = "NA")

  pfile_prefix <- resolve_pfile_prefix(source, chrom)
  status <- if (is.na(pfile_prefix)) "prepared_missing_pfile_prefix" else "prepared"
  msg <- if (is.na(pfile_prefix)) {
    paste0(
      "Could not uniquely resolve chr", chrom,
      " PGEN LD reference. Set eqtl_sources[[", source,
      "]]$pfile_template in config."
    )
  } else NA_character_

  data.frame(
    analysis_group = CFG$analysis_group,
    cell_type = CT,
    gene = gene,
    gene_name = gene_name,
    chromosome = chrom,
    cohort = cohort,
    eqtl_source = source,
    TWAS_p = task$TWAS_p,
    tss = tss,
    region_start = region_start,
    region_end = region_end,
    n_eqtl_region = nrow(eqgene),
    n_overlap_harmonized = n_overlap,
    target_SNP = target$SNP,
    target_pos = target$variant_pos,
    target_eQTL_p = target$P_x,
    target_eQTL_beta = target$BETA_x_aligned,
    target_eQTL_se = target$SE_x,
    target_eQTL_freq = target$EQTL_A1_freq,
    target_GWAS_p = target$P_g,
    target_GWAS_beta = target$BETA_g,
    target_GWAS_se = target$SE_g,
    target_GWAS_freq = target$GWAS_A1_freq,
    gwas_ma = abs_path(gwas_ma),
    eqtl_esd = abs_path(eqtl_esd),
    eqtl_flist = abs_path(eqtl_flist),
    beqtl_prefix = abs_path(beqtl_prefix),
    snplist = abs_path(snplist),
    target_qc = abs_path(target_qc),
    pfile_prefix = ifelse(is.na(pfile_prefix), NA_character_, pfile_prefix),
    task_out_prefix = abs_path(task_out_prefix),
    peqtl_smr = CFG$smr$peqtl_smr,
    peqtl_heidi = CFG$smr$peqtl_heidi,
    heidi_mtd = CFG$smr$heidi_mtd,
    heidi_min_m = CFG$smr$heidi_min_m,
    heidi_max_m = CFG$smr$heidi_max_m,
    ld_lower_limit = CFG$smr$ld_lower_limit,
    ld_upper_limit = CFG$smr$ld_upper_limit,
    maf = CFG$smr$maf,
    diff_freq = CFG$smr$diff_freq,
    diff_freq_prop = CFG$smr$diff_freq_prop,
    thread_num = CFG$smr$thread_num,
    status = status,
    message = msg,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}

# ------------------------------------------------------------------------------
# Run requested cell types
# ------------------------------------------------------------------------------
allowed_ct <- toupper(CFG$cell_types)
if (CT_ARG == "ALL") {
  CTs <- allowed_ct
} else {
  if (!CT_ARG %in% allowed_ct) {
    stop(
      "Cell type ", CT_ARG, " is not configured. Allowed: ",
      paste(c("all", allowed_ct), collapse = ", ")
    )
  }
  CTs <- CT_ARG
}

ensure_dir(CFG$output_dir)
all_manifest <- list()
k <- 0L

for (CT in CTs) {
  message(
    "\n======================================================================\n",
    "PREPARE SMR+HEIDI ", CFG$analysis_group, " ", CT, "\n",
    "======================================================================"
  )

  tasks <- build_tasks_for_celltype(CT)
  if (!nrow(tasks)) {
    warning("No TWAS-significant SMR tasks generated for ", CT)
    next
  }

  task_file <- file.path(
    CFG$output_dir,
    paste0(
      "SMR_HEIDI.", CFG$analysis_group, ".", CT,
      ".significant_gene_tasks.from_TWAS.csv"
    )
  )
  fwrite(tasks, task_file)
  message("Wrote direct-TWAS task table: ", task_file)

  rows <- vector("list", nrow(tasks))
  for (i in seq_len(nrow(tasks))) {
    rows[[i]] <- tryCatch(
      prepare_one_task(tasks[i, , drop = FALSE]),
      error = function(e) {
        message(
          "FAILED PREP: ", CT, " ", tasks$gene[i], " ", tasks$cohort[i],
          " ", tasks$eqtl_source[i], " :: ", conditionMessage(e)
        )
        empty_manifest_row(
          tasks[i, , drop = FALSE],
          "ERROR",
          conditionMessage(e)
        )
      }
    )
  }

  m <- bind_rows(rows)
  ct_manifest <- file.path(
    CFG$output_dir,
    paste0("SMR_HEIDI.", CFG$analysis_group, ".", CT, ".manifest.tsv")
  )
  fwrite(m, ct_manifest, sep = "\t", quote = FALSE, na = "NA")
  message("Wrote manifest: ", ct_manifest)

  k <- k + 1L
  all_manifest[[k]] <- m
}

if (length(all_manifest)) {
  combined <- bind_rows(all_manifest)
  combined_file <- file.path(
    CFG$output_dir,
    paste0("SMR_HEIDI.", CFG$analysis_group, ".all_requested_tasks.manifest.tsv")
  )
  fwrite(combined, combined_file, sep = "\t", quote = FALSE, na = "NA")

  message("\nWrote combined manifest: ", combined_file)
  message("Prepared: ", sum(combined$status == "prepared", na.rm = TRUE), " / ", nrow(combined))
  message(
    "Prepared but missing PGEN prefix: ",
    sum(combined$status == "prepared_missing_pfile_prefix", na.rm = TRUE)
  )
  message(
    "Too few harmonized SNPs: ",
    sum(combined$status == "too_few_harmonized_overlap_snps", na.rm = TRUE)
  )
}
