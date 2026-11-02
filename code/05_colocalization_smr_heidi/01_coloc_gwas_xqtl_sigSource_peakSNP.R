#!/usr/bin/env Rscript

suppressMessages({
  library(data.table)
  library(dplyr)
  library(arrow)
  library(coloc)
})

# ==============================================================================
# Colocalization: significant peakSNP-TWAS genes x CRC GWAS x source-matched eQTL
#
# Usage:
#   Rscript coloc_gwas_xqtl_sigSource_peakSNP.R ABS coloc_peakSNP_NL.config.txt
#   Rscript coloc_gwas_xqtl_sigSource_peakSNP.R all coloc_peakSNP_NL.config.txt
#   Rscript coloc_gwas_xqtl_sigSource_peakSNP.R SSC coloc_peakSNP_ADSER.config.txt
#
# NL:
#   - read EUR/Multi significant genes from the CT-specific Cauchy-meta table
#   - EUR eQTL source comes from "P driven_EUR"
#   - Multi eQTL source comes from "P driven_Multi"
#
# ADSER:
#   - read EUR/Multi significant genes from the CT-specific ADSER table
#   - eQTL source is fixed to VUMC_HTAN_GSE266546
#
# TensorQTL variant convention expected:
#   rsID;chrom_pos_REF_ALT
# with slope / af oriented to ALT.
# ==============================================================================

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2L) {
  stop(
    "Usage: Rscript coloc_gwas_xqtl_sigSource_peakSNP.R ",
    "<CT|all> <config.R>"
  )
}

CT_ARG <- toupper(trimws(args[1]))
CONFIG_FILE <- args[2]

if (!file.exists(CONFIG_FILE)) {
  stop("Config file not found: ", CONFIG_FILE)
}

source(CONFIG_FILE)

if (!exists("COLOC_CONFIG")) {
  stop("COLOC_CONFIG was not defined in: ", CONFIG_FILE)
}

CFG <- COLOC_CONFIG

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

is_yes <- function(x) {
  tolower(trimws(as.character(x))) %in% c(
    "yes", "y", "true", "t", "1", "significant", "sig"
  )
}

first_nonempty <- function(x) {
  x <- as.character(x)
  x <- x[!is.na(x) & nzchar(trimws(x))]
  if (length(x)) x[1] else NA_character_
}

normalize_chr <- function(x) {
  x <- sub("^chr", "", as.character(x), ignore.case = TRUE)
  x
}

normalize_source_name <- function(x) {
  x0 <- trimws(as.character(x))
  out <- rep(NA_character_, length(x0))

  out[grepl("wellcome", x0, ignore.case = TRUE)] <- "WellcomeSanger"
  out[grepl("HTAN|VUMC", x0, ignore.case = TRUE)] <- "VUMC_HTAN_GSE266546"

  out
}

safe_numeric <- function(x) suppressWarnings(as.numeric(x))

ensure_dir <- function(x) {
  dir.create(x, recursive = TRUE, showWarnings = FALSE)
  x
}

# ------------------------------------------------------------------------------
# Read compact GENCODE gene -> chromosome map.
# Uses awk to avoid loading the full GTF into R.
# ------------------------------------------------------------------------------

read_gene_map <- function(gtf_file) {
  if (!file.exists(gtf_file)) {
    stop("Gene annotation not found: ", gtf_file)
  }

  cmd <- paste0(
    "awk 'BEGIN{FS=\"\\t\"; OFS=\"\\t\"} ",
    "$3==\"gene\" {print $1,$4,$5,$9}' ",
    shQuote(gtf_file)
  )

  g <- fread(
    cmd = cmd,
    header = FALSE,
    sep = "\t",
    data.table = FALSE,
    quote = ""
  )

  if (ncol(g) < 4L) {
    stop("Could not parse gene records from GTF: ", gtf_file)
  }

  names(g)[1:4] <- c("chr", "start", "end", "attribute")

  g$gene <- sub(
    '.*gene_id "([^"]+)".*',
    "\\1",
    g$attribute
  )

  g$gene_name <- ifelse(
    grepl('gene_name "', g$attribute, fixed = TRUE),
    sub('.*gene_name "([^"]+)".*', "\\1", g$attribute),
    NA_character_
  )

  g$gene <- canonical_gene_id(g$gene)
  g$chr <- normalize_chr(g$chr)

  # Convert chromosome once. GENCODE also contains X/Y/MT and alternate
  # contigs; suppressing coercion warnings here is intentional because those
  # records are removed immediately below.
  chr_num <- suppressWarnings(as.integer(g$chr))
  keep_autosome <- (
    !is.na(chr_num) &
      chr_num >= 1L &
      chr_num <= 22L &
      grepl("^ENSG[0-9]+$", g$gene)
  )

  g <- g[
    keep_autosome,
    c("gene", "gene_name", "chr", "start", "end"),
    drop = FALSE
  ]

  g <- g[!duplicated(g$gene), , drop = FALSE]

  message("GENCODE gene map: ", nrow(g), " chr1-22 genes")
  g
}

GENE_MAP <- read_gene_map(CFG$gene_annotation)

# ------------------------------------------------------------------------------
# Significant TWAS table
# ------------------------------------------------------------------------------

get_gene_name_from_sig_row <- function(row_df) {
  gene_name_cols <- grep(
    "^gene_name",
    names(row_df),
    value = TRUE
  )

  if (!length(gene_name_cols)) return(NA_character_)

  vals <- unlist(
    row_df[1, gene_name_cols, drop = FALSE],
    use.names = FALSE
  )

  first_nonempty(vals)
}

resolve_nl_source <- function(row_df, cohort) {
  source_col <- CFG$source_columns[[cohort]]

  if (!is.null(source_col) && source_col %in% names(row_df)) {
    src <- normalize_source_name(row_df[[source_col]][1])
    if (!is.na(src) && nzchar(src)) return(src)
  }

  # Fallback 1: META_SOURCE explicitly says one source only.
  meta_col <- paste0("META_SOURCE_", cohort)
  if (meta_col %in% names(row_df)) {
    meta <- as.character(row_df[[meta_col]][1])

    if (grepl("WellcomeSanger.*only", meta, ignore.case = TRUE)) {
      return("WellcomeSanger")
    }
    if (grepl("VUMC_HTAN_GSE266546.*only", meta, ignore.case = TRUE)) {
      return("VUMC_HTAN_GSE266546")
    }
  }

  # Fallback 2: if both datasets contributed but "P driven" is absent,
  # choose the source with the smaller source-specific TWAS p value.
  p_htan <- paste0("pvalue_HTAN_NL_", cohort)
  p_ws   <- paste0("pvalue_WellcomeSanger_NL_", cohort)

  if (all(c(p_htan, p_ws) %in% names(row_df))) {
    ph <- safe_numeric(row_df[[p_htan]][1])
    pw <- safe_numeric(row_df[[p_ws]][1])

    if (is.finite(ph) && is.finite(pw)) {
      if (ph <= pw) return("VUMC_HTAN_GSE266546")
      return("WellcomeSanger")
    }
    if (is.finite(ph)) return("VUMC_HTAN_GSE266546")
    if (is.finite(pw)) return("WellcomeSanger")
  }

  NA_character_
}

build_tasks_for_celltype <- function(CT) {
  sig_file <- fill_template(
    CFG$sig_file_template,
    list(CT = CT)
  )

  if (!file.exists(sig_file)) {
    stop("Significant-gene file not found for ", CT, ": ", sig_file)
  }

  sig <- fread(
    sig_file,
    data.table = FALSE,
    check.names = FALSE
  )

  required <- c(
    CFG$gene_column,
    CFG$sig_flag_columns[["EUR"]],
    CFG$sig_flag_columns[["Multi"]]
  )

  missing <- setdiff(required, names(sig))
  if (length(missing)) {
    stop(
      "Significant-gene file is missing required column(s): ",
      paste(missing, collapse = ", "),
      "\nFile: ", sig_file
    )
  }

  sig$gene_key <- canonical_gene_id(sig[[CFG$gene_column]])

  # One row per gene. If duplicates exist, preserve the first table row.
  sig <- sig[!duplicated(sig$gene_key), , drop = FALSE]

  map_idx <- match(sig$gene_key, GENE_MAP$gene)
  sig$chrom <- GENE_MAP$chr[map_idx]
  sig$gene_name_gtf <- GENE_MAP$gene_name[map_idx]

  tasks <- list()
  k <- 0L

  for (i in seq_len(nrow(sig))) {
    row_df <- sig[i, , drop = FALSE]
    gene <- row_df$gene_key[1]
    chrom <- row_df$chrom[1]

    if (is.na(chrom) || !chrom %in% as.character(1:22)) {
      message("Skipping gene without chr1-22 GENCODE mapping: ", gene)
      next
    }

    gene_name <- get_gene_name_from_sig_row(row_df)
    if (is.na(gene_name) || !nzchar(gene_name)) {
      gene_name <- row_df$gene_name_gtf[1]
    }

    for (cohort in c("EUR", "Multi")) {
      flag_col <- CFG$sig_flag_columns[[cohort]]

      if (!is_yes(row_df[[flag_col]][1])) next

      source <- if (identical(CFG$source_mode, "fixed")) {
        CFG$fixed_eqtl_source
      } else if (identical(CFG$source_mode, "p_driven")) {
        resolve_nl_source(row_df, cohort)
      } else {
        stop("Unknown source_mode: ", CFG$source_mode)
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

      k <- k + 1L

      p_col <- if (identical(CFG$analysis_group, "NL")) {
        paste0("META_P_", cohort)
      } else {
        paste0("pvalue_", cohort)
      }

      tasks[[k]] <- data.frame(
        cell_type = CT,
        gene = gene,
        gene_name = gene_name,
        chromosome = chrom,
        cohort = cohort,
        eqtl_source = source,
        TWAS_p = if (p_col %in% names(row_df)) {
          safe_numeric(row_df[[p_col]][1])
        } else {
          NA_real_
        },
        sig_file = sig_file,
        stringsAsFactors = FALSE
      )
    }
  }

  if (!length(tasks)) {
    return(data.frame())
  }

  bind_rows(tasks)
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

  cohort_cfg <- CFG$gwas[[cohort]]

  path <- fill_template(
    cohort_cfg$file_template,
    list(chrom = chrom)
  )

  if (!file.exists(path)) {
    stop("GWAS file not found: ", path)
  }

  g <- fread(
    path,
    data.table = FALSE,
    check.names = FALSE
  )

  req <- c(
    cohort_cfg$rsid_col,
    cohort_cfg$a1_col,
    cohort_cfg$a2_col,
    cohort_cfg$beta_col,
    cohort_cfg$se_col,
    cohort_cfg$p_col
  )

  miss <- setdiff(req, names(g))
  if (length(miss)) {
    stop(
      "GWAS file missing required column(s): ",
      paste(miss, collapse = ", "),
      "\nFile: ", path
    )
  }

  out <- data.frame(
    SNP = trimws(as.character(g[[cohort_cfg$rsid_col]])),
    A1 = toupper(trimws(as.character(g[[cohort_cfg$a1_col]]))),
    A2 = toupper(trimws(as.character(g[[cohort_cfg$a2_col]]))),
    BETA = safe_numeric(g[[cohort_cfg$beta_col]]),
    SE = safe_numeric(g[[cohort_cfg$se_col]]),
    P = safe_numeric(g[[cohort_cfg$p_col]]),
    FRQ = if (!is.null(cohort_cfg$freq_col) &&
              cohort_cfg$freq_col %in% names(g)) {
      safe_numeric(g[[cohort_cfg$freq_col]])
    } else {
      NA_real_
    },
    stringsAsFactors = FALSE
  )

  keep <- (
    grepl("^rs[0-9]+$", out$SNP, ignore.case = TRUE) &
      out$A1 %in% c("A", "C", "G", "T") &
      out$A2 %in% c("A", "C", "G", "T") &
      out$A1 != out$A2 &
      is.finite(out$BETA) &
      is.finite(out$SE) & out$SE > 0 &
      is.finite(out$P) & out$P > 0 & out$P <= 1
  )

  out <- out[keep, , drop = FALSE]
  out <- out[order(out$P), , drop = FALSE]
  out <- out[!duplicated(tolower(out$SNP)), , drop = FALSE]
  out$SNP <- tolower(out$SNP)

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

  chrom <- vapply(
    parts,
    function(z) if (length(z) >= 4L) z[1] else NA_character_,
    character(1)
  )
  pos <- vapply(
    parts,
    function(z) if (length(z) >= 4L) z[2] else NA_character_,
    character(1)
  )
  ref <- vapply(
    parts,
    function(z) if (length(z) >= 4L) z[3] else NA_character_,
    character(1)
  )
  alt <- vapply(
    parts,
    function(z) if (length(z) >= 4L) z[4] else NA_character_,
    character(1)
  )

  data.frame(
    SNP = tolower(trimws(rsid)),
    variant_chr = normalize_chr(chrom),
    variant_pos = safe_numeric(pos),
    REF = toupper(trimws(ref)),
    ALT = toupper(trimws(alt)),
    stringsAsFactors = FALSE
  )
}

read_eqtl_chr <- function(source, CT, chrom) {
  key <- paste(source, CT, chrom, sep = "|")
  if (exists(key, envir = EQTL_CACHE, inherits = FALSE)) {
    return(get(key, envir = EQTL_CACHE, inherits = FALSE))
  }

  source_cfg <- CFG$eqtl_sources[[source]]

  path <- fill_template(
    source_cfg$file_template,
    list(
      CT = CT,
      chrom = chrom
    )
  )

  if (!file.exists(path)) {
    stop(
      "eQTL parquet not found for source=", source,
      ", CT=", CT, ", chr=", chrom, ":\n", path
    )
  }

  # arrow::read_schema() is for Arrow IPC streams/Buffers, not a Parquet
  # filename.  Infer the Parquet schema through Dataset instead; this works
  # with a character file path and avoids reading the full table just to learn
  # its column names.
  ds <- arrow::open_dataset(path, format = "parquet")
  schema_names <- as.character(ds$schema$names)
  rm(ds)

  required <- c(
    "phenotype_id",
    "variant_id",
    "slope",
    "slope_se",
    "pval_nominal"
  )

  miss <- setdiff(required, schema_names)
  if (length(miss)) {
    stop(
      "eQTL parquet missing required column(s): ",
      paste(miss, collapse = ", "),
      "\nFile: ", path,
      "\nAvailable: ", paste(schema_names, collapse = ", ")
    )
  }

  optional <- intersect(
    c("tss_distance", "af", "ma_samples", "ma_count"),
    schema_names
  )

  eq <- as.data.frame(
    arrow::read_parquet(
      path,
      col_select = c(required, optional)
    )
  )

  eq$phenotype_key <- canonical_gene_id(eq$phenotype_id)

  parsed <- parse_tensorqtl_variant(eq$variant_id)
  eq <- cbind(eq, parsed)

  eq$BETA <- safe_numeric(eq$slope)
  eq$SE <- safe_numeric(eq$slope_se)
  eq$P <- safe_numeric(eq$pval_nominal)
  eq$AF <- if ("af" %in% names(eq)) safe_numeric(eq$af) else NA_real_

  keep <- (
    grepl("^rs[0-9]+$", eq$SNP, ignore.case = TRUE) &
      eq$REF %in% c("A", "C", "G", "T") &
      eq$ALT %in% c("A", "C", "G", "T") &
      eq$REF != eq$ALT &
      is.finite(eq$BETA) &
      is.finite(eq$SE) & eq$SE > 0 &
      is.finite(eq$P) & eq$P > 0 & eq$P <= 1
  )

  if ("tss_distance" %in% names(eq) &&
      is.finite(CFG$flank_length) &&
      CFG$flank_length > 0) {
    keep <- keep &
      is.finite(safe_numeric(eq$tss_distance)) &
      abs(safe_numeric(eq$tss_distance)) <= CFG$flank_length
  }

  eq <- eq[keep, , drop = FALSE]
  eq <- eq[order(eq$P), , drop = FALSE]
  eq <- eq[!duplicated(paste(eq$phenotype_key, eq$SNP, sep = "|")), , drop = FALSE]

  assign(key, eq, envir = EQTL_CACHE)

  message(
    source, " ", CT, " chr", chrom,
    " eQTL rows loaded: ", nrow(eq)
  )

  eq
}

resolve_eqtl_N <- function(eq_gene, source_cfg) {
  if (!is.null(source_cfg$N) &&
      length(source_cfg$N) == 1L &&
      is.finite(source_cfg$N) &&
      source_cfg$N > 2) {
    return(as.integer(round(source_cfg$N)))
  }

  # Auto-estimate N from tensorQTL ma_count + af when an explicit N is not
  # configured. ma_count is minor-allele count.
  if (all(c("ma_count", "af") %in% names(eq_gene))) {
    af <- safe_numeric(eq_gene$af)
    mac <- safe_numeric(eq_gene$ma_count)
    maf <- pmin(af, 1 - af)

    n_est <- mac / (2 * maf)
    n_est <- n_est[
      is.finite(n_est) &
        n_est > 2 &
        abs(n_est - round(n_est)) < 0.25
    ]

    if (length(n_est)) {
      return(as.integer(round(stats::median(n_est))))
    }
  }

  stop(
    "Cannot determine eQTL sample size. ",
    "Set eqtl_sources[[source]]$N in the config, or retain tensorQTL ",
    "'af' + 'ma_count' columns for automatic inference."
  )
}

# ------------------------------------------------------------------------------
# Allele harmonization
# ------------------------------------------------------------------------------

complement_allele <- function(x) {
  chartr("ACGT", "TGCA", toupper(x))
}

is_ambiguous_pair <- function(a1, a2) {
  (a1 == "A" & a2 == "T") |
    (a1 == "T" & a2 == "A") |
    (a1 == "C" & a2 == "G") |
    (a1 == "G" & a2 == "C")
}

harmonize_gwas_eqtl <- function(gwas, eqtl) {
  x <- merge(
    gwas,
    eqtl,
    by = "SNP",
    suffixes = c("_g", "_x")
  )

  if (!nrow(x)) return(x)

  # GWAS effect allele = A1_g
  # tensorQTL eQTL effect allele = ALT
  a1g <- toupper(x$A1)
  a2g <- toupper(x$A2)
  a1x <- toupper(x$ALT)
  a2x <- toupper(x$REF)

  ambiguous <- is_ambiguous_pair(a1g, a2g)

  same <- a1x == a1g & a2x == a2g
  flip <- a1x == a2g & a2x == a1g

  comp_same <- (
    a1x == complement_allele(a1g) &
      a2x == complement_allele(a2g)
  )

  comp_flip <- (
    a1x == complement_allele(a2g) &
      a2x == complement_allele(a1g)
  )

  matched <- same | flip | comp_same | comp_flip

  keep <- matched
  if (isTRUE(CFG$drop_strand_ambiguous)) {
    keep <- keep & !ambiguous
  }

  x <- x[keep, , drop = FALSE]

  if (!nrow(x)) return(x)

  # Recalculate on retained rows.
  a1g <- toupper(x$A1)
  a2g <- toupper(x$A2)
  a1x <- toupper(x$ALT)
  a2x <- toupper(x$REF)

  flip <- a1x == a2g & a2x == a1g
  comp_flip <- (
    a1x == complement_allele(a2g) &
      a2x == complement_allele(a1g)
  )

  x$BETA_x_aligned <- x$BETA_x
  x$BETA_x_aligned[flip | comp_flip] <-
    -x$BETA_x_aligned[flip | comp_flip]

  # Prefer tensorQTL AF for MAF; fall back to numeric GWAS FRQ if needed.
  maf_x <- pmin(x$AF, 1 - x$AF)
  maf_g <- pmin(x$FRQ, 1 - x$FRQ)

  x$maf <- ifelse(
    is.finite(maf_x) & maf_x > 0 & maf_x < 0.5,
    maf_x,
    maf_g
  )

  x <- x[
    is.finite(x$maf) &
      x$maf > 0 &
      x$maf <= 0.5 &
      is.finite(x$BETA_g) &
      is.finite(x$SE_g) &
      is.finite(x$P_g) &
      is.finite(x$BETA_x_aligned) &
      is.finite(x$SE_x) &
      is.finite(x$P_x),
    ,
    drop = FALSE
  ]

  if (!nrow(x)) return(x)

  x <- x[order(x$P_g, x$P_x), , drop = FALSE]
  x <- x[!duplicated(x$SNP), , drop = FALSE]

  x
}

# ------------------------------------------------------------------------------
# One coloc task
# ------------------------------------------------------------------------------

run_one_coloc <- function(task_row, CT_OUTDIR) {
  gene <- task_row$gene
  cohort <- task_row$cohort
  source <- task_row$eqtl_source
  chrom <- task_row$chromosome

  message(
    "\n------------------------------------------------------------\n",
    "CT=", task_row$cell_type,
    " gene=", gene,
    " cohort=", cohort,
    " eQTL=", source,
    " chr=", chrom
  )

  gwas <- read_gwas_chr(cohort, chrom)
  eq_chr <- read_eqtl_chr(source, task_row$cell_type, chrom)

  eq_gene <- eq_chr[
    eq_chr$phenotype_key == gene,
    ,
    drop = FALSE
  ]

  if (!nrow(eq_gene)) {
    return(data.frame(
      task_row,
      qtl_N = NA_integer_,
      n_eqtl_variants = 0L,
      n_overlap_harmonized = 0L,
      nsnps = NA_real_,
      PP.H0.abf = NA_real_,
      PP.H1.abf = NA_real_,
      PP.H2.abf = NA_real_,
      PP.H3.abf = NA_real_,
      PP.H4.abf = NA_real_,
      status = "gene_not_found_in_eqtl",
      stringsAsFactors = FALSE,
      check.names = FALSE
    ))
  }

  source_cfg <- CFG$eqtl_sources[[source]]
  qtl_N <- resolve_eqtl_N(eq_gene, source_cfg)

  matched <- harmonize_gwas_eqtl(gwas, eq_gene)
  n_overlap <- nrow(matched)

  if (n_overlap < CFG$min_overlap_snps) {
    message(
      "Skipping coloc: ", n_overlap,
      " harmonized SNPs < min_overlap_snps=", CFG$min_overlap_snps
    )

    return(data.frame(
      task_row,
      qtl_N = qtl_N,
      n_eqtl_variants = nrow(eq_gene),
      n_overlap_harmonized = n_overlap,
      nsnps = NA_real_,
      PP.H0.abf = NA_real_,
      PP.H1.abf = NA_real_,
      PP.H2.abf = NA_real_,
      PP.H3.abf = NA_real_,
      PP.H4.abf = NA_real_,
      status = "too_few_overlap_snps",
      stringsAsFactors = FALSE,
      check.names = FALSE
    ))
  }

  gwas_cfg <- CFG$gwas[[cohort]]

  ds1 <- list(
    beta = matched$BETA_g,
    varbeta = matched$SE_g^2,
    pvalues = matched$P_g,
    snp = matched$SNP,
    type = "cc",
    s = gwas_cfg$case_fraction
  )

  if (!is.null(gwas_cfg$N) &&
      is.finite(gwas_cfg$N) &&
      gwas_cfg$N > 0) {
    ds1$N <- as.integer(round(gwas_cfg$N))
  }

  ds2 <- list(
    beta = matched$BETA_x_aligned,
    varbeta = matched$SE_x^2,
    pvalues = matched$P_x,
    snp = matched$SNP,
    type = "quant",
    N = qtl_N
  )

  coloc_fit <- coloc::coloc.abf(
    dataset1 = ds1,
    dataset2 = ds2,
    MAF = matched$maf
  )

  s <- coloc_fit$summary

  summary_row <- data.frame(
    task_row,
    qtl_N = qtl_N,
    n_eqtl_variants = nrow(eq_gene),
    n_overlap_harmonized = n_overlap,
    nsnps = unname(s["nsnps"]),
    PP.H0.abf = unname(s["PP.H0.abf"]),
    PP.H1.abf = unname(s["PP.H1.abf"]),
    PP.H2.abf = unname(s["PP.H2.abf"]),
    PP.H3.abf = unname(s["PP.H3.abf"]),
    PP.H4.abf = unname(s["PP.H4.abf"]),
    status = "coloc_completed",
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  # Save variant-level details only for reasonably supported H4 results.
  if (is.finite(summary_row$PP.H4.abf) &&
      summary_row$PP.H4.abf >= CFG$detail_pph4_threshold) {

    detail_dir <- ensure_dir(
      file.path(
        CT_OUTDIR,
        "details",
        cohort,
        source
      )
    )

    detail <- as.data.frame(coloc_fit$results)
    detail$cell_type <- task_row$cell_type
    detail$gene <- gene
    detail$gene_name <- task_row$gene_name
    detail$cohort <- cohort
    detail$eqtl_source <- source

    detail_file <- file.path(
      detail_dir,
      paste0(gene, ".coloc.details.csv")
    )

    fwrite(detail, detail_file)
    message("Wrote H4 details: ", detail_file)
  }

  summary_row
}

# ------------------------------------------------------------------------------
# Run one or all cell types
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

all_summaries <- list()
all_tasks <- list()

for (CT in CTs) {
  message(
    "\n======================================================================\n",
    "COLOC ", CFG$analysis_group, " ", CT, "\n",
    "======================================================================"
  )

  CT_OUTDIR <- ensure_dir(
    file.path(CFG$output_dir, CT)
  )

  tasks <- build_tasks_for_celltype(CT)

  task_file <- file.path(
    CT_OUTDIR,
    paste0(
      "coloc.", CFG$analysis_group, ".", CT,
      ".significant_gene_tasks.csv"
    )
  )
  fwrite(tasks, task_file)
  message("Wrote tasks: ", task_file)

  if (!nrow(tasks)) {
    warning("No coloc tasks generated for ", CT)
    next
  }

  summary_list <- vector("list", nrow(tasks))

  for (i in seq_len(nrow(tasks))) {
    summary_list[[i]] <- tryCatch(
      run_one_coloc(tasks[i, , drop = FALSE], CT_OUTDIR),
      error = function(e) {
        message(
          "FAILED: ", CT, " ", tasks$gene[i], " ",
          tasks$cohort[i], " ", tasks$eqtl_source[i],
          " :: ", conditionMessage(e)
        )

        data.frame(
          tasks[i, , drop = FALSE],
          qtl_N = NA_integer_,
          n_eqtl_variants = NA_integer_,
          n_overlap_harmonized = NA_integer_,
          nsnps = NA_real_,
          PP.H0.abf = NA_real_,
          PP.H1.abf = NA_real_,
          PP.H2.abf = NA_real_,
          PP.H3.abf = NA_real_,
          PP.H4.abf = NA_real_,
          status = paste0("ERROR: ", conditionMessage(e)),
          stringsAsFactors = FALSE,
          check.names = FALSE
        )
      }
    )
  }

  summary_df <- bind_rows(summary_list) |>
    arrange(
      cohort,
      desc(PP.H4.abf),
      gene
    )

  summary_file <- file.path(
    CT_OUTDIR,
    paste0(
      "coloc.", CFG$analysis_group, ".", CT,
      ".GWAS_eQTL_sigSource.summary.csv"
    )
  )

  fwrite(summary_df, summary_file)
  message("Wrote summary: ", summary_file)

  all_summaries[[CT]] <- summary_df
  all_tasks[[CT]] <- tasks
}

if (length(all_summaries)) {
  combined <- bind_rows(all_summaries)

  combined_file <- file.path(
    CFG$output_dir,
    paste0(
      "coloc.", CFG$analysis_group,
      ".all_requested_celltypes.GWAS_eQTL_sigSource.summary.csv"
    )
  )

  fwrite(combined, combined_file)
  message("\nWrote combined summary: ", combined_file)

  if (nrow(combined)) {
    message(
      "Completed coloc tasks: ",
      sum(combined$status == "coloc_completed", na.rm = TRUE),
      " / ", nrow(combined)
    )
    message(
      "PP.H4 >= 0.5: ",
      sum(
        combined$status == "coloc_completed" &
          is.finite(combined$PP.H4.abf) &
          combined$PP.H4.abf >= 0.5,
        na.rm = TRUE
      )
    )
  }
}
