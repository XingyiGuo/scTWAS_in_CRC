#!/usr/bin/env Rscript

# ==============================================================================
# Normal colon (NL), 7 cell types.
#
#   1) Every model-build SNP of every significant TWAS gene is read from the
#      cell-type specific prediction weights (one weight file per chromosome).
#      Normal colon has two models; "META_P driven" in the NL_Multi_TWAS sheet
#      selects the one used for each gene x cell type:
#         WellcomeSanger      -> WellcomeSanger_NL_{CT}_peakSNP/...weights.txt
#         VUMC_HTAN_GSE266546 -> VUMC_HTAN_GSE266546_NL_{CT}_peakSNP/...weights.txt
#   2) Each of those SNPs is annotated with the FULL GWAS record from the CRC
#      meta-analysis. Sig_resource of the gene x cell type (SupplementaryTables
#      sheet NL_Multi_TWAS) selects the summary statistics:
#         EUR_sig              -> EUR GWAS
#         EUR_notSig_Multi_Sig -> multi-ancestry GWAS
#   3) ONE table is written holding every model SNP with all of its GWAS
#      columns; the column is_smallest_GWAS_P flags, within each cell type x
#      gene, the SNP with the smallest GWAS P value (the lead variant).
#   4) Those lead variants drive the two TF summaries that are plotted:
#      number of TWAS genes whose lead variant is occupied by each CRC TF, and
#      number of TWAS genes correlated with each CRC TF above the threshold.
#
# One plot per cell type, plus one combined 7-panel plot.
# ==============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(grid)
})

# ------------------------------------------------------------------------------
# 0. CONFIGURATION
# ------------------------------------------------------------------------------
DATASET_LABEL <- "NL"
CELL_TYPES <- c("ABS", "GOB", "CT", "EE", "TUF", "TAC", "STM")

# --- Supplementary workbook: significant gene x cell type list, Sig_resource,
# --- and the prediction model that drove the meta-analysis P value ------------
SUPP_XLSX      <- "/data/l2_bioinfo1/liq17/scTWAS/ms_figures_tables/Tables/SupplementaryTables.xlsx"
SUPP_SHEET     <- "NL_Multi_TWAS"
SUPP_SKIP      <- 2L                 # two grouping rows sit above the real header
SUPP_MODEL_COL <- "META_P driven"    # WellcomeSanger | VUMC_HTAN_GSE266546

# --- Cell-type specific prediction weights: source of the model-build SNPs ----
# Normal colon has two models; the one named in "META_P driven" is used for each
# gene x cell type. {CT} = cell type, {CHR} = chromosome.
WEIGHT_ROOT <- "/data/l2_bioinfo1/liq17/scTWAS/TWAS_HTAN_GSE266546_WellcomeSanger_res"

MODEL_SOURCE_SINGLE <- NA_character_

WEIGHT_TMPL <- list(
  "VUMC_HTAN_GSE266546" = "VUMC_HTAN_GSE266546_NL_{CT}_peakSNP/VUMC_HTAN_GSE266546_NL_{CT}_peakSNP_chr{CHR}_weights.txt",
  "WellcomeSanger"      = "WellcomeSanger_NL_{CT}_peakSNP/WellcomeSanger_NL_{CT}_peakSNP_chr{CHR}_weights.txt"
)

# --- Cistrome peak overlaps of the model SNPs (TF occupancy) ------------------
OVERLAP_FILE <- "/v5000/data/l2_bioinfo1/liq17/scTWAS/TWAS_HTAN_GSE266546_WellcomeSanger_res/SPrediXcan_peakSNP_association/Cauchy_Meta_NL_WellcomeSanger_NL_EUR_Multi/CRC51TFs_occupy_model_snps/HTAN_NL_WellcomeSanger_NL_7CTs_sigGenes_modelSNPs_CistromePeakOverlaps_long.csv"

# --- GWAS summary statistics --------------------------------------------------
# Sig_resource of the gene x cell type decides which meta-analysis is queried:
#   EUR_sig              -> EUR GWAS
#   EUR_notSig_Multi_Sig -> multi-ancestry GWAS
GWAS_CFG <- list(
  EUR = list(
    dir  = "/data/l2_bioinfo1/GWAS_SS/BIOVU/CRC/metal/GC_Final/processed/EUR/rsid_annotated_CHRPOS",
    tmpl = "final_gc_maf0.01_CRC_EUR_META_1.chr{CHR}.chrpos.rsid.assoc.dosage.gz"
  ),
  Multi = list(
    dir  = "/data/l2_bioinfo1/GWAS_SS/BIOVU/CRC/metal/GC_Final/processed/Multi/rsid_annotated_CHRPOS",
    tmpl = "final_gc_maf0.01_CRC_Multi_META_1.chr{CHR}.chrpos.rsid.assoc.dosage.gz"
  )
)

SIG_RESOURCE_TO_ANCESTRY <- c(
  "EUR_sig"              = "EUR",
  "EUR_notSig_Multi_Sig" = "Multi"
)

# Chromosomes with weight / GWAS files. Add "X" if the models cover it.
CHROMS <- as.character(1:22)

# Directory containing files such as:
#   ABS_CRC51TF_x_NL_TWASriskGene_PearsonR_matrix.csv, ...
COR_DIR <- "/v5000/data/l2_bioinfo1/liq17/scTWAS/CellOracleR_res/WellcomeSanger_NL_TF_TWASgene_pseudobulk_correlation/"

OUTDIR <- "/v5000/data/l2_bioinfo1/liq17/scTWAS/CellOracleR_res/WellcomeSanger_TF_leadVariant_vs_correlation_summary/"
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

COR_THRESHOLD <- 0.5

# If FALSE, TF == risk gene correlations (e.g. CEBPB vs CEBPB = 1)
# are retained, exactly following "corr > 0.5".
EXCLUDE_SELF_CORRELATION <- FALSE

# ------------------------------------------------------------------------------
# 1. HELPERS
# ------------------------------------------------------------------------------

cor_file_for_ct <- function(ct) {

  # Preferred filename based on the analysis design.
  preferred <- file.path(
    COR_DIR,
    paste0(
      ct,"/", ct, "_CRC51TF_x_NL_TWASriskGene_PearsonR_matrix.csv"
    )
  )

  # Fallback path using the same NL correlation-matrix naming convention.
  fallback <- file.path(
    COR_DIR,
    paste0(
      ct,"/", ct, "_CRC51TF_x_NL_TWASriskGene_PearsonR_matrix.csv"
    )
  )

  if (file.exists(preferred)) {
    return(preferred)
  }

  if (file.exists(fallback)) {
    return(fallback)
  }

  stop(
    "Correlation matrix not found for ", ct, ". Tried:\n  ",
    preferred, "\n  ", fallback
  )
}


complete_tf_grid <- function(count_dt, all_tfs) {

  grid <- CJ(
    Cells = CELL_TYPES,
    TFs = all_tfs,
    unique = TRUE
  )

  out <- merge(
    grid,
    count_dt,
    by = c("Cells", "TFs"),
    all.x = TRUE,
    sort = FALSE
  )

  out
}


# ------------------------------------------------------------------------------
# 1b. HELPERS: MODEL WEIGHTS, SUPPLEMENTARY SHEET, GWAS SUMMARY STATISTICS
# ------------------------------------------------------------------------------

# First name in `nms` matching one of `exact` (case-insensitive), else the first
# matching `regex`, else NA.
pick_col <- function(nms, exact, regex = NULL) {

  hit <- which(tolower(trimws(nms)) %in% tolower(exact))

  if (length(hit)) {
    return(nms[hit[1L]])
  }

  if (!is.null(regex)) {
    hit <- grep(regex, nms, ignore.case = TRUE)
    if (length(hit)) {
      return(nms[hit[1L]])
    }
  }

  NA_character_
}


norm_chr <- function(x) {

  x <- toupper(trimws(as.character(x)))
  x <- sub("^CHR", "", x)
  x[x == "23"] <- "X"
  x[x == "24"] <- "Y"
  x[is.na(x) | !nzchar(x) | x %in% c("NA", "NULL", ".")] <- NA_character_

  x
}


clean_pos <- function(x) {

  x <- trimws(as.character(x))
  x <- sub("\\.0+$", "", x)
  x[is.na(x) | !grepl("^[0-9]+$", x)] <- NA_character_

  x
}


clean_rs <- function(x) {

  x <- tolower(trimws(as.character(x)))
  x[is.na(x) | !nzchar(x) | x %in% c("na", ".", "none", "null")] <- NA_character_

  x
}


add_variant_keys <- function(d) {

  d[, key_rs := clean_rs(rs)]
  d[, key_cp := fifelse(
    !is.na(chr) & !is.na(pos),
    paste0(chr, ":", pos),
    NA_character_
  )]
  d[, variant_key := fifelse(!is.na(key_cp), key_cp, key_rs)]

  d[]
}


# fread reads .gz directly when R.utils is installed; otherwise pipe through
# gzip, and fall back to a gzfile() connection if neither is available.
read_maybe_gz <- function(f) {

  attempt <- function(expr) tryCatch(expr, error = function(e) NULL)

  out <- attempt(
    fread(f, data.table = TRUE, check.names = FALSE, showProgress = FALSE)
  )

  if ((is.null(out) || !nrow(out)) && nzchar(Sys.which("gzip"))) {
    out <- attempt(
      fread(
        cmd = paste("gzip -dc", shQuote(f)),
        data.table = TRUE,
        check.names = FALSE,
        showProgress = FALSE
      )
    )
  }

  if (is.null(out) || !nrow(out)) {
    con <- gzfile(f, "rt")
    on.exit(close(con), add = TRUE)
    out <- attempt(
      fread(
        text = readLines(con),
        data.table = TRUE,
        check.names = FALSE,
        showProgress = FALSE
      )
    )
  }

  if (is.null(out)) {
    stop("Could not read ", f)
  }

  out
}


# ---- supplementary workbook ---------------------------------------------------

read_supp_sheet <- function(path, sheet, skip) {

  if (!file.exists(path)) {
    stop("SupplementaryTables workbook not found: ", path)
  }

  if (requireNamespace("readxl", quietly = TRUE)) {

    x <- as.data.table(
      readxl::read_excel(
        path,
        sheet = sheet,
        skip = skip,
        .name_repair = "minimal"
      )
    )

  } else if (requireNamespace("openxlsx", quietly = TRUE)) {

    x <- as.data.table(
      openxlsx::read.xlsx(
        path,
        sheet = sheet,
        startRow = skip + 1L,
        check.names = FALSE
      )
    )

  } else {
    stop("Neither readxl nor openxlsx is available; cannot read ", path)
  }

  x
}


# One row per significant gene x cell type, carrying the GWAS ancestry implied
# by Sig_resource and (optionally) the prediction model implied by `model_col`.
build_gene_table <- function(path, sheet, skip, cell_types,
                             model_col = NA_character_,
                             single_model_source = NA_character_) {

  x <- read_supp_sheet(path, sheet, skip)

  need <- c("cell_type", "gene", "gene_name", "Sig_resource")

  if (!is.na(model_col)) {
    need <- c(need, model_col)
  }

  miss <- setdiff(need, names(x))

  if (length(miss)) {
    stop(
      "Sheet ", sheet, " is missing column(s): ", paste(miss, collapse = ", "),
      "\n  columns seen: ", paste(names(x), collapse = ", "),
      "\n  (SUPP_SKIP is currently ", skip, ")"
    )
  }

  x <- x[!is.na(gene) & nzchar(trimws(as.character(gene)))]

  out <- data.table(
    CellType = trimws(as.character(x$cell_type)),
    gene = sub("\\..*$", "", trimws(as.character(x$gene))),
    gene_name = trimws(as.character(x$gene_name)),
    Sig_resource = trimws(as.character(x$Sig_resource)),
    model_source = if (!is.na(model_col)) {
      trimws(as.character(x[[model_col]]))
    } else {
      single_model_source
    }
  )

  out[, GWAS_ancestry := unname(SIG_RESOURCE_TO_ANCESTRY[Sig_resource])]

  bad <- unique(out$Sig_resource[is.na(out$GWAS_ancestry)])

  if (length(bad)) {
    stop(
      "Unrecognised Sig_resource value(s) in ", sheet, ": ",
      paste(bad, collapse = ", "),
      "\n  extend SIG_RESOURCE_TO_ANCESTRY to cover them."
    )
  }

  bad_model <- setdiff(unique(out$model_source), names(WEIGHT_TMPL))

  if (length(bad_model)) {
    stop(
      "No weight-file template for model source(s): ",
      paste(bad_model, collapse = ", "),
      "\n  extend WEIGHT_TMPL to cover them."
    )
  }

  out <- out[CellType %chin% cell_types]

  clash <- out[, .(n = uniqueN(paste(Sig_resource, model_source))),
               by = .(CellType, gene)][n > 1L]

  if (nrow(clash)) {
    stop(
      "Conflicting Sig_resource / model source for cell type x gene in ", sheet,
      ": ", paste(clash$CellType, clash$gene, sep = "/", collapse = ", ")
    )
  }

  unique(out, by = c("CellType", "gene"))
}


# ---- cell-type-specific prediction weights ------------------------------------

WEIGHT_COL_CANDIDATES <- list(
  gene       = c("gene", "gene_id", "geneid", "ensembl_gene_id"),
  rsid       = c("rsid", "rs_id", "rs", "RSID", "rsID"),
  varID      = c("varID", "var_id", "varid", "variant_id", "snpID", "snp_id",
                 "SNP", "snp", "marker"),
  weight     = c("weight", "beta", "w", "effect_size"),
  ref_allele = c("ref_allele", "ref", "a2", "A2", "non_effect_allele"),
  eff_allele = c("eff_allele", "effect_allele", "alt", "a1", "A1")
)


weight_path <- function(model_source, ct, chr) {

  tmpl <- WEIGHT_TMPL[[model_source]]

  if (is.null(tmpl)) {
    stop("No weight-file template for model source: ", model_source)
  }

  file.path(
    WEIGHT_ROOT,
    gsub("{CHR}", chr, gsub("{CT}", ct, tmpl, fixed = TRUE), fixed = TRUE)
  )
}


# Tolerate small deviations in the weight-file naming.
find_weight_file <- function(model_source, ct, chr) {

  f <- weight_path(model_source, ct, chr)

  if (file.exists(f)) {
    return(f)
  }

  d <- dirname(f)

  if (!dir.exists(d)) {
    return(NA_character_)
  }

  hit <- list.files(
    d,
    pattern = paste0("(^|[^0-9])chr", chr, "([^0-9].*)?weight.*\\.(txt|tsv|csv)(\\.gz)?$"),
    full.names = TRUE,
    ignore.case = TRUE
  )

  if (length(hit)) hit[1L] else NA_character_
}


# Every model-build SNP of every significant gene, read from the cell-type
# specific weight files (one file per chromosome).
read_model_snps <- function(gene_tbl) {

  jobs <- unique(gene_tbl[, .(CellType, model_source)])
  out <- list()

  for (j in seq_len(nrow(jobs))) {

    ct <- jobs$CellType[j]
    ms <- jobs$model_source[j]

    want <- gene_tbl[CellType == ct & model_source == ms, unique(gene)]
    n_files <- 0L

    for (cc in CHROMS) {

      f <- find_weight_file(ms, ct, cc)

      if (is.na(f)) {
        next
      }

      n_files <- n_files + 1L

      w <- fread(f, data.table = TRUE, check.names = FALSE, showProgress = FALSE)

      if (!nrow(w)) {
        next
      }

      sel <- vapply(
        WEIGHT_COL_CANDIDATES,
        function(cand) pick_col(names(w), cand),
        character(1)
      )

      if (is.na(sel[["gene"]]) || (is.na(sel[["rsid"]]) && is.na(sel[["varID"]]))) {
        stop(
          "Weight file lacks a gene and/or rsID column:\n  ", f,
          "\n  columns seen: ", paste(names(w), collapse = ", ")
        )
      }

      d <- data.table(
        CellType = ct,
        model_source = ms,
        gene = sub("\\..*$", "", trimws(as.character(w[[sel[["gene"]]]]))),
        rs = if (!is.na(sel[["rsid"]])) {
          trimws(as.character(w[[sel[["rsid"]]]]))
        } else {
          NA_character_
        },
        varID = if (!is.na(sel[["varID"]])) {
          trimws(as.character(w[[sel[["varID"]]]]))
        } else {
          NA_character_
        },
        model_weight = if (!is.na(sel[["weight"]])) {
          suppressWarnings(as.numeric(w[[sel[["weight"]]]]))
        } else {
          NA_real_
        },
        model_ref_allele = if (!is.na(sel[["ref_allele"]])) {
          as.character(w[[sel[["ref_allele"]]]])
        } else {
          NA_character_
        },
        model_eff_allele = if (!is.na(sel[["eff_allele"]])) {
          as.character(w[[sel[["eff_allele"]]]])
        } else {
          NA_character_
        },
        weight_file = basename(f),
        chr_from_file = as.character(cc)
      )

      d <- d[gene %chin% want]

      if (nrow(d)) {
        out[[length(out) + 1L]] <- d
      }
    }

    if (!n_files) {
      stop(
        "No weight files found for model source '", ms, "', cell type ", ct,
        "\n  expected e.g. ", weight_path(ms, ct, CHROMS[1L]),
        "\n  adjust WEIGHT_ROOT / WEIGHT_TMPL."
      )
    }

    message(
      "Model SNPs read: ", ms, " / ", ct, " - ", n_files, " weight file(s)"
    )
  }

  res <- rbindlist(out, use.names = TRUE, fill = TRUE)

  if (!nrow(res)) {
    stop("No model SNPs were recovered from the weight files.")
  }

  has_cp <- grepl("^(chr)?[0-9XYMT]+[_:][0-9]+", res$varID, ignore.case = TRUE)

  res[, chr := norm_chr(
    fifelse(
      has_cp,
      sub("^((chr)?[0-9XYMT]+)[_:].*$", "\\1", varID, ignore.case = TRUE),
      chr_from_file
    )
  )]

  res[, pos := clean_pos(
    fifelse(
      has_cp,
      sub("^(chr)?[0-9XYMT]+[_:]([0-9]+).*$", "\\2", varID, ignore.case = TRUE),
      NA_character_
    )
  )]

  res <- add_variant_keys(res)

  dropped <- sum(is.na(res$variant_key))

  if (dropped) {
    warning(dropped, " model SNP row(s) had neither rsID nor chr:pos and were dropped.",
            call. = FALSE)
  }

  unique(
    res[!is.na(variant_key)],
    by = c("CellType", "model_source", "gene", "variant_key")
  )
}


# ---- GWAS summary statistics ---------------------------------------------------

# Every column of the GWAS file is kept; only the columns needed to build the
# match keys and to rank by significance are identified by name.
detect_gwas_cols <- function(nms) {

  list(
    rs = pick_col(
      nms,
      c("rsid", "rsID", "rs", "rs_id", "snp_rsid", "rsids", "RS"),
      "^rs[._]?id$"
    ),
    chr = pick_col(
      nms,
      c("chr", "CHR", "chrom", "chromosome", "#CHROM", "CHROM", "chr_name")
    ),
    pos = pick_col(
      nms,
      c("pos", "POS", "bp", "BP", "position", "base_pair_location", "bp_hg38",
        "bp_hg19")
    ),
    chrpos = pick_col(
      nms,
      c("chrpos", "chr_pos", "CHRPOS", "snpid", "SNPID", "MarkerName",
        "markername", "ID", "variant_id", "SNP", "snp")
    ),
    p = pick_col(
      nms,
      c("P", "p", "pval", "p_value", "P_value", "PVALUE", "P.value", "P-value",
        "Pvalue", "pvalue", "P_META"),
      "^p[-_.]?val"
    )
  )
}


gwas_path <- function(ancestry, chr) {

  cfg <- GWAS_CFG[[ancestry]]

  if (is.null(cfg)) {
    stop("No GWAS configuration for ancestry: ", ancestry)
  }

  file.path(cfg$dir, gsub("{CHR}", chr, cfg$tmpl, fixed = TRUE))
}


collect_gwas_hits <- function(need_rs, need_cp, chrs, ancestry) {

  if (!length(chrs)) {
    chrs <- CHROMS
  }

  out <- vector("list", length(chrs))

  for (i in seq_along(chrs)) {

    f <- gwas_path(ancestry, chrs[i])

    if (!file.exists(f)) {
      warning("GWAS file not found, skipped: ", f, call. = FALSE)
      next
    }

    message("  reading ", ancestry, " chr", chrs[i], ": ", basename(f))

    g <- read_maybe_gz(f)
    cols <- detect_gwas_cols(names(g))

    if (is.na(cols$p)) {
      stop(
        "Could not identify a P-value column in\n  ", f,
        "\n  columns seen: ", paste(names(g), collapse = ", ")
      )
    }

    key_rs <- if (!is.na(cols$rs)) {
      clean_rs(g[[cols$rs]])
    } else {
      rep(NA_character_, nrow(g))
    }

    if (!is.na(cols$chr) && !is.na(cols$pos)) {

      key_cp <- paste0(norm_chr(g[[cols$chr]]), ":", clean_pos(g[[cols$pos]]))

    } else if (!is.na(cols$chrpos)) {

      raw <- trimws(as.character(g[[cols$chrpos]]))
      ok <- grepl("^(chr)?[0-9XYMT]+[:_-][0-9]+", raw, ignore.case = TRUE)

      key_cp <- ifelse(
        ok,
        paste0(
          norm_chr(sub("^((chr)?[0-9XYMT]+)[:_-].*$", "\\1", raw, ignore.case = TRUE)),
          ":",
          sub("^(chr)?[0-9XYMT]+[:_-]([0-9]+).*$", "\\2", raw, ignore.case = TRUE)
        ),
        NA_character_
      )

    } else {
      key_cp <- rep(NA_character_, nrow(g))
    }

    key_cp[grepl("NA", key_cp, fixed = TRUE)] <- NA_character_

    pv <- suppressWarnings(as.numeric(g[[cols$p]]))

    setnames(g, names(g), paste0("GWAS_", names(g)))

    g <- cbind(
      data.table(
        key_rs = key_rs,
        key_cp = key_cp,
        GWAS_P_used = pv,
        GWAS_file = basename(f)
      ),
      g
    )

    keep <- g[
      (!is.na(key_rs) & key_rs %chin% need_rs) |
        (!is.na(key_cp) & key_cp %chin% need_cp)
    ]

    rm(g)
    invisible(gc(verbose = FALSE))

    out[[i]] <- keep
  }

  rbindlist(out, use.names = TRUE, fill = TRUE)
}


# Attach the full GWAS record to every model SNP. Matching uses chr:pos when
# both sides carry coordinates and rsID for whatever is left over.
annotate_variants_with_gwas <- function(v) {

  res <- list()

  for (anc in sort(unique(v$GWAS_ancestry))) {

    vv <- copy(v[GWAS_ancestry == anc])

    message(
      "GWAS look-up for ", anc, ": ", nrow(vv), " model SNPs, ",
      uniqueN(vv$gene), " genes"
    )

    hits <- collect_gwas_hits(
      need_rs = unique(vv$key_rs[!is.na(vv$key_rs)]),
      need_cp = unique(vv$key_cp[!is.na(vv$key_cp)]),
      chrs = sort(unique(vv$chr[!is.na(vv$chr)])),
      ancestry = anc
    )

    if (!nrow(hits)) {
      vv[, matched_by := "none"]
      res[[anc]] <- vv
      next
    }

    gw_cols <- setdiff(names(hits), c("key_rs", "key_cp"))

    best_by <- function(h, k) {

      h <- h[!is.na(get(k))]

      if (!nrow(h)) {
        return(NULL)
      }

      setorderv(h, c(k, "GWAS_P_used"), na.last = TRUE)
      n <- h[, .(n_GWAS_rows_matched = .N), by = k]
      h <- unique(h, by = k)

      merge(h, n, by = k, sort = FALSE)
    }

    cp <- best_by(copy(hits), "key_cp")
    rs <- best_by(copy(hits), "key_rs")

    take <- c(gw_cols, "n_GWAS_rows_matched")

    got <- NULL
    miss <- vv

    if (!is.null(cp)) {

      a <- merge(
        vv,
        cp[, c("key_cp", take), with = FALSE],
        by = "key_cp",
        all.x = TRUE,
        sort = FALSE
      )

      got <- a[!is.na(GWAS_file)]

      if (nrow(got)) {
        got[, matched_by := "chrpos"]
      }

      miss <- a[is.na(GWAS_file)]

      if (nrow(miss)) {
        miss[, (take) := NULL]
      }
    }

    b <- NULL

    if (nrow(miss)) {

      if (!is.null(rs)) {
        b <- merge(
          miss,
          rs[, c("key_rs", take), with = FALSE],
          by = "key_rs",
          all.x = TRUE,
          sort = FALSE
        )
      } else {
        b <- copy(miss)
      }

      b[, matched_by := if ("GWAS_file" %in% names(b)) {
        fifelse(!is.na(GWAS_file), "rsid", "none")
      } else {
        "none"
      }]
    }

    res[[anc]] <- rbindlist(list(got, b), use.names = TRUE, fill = TRUE)
  }

  rbindlist(res, use.names = TRUE, fill = TRUE)
}


# ---- Cistrome peak-overlap table (TF occupancy) --------------------------------

OVERLAP_COL_CANDIDATES <- list(
  CellType   = c("CellType", "cell_type", "celltype", "Cells"),
  gene       = c("gene", "gene_id", "geneid", "ensembl_gene_id"),
  genename   = c("genename", "gene_name", "gene_symbol", "symbol", "GeneName"),
  SNP_row_id = c("SNP_row_id", "snp_row_id", "row_id", "SNProwid"),
  rs         = c("rs", "rsid", "rsID", "rs_id", "SNP_rs", "snp_rs"),
  chr        = c("snp_chr", "SNP_chr", "variant_chr", "chr_snp", "chr", "CHR",
                 "chrom", "chromosome"),
  pos        = c("snp_pos", "SNP_pos", "variant_pos", "snp_bp", "pos", "POS",
                 "bp", "BP", "position"),
  TF         = c("Transcription factor", "Transcription_factor", "TF", "TFs",
                 "factor")
)


harmonize_overlap <- function(x, ct_override = NA_character_, src = NA_character_) {

  nms <- names(x)
  sel <- vapply(
    OVERLAP_COL_CANDIDATES,
    function(cand) pick_col(nms, cand),
    character(1)
  )

  for (must in c("gene", "rs", "TF")) {
    if (is.na(sel[[must]])) {
      stop(
        "No '", must, "' column found in ", src,
        "\n  columns seen: ", paste(nms, collapse = ", ")
      )
    }
  }

  out <- data.table(
    CellType = if (!is.na(ct_override)) {
      ct_override
    } else if (!is.na(sel[["CellType"]])) {
      trimws(as.character(x[[sel[["CellType"]]]]))
    } else {
      NA_character_
    },
    gene = sub("\\..*$", "", trimws(as.character(x[[sel[["gene"]]]]))),
    genename = if (!is.na(sel[["genename"]])) {
      trimws(as.character(x[[sel[["genename"]]]]))
    } else {
      NA_character_
    },
    SNP_row_id = if (!is.na(sel[["SNP_row_id"]])) {
      as.character(x[[sel[["SNP_row_id"]]]])
    } else {
      NA_character_
    },
    rs = trimws(as.character(x[[sel[["rs"]]]])),
    chr = if (!is.na(sel[["chr"]])) norm_chr(x[[sel[["chr"]]]]) else NA_character_,
    pos = if (!is.na(sel[["pos"]])) clean_pos(x[[sel[["pos"]]]]) else NA_character_,
    TF = trimws(as.character(x[[sel[["TF"]]]])),
    overlap_file = src
  )

  out <- add_variant_keys(out)

  attr(out, "column_map") <- sel

  out
}


# ------------------------------------------------------------------------------
# 2. SIGNIFICANT GENES AND EVERY ONE OF THEIR MODEL-BUILD SNPs
# ------------------------------------------------------------------------------

gene_tbl <- build_gene_table(
  path = SUPP_XLSX,
  sheet = SUPP_SHEET,
  skip = SUPP_SKIP,
  cell_types = CELL_TYPES,
  model_col = SUPP_MODEL_COL,
  single_model_source = MODEL_SOURCE_SINGLE
)

message(
  "Significant gene x cell type pairs in ", SUPP_SHEET, ": ", nrow(gene_tbl),
  " (", uniqueN(gene_tbl$gene), " genes)"
)

print(
  gene_tbl[
    , .N, by = .(CellType, Sig_resource)
  ][order(match(CellType, CELL_TYPES), Sig_resource)]
)

model_snps <- read_model_snps(gene_tbl)

model_snps <- merge(
  model_snps,
  gene_tbl[, .(CellType, gene, model_source, gene_name, Sig_resource, GWAS_ancestry)],
  by = c("CellType", "gene", "model_source"),
  all.x = TRUE,
  sort = FALSE
)

message(
  "Model-build SNPs collected: ", nrow(model_snps),
  " across ", uniqueN(paste(model_snps$CellType, model_snps$gene)),
  " gene x cell type pairs"
)

no_snps <- gene_tbl[
  !paste(CellType, gene) %chin% unique(paste(model_snps$CellType, model_snps$gene))
]

if (nrow(no_snps)) {
  message("\nGene x cell type pairs with no SNP in the weight files:")
  print(no_snps[, .(CellType, gene, gene_name, model_source)])
}


# ------------------------------------------------------------------------------
# 3. FULL GWAS RECORD FOR EVERY MODEL SNP; SMALLEST-P VARIANT FLAGGED
# ------------------------------------------------------------------------------

variant_gwas <- annotate_variants_with_gwas(model_snps)

# Rank the model SNPs of each gene x cell type by GWAS P. SNPs without a GWAS
# record sort last, so rank 1 is always the smallest observed P value.
variant_gwas[, rs_tiebreak := fifelse(is.na(key_rs), "~", key_rs)]

setorderv(
  variant_gwas,
  c("CellType", "gene", "GWAS_P_used", "rs_tiebreak", "variant_key"),
  na.last = TRUE
)

variant_gwas[, GWAS_P_rank_in_gene := seq_len(.N), by = .(CellType, gene)]

variant_gwas[
  ,
  is_smallest_GWAS_P := (GWAS_P_rank_in_gene == 1L) & !is.na(GWAS_P_used)
]

variant_gwas[, rs_tiebreak := NULL]


# ------------------------------------------------------------------------------
# 4. TF OCCUPANCY OF THE MODEL SNPs AND OF THE LEAD (SMALLEST-P) VARIANT
# ------------------------------------------------------------------------------

overlap_raw <- fread(OVERLAP_FILE, data.table = TRUE, check.names = FALSE)

ov <- harmonize_overlap(overlap_raw, src = basename(OVERLAP_FILE))

message(
  "Overlap column mapping: ",
  paste(
    names(attr(ov, "column_map")), attr(ov, "column_map"),
    sep = " -> ", collapse = " | "
  )
)

ov <- ov[
  CellType %chin% CELL_TYPES &
    !is.na(TF) & nzchar(TF) &
    !is.na(variant_key)
]

# The supplied overlap file contains the 51 CRC TFs.
ALL_TFS <- sort(unique(trimws(ov$TF)))
ALL_TFS <- ALL_TFS[nzchar(ALL_TFS)]

message("TFs found in overlap file: ", length(ALL_TFS))

# The overlap table and the weight files do not necessarily identify a SNP the
# same way (rsID vs chr:pos), so map the overlap rows onto the model SNPs on
# either key before counting TFs.
key_map <- unique(
  rbind(
    variant_gwas[!is.na(key_cp), .(CellType, gene, ov_key = key_cp, variant_key)],
    variant_gwas[!is.na(key_rs), .(CellType, gene, ov_key = key_rs, variant_key)]
  )
)

ov_long <- unique(
  rbind(
    ov[!is.na(key_cp), .(CellType, gene, ov_key = key_cp, TF)],
    ov[!is.na(key_rs), .(CellType, gene, ov_key = key_rs, TF)]
  )
)

ov_mapped <- unique(
  merge(
    ov_long,
    key_map,
    by = c("CellType", "gene", "ov_key"),
    all = FALSE,
    allow.cartesian = TRUE
  )[, .(CellType, gene, variant_key, TF)]
)

message(
  "Overlap gene x SNP x TF rows mapped onto a model SNP: ",
  uniqueN(ov_mapped[, .(CellType, gene, variant_key)]), " SNPs / ",
  nrow(ov_mapped), " TF assignments"
)

if (!nrow(ov_mapped)) {
  warning(
    "No overlap row could be matched to a model SNP. Check that the overlap ",
    "file identifies SNPs by rsID or by chromosome and position.",
    call. = FALSE
  )
}

tf_per_variant <- ov_mapped[
  ,
  .(
    n_TFs_on_variant = uniqueN(TF),
    TFs_on_variant = paste(sort(unique(TF)), collapse = ";")
  ),
  by = .(CellType, gene, variant_key)
]

variant_gwas <- merge(
  variant_gwas,
  tf_per_variant,
  by = c("CellType", "gene", "variant_key"),
  all.x = TRUE,
  sort = FALSE
)

variant_gwas[is.na(n_TFs_on_variant), n_TFs_on_variant := 0L]
variant_gwas[is.na(TFs_on_variant), TFs_on_variant := ""]

variant_gwas[, Dataset := DATASET_LABEL]

drop_cols <- intersect(c("key_rs", "key_cp", "chr_from_file"), names(variant_gwas))

if (length(drop_cols)) {
  variant_gwas[, (drop_cols) := NULL]
}

front_cols <- intersect(
  c(
    "Dataset", "CellType", "gene", "gene_name", "Sig_resource", "model_source",
    "GWAS_ancestry", "variant_key", "rs", "chr", "pos", "varID",
    "model_weight", "model_ref_allele", "model_eff_allele", "weight_file",
    "matched_by", "GWAS_file", "GWAS_P_used", "GWAS_P_rank_in_gene",
    "is_smallest_GWAS_P", "n_GWAS_rows_matched",
    "n_TFs_on_variant", "TFs_on_variant"
  ),
  names(variant_gwas)
)

setcolorder(
  variant_gwas,
  c(front_cols, setdiff(names(variant_gwas), front_cols))
)

variant_gwas[, ct_order := match(CellType, CELL_TYPES)]
setorderv(variant_gwas, c("ct_order", "gene_name", "GWAS_P_rank_in_gene"))
variant_gwas[, ct_order := NULL]

# The lead variant of a gene x cell type is the model SNP flagged above.
lead_keys <- variant_gwas[
  is_smallest_GWAS_P == TRUE,
  .(
    CellType,
    gene,
    variant_key,
    gene_name,
    lead_rs = rs,
    lead_GWAS_P = GWAS_P_used
  )
]

lead_overlap <- merge(
  ov_mapped,
  lead_keys,
  by = c("CellType", "gene", "variant_key"),
  all = FALSE,
  allow.cartesian = TRUE
)

# One gene can have repeated peaks/experiments for the same TF.
# Count each CellType x gene x TF only once.
lead_gene_tf <- unique(
  lead_overlap[
    ,
    .(
      Cells = CellType,
      gene,
      gene_name,
      lead_rs,
      lead_GWAS_P,
      TFs = TF
    )
  ]
)

# Per-gene view of the lead variant and the TFs occupying it.
lead_variant_detail <- lead_gene_tf[
  ,
  .(
    lead_rs = first(lead_rs),
    lead_GWAS_P = first(lead_GWAS_P),
    n_TFs_on_lead_variant = uniqueN(TFs),
    TFs_on_lead_variant = paste(
      sort(unique(TFs)),
      collapse = ";"
    )
  ),
  by = .(
    Cells,
    gene,
    gene_name
  )
]

lead_variant_detail[, cell_order := match(Cells, CELL_TYPES)]
setorder(
  lead_variant_detail,
  cell_order,
  gene_name,
  gene
)
lead_variant_detail[, cell_order := NULL]

# Requested summary:
# [Cells, TFs, #Genes lead variants be occupied]
lead_tf_count <- lead_gene_tf[
  ,
  .(
    n_genes_lead_variant_occupied = uniqueN(gene)
  ),
  by = .(Cells, TFs)
]

lead_tf_count <- complete_tf_grid(
  lead_tf_count,
  ALL_TFS
)

lead_tf_count[
  is.na(n_genes_lead_variant_occupied),
  n_genes_lead_variant_occupied := 0L
]

lead_tf_count[
  ,
  n_genes_lead_variant_occupied :=
    as.integer(n_genes_lead_variant_occupied)
]


# ------------------------------------------------------------------------------
# 5. READ EACH CORRELATION MATRIX AND COUNT GENES WITH r > 0.5
# ------------------------------------------------------------------------------

corr_count_list <- vector(
  "list",
  length(CELL_TYPES)
)

names(corr_count_list) <- CELL_TYPES

for (ct in CELL_TYPES) {

  cor_file <- cor_file_for_ct(ct)

  message(
    "Reading correlation matrix for ",
    ct, ": ", cor_file
  )

  x <- fread(
    cor_file,
    data.table = TRUE,
    check.names = FALSE
  )

  if (!"TF" %in% names(x)) {
    stop(
      "Correlation matrix for ", ct,
      " does not contain a 'TF' column: ",
      cor_file
    )
  }

  x[, TF := trimws(as.character(TF))]

  gene_cols <- setdiff(
    names(x),
    "TF"
  )

  if (!length(gene_cols)) {
    stop(
      "No TWAS risk-gene columns found in ",
      cor_file
    )
  }

  # Convert TF x gene matrix to long format.
  cor_long <- melt(
    x,
    id.vars = "TF",
    measure.vars = gene_cols,
    variable.name = "gene",
    value.name = "r",
    variable.factor = FALSE
  )

  cor_long[, r := suppressWarnings(as.numeric(r))]

  if (EXCLUDE_SELF_CORRELATION) {
    cor_long <- cor_long[
      toupper(trimws(TF)) !=
        toupper(trimws(gene))
    ]
  }

  # Count each risk gene once per TF.
  # In a matrix every TF-gene pair should occur once, but uniqueN()
  # keeps the logic robust.
  tmp <- cor_long[
    is.finite(r) &
      r > COR_THRESHOLD,
    .(
      n_genes_corr_gt_0.5 = uniqueN(gene)
    ),
    by = TF
  ]

  setnames(tmp, "TF", "TFs")
  tmp[, Cells := ct]

  corr_count_list[[ct]] <- tmp[
    ,
    .(
      Cells,
      TFs,
      n_genes_corr_gt_0.5
    )
  ]
}

corr_tf_count <- rbindlist(
  corr_count_list,
  use.names = TRUE,
  fill = TRUE
)

# In case a TF exists in a correlation matrix but not the overlap file,
# retain it in the common 51-TF grid.
TF_FROM_COR <- sort(
  unique(
    corr_tf_count$TFs
  )
)

ALL_TFS <- sort(
  union(
    ALL_TFS,
    TF_FROM_COR
  )
)

lead_tf_count <- complete_tf_grid(
  lead_tf_count[
    ,
    .(
      Cells,
      TFs,
      n_genes_lead_variant_occupied
    )
  ],
  ALL_TFS
)

lead_tf_count[
  is.na(n_genes_lead_variant_occupied),
  n_genes_lead_variant_occupied := 0L
]

corr_tf_count <- complete_tf_grid(
  corr_tf_count,
  ALL_TFS
)

corr_tf_count[
  is.na(n_genes_corr_gt_0.5),
  n_genes_corr_gt_0.5 := 0L
]

corr_tf_count[
  ,
  n_genes_corr_gt_0.5 :=
    as.integer(n_genes_corr_gt_0.5)
]


# ------------------------------------------------------------------------------
# 6. MERGE THE TWO SUMMARIES
# ------------------------------------------------------------------------------

summary_wide <- merge(
  lead_tf_count,
  corr_tf_count,
  by = c("Cells", "TFs"),
  all = TRUE,
  sort = FALSE
)

summary_wide[
  is.na(n_genes_lead_variant_occupied),
  n_genes_lead_variant_occupied := 0L
]

summary_wide[
  is.na(n_genes_corr_gt_0.5),
  n_genes_corr_gt_0.5 := 0L
]

summary_long <- melt(
  summary_wide,
  id.vars = c("Cells", "TFs"),
  measure.vars = c(
    "n_genes_lead_variant_occupied",
    "n_genes_corr_gt_0.5"
  ),
  variable.name = "Metric",
  value.name = "n_genes",
  variable.factor = FALSE
)

summary_long[
  Metric == "n_genes_lead_variant_occupied",
  Metric := "TF-occupied lead variant"
]

summary_long[
  Metric == "n_genes_corr_gt_0.5",
  Metric := paste0("Pearson r > ", COR_THRESHOLD)
]

summary_long[, Cells := factor(
  Cells,
  levels = CELL_TYPES
)]

summary_long[, TFs := factor(
  TFs,
  levels = ALL_TFS
)]


# ------------------------------------------------------------------------------
# 7. WRITE TABLES
# ------------------------------------------------------------------------------

fwrite(
  lead_variant_detail,
  file.path(
    OUTDIR,
    paste0(DATASET_LABEL, "_lead_variant_per_celltype_gene_with_occupying_TFs.csv")
  )
)

fwrite(
  variant_gwas,
  file.path(
    OUTDIR,
    paste0(DATASET_LABEL, "_modelSNPs_allGWASinfo_smallestP_flagged.csv")
  )
)

fwrite(
  lead_tf_count,
  file.path(
    OUTDIR,  paste0(DATASET_LABEL, "_", "TF_counts_genes_with_occupied_lead_variant.csv")
  )
)

fwrite(
  corr_tf_count,
  file.path(
    OUTDIR,
    paste0(
      "TF_counts_genes_with_correlation_gt_", DATASET_LABEL, "_", 
      COR_THRESHOLD,
      ".csv"
    )
  )
)

fwrite(
  summary_wide,
  file.path(
    OUTDIR,
    paste0(DATASET_LABEL, "_", "TF_leadVariantOccupation_vs_geneCorrelation_wide.csv")
  )
)

fwrite(
  summary_long,
  file.path(
    OUTDIR,
    paste0(DATASET_LABEL, "_", "TF_leadVariantOccupation_vs_geneCorrelation_plot_long.csv")
  )
)


# ------------------------------------------------------------------------------
# 8. ONE BAR PLOT PER CELL TYPE
# ------------------------------------------------------------------------------

plot_list <- vector("list", length(CELL_TYPES))
names(plot_list) <- CELL_TYPES

for (ct in CELL_TYPES) {

  # Order TFs within this cell type by the number of TWAS risk genes whose
  # lead variants are occupied by that TF. Ties are broken by the number of
  # correlated genes, then alphabetically by TF name.
  order_dt <- summary_wide[
    Cells == ct,
    .(
      TFs,
      n_genes_lead_variant_occupied,
      n_genes_corr_gt_0.5
    )
  ]

  setorder(
    order_dt,
    -n_genes_lead_variant_occupied,
    -n_genes_corr_gt_0.5,
    TFs
  )

  tf_order_ct <- order_dt$TFs

  plot_dt <- copy(
    summary_long[
      as.character(Cells) == ct
    ]
  )

  plot_dt[, TFs := factor(TFs, levels = tf_order_ct)]

  p <- ggplot(
    plot_dt,
    aes(
      x = TFs,
      y = n_genes,
      fill = Metric
    )
  ) +
    geom_col(
      position = position_dodge(width = 0.82),
      width = 0.72
    ) +
    scale_fill_manual(
      values = setNames(
        c("#4C78A8", "#E45756"),
        c(
          "TF-occupied lead variant",
          paste0("Pearson r > ", COR_THRESHOLD)
        )
      )
    ) +
    labs(
      title = paste0(
        ct,
        ": TF occupancy and co-expression across TWAS risk genes"
      ),
      x = "Transcription factor",
      y = "Number of TWAS risk genes",
      fill = NULL
    ) +
    theme_classic(base_size = 11) +
    theme(
      axis.text.x = element_text(
        angle = 60,
        hjust = 1,
        vjust = 1,
        size = 7
      ),
      legend.position = "top",
      plot.title = element_text(
        face = "bold"
      )
    )

  plot_list[[ct]] <- p

  ggsave(
    filename = file.path(
      OUTDIR,
      paste0(
        DATASET_LABEL, "_", ct,
        "_TF_leadVariantOccupation_vs_corrGT",
        COR_THRESHOLD,
        "_orderedByLeadVariant.pdf"
      )
    ),
    plot = p,
    width = 15,
    height = 6,
    units = "in",
    useDingbats = FALSE
  )

  ggsave(
    filename = file.path(
      OUTDIR,
      paste0(
        DATASET_LABEL, "_", ct,
        "_TF_leadVariantOccupation_vs_corrGT",
        COR_THRESHOLD,
        "_orderedByLeadVariant.png"
      )
    ),
    plot = p,
    width = 15,
    height = 6,
    units = "in",
    dpi = 300
  )
}


# ------------------------------------------------------------------------------
# 9. LARGE COMBINED PNG/PDF WITH ALL 7 CELL TYPES
# ------------------------------------------------------------------------------

save_plot_grid <- function(plot_list, filename, ncol = 3, width = 24, height = 24, res = 300) {

  n <- length(plot_list)
  nrow <- ceiling(n / ncol)

  if (grepl("\\.png$", filename, ignore.case = TRUE)) {
    grDevices::png(
      filename = filename,
      width = width,
      height = height,
      units = "in",
      res = res
    )
  } else if (grepl("\\.pdf$", filename, ignore.case = TRUE)) {
    grDevices::pdf(
      file = filename,
      width = width,
      height = height,
      useDingbats = FALSE
    )
  } else {
    stop("Unsupported output format for: ", filename)
  }

  grid::grid.newpage()
  pushViewport(viewport(layout = grid.layout(nrow, ncol)))

  for (i in seq_along(plot_list)) {
    row_i <- ceiling(i / ncol)
    col_i <- i - (row_i - 1L) * ncol

    print(
      plot_list[[i]],
      vp = viewport(layout.pos.row = row_i, layout.pos.col = col_i)
    )
  }

  invisible(dev.off())
}

save_plot_grid(
  plot_list = plot_list,
  filename = file.path(
    OUTDIR,
    paste0(
      DATASET_LABEL, "_", "ALL7CT_TF_leadVariantOccupation_vs_corrGT",
      COR_THRESHOLD,
      "_orderedByLeadVariant_large.png"
    )
  ),
  ncol = 3,
  width = 24,
  height = 24,
  res = 300
)

save_plot_grid(
  plot_list = plot_list,
  filename = file.path(
    OUTDIR,
    paste0(
      DATASET_LABEL, "_",
      "ALL7CT_TF_leadVariantOccupation_vs_corrGT",
      COR_THRESHOLD,
      "_orderedByLeadVariant_large.pdf"
    )
  ),
  ncol = 3,
  width = 24,
  height = 24
)


# ------------------------------------------------------------------------------
# 10. CONSOLE QC
# ------------------------------------------------------------------------------

message("\n============================================================")
message("DONE")
message("============================================================")
message("Output directory: ", OUTDIR)
message("Cell types: ", paste(CELL_TYPES, collapse = ", "))
message("TFs in final grid: ", length(ALL_TFS))
message("Correlation threshold: r > ", COR_THRESHOLD)
message(
  "Self-correlation excluded: ",
  EXCLUDE_SELF_CORRELATION
)

message("\nModel SNPs, GWAS annotation and lead variants per cell type:")
print(
  variant_gwas[
    ,
    .(
      n_genes = uniqueN(gene),
      n_model_SNPs = .N,
      n_SNPs_with_GWAS_record = sum(matched_by != "none"),
      n_lead_variants = sum(is_smallest_GWAS_P)
    ),
    by = .(Cells = CellType)
  ][order(match(Cells, CELL_TYPES))]
)

message("\nGWAS source and prediction model used (gene x cell type pairs):")
print(
  unique(variant_gwas[, .(CellType, gene, Sig_resource, GWAS_ancestry, model_source)])[
    , .N, by = .(Sig_resource, GWAS_ancestry, model_source)
  ]
)

message("\nHow the model SNPs were matched to the GWAS files:")
print(variant_gwas[, .N, by = matched_by][order(-N)])

message("\nRequested lead-variant TF count table:")
print(
  lead_tf_count[
    order(
      match(Cells, CELL_TYPES),
      TFs
    )
  ]
)

message("\nRequested correlation TF count table:")
print(
  corr_tf_count[
    order(
      match(Cells, CELL_TYPES),
      TFs
    )
  ]
)
