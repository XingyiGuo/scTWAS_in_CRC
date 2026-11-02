#!/usr/bin/env Rscript

# ==============================================================================
# WellcomeSanger healthy/NL:
# cell-type-specific CRC TF <-> NL TWAS-risk-gene pseudobulk correlations
#
# Biological question
# -------------------
# Within each healthy epithelial cell type, are the 51 CRC-risk TFs
# co-expressed across individuals with cell-type-specific NL TWAS-significant
# CRC-risk genes?
#
# Inputs
# ------
# 1. One WellcomeSanger pseudobulk RDS per cell type:
#    pseudobulk_expression_{CT}.rds
#
# 2. Cell-type-specific NL EUR/Multi significant TWAS genes:
#    HTAN_NL_WellcomeSanger_NL_{CT}_peakSNP.EUR_Multi.CauchyMeta.sig.csv
#
# 3. The same 51 CRC TF definition used in the CellOracleR analysis.
#
# For each of 7 healthy epithelial cell types:
#   1) read individual-level pseudobulk expression;
#   2) orient matrix as genes x individuals;
#   3) if the matrix is raw count-like:
#        TMM -> log2 CPM -> quantile normalization;
#      otherwise:
#        retain the existing transformed values and quantile normalize them;
#   4) calculate Pearson correlation:
#        51 TFs (rows) x CT-specific NL TWAS risk genes (columns);
#   5) calculate two-sided P values and BH FDR;
#   6) save correlation/P/FDR matrices and long-format results;
#   7) draw one correlation heatmap per cell type.
#
# IMPORTANT
# ---------
# * The pseudobulk RDS files are already individual-level pseudobulks, so this
#   script does NOT re-aggregate single cells.
# * Correlation is across independent individuals, not cells.
# * Correlation is complementary to CellOracleR; it does not itself establish
#   TF -> target regulation.
# * Self-correlations (TF == TWAS risk gene) are set to NA because r = 1 is
#   trivial and not informative.
#
# Run:
#   Rscript WellcomeSanger_NL_TF_TWASgene_pseudobulk_correlation.R
# ==============================================================================

RENV_ACTIVATE <- "/data/l2_bioinfo1/liq17/scTWAS/src/renv/activate.R"
if (file.exists(RENV_ACTIVATE)) {
  source(RENV_ACTIVATE)
}

suppressPackageStartupMessages({
  library(Matrix)
  library(data.table)
  library(dplyr)
  library(edgeR)
  library(preprocessCore)
  library(pheatmap)
})

options(stringsAsFactors = FALSE)
set.seed(20260903)

# ------------------------------------------------------------------------------
# 0. CONFIGURATION
# ------------------------------------------------------------------------------

CELL_TYPES <- c(
  "ABS", "CT", "EE", "GOB", "STM", "TAC", "TUF"
)

PB_RDS_TEMPLATE <- paste0(
  "/data/l2_bioinfo1/scRNA-seq/WellcomeSanger/processed/E-MTAB-16986/",
  "Healthy_TI_then_rectum_278samples_outputs_LauLabel/",
  "pseudobulk_Epithelial_LauLabel_x_individual/",
  "pseudobulk_expression_{CT}.rds"
)

TWAS_SIG_TEMPLATE <- paste0(
  "/data/l2_bioinfo1/liq17/scTWAS/",
  "TWAS_HTAN_GSE266546_WellcomeSanger_res/",
  "SPrediXcan_peakSNP_association/",
  "Cauchy_Meta_NL_WellcomeSanger_NL_EUR_Multi/",
  "HTAN_NL_WellcomeSanger_NL_{CT}_peakSNP.",
  "EUR_Multi.CauchyMeta.sig.csv"
)

TF_CSV <- Sys.getenv(
  "CRC_TF_CSV",
  unset = paste0(
    "/data/l2_bioinfo1/liq17/scTWAS/support_files/",
    "CRC_51TFs_TWAS_zhishan2024_with_ENSG.csv"
  )
)

OUTDIR <- paste0(
  "/data/l2_bioinfo1/liq17/scTWAS/CellOracleR_res/",
  "WellcomeSanger_NL_TF_TWASgene_pseudobulk_correlation"
)
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

COR_METHOD <- "pearson"
MIN_PB_SAMPLES <- 5L
EXCLUDE_SELF_CORRELATION <- TRUE

# Input-mode behavior:
#   "auto"   -> detect raw count-like vs already transformed expression.
#   "counts" -> force TMM -> log2CPM -> QN.
#   "expr"   -> treat RDS as already transformed expression -> QN only.
PB_INPUT_MODE <- "auto"

# Heatmap
# Plot only Pearson-r values: no EUR/Multi/Both annotation.
# TFs and risk genes are shown alphabetically with compact 6 x 6 point cells.
HEATMAP_LIMIT <- 1
HEATMAP_BREAKS <- seq(
  -HEATMAP_LIMIT,
  HEATMAP_LIMIT,
  length.out = 101
)

HEATMAP_COLORS <- grDevices::colorRampPalette(
  c("#2166AC", "white", "#B2182B")
)(100)

message("Output directory : ", OUTDIR)
message("Cell types       : ", paste(CELL_TYPES, collapse = ", "))
message("Correlation      : ", COR_METHOD)
message("Pseudobulk mode  : ", PB_INPUT_MODE)


# ------------------------------------------------------------------------------
# 1. HELPERS
# ------------------------------------------------------------------------------

replace_ct <- function(x, ct) {
  x <- gsub("\\{CT\\}", ct, x)
  gsub("\\{ct\\}", ct, x)
}

canonical_gene_id <- function(x) {
  x <- trimws(as.character(x))
  is_ensg <- grepl("^ENSG[0-9]+(\\.[0-9]+)?$", x)
  x[is_ensg] <- sub("\\.[0-9]+$", "", x[is_ensg])
  x
}

flag_yes <- function(x) {
  toupper(trimws(as.character(x))) %in%
    c("YES", "TRUE", "1", "Y")
}

first_nonempty_column <- function(df, candidates) {

  candidates <- candidates[
    candidates %in% colnames(df)
  ]

  out <- rep(
    NA_character_,
    nrow(df)
  )

  if (!length(candidates)) {
    return(out)
  }

  for (nm in candidates) {

    v <- trimws(
      as.character(df[[nm]])
    )

    use <- (
      is.na(out) |
        !nzchar(out)
    ) &
      !is.na(v) &
      nzchar(v)

    out[use] <- v[use]
  }

  out
}

# Extract a numeric expression matrix from common RDS object structures.
extract_numeric_matrix <- function(
    x,
    object_name = "RDS") {

  if (
    inherits(x, "Matrix") ||
      is.matrix(x)
  ) {
    return(x)
  }

  if (is.data.frame(x)) {

    y <- x

    # If one non-numeric column contains unique gene IDs, use it as row names.
    nonnum <- which(
      !vapply(
        y,
        is.numeric,
        logical(1)
      )
    )

    default_rn <- is.null(rownames(y)) ||
      all(
        rownames(y) ==
          as.character(seq_len(nrow(y)))
      )

    if (
      length(nonnum) &&
        default_rn
    ) {

      candidate <- y[[nonnum[1]]]

      if (
        length(unique(candidate)) ==
          nrow(y)
      ) {
        rownames(y) <- as.character(
          candidate
        )
        y[[nonnum[1]]] <- NULL
      }
    }

    keep_num <- vapply(
      y,
      is.numeric,
      logical(1)
    )

    if (!any(keep_num)) {
      stop(
        "No numeric columns found in ",
        object_name
      )
    }

    return(
      as.matrix(
        y[, keep_num, drop = FALSE]
      )
    )
  }

  if (
    inherits(
      x,
      "SummarizedExperiment"
    )
  ) {
    return(
      SummarizedExperiment::assay(
        x,
        1
      )
    )
  }

  if (is.list(x)) {

    preferred <- c(
      "expression",
      "expr",
      "pseudobulk",
      "matrix",
      "mat",
      "counts",
      "data",
      "assay"
    )

    for (nm in preferred) {

      if (!is.null(x[[nm]])) {

        out <- try(
          extract_numeric_matrix(
            x[[nm]],
            paste0(
              object_name,
              "$",
              nm
            )
          ),
          silent = TRUE
        )

        if (
          !inherits(
            out,
            "try-error"
          )
        ) {
          return(out)
        }
      }
    }

    # Last-resort recursive search.
    for (i in seq_along(x)) {

      out <- try(
        extract_numeric_matrix(
          x[[i]],
          paste0(
            object_name,
            "[[",
            i,
            "]]"
          )
        ),
        silent = TRUE
      )

      if (
        !inherits(
          out,
          "try-error"
        )
      ) {
        return(out)
      }
    }
  }

  stop(
    "Could not identify a numeric expression matrix in ",
    object_name
  )
}

orient_gene_by_sample <- function(
    m,
    reference_genes) {

  m <- as.matrix(m)

  if (
    is.null(rownames(m)) &&
      is.null(colnames(m))
  ) {
    stop(
      "Pseudobulk matrix has neither rownames nor colnames."
    )
  }

  reference_genes <- canonical_gene_id(
    reference_genes
  )

  rn <- if (
    is.null(rownames(m))
  ) {
    character()
  } else {
    canonical_gene_id(
      rownames(m)
    )
  }

  cn <- if (
    is.null(colnames(m))
  ) {
    character()
  } else {
    canonical_gene_id(
      colnames(m)
    )
  }

  nr <- sum(
    reference_genes %in% rn
  )

  nc <- sum(
    reference_genes %in% cn
  )

  if (nc > nr) {
    message(
      "  Transposing pseudobulk matrix to genes x individuals."
    )
    m <- t(m)
  }

  m
}

is_count_like <- function(m) {

  x <- as.numeric(m)
  x <- x[is.finite(x)]

  if (!length(x)) {
    return(FALSE)
  }

  if (min(x) < 0) {
    return(FALSE)
  }

  frac_integer <- mean(
    abs(x - round(x)) < 1e-8
  )

  # Raw pseudobulk counts should be overwhelmingly integer-valued.
  frac_integer >= 0.98
}

quantile_normalize <- function(m) {

  m <- as.matrix(m)
  storage.mode(m) <- "numeric"

  qn <- preprocessCore::normalize.quantiles(
    m,
    copy = TRUE
  )

  rownames(qn) <- rownames(m)
  colnames(qn) <- colnames(m)

  qn
}

normalize_pseudobulk <- function(
    m,
    input_mode = "auto") {

  m <- as.matrix(m)
  storage.mode(m) <- "numeric"

  if (
    any(!is.finite(m))
  ) {
    stop(
      "Pseudobulk matrix contains non-finite values before normalization."
    )
  }

  if (
    any(m < 0) &&
      input_mode == "counts"
  ) {
    stop(
      "PB_INPUT_MODE='counts' but negative values were detected."
    )
  }

  detected_count_like <- is_count_like(m)

  mode_used <- input_mode

  if (
    input_mode == "auto"
  ) {

    mode_used <- if (
      detected_count_like
    ) {
      "counts"
    } else {
      "expr"
    }
  }

  message(
    "  Input detected/used as: ",
    mode_used,
    if (
      input_mode == "auto"
    ) {
      paste0(
        " (count-like=",
        detected_count_like,
        ")"
      )
    } else {
      ""
    }
  )

  if (
    mode_used == "counts"
  ) {

    if (
      any(m < 0)
    ) {
      stop(
        "Negative values found in count-mode pseudobulk matrix."
      )
    }

    dge <- edgeR::DGEList(
      counts = m
    )

    dge <- edgeR::calcNormFactors(
      dge,
      method = "TMM"
    )

    log_cpm <- edgeR::cpm(
      dge,
      log = TRUE,
      prior.count = 1
    )

    qn <- quantile_normalize(
      log_cpm
    )

    return(
      list(
        input_mode = mode_used,
        raw_or_input = m,
        pre_QN = log_cpm,
        QN = qn,
        norm_factors = dge$samples
      )
    )
  }

  # Already transformed/normalized expression:
  # do not apply TMM to non-count data.
  warning(
    "Pseudobulk RDS does not look like raw integer counts. ",
    "TMM/logCPM will NOT be applied. ",
    "The existing expression values will be quantile normalized directly."
  )

  qn <- quantile_normalize(
    m
  )

  list(
    input_mode = mode_used,
    raw_or_input = m,
    pre_QN = m,
    QN = qn,
    norm_factors = NULL
  )
}

read_tf_table <- function(path) {

  if (!file.exists(path)) {
    stop(
      "51-TF file not found: ",
      path
    )
  }

  x <- data.table::fread(
    path,
    data.table = FALSE,
    check.names = FALSE
  )

  req <- c(
    "Transcription factor",
    "ENSG",
    "GLMM P values"
  )

  miss <- setdiff(
    req,
    colnames(x)
  )

  if (length(miss)) {
    stop(
      "TF file is missing required columns: ",
      paste(
        miss,
        collapse = ", "
      )
    )
  }

  x$`GLMM P values` <- as.numeric(
    x$`GLMM P values`
  )

  x$TF_ENSG <- canonical_gene_id(
    x$ENSG
  )

  x$TF_symbol <- trimws(
    as.character(
      x$`Transcription factor`
    )
  )

  x <- x[
    !is.na(x$`GLMM P values`) &
      x$`GLMM P values` < 0.001 &
      !is.na(x$TF_ENSG) &
      nzchar(x$TF_ENSG) &
      !is.na(x$TF_symbol) &
      nzchar(x$TF_symbol),
    ,
    drop = FALSE
  ]

  x <- x[
    !duplicated(x$TF_ENSG),
    ,
    drop = FALSE
  ]

  if (nrow(x) != 51L) {
    warning(
      "Expected 51 CRC TFs after GLMM P < 0.001; found ",
      nrow(x),
      "."
    )
  }

  x
}

read_ct_risk_genes <- function(
    path,
    ct) {

  if (!file.exists(path)) {
    stop(
      "NL TWAS significant-gene file not found for ",
      ct,
      ": ",
      path
    )
  }

  x <- data.table::fread(
    path,
    data.table = FALSE,
    check.names = FALSE
  )

  req <- c(
    "gene",
    "EUR_sig",
    "Multi_sig"
  )

  miss <- setdiff(
    req,
    colnames(x)
  )

  if (length(miss)) {
    stop(
      "TWAS file for ",
      ct,
      " is missing required columns: ",
      paste(
        miss,
        collapse = ", "
      )
    )
  }

  if (
    "CellType" %in%
      colnames(x)
  ) {

    input_ct <- unique(
      trimws(
        as.character(
          x$CellType
        )
      )
    )

    input_ct <- input_ct[
      !is.na(input_ct) &
        nzchar(input_ct)
    ]

    if (
      length(input_ct) &&
        any(input_ct != ct)
    ) {
      warning(
        "CellType values in ",
        basename(path),
        " do not all match ",
        ct,
        "; filename-derived CT is used."
      )
    }
  }

  x$gene <- canonical_gene_id(
    x$gene
  )

  gene_symbol <- first_nonempty_column(
    x,
    c(
      "gene_name_WellcomeSanger_NL_EUR",
      "gene_name_WellcomeSanger_NL_Multi",
      "gene_name_HTAN_NL_EUR",
      "gene_name_HTAN_NL_Multi",
      "gene_name_EUR",
      "gene_name_Multi",
      "gene_name"
    )
  )

  gene_symbol[
    is.na(gene_symbol) |
      !nzchar(gene_symbol)
  ] <- x$gene[
    is.na(gene_symbol) |
      !nzchar(gene_symbol)
  ]

  eur_yes <- flag_yes(
    x$EUR_sig
  )

  multi_yes <- flag_yes(
    x$Multi_sig
  )

  support <- ifelse(
    eur_yes &
      multi_yes,
    "Both_EUR_Multi",
    ifelse(
      eur_yes,
      "EUR_only",
      ifelse(
        multi_yes,
        "Multi_only",
        "Not_significant"
      )
    )
  )

  keep <- !is.na(x$gene) &
    nzchar(x$gene) &
    grepl(
      "^ENSG[0-9]+$",
      x$gene
    ) &
    (
      eur_yes |
        multi_yes
    )

  out <- data.frame(
    risk_ENSG = x$gene[keep],
    risk_gene_symbol = gene_symbol[keep],
    EUR_sig = as.character(
      x$EUR_sig[keep]
    ),
    Multi_sig = as.character(
      x$Multi_sig[keep]
    ),
    GWAS_support = support[keep],
    stringsAsFactors = FALSE
  )

  support_rank <- match(
    out$GWAS_support,
    c(
      "Both_EUR_Multi",
      "EUR_only",
      "Multi_only"
    )
  )

  out$.support_rank <- support_rank

  out <- out |>
    dplyr::arrange(
      .support_rank,
      risk_ENSG
    ) |>
    dplyr::distinct(
      risk_ENSG,
      .keep_all = TRUE
    ) |>
    dplyr::select(
      -.support_rank
    )

  out
}

# Map ENSG and symbol identifiers onto exact row names in the pseudobulk matrix.
make_feature_mapper <- function(
    feature_names) {

  feature_names <- as.character(
    feature_names
  )

  feature_keys <- canonical_gene_id(
    feature_names
  )

  ensg_lookup <- setNames(
    feature_names,
    feature_keys
  )

  ensg_lookup <- ensg_lookup[
    !duplicated(
      names(ensg_lookup)
    )
  ]

  symbol_lookup <- setNames(
    feature_names,
    toupper(feature_names)
  )

  symbol_lookup <- symbol_lookup[
    !duplicated(
      names(symbol_lookup)
    )
  ]

  function(
      ensg,
      symbol = NA_character_) {

    out <- unname(
      ensg_lookup[
        canonical_gene_id(ensg)
      ]
    )

    miss <- is.na(out) |
      !nzchar(out)

    if (any(miss)) {

      sym <- toupper(
        trimws(
          as.character(
            symbol[miss]
          )
        )
      )

      out[miss] <- unname(
        symbol_lookup[sym]
      )
    }

    out
  }
}

calculate_correlations <- function(
    expr,
    tf_meta,
    tf_features,
    risk_meta,
    risk_features,
    risk_labels,
    ct) {

  n <- ncol(expr)

  cor_mat <- matrix(
    NA_real_,
    nrow = nrow(tf_meta),
    ncol = nrow(risk_meta),
    dimnames = list(
      tf_meta$TF_symbol,
      risk_labels
    )
  )

  p_mat <- cor_mat

  if (
    n < MIN_PB_SAMPLES
  ) {
    warning(
      ct,
      ": only ",
      n,
      " individual pseudobulk samples; ",
      "correlations left as NA."
    )

    return(
      list(
        cor = cor_mat,
        p = p_mat,
        n = n
      )
    )
  }

  tf_present <- !is.na(tf_features) &
    tf_features %in%
      rownames(expr)

  risk_present <- !is.na(risk_features) &
    risk_features %in%
      rownames(expr)

  if (
    !any(tf_present) ||
      !any(risk_present)
  ) {
    return(
      list(
        cor = cor_mat,
        p = p_mat,
        n = n
      )
    )
  }

  tf_expr <- expr[
    tf_features[tf_present],
    ,
    drop = FALSE
  ]

  risk_expr <- expr[
    risk_features[risk_present],
    ,
    drop = FALSE
  ]

  r_sub <- stats::cor(
    t(tf_expr),
    t(risk_expr),
    method = COR_METHOD,
    use = "pairwise.complete.obs"
  )

  cor_mat[
    which(tf_present),
    which(risk_present)
  ] <- r_sub

  # Pearson correlation two-sided P value.
  denom <- pmax(
    1 - r_sub^2,
    .Machine$double.eps
  )

  t_stat <- r_sub *
    sqrt(
      (n - 2) /
        denom
    )

  p_sub <- 2 *
    stats::pt(
      -abs(t_stat),
      df = n - 2
    )

  p_sub[
    is.na(r_sub)
  ] <- NA_real_

  p_mat[
    which(tf_present),
    which(risk_present)
  ] <- p_sub

  # Remove trivial TF self-correlations.
  if (EXCLUDE_SELF_CORRELATION) {

    tf_ids <- canonical_gene_id(
      tf_meta$TF_ENSG
    )

    risk_ids <- canonical_gene_id(
      risk_meta$risk_ENSG
    )

    for (
      i in seq_along(tf_ids)
    ) {

      j <- which(
        risk_ids ==
          tf_ids[i]
      )

      if (length(j)) {
        cor_mat[i, j] <- NA_real_
        p_mat[i, j] <- NA_real_
      }
    }
  }

  list(
    cor = cor_mat,
    p = p_mat,
    n = n
  )
}

matrix_to_long <- function(
    cor_mat,
    p_mat,
    ct,
    tf_meta,
    risk_meta,
    n_samples,
    tf_features,
    risk_features) {

  idx <- expand.grid(
    TF_row = seq_len(
      nrow(cor_mat)
    ),
    gene_col = seq_len(
      ncol(cor_mat)
    ),
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )

  out <- data.frame(
    cell_type = ct,
    TF = rownames(cor_mat)[
      idx$TF_row
    ],
    TF_ENSG = tf_meta$TF_ENSG[
      idx$TF_row
    ],
    TF_feature_in_RDS = tf_features[
      idx$TF_row
    ],
    risk_gene = colnames(cor_mat)[
      idx$gene_col
    ],
    risk_gene_symbol = risk_meta$risk_gene_symbol[
      idx$gene_col
    ],
    risk_gene_ENSG = risk_meta$risk_ENSG[
      idx$gene_col
    ],
    risk_feature_in_RDS = risk_features[
      idx$gene_col
    ],
    EUR_sig = risk_meta$EUR_sig[
      idx$gene_col
    ],
    Multi_sig = risk_meta$Multi_sig[
      idx$gene_col
    ],
    GWAS_support = risk_meta$GWAS_support[
      idx$gene_col
    ],
    r = cor_mat[
      cbind(
        idx$TF_row,
        idx$gene_col
      )
    ],
    P = p_mat[
      cbind(
        idx$TF_row,
        idx$gene_col
      )
    ],
    n_individuals = n_samples,
    stringsAsFactors = FALSE
  )

  out$self_pair <- canonical_gene_id(
    out$TF_ENSG
  ) ==
    canonical_gene_id(
      out$risk_gene_ENSG
    )

  out$FDR <- stats::p.adjust(
    out$P,
    method = "BH"
  )

  out
}

# Draw one heatmap per cell type, preserving alphabetical TF and gene order.
draw_single_heatmap <- function(
    cor_mat,
    ct,
    source_label,
    n_samples,
    outfile) {

  stopifnot(
    identical(rownames(cor_mat), sort(rownames(cor_mat))),
    identical(colnames(cor_mat), sort(colnames(cor_mat)))
  )

  # --------------------------------------------------------------------------
  # Draw heatmap first without writing directly to PDF.
  # --------------------------------------------------------------------------
  ph <- pheatmap::pheatmap(
    cor_mat,
    color = HEATMAP_COLORS,
    breaks = HEATMAP_BREAKS,

    cluster_rows = TRUE,
    cluster_cols = TRUE,

    main = paste0(
      ct, ": 51 CRC TFs vs cell-type-specific TWAS risk genes\n",
      source_label,
      "; patient-level pseudobulk; n = ",
      n_samples
    ),

    border_color = NA,
    na_col = "grey90",

    # Heatmap cells
    cellwidth = 9,
    cellheight = 9,

    # Global fontsize mainly affects legend + main
    fontsize = 9,

    # Explicit row/column label sizes
    fontsize_row = 6,
    fontsize_col = 6,

    angle_col = 45,
    legend = TRUE,

    silent = TRUE
  )

  # --------------------------------------------------------------------------
  # Add outer margins around the complete pheatmap. Units are inches.
  # --------------------------------------------------------------------------
  margin_left   <- 0.7
  margin_right  <- 0.8
  margin_top    <- 0.6
  margin_bottom <- 0.7

  # Native pheatmap gtable dimensions
  heatmap_width <- grid::convertWidth(
    sum(ph$gtable$widths),
    unitTo = "in",
    valueOnly = TRUE
  )

  heatmap_height <- grid::convertHeight(
    sum(ph$gtable$heights),
    unitTo = "in",
    valueOnly = TRUE
  )

  pdf_width <- heatmap_width +
    margin_left +
    margin_right

  pdf_height <- heatmap_height +
    margin_top +
    margin_bottom

  # --------------------------------------------------------------------------
  # Save manually so margins can be controlled.
  # --------------------------------------------------------------------------
  grDevices::pdf(
    outfile,
    width = pdf_width,
    height = pdf_height,
    useDingbats = FALSE
  )

  grid::grid.newpage()

  # Position the complete heatmap inside the larger PDF canvas
  vp <- grid::viewport(
    x = grid::unit(
      margin_left + heatmap_width / 2,
      "in"
    ),
    y = grid::unit(
      margin_bottom + heatmap_height / 2,
      "in"
    ),
    width = grid::unit(
      heatmap_width,
      "in"
    ),
    height = grid::unit(
      heatmap_height,
      "in"
    )
  )

  grid::pushViewport(vp)
  grid::grid.draw(ph$gtable)
  grid::popViewport()

  grDevices::dev.off()

  # Preserve current downstream files that record the heatmap ordering.
  list(
    row_order = rownames(cor_mat),
    col_order = colnames(cor_mat)
  )
}

# ------------------------------------------------------------------------------
# 2. READ 51 TFs
# ------------------------------------------------------------------------------

tf_meta <- read_tf_table(
  TF_CSV
)

# Fix TF order alphabetically for all downstream matrices and heatmaps.
tf_meta <- tf_meta |>
  dplyr::arrange(TF_symbol)

message(
  "CRC TFs retained in alphabetical order: ",
  nrow(tf_meta)
)

data.table::fwrite(
  tf_meta,
  file.path(
    OUTDIR,
    "CRC51TFs_used.csv"
  )
)


# ------------------------------------------------------------------------------
# 3. MAIN LOOP
# ------------------------------------------------------------------------------

all_long <- list()
all_qc <- list()

for (
  CT in CELL_TYPES
) {

  message("\n")
  message(
    "======================================================================"
  )
  message(
    "CELL TYPE: ",
    CT
  )
  message(
    "======================================================================"
  )

  ct_outdir <- file.path(
    OUTDIR,
    CT
  )

  dir.create(
    ct_outdir,
    recursive = TRUE,
    showWarnings = FALSE
  )

  pb_file <- replace_ct(
    PB_RDS_TEMPLATE,
    CT
  )

  twas_file <- replace_ct(
    TWAS_SIG_TEMPLATE,
    CT
  )

  message(
    "Pseudobulk RDS: ",
    pb_file
  )

  message(
    "TWAS file     : ",
    twas_file
  )

  if (
    !file.exists(pb_file)
  ) {
    warning(
      "Missing pseudobulk RDS for ",
      CT,
      ": ",
      pb_file
    )
    next
  }

  risk_meta <- read_ct_risk_genes(
    twas_file,
    CT
  )

  if (
    !nrow(risk_meta)
  ) {
    warning(
      "No significant NL TWAS genes for ",
      CT,
      "; skipping."
    )
    next
  }

  message(
    "CT-specific NL TWAS risk genes: ",
    nrow(risk_meta)
  )

  risk_labels <- ifelse(
    !is.na(
      risk_meta$risk_gene_symbol
    ) &
      nzchar(
        risk_meta$risk_gene_symbol
      ),
    risk_meta$risk_gene_symbol,
    risk_meta$risk_ENSG
  )

  risk_meta$risk_gene_label <-
    risk_labels

  # Fix risk-gene order alphabetically; ENSG is a deterministic tie-breaker.
  risk_meta <- risk_meta |>
    dplyr::arrange(
      risk_gene_label,
      risk_ENSG
    )

  # Keep both genes if duplicate display symbols occur.
  risk_meta$risk_gene_label <- make.unique(
    as.character(risk_meta$risk_gene_label)
  )

  risk_labels <- risk_meta$risk_gene_label

  data.table::fwrite(
    risk_meta,
    file.path(
      ct_outdir,
      paste0(
        CT,
        "_NL_TWAS_risk_genes_used.csv"
      )
    )
  )

  pb_obj <- readRDS(
    pb_file
  )

  pb <- extract_numeric_matrix(
    pb_obj,
    basename(pb_file)
  )

  # Use TF + TWAS ENSGs to determine orientation.
  orientation_reference <- unique(
    c(
      tf_meta$TF_ENSG,
      risk_meta$risk_ENSG
    )
  )

  pb <- orient_gene_by_sample(
    pb,
    orientation_reference
  )

  if (
    is.null(rownames(pb)) ||
      is.null(colnames(pb))
  ) {
    stop(
      CT,
      ": pseudobulk matrix must have both gene row names and individual column names."
    )
  }

  # Canonicalize ENSG row names only if that does not create duplicates.
  pb_rn_canonical <- canonical_gene_id(
    rownames(pb)
  )

  if (
    !anyDuplicated(
      pb_rn_canonical
    )
  ) {
    rownames(pb) <-
      pb_rn_canonical
  }

  message(
    "Pseudobulk matrix: ",
    nrow(pb),
    " genes x ",
    ncol(pb),
    " individuals"
  )

  if (
    ncol(pb) < MIN_PB_SAMPLES
  ) {
    warning(
      CT,
      ": fewer than ",
      MIN_PB_SAMPLES,
      " individuals."
    )
  }

  # Basic input diagnostics before normalization.
  finite_pb <- pb[
    is.finite(pb)
  ]

  pb_input_qc <- data.frame(
    cell_type = CT,
    n_genes = nrow(pb),
    n_individuals = ncol(pb),
    min_value = min(
      finite_pb,
      na.rm = TRUE
    ),
    median_value = stats::median(
      finite_pb,
      na.rm = TRUE
    ),
    q99_value = unname(
      stats::quantile(
        finite_pb,
        0.99,
        na.rm = TRUE
      )
    ),
    max_value = max(
      finite_pb,
      na.rm = TRUE
    ),
    fraction_integer_like = mean(
      abs(
        finite_pb -
          round(finite_pb)
      ) < 1e-8
    ),
    count_like = is_count_like(pb),
    stringsAsFactors = FALSE
  )

  data.table::fwrite(
    pb_input_qc,
    file.path(
      ct_outdir,
      paste0(
        CT,
        "_pseudobulk_input_QC.csv"
      )
    )
  )

  norm <- normalize_pseudobulk(
    pb,
    input_mode = PB_INPUT_MODE
  )

  expr <- norm$QN

  # Save normalization-factor information when count mode was used.
  if (
    !is.null(
      norm$norm_factors
    )
  ) {

    norm_info <- data.frame(
      individual_id =
        rownames(
          norm$norm_factors
        ),
      norm$norm_factors,
      row.names = NULL,
      check.names = FALSE
    )

    data.table::fwrite(
      norm_info,
      file.path(
        ct_outdir,
        paste0(
          CT,
          "_TMM_normalization_factors.csv"
        )
      )
    )
  }

  mapper <- make_feature_mapper(
    rownames(expr)
  )

  tf_features <- mapper(
    tf_meta$TF_ENSG,
    tf_meta$TF_symbol
  )

  risk_features <- mapper(
    risk_meta$risk_ENSG,
    risk_meta$risk_gene_symbol
  )

  tf_presence <- data.frame(
    TF = tf_meta$TF_symbol,
    TF_ENSG = tf_meta$TF_ENSG,
    feature_in_RDS =
      tf_features,
    present =
      !is.na(tf_features),
    stringsAsFactors = FALSE
  )

  risk_presence <- data.frame(
    risk_gene = risk_labels,
    risk_gene_symbol =
      risk_meta$risk_gene_symbol,
    risk_gene_ENSG =
      risk_meta$risk_ENSG,
    GWAS_support =
      risk_meta$GWAS_support,
    feature_in_RDS =
      risk_features,
    present =
      !is.na(risk_features),
    stringsAsFactors = FALSE
  )

  data.table::fwrite(
    tf_presence,
    file.path(
      ct_outdir,
      paste0(
        CT,
        "_CRC51TF_feature_presence.csv"
      )
    )
  )

  data.table::fwrite(
    risk_presence,
    file.path(
      ct_outdir,
      paste0(
        CT,
        "_NL_TWASriskGene_feature_presence.csv"
      )
    )
  )

  selected_features <- unique(
    na.omit(
      c(
        tf_features,
        risk_features
      )
    )
  )

  # Save the expression used for correlation for TF + risk genes.
  selected_expr <- expr[
    selected_features,
    ,
    drop = FALSE
  ]

  data.table::fwrite(
    data.frame(
      gene_feature =
        rownames(
          selected_expr
        ),
      selected_expr,
      check.names = FALSE
    ),
    file.path(
      ct_outdir,
      paste0(
        CT,
        "_selectedGenes_expression_used_for_correlation.csv.gz"
      )
    )
  )

  cor_res <- calculate_correlations(
    expr = expr,
    tf_meta = tf_meta,
    tf_features = tf_features,
    risk_meta = risk_meta,
    risk_features = risk_features,
    risk_labels = risk_labels,
    ct = CT
  )

  # Correlation matrix.
  data.table::fwrite(
    data.frame(
      TF =
        rownames(
          cor_res$cor
        ),
      cor_res$cor,
      check.names = FALSE
    ),
    file.path(
      ct_outdir,
      paste0(
        CT,
        "_CRC51TF_x_NL_TWASriskGene_PearsonR_matrix.csv"
      )
    )
  )

  # P-value matrix.
  data.table::fwrite(
    data.frame(
      TF =
        rownames(
          cor_res$p
        ),
      cor_res$p,
      check.names = FALSE
    ),
    file.path(
      ct_outdir,
      paste0(
        CT,
        "_CRC51TF_x_NL_TWASriskGene_P_matrix.csv"
      )
    )
  )

  long_res <- matrix_to_long(
    cor_mat = cor_res$cor,
    p_mat = cor_res$p,
    ct = CT,
    tf_meta = tf_meta,
    risk_meta = risk_meta,
    n_samples = cor_res$n,
    tf_features = tf_features,
    risk_features = risk_features
  )

  data.table::fwrite(
    long_res,
    file.path(
      ct_outdir,
      paste0(
        CT,
        "_CRC51TF_x_NL_TWASriskGene_correlation_long.csv"
      )
    )
  )

  # FDR matrix.
  fdr_mat <- matrix(
    long_res$FDR,
    nrow = nrow(cor_res$cor),
    ncol = ncol(cor_res$cor),
    byrow = FALSE,
    dimnames = dimnames(
      cor_res$cor
    )
  )

  data.table::fwrite(
    data.frame(
      TF =
        rownames(
          fdr_mat
        ),
      fdr_mat,
      check.names = FALSE
    ),
    file.path(
      ct_outdir,
      paste0(
        CT,
        "_CRC51TF_x_NL_TWASriskGene_FDR_matrix.csv"
      )
    )
  )

  sig_res <- long_res |>
    dplyr::filter(
      !is.na(FDR),
      FDR < 0.05
    ) |>
    dplyr::arrange(
      FDR,
      dplyr::desc(
        abs(r)
      )
    )

  data.table::fwrite(
    sig_res,
    file.path(
      ct_outdir,
      paste0(
        CT,
        "_CRC51TF_x_NL_TWASriskGene_FDR005.csv"
      )
    )
  )

  heatmap_order <- draw_single_heatmap(
    cor_mat = cor_res$cor,
    ct = CT,
    source_label = "WellcomeSanger healthy",
    n_samples = cor_res$n,
    outfile = file.path(
      ct_outdir,
      paste0(
        CT,
        "_CRC51TF_x_NL_TWASriskGene_correlation_heatmap.pdf"
      )
    )
  )

  data.table::fwrite(
    data.frame(
      order =
        seq_along(
          heatmap_order$row_order
        ),
      TF =
        heatmap_order$row_order
    ),
    file.path(
      ct_outdir,
      paste0(
        CT,
        "_heatmap_TF_order.csv"
      )
    )
  )

  data.table::fwrite(
    data.frame(
      order =
        seq_along(
          heatmap_order$col_order
        ),
      risk_gene =
        heatmap_order$col_order
    ),
    file.path(
      ct_outdir,
      paste0(
        CT,
        "_heatmap_riskGene_order.csv"
      )
    )
  )

  all_long[[CT]] <- long_res

  all_qc[[CT]] <- data.frame(
    cell_type = CT,
    pseudobulk_input_mode =
      norm$input_mode,
    n_individuals =
      ncol(expr),
    n_CRC51_TFs =
      nrow(tf_meta),
    n_CRC51_TFs_present =
      sum(
        !is.na(
          tf_features
        )
      ),
    n_NL_TWAS_risk_genes =
      nrow(risk_meta),
    n_NL_TWAS_risk_genes_present =
      sum(
        !is.na(
          risk_features
        )
      ),
    n_tested_pairs =
      sum(
        !is.na(
          long_res$r
        )
      ),
    n_FDR005_pairs =
      sum(
        !is.na(
          long_res$FDR
        ) &
          long_res$FDR < 0.05
      ),
    stringsAsFactors = FALSE
  )

  rm(
    pb_obj,
    pb,
    expr,
    norm
  )

  gc(
    verbose = FALSE
  )
}


# ------------------------------------------------------------------------------
# 4. COMBINE RESULTS ACROSS 7 CELL TYPES
# ------------------------------------------------------------------------------

if (
  length(all_long)
) {

  all_long_df <- dplyr::bind_rows(
    all_long
  )

  data.table::fwrite(
    all_long_df,
    file.path(
      OUTDIR,
      "ALL7CT_WellcomeSanger_NL_CRC51TF_x_TWASriskGene_correlation_long.csv.gz"
    )
  )

  all_sig <- all_long_df |>
    dplyr::filter(
      !is.na(FDR),
      FDR < 0.05
    ) |>
    dplyr::arrange(
      cell_type,
      FDR,
      dplyr::desc(
        abs(r)
      )
    )

  data.table::fwrite(
    all_sig,
    file.path(
      OUTDIR,
      "ALL7CT_WellcomeSanger_NL_CRC51TF_x_TWASriskGene_FDR005.csv"
    )
  )

  # TF-level summary: which TFs have the most significant NL risk-gene
  # correlations in each cell type?
  tf_summary <- all_long_df |>
    dplyr::group_by(
      cell_type,
      TF,
      TF_ENSG
    ) |>
    dplyr::summarise(
      n_risk_genes_tested =
        sum(
          !is.na(r)
        ),
      n_FDR005_risk_genes =
        sum(
          !is.na(FDR) &
            FDR < 0.05
        ),
      n_positive_FDR005 =
        sum(
          !is.na(FDR) &
            FDR < 0.05 &
            r > 0
        ),
      n_negative_FDR005 =
        sum(
          !is.na(FDR) &
            FDR < 0.05 &
            r < 0
        ),
      max_abs_r =
        if (
          any(
            is.finite(r)
          )
        ) {
          max(
            abs(r),
            na.rm = TRUE
          )
        } else {
          NA_real_
        },
      .groups = "drop"
    ) |>
    dplyr::arrange(
      cell_type,
      dplyr::desc(
        n_FDR005_risk_genes
      ),
      dplyr::desc(
        max_abs_r
      )
    )

  data.table::fwrite(
    tf_summary,
    file.path(
      OUTDIR,
      "ALL7CT_WellcomeSanger_NL_TF_correlation_summary.csv"
    )
  )
}

if (
  length(all_qc)
) {

  qc_df <- dplyr::bind_rows(
    all_qc
  )

  data.table::fwrite(
    qc_df,
    file.path(
      OUTDIR,
      "ALL7CT_WellcomeSanger_NL_correlation_QC_summary.csv"
    )
  )

  print(
    qc_df
  )
}


# ------------------------------------------------------------------------------
# 5. SESSION INFO
# ------------------------------------------------------------------------------

capture.output(
  sessionInfo(),
  file = file.path(
    OUTDIR,
    "sessionInfo.txt"
  )
)

message("\n")
message(
  "======================================================================"
)
message("DONE")
message(
  "======================================================================"
)
message(
  "Results: ",
  OUTDIR
)
message("")
message(
  "Key output per cell type:"
)
message(
  "  *_CRC51TF_x_NL_TWASriskGene_PearsonR_matrix.csv"
)
message(
  "  *_CRC51TF_x_NL_TWASriskGene_P_matrix.csv"
)
message(
  "  *_CRC51TF_x_NL_TWASriskGene_FDR_matrix.csv"
)
message(
  "  *_CRC51TF_x_NL_TWASriskGene_correlation_long.csv"
)
message(
  "  *_CRC51TF_x_NL_TWASriskGene_FDR005.csv"
)
message(
  "  *_CRC51TF_x_NL_TWASriskGene_correlation_heatmap.pdf"
)
