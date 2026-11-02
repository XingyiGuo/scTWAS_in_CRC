#!/usr/bin/env Rscript

# ==============================================================================
# HTAN ADSER: one cell-type-specific TF <-> TWAS-risk-gene pseudobulk correlation
#
# Biological question
# -------------------
# Within each epithelial cell type, are the 51 CRC-risk TFs co-expressed with
# cell-type-specific ADSER TWAS-significant genes across independent patients?
#
# Analysis design
# ---------------
# * ABS, GOB, CT, EE, TUF, TAC, STM:
#     pool AD + SER cells FIRST, then aggregate to PatientID-level pseudobulk,
#     normalize expression, and calculate ONE TF-risk-gene correlation matrix.
# * ASC:
#     use AD cells only, then calculate ONE correlation matrix.
# * SSC:
#     use SER cells only, then calculate ONE correlation matrix.
# * Full correlation matrices retain all TFs/genes; undefined correlations are NA.
# * Heatmaps cluster TFs and TWAS risk genes by correlation pattern.
# * Features with no finite correlations are excluded from heatmap clustering only.
#
# Correlation is calculated across PATIENT-LEVEL pseudobulk samples, NOT cells.
# If the same PatientID contributes both AD and SER cells in a shared lineage,
# those cells are combined into that patient's single pooled ADSER pseudobulk.
#
# Correlation does not itself establish TF -> target regulation. It provides
# expression-level support complementary to CellOracleR TF -> target predictions.
#
# Run:
#   Rscript HTAN_ADSER_TF_TWASgene_pseudobulk_correlation_ADSER_merged.R
# ==============================================================================

RENV_ACTIVATE <- "/data/l2_bioinfo1/liq17/scTWAS/src/renv/activate.R"
if (file.exists(RENV_ACTIVATE)) {
  source(RENV_ACTIVATE)
}

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(data.table)
  library(dplyr)
  library(edgeR)
  library(preprocessCore)
  library(pheatmap)
})

options(stringsAsFactors = FALSE)
set.seed(20260910)

# ------------------------------------------------------------------------------
# 0. CONFIGURATION
# ------------------------------------------------------------------------------

CELL_TYPES <- c(
  "ABS", "GOB", "CT", "EE", "TUF", "TAC", "STM", "ASC", "SSC"
)

SHARED_CELL_TYPES <- c(
  "ABS", "GOB", "CT", "EE", "TUF", "TAC", "STM"
)

RDS_TEMPLATE <- paste0(
  "/data/l2_bioinfo1/liq17/scTWAS/data/VUMC_HTAN/VUMC_HTAN_DIS_VAL_ADSER/Results/by_celltype_DISVAL_PCA83_ADSER_HPtoSER_ASCSSCfiltered/",
  "HTAN_DISVAL_PCA83_ADSER_HPtoSER_{ct}_CellOracleR_input.rds"
)

TWAS_SIG_TEMPLATE <- paste0(
  "/v5000/data/l2_bioinfo1/liq17/scTWAS/",
  "TWAS_HTAN_GSE266546_WellcomeSanger_res/",
  "SPrediXcan_peakSNP_association/",
  "VUMC_HTAN_ADSER_EUR_Multi_sig/",
  "VUMC_HTAN_ADSER_{ct}_peakSNP.EUR_Multi.sig.csv"
)

TF_CSV <- paste0(
  "/data/l2_bioinfo1/liq17/scTWAS/support_files/",
  "CRC_51TFs_TWAS_zhishan2024_with_ENSG.csv"
)

OUTDIR <- paste0(
  "/data/l2_bioinfo1/liq17/scTWAS/CellOracleR_res/",
  "HTAN_ADSER_TF_TWASgene_pseudobulk_correlation_merged"
)
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

PSEUDOBULK_ID_COL <- "PatientID"
CONDITION_COL <- "sample_condition"

# Earlier preparation scripts sometimes used "SET" for the serrated group.
SER_ALIASES <- c("SER", "SET")

MIN_PB_SAMPLES <- 5L
COR_METHOD <- "pearson"

HEATMAP_LIMIT <- 1
HEATMAP_BREAKS <- seq(-HEATMAP_LIMIT, HEATMAP_LIMIT, length.out = 101)
HEATMAP_COLORS <- grDevices::colorRampPalette(
  c("#2166AC", "white", "#B2182B")
)(100)

message("Output directory: ", OUTDIR)
message("Cell types: ", paste(CELL_TYPES, collapse = ", "))
message("Pseudobulk unit: ", PSEUDOBULK_ID_COL)
message(
  "Design: shared 7 CTs = pooled AD+SER; ASC = AD only; SSC = SER only."
)

# ------------------------------------------------------------------------------
# 1. HELPERS
# ------------------------------------------------------------------------------

replace_ct <- function(x, ct) {
  x <- gsub("\\{ct\\}", ct, x)
  gsub("\\{ct\\}", ct, x)
}

canonical_gene_id <- function(x) {
  x <- trimws(as.character(x))
  is_ensg <- grepl("^ENSG[0-9]+(\\.[0-9]+)?$", x)
  x[is_ensg] <- sub("\\.[0-9]+$", "", x[is_ensg])
  x
}

flag_yes <- function(x) {
  toupper(trimws(as.character(x))) %in% c("YES", "TRUE", "1", "Y")
}

first_nonempty_column <- function(df, candidates) {
  candidates <- candidates[candidates %in% colnames(df)]
  out <- rep(NA_character_, nrow(df))

  if (!length(candidates)) {
    return(out)
  }

  for (nm in candidates) {
    v <- trimws(as.character(df[[nm]]))
    use <- (is.na(out) | !nzchar(out)) & !is.na(v) & nzchar(v)
    out[use] <- v[use]
  }

  out
}

canonical_condition <- function(x) {
  y <- toupper(trimws(as.character(x)))

  out <- rep(NA_character_, length(y))
  out[y == "AD"] <- "AD"
  out[y %in% SER_ALIASES] <- "SER"

  out
}

ensure_joined_counts <- function(obj) {

  DefaultAssay(obj) <- "RNA"
  layers <- SeuratObject::Layers(obj[["RNA"]])

  count_layers <- grep("^counts(\\.|$)", layers, value = TRUE)

  if (length(count_layers) > 1L) {
    message("  Joining ", length(count_layers), " RNA count layers.")
    obj <- SeuratObject::JoinLayers(obj, assay = "RNA")
    layers <- SeuratObject::Layers(obj[["RNA"]])
  }

  if (!"counts" %in% layers) {
    stop(
      "RNA counts layer is unavailable. Raw counts are required for pseudobulk."
    )
  }

  obj
}

# Aggregate a genes x cells sparse count matrix to genes x pseudobulk samples.
aggregate_sparse_counts <- function(counts, group) {

  group <- as.character(group)

  if (length(group) != ncol(counts)) {
    stop("Length of pseudobulk group vector does not match number of cells.")
  }

  if (anyNA(group) || any(!nzchar(group))) {
    stop("Missing/empty pseudobulk IDs reached aggregate_sparse_counts().")
  }

  f <- factor(group, levels = unique(group))

  design <- Matrix::sparseMatrix(
    i = seq_along(f),
    j = as.integer(f),
    x = 1,
    dims = c(length(f), nlevels(f)),
    dimnames = list(
      colnames(counts),
      levels(f)
    )
  )

  pb <- counts %*% design
  pb <- as(pb, "dgCMatrix")

  colnames(pb) <- levels(f)
  pb
}

normalize_pseudobulk <- function(pb_counts) {

  if (ncol(pb_counts) < 2L) {
    stop("At least 2 pseudobulk samples are needed for normalization.")
  }

  dge <- edgeR::DGEList(
    counts = as.matrix(pb_counts)
  )

  # TMM normalization corrects composition/library-size differences.
  dge <- edgeR::calcNormFactors(
    dge,
    method = "TMM"
  )

  # edgeR cpm(log=TRUE) is log2 CPM.
  log_cpm <- edgeR::cpm(
    dge,
    log = TRUE,
    prior.count = 1
  )

  # Quantile normalize COLUMNS = pseudobulk samples.
  qn <- preprocessCore::normalize.quantiles(
    as.matrix(log_cpm),
    copy = TRUE
  )

  rownames(qn) <- rownames(log_cpm)
  colnames(qn) <- colnames(log_cpm)

  list(
    logCPM = log_cpm,
    QN = qn,
    norm_factors = dge$samples
  )
}

read_tf_table <- function(path) {

  if (!file.exists(path)) {
    stop("51-TF CSV not found: ", path)
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

  miss <- setdiff(req, colnames(x))
  if (length(miss)) {
    stop(
      "TF file is missing required columns: ",
      paste(miss, collapse = ", ")
    )
  }

  x$`GLMM P values` <- as.numeric(x$`GLMM P values`)
  x$TF_ENSG <- canonical_gene_id(x$ENSG)
  x$TF_symbol <- trimws(as.character(x$`Transcription factor`))

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

  x <- x[!duplicated(x$TF_ENSG), , drop = FALSE]

  if (nrow(x) != 51L) {
    warning(
      "Expected 51 CRC TFs after GLMM P < 0.001, but found ",
      nrow(x), "."
    )
  }

  x
}

read_ct_risk_genes <- function(path, ct) {

  if (!file.exists(path)) {
    stop("TWAS significant-gene file not found for ", ct, ": ", path)
  }

  x <- data.table::fread(
    path,
    data.table = FALSE,
    check.names = FALSE
  )

  req <- c("gene", "EUR_sig", "Multi_sig")
  miss <- setdiff(req, colnames(x))

  if (length(miss)) {
    stop(
      "TWAS file for ", ct,
      " is missing required columns: ",
      paste(miss, collapse = ", ")
    )
  }

  x$gene <- canonical_gene_id(x$gene)

  symbol_candidates <- c(
    "gene_name_HTAN_ADSER_EUR",
    "gene_name_HTAN_ADSER_Multi",
    "gene_name_EUR",
    "gene_name_Multi",
    "gene_name",
    "Gene"
  )

  gene_symbol <- first_nonempty_column(
    x,
    symbol_candidates
  )

  gene_symbol[
    is.na(gene_symbol) | !nzchar(gene_symbol)
  ] <- x$gene[
    is.na(gene_symbol) | !nzchar(gene_symbol)
  ]

  eur_yes <- flag_yes(x$EUR_sig)
  multi_yes <- flag_yes(x$Multi_sig)

  support <- ifelse(
    eur_yes & multi_yes,
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
    grepl("^ENSG[0-9]+$", x$gene) &
    (eur_yes | multi_yes)

  out <- data.frame(
    risk_ENSG = x$gene[keep],
    risk_gene_symbol = gene_symbol[keep],
    EUR_sig = as.character(x$EUR_sig[keep]),
    Multi_sig = as.character(x$Multi_sig[keep]),
    GWAS_support = support[keep],
    stringsAsFactors = FALSE
  )

  # If duplicate rows exist, prioritize cross-GWAS support.
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
    dplyr::select(-.support_rank)

  out
}

# Map ENSG or symbol IDs to exact Seurat feature names.
make_feature_mapper <- function(feature_names) {

  feature_names <- as.character(feature_names)
  feature_keys <- canonical_gene_id(feature_names)

  ensg_lookup <- setNames(
    feature_names,
    feature_keys
  )
  ensg_lookup <- ensg_lookup[
    !duplicated(names(ensg_lookup))
  ]

  symbol_lookup <- setNames(
    feature_names,
    toupper(feature_names)
  )
  symbol_lookup <- symbol_lookup[
    !duplicated(names(symbol_lookup))
  ]

  function(ensg, symbol = NA_character_) {

    out <- unname(
      ensg_lookup[canonical_gene_id(ensg)]
    )

    miss <- is.na(out) | !nzchar(out)

    if (any(miss)) {
      sym <- toupper(trimws(as.character(symbol[miss])))
      out[miss] <- unname(
        symbol_lookup[sym]
      )
    }

    out
  }
}

correlation_from_expression <- function(
    expr,
    tf_features,
    tf_symbols,
    risk_features,
    risk_labels,
    ct,
    condition,
    min_samples = 5L) {

  n <- ncol(expr)

  # Always return a full 51 x risk-gene matrix, including NA when unavailable.
  cor_mat <- matrix(
    NA_real_,
    nrow = length(tf_features),
    ncol = length(risk_features),
    dimnames = list(
      tf_symbols,
      risk_labels
    )
  )

  p_mat <- cor_mat

  if (n < min_samples) {
    warning(
      ct, " ", condition, ": only ", n,
      " pseudobulk samples; correlations left as NA."
    )

    return(list(
      cor = cor_mat,
      p = p_mat,
      n = n
    ))
  }

  tf_present <- !is.na(tf_features) & tf_features %in% rownames(expr)
  risk_present <- !is.na(risk_features) & risk_features %in% rownames(expr)

  if (!any(tf_present) || !any(risk_present)) {
    return(list(
      cor = cor_mat,
      p = p_mat,
      n = n
    ))
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

  # Pearson correlation is undefined for genes with zero variance across
  # pseudobulk patients. Detect those genes explicitly instead of allowing
  # cor() to emit "the standard deviation is zero" warnings.
  tf_sd <- apply(tf_expr, 1L, stats::sd, na.rm = TRUE)
  risk_sd <- apply(risk_expr, 1L, stats::sd, na.rm = TRUE)

  tf_variable <- is.finite(tf_sd) & tf_sd > 0
  risk_variable <- is.finite(risk_sd) & risk_sd > 0

  if (any(!tf_variable)) {
    warning(
      ct, " ", condition, ": ",
      sum(!tf_variable),
      " TF(s) have zero/non-finite SD across pseudobulk patients; ",
      "their correlations are left as NA: ",
      paste(rownames(tf_expr)[!tf_variable], collapse = ", ")
    )
  }

  if (any(!risk_variable)) {
    warning(
      ct, " ", condition, ": ",
      sum(!risk_variable),
      " risk gene(s) have zero/non-finite SD across pseudobulk patients; ",
      "their correlations are left as NA: ",
      paste(rownames(risk_expr)[!risk_variable], collapse = ", ")
    )
  }

  tf_full_idx <- which(tf_present)[tf_variable]
  risk_full_idx <- which(risk_present)[risk_variable]

  if (length(tf_full_idx) > 0L && length(risk_full_idx) > 0L) {

    # cor() expects observations in rows, variables in columns.
    r_sub <- stats::cor(
      t(tf_expr[tf_variable, , drop = FALSE]),
      t(risk_expr[risk_variable, , drop = FALSE]),
      method = COR_METHOD,
      use = "pairwise.complete.obs"
    )

    # Insert valid correlations into the full 51 x risk-gene matrix.
    cor_mat[
      tf_full_idx,
      risk_full_idx
    ] <- r_sub

    # For complete normalized pseudobulk expression, n is common to all pairs.
    # Compute two-sided Pearson P values analytically from r.
    r_for_p <- r_sub

    denom <- pmax(
      1 - r_for_p^2,
      .Machine$double.eps
    )

    t_stat <- r_for_p * sqrt(
      (n - 2) / denom
    )

    p_sub <- 2 * stats::pt(
      -abs(t_stat),
      df = n - 2
    )

    p_sub[!is.finite(r_sub)] <- NA_real_

    p_mat[
      tf_full_idx,
      risk_full_idx
    ] <- p_sub
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
    condition,
    tf_meta,
    risk_meta,
    n_samples) {

  idx <- expand.grid(
    TF_row = seq_len(nrow(cor_mat)),
    gene_col = seq_len(ncol(cor_mat)),
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )

  out <- data.frame(
    cell_type = ct,
    condition = condition,
    TF = rownames(cor_mat)[idx$TF_row],
    TF_ENSG = tf_meta$TF_ENSG[idx$TF_row],
    risk_gene = colnames(cor_mat)[idx$gene_col],
    risk_gene_symbol = risk_meta$risk_gene_symbol[idx$gene_col],
    risk_gene_ENSG = risk_meta$risk_ENSG[idx$gene_col],
    EUR_sig = risk_meta$EUR_sig[idx$gene_col],
    Multi_sig = risk_meta$Multi_sig[idx$gene_col],
    GWAS_support = risk_meta$GWAS_support[idx$gene_col],
    r = cor_mat[
      cbind(idx$TF_row, idx$gene_col)
    ],
    P = p_mat[
      cbind(idx$TF_row, idx$gene_col)
    ],
    n_pseudobulk_samples = n_samples,
    stringsAsFactors = FALSE
  )

  out$FDR <- stats::p.adjust(
    out$P,
    method = "BH"
  )

  out
}

# Draw one clustered heatmap per cell type.
# Genes/TFs with completely undefined correlations are retained in the CSV
# outputs but excluded from the clustered heatmap because hclust() cannot use NA.
draw_single_heatmap <- function(
    cor_mat,
    ct,
    source_label,
    n_samples,
    outfile) {

  # --------------------------------------------------------------------------
  # Remove only rows/columns that contain no finite correlation at all.
  # This commonly occurs when a feature is absent from the RDS or has zero
  # variance across pseudobulk patients.
  # --------------------------------------------------------------------------
  row_has_finite <- apply(
    cor_mat,
    1L,
    function(x) any(is.finite(x))
  )

  col_has_finite <- apply(
    cor_mat,
    2L,
    function(x) any(is.finite(x))
  )

  dropped_rows <- rownames(cor_mat)[!row_has_finite]
  dropped_cols <- colnames(cor_mat)[!col_has_finite]

  if (length(dropped_rows)) {
    message(
      "  Heatmap ", ct, ": excluding ",
      length(dropped_rows),
      " TF row(s) with no finite correlations: ",
      paste(dropped_rows, collapse = ", ")
    )
  }

  if (length(dropped_cols)) {
    message(
      "  Heatmap ", ct, ": excluding ",
      length(dropped_cols),
      " risk-gene column(s) with no finite correlations: ",
      paste(dropped_cols, collapse = ", ")
    )
  }

  plot_mat <- cor_mat[
    row_has_finite,
    col_has_finite,
    drop = FALSE
  ]

  if (nrow(plot_mat) == 0L || ncol(plot_mat) == 0L) {
    warning(
      ct,
      ": no finite correlations are available for the heatmap; skipping PDF."
    )
    return(invisible(NULL))
  }

  # After removing completely undefined rows/columns, there normally should be
  # no NA values because normalized pseudobulk expression is complete.
  # If isolated non-finite cells remain, use 0 ONLY to construct dendrograms;
  # the displayed matrix itself keeps those cells as NA/grey.
  clustering_mat <- plot_mat

  if (any(!is.finite(clustering_mat))) {
    warning(
      ct,
      ": isolated NA/NaN/Inf correlations remain after removing all-NA ",
      "rows/columns. They will be represented as 0 for dendrogram calculation ",
      "only; displayed cells remain NA."
    )
    clustering_mat[!is.finite(clustering_mat)] <- 0
  }

  row_cluster <- if (nrow(plot_mat) > 1L) {
    stats::hclust(
      stats::dist(clustering_mat),
      method = "complete"
    )
  } else {
    FALSE
  }

  col_cluster <- if (ncol(plot_mat) > 1L) {
    stats::hclust(
      stats::dist(t(clustering_mat)),
      method = "complete"
    )
  } else {
    FALSE
  }

  # --------------------------------------------------------------------------
  # Draw heatmap first without writing directly to PDF
  # --------------------------------------------------------------------------
  ph <- pheatmap::pheatmap(
    plot_mat,
    color = HEATMAP_COLORS,
    breaks = HEATMAP_BREAKS,

    cluster_rows = row_cluster,
    cluster_cols = col_cluster,

    main = paste0(
      ct, ": 51 CRC TFs vs cell-type-specific risk genes\n",
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
  # Add outer margins around the complete pheatmap
  #
  # Increase these values if you want more white space.
  # Units are inches.
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
  # Save manually so margins can be controlled
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

  invisible(
    list(
      n_rows_plotted = nrow(plot_mat),
      n_cols_plotted = ncol(plot_mat),
      dropped_rows = dropped_rows,
      dropped_cols = dropped_cols
    )
  )
}


# Return the requested source conditions for each cell type.
conditions_for_celltype <- function(ct) {
  if (ct %in% SHARED_CELL_TYPES) {
    return(c("AD", "SER"))
  }
  if (ct == "ASC") {
    return("AD")
  }
  if (ct == "SSC") {
    return("SER")
  }
  stop("Unexpected cell type: ", ct)
}


analysis_label_for_celltype <- function(ct) {
  if (ct %in% SHARED_CELL_TYPES) {
    return("ADSER")
  }
  if (ct == "ASC") {
    return("AD")
  }
  if (ct == "SSC") {
    return("SER")
  }
  stop("Unexpected cell type: ", ct)
}


source_description_for_celltype <- function(ct) {
  if (ct %in% SHARED_CELL_TYPES) {
    return("AD + SER pooled before pseudobulk")
  }
  if (ct == "ASC") {
    return("AD only")
  }
  if (ct == "SSC") {
    return("SER only")
  }
  stop("Unexpected cell type: ", ct)
}


# ------------------------------------------------------------------------------
# 2. READ 51 CRC TFs
# ------------------------------------------------------------------------------

tf_meta <- read_tf_table(TF_CSV)

# Alphabetical TF order is fixed here and retained in every matrix/heatmap.
tf_meta <- tf_meta |>
  dplyr::arrange(TF_symbol)

message(
  "51-TF table: ",
  nrow(tf_meta),
  " TFs retained in alphabetical order."
)

data.table::fwrite(
  tf_meta,
  file.path(
    OUTDIR,
    "CRC51TFs_used_alphabetical.csv"
  )
)


# ------------------------------------------------------------------------------
# 3. MAIN LOOP OVER CELL TYPES
# ------------------------------------------------------------------------------

all_long <- list()
all_qc <- list()

for (ct in CELL_TYPES) {

  message("\n")
  message("======================================================================")
  message("CELL TYPE: ", ct)
  message("======================================================================")

  ct_outdir <- file.path(OUTDIR, ct)
  dir.create(
    ct_outdir,
    recursive = TRUE,
    showWarnings = FALSE
  )

  rds_file <- replace_ct(RDS_TEMPLATE, ct)
  twas_file <- replace_ct(TWAS_SIG_TEMPLATE, ct)

  message("RDS input: ", rds_file)
  message("TWAS input: ", twas_file)

  if (!file.exists(rds_file)) {
    warning("Missing cell-type RDS; skipping ", ct, ": ", rds_file)
    next
  }

  risk_meta <- read_ct_risk_genes(twas_file, ct)

  message(
    "TWAS risk genes for ", ct, ": ",
    nrow(risk_meta)
  )

  if (!nrow(risk_meta)) {
    warning("No significant TWAS risk genes for ", ct, "; skipping.")
    next
  }

  # Stable display labels.
  risk_labels <- ifelse(
    !is.na(risk_meta$risk_gene_symbol) &
      nzchar(risk_meta$risk_gene_symbol),
    risk_meta$risk_gene_symbol,
    risk_meta$risk_ENSG
  )

  risk_meta$risk_gene_label <- risk_labels

  # Alphabetical gene order. Use ENSG as a deterministic tie-breaker.
  risk_meta <- risk_meta |>
    dplyr::arrange(risk_gene_label, risk_ENSG)

  # If two ENSGs share the same displayed symbol, keep both labels unique.
  risk_meta$risk_gene_label <- make.unique(
    as.character(risk_meta$risk_gene_label)
  )

  risk_labels <- risk_meta$risk_gene_label

  data.table::fwrite(
    risk_meta,
    file.path(
      ct_outdir,
      paste0(ct, "_TWAS_risk_genes_used_alphabetical.csv")
    )
  )

  obj <- readRDS(rds_file)

  if (!inherits(obj, "Seurat")) {
    stop(ct, " RDS is not a Seurat object: ", rds_file)
  }

  obj <- ensure_joined_counts(obj)

  if (!CONDITION_COL %in% colnames(obj@meta.data)) {
    stop(
      CONDITION_COL,
      " is absent from metadata in ",
      rds_file
    )
  }

  if (!PSEUDOBULK_ID_COL %in% colnames(obj@meta.data)) {
    stop(
      PSEUDOBULK_ID_COL,
      " is absent from metadata in ",
      rds_file
    )
  }

  raw_condition <- as.character(
    obj@meta.data[[CONDITION_COL]]
  )
  canonical_cond <- canonical_condition(raw_condition)
  obj$condition_for_correlation <- canonical_cond

  if (any(
    toupper(trimws(raw_condition)) == "SET",
    na.rm = TRUE
  )) {
    message("  NOTE: sample_condition='SET' detected and treated as SER.")
  }

  message("Canonical AD/SER cell counts in input:")
  print(
    table(
      obj$condition_for_correlation,
      useNA = "ifany"
    )
  )

  wanted_conditions <- conditions_for_celltype(ct)
  analysis_label <- analysis_label_for_celltype(ct)
  source_description <- source_description_for_celltype(ct)

  # Shared lineages must actually have both AD and SER cells before pooling.
  if (ct %in% SHARED_CELL_TYPES) {
    n_ad_cells <- sum(
      obj$condition_for_correlation == "AD",
      na.rm = TRUE
    )
    n_ser_cells <- sum(
      obj$condition_for_correlation == "SER",
      na.rm = TRUE
    )

    if (n_ad_cells == 0L || n_ser_cells == 0L) {
      stop(
        "\nINPUT RDS CONDITION ERROR for ", ct, "\n",
        "  file: ", rds_file, "\n",
        "  AD cells : ", n_ad_cells, "\n",
        "  SER cells: ", n_ser_cells, "\n",
        "This analysis requires both AD and SER cells for the seven shared ",
        "cell types before pooling."
      )
    }
  }

  keep_cell <- !is.na(obj$condition_for_correlation) &
    obj$condition_for_correlation %in% wanted_conditions

  cells_use <- colnames(obj)[keep_cell]

  if (!length(cells_use)) {
    warning(
      ct, ": no cells available for requested source: ",
      paste(wanted_conditions, collapse = "+"),
      "; skipping."
    )
    next
  }

  md_use <- obj@meta.data[
    cells_use,
    ,
    drop = FALSE
  ]

  cond_use <- as.character(
    md_use$condition_for_correlation
  )

  pb_id <- trimws(
    as.character(
      md_use[[PSEUDOBULK_ID_COL]]
    )
  )

  valid_pb <- !is.na(pb_id) & nzchar(pb_id)

  cells_use <- cells_use[valid_pb]
  cond_use <- cond_use[valid_pb]
  pb_id <- pb_id[valid_pb]

  if (!length(cells_use)) {
    warning(ct, ": all selected cells had missing PatientID; skipping.")
    next
  }

  counts_all <- SeuratObject::GetAssayData(
    obj,
    assay = "RNA",
    layer = "counts"
  )

  mapper <- make_feature_mapper(
    rownames(counts_all)
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
    feature_in_RDS = tf_features,
    present = !is.na(tf_features),
    stringsAsFactors = FALSE
  )

  risk_presence <- data.frame(
    risk_gene = risk_labels,
    risk_gene_symbol = risk_meta$risk_gene_symbol,
    risk_gene_ENSG = risk_meta$risk_ENSG,
    GWAS_support = risk_meta$GWAS_support,
    feature_in_RDS = risk_features,
    present = !is.na(risk_features),
    stringsAsFactors = FALSE
  )

  data.table::fwrite(
    tf_presence,
    file.path(
      ct_outdir,
      paste0(ct, "_CRC51TF_feature_presence_alphabetical.csv")
    )
  )

  data.table::fwrite(
    risk_presence,
    file.path(
      ct_outdir,
      paste0(ct, "_TWASriskGene_feature_presence_alphabetical.csv")
    )
  )

  counts_use <- counts_all[
    ,
    cells_use,
    drop = FALSE
  ]

  # IMPORTANT:
  # For the seven shared lineages, this groups ALL selected AD+SER cells by
  # PatientID in a single operation. Therefore AD and SER are pooled BEFORE
  # pseudobulk normalization and correlation.
  pb_counts <- aggregate_sparse_counts(
    counts_use,
    pb_id
  )

  message(
    "  source: ", source_description,
    "; cells: ", length(cells_use),
    "; pseudobulk patients: ", ncol(pb_counts)
  )

  # Per-patient composition QC. For shared lineages this shows how many AD and
  # SER cells contributed to each pooled PatientID pseudobulk.
  cond_table <- table(
    factor(pb_id, levels = colnames(pb_counts)),
    factor(cond_use, levels = c("AD", "SER"))
  )

  pb_qc <- data.frame(
    cell_type = ct,
    analysis_group = analysis_label,
    pseudobulk_id = colnames(pb_counts),
    n_cells = as.integer(
      table(
        factor(
          pb_id,
          levels = colnames(pb_counts)
        )
      )
    ),
    n_AD_cells = as.integer(cond_table[, "AD"]),
    n_SER_cells = as.integer(cond_table[, "SER"]),
    library_size = Matrix::colSums(pb_counts),
    detected_genes = Matrix::colSums(pb_counts > 0),
    stringsAsFactors = FALSE
  )

  data.table::fwrite(
    pb_qc,
    file.path(
      ct_outdir,
      paste0(
        ct, "_", analysis_label,
        "_pseudobulk_QC.csv"
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

  pb_selected <- as.matrix(
    pb_counts[
      selected_features,
      ,
      drop = FALSE
    ]
  )

  data.table::fwrite(
    data.frame(
      gene_feature = rownames(pb_selected),
      pb_selected,
      check.names = FALSE
    ),
    file.path(
      ct_outdir,
      paste0(
        ct, "_", analysis_label,
        "_selectedGenes_pseudobulk_counts.csv.gz"
      )
    )
  )

  if (ncol(pb_counts) < 2L) {
    warning(
      ct, ": fewer than 2 pseudobulk samples; normalization/correlation skipped."
    )
    next
  }

  norm <- normalize_pseudobulk(
    pb_counts
  )

  norm_info <- data.frame(
    pseudobulk_id = rownames(norm$norm_factors),
    norm$norm_factors,
    row.names = NULL,
    check.names = FALSE
  )

  data.table::fwrite(
    norm_info,
    file.path(
      ct_outdir,
      paste0(
        ct, "_", analysis_label,
        "_TMM_normalization_factors.csv"
      )
    )
  )

  log_selected <- norm$logCPM[
    selected_features,
    ,
    drop = FALSE
  ]

  qn_selected <- norm$QN[
    selected_features,
    ,
    drop = FALSE
  ]

  data.table::fwrite(
    data.frame(
      gene_feature = rownames(log_selected),
      log_selected,
      check.names = FALSE
    ),
    file.path(
      ct_outdir,
      paste0(
        ct, "_", analysis_label,
        "_selectedGenes_TMM_log2CPM.csv.gz"
      )
    )
  )

  data.table::fwrite(
    data.frame(
      gene_feature = rownames(qn_selected),
      qn_selected,
      check.names = FALSE
    ),
    file.path(
      ct_outdir,
      paste0(
        ct, "_", analysis_label,
        "_selectedGenes_TMM_log2CPM_QN.csv.gz"
      )
    )
  )

  cor_res <- correlation_from_expression(
    expr = norm$QN,
    tf_features = tf_features,
    tf_symbols = tf_meta$TF_symbol,
    risk_features = risk_features,
    risk_labels = risk_labels,
    ct = ct,
    condition = analysis_label,
    min_samples = MIN_PB_SAMPLES
  )

  # The metadata tables were alphabetically sorted before mapping, so these
  # matrices are already alphabetically ordered. Enforce this as a check.
  stopifnot(
    identical(rownames(cor_res$cor), sort(rownames(cor_res$cor))),
    identical(colnames(cor_res$cor), sort(colnames(cor_res$cor)))
  )

  cor_file <- file.path(
    ct_outdir,
    paste0(
      ct, "_", analysis_label,
      "_CRC51TF_x_TWASriskGene_PearsonR_matrix_alphabetical.csv"
    )
  )

  data.table::fwrite(
    data.frame(
      TF = rownames(cor_res$cor),
      cor_res$cor,
      check.names = FALSE
    ),
    cor_file
  )

  data.table::fwrite(
    data.frame(
      TF = rownames(cor_res$p),
      cor_res$p,
      check.names = FALSE
    ),
    file.path(
      ct_outdir,
      paste0(
        ct, "_", analysis_label,
        "_CRC51TF_x_TWASriskGene_P_matrix_alphabetical.csv"
      )
    )
  )

  long_res <- matrix_to_long(
    cor_mat = cor_res$cor,
    p_mat = cor_res$p,
    ct = ct,
    condition = analysis_label,
    tf_meta = tf_meta,
    risk_meta = risk_meta,
    n_samples = cor_res$n
  )

  names(long_res)[names(long_res) == "condition"] <- "analysis_group"
  long_res$source_cells <- source_description

  # Keep output easy to inspect: alphabetical TF, then alphabetical risk gene.
  long_res <- long_res |>
    dplyr::arrange(TF, risk_gene)

  data.table::fwrite(
    long_res,
    file.path(
      ct_outdir,
      paste0(
        ct, "_", analysis_label,
        "_CRC51TF_x_TWASriskGene_correlation_long.csv"
      )
    )
  )

  all_long[[ct]] <- long_res

  draw_single_heatmap(
    cor_mat = cor_res$cor,
    ct = ct,
    source_label = source_description,
    n_samples = cor_res$n,
    outfile = file.path(
      ct_outdir,
      paste0(
        ct, "_", analysis_label,
        "_TF_TWASriskGene_correlation_heatmap_clustered.pdf"
      )
    )
  )

  all_qc[[ct]] <- data.frame(
    cell_type = ct,
    analysis_group = analysis_label,
    source_cells = source_description,
    n_input_cells = length(cells_use),
    n_AD_cells = sum(cond_use == "AD", na.rm = TRUE),
    n_SER_cells = sum(cond_use == "SER", na.rm = TRUE),
    n_pseudobulk = ncol(pb_counts),
    n_TWAS_risk_genes = nrow(risk_meta),
    n_TWAS_risk_genes_present = sum(!is.na(risk_features)),
    n_CRC51_TFs_present = sum(!is.na(tf_features)),
    correlation_computed = cor_res$n >= MIN_PB_SAMPLES,
    stringsAsFactors = FALSE
  )

  rm(
    obj,
    counts_all,
    counts_use,
    pb_counts
  )
  gc(verbose = FALSE)
}


# ------------------------------------------------------------------------------
# 4. COMBINE RESULTS ACROSS ALL CELL TYPES
# ------------------------------------------------------------------------------

if (length(all_long)) {

  all_long_df <- dplyr::bind_rows(
    all_long
  )

  data.table::fwrite(
    all_long_df,
    file.path(
      OUTDIR,
      "ALL9CT_oneCorrelation_CRC51TF_x_TWASriskGene_long.csv.gz"
    )
  )

  significant_correlations <- all_long_df |>
    dplyr::filter(
      !is.na(FDR),
      FDR < 0.05
    ) |>
    dplyr::arrange(
      cell_type,
      TF,
      risk_gene
    )

  data.table::fwrite(
    significant_correlations,
    file.path(
      OUTDIR,
      "ALL9CT_oneCorrelation_CRC51TF_x_TWASriskGene_FDR005.csv"
    )
  )
}


if (length(all_qc)) {

  all_qc_df <- dplyr::bind_rows(
    all_qc
  )

  data.table::fwrite(
    all_qc_df,
    file.path(
      OUTDIR,
      "ALL9CT_oneCorrelation_QC_summary.csv"
    )
  )

  print(all_qc_df)
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

message("\n======================================================================")
message("DONE")
message("======================================================================")
message("Results: ", OUTDIR)
message("")
message("Per cell type:")
message("  Shared 7 CTs: one AD+SER pooled correlation")
message("  ASC: one AD-only correlation")
message("  SSC: one SER-only correlation")
message("  TFs and risk genes: alphabetical order; no clustering")
