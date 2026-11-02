# =============================================================================
# GenoPhenoBasic.R
# Generic helper functions for TWAS model building
#
# Cleaned for the current PLINK .traw + PEER residual workflow:
#   - no VCF/vcf.head dependency
#   - no dataset-specific sample ID mapping
#   - no duplicated helper definitions
#   - annotation functions use their function arguments (no hidden globals)
# =============================================================================

suppressPackageStartupMessages({
  library(methods)
  library(stringr)
  library(data.table)
  library(dplyr)
})


# =============================================================================
# 1. Gene annotation
# =============================================================================

get_gene_annotation <- function(gene_annot_file_name, chrom) {

  gene_df <- read.table(
    gene_annot_file_name,
    header = FALSE,
    stringsAsFactors = FALSE,
    sep = "\t",
    fill = TRUE,
    quote = ""
  )

  gene_df1 <- gene_df[gene_df$V3 == "gene", , drop = FALSE]

  geneid <- stringr::str_extract(
    gene_df1[, 9],
    "ENSG\\d+(?:\\.\\d+)?"
  )

  genename <- stringr::str_match(
    gene_df1[, 9],
    'gene_name "([^"]+)"'
  )[, 2]

  gene_used <- data.frame(
    geneid = geneid,
    genename = genename,
    chr = gene_df1$V1,
    start = as.numeric(gene_df1$V4),
    end = as.numeric(gene_df1$V5),
    anno = gene_df1$V3,
    direction = gene_df1$V7,
    stringsAsFactors = FALSE
  )

  gene_used$gene_length <- gene_used$end - gene_used$start

  # Remove ENSG version suffix.
  gene_used$geneid <- sub(
    "^(ENSG[0-9]+)\\.[0-9]+$",
    "\\1",
    gene_used$geneid
  )

  # Keep requested chromosome.
  gene_used <- gene_used[
    gene_used$chr == paste0("chr", chrom),
    ,
    drop = FALSE
  ]

  # For duplicated gene symbols, keep the longest gene record.
  gene_used <- gene_used[
    order(gene_used$genename, -gene_used$gene_length),
    ,
    drop = FALSE
  ]

  gene_used <- gene_used[
    !duplicated(gene_used$genename),
    ,
    drop = FALSE
  ]

  rownames(gene_used) <- NULL
  gene_used
}


get_gene_annotation_all <- function(gene_annot_file_name) {

  gene_df <- read.table(
    gene_annot_file_name,
    header = FALSE,
    stringsAsFactors = FALSE,
    sep = "\t",
    fill = TRUE,
    quote = ""
  )

  gene_df1 <- gene_df[gene_df$V3 == "gene", , drop = FALSE]

  geneid <- stringr::str_extract(
    gene_df1[, 9],
    "ENSG\\d+(?:\\.\\d+)?"
  )

  genename <- stringr::str_match(
    gene_df1[, 9],
    'gene_name "([^"]+)"'
  )[, 2]

  gene_used <- data.frame(
    geneid = geneid,
    genename = genename,
    chr = gene_df1$V1,
    start = as.numeric(gene_df1$V4),
    end = as.numeric(gene_df1$V5),
    anno = gene_df1$V3,
    direction = gene_df1$V7,
    stringsAsFactors = FALSE
  )

  gene_used$gene_length <- gene_used$end - gene_used$start

  gene_used$geneid <- sub(
    "^(ENSG[0-9]+)\\.[0-9]+$",
    "\\1",
    gene_used$geneid
  )

  gene_used <- gene_used[
    order(gene_used$genename, -gene_used$gene_length),
    ,
    drop = FALSE
  ]

  gene_used <- gene_used[
    !duplicated(gene_used$genename),
    ,
    drop = FALSE
  ]

  rownames(gene_used) <- NULL
  gene_used
}


# =============================================================================
# 2. Generic expression reader
#
# Expected input:
#   rows    = genes
#   columns = samples
#   first column = gene ID or gene symbol
#
# Returns:
#   rows    = samples
#   columns = GENCODE ENSG IDs
#
# NOTE:
#   expression_id -> genotype_id conversion is intentionally NOT done here.
#   That mapping belongs in the dataset-specific TWAS_modelbuilding.R script.
# =============================================================================

get_gene_expression <- function(expression_file, gene_annot) {

  x <- data.table::fread(
    expression_file,
    data.table = FALSE,
    check.names = FALSE
  )

  if (ncol(x) < 2) {
    stop(
      "Expression file must contain one gene column and at least one sample column."
    )
  }

  input_gene <- trimws(as.character(x[[1]]))
  sample_ids <- colnames(x)[-1]

  if (anyDuplicated(sample_ids)) {
    dup <- unique(sample_ids[duplicated(sample_ids)])
    stop(
      "Duplicated expression sample IDs: ",
      paste(head(dup, 20), collapse = ", ")
    )
  }

  if (!all(c("geneid", "genename") %in% colnames(gene_annot))) {
    stop("gene_annot must contain columns: geneid and genename.")
  }

  annot_geneid <- as.character(gene_annot$geneid)
  annot_genename <- as.character(gene_annot$genename)

  # 1. Exact ENSG match.
  idx <- match(input_gene, annot_geneid)

  # 2. ENSG match ignoring version suffix.
  need <- is.na(idx) & grepl("^ENSG[0-9]+", input_gene)

  if (any(need)) {
    input_base <- sub("\\.[0-9]+$", "", input_gene[need])
    annot_base <- sub("\\.[0-9]+$", "", annot_geneid)
    idx[need] <- match(input_base, annot_base)
  }

  # 3. Gene-symbol match.
  need <- is.na(idx)

  if (any(need)) {
    idx[need] <- match(input_gene[need], annot_genename)
  }

  keep <- !is.na(idx)

  if (!any(keep)) {
    stop(
      "No expression genes matched the supplied gene annotation."
    )
  }

  x <- x[keep, , drop = FALSE]
  canonical_geneid <- annot_geneid[idx[keep]]

  # Keep one row per canonical gene ID.
  nondup <- !duplicated(canonical_geneid)
  x <- x[nondup, , drop = FALSE]
  canonical_geneid <- canonical_geneid[nondup]

  expr_mat <- as.matrix(x[, -1, drop = FALSE])
  storage.mode(expr_mat) <- "numeric"

  if (any(!is.finite(expr_mat))) {
    stop(
      "Expression matrix contains NA/NaN/Inf."
    )
  }

  # genes x samples -> samples x genes
  expr_df <- as.data.frame(
    t(expr_mat),
    check.names = FALSE
  )

  rownames(expr_df) <- sample_ids
  colnames(expr_df) <- canonical_geneid

  expr_df
}


# =============================================================================
# 3. Gene coordinates
# =============================================================================

get_gene_coords <- function(gene_annot, gene) {

  if (grepl("^ENSG", gene)) {
    row <- gene_annot[
      gene_annot$geneid == gene,
      ,
      drop = FALSE
    ]
  } else {
    row <- gene_annot[
      gene_annot$genename == gene,
      ,
      drop = FALSE
    ]
  }

  if (nrow(row) == 0) {
    warning("Gene not found in annotation: ", gene)
    return(c(NA_real_, NA_real_))
  }

  # TSS according to strand.
  if (row$direction[1] == "+") {
    tss <- row$start[1]
  } else if (row$direction[1] == "-") {
    tss <- row$end[1]
  } else {
    warning(
      "Unknown strand for gene ", gene,
      "; using start coordinate as TSS."
    )
    tss <- row$start[1]
  }

  c(tss, tss)
}


# =============================================================================
# 4. Cis-genotype extraction
# =============================================================================

get_cis_genotype <- function(
  gt_df,
  snp_annot,
  pos_col_name,
  coords,
  cis_window
) {

  if (anyNA(coords)) {
    return(
      matrix(
        numeric(0),
        nrow = nrow(gt_df),
        ncol = 0,
        dimnames = list(rownames(gt_df), NULL)
      )
    )
  }

  if (!pos_col_name %in% colnames(snp_annot)) {
    stop(
      "Position column not found in SNP annotation: ",
      pos_col_name
    )
  }

  if (!"SNP" %in% colnames(snp_annot)) {
    stop("SNP annotation does not contain column 'SNP'.")
  }

  snp_info <- snp_annot %>%
    filter(
      .data[[pos_col_name]] >= (coords[1] - cis_window),
      .data[[pos_col_name]] <= (coords[2] + cis_window)
    )

  keep <- intersect(colnames(gt_df), snp_info$SNP)

  as.matrix(
    gt_df[, keep, drop = FALSE]
  )
}


get_cis_genotype_varID <- function(
  gt_df,
  snp_annot,
  pos_col_name,
  coords,
  cis_window
) {

  if (anyNA(coords)) {
    return(
      matrix(
        numeric(0),
        nrow = nrow(gt_df),
        ncol = 0,
        dimnames = list(rownames(gt_df), NULL)
      )
    )
  }

  if (!pos_col_name %in% colnames(snp_annot)) {
    stop(
      "Position column not found in SNP annotation: ",
      pos_col_name
    )
  }

  if (!"varID" %in% colnames(snp_annot)) {
    stop("SNP annotation does not contain column 'varID'.")
  }

  snp_info <- snp_annot %>%
    filter(
      .data[[pos_col_name]] >= (coords[1] - cis_window),
      .data[[pos_col_name]] <= (coords[2] + cis_window)
    )

  keep <- intersect(colnames(gt_df), snp_info$varID)

  as.matrix(
    gt_df[, keep, drop = FALSE]
  )
}


# =============================================================================
# 5. Genotype covariance for prediction models
# =============================================================================

do_covariance <- function(
  gene_id,
  cis_gt,
  rsids,
  varIDs,
  out_file
) {

  if (length(rsids) != length(varIDs)) {
    stop("rsids and varIDs must have the same length.")
  }

  if (length(varIDs) == 0) {
    return(invisible(NULL))
  }

  missing_varids <- setdiff(varIDs, colnames(cis_gt))

  if (length(missing_varids) > 0) {
    stop(
      "Model varID(s) not found in cis genotype matrix: ",
      paste(head(missing_varids, 20), collapse = ", ")
    )
  }

  model_gt <- cis_gt[
    ,
    varIDs,
    drop = FALSE
  ]

  geno_cov <- stats::cov(
    model_gt,
    use = "pairwise.complete.obs"
  )

  out <- vector(
    "list",
    length(varIDs) * (length(varIDs) + 1) / 2
  )

  k <- 1L

  for (i in seq_along(varIDs)) {
    for (j in i:length(varIDs)) {

      out[[k]] <- data.frame(
        gene = gene_id,
        rsid1 = rsids[i],
        rsid2 = rsids[j],
        covariance = geno_cov[i, j],
        stringsAsFactors = FALSE
      )

      k <- k + 1L
    }
  }

  cov_df <- do.call(
    rbind,
    out[seq_len(k - 1L)]
  )

  write.table(
    cov_df,
    file = out_file,
    append = TRUE,
    quote = FALSE,
    col.names = FALSE,
    row.names = FALSE,
    sep = "\t"
  )

  invisible(cov_df)
}


# =============================================================================
# 6. Small utility
# =============================================================================

check_unique <- function(column) {
  length(unique(column)) == 1L
}
