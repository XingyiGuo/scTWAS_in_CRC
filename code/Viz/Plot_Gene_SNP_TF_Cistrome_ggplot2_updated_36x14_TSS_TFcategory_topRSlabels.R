#!/usr/bin/env Rscript

# ============================================================================
# Gene–variant–TF occupancy illustration using ggplot2
#
# This script plots, for one gene in one cell type (or for all CellType-gene
# pairs):
#   - the genomic axis
#   - the gene body as a rectangle
#   - model SNPs at their hg38 genomic positions
#   - stacked TF symbols for Cistrome occupancy
#     (same colour = same functional program; different shape = different TF)
#   - separated rsID labels connected back to SNP positions by line segments
#   - TSS annotation, when gene_direction/strand is available
#
# TSS definition:
#   gene_direction == "+"  -> TSS = gene_start
#   gene_direction == "-"  -> TSS = gene_end
#
# Note:
#   The input file should ideally contain:
#     CellType, gene, genename, gene_chr, gene_start, gene_end,
#     gene_direction, rs, GWAS_SNP_CHR_POS,
#     Cistrome_overlap, N_Cistrome_TFs, Cistrome_TFs
#
#   If gene_direction is missing, the script will still run, but TSS will not
#   be displayed.
# ============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(scales)
})

set.seed(20260909)

# ============================================================================
# 1. Settings
# ============================================================================

INPUT_FILE <- paste0(
  "/data/l2_bioinfo1/liq17/scTWAS/",
  "TWAS_HTAN_GSE266546_WellcomeSanger_res/",
  "SPrediXcan_peakSNP_association/",
  "VUMC_HTAN_ADSER_EUR_Multi_sig/",
  "CRC51TFs_occupy_model_snps/",
  "VUMC_HTAN_ADSER_9CTs_sigGenes_modelSNPs_CistromeAnnotated.csv"
)

OUT_DIR <- file.path(
  dirname(INPUT_FILE),
  "Gene_SNP_TF_illustrations"
)

dir.create(
  OUT_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)

# --------------------------------------------------------------------------
# Run mode:
#   "single" = plot one selected gene in one selected cell type
#   "all"    = plot every CellType-gene pair in the input file
# --------------------------------------------------------------------------
RUN_MODE <- "single"

# Used only when RUN_MODE == "single"
TARGET_CELLTYPE <- "ABS"
TARGET_GENE <- "LAMC1"   # gene symbol or ENSG ID

# --------------------------------------------------------------------------
# Plotting options
# --------------------------------------------------------------------------
SHOW_ALL_MODEL_SNPS <- TRUE
SHOW_SNP_LABELS <- TRUE
# rs labels are placed directly above the occupied SNP TF stacks, so connectors
# are no longer needed.
SHOW_SNP_LABEL_CONNECTORS <- FALSE
SHOW_TSS <- TRUE

# Optional vector of TFs to display. NULL = display all TFs in Cistrome_TFs.
TF_KEEP <- NULL

# Geometry / layout
# Raise the first TF symbol above the genomic axis/gene annotation region.
TF_BASE_Y <- 0.86
TF_VERTICAL_SPACING <- 0.30
TF_POINT_SIZE <- 2.7
SNP_POINT_SIZE <- 1.7

# Display occupied SNPs with LOCAL spreading rather than distributing them
# uniformly across the whole locus. Isolated SNPs remain near their true
# genomic positions; only crowded neighboring SNPs are separated.
SPREAD_OCCUPIED_SNPS <- TRUE
# Slightly smaller minimum gap keeps neighboring occupied SNP stacks closer
# together while remaining readable.
OCCUPIED_SNP_MIN_GAP_FRAC <- 0.026
# Keep stacks away from the extreme plot edges.
OCCUPIED_SNP_BOUNDARY_FRAC <- 0.030

# Extra padding at the locus edges helps prevent clipping of the right-most
# TF stacks and labels.
AXIS_PADDING_FRAC <- 0.06
AXIS_PADDING_MIN_BP <- 100

# rs labels are now placed above the TF stacks.
SNP_LABEL_STACK_LEVELS <- 6L
SNP_LABEL_ANGLE <- 90
SNP_LABEL_SIZE <- 2.7
SNP_LABEL_TOP_PADDING <- 0.24
SNP_LABEL_CONNECTOR_X_OFFSET <- 9000

# Custom TF legend:
# show all TF names that occupy SNPs in the current plot, grouped by category,
# with the actual plotting symbol displayed before each TF name.
SHOW_CUSTOM_TF_TEXT_LEGEND <- TRUE
CUSTOM_TF_LEGEND_TITLE_SIZE <- 3.6
CUSTOM_TF_LEGEND_CATEGORY_LABEL_SIZE <- 3.0
CUSTOM_TF_LEGEND_TEXT_SIZE <- 2.8
CUSTOM_TF_LEGEND_POINT_SIZE <- 2.5
CUSTOM_TF_LEGEND_LINE_SPACING <- 0.62
CUSTOM_TF_LEGEND_Y_PADDING <- 0.42
CUSTOM_TF_LEGEND_LEFT_FRAC <- 0.31
CUSTOM_TF_LEGEND_RIGHT_FRAC <- 0.99
CUSTOM_TF_LEGEND_CATEGORY_BLOCK_FRAC <- 0.19
CUSTOM_TF_LEGEND_CHAR_WIDTH_FRAC <- 0.0038
CUSTOM_TF_LEGEND_POINT_TEXT_GAP_FRAC <- 0.0052
CUSTOM_TF_LEGEND_ITEM_GAP_FRAC <- 0.011

# Gene/TSS annotation sizes and vertical positions.
# These are intentionally below the gene box/baseline.
TSS_LABEL_SIZE <- 3.2
GENE_LABEL_SIZE <- 5.8
GENE_COORD_LABEL_SIZE <- 3.1
GENE_COORD_Y <- -0.31
TSS_LABEL_Y <- -0.62

# Fixed extra-wide manuscript layout
PLOT_WIDTH <- 36
PLOT_HEIGHT <- 14

# The default ggplot legend is suppressed.
# A custom symbol + text legend is drawn inside the panel:
# one row per functional category, with each TF shown by its symbol and name.
LEGEND_POSITION <- "none"
LEGEND_NCOL <- 4L

# --------------------------------------------------------------------------
# CRC 51-TF/regulator functional-program classification
#
# Each regulator is assigned to ONE primary category for visualization.
# This makes the colour mapping mutually exclusive and stable across genes.
# Biological functions can overlap in reality; these are primary labels chosen
# for the locus-level visualization.
# --------------------------------------------------------------------------

TF_CATEGORY_LEVELS <- c(
  "Intestinal lineage / WNT / differentiation",
  "Proliferation / oncogenic growth",
  "Inflammation / stress / AP-1",
  "Metabolic / nuclear-receptor signaling",
  "Development / lineage plasticity",
  "Chromatin / enhancer / 3D genome",
  "General transcription / regulation",
  "Other sequence-specific regulators"
)

TF_CLASSIFICATION <- rbindlist(list(

  data.table(
    TF = c(
      "CDX2", "FOXA2", "GATA6", "HNF4A", "KLF5",
      "PROX1", "SOX9", "TCF7L1", "TCF7L2"
    ),
    TF_category = "Intestinal lineage / WNT / differentiation"
  ),

  data.table(
    TF = c(
      "E2F3", "E2F8", "ETS2", "ETV5",
      "MAX", "MYC", "TFAP4"
    ),
    TF_category = "Proliferation / oncogenic growth"
  ),

  data.table(
    TF = c(
      "ATF2", "ATF3", "CEBPB", "EGR1", "ELF1",
      "FOSL1", "IRF1", "JUN", "JUND", "RELA", "ZBTB7B"
    ),
    TF_category = "Inflammation / stress / AP-1"
  ),

  data.table(
    TF = c("NR1H2", "NR1H3", "PPARG", "VDR"),
    TF_category = "Metabolic / nuclear-receptor signaling"
  ),

  data.table(
    TF = c("HOXA6", "SOX2", "SRF"),
    TF_category = "Development / lineage plasticity"
  ),

  data.table(
    TF = c(
      "CBX3", "EP300", "KMT2B", "NIPBL",
      "RAD21", "SIN3A", "SMC3", "YY1"
    ),
    TF_category = "Chromatin / enhancer / 3D genome"
  ),

  data.table(
    TF = c(
      "CNOT3", "ICE1", "MED12", "PAF1",
      "POLR2A", "SP1", "USF1"
    ),
    TF_category = "General transcription / regulation"
  ),

  data.table(
    TF = c("REPIN1", "ZNF83"),
    TF_category = "Other sequence-specific regulators"
  )
))

# Expected complete set of 51 CRC TFs/regulators.
CRC51_TFS <- c(
  "ATF2", "ATF3", "CBX3", "CDX2", "CEBPB", "CNOT3", "E2F3", "E2F8",
  "EGR1", "ELF1", "EP300", "ETS2", "ETV5", "FOSL1", "FOXA2", "GATA6",
  "HNF4A", "HOXA6", "ICE1", "IRF1", "JUN", "JUND", "KLF5", "KMT2B",
  "MAX", "MED12", "MYC", "NIPBL", "NR1H2", "NR1H3", "PAF1", "POLR2A",
  "PPARG", "PROX1", "RAD21", "RELA", "REPIN1", "SIN3A", "SMC3", "SOX2",
  "SOX9", "SP1", "SRF", "TCF7L1", "TCF7L2", "TFAP4", "USF1", "VDR",
  "YY1", "ZBTB7B", "ZNF83"
)

# Validate that every one of the 51 regulators occurs exactly once.
if (anyDuplicated(TF_CLASSIFICATION$TF)) {
  stop(
    "Duplicate TF(s) in TF_CLASSIFICATION: ",
    paste(
      unique(TF_CLASSIFICATION$TF[duplicated(TF_CLASSIFICATION$TF)]),
      collapse = ", "
    )
  )
}

missing_from_classification <- setdiff(CRC51_TFS, TF_CLASSIFICATION$TF)
extra_in_classification <- setdiff(TF_CLASSIFICATION$TF, CRC51_TFS)

if (length(missing_from_classification) > 0L) {
  stop(
    "TF_CLASSIFICATION is missing CRC51 TF(s): ",
    paste(missing_from_classification, collapse = ", ")
  )
}

if (length(extra_in_classification) > 0L) {
  stop(
    "TF_CLASSIFICATION contains unexpected TF(s): ",
    paste(extra_in_classification, collapse = ", ")
  )
}

if (nrow(TF_CLASSIFICATION) != 51L) {
  stop(
    "TF_CLASSIFICATION should contain 51 rows, but contains ",
    nrow(TF_CLASSIFICATION),
    "."
  )
}

TF_CLASSIFICATION[
  ,
  TF_category := factor(
    TF_category,
    levels = TF_CATEGORY_LEVELS
  )
]

# Stable category colours across all genes/cell types.
TF_CATEGORY_COLORS <- c(
  "Intestinal lineage / WNT / differentiation" = "#0072B2",
  "Proliferation / oncogenic growth" = "#D55E00",
  "Inflammation / stress / AP-1" = "#CC79A7",
  "Metabolic / nuclear-receptor signaling" = "#009E73",
  "Development / lineage plasticity" = "#E69F00",
  "Chromatin / enhancer / 3D genome" = "#6A3D9A",
  "General transcription / regulation" = "#56B4E9",
  "Other sequence-specific regulators" = "#7F7F7F"
)

# Shapes repeat ACROSS categories, but are unique WITHIN each category.
# This keeps the figure readable while allowing individual TFs in the same
# functional program to be visually distinguished.
TF_SHAPE_POOL <- c(
  16, 17, 15, 18, 3, 4, 8, 7, 9, 10, 0, 1, 2, 5, 6, 11, 12, 13, 14
)

TF_CLASSIFICATION[
  ,
  TF_shape := TF_SHAPE_POOL[seq_len(.N)],
  by = TF_category
]

TF_SHAPE_VALUES <- setNames(
  TF_CLASSIFICATION$TF_shape,
  TF_CLASSIFICATION$TF
)

# Write the classification used by the plots for reproducibility.
fwrite(
  TF_CLASSIFICATION[, .(
    TF,
    TF_category = as.character(TF_category),
    TF_shape
  )],
  file.path(
    OUT_DIR,
    "CRC51_TF_functional_program_classification.csv"
  )
)

SAVE_PDF <- TRUE
SAVE_PNG <- TRUE
PNG_DPI <- 400

# ============================================================================
# 2. Read and validate input
# ============================================================================

if (!file.exists(INPUT_FILE)) {
  stop("Input file not found: ", INPUT_FILE)
}

dat <- fread(
  INPUT_FILE,
  na.strings = c("", "NA", "NAN"),
  showProgress = FALSE
)

required_cols <- c(
  "CellType",
  "gene",
  "genename",
  "gene_chr",
  "gene_start",
  "gene_end",
  "rs",
  "GWAS_SNP_CHR_POS",
  "Cistrome_overlap",
  "N_Cistrome_TFs",
  "Cistrome_TFs"
)

missing_cols <- setdiff(required_cols, names(dat))
if (length(missing_cols) > 0L) {
  stop(
    "Input file is missing required column(s): ",
    paste(missing_cols, collapse = ", ")
  )
}

# gene_direction is optional; create it if absent
if (!"gene_direction" %in% names(dat)) {
  dat[, gene_direction := NA_character_]
}

# ============================================================================
# 3. Standardize fields
# ============================================================================

dat[, CellType := toupper(trimws(as.character(CellType)))]
dat[, gene := trimws(as.character(gene))]
dat[, genename := trimws(as.character(genename))]
dat[, rs := trimws(as.character(rs))]

# Gene chromosome
dat[, gene_chr := as.character(gene_chr)]
dat[, gene_chr := sub("^chr", "", gene_chr, ignore.case = TRUE)]
dat[, gene_chr := paste0("chr", gene_chr)]

# Gene positions
dat[, gene_start := as.numeric(gene_start)]
dat[, gene_end := as.numeric(gene_end)]

# Gene direction / strand
dat[, gene_direction := trimws(as.character(gene_direction))]
dat[
  gene_direction %in% c("1", "plus", "PLUS"),
  gene_direction := "+"
]
dat[
  gene_direction %in% c("-1", "minus", "MINUS"),
  gene_direction := "-"
]

# Cistrome overlap
dat[, Cistrome_overlap := toupper(
  trimws(as.character(Cistrome_overlap))
) == "TRUE"]

dat[, N_Cistrome_TFs := as.integer(N_Cistrome_TFs)]

# Parse SNP coordinates from GWAS_SNP_CHR_POS
dat[, SNP_chr := sub(
  ":.*$",
  "",
  as.character(GWAS_SNP_CHR_POS)
)]
dat[, SNP_chr := sub("^chr", "", SNP_chr, ignore.case = TRUE)]
dat[, SNP_chr := paste0("chr", SNP_chr)]

dat[, SNP_pos := suppressWarnings(
  as.numeric(
    sub(
      "^[^:]+:",
      "",
      as.character(GWAS_SNP_CHR_POS)
    )
  )
)]

# ============================================================================
# 4. Expand semicolon-separated TF annotations
# ============================================================================

expand_tf_annotations <- function(d) {

  d_tf <- d[
    Cistrome_overlap == TRUE &
      !is.na(Cistrome_TFs) &
      trimws(Cistrome_TFs) != ""
  ]

  if (nrow(d_tf) == 0L) {
    return(
      data.table(
        rs = character(),
        SNP_chr = character(),
        SNP_pos = numeric(),
        TF = character()
      )
    )
  }

  tf_long <- d_tf[
    ,
    {
      tf_vec <- trimws(
        unlist(
          strsplit(
            as.character(Cistrome_TFs[1]),
            ";",
            fixed = TRUE
          )
        )
      )

      tf_vec <- tf_vec[
        !is.na(tf_vec) &
          tf_vec != ""
      ]

      .(TF = unique(tf_vec))
    },
    by = .(
      rs,
      SNP_chr,
      SNP_pos
    )
  ]

  unique(tf_long)
}

# ============================================================================
# 5. Main plotting function
# ============================================================================

plot_gene_snp_tf <- function(
  x,
  target_gene,
  target_celltype,
  show_all_model_snps = TRUE,
  tf_keep = NULL
) {

  target_celltype <- toupper(target_celltype)

  # ------------------------------------------------------------------------
  # Select one gene / one cell type.
  # target_gene can be either ENSG or gene symbol.
  # ------------------------------------------------------------------------
  d <- x[
    CellType == target_celltype &
      (gene == target_gene | genename == target_gene)
  ]

  if (nrow(d) == 0L) {
    stop(
      "No rows found for CellType=",
      target_celltype,
      " and gene=",
      target_gene
    )
  }

  # ------------------------------------------------------------------------
  # Gene identity / coordinates / strand
  # ------------------------------------------------------------------------
  gene_ids <- unique(
    d[!is.na(gene) & gene != "", gene]
  )

  gene_symbols <- unique(
    d[!is.na(genename) & genename != "", genename]
  )

  gene_chrs <- unique(
    d[!is.na(gene_chr) & gene_chr != "", gene_chr]
  )

  gene_starts <- unique(
    d[!is.na(gene_start), gene_start]
  )

  gene_ends <- unique(
    d[!is.na(gene_end), gene_end]
  )

  gene_directions <- unique(
    d[
      !is.na(gene_direction) &
        gene_direction %in% c("+", "-"),
      gene_direction
    ]
  )

  if (
    length(gene_chrs) != 1L ||
      length(gene_starts) == 0L ||
      length(gene_ends) == 0L
  ) {
    stop(
      "Could not determine a unique gene chromosome/start/end for ",
      target_gene,
      " in ",
      target_celltype,
      "."
    )
  }

  gene_id <- if (length(gene_ids) > 0L) gene_ids[1] else target_gene
  gene_symbol <- if (length(gene_symbols) > 0L) gene_symbols[1] else target_gene
  gene_chr <- gene_chrs[1]
  gene_start <- gene_starts[1]
  gene_end <- gene_ends[1]
  gene_direction <- if (length(gene_directions) == 1L) gene_directions[1] else NA_character_

  # TSS
  tss_pos <- NA_real_
  if (!is.na(gene_direction) && gene_direction == "+") {
    tss_pos <- gene_start
  } else if (!is.na(gene_direction) && gene_direction == "-") {
    tss_pos <- gene_end
  }

  # ------------------------------------------------------------------------
  # One row per SNP
  # ------------------------------------------------------------------------
  snps <- unique(
    d[
      !is.na(SNP_pos) &
        !is.na(rs) &
        rs != "",
      .(
        rs,
        SNP_chr,
        SNP_pos,
        weight,
        Cistrome_overlap,
        N_Cistrome_TFs
      )
    ],
    by = c("rs", "SNP_chr", "SNP_pos")
  )

  if (nrow(snps) == 0L) {
    stop(
      "No SNPs with usable GWAS_SNP_CHR_POS coordinates were found for ",
      gene_symbol,
      " / ",
      target_celltype,
      "."
    )
  }

  # Optional: keep all model SNPs (default behavior)
  if (!show_all_model_snps) {
    snps <- snps[
      Cistrome_overlap == TRUE
    ]
  }

  # ------------------------------------------------------------------------
  # Expand TF occupancy annotations
  # ------------------------------------------------------------------------
  tf_long <- expand_tf_annotations(d)

  if (!is.null(tf_keep)) {
    tf_keep <- unique(trimws(as.character(tf_keep)))
    tf_long <- tf_long[TF %chin% tf_keep]
  }

  # Keep TF rows only for SNPs that are displayed
  tf_long <- tf_long[
    rs %chin% snps$rs
  ]

  # ------------------------------------------------------------------------
  # Add mutually exclusive TF functional-program category and TF-specific shape
  # ------------------------------------------------------------------------
  if (nrow(tf_long) > 0L) {

    unclassified_tfs <- setdiff(
      unique(tf_long$TF),
      TF_CLASSIFICATION$TF
    )

    if (length(unclassified_tfs) > 0L) {
      warning(
        "The following Cistrome TFs are not in the predefined CRC51 ",
        "classification and will be omitted from the plot: ",
        paste(unclassified_tfs, collapse = ", ")
      )

      tf_long <- tf_long[
        TF %chin% TF_CLASSIFICATION$TF
      ]
    }

    tf_long <- merge(
      tf_long,
      TF_CLASSIFICATION[, .(
        TF,
        TF_category,
        TF_shape
      )],
      by = "TF",
      all.x = TRUE,
      sort = FALSE
    )

    tf_long[
      ,
      TF_category := factor(
        TF_category,
        levels = TF_CATEGORY_LEVELS
      )
    ]
  }

  occupied_rs <- unique(tf_long$rs)

  occupied_snps <- snps[
    rs %chin% occupied_rs
  ]

  # ------------------------------------------------------------------------
  # Stack TF symbols vertically at each displayed SNP.
  # Sort by functional program first so same-colour TFs appear together.
  # ------------------------------------------------------------------------
  if (nrow(tf_long) > 0L) {

    setorder(tf_long, SNP_pos, TF_category, TF)

    tf_long[
      ,
      TF_stack_index := seq_len(.N),
      by = .(rs, SNP_pos)
    ]

    tf_long[
      ,
      TF_y := TF_BASE_Y + (TF_stack_index - 1L) * TF_VERTICAL_SPACING
    ]

    max_tf_stack <- max(tf_long$TF_stack_index, na.rm = TRUE)

  } else {

    tf_long[, `:=`(
      TF_stack_index = integer(),
      TF_y = numeric()
    )]

    max_tf_stack <- 0L
  }

  # ------------------------------------------------------------------------
  # X range
  # ------------------------------------------------------------------------
  x_min <- min(
    c(gene_start, snps$SNP_pos),
    na.rm = TRUE
  )

  x_max <- max(
    c(gene_end, snps$SNP_pos),
    na.rm = TRUE
  )

  locus_width <- x_max - x_min
  if (!is.finite(locus_width) || locus_width <= 0) {
    locus_width <- 1000
  }

  padding <- max(locus_width * AXIS_PADDING_FRAC, AXIS_PADDING_MIN_BP)
  axis_start <- x_min - padding
  axis_end <- x_max + padding

  # ------------------------------------------------------------------------
  # Occupied SNP display positions:
  # keep the true SNP genomic positions on the baseline, but optionally spread
  # the displayed occupied-SNP stacks more sparsely across the x-range.
  # ------------------------------------------------------------------------
  occupied_snps <- copy(occupied_snps)
  setorder(occupied_snps, SNP_pos, rs)

  if (nrow(occupied_snps) > 0L) {

    occupied_snps[, SNP_plot_x := SNP_pos]

    if (SPREAD_OCCUPIED_SNPS && nrow(occupied_snps) > 1L) {

      locus_span <- axis_end - axis_start
      min_gap <- OCCUPIED_SNP_MIN_GAP_FRAC * locus_span
      display_left <- axis_start + OCCUPIED_SNP_BOUNDARY_FRAC * locus_span
      display_right <- axis_end - OCCUPIED_SNP_BOUNDARY_FRAC * locus_span

      # Forward pass: resolve local overlaps while preserving order.
      for (ii in 2:nrow(occupied_snps)) {
        occupied_snps$SNP_plot_x[ii] <- max(
          occupied_snps$SNP_pos[ii],
          occupied_snps$SNP_plot_x[ii - 1L] + min_gap
        )
      }

      # Shift back if the right-most displayed SNP exceeds the plotting range.
      right_excess <- max(occupied_snps$SNP_plot_x) - display_right
      if (is.finite(right_excess) && right_excess > 0) {
        occupied_snps[, SNP_plot_x := SNP_plot_x - right_excess]
      }

      # Backward pass: maintain the minimum gap after the right-boundary shift.
      if (nrow(occupied_snps) >= 2L) {
        for (ii in (nrow(occupied_snps) - 1L):1L) {
          occupied_snps$SNP_plot_x[ii] <- min(
            occupied_snps$SNP_plot_x[ii],
            occupied_snps$SNP_plot_x[ii + 1L] - min_gap
          )
        }
      }

      # If the left-most SNP is pushed too far left, shift the full set right.
      left_excess <- display_left - min(occupied_snps$SNP_plot_x)
      if (is.finite(left_excess) && left_excess > 0) {
        occupied_snps[, SNP_plot_x := SNP_plot_x + left_excess]
      }

      # Keep isolated SNPs as close as possible to their true locations.
      # The above passes only move SNPs when needed to prevent overlap.
    }

  } else {
    occupied_snps[, SNP_plot_x := numeric()]
  }

  snps <- merge(
    snps,
    occupied_snps[, .(rs, SNP_plot_x)],
    by = "rs",
    all.x = TRUE,
    sort = FALSE
  )
  snps[is.na(SNP_plot_x), SNP_plot_x := SNP_pos]

  if (nrow(tf_long) > 0L) {
    tf_long <- merge(
      tf_long,
      occupied_snps[, .(rs, SNP_plot_x)],
      by = "rs",
      all.x = TRUE,
      sort = FALSE
    )
    tf_long[is.na(SNP_plot_x), SNP_plot_x := SNP_pos]
  } else {
    tf_long[, SNP_plot_x := numeric()]
  }

  # Gene label placement
  gene_mid <- (gene_start + gene_end) / 2
  gene_width <- gene_end - gene_start
  gene_fraction <- gene_width / locus_width

  if (is.finite(gene_fraction) && gene_fraction >= 0.08) {
    gene_label_y <- 0
  } else {
    gene_label_y <- -0.46
  }

  # ------------------------------------------------------------------------
  # rs labels: place labels directly above the occupied SNP TF stacks.
  # Vertical text (90°) keeps labels compact and avoids the bottom clutter.
  # ------------------------------------------------------------------------
  setorder(snps, SNP_pos, rs)

  snp_label_df <- data.table(
    rs = character(),
    SNP_pos = numeric(),
    SNP_plot_x = numeric(),
    SNP_label_x = numeric(),
    SNP_label_y = numeric()
  )

  if (SHOW_SNP_LABELS && nrow(occupied_snps) > 0L) {

    snp_label_df <- copy(occupied_snps)

    if (nrow(tf_long) > 0L) {
      stack_top_df <- tf_long[
        ,
        .(stack_top_y = max(TF_y, na.rm = TRUE)),
        by = .(rs, SNP_pos, SNP_plot_x)
      ]

      snp_label_df <- merge(
        snp_label_df,
        stack_top_df,
        by = c("rs", "SNP_pos", "SNP_plot_x"),
        all.x = TRUE,
        sort = FALSE
      )
    } else {
      snp_label_df[, stack_top_y := TF_BASE_Y]
    }

    snp_label_df[is.na(stack_top_y), stack_top_y := TF_BASE_Y]

    setorder(snp_label_df, SNP_plot_x, rs)

    snp_label_df[, SNP_label_x := SNP_plot_x]
    snp_label_df[, SNP_label_y := stack_top_y + SNP_LABEL_TOP_PADDING]

    # Apply a very small vertical stagger for immediately neighboring stacks
    # so that top labels remain clear when two displayed SNP stacks are close.
    local_gap_for_stagger <- 0.65 * OCCUPIED_SNP_MIN_GAP_FRAC * (axis_end - axis_start)
    if (!is.finite(local_gap_for_stagger) || local_gap_for_stagger <= 0) {
      local_gap_for_stagger <- Inf
    }

    snp_label_df[
      ,
      close_pair := c(FALSE, diff(SNP_plot_x) < local_gap_for_stagger)
    ]

    if (nrow(snp_label_df) >= 2L) {
      snp_label_df[
        close_pair == TRUE,
        SNP_label_y := SNP_label_y + 0.12
      ]
    }

  }

  snp_label_angle <- SNP_LABEL_ANGLE

  # ------------------------------------------------------------------------
  # Y range and output dimensions
  # ------------------------------------------------------------------------
  base_y_max <- max(
    1.30,
    TF_BASE_Y + max(1L, max_tf_stack - 1L) * TF_VERTICAL_SPACING + 0.70
  )

  if (exists("snp_label_df") && nrow(snp_label_df) > 0L) {
    base_y_max <- max(
      base_y_max,
      max(snp_label_df$SNP_label_y, na.rm = TRUE) + 0.35
    )
  }

  # Custom legend showing all occupied-SNP TF names grouped by category,
  # with each TF displayed using its plot symbol.
  legend_category_df <- data.table(
    TF_category = factor(),
    category_x = numeric(),
    category_y = numeric(),
    category_label = character()
  )

  legend_item_df <- data.table(
    TF_category = factor(),
    TF = character(),
    TF_shape = numeric(),
    point_x = numeric(),
    point_y = numeric(),
    text_x = numeric(),
    text_y = numeric()
  )

  if (SHOW_CUSTOM_TF_TEXT_LEGEND && nrow(tf_long) > 0L) {

    legend_item_df <- unique(
      tf_long[
        !is.na(TF_category) & !is.na(TF) & TF != "",
        .(TF_category, TF, TF_shape)
      ]
    )

    setorder(legend_item_df, TF_category, TF)

    displayed_categories <- unique(legend_item_df$TF_category)
    legend_n_rows <- length(displayed_categories)

    legend_extra_height <- CUSTOM_TF_LEGEND_Y_PADDING +
      (legend_n_rows + 1L) * CUSTOM_TF_LEGEND_LINE_SPACING + 0.35

  } else {
    displayed_categories <- factor()
    legend_n_rows <- 0L
    legend_extra_height <- 0
  }

  y_max <- base_y_max + legend_extra_height

  if (legend_n_rows > 0L) {
    legend_title_y <- y_max - CUSTOM_TF_LEGEND_Y_PADDING

    legend_left_x <- axis_start +
      CUSTOM_TF_LEGEND_LEFT_FRAC * (axis_end - axis_start)
    legend_right_x <- axis_start +
      CUSTOM_TF_LEGEND_RIGHT_FRAC * (axis_end - axis_start)

    category_block_width <-
      CUSTOM_TF_LEGEND_CATEGORY_BLOCK_FRAC * (axis_end - axis_start)
    est_char_width <-
      CUSTOM_TF_LEGEND_CHAR_WIDTH_FRAC * (axis_end - axis_start)
    point_text_gap <-
      CUSTOM_TF_LEGEND_POINT_TEXT_GAP_FRAC * (axis_end - axis_start)
    item_gap <-
      CUSTOM_TF_LEGEND_ITEM_GAP_FRAC * (axis_end - axis_start)

    legend_category_df <- data.table(
      TF_category = displayed_categories,
      category_x = legend_left_x,
      category_y = legend_title_y -
        seq_len(legend_n_rows) * CUSTOM_TF_LEGEND_LINE_SPACING,
      category_label = paste0(as.character(displayed_categories), ":")
    )

    legend_item_df <- merge(
      legend_item_df,
      legend_category_df[, .(TF_category, category_y)],
      by = "TF_category",
      all.x = TRUE,
      sort = FALSE
    )

    legend_item_df[
      ,
      `:=`(
        point_y = category_y,
        text_y = category_y,
        TF_nchar = nchar(TF)
      )
    ]

    legend_item_df[
      ,
      est_item_width := point_text_gap + TF_nchar * est_char_width + item_gap
    ]

    legend_item_df[
      ,
      point_x := legend_left_x + category_block_width +
        cumsum(shift(est_item_width, fill = 0)),
      by = TF_category
    ]

    legend_item_df[, text_x := point_x + point_text_gap]

  } else {
    legend_title_y <- NA_real_
  }

  y_min <- -0.98

  n_snps <- nrow(snps)
  n_unique_tfs <- uniqueN(tf_long$TF)
  n_unique_tf_categories <- uniqueN(
    tf_long[!is.na(TF_category), TF_category]
  )

  plot_width <- PLOT_WIDTH
  plot_height <- PLOT_HEIGHT

  # ------------------------------------------------------------------------
  # Build figure
  # ------------------------------------------------------------------------
  p <- ggplot() +

    # Genomic baseline with arrow
    geom_segment(
      aes(
        x = axis_start,
        xend = axis_end,
        y = 0,
        yend = 0
      ),
      linewidth = 0.65,
      arrow = grid::arrow(
        type = "closed",
        length = grid::unit(0.12, "inches")
      )
    ) +

    # Gene body rectangle
    annotate(
      "rect",
      xmin = gene_start,
      xmax = gene_end,
      ymin = -0.14,
      ymax = 0.14,
      fill = "grey92",
      colour = "steelblue4",
      linewidth = 0.9
    ) +

    # Gene symbol
    annotate(
      "text",
      x = gene_mid,
      y = gene_label_y,
      label = gene_symbol,
      fontface = "italic",
      size = GENE_LABEL_SIZE
    ) +

    # Start / end coordinate labels
    annotate(
      "text",
      x = gene_start,
      y = GENE_COORD_Y,
      label = paste0("start: ", comma(gene_start, accuracy = 1)),
      hjust = 1,
      vjust = 1,
      size = GENE_COORD_LABEL_SIZE
    ) +

    annotate(
      "text",
      x = gene_end,
      y = GENE_COORD_Y,
      label = paste0("end: ", comma(gene_end, accuracy = 1)),
      hjust = 0,
      vjust = 1,
      size = GENE_COORD_LABEL_SIZE
    ) +

    # Guide lines from gene ends to baseline
    geom_segment(
      data = data.frame(
        x = c(gene_start, gene_end),
        xend = c(gene_start, gene_end),
        y = c(-0.14, -0.14),
        yend = c(0, 0)
      ),
      aes(x = x, xend = xend, y = y, yend = yend),
      linewidth = 0.35,
      linetype = "dashed"
    )

  # TSS
  if (SHOW_TSS && !is.na(tss_pos)) {
    p <- p +
      geom_segment(
        data = data.frame(
          x = tss_pos,
          xend = tss_pos,
          y = -0.14,
          yend = TSS_LABEL_Y + 0.10
        ),
        aes(x = x, xend = xend, y = y, yend = yend),
        linewidth = 0.45,
        linetype = "dotted"
      ) +
      annotate(
        "text",
        x = tss_pos,
        y = TSS_LABEL_Y,
        label = paste0(
          "TSS (", gene_direction, "): ",
          comma(tss_pos, accuracy = 1)
        ),
        hjust = 0.5,
        vjust = 1,
        size = TSS_LABEL_SIZE
      )
  }

  p <- p +

    # All model SNPs at their true genomic positions
    geom_point(
      data = snps,
      aes(
        x = SNP_pos,
        y = 0.22
      ),
      size = SNP_POINT_SIZE,
      shape = 16
    ) +

    # Connect true genomic positions to the spread-out displayed occupied SNPs
    geom_segment(
      data = occupied_snps,
      aes(
        x = SNP_pos,
        xend = SNP_plot_x,
        y = 0.22,
        yend = 0.70
      ),
      linewidth = 0.28,
      colour = "grey45"
    ) +

    # Stems for the displayed occupied SNP stacks
    geom_segment(
      data = occupied_snps,
      aes(
        x = SNP_plot_x,
        xend = SNP_plot_x,
        y = 0.70,
        yend = TF_BASE_Y - 0.08
      ),
      linewidth = 0.40
    )

  # TF occupancy symbols:
  #   colour = functional program
  #   shape  = individual TF within that program
  if (nrow(tf_long) > 0L) {
    p <- p +
      geom_point(
        data = tf_long,
        aes(
          x = SNP_plot_x,
          y = TF_y,
          colour = TF_category,
          shape = TF
        ),
        size = TF_POINT_SIZE,
        stroke = 0.6
      )
  }

  # rsID labels above occupied SNP TF stacks
  if (SHOW_SNP_LABELS && exists("snp_label_df") && nrow(snp_label_df) > 0L) {

    if (SHOW_SNP_LABEL_CONNECTORS) {
      p <- p +
        geom_segment(
          data = snp_label_df,
          aes(
            x = SNP_plot_x,
            xend = SNP_label_x,
            y = stack_top_y + 0.03,
            yend = SNP_label_y - 0.03
          ),
          linewidth = 0.28
        )
    }

    p <- p +
      geom_text(
        data = snp_label_df,
        aes(
          x = SNP_label_x,
          y = SNP_label_y,
          label = rs
        ),
        angle = snp_label_angle,
        hjust = 0,
        vjust = 0.5,
        size = SNP_LABEL_SIZE
      )
  }

  subtitle_text <- paste0(
    target_celltype,
    " | ",
    gene_chr,
    " | ",
    gene_id,
    " | gene: ",
    comma(gene_start, accuracy = 1),
    "-",
    comma(gene_end, accuracy = 1),
    ifelse(
      is.na(gene_direction),
      "",
      paste0(" | strand: ", gene_direction)
    ),
    ifelse(
      is.na(tss_pos),
      "",
      paste0(" | TSS: ", comma(tss_pos, accuracy = 1))
    ),
    " | ",
    n_snps,
    " model SNPs shown | ",
    sum(snps$Cistrome_overlap, na.rm = TRUE),
    " SNPs with Cistrome overlap | ",
    n_unique_tfs,
    " unique TFs displayed across ",
    n_unique_tf_categories,
    " functional programs"
  )

  p <- p +
    scale_x_continuous(
      labels = label_number(big.mark = ",", accuracy = 1),
      breaks = pretty_breaks(n = 8),
      expand = expansion(mult = c(0.01, 0.01))
    ) +
    scale_colour_manual(
      values = TF_CATEGORY_COLORS,
      drop = TRUE,
      guide = "none"
    ) +
    scale_shape_manual(
      values = TF_SHAPE_VALUES,
      drop = TRUE,
      guide = "none"
    ) +
    labs(
      title = paste0(
        gene_symbol,
        ": model variants and TF occupancy"
      ),
      subtitle = subtitle_text,
      x = paste0(gene_chr, " genomic position (hg38)"),
      y = NULL,
      colour = "TF functional program",
      shape = "TF"
    ) +
    coord_cartesian(
      xlim = c(axis_start, axis_end),
      ylim = c(y_min, y_max),
      clip = "off"
    ) +
    theme_classic(base_size = 16) +
    theme(
      axis.line.x = element_blank(),
      axis.line.y = element_blank(),
      axis.text.x = element_text(size = 12),
      axis.text.y = element_blank(),
      axis.title.x = element_text(size = 14),
      axis.title.y = element_blank(),
      axis.ticks.y = element_blank(),
      plot.title = element_text(face = "bold", size = 19),
      plot.subtitle = element_text(size = 13),
      legend.position = "none",
      plot.margin = margin(t = 18, r = 34, b = 28, l = 14)
    )

  if (nrow(legend_category_df) > 0L) {
    p <- p +
      annotate(
        "text",
        x = axis_start + CUSTOM_TF_LEGEND_LEFT_FRAC * (axis_end - axis_start),
        y = legend_title_y,
        label = "TFs occupying SNPs (grouped by functional program)",
        hjust = 0,
        vjust = 0.5,
        size = CUSTOM_TF_LEGEND_TITLE_SIZE,
        fontface = "bold"
      ) +
      geom_text(
        data = legend_category_df,
        aes(
          x = category_x,
          y = category_y,
          label = category_label,
          colour = TF_category
        ),
        hjust = 0,
        vjust = 0.5,
        size = CUSTOM_TF_LEGEND_CATEGORY_LABEL_SIZE,
        fontface = "bold",
        show.legend = FALSE
      ) +
      geom_point(
        data = legend_item_df,
        aes(
          x = point_x,
          y = point_y,
          colour = TF_category,
          shape = TF
        ),
        size = CUSTOM_TF_LEGEND_POINT_SIZE,
        stroke = 0.6,
        show.legend = FALSE
      ) +
      geom_text(
        data = legend_item_df,
        aes(
          x = text_x,
          y = text_y,
          label = TF,
          colour = TF_category
        ),
        hjust = 0,
        vjust = 0.5,
        size = CUSTOM_TF_LEGEND_TEXT_SIZE,
        show.legend = FALSE
      )
  }

  list(
    plot = p,
    snps = snps,
    tf_long = tf_long,
    gene_info = data.table(
      CellType = target_celltype,
      gene = gene_id,
      genename = gene_symbol,
      gene_chr = gene_chr,
      gene_start = gene_start,
      gene_end = gene_end,
      gene_direction = gene_direction,
      gene_TSS = tss_pos,
      N_model_SNPs_plotted = n_snps,
      N_TF_occupied_SNPs = sum(snps$Cistrome_overlap, na.rm = TRUE),
      N_unique_TFs_displayed = n_unique_tfs,
      N_TF_functional_programs_displayed = n_unique_tf_categories
    ),
    width = plot_width,
    height = plot_height
  )
}

# ============================================================================
# 6. Save outputs
# ============================================================================

save_gene_snp_tf <- function(
  result,
  out_dir,
  prefix
) {

  dir.create(
    out_dir,
    recursive = TRUE,
    showWarnings = FALSE
  )

  safe_prefix <- gsub(
    "[^A-Za-z0-9_.-]+",
    "_",
    prefix
  )

  if (SAVE_PDF) {
    ggsave(
      filename = file.path(
        out_dir,
        paste0(safe_prefix, "_Gene_SNP_TF_occupancy.pdf")
      ),
      plot = result$plot,
      width = result$width,
      height = result$height,
      units = "in",
      device = cairo_pdf,
      limitsize = FALSE
    )
  }

  if (SAVE_PNG) {
    ggsave(
      filename = file.path(
        out_dir,
        paste0(safe_prefix, "_Gene_SNP_TF_occupancy.png")
      ),
      plot = result$plot,
      width = result$width,
      height = result$height,
      units = "in",
      dpi = PNG_DPI,
      limitsize = FALSE
    )
  }

  fwrite(
    result$snps,
    file.path(
      out_dir,
      paste0(safe_prefix, "_plotted_SNPs.csv")
    )
  )

  fwrite(
    result$tf_long,
    file.path(
      out_dir,
      paste0(safe_prefix, "_plotted_SNP_TF_pairs.csv")
    )
  )

  fwrite(
    result$gene_info,
    file.path(
      out_dir,
      paste0(safe_prefix, "_plot_summary.csv")
    )
  )
}

# ============================================================================
# 7. Run
# ============================================================================

if (RUN_MODE == "single") {

  message(
    "Plotting ",
    TARGET_CELLTYPE,
    " / ",
    TARGET_GENE
  )

  result <- plot_gene_snp_tf(
    x = dat,
    target_gene = TARGET_GENE,
    target_celltype = TARGET_CELLTYPE,
    show_all_model_snps = SHOW_ALL_MODEL_SNPS,
    tf_keep = TF_KEEP
  )

  print(result$plot)

  gene_symbol_out <- result$gene_info$genename[1]

  single_out_dir <- file.path(
    OUT_DIR,
    TARGET_CELLTYPE
  )

  save_gene_snp_tf(
    result = result,
    out_dir = single_out_dir,
    prefix = paste0(
      TARGET_CELLTYPE,
      "_",
      gene_symbol_out
    )
  )

  message(
    "Finished: ",
    TARGET_CELLTYPE,
    " / ",
    gene_symbol_out
  )

} else if (RUN_MODE == "all") {

  plot_targets <- unique(
    dat[
      !is.na(CellType) &
        CellType != "" &
        !is.na(gene) &
        gene != "",
      .(
        CellType,
        gene,
        genename
      )
    ]
  )

  setorder(
    plot_targets,
    CellType,
    genename,
    gene
  )

  message(
    "Generating figures for ",
    nrow(plot_targets),
    " CellType-gene pairs."
  )

  for (i in seq_len(nrow(plot_targets))) {

    ct_i <- plot_targets$CellType[i]
    gene_i <- plot_targets$gene[i]
    gene_symbol_i <- plot_targets$genename[i]

    message(
      "[", i, "/", nrow(plot_targets), "] ",
      ct_i, " / ", gene_symbol_i, " (", gene_i, ")"
    )

    result_i <- tryCatch(
      plot_gene_snp_tf(
        x = dat,
        target_gene = gene_i,
        target_celltype = ct_i,
        show_all_model_snps = SHOW_ALL_MODEL_SNPS,
        tf_keep = TF_KEEP
      ),
      error = function(e) {
        warning(
          "Skipping ",
          ct_i,
          " / ",
          gene_i,
          ": ",
          conditionMessage(e)
        )
        NULL
      }
    )

    if (is.null(result_i)) {
      next
    }

    ct_out_dir <- file.path(
      OUT_DIR,
      ct_i
    )

    save_gene_snp_tf(
      result = result_i,
      out_dir = ct_out_dir,
      prefix = paste0(
        ct_i, "_", gene_symbol_i, "_", gene_i
      )
    )
  }

} else {

  stop(
    "Unknown RUN_MODE: ",
    RUN_MODE,
    ". Use 'single' or 'all'."
  )
}

message("\n============================================================")
message("Gene–variant–TF plotting finished.")
message("Output directory: ", OUT_DIR)
message("============================================================")
