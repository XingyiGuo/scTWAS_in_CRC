#!/usr/bin/env Rscript

# ==============================================================================
# Manhattan plots for HTAN_NL + WellcomeSanger_NL Cauchy meta-TWAS (EUR)
#
# Input files:
#   HTAN_NL_WellcomeSanger_NL_{ct}_peakSNP.EUR.CauchyMeta.csv
#
# Cell types:
#   ABS, GOB, CT, STM, EE, TUF, TAC
#
# Plotting/testing rules:
#   1. Use META_P as the EUR Cauchy meta-analysis P value.
#   2. Bonferroni denominator = number of eligible meta-TWAS tests.
#      Eligibility is based on source-specific R2_global > 0.01:
#        - HTAN_NL_only:                  R2_global_HTAN_NL > 0.01
#        - WellcomeSanger_NL_only:         R2_global_WellcomeSanger_NL > 0.01
#        - HTAN_NL + WellcomeSanger_NL:    BOTH R2_global values > 0.01
#      For an unrecognized META_SOURCE, require at least one R2_global > 0.01.
#   3. Bonferroni cutoff = 0.05 / number of eligible tests, separately by cell type.
#   4. Plot only eligible tests.
#   5. Non-Bonferroni genes:
#        odd chromosomes  = blue   (#2B57C4)
#        even chromosomes = orange (#F29400)
#   6. Bonferroni-significant genes are drawn as cell-type-colored circles with
#      a black rim, and their gene symbols are labeled in regular-weight black
#      Calibri. Panel titles are black.
#   7. DISPLAY caps are cell-type-specific:
#        - ABS, GOB, CT, STM: META_P < 1e-40 is plotted at 1e-40 (-log10 P = 40).
#        - EE, TUF, TAC:      META_P < 1e-20 is plotted at 1e-20 (-log10 P = 20).
#      Significance always uses the original uncapped META_P.
#   8. Panel-height ratio is 2:1:
#        - ABS, GOB, CT, STM have height 2.
#        - EE, TUF, TAC have height 1.
#   9. Output layout: one combined PNG and one combined PDF for all 7 cell types.
# ==============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(ggplot2)
  library(ggrepel)
  library(patchwork)
  library(ensembldb)
  library(EnsDb.Hsapiens.v86)
})

# ------------------------------------------------------------------------------
# User settings
# ------------------------------------------------------------------------------

ASSOC_ROOT <- "/v5000/data/l2_bioinfo1/liq17/scTWAS/TWAS_HTAN_GSE266546_WellcomeSanger_res/SPrediXcan_peakSNP_association/Cauchy_Meta_NL_WellcomeSanger_NL_EUR/"

INPUT_PATTERN <- file.path(
  ASSOC_ROOT,
  "HTAN_NL_WellcomeSanger_NL_%s_peakSNP.EUR.CauchyMeta.csv"
)

CTS <- c("ABS", "GOB", "CT", "STM", "EE", "TUF", "TAC")

OUTDIR <- file.path(
  ASSOC_ROOT,
  "Manhattan_HTAN_NL_WellcomeSanger_NL_EUR_CauchyMeta_7CT_R2gt0.01"
)
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

ALPHA <- 0.05
R2_THRESHOLD <- 0.01

# Cell-type-specific display caps.
# The Manhattan y-axis is -log10(P), so these correspond to maxima of 40 and 20.
PLOT_P_FLOOR_BY_CT <- c(
  "ABS" = 1e-40,
  "GOB" = 1e-40,
  "CT"  = 1e-40,
  "STM" = 1e-40,
  "EE"  = 1e-20,
  "TUF" = 1e-20,
  "TAC" = 1e-20
)

Y_MAX_BY_CT <- -log10(PLOT_P_FLOOR_BY_CT)

# Make the physical y-axis/panel height for the 40-unit panels twice that
# of the 20-unit panels, preserving the same vertical scale per -log10(P) unit.
PANEL_HEIGHTS <- c(
  "ABS" = 2,
  "GOB" = 2,
  "CT"  = 2,
  "STM" = 2,
  "EE"  = 1,
  "TUF" = 1,
  "TAC" = 1
)

stopifnot(
  all(CTS %in% names(PLOT_P_FLOOR_BY_CT)),
  all(CTS %in% names(Y_MAX_BY_CT)),
  all(CTS %in% names(PANEL_HEIGHTS))
)

# A4 output sizes (inches)
A4_WIDTH <- 9
A4_HEIGHT <- 15
A4_DPI <- 300

LABEL_META_SIGNIFICANT <- TRUE
POINT_SIZE <- 0.55
POINT_ALPHA <- 0.85
SIG_POINT_SIZE <- 1.7
SIG_POINT_STROKE <- 0.35

# Odd/even chromosome colors for non-significant genes.
# Requested 2026-09-28: blue for odd chromosomes, orange for even chromosomes
# (the alternating scheme used by the published pancreatic-cancer TWAS Manhattan).
# Both clear 3:1 contrast against white and are distinguishable in greyscale
# (blue prints darker than orange).
CHR_COLORS <- c(
  "odd"  = "#2B57C4",
  "even" = "#F29400"
)

# ------------------------------------------------------------------------------
# Typography and gene-label style
# ------------------------------------------------------------------------------
# All text in the figure uses ONE family. Naming a family that the graphics
# device cannot actually draw makes EVERY glyph come out as fragments, which is
# how the 2026-09-28 cluster run lost its gene symbols, axis numbers and panel
# titles while the points drew normally. The family is therefore VERIFIED here
# before it is used:
#   1. fontconfig must resolve the requested name to that SAME family
#      (it silently substitutes a different font otherwise),
#   2. the resolved file must exist and be a scalable outline font, and
#   3. it must really contain Latin glyphs (glyph index 0 = missing glyph).
# If no preference passes, PLOT_FONT stays "" = the device's built-in sans,
# which is what these scripts used before any family was requested, so text
# always renders even when no named font is available.
#
# PLOT_FONT_OVERRIDE short-circuits the check:
#   PLOT_FONT_OVERRIDE <- ""        # force the device default (always safe)
#   PLOT_FONT_OVERRIDE <- "Arial"   # force a family, no verification
PLOT_FONT_OVERRIDE <- NA_character_
PLOT_FONT_PREFS <- c("Calibri", "Arial", "Helvetica", "Liberation Sans", "DejaVu Sans")

font_is_usable <- function(fam) {
  if (!requireNamespace("systemfonts", quietly = TRUE)) return(FALSE)
  info <- try(systemfonts::font_info(family = fam, size = 10), silent = TRUE)
  if (inherits(info, "try-error")) return(FALSE)
  if (!isTRUE(tolower(as.character(info$family[1])) == tolower(fam))) return(FALSE)
  if (!isTRUE(file.exists(as.character(info$path[1])))) return(FALSE)
  if (!isTRUE(as.logical(info$scalable[1]))) return(FALSE)
  gl <- try(
    systemfonts::glyph_info(c("A", "g", "1"), family = fam, size = 10),
    silent = TRUE
  )
  if (inherits(gl, "try-error")) return(FALSE)
  isTRUE(all(!is.na(gl$index)) && all(gl$index > 0))
}

resolve_plot_font <- function(prefs = PLOT_FONT_PREFS) {
  if (!is.na(PLOT_FONT_OVERRIDE)) {
    message(
      "Font: PLOT_FONT_OVERRIDE is set to '", PLOT_FONT_OVERRIDE,
      "'; using it without verification."
    )
    return(PLOT_FONT_OVERRIDE)
  }
  for (fam in prefs) {
    if (font_is_usable(fam)) {
      if (.Platform$OS.type == "windows") {
        try(
          do.call(
            grDevices::windowsFonts,
            stats::setNames(list(grDevices::windowsFont(fam)), fam)
          ),
          silent = TRUE
        )
      }
      info <- systemfonts::font_info(family = fam, size = 10)
      message("Font: using ", fam, " (", basename(as.character(info$path[1])), ").")
      return(fam)
    }
    message("Font: ", fam, " is not usable here; trying the next preference.")
  }
  message(
    "Font: none of ", paste(prefs, collapse = ", "),
    " could be verified; falling back to the device default sans family. ",
    "Text will be legible but it will not be Arial."
  )
  ""
}

PLOT_FONT <- resolve_plot_font()

# Bonferroni-significant gene symbols are drawn in BLACK rather than in the
# cell-type color, so they read as annotation instead of as data, at REGULAR
# weight (requested 2026-09-28), with a thin white halo (ggrepel >= 0.9.0) so
# they stay legible where a label overlaps the point cloud. Significant POINTS
# keep the cell-type color, now as a filled circle with a black rim, which is
# what makes them obvious against the blue/orange background.
#
# On "Calibri Light": it is NOT a family of its own, it is the light WEIGHT of
# the Calibri family (systemfonts reports family = "Calibri",
# style = "Calibri Light"). Asking a graphics device for family
# "Calibri Light" therefore does not give you that face - measured on the
# Windows workstation, the cairo device renders family "Calibri Light"
# pixel-identically to Arial, i.e. it silently substitutes. Addressing the
# light face needs systemfonts::register_variant() plus a systemfonts-backed
# device (ragg for PNG); grDevices::cairo_pdf does not consult systemfonts, so
# a registered alias would leave the PDF and the PNG in different fonts. The
# labels are consequently set in the verified base family (Calibri where
# available) at regular weight, which is the closest faithful rendering that
# both output devices agree on.
GENE_LABEL_COLOR <- "black"
GENE_LABEL_FACE <- "plain"
GENE_LABEL_SIZE <- 2.6
GENE_LABEL_HALO <- "white"
GENE_LABEL_HALO_R <- 0.12
GENE_LABEL_SEGMENT <- "grey30"

# Cell-type-specific colors for Bonferroni-significant genes.
# Cell-type highlight colors, shared verbatim by all four Manhattan scripts.
# Nature-journal (NPG) palette, darkened so every color clears a 4.5:1 contrast
# ratio against white (WCAG AA for text) -- point colors, threshold lines and
# ggrepel gene symbols all use these, so gene labels stay legible in print.
# Contrast vs white 4.5-9.2:1; minimum CIELAB dE between any two co-occurring
# cell types is 30 (NL 7CT) / 30 (ADSER 6CT).
CT_SIG_COLORS <- c(
  "ABS" = "#B5382A",  # red        (Nature/NPG #E64B35, darkened)
  "GOB" = "#17557F",  # blue
  "CT"  = "#A76A00",  # amber
  "EE"  = "#6E4B3A",  # brown      (Nature/NPG #7E6148, darkened)
  "STM" = "#4433A8",  # indigo
  "TUF" = "#6E6E22",  # olive
  "TAC" = "#17695E",  # teal       (Nature/NPG #00A087 / #4DBBD5 family, darkened)
  "ASC" = "#8A5015",  # ochre
  "SSC" = "#882255"   # wine
)

# Full cell-type display names. These are the ONLY text drawn inside a panel.
CT_FULL_NAMES <- c(
  "ABS" = "Absorptive Cell Type (ABS)",
  "CT"  = "Crypt-Top Cell Type (CT)",
  "EE"  = "Enteroendocrine Cell Type (EE)",
  "GOB" = "Goblet Cell Type (GOB)",
  "STM" = "Stem Cell Type (STM)",
  "TAC" = "Transit-Amplifying Cell Type (TAC)",
  "TUF" = "Tuft Cell Type (TUF)",
  "ASC" = "Adenoma-Specific Cell Type (ASC)",
  "SSC" = "Serrated-Polyp-Specific Cell Type (SSC)"
)

# Size (pt) and vertical placement of the in-panel name.
# CT_NAME_Y_FRAC is a fraction of the y range; 0.97 with vjust = 1 puts the
# name top-centre inside the panel. Lower it if a tall peak sits mid-panel.
CT_NAME_SIZE <- 3.6
# Panel title color. "black" (requested 2026-09-28) keeps the cell-type coding
# in the points and the threshold line only; set CT_NAME_COLOR <- NA to go back
# to drawing the name in the cell-type color.
CT_NAME_COLOR <- "black"
# 1.035 puts the name in the 7% head-room added above the data range by the
# y-scale expansion, so it can no longer collide with a gene label (gene labels
# are separately constrained to stay at or below the top of the data range).
# Values below 1 move it back inside the panel.
CT_NAME_Y_FRAC <- 1.035

stopifnot(all(CTS %in% names(CT_FULL_NAMES)))

stopifnot(all(CTS %in% names(CT_SIG_COLORS)))

# ------------------------------------------------------------------------------
# GRCh38 autosome lengths
# ------------------------------------------------------------------------------

CHR_LENGTHS <- tibble::tibble(
  chrom = 1:22,
  chr_len = c(
    248956422, 242193529, 198295559, 190214555, 181538259,
    170805979, 159345973, 145138636, 138394717, 133797422,
    135086622, 133275309, 114364328, 107043718, 101991189,
     90338345,  83257441,  80373285,  58617616,  64444167,
     46709983,  50818468
  )
) %>%
  dplyr::mutate(
    chr_offset = cumsum(chr_len) - chr_len,
    chr_center = chr_offset + chr_len / 2
  )

X_MAX <- sum(CHR_LENGTHS$chr_len)

# ------------------------------------------------------------------------------
# Gene annotation: Ensembl gene ID -> symbol, chromosome, midpoint (GRCh38)
# Explicit dplyr:: calls avoid namespace conflicts with Bioconductor packages.
# ------------------------------------------------------------------------------

message("Loading GRCh38 gene coordinates from EnsDb.Hsapiens.v86 ...")

GENE_ANNOT <- ensembldb::genes(
  EnsDb.Hsapiens.v86,
  return.type = "data.frame"
) %>%
  dplyr::transmute(
    gene = sub("\\..*$", "", as.character(gene_id)),
    annot_gene_name = as.character(gene_name),
    chrom = suppressWarnings(as.integer(as.character(seq_name))),
    gene_start = as.numeric(gene_seq_start),
    gene_end = as.numeric(gene_seq_end),
    bp = floor((gene_start + gene_end) / 2)
  ) %>%
  dplyr::filter(chrom %in% 1:22, !is.na(bp)) %>%
  dplyr::distinct(gene, .keep_all = TRUE)

# ------------------------------------------------------------------------------
# Read one cell type
# ------------------------------------------------------------------------------

read_one_ct <- function(ct) {

  plot_p_floor <- unname(PLOT_P_FLOOR_BY_CT[ct])
  y_max <- unname(Y_MAX_BY_CT[ct])

  infile <- sprintf(INPUT_PATTERN, ct)

  if (!file.exists(infile)) {
    stop(ct, ": input file not found:\n  ", infile)
  }

  message("Reading ", ct, ": ", infile)
  df <- data.table::fread(infile, data.table = FALSE)

  required <- c(
    "gene",
    "gene_name_HTAN_NL",
    "gene_name_WellcomeSanger_NL",
    "R2_global_HTAN_NL",
    "R2_global_WellcomeSanger_NL",
    "META_P",
    "META_SOURCE"
  )

  missing_cols <- setdiff(required, names(df))
  if (length(missing_cols) > 0) {
    stop(
      ct, ": missing required column(s): ",
      paste(missing_cols, collapse = ", "),
      "\nFile: ", infile
    )
  }

  # --------------------------------------------------------------------------
  # Define the Cauchy meta-TWAS testing universe.
  # --------------------------------------------------------------------------

  meta_p_all <- suppressWarnings(as.numeric(df$META_P))
  r2_htan_all <- suppressWarnings(as.numeric(df$R2_global_HTAN_NL))
  r2_ws_all <- suppressWarnings(as.numeric(df$R2_global_WellcomeSanger_NL))
  source_all <- as.character(df$META_SOURCE)

  valid_meta_p <-
    !is.na(meta_p_all) &
    is.finite(meta_p_all) &
    meta_p_all >= 0 &
    meta_p_all <= 1

  htan_r2_pass <-
    !is.na(r2_htan_all) &
    is.finite(r2_htan_all) &
    r2_htan_all > R2_THRESHOLD

  ws_r2_pass <-
    !is.na(r2_ws_all) &
    is.finite(r2_ws_all) &
    r2_ws_all > R2_THRESHOLD

  # Source-specific R2 rule.
  r2_eligible <- dplyr::case_when(
    source_all == "VUMC_HTAN_GSE266546_NL_only" ~ htan_r2_pass,
    source_all == "WellcomeSanger_NL_only" ~ ws_r2_pass,
    source_all == "VUMC_HTAN_GSE266546_NL+WellcomeSanger_NL" ~
      (htan_r2_pass & ws_r2_pass),
    TRUE ~ (htan_r2_pass | ws_r2_pass)
  )

  eligible_test <- valid_meta_p & r2_eligible

  n_total_rows <- nrow(df)
  n_valid_meta_p <- sum(valid_meta_p)
  n_tests <- sum(eligible_test)

  if (n_tests == 0) {
    stop(ct, ": no eligible Cauchy meta-TWAS tests after R2_global > ", R2_THRESHOLD)
  }

  bonf_cutoff <- ALPHA / n_tests

  # Diagnostic: eligible rows vs unique genes.
  eligible_gene_ids <- sub("\\..*$", "", as.character(df$gene[eligible_test]))
  n_unique_genes <- dplyr::n_distinct(eligible_gene_ids)

  if (n_unique_genes != n_tests) {
    warning(
      ct, ": ", n_tests, " eligible meta-TWAS rows but ", n_unique_genes,
      " unique gene IDs. Bonferroni denominator remains the number of eligible tests (rows)."
    )
  }

  # Source counts among eligible tests for QC.
  source_counts <- table(source_all[eligible_test], useNA = "ifany")
  source_msg <- paste(
    paste0(names(source_counts), "=", as.integer(source_counts)),
    collapse = "; "
  )

  out <- df %>%
    dplyr::transmute(
      CellType = ct,
      gene = sub("\\..*$", "", as.character(gene)),
      gene_name_HTAN_NL = as.character(gene_name_HTAN_NL),
      gene_name_WellcomeSanger_NL = as.character(gene_name_WellcomeSanger_NL),
      META_SOURCE = as.character(META_SOURCE),
      META_P_raw = suppressWarnings(as.numeric(META_P)),
      R2_global_HTAN_NL = suppressWarnings(as.numeric(R2_global_HTAN_NL)),
      R2_global_WellcomeSanger_NL = suppressWarnings(as.numeric(R2_global_WellcomeSanger_NL))
    ) %>%
    dplyr::mutate(
      HTAN_R2_pass =
        !is.na(R2_global_HTAN_NL) &
        is.finite(R2_global_HTAN_NL) &
        R2_global_HTAN_NL > R2_THRESHOLD,

      WellcomeSanger_R2_pass =
        !is.na(R2_global_WellcomeSanger_NL) &
        is.finite(R2_global_WellcomeSanger_NL) &
        R2_global_WellcomeSanger_NL > R2_THRESHOLD,

      R2_eligible = dplyr::case_when(
        META_SOURCE == "VUMC_HTAN_GSE266546_NL_only" ~ HTAN_R2_pass,
        META_SOURCE == "WellcomeSanger_NL_only" ~ WellcomeSanger_R2_pass,
        META_SOURCE == "VUMC_HTAN_GSE266546_NL+WellcomeSanger_NL" ~
          (HTAN_R2_pass & WellcomeSanger_R2_pass),
        TRUE ~ (HTAN_R2_pass | WellcomeSanger_R2_pass)
      )
    ) %>%
    dplyr::filter(
      !is.na(META_P_raw),
      is.finite(META_P_raw),
      META_P_raw >= 0,
      META_P_raw <= 1,
      R2_eligible
    ) %>%
    dplyr::left_join(GENE_ANNOT, by = "gene") %>%
    dplyr::left_join(
      CHR_LENGTHS %>% dplyr::select(chrom, chr_offset),
      by = "chrom"
    ) %>%
    dplyr::mutate(
      # Prefer HTAN symbol, then WellcomeSanger symbol, then EnsDb symbol, then ENSG.
      gene_name = dplyr::coalesce(
        dplyr::na_if(gene_name_HTAN_NL, ""),
        dplyr::na_if(gene_name_WellcomeSanger_NL, ""),
        dplyr::na_if(annot_gene_name, ""),
        gene
      ),

      N_total_rows = n_total_rows,
      N_valid_META_P = n_valid_meta_p,
      N_META_tests = n_tests,
      R2_threshold = R2_THRESHOLD,
      Bonferroni_cutoff = bonf_cutoff,

      # Significance uses the original uncapped meta P value.
      META_sig = META_P_raw < Bonferroni_cutoff,

      # Cell-type-specific DISPLAY cap only.
      # Significance continues to use META_P_raw above.
      Plot_P_floor = plot_p_floor,
      Y_max = y_max,
      META_P_plot = pmax(META_P_raw, Plot_P_floor),
      minus_log10P = -log10(META_P_plot),

      BPcum = bp + chr_offset,
      chr_parity = dplyr::if_else(chrom %% 2 == 0, "even", "odd")
    )

  n_missing_pos <- sum(is.na(out$chrom) | is.na(out$bp))
  if (n_missing_pos > 0) {
    warning(
      ct, ": ", n_missing_pos,
      " eligible meta-TWAS gene(s) lack an autosomal GRCh38 coordinate and will be omitted from the plot."
    )
  }

  out <- out %>%
    dplyr::filter(!is.na(chrom), !is.na(bp), !is.na(BPcum))

  message(
    "  Total meta rows = ", n_total_rows,
    " | eligible tests = ", n_tests,
    " | Bonferroni cutoff = ", format(bonf_cutoff, scientific = TRUE, digits = 4),
    " | plottable genes = ", nrow(out),
    " | Bonferroni-significant = ", sum(out$META_sig, na.rm = TRUE),
    " | display P floor = ", format(plot_p_floor, scientific = TRUE),
    " | y max = ", y_max,
    " | P values capped = ", sum(out$META_P_raw < plot_p_floor, na.rm = TRUE)
  )
  message("  Eligible META_SOURCE counts: ", source_msg)

  out
}

# ------------------------------------------------------------------------------
# Read all 7 normal-colon cell types
# ------------------------------------------------------------------------------

ALL <- lapply(CTS, read_one_ct) %>%
  dplyr::bind_rows()

# ------------------------------------------------------------------------------
# Summary table
# ------------------------------------------------------------------------------

SUMMARY <- ALL %>%
  dplyr::group_by(CellType) %>%
  dplyr::summarise(
    N_total_rows = dplyr::first(N_total_rows),
    N_valid_META_P = dplyr::first(N_valid_META_P),
    N_META_tests_R2gt0.01 = dplyr::first(N_META_tests),
    R2_threshold = dplyr::first(R2_threshold),
    Bonferroni_cutoff = dplyr::first(Bonferroni_cutoff),
    Plot_P_floor = dplyr::first(Plot_P_floor),
    Y_max = dplyr::first(Y_max),
    Panel_height_ratio = unname(PANEL_HEIGHTS[dplyr::first(CellType)]),
    N_rows_plotted = dplyr::n(),
    N_META_Bonferroni_sig = sum(META_sig, na.rm = TRUE),
    N_P_capped_at_display_floor = sum(META_P_raw < Plot_P_floor, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::mutate(CellType = factor(CellType, levels = CTS)) %>%
  dplyr::arrange(CellType)

data.table::fwrite(
  as.data.frame(SUMMARY),
  file.path(
    OUTDIR,
    "HTAN_NL_WellcomeSanger_NL_EUR_CauchyMeta_7CT_Manhattan_Bonferroni_summary.csv"
  )
)

print(SUMMARY)

# Save plotted gene-level data for QC/reproducibility.
data.table::fwrite(
  as.data.frame(ALL),
  file.path(
    OUTDIR,
    "HTAN_NL_WellcomeSanger_NL_EUR_CauchyMeta_7CT_Manhattan_plotted_genes.csv"
  )
)

# ------------------------------------------------------------------------------
# Plot one cell type
# ------------------------------------------------------------------------------

plot_one_ct <- function(ct) {

  df <- ALL %>%
    dplyr::filter(CellType == ct)

  if (nrow(df) == 0) {
    stop(ct, ": no plottable eligible genes.")
  }

  cutoff <- unique(df$Bonferroni_cutoff)
  n_tests <- unique(df$N_META_tests)
  y_max <- unique(df$Y_max)

  stopifnot(
    length(cutoff) == 1,
    length(n_tests) == 1,
    length(y_max) == 1
  )

  sig_color <- unname(CT_SIG_COLORS[ct])

  # Label each Bonferroni-significant gene once.
  label_df <- df %>%
    dplyr::filter(META_sig, !is.na(gene_name), gene_name != "") %>%
    dplyr::group_by(gene_name) %>%
    dplyr::slice_min(order_by = META_P_raw, n = 1, with_ties = FALSE) %>%
    dplyr::ungroup()

  p <- ggplot2::ggplot(
    df,
    ggplot2::aes(x = BPcum, y = minus_log10P)
  ) +

    # Non-Bonferroni genes: alternating chromosome colors (blue odd / orange even).
    ggplot2::geom_point(
      data = df %>% dplyr::filter(!META_sig),
      ggplot2::aes(color = chr_parity),
      size = POINT_SIZE,
      alpha = POINT_ALPHA
    ) +

    # Bonferroni-significant genes: cell-type-specific color.
    ggplot2::geom_point(
      data = df %>% dplyr::filter(META_sig),
      shape = 21,
      fill = sig_color,
      color = "black",
      stroke = SIG_POINT_STROKE,
      size = SIG_POINT_SIZE,
      alpha = 0.98
    ) +

    # Cell-type-specific Bonferroni line.
    ggplot2::geom_hline(
      yintercept = -log10(cutoff),
      linetype = "dashed",
      linewidth = 0.5,
      color = sig_color
    ) +

    ggplot2::scale_color_manual(
      values = CHR_COLORS,
      guide = "none"
    ) +

    ggplot2::scale_x_continuous(
      limits = c(0, X_MAX),
      breaks = CHR_LENGTHS$chr_center,
      labels = CHR_LENGTHS$chrom,
      expand = ggplot2::expansion(mult = c(0.002, 0.002))
    ) +

    # Cell-type-specific display range:
    # ABS/GOB/CT/STM: 0-40; EE/TUF/TAC: 0-20.
    ggplot2::scale_y_continuous(
      breaks = seq(0, y_max, by = 5),
      expand = ggplot2::expansion(mult = c(0, 0.07))
    ) +
    ggplot2::coord_cartesian(
      ylim = c(0, y_max),
      clip = "off"
    ) +

    # Cell-type name centred in the panel; no other text is drawn here.
    ggplot2::annotate(
      "text",
      x = X_MAX / 2,
      y = y_max * CT_NAME_Y_FRAC,
      label = unname(CT_FULL_NAMES[ct]),
      color = if (is.na(CT_NAME_COLOR)) sig_color else CT_NAME_COLOR,
      fontface = "bold",
      family = PLOT_FONT,
      size = CT_NAME_SIZE,
      hjust = 0.5,
      vjust = 0.5
    ) +

    ggplot2::labs(
      x = "Chromosome",
      y = expression(-log[10](italic(P)[Meta]))
    ) +

    ggplot2::theme_classic(base_size = 9, base_family = PLOT_FONT) +
    # Tightened format: smaller type, shorter ticks, axis titles and tick
    # labels pulled in towards the panel, and a narrower plot margin, so the
    # seven stacked panels carry less white space and more data area.
    ggplot2::theme(
      text = ggplot2::element_text(family = PLOT_FONT, color = "black"),
      legend.position = "none",
      axis.title = ggplot2::element_text(size = 8.5, color = "black"),
      axis.title.x = ggplot2::element_text(margin = ggplot2::margin(t = 1.5)),
      axis.title.y = ggplot2::element_text(margin = ggplot2::margin(r = 1.5)),
      axis.text.x = ggplot2::element_text(
        size = 6.2, color = "black", margin = ggplot2::margin(t = 1)
      ),
      axis.text.y = ggplot2::element_text(
        size = 6.8, color = "black", margin = ggplot2::margin(r = 1)
      ),
      axis.ticks = ggplot2::element_line(linewidth = 0.2, color = "black"),
      axis.ticks.length = grid::unit(1.2, "pt"),
      axis.line = ggplot2::element_line(linewidth = 0.3, color = "black"),
      plot.margin = ggplot2::margin(1, 4, 1, 2)
    )

  if (LABEL_META_SIGNIFICANT && nrow(label_df) > 0) {
    # Black, bold, haloed gene symbols. bg.color/bg.r need ggrepel >= 0.9.0,
    # so they are added conditionally and the layer still builds with older
    # versions installed on the cluster.
    repel_args <- list(
      data = label_df,
      mapping = ggplot2::aes(label = gene_name),
      color = GENE_LABEL_COLOR,
      fontface = GENE_LABEL_FACE,
      family = PLOT_FONT,
      size = GENE_LABEL_SIZE,
      max.overlaps = Inf,
      box.padding = 0.16,
      point.padding = 0.10,
      min.segment.length = 0,
      segment.color = GENE_LABEL_SEGMENT,
      segment.linewidth = 0.18,
      segment.alpha = 0.9,
      seed = 123,
      direction = "both",
      ylim = c(NA, y_max)
    )
    if (utils::packageVersion("ggrepel") >= "0.9.0") {
      repel_args$bg.color <- GENE_LABEL_HALO
      repel_args$bg.r <- GENE_LABEL_HALO_R
    }
    p <- p + do.call(ggrepel::geom_text_repel, repel_args)
  }

  p
}

# ------------------------------------------------------------------------------
# Generate plots
# ------------------------------------------------------------------------------

PLOTS <- lapply(CTS, plot_one_ct)
names(PLOTS) <- CTS

# No separate per-cell-type files are written here.
# The requested deliverable is one combined A4 PNG and one combined A4 PDF
# containing all 7 cell-type Manhattan plots.

# ------------------------------------------------------------------------------
# Combined outputs: one A4 PNG and one A4 PDF for all 7 cell types
# ------------------------------------------------------------------------------

COMBINED <- patchwork::wrap_plots(
  PLOTS,
  ncol = 1,
  heights = unname(PANEL_HEIGHTS[CTS])
) &
  ggplot2::theme(
    plot.margin = ggplot2::margin(1, 4, 1, 3)
  )

# Combined figure: all 7 cell types in one column.
# The first four panels use twice the height of the last three panels.
ggplot2::ggsave(
  filename = file.path(
    OUTDIR,
    "HTAN_NL_WellcomeSanger_NL_7CT_EUR_CauchyMeta_R2gt0.01_Manhattan_A4.png"
  ),
  plot = COMBINED,
  width = A4_WIDTH,
  height = A4_HEIGHT,
  dpi = A4_DPI,
  bg = "white"
)

ggplot2::ggsave(
  filename = file.path(
    OUTDIR,
    "HTAN_NL_WellcomeSanger_NL_7CT_EUR_CauchyMeta_R2gt0.01_Manhattan_A4.pdf"
  ),
  plot = COMBINED,
  width = A4_WIDTH,
  height = A4_HEIGHT,
  device = grDevices::cairo_pdf,
  bg = "white"
)

message("Done. Outputs written to: ", OUTDIR)
message("Generated:")
message("  - One combined A4 PNG for all 7 cell types")
message("  - One combined A4 PDF for all 7 cell types")
