#!/usr/bin/env Rscript
###############################################################################
# Fig_peaks_R2_violin.R
#
# PURPOSE
#   TWO figures, one per model group, each one row x two columns:
#
#     Figure 1  Normal colon   a) scATAC-seq peaks per cell type (bar)
#                              b) global R^2 of the scTWAS prediction models
#                                 (violin), TWO resources side by side
#     Figure 2  Colon polyps   a) scATAC-seq peaks per cell type (bar)
#                              b) global R^2 (violin), ONE resource
#
#   Panel b plots one violin per cell type per resource over the models that
#   pass R2_global > R2_MIN (0.01).
#
#   ON THE Y CAP.  R2_global has a long upper tail (single models reach ~0.8),
#   so a violin built on the full tail and then clipped at 0.10 shows only the
#   bottom eighth of its own shape - a sliver that reads as a column and says
#   nothing about the distribution.  The violins are therefore built on the
#   models INSIDE the displayed window (R2_MIN < R2_global <= Y_CAP), which is
#   what a reader of a capped axis interprets them as, and the models above the
#   cap are counted rather than hidden: the count and percentage per cell type
#   per resource go to the *_R2_sourcedata.csv, and the run prints the totals.
#   Set VIOLIN_WITHIN_CAP <- FALSE to keep the old clip-the-tail behaviour.
#
# INPUTS (all absolute; nothing is written to these locations)
#   scATAC peaks
#     F:/Dropbox/GuoLab/COREAD_ToDrGuo/scTWAS/ms-guo/SupplementaryTables.xlsx
#     sheet "scATAC-seq peaks" - Supplementary Table S5.  Its header is on the
#     SECOND row and the third row is a sub-header (Sample / Num scATAC peaks
#     under two study columns), so the reader below takes header = 2 and drops
#     rows with no CellType.  Peaks are taken from Num_regions_after_chain_merge
#     and the rows are selected with Group == "Normal colon" / "Colon polyps".
#
#   scTWAS model summaries (one CSV per cell type; column R2_global)
#     <TWAS_ROOT>/VUMC_NL/VUMC_HTAN_GSE266546_NL_<CT>_peakSNP_model_summaries.csv
#     <TWAS_ROOT>/WellcomeSanger_NL/WellcomeSanger_NL_<CT>_peakSNP_model_summaries.csv
#     <TWAS_ROOT>/VUMC_ADSER/VUMC_HTAN_ADSER_<CT>_peakSNP_model_summaries.csv
#   with TWAS_ROOT =
#     F:/Dropbox/GuoLab/COREAD_ToDrGuo/scTWAS/ms_figures_tables/res/TWAS/SPrediXcan_peakSNP_association
#
#   NOTE ON THE TWO NORMAL-COLON RESOURCES.  The folder names and the file
#   prefixes disagree with the way the request described them: VUMC_NL/ holds
#   VUMC_HTAN_GSE266546_NL_* files and WellcomeSanger_NL/ holds
#   WellcomeSanger_NL_* files.  Every resource label in these figures is read
#   from the FILE PREFIX, not from the folder name, so the legend cannot be
#   crossed by a folder that is named after the other cohort.
#
# OUTPUTS (OUT_DIR)
#   Fig_NormalColon_peaks_R2violin.pdf / .png
#   Fig_ColonPolyps_peaks_R2violin.pdf / .png
#   Fig_<group>_peaks_sourcedata.csv    <- the bars as drawn
#   Fig_<group>_R2_sourcedata.csv       <- per cell type x resource: n models,
#                                          five-number summary, n above the cap
#   Fig_peaks_R2_violin_sessionInfo.txt
#
# COLOURS
#   Cell types carry the colours of the UMAP, which is ggplot's default hue
#   wheel assigned in ALPHABETICAL cell-type order - verified here:
#   scales::hue_pal()(7) returns exactly the seven hexes hard-coded below.
#   They are then DESATURATED by SATURATION (0.6 of the original HSV
#   saturation) as requested, which keeps the hue identity of each type while
#   taking the glare off large filled areas.  The two polyp-specific states
#   (ASC, SSC) are not in the UMAP, so they use this project's existing ochre
#   and wine, desaturated by the same factor.
#
# HOW TO RUN
#   Rscript "W:/VUMC_projects/scTWAS/Viz/Fig_peaks_R2_violin.R"
###############################################################################

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(data.table)
  library(readxl)
})

## ===========================================================================
## CONFIG
## ===========================================================================

# Primary output folder: the working directory for these two figures, beside
# the model summaries they are built from (set 2026-09-29).  If it is ever
# unreachable the script says so loudly and writes beside itself rather than
# dying half way through.
OUT_DIR <- paste0("F:/Dropbox/GuoLab/COREAD_ToDrGuo/scTWAS/ms_figures_tables/",
                  "res/TWAS/SPrediXcan_peakSNP_association/Figs_peaks_R2_violin")
if (!dir.exists(OUT_DIR)) {
  alt <- getwd()
  message("OUT_DIR is unreachable: ", OUT_DIR, "\n  writing to ", alt, " instead")
  OUT_DIR <- alt
}

PEAKS_XLSX  <- "F:/Dropbox/GuoLab/COREAD_ToDrGuo/scTWAS/ms-guo/SupplementaryTables.xlsx"
PEAKS_SHEET <- "scATAC-seq peaks"
PEAKS_COL   <- "Num_regions_after_chain_merge"

TWAS_ROOT <- paste0("F:/Dropbox/GuoLab/COREAD_ToDrGuo/scTWAS/ms_figures_tables/",
                    "res/TWAS/SPrediXcan_peakSNP_association")

# group -> the model folders that belong to it.  Labels come from the file
# prefixes inside these folders, not from the folder names (see header).
GROUPS <- list(
  "Normal colon" = list(
    file  = "Fig_NormalColon_peaks_R2violin",
    dirs  = c("VUMC_NL", "WellcomeSanger_NL"),
    cts   = c("STM", "TAC", "ABS", "CT", "GOB", "EE", "TUF")
  ),
  "Colon polyps" = list(
    file  = "Fig_ColonPolyps_peaks_R2violin",
    dirs  = "VUMC_ADSER",
    cts   = c("STM", "TAC", "ABS", "CT", "GOB", "EE", "TUF", "ASC", "SSC")
  )
)

# Which performance statistic panel b shows (requested 2026-09-29: the
# cross-validated one).
#   "pred.perf.r2" - CROSS-VALIDATION R^2 of the elastic-net model: a squared
#                    correlation of held-out predictions, so it is never
#                    negative in these files.
#   "R2_global"    - the whole-sample global R^2 used in the earlier version of
#                    this figure.  It IS frequently negative (53% of ADSER, 52%
#                    of VUMC_HTAN_GSE266546_NL, 43% of WellcomeSanger_NL
#                    models), which is why that version filtered much harder.
# The run log prints the counts passing under BOTH, so switching this constant
# never changes the figure silently.
# DIRECTIVE 2026-09-29: use R2_global and NEVER pred.perf.r2.  The reason is
# in the files: pred.perf.r2 is exactly r.global^2 (verified to 1e-15 over all
# 144,713 models) and r_squared_global is bit-identical to r.global, so
# pred.perf.r2 is a squared GLOBAL correlation, not a cross-validated
# statistic; squaring it also hides sign, admitting models whose predictions
# run opposite to measured expression (r.global < 0 for 19.6% of the VUMC and
# 7.6% of the WellcomeSanger models that pass pred.perf.r2 > 0.01).
METRIC_COL <- "R2_global"
METRIC_LAB <- expression(global~italic(R)^2)

R2_MIN <- 0.01     # models are kept when METRIC_COL is GREATER than this
# Panel b starts AT the inclusion threshold (requested 2026-09-29), so the
# axis floor and the filter are the same number and no empty band is drawn
# below the data.  When Y_LO equals R2_MIN the dashed threshold line would fall
# exactly on the axis line, so it is suppressed.
Y_LO   <- 0.01     # y axis floor of panel b
Y_CAP  <- 0.10     # y axis ceiling of panel b
Y_BREAKS <- c(0.01, 0.02, 0.04, 0.06, 0.08, 0.10)

FIG_W <- 7.2         # inches, single-resource figure (polyps)
FIG_W_FACET <- 9.2   # inches when panel b is split into two resource facets
FIG_H <- 3.0
OUT_DPI <- 600
WIDTHS <- c(1, 1.2)        # peaks panel : violin panel, single resource
WIDTHS_FACET <- c(1, 2.0)  # ... when panel b carries two facets of 7 violins
WIDTHS_2x2   <- c(1, 1)    # ... in the 2 x 2 split layout (a|b over c|d)

BASE_PT  <- 7      # axis titles, legend text
TICK_PT  <- 6      # tick labels
VALUE_PT <- 5      # numbers above the bars
TAG_PT   <- 9      # panel letters

BAR_WIDTH  <- 0.72
VIOLIN_W   <- 0.80
DODGE      <- 0.82
SHOW_MEDIAN <- TRUE     # white median crossbar inside each violin
SHOW_N      <- TRUE     # n models above each violin (the two cohorts differ
                        # several-fold in how many genes clear R2_MIN, so the
                        # shapes are not comparable without the counts)

# Cross-cohort concordance scatter, written as its own figure for any group
# with two resources: one point per gene x cell type modelled in BOTH cohorts.
MAKE_CONCORDANCE  <- TRUE
CONC_REQUIRE_BOTH <- TRUE    # TRUE: both cohorts > R2_MIN; FALSE: either one
# Cell types drawn in the in-figure scatter (panel d).  ABS and GOB are the two
# with the most shared genes; the standalone concordance figure keeps all of
# them.  NULL = all.
CONC_CTS_IN_FIG <- c("ABS", "GOB")
# With few cell types in the panel, the pooled rho hides that each lineage has
# its own agreement; TRUE writes one coloured line per cell type instead.
CONC_NOTE_BY_CT <- TRUE
CONC_NOTE_MAX_CT <- 3     # fall back to the pooled note above this many
CONC_W <- 3.9
CONC_H <- 3.6
CONC_BREAKS <- c(0.01, 0.02, 0.05, 0.1, 0.2, 0.5)
CONC_PT <- 0.5
CONC_ALPHA <- 0.55
# Requested 2026-10-05: shape each violin on ALL models passing R2_MIN, not
# only those inside the displayed window.  The y axis still stops at Y_CAP, so
# the thin upper tail of every violin is clipped at the panel edge; the count
# and percentage above the cap stay in the *_R2_sourcedata.csv and in the run
# log, and median_pass (not median_plotted) is now the bar you see.
VIOLIN_WITHIN_CAP <- FALSE

# How panel b arranges two resources:
#   "split" - each resource gets its OWN lettered panel (b and c) on identical
#             y axes, and the figure becomes a 2 x 2 grid a/b/c/d with the
#             cross-cohort scatter as d (requested 2026-10-05).
#   "facet" - one block of 7 violins per resource side by side inside a single
#             panel b, resource named on the strip above its block.
#   "dodge" - the earliest layout: the two resources interleaved within each
#             cell type, distinguished by alpha and outline with a legend.
#   "combined" - ONE panel b holding both cohorts together, dodged within each
#                cell type (requested 2026-10-05).  The figure is then three
#                panels on the same 2-row canvas: a peaks and c scatter on the
#                top row, b spanning the full bottom row so 14 boxes fit.
PANEL_B_LAYOUT <- "combined"

# Geometry for panel b: "violin" or "box".
PANEL_B_GEOM <- "violin"
BOX_WIDTH <- 0.70

# Requested 2026-10-05: panel b shows ONE distribution per cell type with the
# two normal-colon cohorts pooled, no VUMC/WellcomeSanger distinction.  The
# pooled set is dominated by WellcomeSanger (34,626 of 42,625 models, 81%), so
# the per-cohort numbers stay in the source-data CSV and in panel c.  Set FALSE
# to go back to a cohort-resolved panel b (see PANEL_B_LAYOUT).
POOL_RESOURCES <- TRUE

# Only the "split" run writes the shipped file name; the alternative layouts
# get a suffix so re-running them cannot overwrite the 2 x 2 deliverable.
# "facet" therefore lands as Fig_NormalColon_peaks_R2violin_facet2panel.*
LAYOUT_SUFFIX <- c(combined = "", split = "_4panel", facet = "_facet2panel",
                   dodge = "_dodge2panel")
ALPHA_BY_RESOURCE <- c(1, 0.55)   # "dodge" layout only: first solid, second lighter

# Fill opacity per resource in the "facet" layout, matched by name (anything
# not listed is drawn solid).  Requested 2026-10-05: WellcomeSanger slightly
# more transparent than VUMC_HTAN.
RESOURCE_ALPHA <- c(VUMC_HTAN_GSE266546 = 1, VUMC_HTAN_ADSER = 1,
                    WellcomeSanger = 0.8)

# Strip labels for panel b.  The directory prefixes are the dataset ids; the
# figure shows the shorter cohort names the manuscript uses.
RESOURCE_LABELS <- c(VUMC_HTAN_GSE266546 = "VUMC_HTAN",
                     VUMC_HTAN_ADSER = "VUMC_HTAN",
                     WellcomeSanger = "WellcomeSanger")

# UMAP colours = scales::hue_pal()(7) in alphabetical cell-type order.
UMAP_HUES <- c(ABS = "#F8766D", CT = "#C49A00", EE = "#53B400", GOB = "#00C094",
               STM = "#00B6EB", TAC = "#A58AFF", TUF = "#FB61D7")
# Polyp-specific states.  They are absent from the normal-colon UMAP, and the
# polyp UMAP these panels are assembled next to (FigureS_UMAP_Peaks_S2, panel
# A) draws ASC in navy and SSC in dark olive.  Those two hexes were read off
# that figure's legend keys, so all three panels agree on every cell type.
# After desaturation the closest pair among the nine is TAC-TUF at dE 27
# (CIELAB), i.e. nothing collides.
EXTRA_HUES <- c(ASC = "#1F3B73", SSC = "#4F5311")
SATURATION <- 0.6   # multiply HSV saturation by this ("lower saturation")

## ===========================================================================
## HELPERS
## ===========================================================================

# Lower the saturation of a colour without moving its hue or its value, so a
# desaturated palette still reads as "the same colours as the UMAP".
desaturate <- function(cols, f = SATURATION) {
  hsv_mat <- grDevices::rgb2hsv(grDevices::col2rgb(cols))
  out <- grDevices::hsv(hsv_mat["h", ], pmin(1, hsv_mat["s", ] * f), hsv_mat["v", ])
  names(out) <- names(cols)
  out
}

PAL <- c(desaturate(UMAP_HUES), desaturate(EXTRA_HUES))

# Sanity check rather than a promise: say so if this ggplot2/scales build no
# longer produces the hue wheel the UMAP was drawn with.
if (requireNamespace("scales", quietly = TRUE)) {
  hp <- scales::hue_pal()(7)
  if (!identical(unname(UMAP_HUES), hp))
    message("NOTE: scales::hue_pal()(7) is ", paste(hp, collapse = " "),
            "\n      but the UMAP hues hard-coded here are ",
            paste(unname(UMAP_HUES), collapse = " "),
            "\n      Panel colours still follow the hard-coded UMAP hues.")
}

theme_panel <- function() {
  theme_classic(base_size = BASE_PT) +
    theme(
      axis.text        = element_text(size = TICK_PT, colour = "black"),
      axis.title       = element_text(size = BASE_PT, colour = "black"),
      axis.ticks       = element_line(linewidth = 0.25, colour = "black"),
      axis.ticks.length = unit(1.4, "pt"),
      axis.line        = element_line(linewidth = 0.3, colour = "black"),
      legend.key.size  = unit(6, "pt"),
      legend.text      = element_text(size = TICK_PT),
      legend.title     = element_blank(),
      legend.margin    = margin(0, 0, 0, 0),
      legend.box.spacing = unit(2, "pt"),
      legend.background = element_blank(),
      plot.margin      = margin(2, 3, 1, 2)
    )
}

## --- peaks -------------------------------------------------------------------
read_peaks <- function() {
  d <- suppressMessages(
    readxl::read_excel(PEAKS_XLSX, sheet = PEAKS_SHEET, skip = 1,
                       .name_repair = "unique_quiet"))
  d <- as.data.table(d)
  stopifnot(all(c("Group", "CellType", PEAKS_COL) %in% names(d)))
  d <- d[!is.na(CellType)]                       # drops the sub-header row
  d[, .(group = as.character(Group), cell_type = as.character(CellType),
        peaks = suppressWarnings(as.numeric(get(PEAKS_COL))))][
    !is.na(peaks)]
}

## --- model summaries ---------------------------------------------------------
# The resource label is derived from the FILE name, so a folder named after the
# wrong cohort cannot mislabel a violin.
resource_of <- function(path) {
  b <- basename(path)
  if (grepl("^WellcomeSanger", b)) "WellcomeSanger"
  else if (grepl("^VUMC_HTAN_GSE266546", b)) "VUMC_HTAN_GSE266546"
  else if (grepl("^VUMC_HTAN_ADSER", b)) "VUMC_HTAN_ADSER"
  else sub("_[A-Za-z0-9]+_peakSNP.*$", "", b)
}

read_models <- function(dirs) {
  files <- unlist(lapply(dirs, function(d)
    list.files(file.path(TWAS_ROOT, d), pattern = "_peakSNP_model_summaries\\.csv$",
               full.names = TRUE)))
  if (!length(files)) stop("no model summaries under: ",
                           paste(file.path(TWAS_ROOT, dirs), collapse = ", "))
  rbindlist(lapply(files, function(f) {
    d <- fread(f, showProgress = FALSE)
    stopifnot(METRIC_COL %in% names(d), "R2_global" %in% names(d))
    gene_col <- intersect(c("genename", "gene_name", "gene"), names(d))[1]
    data.table(
      resource  = resource_of(f),
      cell_type = sub("^.*_([A-Za-z0-9]+)_peakSNP_model_summaries\\.csv$", "\\1",
                      basename(f)),
      gene_name = as.character(d[[gene_col]]),
      metric    = suppressWarnings(as.numeric(d[[METRIC_COL]])),
      R2_global = suppressWarnings(as.numeric(d$R2_global))
    )
  }))
}

## ===========================================================================
## PANELS
## ===========================================================================

peaks_panel <- function(pk, cts, group_label) {
  pk <- pk[cell_type %in% cts]
  missing <- setdiff(cts, pk$cell_type)
  if (length(missing))
    warning(group_label, ": no peak count for ", paste(missing, collapse = ", "))
  pk[, cell_type := factor(cell_type, levels = cts)]
  # Peaks are drawn in HUNDREDS so the number printed on each bar is a whole
  # number (requested 2026-09-29); in thousands it carried a decimal.
  ggplot(pk, aes(cell_type, peaks / 100, fill = cell_type)) +
    geom_col(width = BAR_WIDTH, colour = NA, show.legend = FALSE) +
    geom_text(aes(label = format(round(peaks / 100), big.mark = ",", trim = TRUE)),
              vjust = -0.45, size = VALUE_PT / .pt, colour = "black") +
    scale_fill_manual(values = PAL) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
    labs(x = group_label, y = expression(scATAC~peaks~(10^2))) +
    theme_panel()
}

r2_panel <- function(md, cts, group_label, title = NULL) {
  md <- md[cell_type %in% cts & is.finite(metric) & metric > R2_MIN]
  if (VIOLIN_WITHIN_CAP) {
    n_over <- md[metric > Y_CAP, .N]
    message(sprintf("  violins shaped on %s models in (%g, %g]; %s (%.1f%%) above the cap, max %.3f",
                    format(md[metric <= Y_CAP, .N], big.mark = ","), R2_MIN, Y_CAP,
                    format(n_over, big.mark = ","), 100 * n_over / nrow(md),
                    max(md$metric)))
    md <- md[metric <= Y_CAP]
  }
  md[, cell_type := factor(cell_type, levels = cts)]
  if (POOL_RESOURCES) md[, resource := "pooled"]
  res <- sort(unique(md$resource))
  md[, resource := factor(resource, levels = res)]
  multi <- length(res) > 1

  facet_mode <- multi && identical(PANEL_B_LAYOUT, "facet")

  # The two cohorts clear R2_MIN with several-fold different numbers of genes,
  # and "scale = width" hides that by drawing every violin the same width.  The
  # count therefore goes on the axis, as a second line of the tick label: a
  # floating label at the top of the panel would sit on the clipped violin
  # tails.  Each facet needs its OWN labels, so x is a resource x cell-type key
  # with free x scales rather than the bare cell type.
  md[, x_key := if (facet_mode) paste(resource, cell_type, sep = "|")
                else as.character(cell_type)]
  key_lv <- if (facet_mode) as.vector(outer(cts, res, function(a, b) paste(b, a, sep = "|")))
            else cts
  key_lv <- key_lv[key_lv %in% md$x_key]
  md[, x_key := factor(x_key, levels = key_lv)]
  x_labs <- setNames(sub("^.*\\|", "", key_lv), key_lv)
  if (SHOW_N) {
    if (facet_mode || !multi) {
      n_by_key <- md[, .N, keyby = x_key]
      n_txt <- format(n_by_key[match(key_lv, x_key), N], big.mark = ",", trim = TRUE)
    } else {
      # Both cohorts share one tick, so the tick carries both counts in the
      # order the boxes are dodged; the axis title names that order.
      n_wide <- dcast(md[, .N, by = .(cell_type, resource)],
                      cell_type ~ resource, value.var = "N", fill = 0)
      n_txt <- apply(n_wide[match(key_lv, as.character(cell_type)), ..res], 1,
                     function(v) paste(format(v, big.mark = ",", trim = TRUE),
                                       collapse = " / "))
    }
    x_labs <- setNames(paste0(x_labs, "\n", n_txt), key_lv)
  }

  p <- ggplot(md, aes(x_key, metric, fill = cell_type))
  if (facet_mode) {
    # Two blocks of seven violins, one per resource, sharing the y scale.
    alph <- RESOURCE_ALPHA[res]
    alph[is.na(alph)] <- 1
    names(alph) <- res
    labs_res <- RESOURCE_LABELS[res]
    labs_res[is.na(labs_res)] <- res[is.na(labs_res)]
    names(labs_res) <- res
    p <- p +
      geom_violin(aes(alpha = resource), width = VIOLIN_W, scale = "width",
                  trim = TRUE, linewidth = 0.2, colour = "grey25",
                  show.legend = FALSE) +
      scale_alpha_manual(values = alph, guide = "none") +
      facet_wrap(~ resource, nrow = 1, scales = "free_x",
                 labeller = labeller(resource = labs_res))
    if (SHOW_MEDIAN && !identical(PANEL_B_GEOM, "box"))
      p <- p + stat_summary(fun = median, geom = "crossbar", width = 0.32,
                            linewidth = 0.25, colour = "white", show.legend = FALSE)
  } else if (multi) {
    p <- p +
      {if (identical(PANEL_B_GEOM, "box"))
         # Outliers are suppressed: with thousands of models per box the cloud
         # would hide the boxes, and the axis stops at Y_CAP anyway.  The box
         # already draws its own median, so no crossbar is added below.
         geom_boxplot(aes(alpha = resource, colour = resource),
                      position = position_dodge(width = DODGE), width = BOX_WIDTH,
                      outlier.shape = NA, linewidth = 0.2)
       else
         geom_violin(aes(alpha = resource, colour = resource),
                     position = position_dodge(width = DODGE), width = VIOLIN_W,
                     scale = "width", trim = TRUE, linewidth = 0.2)} +
      # The keys must carry the alpha difference on a visible swatch: fill is
      # mapped to cell type, so without override.aes the resource keys come out
      # as empty outlined boxes.
      scale_alpha_manual(values = setNames(RESOURCE_ALPHA[res] %|0% ALPHA_BY_RESOURCE[seq_along(res)], res),
                         labels = unname(RESOURCE_LABELS[res] %|0% res),
                         guide = guide_legend(
                           override.aes = list(fill = "grey25", colour = "grey25"))) +
      scale_colour_manual(values = setNames(c("grey15", "grey45")[seq_along(res)], res),
                          labels = unname(RESOURCE_LABELS[res] %|0% res),
                          guide = guide_legend(
                            override.aes = list(fill = c("grey25", "grey25")[seq_along(res)],
                                                alpha = RESOURCE_ALPHA[res] %|0%
                                                        ALPHA_BY_RESOURCE[seq_along(res)])))
    if (SHOW_MEDIAN && !identical(PANEL_B_GEOM, "box"))
      p <- p + stat_summary(aes(group = interaction(x_key, resource)),
                            fun = median, geom = "crossbar", width = 0.32,
                            linewidth = 0.25, colour = "white",
                            position = position_dodge(width = DODGE),
                            show.legend = FALSE)
  } else {
    alph <- RESOURCE_ALPHA[res]
    alph[is.na(alph)] <- 1
    names(alph) <- res
    p <- p +
      {if (identical(PANEL_B_GEOM, "box"))
         geom_boxplot(aes(alpha = resource), width = BOX_WIDTH,
                      outlier.shape = NA, linewidth = 0.2, colour = "grey25",
                      show.legend = FALSE)
       else
         geom_violin(aes(alpha = resource), width = VIOLIN_W, scale = "width",
                     trim = TRUE, linewidth = 0.2, colour = "grey25",
                     show.legend = FALSE)} +
      scale_alpha_manual(values = alph, guide = "none")
    if (SHOW_MEDIAN && !identical(PANEL_B_GEOM, "box"))
      p <- p + stat_summary(fun = median, geom = "crossbar", width = 0.32,
                            linewidth = 0.25, colour = "white", show.legend = FALSE)
  }

  p +
    {if (!isTRUE(all.equal(Y_LO, R2_MIN)))
       geom_hline(yintercept = R2_MIN, linetype = "dashed", linewidth = 0.3,
                  colour = "grey30")} +
    scale_fill_manual(values = PAL, guide = "none") +
    scale_x_discrete(labels = x_labs) +
    scale_y_continuous(breaks = Y_BREAKS,
                       labels = function(b) sub("^0\\.10$", "0.10", sprintf("%.2f", b)),
                       expand = expansion(mult = c(0.02, 0.03))) +
    coord_cartesian(ylim = c(Y_LO, Y_CAP)) +   # clip the tail, keep the rows
    labs(x = if (multi && !facet_mode && SHOW_N)
              sprintf("%s  (n models: %s)", group_label,
                      paste(unname(RESOURCE_LABELS[res] %|0% res), collapse = " / "))
            else group_label,
         y = METRIC_LAB, title = title) +
    theme_panel() +
    theme(plot.title = element_text(size = BASE_PT, colour = "black",
                                    hjust = 0.5, margin = margin(0, 0, 2, 0))) +
    theme(legend.position = if (multi && !facet_mode) "top" else "none",
          legend.direction = "horizontal",
          legend.key.size = unit(9, "pt"),
          strip.background = element_blank(),
          strip.text = element_text(size = BASE_PT, colour = "black",
                                    margin = margin(0, 0, 1.5, 0)),
          panel.spacing.x = unit(5, "pt"))
}

## One point per gene x cell type that BOTH cohorts modelled.  The axes are
## log10 because the metric spans 0.01 to ~0.8 and a linear pair of axes puts
## nine tenths of the points into the bottom corner; coord_fixed keeps the
## y = x reference at 45 degrees, so distance from it reads directly as the
## gap between the two cohorts for that gene.
`%|%` <- function(a, b) if (is.na(a)) b else a
`%|0%` <- function(a, b) { a <- unname(a); b <- unname(b); a[is.na(a)] <- b[is.na(a)]; a }

# Point fills are desaturated for the violins; the same hue as text on white
# is too light to read, so annotation text uses a darkened version.
darken <- function(cols, f = 0.65) {
  m <- grDevices::col2rgb(cols) * f
  grDevices::rgb(m[1, ], m[2, ], m[3, ], maxColorValue = 255)
}

concordance_figure <- function(md, cts, x_res, y_res, group_label) {
  m <- md[cell_type %in% cts & is.finite(metric)]
  w <- dcast(m, cell_type + gene_name ~ resource, value.var = "metric",
             fun.aggregate = function(v) if (length(v)) max(v) else NA_real_,
             fill = NA_real_)
  setnames(w, c(x_res, y_res), c("x", "y"))
  w <- w[is.finite(x) & is.finite(y)]
  w <- if (CONC_REQUIRE_BOTH) w[x > R2_MIN & y > R2_MIN] else w[x > R2_MIN | y > R2_MIN]
  if (!nrow(w)) return(NULL)
  w[, cell_type := factor(cell_type, levels = cts)]

  rho <- suppressWarnings(cor(w$x, w$y, method = "spearman"))
  r   <- suppressWarnings(cor(w$x, w$y))
  lim <- range(c(w$x, w$y))
  lim <- c(10^floor(log10(lim[1])), 10^ceiling(log10(lim[2])))
  by_ct <- w[, .(n = .N,
                 rho = suppressWarnings(cor(x, y, method = "spearman")),
                 r = suppressWarnings(cor(x, y)),
                 above = 100 * mean(y > x)), by = cell_type][order(cell_type)]
  per_ct_note <- CONC_NOTE_BY_CT && nrow(by_ct) <= CONC_NOTE_MAX_CT && nrow(by_ct) > 1
  note <- sprintf("n = %s gene \u00d7 cell type\nSpearman \u03c1 = %.2f, r = %.2f\n%.0f%% above the diagonal",
                  format(nrow(w), big.mark = ","), rho, r, 100 * mean(w$y > w$x))
  by_ct[, lab := sprintf("%s  n = %s, \u03c1 = %.2f, r = %.2f",
                         cell_type, format(n, big.mark = ","), rho, r)]
  xl <- unname(RESOURCE_LABELS[x_res] %|% x_res)
  yl <- unname(RESOURCE_LABELS[y_res] %|% y_res)

  p <- ggplot(w, aes(x, y, colour = cell_type)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed",
                linewidth = 0.3, colour = "grey30") +
    geom_point(size = CONC_PT, alpha = CONC_ALPHA, stroke = 0) +
    {if (per_ct_note)
       annotate("text", x = lim[1], y = lim[2] * 10^(-0.07 * (seq_len(nrow(by_ct)) - 1)),
                hjust = 0, vjust = 1, label = by_ct$lab, size = TICK_PT / .pt,
                colour = darken(PAL[as.character(by_ct$cell_type)], 0.65))
     else
       annotate("text", x = lim[1], y = lim[2], hjust = 0, vjust = 1,
                label = note, size = TICK_PT / .pt, colour = "grey20",
                lineheight = 1.15)} +
    scale_colour_manual(values = PAL, name = NULL,
                        guide = guide_legend(nrow = 1,
                          override.aes = list(size = 1.5, alpha = 1))) +
    scale_x_log10(breaks = CONC_BREAKS, labels = CONC_BREAKS, limits = lim) +
    scale_y_log10(breaks = CONC_BREAKS, labels = CONC_BREAKS, limits = lim) +
    coord_fixed() +
    labs(x = bquote(.(xl) ~ global ~ italic(R)^2),
         y = bquote(.(yl) ~ global ~ italic(R)^2),
         title = if (is.null(group_label)) NULL else group_label) +
    theme_panel() +
    theme(legend.position = "top", legend.direction = "horizontal",
          legend.key.size = unit(7, "pt"),
          plot.title = element_text(size = BASE_PT, colour = "black",
                                    hjust = 0.5, margin = margin(0, 0, 2, 0)))
  list(plot = p, data = w, rho = rho, r = r)
}

## ===========================================================================
## BUILD
## ===========================================================================

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
PEAKS <- read_peaks()
message("peaks sheet: ", nrow(PEAKS), " rows (",
        paste(sprintf("%s %d", names(table(PEAKS$group)), table(PEAKS$group)),
              collapse = "; "), ")")

for (g in names(GROUPS)) {
  spec <- GROUPS[[g]]
  message("\n== ", g, " ==")

  pk <- PEAKS[group == g]
  if (!nrow(pk)) stop("no peak rows with Group == '", g, "'")

  md <- read_models(spec$dirs)
  kept <- md[is.finite(metric) & metric > R2_MIN]
  message(sprintf("  models: %s read, %s with %s > %g (%.1f%%)",
                  format(nrow(md), big.mark = ","), format(nrow(kept), big.mark = ","),
                  METRIC_COL, R2_MIN, 100 * nrow(kept) / nrow(md)))
  message(sprintf("  for comparison: %s pass R2_global > %g (%.1f%%); %s of all models have R2_global < 0",
                  format(md[is.finite(R2_global) & R2_global > R2_MIN, .N], big.mark = ","), R2_MIN,
                  100 * md[is.finite(R2_global) & R2_global > R2_MIN, .N] / nrow(md),
                  format(md[is.finite(R2_global) & R2_global < 0, .N], big.mark = ",")))
  for (r in sort(unique(md$resource)))
    message("    ", r, ": ",
            paste(sprintf("%s %s", spec$cts,
                          format(sapply(spec$cts, function(ct)
                            kept[resource == r & cell_type == ct, .N]),
                            big.mark = ",", trim = TRUE)), collapse = " / "))

  res_fig <- sort(unique(kept$resource))
  # With the cohorts pooled there is no cohort-resolved panel and no in-figure
  # scatter: the figure is a + b only, the same two-panel shape as the polyp
  # figure (requested 2026-10-06).  The cross-cohort scatter is still written
  # as its own file below.
  combined <- length(res_fig) == 2 && identical(PANEL_B_LAYOUT, "combined") &&
              !POOL_RESOURCES
  split_2x2 <- length(res_fig) == 2 && identical(PANEL_B_LAYOUT, "split") &&
               !POOL_RESOURCES
  faceted   <- length(res_fig) > 1 && identical(PANEL_B_LAYOUT, "facet")

  if (combined) {
    # a peaks (top left), c scatter (top right), b across the bottom row: the
    # 14 dodged boxes need the full width, the square scatter does not.
    cts_d <- if (is.null(CONC_CTS_IN_FIG)) spec$cts
             else intersect(spec$cts, CONC_CTS_IN_FIG)
    cd <- concordance_figure(md, cts_d, res_fig[1], res_fig[2], NULL)
    if (is.null(cd)) stop("no gene x cell type pairs for the scatter panel")
    fig <- peaks_panel(pk, spec$cts, g) + r2_panel(md, spec$cts, g) + cd$plot +
      plot_layout(design = "AC\nBB", widths = c(1, 1), heights = c(1, 1)) +
      plot_annotation(tag_levels = "a",
                      theme = theme(plot.margin = margin(1, 2, 1, 2))) &
      theme(plot.tag = element_text(size = TAG_PT, face = "bold"))
    fig_w <- FIG_W
    fig_h <- FIG_H * 2
  } else if (split_2x2) {
    # 2 x 2: a peaks, b first cohort, c second cohort, d cross-cohort scatter.
    # Same width as the single-resource figure, twice the height.
    lab_of <- function(r) unname(RESOURCE_LABELS[r] %|% r)
    cts_d  <- if (is.null(CONC_CTS_IN_FIG)) spec$cts
              else intersect(spec$cts, CONC_CTS_IN_FIG)
    cd <- concordance_figure(md, cts_d, res_fig[1], res_fig[2], NULL)
    if (is.null(cd)) stop("no gene x cell type pairs for the scatter panel")
    fig <- peaks_panel(pk, spec$cts, g) +
      r2_panel(md[resource == res_fig[1]], spec$cts, g, title = lab_of(res_fig[1])) +
      r2_panel(md[resource == res_fig[2]], spec$cts, g, title = lab_of(res_fig[2])) +
      cd$plot +
      plot_layout(ncol = 2, nrow = 2, widths = WIDTHS_2x2) +
      plot_annotation(tag_levels = "a",
                      theme = theme(plot.margin = margin(1, 2, 1, 2))) &
      theme(plot.tag = element_text(size = TAG_PT, face = "bold"))
    fig_w <- FIG_W
    fig_h <- FIG_H * 2
  } else {
    fig_w  <- if (faceted) FIG_W_FACET else FIG_W
    fig_h  <- FIG_H
    widths <- if (faceted) WIDTHS_FACET else WIDTHS
    fig <- peaks_panel(pk, spec$cts, g) + r2_panel(md, spec$cts, g) +
      plot_layout(ncol = 2, widths = widths) +
      plot_annotation(tag_levels = "a",
                      theme = theme(plot.margin = margin(1, 2, 1, 2))) &
      theme(plot.tag = element_text(size = TAG_PT, face = "bold"))
  }

  # Single-resource groups (the polyps) have only one possible layout, so they
  # keep the plain name whatever PANEL_B_LAYOUT says.
  sfx <- if (length(res_fig) < 2 || POOL_RESOURCES) ""
         else unname(LAYOUT_SUFFIX[PANEL_B_LAYOUT] %|% "")
  stem <- paste0(spec$file, sfx)
  ggsave(file.path(OUT_DIR, paste0(stem, ".pdf")), fig,
         width = fig_w, height = fig_h, device = grDevices::cairo_pdf, bg = "white")
  ggsave(file.path(OUT_DIR, paste0(stem, ".png")), fig,
         width = fig_w, height = fig_h, dpi = OUT_DPI, bg = "white")

  ## source data -------------------------------------------------------------
  fwrite(pk[cell_type %in% spec$cts][
    , .(group, cell_type, peaks_raw = peaks, peaks_hundreds = round(peaks / 100),
          peaks_thousands = peaks / 1000)][
    order(match(cell_type, spec$cts))],
    file.path(OUT_DIR, paste0(spec$file, "_peaks_sourcedata.csv")))

  # Two summaries per cell type x resource, because the violin and the filter
  # describe different sets and conflating them would misreport the figure:
  #   *_pass    = every model with R2_global > R2_MIN
  #   *_plotted = the subset the violin is actually shaped on, i.e. also
  #               <= Y_CAP when VIOLIN_WITHIN_CAP is TRUE.  median_plotted is
  #               the white crossbar drawn in the panel.
  kept_ct <- kept[cell_type %in% spec$cts]
  drawn <- if (VIOLIN_WITHIN_CAP) kept_ct[metric <= Y_CAP] else kept_ct
  stat <- merge(
    kept_ct[, .(n_models_pass = .N,
                min_pass = min(metric), q1_pass = quantile(metric, 0.25),
                median_pass = median(metric), q3_pass = quantile(metric, 0.75),
                max_pass = max(metric),
                n_above_y_cap = sum(metric > Y_CAP),
                pct_above_y_cap = round(100 * sum(metric > Y_CAP) / .N, 2)),
            by = .(resource, cell_type)],
    drawn[, .(n_models_plotted = .N,
              min_plotted = min(metric), q1_plotted = quantile(metric, 0.25),
              median_plotted = median(metric), q3_plotted = quantile(metric, 0.75),
              max_plotted = max(metric)),
          by = .(resource, cell_type)],
    by = c("resource", "cell_type"), all.x = TRUE
  )[order(resource, match(cell_type, spec$cts))]
  # When panel b pools the cohorts, the violin and its median belong to the
  # pooled set, so that set needs its own rows in the source data - otherwise
  # nothing in the file reproduces the crossbar the reader sees.
  if (POOL_RESOURCES && length(res_fig) > 1) {
    pooled <- merge(
      kept_ct[, .(resource = "pooled", n_models_pass = .N,
                  min_pass = min(metric), q1_pass = quantile(metric, 0.25),
                  median_pass = median(metric), q3_pass = quantile(metric, 0.75),
                  max_pass = max(metric),
                  n_above_y_cap = sum(metric > Y_CAP),
                  pct_above_y_cap = round(100 * sum(metric > Y_CAP) / .N, 2)),
              by = .(cell_type)],
      drawn[, .(resource = "pooled", n_models_plotted = .N,
                min_plotted = min(metric), q1_plotted = quantile(metric, 0.25),
                median_plotted = median(metric), q3_plotted = quantile(metric, 0.75),
                max_plotted = max(metric)),
            by = .(cell_type)],
      by = c("resource", "cell_type"), all.x = TRUE
    )[order(match(cell_type, spec$cts))]
    stat <- rbind(pooled, stat)
    stat[, pooled_from := paste(res_fig, collapse = " + ")]
  }
  stat[, `:=`(group = g, metric = METRIC_COL, metric_min_filter = R2_MIN,
              y_axis_from = Y_LO, y_cap = Y_CAP,
              violin_within_cap = VIOLIN_WITHIN_CAP)]
  fwrite(stat, file.path(OUT_DIR, paste0(spec$file, "_R2_sourcedata.csv")))

  ## cross-cohort concordance ------------------------------------------------
  res_all <- sort(unique(kept$resource))
  if (MAKE_CONCORDANCE && length(res_all) == 2) {
    cf <- concordance_figure(md, spec$cts, res_all[1], res_all[2], g)
    if (!is.null(cf)) {
      f <- paste0(spec$file, "_concordance")
      ggsave(file.path(OUT_DIR, paste0(f, ".pdf")), cf$plot,
             width = CONC_W, height = CONC_H, device = grDevices::cairo_pdf, bg = "white")
      ggsave(file.path(OUT_DIR, paste0(f, ".png")), cf$plot,
             width = CONC_W, height = CONC_H, dpi = OUT_DPI, bg = "white")
      out <- copy(cf$data)
      setnames(out, c("x", "y"), paste0("R2_global_", res_all))
      fwrite(out[order(cell_type, gene_name)],
             file.path(OUT_DIR, paste0(f, "_sourcedata.csv")))
      summ <- function(d) d[, .(n_pairs = .N,
                                spearman_rho = round(cor(x, y, method = "spearman"), 3),
                                pearson_r = round(cor(x, y), 3),
                                median_x = round(median(x), 4),
                                median_y = round(median(y), 4),
                                pct_y_gt_x = round(100 * mean(y > x), 1))]
      per_ct <- rbind(cf$data[, summ(.SD), by = cell_type],
                      cbind(cell_type = "all", summ(cf$data)))
      per_ct[, `:=`(x_resource = res_all[1], y_resource = res_all[2],
                    metric = METRIC_COL, require_both_above_R2_MIN = CONC_REQUIRE_BOTH,
                    R2_min_filter = R2_MIN)]
      fwrite(per_ct, file.path(OUT_DIR, paste0(f, "_percelltype.csv")))
      message(sprintf("  wrote %s.pdf/.png (%s pairs, Spearman %.2f, Pearson %.2f)",
                      f, format(nrow(cf$data), big.mark = ","), cf$rho, cf$r))
    }
  }

  message("  wrote ", stem, ".pdf/.png + 2 source-data CSVs")
}

writeLines(capture.output(sessionInfo()),
           file.path(OUT_DIR, "Fig_peaks_R2_violin_sessionInfo.txt"))
message("\ndone: ", OUT_DIR)
