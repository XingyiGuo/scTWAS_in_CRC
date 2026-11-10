## 04_de_support_grouped_bars.R -------------------------------------------
##
## Grouped bar chart of Up/Down support across stage contrasts.
##
##   x            the 8 stage contrasts (contrast_short), grouped into the three
##                routes from normal epithelium to carcinoma
##   within each  14 cell-type models, dodged: all 7 cell types x NL and Polyps
##   above zero   count of Up-supported genes, warm ramps
##   below zero   count of Down-supported genes, cool ramps
##   hue family   red/blue for the normal-colon models, purple/teal for polyp
##
## Reads the long count table from 03_summarize_scrna_support_counts.R.
##
## Three design decisions worth knowing
##
## 1. Two ramps, one fill scale. ggplot allows a single fill scale per plot, so
##    the 18 combinations of direction x cell type are encoded as one factor
##    (fill_key) with 18 manual colours: a warm ramp for Up and a cool ramp for
##    Down, both walking the SAME cell-type order. A given position in either
##    ramp therefore means the same cell type, which is what makes the halves
##    readable against each other.
##
## 2. Down counts are plotted negative but the axis is labelled with absolute
##    values, because the bars are counts and a "-14 genes" tick would be wrong.
##    The axis title names which side is which.
##
## 3. Denominators live in the legend. Each model tests a different number of
##    gene x cell-type pairs (NL ABS 93 down to Polyps GOB 27), so a taller bar
##    is not necessarily a higher rate. The legend label carries n for each
##    model, the same convention the count grid used on its y axis.
##
## Output: figure/Fig_scTWAS_DEsupport_UpDown_bars.{png,pdf} at 300 dpi
##         + _sourcedata.csv (every plotted count)
## -----------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(data.table); library(ggplot2)
})

PROJ <- Sys.getenv("SC_TWAS_SUPPORT_ROOT", "W:/VUMC_projects/scTWAS/MultiOmics_DE_support")
LONG_IN <- file.path(PROJ, "table", "DE_key_results",
                     "scTWAS_DEsupport_summary_UpDown_long.csv")
DIR_FIG <- file.path(PROJ, "figure")
dir.create(DIR_FIG, showWarnings = FALSE, recursive = TRUE)

## ---- the 14 models, in the project's cell-type order -----------------------
## All seven cell types in both datasets. The five polyp models added here have
## small denominators (TAC 18, STM 15, TUF 2, EE 4, CT 3 pairs), so their bars
## are shown with n in the legend and flagged in the caption: a zero there is
## sample size, not biology.
CT_ORDER <- c("ABS", "GOB", "TAC", "STM", "TUF", "EE", "CT")
MODELS <- rbind(
  data.table(Dataset = "NL",    CellType = CT_ORDER),
  data.table(Dataset = "ADSER", CellType = CT_ORDER))
MODELS[, model := fifelse(Dataset == "NL", paste("NL", CellType),
                          paste("Polyp", CellType))]

## ---- input -----------------------------------------------------------------
lg <- fread(LONG_IN, encoding = "UTF-8")
d <- merge(lg, MODELS, by = c("Dataset", "CellType"))
stopifnot(nrow(d) == nrow(MODELS) * 8L * 2L)         # 14 models x 8 contrasts x 2

## contrast order as exported, and the from-stage grouping that makes the
## carcinoma contrasts read separately from the polyp ones
CON_ORDER <- c("Normal \u2192 MSI-H", "Normal \u2192 MSS",
               "Normal \u2192 ASC", "Normal \u2192 SSC",
               "ASC \u2192 MSI-H", "ASC \u2192 MSS",
               "SSC \u2192 MSI-H", "SSC \u2192 MSS")
stopifnot(setequal(unique(d$contrast_short), CON_ORDER))
d[, contrast_short := factor(contrast_short, levels = CON_ORDER)]
## Three sections, each a route from normal epithelium to carcinoma: the direct
## one, and the two that pass through a polyp. Grouping on the route rather than
## on the from-stage keeps each lesion's entry and exit contrasts adjacent.
SECTION <- data.table(
  contrast_short = CON_ORDER,
  section = c(rep("Normal \u2192 carcinoma", 2),
              "Normal \u2192 ASC \u2192 carcinoma", "Normal \u2192 SSC \u2192 carcinoma",
              rep("Normal \u2192 ASC \u2192 carcinoma", 2),
              rep("Normal \u2192 SSC \u2192 carcinoma", 2)))
SEC_ORDER <- c("Normal \u2192 carcinoma", "Normal \u2192 ASC \u2192 carcinoma",
               "Normal \u2192 SSC \u2192 carcinoma")
d <- merge(d, SECTION, by = "contrast_short")
d[, section := factor(section, levels = SEC_ORDER)]
## merge() reorders rows, so check by section name rather than by position
.chk <- d[, uniqueN(contrast_short), by = section]
stopifnot(identical(.chk[match(SEC_ORDER, as.character(section))]$V1,
                    c(2L, 3L, 3L)))
## within a section, contrasts run in route order: entry first, then the exits
d[, contrast_short := factor(as.character(contrast_short), levels = CON_ORDER)]

d[, model := factor(model, levels = MODELS$model)]
d[, direction := factor(direction, levels = c("Up", "Down"))]

## signed height: Down below the axis
d[, y := fifelse(direction == "Up", n_supported, -n_supported)]
d[, fill_key := paste(direction, model)]

## ---- palettes: warm above, cool below, one gradient family per dataset ------
## Four ramps, not two. Direction is carried by temperature (warm above the
## axis, cool below) and DATASET by hue family: the normal-colon models use
## red/blue, the polyp models purple/teal. Within each ramp the seven cell types
## run light to dark in the project order, so the same depth means the same cell
## type in all four ramps.
ramp <- function(cols) colorRampPalette(cols)(length(CT_ORDER))
NL_UP    <- ramp(c("#FFEBB4", "#FDC161", "#F59B42", "#E3603F", "#B2182B", "#6E0D1B"))
NL_DOWN  <- ramp(c("#DCEBF7", "#A8CDE6", "#6BA8D0", "#3B7CB3", "#21548F", "#0B2F5E"))
PY_UP    <- ramp(c("#FBE3EF", "#F3B6D6", "#E084B4", "#C4568E", "#98306A", "#611443"))
PY_DOWN  <- ramp(c("#DDF0EC", "#A9DBD2", "#6DBFB3", "#3C9B92", "#1F7169", "#0A4641"))

nl_models <- MODELS[Dataset == "NL",    model]
py_models <- MODELS[Dataset == "ADSER", model]
FILL <- c(setNames(NL_UP,   paste("Up",   nl_models)),
          setNames(NL_DOWN, paste("Down", nl_models)),
          setNames(PY_UP,   paste("Up",   py_models)),
          setNames(PY_DOWN, paste("Down", py_models)))
stopifnot(length(FILL) == 2L * nrow(MODELS), !anyDuplicated(names(FILL)))

## legend: two rows, normal-colon models above polyp models, so the two gradient
## families read as families. Keys use the Up colour; the cool ramps are keyed by
## position, each model's Down bar sitting directly beneath its Up bar.
DENOM <- d[direction == "Up", .(n = first(n_genes)), by = model]
DENOM <- setNames(DENOM$n, DENOM$model)
BREAKS <- paste("Up", MODELS$model)
LABELS_WARM <- sprintf("%s (n=%d)", MODELS$model, DENOM[MODELS$model])

## ---- figure ----------------------------------------------------------------
YMAX <- max(d$n_supported)
DIRLAB <- data.table(
  section = factor(SEC_ORDER[1], levels = SEC_ORDER),
  x = 0.52,
  y = c(YMAX * 0.99, -YMAX * 0.70),
  label = c("Up (warm)", "Down (cool)"))
INK <- "#0b0b0b"; INK2 <- "#52514e"; INK3 <- "#8a8880"

p <- ggplot(d, aes(x = contrast_short, y = y, fill = fill_key)) +
  geom_hline(yintercept = 0, colour = INK, linewidth = 0.4) +
  ## group = model is load-bearing: dodging on fill_key would give each of the
  ## 18 direction x model combinations its own slot, so a model's Up bar and its
  ## Down bar would sit at different x positions instead of mirroring.
  geom_col(aes(group = model),
           position = position_dodge(width = 0.86, preserve = "single"),
           width = 0.78, colour = "white", linewidth = 0.08) +
  facet_grid(. ~ section, scales = "free_x", space = "free_x") +
  scale_fill_manual(values = FILL, breaks = BREAKS, labels = LABELS_WARM,
                    name = NULL,
                    guide = guide_legend(nrow = 2, byrow = TRUE,
                                         keyheight = unit(3.1, "mm"),
                                         keywidth  = unit(3.1, "mm"))) +
  ## the cool half is keyed by position, not by a second set of swatches.
  ## annotate() would repeat these in every facet, so they are drawn from a
  ## one-row-per-label data frame pinned to the left-hand facet only.
  geom_text(data = DIRLAB, aes(x = x, y = y, label = label, colour = label),
            inherit.aes = FALSE, hjust = 0, size = 2.5, fontface = "bold",
            show.legend = FALSE) +
  scale_colour_manual(values = c("Up (warm)" = "#C1272D",
                                 "Down (cool)" = "#21548F"), guide = "none") +
  scale_y_continuous(labels = function(v) abs(v),
                     breaks = scales::pretty_breaks(8),
                     expand = expansion(mult = c(0.05, 0.07))) +
  labs(x = NULL,
       y = "Direction-concordant supported genes\n(Up above axis, Down below)",
       title = paste0("Up- and down-supported scTWAS genes per stage contrast, ",
                      "split by cell-type model"),
       subtitle = paste0("14 models: all 7 epithelial cell types in the ",
                         "normal-colon and polyp datasets"),
       caption = paste0(
         "A gene counts once per model x contrast when it is significant at FDR < 0.05, |log2FC| > 0.5, and sign(log2FC) matches the sign of its scTWAS Z in that\n",
         "model - so Up and Down are both support, and record the direction of the expression change. Sections are the three routes from normal epithelium to\n",
         "carcinoma: direct, via a conventional adenoma (ASC), and via a sessile serrated lesion (SSC). Within a section the lesion entry contrast comes first,\n",
         "then its two exits. Warm bars above the axis and cool below; the normal-colon models use the red/blue pair and the polyp models the purple/teal pair,\n",
         "with the seven cell types running light to dark in the same order in all four ramps. Bars are gene counts, not rates: each model tests a different\n",
         "number of gene x cell-type pairs (n in the legend), from NL ABS 93 down to Polyp TUF 2, so a taller bar can reflect the larger denominator rather than\n",
         "a stronger signal. The five smallest polyp models in particular (TAC 18, STM 15, CT 3, EE 4, TUF 2 pairs) cannot show much either way, and a zero there\n",
         "is sample size rather than biology.")) +
  theme_minimal(base_size = 9) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor   = element_blank(),
    panel.grid.major.y = element_line(colour = "grey91", linewidth = 0.25),
    strip.text         = element_text(size = 8.4, face = "bold", colour = INK),
    strip.background   = element_rect(fill = "grey96", colour = NA),
    axis.text.x        = element_text(size = 8, colour = INK, angle = 20,
                                      hjust = 1),
    axis.text.y        = element_text(size = 7.5, colour = INK2),
    axis.title.y       = element_text(size = 8, colour = INK2,
                                      lineheight = 1.15),
    plot.title         = element_text(size = 10, face = "bold", colour = INK,
                                      margin = margin(b = 1)),
    plot.subtitle      = element_text(size = 8, colour = INK2,
                                      margin = margin(b = 5)),
    plot.caption       = element_text(size = 5.4, colour = INK3, hjust = 0,
                                      lineheight = 1.3, margin = margin(t = 5)),
    legend.position    = "top",
    legend.text        = element_text(size = 6.6, colour = INK2),
    legend.margin      = margin(b = 2),
    plot.margin        = margin(6, 8, 4, 6))

W <- 14.4; H <- 6.0
ggsave(file.path(DIR_FIG, "Fig_scTWAS_DEsupport_UpDown_bars.png"), p,
       width = W, height = H, dpi = 300, bg = "white")
ggsave(file.path(DIR_FIG, "Fig_scTWAS_DEsupport_UpDown_bars.pdf"), p,
       width = W, height = H, bg = "white", device = cairo_pdf)

## ---- source data -----------------------------------------------------------
sd <- d[, .(section, contrast, contrast_short, Dataset, CellType, model,
            n_pairs, n_genes, direction, n_supported, plotted_y = y)]
setorder(sd, section, contrast_short, model, direction)
fwrite(sd, file.path(DIR_FIG, "Fig_scTWAS_DEsupport_UpDown_bars_sourcedata.csv"))

cat("bars:", nrow(d), "| tallest:", YMAX,
    "| Up total:", d[direction == "Up", sum(n_supported)],
    "| Down total:", d[direction == "Down", sum(n_supported)], "\n")

## ---- compact 18 x 5 cm version, for placement as a composite panel ---------
## Requested: 18 cm wide x 5 cm tall = 7.09 x 1.97 in, tight, axis ticks 8-11 pt.
##
## What is removed and why. At this height the title, subtitle and the nine-line
## caption cannot be set at a legible size and still leave room for bars, and
## they are not needed on a panel that sits inside a composite whose legend
## text carries them. Removed: title, subtitle, caption, and the "(n=..)"
## denominators from the legend labels - at 8 pt the longest label with its
## denominator is about 1.2 in, so seven of them exceed the 7.09 in width,
## while the bare model names fit in roughly 5.6 in.
##
## TWO THINGS THE COMPOSITE LEGEND MUST NOW STATE, since the figure no longer
## does: (1) bars are gene counts and each model tests a different number of
## gene x cell-type pairs, from NL ABS 93 down to Polyp TUF 2, so a taller bar
## can reflect the larger denominator rather than a stronger signal; (2) the
## five smallest polyp models (TAC 18, STM 15, EE 4, CT 3, TUF 2 pairs) cannot
## show much either way, so a zero there is sample size, not biology.
## The denominators remain in the _sourcedata.csv.
CM <- 1 / 2.54
FS_TICK <- 8          # x and y tick labels, the size asked for
LABELS_C <- MODELS$model

p_c <- p +
  scale_fill_manual(values = FILL, breaks = BREAKS, labels = LABELS_C,
                    name = NULL,
                    guide = guide_legend(nrow = 2, byrow = TRUE,
                                         keyheight = unit(2.6, "mm"),
                                         keywidth  = unit(2.6, "mm"))) +
  labs(title = NULL, subtitle = NULL, caption = NULL,
       y = "Supported genes") +
  theme(
    axis.text.x   = element_text(size = FS_TICK, colour = INK, angle = 20,
                                 hjust = 1, vjust = 1),
    axis.text.y   = element_text(size = FS_TICK, colour = INK2),
    axis.title.y  = element_text(size = FS_TICK, colour = INK2,
                                 margin = margin(r = 1)),
    strip.text    = element_text(size = FS_TICK, face = "bold", colour = INK,
                                 margin = margin(1, 1, 1, 1)),
    legend.text   = element_text(size = FS_TICK, colour = INK2),
    legend.position = "top",
    legend.margin = margin(0, 0, 1, 0),
    legend.box.spacing = unit(1, "mm"),
    legend.spacing.x = unit(1, "mm"),
    panel.spacing = unit(2, "mm"),
    plot.margin   = margin(2, 3, 1, 2))

## the in-panel direction labels shrink with the panel; keep them readable but
## out of the bars' way
p_c$layers[[which(vapply(p_c$layers,
                         function(L) inherits(L$geom, "GeomText"),
                         logical(1)))[1]]]$aes_params$size <- 5.5 / .pt

STEM_C <- file.path(DIR_FIG, "Fig_scTWAS_DEsupport_UpDown_bars_18x5cm")
ggsave(paste0(STEM_C, ".png"), p_c, width = 18 * CM, height = 5 * CM,
       units = "in", dpi = 600)
ggsave(paste0(STEM_C, ".pdf"), p_c, width = 18 * CM, height = 5 * CM,
       units = "in", device = cairo_pdf)
cat("wrote", paste0(STEM_C, ".{png,pdf}"), "at 18 x 5 cm, ticks", FS_TICK, "pt\n")
