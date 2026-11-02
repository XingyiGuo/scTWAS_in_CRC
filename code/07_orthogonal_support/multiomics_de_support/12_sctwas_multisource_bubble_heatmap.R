## 12_sctwas_multisource_bubble_heatmap.R --------------------------------------
##
## Redesign of Fig_scTWAS_multisource_DEsupport_diverging as two bubble
## heatmaps, split by the direction of the genetic effect.
##
##   panel a  oncogene-like        rep_Z > 0
##   panel b  tumour-suppressor-like  rep_Z < 0
##
## Four rows per panel, top to bottom:
##   1  scTWAS          bubble colour = rep_Z,                 size = rep_P
##   2  CosMx NOR->CORE bubble colour = log2FC,                size = padj
##   3  TCGA            bubble colour = log2FC,                size = padj
##   4  CPTAC           bubble colour = log2FC (protein),      size = padj
##
## Genes shown: the 109 with at least one Up/Down support call among CosMx,
## TCGA and CPTAC (57 oncogene-like, 52 suppressor-like). Script 05 is left
## untouched; this reads the same canonical table.
##
## Four decisions the brief left open, each visible in the figure
##
## 1. Two colour scales, not four. Row 1 is a Z score and rows 2-4 are log2
##    fold changes, which are different units, so they cannot share a ramp. But
##    the three fold-change rows DO share one scale, which is the point - it
##    makes CosMx, TCGA and CPTAC directly comparable to each other. ggplot
##    allows one fill scale per plot, so each panel is built as two stacked
##    sub-plots (Z row, then the three log2FC rows) and composed.
##
## 2. Bubble size is a five-level significance bin, not continuous -log10(p).
##    The p-value ranges are wildly unequal: CosMx padj bottoms out at 4.2e-05
##    because this project recomputes its FDR over the risk genes alone, while
##    TCGA reaches 7.4e-119 genome-wide. On a shared continuous size scale every
##    CosMx bubble would render as a dot, implying CosMx is weak when it is
##    merely differently normalised. Bins compress that and stay readable.
##    The consequence to keep: sizes are comparable WITHIN a row, not across
##    rows, and the caption says so.
##
## 3. Colour is capped, not clipped. log2FC limits are +/-2.5 and Z limits
##    +/-13, with values beyond squished to the end colour rather than dropped
##    (7 of 259 fold changes and 5 of 109 Z scores). Nothing is removed.
##
## 4. Untested gene x resource combinations are drawn as a small grey cross, not
##    left blank. CPTAC measures only 57 of these 109 genes and CosMx 93, so
##    blank cells would read as null results rather than absent measurements.
##
## Output: figure/Fig_scTWAS_multisource_bubble.{png,pdf} at 300 dpi
##         + _sourcedata.csv  (every plotted bubble)
##         + table/DE_key_results/scTWAS_bubble_summary.csv  (summary level)
## -----------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(patchwork)
})

PROJ <- "W:/VUMC_projects/scTWAS/MultiOmics_DE_support"
TAB  <- file.path(PROJ, "table", "scTWAS_347genes_multisource_DE_support.csv")
DIR_FIG <- file.path(PROJ, "figure")
DIR_KEY <- file.path(PROJ, "table", "DE_key_results")
dir.create(DIR_FIG, showWarnings = FALSE, recursive = TRUE)

LFC_LIM <- 2.5
Z_LIM   <- 13
INK <- "#0b0b0b"; INK2 <- "#52514e"; INK3 <- "#8a8880"
NA_COL <- "#b9b7b1"

## ---- input and gene selection ----------------------------------------------
ms <- fread(TAB)
## Row labels are the platform name alone; the contrast is the same in all
## three rows (normal epithelium versus carcinoma) and is stated once in the
## figure legend rather than repeated on every row.
##
## The first row is labelled GeoMx, not CosMx. The column it draws from is
## named CosMx_NOR_CORE, but that name is wrong and was verified so earlier:
## its values are identical to four decimals to the GeoMx DSP PanCK+ NOR->CORE
## rows of spatial_in/geomx_DE_compact.tsv (GTF2I 0.9750, TRIM28 1.2973,
## UBE2V1 1.0423). Figure5_legend.docx and the manuscript review were already
## corrected to GeoMx, so labelling this row CosMx would contradict them and
## would also collide with the genuine single-cell CosMx data used elsewhere.
## To revert, change this one string.
RES <- data.table(
  prefix = c("CosMx_NOR_CORE", "TCGA", "CPTAC"),
  row_label = c("GeoMx", "TCGA", "CPTAC"))

lab_cols <- paste0(RES$prefix, "_DEsupport")
## a character matrix has to be compared elementwise; %in% would flatten it
ms[, n_support := rowSums(as.matrix(.SD) == "Up" | as.matrix(.SD) == "Down",
                          na.rm = TRUE), .SDcols = lab_cols]
keep <- ms[n_support >= 1]
keep[, gene_class := fifelse(rep_Z > 0, "Oncogene-like  (scTWAS Z > 0)",
                             "Tumour-suppressor-like  (scTWAS Z < 0)")]
stopifnot(nrow(keep) == 109L, keep[, uniqueN(gene_name)] == nrow(keep))

## gene order: best-supported first, then by strength of the genetic signal
## setorder() takes column names, not expressions, so |Z| becomes a column
keep[, abs_Z := abs(rep_Z)]
setorder(keep, gene_class, -n_support, -abs_Z)
keep[, gene_name := factor(gene_name, levels = unique(gene_name))]

## ---- long form: one row per bubble -----------------------------------------
ROWS <- c("scTWAS Z score", RES$row_label)
mk <- function(dt, row_label, value, pval, kind) {
  data.table(gene_name = dt$gene_name, gene_class = dt$gene_class,
             n_support = dt$n_support, row_label = row_label,
             value = value, pval = pval, kind = kind)
}
bub <- rbindlist(c(
  list(mk(keep, "scTWAS Z score", keep$rep_Z, keep$rep_P, "Z")),
  lapply(seq_len(nrow(RES)), function(i)
    mk(keep, RES$row_label[i],
       keep[[paste0(RES$prefix[i], "_log2FC")]],
       keep[[paste0(RES$prefix[i], "_padj")]], "log2FC"))))
## ggplot draws the first factor level at the BOTTOM, so the display order
## (scTWAS, CosMx, TCGA, CPTAC top-to-bottom) is the reverse of the levels.
ROW_DISPLAY <- c("scTWAS Z score", RES$row_label)
stopifnot(setequal(ROW_DISPLAY, ROWS))
bub[, row_label := factor(row_label, levels = rev(ROW_DISPLAY))]

## significance bins. Values are binned on the p or padj of that row's own test.
BINS <- c(0, 1e-8, 1e-4, 0.01, 0.05, Inf)
BIN_LAB <- c("< 1e-8", "< 1e-4", "< 0.01", "< 0.05", "\u2265 0.05 (n.s.)")
bub[, sig_bin := cut(pval, breaks = BINS, labels = BIN_LAB, right = FALSE)]
bub[, sig_bin := factor(sig_bin, levels = BIN_LAB)]
SIZES <- setNames(c(4.4, 3.4, 2.5, 1.7, 0.9), BIN_LAB)

bub[, tested := !is.na(value)]
stopifnot(bub[tested == TRUE, !any(is.na(sig_bin))])

## ---- shared scales ---------------------------------------------------------
## Both ramps are diverging and centred on zero, which is the semantically
## meaningful centre for a Z score and for a log fold change alike.
DIV_LFC <- c("#2166AC", "#6BAED6", "#F7F7F7", "#F4A582", "#B2182B")
DIV_Z   <- c("#2166AC", "#7FB0D4", "#F7F7F7", "#EF9273", "#8F1D2E")

theme_bub <- function(show_x = FALSE) {
  theme_minimal(base_size = 8) +
    theme(
      panel.grid.major.y = element_line(colour = "grey93", linewidth = 0.25),
      panel.grid.major.x = element_line(colour = "grey96", linewidth = 0.18),
      panel.grid.minor   = element_blank(),
      axis.title         = element_blank(),
      axis.text.y        = element_text(size = 8.5, colour = INK),
      ## rotated 90 deg, so the limiting dimension is the per-gene slot width:
      ## 12.6 in / 2 panels / 57 genes = 0.11 in = about 7.9 pt. 6.5 clears it.
      axis.text.x        = if (show_x)
        element_text(size = 6.5, colour = INK, angle = 90, hjust = 1,
                     vjust = 0.5, face = "italic") else element_blank(),
      axis.ticks         = element_blank(),
      legend.text        = element_text(size = 7.5, colour = INK2),
      legend.title       = element_text(size = 7.5, colour = INK2),
      plot.margin        = margin(1, 3, 1, 3))
}

bubble_layer <- function(d) {
  list(
    ## untested combinations first, so a real bubble always draws over the cross
    geom_point(data = d[tested == FALSE], aes(x = gene_name, y = row_label),
               shape = 4, size = 0.7, stroke = 0.35, colour = NA_COL,
               inherit.aes = FALSE),
    geom_point(data = d[tested == TRUE],
               aes(x = gene_name, y = row_label, fill = value, size = sig_bin),
               shape = 21, colour = "grey30", stroke = 0.16))
}

panel_Z <- function(cls) {
  d <- bub[gene_class == cls & kind == "Z"]
  ggplot() + bubble_layer(d) +
    scale_fill_gradientn(colours = DIV_Z, limits = c(-Z_LIM, Z_LIM),
                         oob = scales::squish, name = "scTWAS Z",
                         guide = guide_colourbar(direction = "horizontal",
                                                 barwidth = unit(20, "mm"),
                                                 barheight = unit(2.4, "mm"),
                                                 title.position = "top",
                                                 ticks.colour = NA, order = 1)) +
    ## one size legend for the whole figure, carried by the log2FC sub-plots
    scale_size_manual(values = SIZES, drop = FALSE, guide = "none") +
    scale_x_discrete(drop = TRUE) +
    scale_y_discrete(limits = "scTWAS Z score") + theme_bub(FALSE)
}

panel_LFC <- function(cls) {
  d <- bub[gene_class == cls & kind == "log2FC"]
  ggplot() + bubble_layer(d) +
    scale_fill_gradientn(colours = DIV_LFC, limits = c(-LFC_LIM, LFC_LIM),
                         oob = scales::squish,
                         name = expression(log[2]~"fold change"),
                         guide = guide_colourbar(direction = "horizontal",
                                                 barwidth = unit(20, "mm"),
                                                 barheight = unit(2.4, "mm"),
                                                 title.position = "top",
                                                 ticks.colour = NA, order = 2)) +
    scale_size_manual(values = SIZES, drop = FALSE, name = "p or FDR",
                      guide = guide_legend(nrow = 1, order = 3,
                                           override.aes = list(fill = "grey75"),
                                           title.position = "top")) +
    scale_x_discrete(drop = TRUE) +
    ## explicit limits: layer-merged discrete levels are not order-stable
    scale_y_discrete(limits = rev(RES$row_label)) + theme_bub(TRUE)
}

CLS <- c("Oncogene-like  (scTWAS Z > 0)",
         "Tumour-suppressor-like  (scTWAS Z < 0)")
subt <- function(cls) {
  d <- keep[gene_class == cls]
  sprintf("%d genes \u00b7 %d supported by one resource, %d by two, %d by all three",
          nrow(d), d[n_support == 1, .N], d[n_support == 2, .N],
          d[n_support == 3, .N])
}

## Composition: one flat patchwork of four sub-plots, NOT two nested panels.
## Nesting each panel in its own plot_annotation forces wrap_elements(), which
## walls the guides off so every sub-plot keeps its own legend - four colourbars
## and four size keys. Flat + guides = "collect" merges the identical scales
## into one Z bar, one log2FC bar and one size key.
paZ <- panel_Z(CLS[1]) +
  labs(title = CLS[1], subtitle = subt(CLS[1]))
paL <- panel_LFC(CLS[1])
pbZ <- panel_Z(CLS[2]) +
  labs(title = CLS[2], subtitle = subt(CLS[2]))
pbL <- panel_LFC(CLS[2])

fig <- paZ / paL / pbZ / pbL +
  plot_layout(heights = c(1, 3, 1, 3), guides = "collect") +
  plot_annotation(
    title = paste0("Multi-resource support for scTWAS risk genes, split by the ",
                   "direction of the genetic effect"),
    caption = paste0(
      "All three fold-change rows are the same contrast, normal colon epithelium versus carcinoma: GeoMx DSP PanCK-positive compartment, TCGA bulk RNA and CPTAC bulk\n",
      "protein. Each bubble is one gene in one resource. Colour is the effect (scTWAS Z in the top row of each panel, log2 fold change in the other three) on a diverging scale\n",
      "centred at zero; colour is capped at |Z| = 13 and |log2FC| = 2.5, with more extreme values drawn at the end colour rather than dropped. Bubble size is a binned\n",
      "significance level for that row's own test. The bins are NOT comparable between rows: the GeoMx FDR is recomputed over the risk genes alone and bottoms out at\n",
      "4.2e-05, whereas TCGA and CPTAC FDRs are genome- and proteome-wide and reach 1e-119 and 1e-43. Compare sizes along a row, not down a column. A grey cross marks\n",
      "a gene the resource does not measure - CPTAC covers 57 of these 109 genes and GeoMx 93 - so an absent bubble is a missing measurement, not a null result.\n",
      "Genes shown are those with at least one direction-concordant Up or Down call (FDR < 0.05, |log2FC| > 0.5, sign matching the scTWAS Z) in GeoMx, TCGA or CPTAC,\n",
      "ordered by how many resources support them and then by |Z|. The three fold-change rows share one colour scale and so are directly comparable to each other."),
    theme = theme(
      plot.title   = element_text(size = 10.5, face = "bold", colour = INK,
                                  margin = margin(b = 4)),
      plot.caption = element_text(size = 5, colour = INK3, hjust = 0,
                                  lineheight = 1.3, margin = margin(t = 5)),
      plot.margin  = margin(6, 6, 4, 6))) &
  theme(legend.position = "top", legend.box = "horizontal",
        legend.justification = "center", legend.box.just = "bottom",
        legend.box.spacing = unit(1.5, "mm"),
        legend.spacing.x = unit(4, "mm"),
        legend.key.size = unit(3, "mm"),
        plot.title = element_text(size = 9.6, face = "bold", colour = INK,
                                  margin = margin(t = 4, b = 1)),
        plot.subtitle = element_text(size = 7.2, colour = INK2,
                                     margin = margin(b = 3)))

W <- 12.6; H <- 7.6
ggsave(file.path(DIR_FIG, "Fig_scTWAS_multisource_bubble.png"), fig,
       width = W, height = H, dpi = 300, bg = "white")
ggsave(file.path(DIR_FIG, "Fig_scTWAS_multisource_bubble.pdf"), fig,
       width = W, height = H, bg = "white")

## ---- source data -----------------------------------------------------------
sd <- bub[, .(gene_class, gene_name, n_support, row_label, kind,
              value, pval, sig_bin, tested)]
setorder(sd, gene_class, gene_name, row_label)
fwrite(sd, file.path(DIR_FIG, "Fig_scTWAS_multisource_bubble_sourcedata.csv"))

## ---- summary level ---------------------------------------------------------
## per class x resource: how many genes are measured, called Up, called Down
res_long <- rbindlist(lapply(seq_len(nrow(RES)), function(i) {
  p <- RES$prefix[i]
  data.table(gene_class = keep$gene_class, resource = RES$row_label[i],
             call = keep[[paste0(p, "_DEsupport")]],
             lfc = keep[[paste0(p, "_log2FC")]])
}))
summ <- res_long[, .(
  genes_measured = sum(!is.na(lfc)),
  supported_Up   = sum(call == "Up", na.rm = TRUE),
  supported_Down = sum(call == "Down", na.rm = TRUE),
  discordant     = sum(call == "DE_discordant", na.rm = TRUE),
  tested_ns      = sum(call == "NS", na.rm = TRUE),
  not_measured   = sum(is.na(lfc)),
  median_abs_lfc = round(median(abs(lfc), na.rm = TRUE), 3)),
  by = .(gene_class, resource)]
summ[, pct_supported_of_measured :=
       round(100 * (supported_Up + supported_Down) / genes_measured, 1)]

overall <- keep[, .(genes = .N,
                    one_resource = sum(n_support == 1),
                    two_resources = sum(n_support == 2),
                    three_resources = sum(n_support == 3),
                    median_abs_Z = round(median(abs(rep_Z)), 2),
                    max_abs_Z = round(max(abs(rep_Z)), 2)), by = gene_class]

fwrite(summ, file.path(DIR_KEY, "scTWAS_bubble_summary.csv"))
fwrite(overall, file.path(DIR_KEY, "scTWAS_bubble_summary_byclass.csv"))

cat("bubbles:", nrow(bub), "| genes:", nrow(keep), "\n")
print(overall)
print(summ)
