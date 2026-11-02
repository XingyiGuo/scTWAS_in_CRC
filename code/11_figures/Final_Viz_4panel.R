###############################################################################
# Final_Viz_4panel.R
#
# PURPOSE
#   Build ONE 6 x 6 inch, four-panel PDF summarising the scTWAS resource:
#     A  Integrated GSE266546 + HTAN-VUMC epithelial UMAP, reviewed cell types,
#        DRAWN HERE from the embedding's cell-level table.  The plotting code is
#        lifted from L:/scRNA-seq/VUMC_HTAN_GSE266546/src/05_annotation.R, which
#        renders this figure as
#          DimPlot(obj, reduction = "umap_harmony",
#                  group.by = "celltype_annot_reviewed", label = TRUE,
#                  raster = FALSE) + NoLegend()
#     B  WellcomeSanger healthy epithelium, seven-type UMAP, DRAWN HERE from the
#        embedding's cell-level table.  Plotting code lifted from
#        L:/scRNA-seq/WellcomeSanger/src/04_umap_epithelial_LauLabel.R, which
#        draws it as sub-panel b ("LauLabel annotation") of a three-panel figure:
#          ggplot(plot_df, aes(UMAP_1, UMAP_2, colour = LauLabel)) +
#            geom_point(size = 0.22, alpha = 0.85, show.legend = FALSE) +
#            geom_text(data = centroids, ...) + scale_colour_manual(...)
#
#        Panels A and B therefore go through ONE builder, umap_panel(), so their
#        point size, palette, label placement and aspect ratio agree by
#        construction instead of by eye.  Both keep a raster fallback:
#        PANEL_A_MODE / PANEL_B_MODE <- "raster" embeds the finished source
#        figure as a rasterised PDF page, which is what this figure shipped with
#        before, and is also what happens automatically when a panel's
#        cell-level table cannot be found.
#     C  Number of scATAC-seq peaks per cell type in normal colon (bar plot,
#        numbers typed in below - they are the raw data supplied for this figure)
#     D  Distribution of global R^2 over ALL scTWAS prediction models, one box
#        per cell type per expression resource, with the predictability
#        threshold R2_global = 0.01 drawn as a dashed reference line
#        (box plot, computed here from the *_model_summaries.csv files)
#
#   Layout is a 2 x 2 grid: A | B on row 1, C | D on row 2.
#
# HOW TO RUN
#   Rscript "W:/VUMC_projects/scTWAS/Viz/Final_Viz_4panel.R"
#   (no arguments; every path is absolute and listed in the CONFIG block below)
#
# R PACKAGES
#   ggplot2, patchwork, data.table, pdftools, png
#   pdftools bundles poppler, so panels A and B are rasterised WITHOUT needing
#   ghostscript or ImageMagick.  If pdftools is unavailable the script falls
#   back to the sibling .png export of each source figure (see read_fig()).
#
# OUTPUTS  (all written next to this script, in W:/VUMC_projects/scTWAS/Viz)
#   Final_Viz_4panel.pdf                    <- the 6 x 6 in deliverable
#   Final_Viz_4panel.png                    <- 600 dpi raster preview
#   Final_Viz_4panel_panelA_sourcedata.csv  <- the cells plotted in panel A
#   Final_Viz_4panel_panelB_sourcedata.csv  <- the cells plotted in panel B
#   Final_Viz_4panel_panelC_sourcedata.csv  <- peak counts as plotted
#   Final_Viz_4panel_panelD_sourcedata.csv  <- box statistics as plotted
#   Final_Viz_4panel_sessionInfo.txt        <- R/package versions of the run
#
# NOTES / CAVEATS worth knowing before editing
#   * Panels A and B are DRAWN HERE (see umap_panel), so the whole figure obeys
#     one type ladder (BASE_PT / TICK_PT) and one palette.  In the raster
#     fallback that is not true: a rasterised source page carries text sized for
#     its own canvas (8 x 6 in for A, 9.5 x 4.6 in for B), which then shrinks
#     when it is squeezed into a 6 x 6 in page.
#   * PANEL A NEEDS A COORDINATE TABLE, and the only copy of those coordinates
#     is inside a 5.8 GB Seurat object (PANEL_A_RDS).  extract_umap_sourcedata()
#     below writes the small table once - cell_id, UMAPHARMONY_1, UMAPHARMONY_2,
#     annotation - and every later run reads that CSV.  Loading the object needs
#     Seurat and roughly 12 GB of free RAM, so if the extraction cannot run the
#     script says so and falls back to the rasterised page.
#   * HOW MANY CELLS PANEL B SHOWS: 21,000 (3,000 per LauLabel) - that is the
#     only WellcomeSanger embedding that exists, and it is what both the source
#     figure and PANEL_B_CSV are built from.  In "native" mode the panel plots
#     every row of the first source-data file in the candidate list, so it
#     follows the embedding automatically.  To show all 358,207 retained epithelial
#     cells the embedding itself has to be recomputed on the cluster - a UMAP of
#     that size cannot run on this Windows box - with
#         WS_CELLS_PER_TYPE=all Rscript 04_umap_epithelial_LauLabel.R
#     (L:/scRNA-seq/WellcomeSanger/src/04_umap_epithelial_LauLabel.R, which
#     writes umap_WellcomeSanger_epithelial_clusters_LauLabel_sourcedata.csv).
#     That filename is already first in the candidate list, so once the file is
#     copied back into the umap folder this script picks it up with no edit, and
#     the message it prints reports the cells actually plotted per type.
#   * The cell-type glosses are not on the figure (STM - stem, TAC - transit-
#     amplifying, ABS - absorptive/enterocyte, CT - colonocyte, GOB - goblet,
#     EE - enteroendocrine, TUF - tuft); state them in the figure legend text.
#     Both UMAPs carry their type codes on the plot and panel C names all seven
#     on its x axis, so the composite still carries a colour key.
#   * Panel A is cropped (CROP_A, fractions of the source page, 0 to keep) to
#     drop its internal title "celltype_annot_reviewed".
#   * Colour harmonisation across panels A, B and C: all three are drawn with
#     CT_COLORS (the Okabe-Ito based LauLabel palette), which is the main reason
#     to draw A and B here rather than embed them.  The source figures do not
#     agree with each other - 05_annotation.R lets ggplot assign its default hue
#     wheel, 04_umap_epithelial_LauLabel.R uses its own LauLabel palette - and in
#     the raster fallback panel A's hues are therefore remapped IN THE PIXELS by
#     recolour_hues() (see the comment there for the verification).  That whole
#     mechanism is bypassed in native mode.
#   * Panel D summarises the models with R2_global > 0 (44% of all models are
#     negative, since R2_global is cross-validated), and marks the 0.01
#     predictability threshold with a dashed line and a labelled tick rather
#     than filtering on it.  y runs 0 to 0.10 via coord_cartesian(): every box
#     fits (largest q3 0.088, WellcomeSanger ABS), 11 of the 14 upper whiskers
#     are clipped (whisker range 0.043-0.207) and the un-drawn outliers reach
#     0.82.  Outliers are suppressed (outlier.shape = NA) because ~64,000
#     models would bury the boxes and bloat the PDF.  Full five-number
#     summaries, whisker limits and the per-group counts above the threshold go
#     to Final_Viz_4panel_panelD_sourcedata.csv.
#
# Author: prepared for the scTWAS manuscript figure set.
###############################################################################

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(data.table)
  library(grid)
  library(png)
})

## ===========================================================================
## CONFIG - every input and output path, in full
## ===========================================================================

OUT_DIR <- "W:/VUMC_projects/scTWAS/Viz"

# --- Panel A: integrated GSE266546 + HTAN-VUMC UMAP ------------------------
# Produced by L:/scRNA-seq/VUMC_HTAN_GSE266546/src/05_annotation.R
# (cluster-level dominant-vote annotation, Harmony on `batch`, res 0.8,
#  cluster 23 removed, ABS-heavy clusters re-clustered in 04_sub_cluster.R).
#
# PANEL_A_MODE, like PANEL_B_MODE below:
#   "native" - redraw the UMAP here from its cell-level table (PANEL_A_CSV), so
#              palette, point size and label style are this script's to set.
#   "raster" - embed the finished source page as a picture (PANEL_A_PDF), with
#              its hues remapped by recolour_hues().  This is also the automatic
#              fallback when no cell-level table is available.
PANEL_A_MODE <- "native"

# The source figure, kept for the raster fallback.
PANEL_A_PDF <- "L:/scRNA-seq/VUMC_HTAN_GSE266546/figure/05e_umap_reviewed_celltype_res0.8.rm23_sub.pdf"
PANEL_A_PNG <- "L:/scRNA-seq/VUMC_HTAN_GSE266546/figure/Fig1_umap_reviewed_celltype.png"  # fallback only

# Cell-level table for native mode: first file that exists wins.  Written by
# extract_umap_sourcedata() from PANEL_A_RDS, which is the object 05_annotation.R
# saves after the reviewed-annotation step (5.76 GB).  The metadata CSV that
# 05_annotation.R writes alongside it is NOT enough: Seurat keeps reduction
# embeddings outside @meta.data, so that file carries celltype_annot_reviewed
# but no coordinates.
PANEL_A_CSV <- c(
  "L:/scRNA-seq/VUMC_HTAN_GSE266546/table/umap_HTAN_GSE266546_reviewed_celltype_sourcedata.csv",
  "W:/VUMC_projects/scTWAS/Viz/Final_Viz_4panel_panelA_sourcedata.csv"
)
PANEL_A_RDS <- "L:/scRNA-seq/VUMC_HTAN_GSE266546/raw/rds2/integrated_seurat_final2_reviewed.res0.8.rm23_sub.rds"
PANEL_A_REDUCTION <- "umap_harmony"          # as in 05_annotation.R
PANEL_A_GROUP <- "celltype_annot_reviewed"   # as in 05_annotation.R (05e)

# Try the one-off extraction when PANEL_A_CSV is missing.  Set FALSE to go
# straight to the raster fallback without touching the 5.8 GB object.
EXTRACT_PANEL_A <- TRUE

# --- Panel B: WellcomeSanger healthy epithelium UMAP ----------------------
# Embedding produced by L:/scRNA-seq/WellcomeSanger/src/04_umap_epithelial_LauLabel.R
# from Healthy_TI_then_rectum_278samples.seurat.rds / Epithelial_expression_BPCELLS
# (358,207 retained epithelial cells, 278 donors): log-normalise 1e4 -> 2,000
# variable genes -> 50 PCs -> Harmony on individual_id -> k=15 -> UMAP.
# Panel B is drawn from the CELL-LEVEL source data of that run.  First file that
# exists wins, so an all-cells re-run is picked up automatically:
UMAP_DIR <- "L:/scRNA-seq/WellcomeSanger/Healthy_TI_then_rectum_278samples_outputs_LauLabel/umap"

# PANEL_B_MODE picks which of the two panel B builds is used:
#   "raster" (what this figure ships with) - embed the LEFT sub-panel of the
#            finished two-part figure as a picture, exactly as the source figure
#            drew it: its own bold white-outlined coloured type codes, its own
#            axis type.  Cropped by CROP_B.
#   "native" - redraw the same embedding from its cell-level table, which puts
#            the palette, type sizes and label style under this script's control
#            and plots every row of PANEL_B_CSV (so an all-cells re-run of the
#            embedding is picked up automatically).  Labels then follow
#            B_LABEL_* below.
# Either way the plotted cells are written to OUT_B, since both builds come from
# the same embedding.
PANEL_B_MODE <- "native"

# Sub-panel b ("LauLabel annotation") of the three-panel clustering figure that
# 04_umap_epithelial_LauLabel.R writes: 28,000 cells plotted (4,000 per LauLabel,
# seed 20260923) of the 358,207 retained, Louvain res 0.5 on 50 PCs, and the
# seven type codes annotated with their glosses in plain black text - the same
# annotation style as panel A.  Sub-panels a (clusters) and c (Final_labels),
# the sub-panel title, the two-column cluster key and the caption are cropped.
PANEL_B_PDF <- file.path(UMAP_DIR, "umap_WellcomeSanger_epithelial_clusters_LauLabel.pdf")
PANEL_B_PNG <- file.path(UMAP_DIR, "umap_WellcomeSanger_epithelial_clusters_LauLabel.png")  # fallback only

# panel B raster crop, c(top, right, bottom, left) as fractions of the source
# page, every edge measured as a blank band on the rendered page (9.5 x 4.6 in):
# the three UMAP panels occupy rows 0.084-0.648 and are separated by blank
# column bands at 0.3179-0.3551 (a|b) and 0.6456-0.6828 (b|c), so sub-panel b is
# the span between them; the sub-panel titles sit in rows 0.026-0.052 and the
# cluster key plus caption below row 0.687.  Re-measure if the source figure is
# re-rendered at a different size or with a different number of sub-panels.
CROP_B <- c(0.070, 0.336, 0.340, 0.336)

PANEL_B_CSV <- c(
  # the table behind PANEL_B_PDF: the clustered run of 04, one row per plotted
  # cell with seurat_cluster and cluster_majority_LauLabel
  file.path(UMAP_DIR, "umap_WellcomeSanger_epithelial_clusters_LauLabel_sourcedata.csv"),
  # all 358,207 retained epithelial cells (the uncapped embedding)
  file.path(UMAP_DIR, "umap_allcells_sourcedata.csv.gz"),
  # earlier 3,000-cells-per-type exports, if they are still around
  file.path(UMAP_DIR, "umap_WellcomeSanger_7types_sourcedata_v2.csv"),
  file.path(UMAP_DIR, "umap_WellcomeSanger_figure_sourcedata.csv")
)

# --- Panel D: scTWAS model summaries --------------------------------------
# One CSV per cell type per expression resource, columns:
#   gene, genename, pred.perf.r2, pred.perf.pval, pred.perf.qval,
#   n.snps.in.model, R2_global, r_squared_global, r.global
TWAS_ROOT   <- "W:/VUMC_projects/scTWAS/TWAS"
TWAS_SETS   <- c(                                   # display label -> folder / file prefix
  "VUMC_HTAN_GSE266546" = "VUMC_HTAN_GSE266546_NL",
  "WellcomeSanger"      = "WellcomeSanger_NL"
)
# file pattern: <TWAS_ROOT>/<prefix>/<prefix>_<celltype>_peakSNP_model_summaries.csv
R2_CUTOFF   <- 0.01                                  # "models with R2_global > 0.01"

# --- Outputs ---------------------------------------------------------------
OUT_PDF <- file.path(OUT_DIR, "Final_Viz_4panel.pdf")
OUT_PNG <- file.path(OUT_DIR, "Final_Viz_4panel.png")
OUT_A   <- file.path(OUT_DIR, "Final_Viz_4panel_panelA_sourcedata.csv")
OUT_B   <- file.path(OUT_DIR, "Final_Viz_4panel_panelB_sourcedata.csv")
OUT_C   <- file.path(OUT_DIR, "Final_Viz_4panel_panelC_sourcedata.csv")
OUT_D   <- file.path(OUT_DIR, "Final_Viz_4panel_panelD_sourcedata.csv")
OUT_SI  <- file.path(OUT_DIR, "Final_Viz_4panel_sessionInfo.txt")

## ===========================================================================
## STYLE - shared type ladder, cell-type order and palettes
## ===========================================================================

FIG_W <- 6      # inches, final page width
FIG_H <- 6      # inches, final page height
OUT_DPI <- 600  # effective resolution of the embedded rasters in panels A/B

BASE_PT  <- 7   # axis titles / legend text in panels C and D
TICK_PT  <- 6   # tick labels
VALUE_PT <- 5   # in-bar value labels (panel C); ~5 pt is the legibility floor
TAG_PT   <- 9   # panel letters A-D

# Crypt order (stem -> transit-amplifying -> absorptive -> colonocyte, then the
# secretory/rare types), the order used throughout this manuscript's figures.
CT_LEVELS <- c("STM", "TAC", "ABS", "CT", "GOB", "EE", "TUF")

# LauLabel palette, copied verbatim from
#   L:/scRNA-seq/WellcomeSanger/src/04_umap_epithelial_LauLabel.R
# so panel C is colour-consistent with panel B.
CT_COLORS <- c(
  STM = "#CC79A7",  # pink
  TAC = "#6A3D9A",  # purple
  ABS = "#0072B2",  # blue
  CT  = "#E69F00",  # amber
  GOB = "#009E73",  # green
  EE  = "#D55E00",  # vermillion
  TUF = "#56B4E9"   # sky blue
)

# Panel D encodes RESOURCE, not cell type, so it deliberately uses two hues
# that are not in CT_COLORS (NPG dark blue / NPG red).
SET_COLORS <- c("VUMC_HTAN_GSE266546" = "#3C5488",
                "WellcomeSanger"      = "#E64B35")

theme_panel <- function() {
  theme_classic(base_size = BASE_PT) +
    theme(
      axis.text       = element_text(size = TICK_PT, colour = "black"),
      axis.title      = element_text(size = BASE_PT, colour = "black"),
      axis.ticks      = element_line(linewidth = 0.25, colour = "black"),
      axis.ticks.length = unit(1.4, "pt"),          # outward ticks
      axis.line       = element_line(linewidth = 0.3, colour = "black"),
      legend.key.size = unit(6, "pt"),
      legend.text     = element_text(size = TICK_PT),
      legend.title    = element_blank(),
      legend.margin   = margin(0, 0, 0, 0),
      legend.box.spacing = unit(2, "pt"),
      legend.background  = element_blank(),          # frameless legend
      plot.margin     = margin(2, 3, 1, 2)
    )
}

## ===========================================================================
## HELPERS - embed an existing figure page as an aspect-preserving panel
## ===========================================================================

# Drop fully white margins so the embedded figure fills its cell.
# `bmp` is the raw array returned by pdftools::pdf_render_page(), dim = c(4, w, h)
# (channels R,G,B,A by x by y).  `tol` is the per-channel value at or above
# which a pixel counts as white; `pad` keeps a few pixels of breathing room.
trim_white_raw <- function(bmp, tol = 250L, pad = 4L) {
  d <- dim(bmp); w <- d[2]; h <- d[3]
  r <- matrix(as.integer(bmp[1, , ]), w, h)
  g <- matrix(as.integer(bmp[2, , ]), w, h)
  b <- matrix(as.integer(bmp[3, , ]), w, h)
  ink <- !(r >= tol & g >= tol & b >= tol)          # TRUE where something is drawn
  xs <- which(rowSums(ink) > 0)                     # occupied x (width) indices
  ys <- which(colSums(ink) > 0)                     # occupied y (height) indices
  if (!length(xs) || !length(ys)) return(bmp)       # blank page: leave as is
  x1 <- max(1L, min(xs) - pad); x2 <- min(w, max(xs) + pad)
  y1 <- max(1L, min(ys) - pad); y2 <- min(h, max(ys) + pad)
  bmp[, x1:x2, y1:y2, drop = FALSE]
}

# Add white margin (in pixels) on each side AFTER trimming: c(top, right, bottom, left).
# Panels A and B are trimmed flush to their ink, which puts the patchwork panel
# letter on top of the first tick label; a little padding gives the tag a home.
pad_raw <- function(bmp, pad_px = c(0, 0, 0, 0)) {
  if (all(pad_px == 0)) return(bmp)
  d <- dim(bmp); w <- d[2]; h <- d[3]
  out <- array(as.raw(255L), c(4, w + pad_px[2] + pad_px[4], h + pad_px[1] + pad_px[3]))
  out[, pad_px[4] + seq_len(w), pad_px[1] + seq_len(h)] <- bmp
  out
}

# Crop a fraction off each side BEFORE trimming: crop = c(top, right, bottom, left).
crop_raw <- function(bmp, crop = c(0, 0, 0, 0)) {
  if (all(crop == 0)) return(bmp)
  d <- dim(bmp); w <- d[2]; h <- d[3]
  x1 <- 1L + floor(crop[4] * w); x2 <- w - floor(crop[2] * w)
  y1 <- 1L + floor(crop[1] * h); y2 <- h - floor(crop[3] * h)
  bmp[, x1:x2, y1:y2, drop = FALSE]
}

# ---------------------------------------------------------------------------
# Remap the point colours of an already-rendered figure.
#
# Panel A's PDF was drawn with ggplot's default hue wheel.  scales::hue_pal()(7)
# is exactly #F8766D #C49A00 #53B400 #00C094 #00B6EB #A58AFF #FB61D7 and those
# seven hexes are precisely the seven colours found in the rendered page, in
# ggplot's alphabetical factor order ABS, CT, EE, GOB, STM, TAC, TUF.  The
# pixel counts corroborate that assignment independently of the ordering rule:
# at 100 dpi the page holds 39,274 #F8766D pixels (ABS, the dominant mass),
# 10,845 #00C094 (GOB, the large adjacent lobe), 4,450 #A58AFF (TAC), 2,446
# #C49A00 (CT, the thin spur), 1,669 #00B6EB (STM), 757 #FB61D7 (TUF) and 464
# #53B400 (EE) - the same size ranking the seven clusters have in the data.
#
# Each pixel is treated as a blend of white and one base hue, p = a*hue +
# (1-a)*white, so anti-aliased edges and light point overlaps are recoloured at
# their own opacity instead of leaving haloes.  A pixel is only touched when
# some source hue explains it to within `resid_tol` and with opacity >= a_min,
# which is what keeps black text, grey axes and white background untouched
# (their residual against every hue is > 0.18).
recolour_hues <- function(bmp, from, to, resid_tol = 0.10, a_min = 0.12) {
  stopifnot(length(from) == length(to))
  d <- dim(bmp)
  v <- 1 - matrix(as.integer(bmp[1:3, , ]), nrow = 3) / 255   # 3 x npix, 0 = white
  best_res <- rep(Inf, ncol(v)); best_a <- numeric(ncol(v)); best_k <- integer(ncol(v))
  U <- 1 - grDevices::col2rgb(from) / 255                     # 3 x K hue directions
  for (k in seq_along(from)) {
    u  <- U[, k]
    a  <- as.vector(crossprod(u, v)) / sum(u * u)             # least-squares opacity
    res <- sqrt(colSums((v - outer(u, a))^2))
    hit <- res < best_res
    best_res[hit] <- res[hit]; best_a[hit] <- a[hit]; best_k[hit] <- k
  }
  ok <- best_res < resid_tol & best_a >= a_min & best_a <= 1.25
  Tm <- 1 - grDevices::col2rgb(to) / 255
  newv <- Tm[, best_k[ok], drop = FALSE] * rep(best_a[ok], each = 3)
  v[, ok] <- newv
  out <- bmp
  out[1:3, , ] <- as.raw(round(pmin(pmax(1 - v, 0), 1) * 255))
  # Diagnostic: the per-hue pixel counts should keep the cluster size ranking
  # (ABS >> GOB > TAC > CT > STM > TUF > EE); if they do not, the from/to
  # mapping is wrong and the panel must not be used.
  message(sprintf("  recoloured %s px (%s)",
                  format(sum(ok), big.mark = ","),
                  paste(sprintf("%s:%s", names(to),
                                format(tabulate(best_k[ok], length(from)), big.mark = ",",
                                       trim = TRUE)), collapse = " ")))
  out
}

# Read a finished figure and return a native raster (integer matrix: one value
# per pixel, so it stays small in memory).  The page is rasterised at whatever
# dpi gives `target_w_in` inches of final width at OUT_DPI - i.e. no more pixels
# than the 6 x 6 in page can show, which keeps the PDF a few MB rather than tens.
read_fig <- function(pdf_path, png_path = NULL, target_w_in = 3,
                     crop = c(0, 0, 0, 0), trim = TRUE,
                     pad_in = c(0, 0, 0, 0), recolour = NULL) {
  tmp <- tempfile(fileext = ".png")
  if (requireNamespace("pdftools", quietly = TRUE) && file.exists(pdf_path)) {
    # dpi is set from the width that SURVIVES the crop, so a panel cropped out
    # of a wide page (panel B) is still rasterised at OUT_DPI of final size.
    src_w_in <- pdftools::pdf_pagesize(pdf_path)$width[1] / 72 * (1 - crop[2] - crop[4])
    dpi <- ceiling(target_w_in * OUT_DPI / src_w_in)
    bmp <- pdftools::pdf_render_page(pdf_path, page = 1, dpi = dpi, numeric = FALSE)
    bmp <- crop_raw(bmp, crop)
    if (trim) bmp <- trim_white_raw(bmp)
    if (!is.null(recolour)) bmp <- recolour_hues(bmp, recolour$from, recolour$to)
    # pad_in is expressed in FINAL figure inches, hence the OUT_DPI scaling
    bmp <- pad_raw(bmp, round(pad_in * OUT_DPI))
    png::writePNG(bmp, tmp)
    message(sprintf("  embedded %s at %d dpi -> %d x %d px",
                    basename(pdf_path), dpi, dim(bmp)[2], dim(bmp)[3]))
  } else {
    # Fallback: the sibling PNG export of the same figure.
    stopifnot(!is.null(png_path), file.exists(png_path))
    a <- png::readPNG(png_path)                      # dim = c(h, w, channels)
    if (length(dim(a)) == 2) a <- array(a, c(dim(a), 1))
    bmp <- aperm(a, c(3, 2, 1))                      # -> c(channels, w, h)
    bmp <- as.raw(round(bmp * 255)); dim(bmp) <- c(dim(a)[3], dim(a)[2], dim(a)[1])
    if (dim(bmp)[1] == 3) bmp <- abind3(bmp)          # add an opaque alpha plane
    bmp <- crop_raw(bmp, crop)
    if (trim) bmp <- trim_white_raw(bmp)
    png::writePNG(bmp, tmp)
    message(sprintf("  embedded PNG fallback %s -> %d x %d px",
                    basename(png_path), dim(bmp)[2], dim(bmp)[3]))
  }
  png::readPNG(tmp, native = TRUE)
}

# helper for the RGB -> RGBA case in the PNG fallback
abind3 <- function(bmp) {
  d <- dim(bmp)
  out <- array(as.raw(255), c(4, d[2], d[3]))
  out[1:3, , ] <- bmp
  out
}

# grid::rasterGrob with width/height left NULL preserves the image's aspect
# ratio and centres it inside whatever cell patchwork hands it, so panels A and
# B can never be stretched by the layout.
as_panel <- function(nr) {
  wrap_elements(full = grid::rasterGrob(nr, interpolate = TRUE))
}

## ===========================================================================
## PANEL A + B - two UMAPs, one builder
## ===========================================================================
## Both panels are the same kind of object: a Harmony-corrected UMAP of colonic
## epithelium, coloured by the seven crypt cell types, with the type code typed
## on each cluster.  The two source scripts draw them differently (Seurat's
## DimPlot with ggplot's default hue wheel in 05_annotation.R; a hand-built
## ggplot with the LauLabel palette in 04_umap_epithelial_LauLabel.R), so the
## plotting code is consolidated here into umap_panel() and both panels are fed
## through it.  What each source script contributes is its DATA and its
## definition of the label column; nothing else is inherited.

# Point transparency, taken from 04_umap_epithelial_LauLabel.R's geom_point.
UMAP_ALPHA <- 0.85

# Type codes are annotated in panel A's style: plain (not bold) black text
# sitting directly on the coloured cluster, one label per cell type, no halo and
# no legend.  UMAP_LABEL_HALO <- TRUE restores the white-outlined coloured
# labels the WellcomeSanger source figure used (the halo is the same text drawn
# 8 times at a small offset, which is how you outline text without an extra
# package).
UMAP_LABEL_HALO   <- FALSE
UMAP_LABEL_COLOUR <- "black"   # set to NA to colour each label by its cell type
UMAP_LABEL_FACE   <- "plain"

# Label nudges in DATA units, per panel, applied after the anchor is found.
A_NUDGE <- list()
B_NUDGE <- list(STM = c(2.2, 1.6), TAC = c(-0.6, -1.2))

## --- one-off: turn a Seurat object into a small coordinate table -------------
# The embeddings behind panel A exist only inside a 5.76 GB Seurat object, so
# they are extracted once into a CSV and every later run reads that instead.
# Needs Seurat and roughly 12 GB of free RAM; returns the path on success and
# NULL on any failure, so the caller can fall back to the rasterised page.
extract_umap_sourcedata <- function(rds, reduction, group, out_csv) {
  if (!file.exists(rds)) {
    message("  no Seurat object at ", rds)
    return(NULL)
  }
  if (!requireNamespace("SeuratObject", quietly = TRUE) &&
      !requireNamespace("Seurat", quietly = TRUE)) {
    message("  Seurat/SeuratObject is not installed, so the embeddings cannot be read")
    return(NULL)
  }
  message("  reading ", basename(rds), " (",
          format(round(file.size(rds) / 1e9, 2)), " GB) - this is the slow step")
  ok <- try({
    obj <- readRDS(rds)
    emb <- as.data.frame(SeuratObject::Embeddings(obj, reduction = reduction)[, 1:2])
    names(emb) <- c("UMAPHARMONY_1", "UMAPHARMONY_2")
    md <- obj@meta.data
    stopifnot(group %in% names(md), nrow(md) == nrow(emb))
    out <- data.table(cell_id = rownames(emb),
                      UMAPHARMONY_1 = emb$UMAPHARMONY_1,
                      UMAPHARMONY_2 = emb$UMAPHARMONY_2,
                      annotation = as.character(md[[group]]))
    dir.create(dirname(out_csv), recursive = TRUE, showWarnings = FALSE)
    fwrite(out, out_csv)
    rm(obj, emb, md, out); invisible(gc())
  }, silent = TRUE)
  if (inherits(ok, "try-error")) {
    message("  extraction failed: ", conditionMessage(attr(ok, "condition")))
    return(NULL)
  }
  message("  wrote ", out_csv)
  out_csv
}

## --- read a cell-level UMAP table and normalise its column names ------------
# The two runs name their columns differently (UMAPHARMONY_1/2 + annotation in
# the panel A extract, UMAP_1/2 + LauLabel in 04's own export), so take the
# first coordinate pair and the first label column that are present.
read_umap_csv <- function(path) {
  d <- if (grepl("\\.gz$", path)) as.data.table(read.csv(gzfile(path))) else fread(path)
  xy <- Filter(function(p) all(p %in% names(d)),
               list(c("UMAPHARMONY_1", "UMAPHARMONY_2"), c("UMAP7_1", "UMAP7_2"),
                    c("UMAP_1", "UMAP_2"), c("umapharmony_1", "umapharmony_2"),
                    c("umap_1", "umap_2")))
  if (!length(xy)) stop(path, ": no UMAP coordinate columns found")
  lab <- intersect(c("cell_type", "annotation", "LauLabel", "celltype_annot_reviewed",
                     "cluster_majority_LauLabel", "Final_labels"), names(d))
  if (!length(lab)) stop(path, ": no cell-type column found")
  out <- d[, .(cell_id = if ("cell_id" %in% names(d)) as.character(get("cell_id"))
                         else if ("V1" %in% names(d)) as.character(get("V1"))
                         else as.character(seq_len(.N)),
               x = as.numeric(get(xy[[1]][1])), y = as.numeric(get(xy[[1]][2])),
               cell_type = as.character(get(lab[1])))]
  dropped <- setdiff(unique(out$cell_type), CT_LEVELS)
  if (length(dropped)) {
    n_drop <- out[cell_type %in% dropped, .N]
    message(sprintf("  %s: dropping %s cell(s) outside the seven crypt types (%s)",
                    basename(path), format(n_drop, big.mark = ","),
                    paste(head(dropped, 8), collapse = ", ")))
  }
  out <- out[is.finite(x) & is.finite(y) & cell_type %in% CT_LEVELS]
  out[, cell_type := factor(cell_type, levels = CT_LEVELS)]
  message(sprintf("  %s: %s cells (%s), label column '%s'",
                  basename(path), format(nrow(out), big.mark = ","),
                  paste(sprintf("%s %s", levels(out$cell_type),
                                format(as.vector(table(out$cell_type)), big.mark = ",",
                                       trim = TRUE)), collapse = " / "),
                  lab[1]))
  out[]
}

## --- the shared panel builder ------------------------------------------------
umap_panel <- function(dt, nudge = list(), pt = NULL, seed = 20260923) {
  stopifnot(all(c("x", "y", "cell_type") %in% names(dt)))
  dt <- copy(dt)

  # Draw order is randomised (fixed seed) so no cell type is systematically
  # painted over another - in panel A ABS alone would otherwise cover the crypt
  # types it overlaps.
  set.seed(seed)
  dt <- dt[sample(.N)]

  # Point size scaled to the number of cells: dense enough to read the clusters,
  # sparse enough not to fuse into solid blocks.
  if (is.null(pt)) pt <- if (nrow(dt) > 2e5) 0.035 else
                         if (nrow(dt) > 5e4) 0.08 else 0.22

  # The label anchor is NOT the per-type median: several types have two or more
  # lobes (EE is a cloud of islands, STM sits either side of TAC) and the median
  # then lands in the gap between them, which is how labels end up on top of
  # each other.  Instead each label goes to the type's DENSEST cell: bin the
  # embedding 40 x 40, take the fullest bin for that type, and use the median of
  # the cells inside it.  `nudge` moves individual labels afterwards, in data
  # units, when two dense lobes are genuinely adjacent.
  dt[, `:=`(bx = cut(x, 40, labels = FALSE), by = cut(y, 40, labels = FALSE))]
  topbin <- dt[, .N, by = .(cell_type, bx, by)][order(-N), .SD[1], by = cell_type]
  lab <- merge(dt, topbin[, .(cell_type, bx, by)], by = c("cell_type", "bx", "by"))[
    , .(x = median(x), y = median(y)), by = cell_type]
  for (nm in names(nudge)) lab[cell_type == nm,
    `:=`(x = x + nudge[[nm]][1], y = y + nudge[[nm]][2])]
  halo <- do.call(rbind, lapply(seq(0, 7), function(i)
    transform(lab, dx = 0.12 * cos(i * pi / 4), dy = 0.12 * sin(i * pi / 4))))

  lab_layer <- if (is.na(UMAP_LABEL_COLOUR)) {
    geom_text(data = lab, aes(x, y, label = cell_type, colour = cell_type),
              inherit.aes = FALSE, fontface = UMAP_LABEL_FACE, size = BASE_PT / .pt,
              show.legend = FALSE)
  } else {
    geom_text(data = lab, aes(x, y, label = cell_type), inherit.aes = FALSE,
              colour = UMAP_LABEL_COLOUR, fontface = UMAP_LABEL_FACE,
              size = BASE_PT / .pt)
  }

  ggplot(dt, aes(x, y, colour = cell_type)) +
    geom_point(size = pt, alpha = UMAP_ALPHA, stroke = 0, shape = 16,
               show.legend = FALSE) +
    {if (UMAP_LABEL_HALO)
       geom_text(data = halo, aes(x + dx, y + dy, label = cell_type),
                 inherit.aes = FALSE, colour = "white", fontface = "bold",
                 size = BASE_PT / .pt)} +
    lab_layer +
    scale_colour_manual(values = CT_COLORS, guide = "none") +
    coord_fixed(ratio = 1) +                 # a UMAP must not be stretched
    labs(x = "UMAP 1", y = "UMAP 2") +
    theme_panel()
}

## --- panel A -----------------------------------------------------------------
message("Panel A: integrated GSE266546 + HTAN-VUMC epithelial UMAP")

# crop the internal figure title ("celltype_annot_reviewed") off the top of A,
# and remap its default ggplot hues - both used by the raster fallback only.
# ggplot assigns hue_pal()(7) in alphabetical factor order, hence this mapping.
CROP_A <- c(0.085, 0, 0, 0)
RECOLOR_A <- TRUE
A_FROM <- c(ABS = "#F8766D", CT = "#C49A00", EE = "#53B400", GOB = "#00C094",
            STM = "#00B6EB", TAC = "#A58AFF", TUF = "#FB61D7")
A_MAP  <- if (RECOLOR_A) list(from = A_FROM, to = CT_COLORS[names(A_FROM)]) else NULL

a_file <- PANEL_A_CSV[file.exists(PANEL_A_CSV)][1]
if (is.na(a_file) && PANEL_A_MODE == "native" && EXTRACT_PANEL_A) {
  message("  no coordinate table yet; extracting one from the Seurat object")
  a_file <- extract_umap_sourcedata(PANEL_A_RDS, PANEL_A_REDUCTION, PANEL_A_GROUP,
                                    PANEL_A_CSV[length(PANEL_A_CSV)])
  if (is.null(a_file)) a_file <- NA_character_
}

if (PANEL_A_MODE == "native" && !is.na(a_file)) {
  ad <- read_umap_csv(a_file)
  fwrite(ad[, .(cell_id, UMAPHARMONY_1 = x, UMAPHARMONY_2 = y, cell_type)], OUT_A)
  pA <- umap_panel(ad, nudge = A_NUDGE)
} else {
  if (PANEL_A_MODE == "native")
    message("  falling back to the rasterised source page (no coordinate table)")
  pA <- as_panel(read_fig(PANEL_A_PDF, PANEL_A_PNG, target_w_in = 2.95, crop = CROP_A,
                          pad_in = c(0.07, 0, 0, 0.20),    # room for the "A" tag
                          recolour = A_MAP))
}

## --- panel B -----------------------------------------------------------------
message("Panel B: WellcomeSanger healthy epithelium UMAP")

b_file <- PANEL_B_CSV[file.exists(PANEL_B_CSV)][1]
if (is.na(b_file) && PANEL_B_MODE == "native") {
  message("  falling back to the rasterised source page (no coordinate table)")
  PANEL_B_MODE <- "raster"
}

if (PANEL_B_MODE == "native") {

  bd <- read_umap_csv(b_file)
  fwrite(bd[, .(cell_id, UMAP_1 = x, UMAP_2 = y, cell_type)], OUT_B)
  pB <- umap_panel(bd, nudge = B_NUDGE)

} else {

  # Panel B's source data is the embedding's own cell table, so it is copied
  # beside the figure rather than re-derived (358,207 rows: re-parsing and
  # re-writing it on every run costs a minute and gains nothing).
  if (!is.na(b_file)) {
    OUT_B <- paste0(sub("\\.csv$", "", OUT_B), if (grepl("\\.gz$", b_file)) ".csv.gz" else ".csv")
    file.copy(b_file, OUT_B, overwrite = TRUE)
    message("  source data: ", basename(b_file), " -> ", basename(OUT_B))
  }

  pB <- as_panel(read_fig(PANEL_B_PDF, PANEL_B_PNG, target_w_in = 2.95,
                          crop = CROP_B,
                          pad_in = c(0.07, 0, 0, 0.20)))   # room for the "B" tag

}   # end PANEL_B_MODE branch

## ===========================================================================
## PANEL C - scATAC peaks per cell type (normal colon)
## ===========================================================================
# Raw data as supplied for this figure (tissue, cell type, number of peaks).
# Keep this block as the single source of truth for panel C; it is written out
# verbatim to Final_Viz_4panel_panelC_sourcedata.csv.

peaks <- data.table(
  tissue    = "Normal colon",
  cell_type = c("ABS",    "CT",     "EE",    "GOB",    "STM",    "TAC",    "TUF"),
  n_peaks   = c( 150086L,  152781L,  68521L,  152548L,  145960L,  149672L,  55454L)
)
peaks[, cell_type := factor(cell_type, levels = CT_LEVELS)]
setorder(peaks, cell_type)
fwrite(peaks, OUT_C)

pC <- ggplot(peaks, aes(cell_type, n_peaks / 1000, fill = cell_type)) +
  geom_col(width = 0.72) +
  geom_text(aes(label = sprintf("%.1f", n_peaks / 1000)),
            vjust = -0.35, size = VALUE_PT / .pt, colour = "black") +
  scale_fill_manual(values = CT_COLORS, guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.13))) +
  labs(x = "Normal colon", y = expression("scATAC peaks (" * 10^3 * ")")) +
  theme_panel()

## ===========================================================================
## PANEL D - global R^2 of ALL scTWAS models, threshold 0.01 marked
## ===========================================================================
message("Panel D: reading model summaries")

read_one <- function(set_label, prefix, ct) {
  f <- file.path(TWAS_ROOT, prefix,
                 sprintf("%s_%s_peakSNP_model_summaries.csv", prefix, ct))
  if (!file.exists(f)) stop("missing model summary file: ", f)
  d <- fread(f, select = c("gene", "genename", "n.snps.in.model", "R2_global"))
  d[, `:=`(cell_type = ct, resource = set_label)]
  d[]
}

mods <- rbindlist(lapply(CT_LEVELS, function(ct)
  rbindlist(Map(read_one, names(TWAS_SETS), TWAS_SETS, ct))))
mods[, cell_type := factor(cell_type, levels = CT_LEVELS)]
mods[, resource  := factor(resource,  levels = names(TWAS_SETS))]

# R2_global is a cross-validated statistic, so 44% of all models are negative -
# the model predicts worse than the mean of the expression it is fitting.  Those
# are dropped here: the boxes summarise the models with R2_global > 0, and the
# 0.01 predictability threshold is drawn as a dashed reference line rather than
# used as a filter, so the panel shows what fraction of the positive-R^2 mass
# sits above it.
allm <- mods[is.finite(R2_global) & R2_global > 0]

# Source data = exactly the statistics the boxes draw (Tukey five-number
# summary plus 1.5 x IQR whiskers), alongside the counts above the threshold.
statD <- allm[, {
  q <- quantile(R2_global, c(0.25, 0.5, 0.75), names = FALSE)
  iqr <- q[3] - q[1]
  .(n_models        = .N,
    n_above_0.01    = sum(R2_global > R2_CUTOFF),
    min             = min(R2_global),
    q1              = q[1],
    median          = q[2],
    q3              = q[3],
    max             = max(R2_global),
    whisker_lower   = min(R2_global[R2_global >= q[1] - 1.5 * iqr]),
    whisker_upper   = max(R2_global[R2_global <= q[3] + 1.5 * iqr]))
}, by = .(cell_type, resource)]
statD[, pct_above_0.01 := round(100 * `n_above_0.01` / n_models, 1)]
setcolorder(statD, c("cell_type", "resource", "n_models",
                     "n_above_0.01", "pct_above_0.01"))
setorder(statD, cell_type, resource)
fwrite(statD, OUT_D)

# Display window.  coord_cartesian() clips the drawing without dropping rows,
# so the quartiles always come from every model with R2_global > 0.  Every box
# is still fully inside this window (the largest q3 is 0.088, WellcomeSanger
# ABS), but 11 of the 14 upper whiskers now run past the cap (whisker range
# 0.043-0.207) and the un-drawn outliers reach 0.82.
Y_LO <- 0
Y_HI <- 0.10

pD <- ggplot(allm, aes(cell_type, R2_global, fill = resource)) +
  geom_hline(yintercept = R2_CUTOFF, linetype = "22", linewidth = 0.3,
             colour = "grey35") +
  geom_boxplot(outlier.shape = NA, linewidth = 0.22, width = 0.68,
               position = position_dodge(width = 0.76), colour = "grey20") +
  # no in-panel label for the threshold: 0.01 is a labelled tick on the y axis
  scale_fill_manual(values = SET_COLORS) +
  scale_y_continuous(breaks = c(0, 0.01, 0.025, 0.05, 0.075, 0.10),
                     labels = c("0", "0.01", "0.025", "0.05", "0.075", "0.10"),
                     expand = expansion(mult = c(0, 0.02))) +
  # clip = "on" (the default) matters here: with clip = "off" the whiskers that
  # run past Y_HI are drawn outside the panel and land on top of panel B.
  coord_cartesian(ylim = c(Y_LO, Y_HI)) +
  labs(x = NULL, y = expression("global " * italic(R)^2)) +
  theme_panel() +
  theme(legend.position = "top", legend.direction = "horizontal",
        legend.text = element_text(size = TICK_PT - 0.5))

## ===========================================================================
## ASSEMBLE - 6 x 6 in page, 2 x 2 grid: A | B over C | D
## ===========================================================================
#   row 1  2.50 in   A (2.95 in wide, aspect 1.45:1 -> 2.04 in tall)
#                    B (2.95 in wide, aspect 1.31:1 -> 2.26 in tall)
#   row 2  3.50 in   C | D, each 2.95 in wide
# The two embedded rasters keep their own aspect ratios (see as_panel), so they
# are centred in their cells rather than stretched to match each other; panel A
# therefore sits about 0.2 in shorter than panel B.  Row 1 is deliberately not
# taller than 2.5 in: nothing is gained by it, since at 2.95 in wide neither
# raster can grow past its width-limited height.
#
# IMPORTANT: patchwork fills the areas of `design` in the order the plots are
# ADDED (1st plot -> area "A", 2nd -> "B", ...), and plot_annotation() letters
# them in that same order, so keep the pA, pB, pC, pD order.

design <- "
AB
CD
"

fig <- pA + pB + pC + pD +
  plot_layout(design = design,
              widths  = c(1, 1),
              heights = c(2.50, 3.50)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = TAG_PT, face = "bold", hjust = 0, vjust = 1))

message("writing ", OUT_PDF)
ggsave(OUT_PDF, fig, width = FIG_W, height = FIG_H, units = "in",
       device = grDevices::cairo_pdf)
ggsave(OUT_PNG, fig, width = FIG_W, height = FIG_H, units = "in", dpi = OUT_DPI)

writeLines(capture.output(sessionInfo()), OUT_SI)

message("done:")
message("  ", OUT_PDF)
message("  ", OUT_PNG)
message("  ", OUT_B)
message("  ", OUT_C)
message("  ", OUT_D)
###############################################################################
