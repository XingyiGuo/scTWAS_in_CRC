## 09_fig_NOR_CORE_highlight_4panel.R -----------------------------------------
##
## One row, four columns: the CosMx slide's cell composition, then three scTWAS
## risk genes up in the normal -> carcinoma transition, each annotated with its
## scTWAS statistics and every independent evidence layer that supports it.
##
##   a  cell composition   (source: fig16 panel a)
##   b  GTF2I              (source: fig18 panel d)
##   c  TRIM28             (source: fig16 panel d)
##   d  UBE2V1             (source: fig18 panel e)
##
## Inputs, written by 08_prep_cosmx_highlight_panels.py
##   cosmx_panelA_composition.csv.gz   x, y, group
##   cosmx_panel_expr.csv.gz           gene, x, y, expr, expressing
##   panel_stats.csv                   per-gene scTWAS / CosMx / other layers
##
## Two things about this figure that are easy to get wrong
##
## 1. Panels b and d come from a different cell subset than a and c. Both
##    archives sample the same slide, but each keeps every cell expressing one
##    of ITS genes plus a random context sample, so the point clouds differ in
##    density and outline. That is why each panel carries its own cell count
##    rather than one shared count in the title.
##
## 2. The layer labelled GeoMx is the compartment-level DE that the scTWAS
##    tables call "CosMx_NOR_CORE". It is GeoMx DSP, PanCK+ compartment, and it
##    is NOT independent of that column - it is the same measurement. The
##    independent line of evidence is the per-cell CosMx delta, which is a
##    paired carcinoma-minus-normal-epithelium difference over 7 patients with
##    a nominal p (smallest attainable 0.016 on a 7-pair sign test). The
##    annotation names both explicitly so the two are never read as one.
##
## Output: figure/Fig_NOR_CORE_highlight_4panel.{png,pdf} + _sourcedata.csv
## -----------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(patchwork); library(grid)
})

DIR_FIG <- "figure"
dir.create(DIR_FIG, showWarnings = FALSE, recursive = TRUE)

## ---- why the point clouds are rasterised by hand ---------------------------
## The four panels carry ~1.2 million cells between them. Drawn as points, the
## PDF would hold that many vector circles and become unopenable, and the PNG
## would take minutes. ggrastr is the usual answer but cannot be used here: it
## is absent from the shared R library, and a CRAN install into the workspace
## cannot be loaded because this sandbox refuses to dyn.load a shared library
## from a writable path.
##
## So each cloud is reduced to a pixel grid in base R before plotting. This is
## not an approximation of the original figures - it reproduces them. Their
## panels overplot points in ascending value order, so the visible colour at
## any pixel is the LARGEST value landing there, which is exactly what taking
## the per-pixel maximum computes. Panel a resolves ties the same way, by
## drawing group priority order.
##
## The consequence to be aware of: point SIZE no longer varies. In the source
## panels the expressing cells are drawn slightly larger than the grey context
## cells (0.34 vs 0.10), which dilates sparse positive signal. Here one cell is
## one pixel, so sparse expression reads slightly thinner. PX_H sets the
## resolution; raising it makes cells smaller, not more numerous.
PX_H <- 820L

rasterise_cloud <- function(x, y, value = NULL, ramp = NULL, vmax = NULL,
                            base_col = NULL, px_h = PX_H) {
  ## Return a matrix of hex colours laid out for annotation_raster, plus the
  ## data-space extent it should be drawn into. NA cells stay transparent.
  xr <- range(x); yr <- range(y)
  asp <- diff(xr) / diff(yr)
  ny <- px_h
  nx <- max(1L, as.integer(round(px_h * asp)))          # square pixels
  ix <- pmin(nx, pmax(1L, as.integer(floor((x - xr[1]) / diff(xr) * (nx - 1))) + 1L))
  iy <- pmin(ny, pmax(1L, as.integer(floor((y - yr[1]) / diff(yr) * (ny - 1))) + 1L))
  cell <- (iy - 1L) * nx + ix                           # flat pixel index

  if (is.null(value)) {
    ## categorical: value carried in base_col as a per-point colour, and the
    ## caller has pre-sorted so the last write wins
    px <- rep(NA_character_, nx * ny)
    px[cell] <- base_col
  } else {
    ## per-pixel maximum, computed without a full sort
    mx <- rep(NA_real_, nx * ny)
    o <- order(value)                                   # ascending: last wins
    mx[cell[o]] <- value[o]
    px <- rep(NA_character_, nx * ny)
    seen <- !is.na(mx)
    px[seen] <- base_col                                # any cell here at all
    hot <- seen & mx > 0
    if (any(hot)) {
      v <- pmin(mx[hot] / vmax, 1)
      px[hot] <- scales::gradient_n_pal(ramp)(v)
    }
  }
  m <- matrix(px, nrow = ny, ncol = nx, byrow = TRUE)
  ## annotation_raster puts row 1 at the TOP, while iy increases upward
  list(raster = m[ny:1L, , drop = FALSE], xr = xr, yr = yr,
       n_px = sum(!is.na(px)))
}

## ---- palettes, taken from the source figures so panels match ---------------
GCOL <- c("Carcinoma"          = "#e34948",
          "Normal epithelium"  = "#1baf7a",
          "Oncofetal"          = "#eda100",
          "Fibroblast / stroma"= "#8a8880",
          "Immune"             = "#9ec5f4",
          "Other"              = "#d8d7d2")
G_ORDER <- c("Other", "Immune", "Fibroblast / stroma", "Oncofetal",
             "Normal epithelium", "Carcinoma")
EXP_RAMP <- c("#e9e7e1", "#ffd9a0", "#f59b42", "#e34948", "#8f1d2e")
GREY_OFF <- "#eceae5"     # non-expressing cells
INK <- "#0b0b0b"; INK2 <- "#52514e"; INK3 <- "#8a8880"

theme_map <- function() {
  theme_void(base_size = 8) +
    ## the panel letter is drawn by patchwork at the top-left of this same
    ## corner, so the title block is indented past it rather than fighting it
    theme(plot.title    = element_text(size = 8.4, face = "bold", colour = INK,
                                       hjust = 0, margin = margin(b = 1, l = 11)),
          plot.subtitle = element_text(size = 5.1, colour = INK2, hjust = 0,
                                       lineheight = 1.25,
                                       margin = margin(b = 2, l = 11)),
          plot.margin   = margin(11, 3, 2, 3),
          legend.title  = element_text(size = 6, colour = INK2),
          legend.text   = element_text(size = 5.8, colour = INK2))
}

## ---- inputs ----------------------------------------------------------------
## fread() only reads .gz when R.utils is installed, which is not guaranteed
## here, so fall back to base R's native gzip connection.
read_any <- function(path) {
  if (grepl("\\.gz$", path) && !requireNamespace("R.utils", quietly = TRUE)) {
    as.data.table(utils::read.csv(gzfile(path), stringsAsFactors = FALSE))
  } else {
    fread(path)
  }
}
comp <- read_any("cosmx_panelA_composition.csv.gz")
expr <- read_any("cosmx_panel_expr.csv.gz")
st   <- read_any("panel_stats.csv")
setkey(st, gene)

GENES <- c("GTF2I", "TRIM28", "UBE2V1")
stopifnot(all(GENES %in% expr$gene), all(GENES %in% st$gene))

## ---- number formatting -----------------------------------------------------
## p values span 1e-20 to 0.05 here, so a fixed number of decimals is useless;
## anything below 1e-3 goes to scientific notation.
fmt_p <- function(p) ifelse(is.na(p), "NA",
                     ifelse(p < 1e-3, formatC(p, format = "e", digits = 1),
                            formatC(p, format = "f", digits = 3)))
fmt_sg <- function(v, d = 2) ifelse(is.na(v), "NA",
                                    sprintf("%+.*f", d, v))     # keep the sign

## One subtitle per gene panel. Order is deliberate: the risk statistic first
## (why the gene is a candidate), then the measurement this panel actually
## shows, then the compartment DE it replicates, then the remaining layers on
## one line. Only supported layers are listed; an unsupported one would invite
## reading a null as evidence.
subtitle_for <- function(g) {
  r <- st[g]
  others <- c()
  for (nm in c("Visium_M1", "Visium_M3", "TCGA_bulk_RNA", "CPTAC_bulk_protein")) {
    if (!is.na(r[[paste0(nm, "_label")]]) &&
        r[[paste0(nm, "_label")]] %in% c("Up", "Down")) {
      pretty <- c(Visium_M1 = "Visium M1", Visium_M3 = "Visium M3",
                  TCGA_bulk_RNA = "TCGA RNA",
                  CPTAC_bulk_protein = "CPTAC protein")[nm]
      others <- c(others, sprintf("%s %s", pretty,
                                  fmt_sg(r[[paste0(nm, "_log2FC")]])))
    }
  }
  paste0(
    sprintf("scTWAS Z %s, P %s  (%s", fmt_sg(r$rep_Z), fmt_p(r$rep_P), r$rep_model),
    if (r$n_models > 1) sprintf(" \u00b7 %d models, %s)", r$n_models,
                                tolower(r$direction_status)) else ")",
    "\n",
    sprintf("CosMx per-cell \u0394 %s, P %s (n = %d patients)",
            fmt_sg(r$cosmx_delta, 3), fmt_p(r$cosmx_p), r$cosmx_n_patients),
    "\n",
    sprintf("GeoMx PanCK+ NOR\u2192CORE log2FC %s, FDR %s",
            fmt_sg(r$GeoMx_PanCKpos_NOR_CORE_log2FC),
            fmt_p(r$GeoMx_PanCKpos_NOR_CORE_padj)),
    if (length(others)) paste0("\n", paste(others, collapse = " \u00b7 ")) else "")
}

## A raster layer carries no colour aesthetic, so nothing would key a legend.
## Each panel therefore gets one invisible point per legend level, placed inside
## the limits at alpha 0, purely to instantiate the guide. The override makes
## the swatches opaque.
ghost <- function(xr, yr, n) data.table(x = rep(xr[1], n), y = rep(yr[1], n))

## ---- panel a: composition --------------------------------------------------
comp[, group := factor(group, levels = G_ORDER)]
setorder(comp, group)               # carcinoma and normal epithelium last
ra <- rasterise_cloud(comp$x, comp$y, base_col = GCOL[as.character(comp$group)])

gh_a <- ghost(ra$xr, ra$yr, length(G_ORDER))[, group := factor(G_ORDER, levels = G_ORDER)]
pa <- ggplot() +
  annotation_raster(ra$raster, xmin = ra$xr[1], xmax = ra$xr[2],
                    ymin = ra$yr[1], ymax = ra$yr[2], interpolate = FALSE) +
  geom_point(data = gh_a, aes(x, y, colour = group), alpha = 0, size = 0.1) +
  scale_colour_manual(values = GCOL, breaks = rev(G_ORDER), name = NULL,
                      guide = guide_legend(
                        override.aes = list(alpha = 1, size = 1.9))) +
  scale_x_continuous(limits = ra$xr) + scale_y_continuous(limits = ra$yr) +
  coord_fixed(clip = "off") +
  labs(title = "Cell composition",
       subtitle = sprintf("Slide 2 \u00b7 %s cells",
                          format(nrow(comp), big.mark = ","))) +
  theme_map() +
  theme(legend.position = "inside", legend.position.inside = c(0.02, 0.30),
        legend.justification = c(0, 0.5), legend.key.height = unit(2.9, "mm"),
        legend.key.width = unit(2.9, "mm"), legend.background = element_blank(),
        plot.title = element_text(size = 8.4, face = "bold", colour = INK,
                                  hjust = 0, margin = margin(b = 1, l = 11)))

## ---- panels b-d: per-gene expression ---------------------------------------
gene_panel <- function(g) {
  d <- expr[gene == g]
  vmax <- st[g]$expr_vmax
  rg <- rasterise_cloud(d$x, d$y, value = d$expr, ramp = EXP_RAMP,
                        vmax = vmax, base_col = GREY_OFF)
  gh <- ghost(rg$xr, rg$yr, 2)[, v := c(0, vmax)]
  ggplot() +
    annotation_raster(rg$raster, xmin = rg$xr[1], xmax = rg$xr[2],
                      ymin = rg$yr[1], ymax = rg$yr[2], interpolate = FALSE) +
    geom_point(data = gh, aes(x, y, colour = v), alpha = 0, size = 0.1) +
    scale_colour_gradientn(colours = EXP_RAMP, limits = c(0, vmax),
                           oob = scales::squish,
                           name = expression(log[2]~"norm. counts"),
                           guide = guide_colourbar(barwidth = unit(2, "mm"),
                                                   barheight = unit(12, "mm"),
                                                   ticks.colour = NA,
                                                   title.position = "left")) +
    scale_x_continuous(limits = rg$xr) + scale_y_continuous(limits = rg$yr) +
    coord_fixed(clip = "off") +
    labs(title = g, subtitle = subtitle_for(g)) +
    theme_map() +
    theme(legend.position = "inside", legend.position.inside = c(0.04, 0.28),
          legend.justification = c(0, 0.5),
          legend.title = element_text(size = 5.6, angle = 90, hjust = 0.5),
          ## gene symbols are italicised per nomenclature convention
          plot.title = element_text(size = 8.4, face = "bold.italic",
                                    colour = INK, hjust = 0,
                                    margin = margin(b = 1, l = 11)))
}

pb <- gene_panel("GTF2I")
pc <- gene_panel("TRIM28")
pd <- gene_panel("UBE2V1")

## ---- compose ---------------------------------------------------------------
fig <- (pa | pb | pc | pd) +
  plot_annotation(
    tag_levels = "a",
    title = paste0("Three scTWAS risk genes are up in carcinoma epithelium ",
                   "across every platform that measures them"),
    caption = paste0(
      "Each map is one CosMx slide; colour is log2 normalised counts in panels b-d and cell-type group in panel a. Non-expressing cells are drawn in pale grey beneath the\n",
      "expressing cells, and the colour scale is capped at the 98th percentile of the expressing values, so the cap differs between panels. Panels a and c share one cell subset\n",
      "(339,739 cells) and panels b and d another (260,000 cells); both sample the same slide but keep every cell expressing one of their own genes plus a random context\n",
      "sample, so outlines differ slightly. Clouds are rendered at one cell per pixel (", PX_H, " px tall), the visible colour at a pixel being the largest value landing there, which is\n",
      "what the source panels' ascending draw order shows. scTWAS Z and P are the gene's representative cell-type model. The CosMx per-cell \u0394 is a paired carcinoma-minus-normal-epithelium\n",
      "difference across 7 patients with a nominal, unadjusted P - the smallest value a 7-pair sign test can return is 0.016, so this is corroboration rather than independent\n",
      "significance. The GeoMx PanCK+ NOR\u2192CORE statistic is the compartment-level DE that the scTWAS support tables label CosMx_NOR_CORE; it is the same measurement as\n",
      "that column, not a second one. Remaining layers are listed only where they support the gene in the direction its risk Z predicts."),
    theme = theme(
      plot.title   = element_text(size = 9.2, face = "bold", colour = INK,
                                  margin = margin(b = 3)),
      plot.caption = element_text(size = 4.9, colour = INK3, hjust = 0,
                                  lineheight = 1.3, margin = margin(t = 4)),
      plot.tag     = element_text(size = 9.6, face = "bold", colour = INK),
      plot.margin  = margin(5, 5, 4, 5))) &
  theme(plot.tag.position = c(0, 1))

W <- 12.4; H <- 5.2
ggsave(file.path(DIR_FIG, "Fig_NOR_CORE_highlight_4panel.png"), fig,
       width = W, height = H, dpi = 400, bg = "white")
## the clouds are already rasters, so the PDF stays small and its text vector
ggsave(file.path(DIR_FIG, "Fig_NOR_CORE_highlight_4panel.pdf"), fig,
       width = W, height = H, bg = "white")

## ---- source data -----------------------------------------------------------
## Cell coordinates are not repeated here - they are the two input CSVs. What
## ships with the figure is every number printed on it.
sd <- st[, .(gene, panel_source, rep_model, rep_Z, rep_P, n_models,
             direction_status, cosmx_delta, cosmx_p, cosmx_fdr,
             cosmx_n_patients, expr_vmax,
             GeoMx_PanCKpos_NOR_CORE_log2FC, GeoMx_PanCKpos_NOR_CORE_padj,
             Visium_M1_log2FC, Visium_M1_label,
             Visium_M3_log2FC, Visium_M3_label,
             TCGA_bulk_RNA_log2FC, TCGA_bulk_RNA_label,
             CPTAC_bulk_protein_log2FC, CPTAC_bulk_protein_label, n_layers)]
sd <- rbind(sd, data.table(gene = "cells_plotted", panel_source =
  sprintf("a,c: %d | b,d: %d", nrow(comp), expr[gene == "GTF2I", .N])),
  fill = TRUE)
fwrite(sd, file.path(DIR_FIG, "Fig_NOR_CORE_highlight_4panel_sourcedata.csv"),
       na = "")

cat("wrote", file.path(DIR_FIG, "Fig_NOR_CORE_highlight_4panel.png"),
    "| raster px:", PX_H, "\n")
