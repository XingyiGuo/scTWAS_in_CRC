#!/usr/bin/env Rscript
# 05_annotation.R
# ----------------------------------------------------------------------------
# Annotate clusters using the dominant ORIGINAL cell type in each cluster.
# Workflow:
#   1) Load the final clustered object from 04_sub_cluster.R
#   2) Harmonize original CellType_cb labels
#   3) Calculate original cell-type composition within each cluster
#   4) Identify the dominant cell type and its percentage in each cluster
#   5) Automatically annotate only clusters with sufficient dominance
#      (default >= 60%); otherwise label as "Review"
#   6) Save composition tables, annotated object, metadata, and UMAP
# ----------------------------------------------------------------------------

options(bitmapType = "cairo")
set.seed(123)

# ---- project / renv ---------------------------------------------------------
PROJECT <- "/data/l2_bioinfo1/liq17/scTWAS/data/VUMC_HTAN_GSE266546/src"
source(file.path("/data/l2_bioinfo1/liq17/scTWAS/src/renv/activate.R"))

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
})

# ---- paths ------------------------------------------------------------------
RESULTS_DIR <- "/data/l2_bioinfo1/liq17/scTWAS/data/VUMC_HTAN_GSE266546/raw/rds2"
FIG_DIR <- file.path(RESULTS_DIR, "figures")
dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)

rds_in <- file.path(
  RESULTS_DIR,
  "integrated_seurat_final1.res0.8.rm23_sub.rds"
)

stopifnot(file.exists(rds_in))
message("Loading: ", rds_in)
obj <- readRDS(rds_in)

# ---- input checks ------------------------------------------------------------
required_meta <- c("cluster", "CellType_cb")
missing_meta <- setdiff(required_meta, colnames(obj@meta.data))
if (length(missing_meta) > 0) {
  stop("Missing metadata columns: ", paste(missing_meta, collapse = ", "))
}

if (!"umap_harmony" %in% Reductions(obj)) {
  stop("Reduction 'umap_harmony' is missing from the object.")
}

message("Clusters from 04:")
print(table(obj$cluster, useNA = "ifany"))

# =============================================================================
# 1. Harmonize the ORIGINAL annotations
# =============================================================================
obj$CellType_cb2 <- as.character(obj$CellType_cb)

obj$CellType_cb2[
  grepl("Goblet", obj$CellType_cb2, ignore.case = TRUE)
] <- "GOB"

obj$CellType_cb2[
  grepl("Stem", obj$CellType_cb2, ignore.case = TRUE)
] <- "STM"

obj$CellType_cb2[
  grepl("Tuft", obj$CellType_cb2, ignore.case = TRUE)
] <- "TUF"

obj$CellType_cb2[
  grepl("EEC|Enteroendocrine", obj$CellType_cb2, ignore.case = TRUE)
] <- "EE"

obj$CellType_cb3 <- as.character(obj$CellType_cb2)

obj$CellType_cb3[
  grepl("Colonocyte|Enterocyte", obj$CellType_cb3, ignore.case = TRUE)
] <- "ABS"

obj$CellType_cb3[
  grepl("Cycling TA", obj$CellType_cb3, ignore.case = TRUE)
] <- "TAC"

# Save original-label cross-tabs for inspection
original_by_cluster <- as.data.frame(
  table(obj$CellType_cb, obj$cluster, useNA = "ifany")
) %>%
  rename(CellType_cb = Var1, cluster = Var2, n = Freq) %>%
  filter(n > 0) %>%
  group_by(cluster) %>%
  arrange(desc(n), .by_group = TRUE) %>%
  ungroup()

write.csv(
  original_by_cluster,
  file.path(
    RESULTS_DIR,
    "cluster_by_original_CellType_cb_res0.8.rm23_sub.csv"
  ),
  row.names = FALSE
)

# =============================================================================
# 2. Cell-type composition within each cluster
# =============================================================================
# Cluster total includes all cells.
cluster_totals <- obj@meta.data %>%
  transmute(cluster = as.character(cluster)) %>%
  count(cluster, name = "cluster_total")

# Composition is calculated from cells that have a usable original annotation.
cluster_comp <- obj@meta.data %>%
  transmute(
    cluster = as.character(cluster),
    celltype = as.character(CellType_cb3)
  ) %>%
  filter(!is.na(celltype), celltype != "") %>%
  count(cluster, celltype, name = "n_cells") %>%
  group_by(cluster) %>%
  mutate(
    annotated_total = sum(n_cells),
    percent_of_annotated = round(100 * n_cells / annotated_total, 2)
  ) %>%
  ungroup() %>%
  left_join(cluster_totals, by = "cluster") %>%
  mutate(
    percent_of_cluster = round(100 * n_cells / cluster_total, 2)
  ) %>%
  arrange(as.numeric(cluster), desc(percent_of_annotated))

write.csv(
  cluster_comp,
  file.path(
    RESULTS_DIR,
    "cluster_celltype_composition_res0.8.rm23_sub.csv"
  ),
  row.names = FALSE
)

message("\nCell-type composition by cluster:")
print(cluster_comp, n = Inf)

# =============================================================================
# 3. Dominant original cell type for each cluster
# =============================================================================
dominant_celltype <- cluster_comp %>%
  group_by(cluster) %>%
  slice_max(
    order_by = percent_of_annotated,
    n = 1,
    with_ties = FALSE
  ) %>%
  ungroup() %>%
  arrange(as.numeric(cluster)) %>%
  rename(
    dominant_celltype = celltype,
    dominant_n = n_cells,
    dominant_percent = percent_of_annotated
  )

# ---- confidence threshold ----------------------------------------------------
# A cluster is automatically assigned to its dominant original cell type only
# if that type accounts for at least this percentage of annotated cells.
# Change this value if you want a stricter or more permissive criterion.
DOMINANT_MIN_PCT <- 60

dominant_celltype <- dominant_celltype %>%
  mutate(
    annotation_status = ifelse(
      dominant_percent >= DOMINANT_MIN_PCT,
      "auto",
      "Review"
    ),
    proposed_annotation = ifelse(
      dominant_percent >= DOMINANT_MIN_PCT,
      dominant_celltype,
      "Review"
    )
  )

write.csv(
  dominant_celltype,
  file.path(
    RESULTS_DIR,
    "cluster_dominant_celltype_res0.8.rm23_sub.csv"
  ),
  row.names = FALSE
)

message("\nDominant cell type for each cluster:")
print(
  dominant_celltype %>%
    select(
      cluster,
      dominant_celltype,
      dominant_n,
      cluster_total,
      dominant_percent,
      annotation_status,
      proposed_annotation
    ),
  n = Inf
)

# =============================================================================
# 4. Assign cluster-level annotation
# =============================================================================
annotation_map <- setNames(
  dominant_celltype$proposed_annotation,
  dominant_celltype$cluster
)

obj$celltype_annot <- unname(
  annotation_map[as.character(obj$cluster)]
)

# Any cluster absent from the composition table should be reviewed manually.
obj$celltype_annot[is.na(obj$celltype_annot)] <- "Review"

message("\nCluster -> final annotation:")
print(table(obj$cluster, obj$celltype_annot, useNA = "ifany"))

# =============================================================================
# 5. Detailed comparison table
# =============================================================================
a4 <- obj@meta.data %>%
  mutate(
    cluster = as.character(cluster),
    CellType_cb = as.character(CellType_cb),
    CellType_cb3 = as.character(CellType_cb3),
    celltype_annot = as.character(celltype_annot)
  ) %>%
  group_by(cluster, CellType_cb, CellType_cb3, celltype_annot) %>%
  summarise(n_cells = n(), .groups = "drop") %>%
  group_by(cluster) %>%
  mutate(
    cluster_total = sum(n_cells),
    percent_within_cluster = round(100 * n_cells / cluster_total, 2)
  ) %>%
  arrange(as.numeric(cluster), desc(percent_within_cluster)) %>%
  ungroup()

write.csv(
  a4,
  file.path(
    RESULTS_DIR,
    "cluster_by_original_CellType_cb_annot_pct_res0.8.rm23_sub.csv"
  ),
  row.names = FALSE
)

# =============================================================================
# 6. Save annotated object and metadata
# =============================================================================
message(
  "    ", nrow(obj), " genes x ",
  ncol(obj), " cells"
)

message(
  "    ", length(unique(obj$SampleID)), " samples, ",
  length(unique(obj$PatientID)), " patients"
)

message(
  "    ", length(unique(obj$cluster)), " clusters; ",
  length(unique(obj$celltype_annot)), " annotation labels"
)

saveRDS(
  obj,
  file.path(
    RESULTS_DIR,
    "integrated_seurat_final2.res0.8.rm23_sub.rds"
  )
)

# write.csv(
  # obj@meta.data,
  # file.path(
    # RESULTS_DIR,
    # "integrated_seurat_final2.res0.8.rm23_sub_metadata.csv"
  # ),
  # row.names = FALSE
# )

# barcode_both <- data.frame(
  # barcode = rownames(obj@meta.data),
  # CellType_cb = obj$CellType_cb,
  # CellType_cb3 = obj$CellType_cb3,
  # cluster = as.character(obj$cluster),
  # celltype_annot = obj$celltype_annot,
  # stringsAsFactors = FALSE
# )

# write.csv(
  # barcode_both,
  # file.path(
    # RESULTS_DIR,
    # "barcode_original_annotation2_res0.8.rm23_sub.csv"
  # ),
  # row.names = FALSE
# )

# =============================================================================
# 7. UMAP plots
# =============================================================================
# ---- Cell-type label geometry: side columns hugging the cloud ---------------
# Ported VERBATIM from 04_umap_epithelial_LauLabel_04c_filtered.R (densest-bin
# anchor rule from Final_Viz_4panel.R) so every UMAP panel in the project
# places its labels the same way. Both functions read UMAP_1/UMAP_2 off the
# data frame handed to them, so this object's UMAPHARMONY_1/2 columns are
# copied under those names at the call site rather than forking the helpers.
# Anchors are each group's DENSEST 40 x 40 bin, not its median: several fine
# groups have two or more lobes and the median then lands in the gap between
# them, which is how a leader line ends up pointing at another group's cells.
# (Same rule Final_Viz_4panel.R uses.)
label_anchors <- function(df, group) {
  d <- df
  d$.bx <- cut(d$UMAP_1, 40, labels = FALSE)
  d$.by <- cut(d$UMAP_2, 40, labels = FALSE)
  top <- d %>% count(.data[[group]], .bx, .by, name = ".n") %>%
    group_by(.data[[group]]) %>% slice_max(.n, n = 1, with_ties = FALSE) %>%
    ungroup() %>% select(all_of(group), .bx, .by)
  cen <- d %>% inner_join(top, by = c(group, ".bx", ".by")) %>%
    group_by(.data[[group]]) %>%
    summarise(x = median(UMAP_1), y = median(UMAP_2), .groups = "drop")
  names(cen)[1] <- "label"
  cen$label <- as.character(cen$label)
  cen
}

# Labels sit in a LEFT and a RIGHT column, each just clear of the point cloud
# rather than on a fixed ring: for a given label the column steps in to the
# cloud's own outline at that height, so a group near the edge gets a short
# leader instead of one stretched to a common radius. Within a column the
# labels keep their anchors' vertical order and are pushed apart only as far as
# the text height requires, so they read as an ordered list, not a scatter.
#
# The horizontal margin is COMPUTED from the longest label on each side, its
# font size and the figure width, then the panel limits are solved so the text
# ends inside the panel: x_total = Rx / (1 - (cL + cR) / Uw). The previous
# fixed 10% expansion is what clipped the long ABS/EE/TAC names.
side_leader_labels <- function(p, df, group, colours, label.size = 2.6,
                               label.family = "",
                               plot_w = 11, plot_h = 7, gap_in = 0.07,
                               line.width = 0.25, dot.size = 1.2,
                               axis_w_in = 1.15, axis_h_in = 0.95) {

  an <- label_anchors(df, group)

  xr <- range(df$UMAP_1); yr <- range(df$UMAP_2)
  Rx <- diff(xr); Ry <- diff(yr)
  cx <- mean(xr)

  Uw <- max(plot_w - axis_w_in, 1)   # inches of panel actually available
  Uh <- max(plot_h - axis_h_in, 1)

  fs_pt   <- label.size * .pt        # geom_text size is mm; .pt converts to pt
  char_in <- 0.5 * fs_pt / 72        # mean glyph width
  line_in <- 1.25 * fs_pt / 72       # line height

  an$side <- ifelse(an$x >= cx, "R", "L")
  an$w_in <- nchar(an$label) * char_in + gap_in

  cL <- if (any(an$side == "L")) max(an$w_in[an$side == "L"]) else 0
  cR <- if (any(an$side == "R")) max(an$w_in[an$side == "R"]) else 0
  if ((cL + cR) / Uw > 0.8) {
    warning("04c labels need ", round(100 * (cL + cR) / Uw),
            "% of the panel width; capping at 80%. Widen the figure or shorten",
            " simplified_annotation in the mapping CSV.")
    sc <- 0.8 * Uw / (cL + cR); cL <- cL * sc; cR <- cR * sc
  }
  x_total <- Rx / (1 - (cL + cR) / Uw)
  xlim <- c(xr[1] - cL * x_total / Uw, xr[2] + cR * x_total / Uw)
  gap_data <- gap_in * x_total / Uw

  # vertical room: give each column at least one text line per label
  dy <- line_in / Uh * Ry * 1.05
  need <- max(tapply(an$label, an$side, length)) * dy
  ylim <- if (need > Ry) {
    c(mean(yr) - need / 2, mean(yr) + need / 2)
  } else {
    yr
  }
  dy <- line_in / Uh * diff(ylim) * 1.05

  place <- function(s) {
    d <- an[an$side == s, , drop = FALSE]
    if (!nrow(d)) return(NULL)
    d <- d[order(-d$y), , drop = FALSE]
    ly <- d$y
    for (it in seq_len(400)) {
      moved <- FALSE
      if (nrow(d) > 1) for (i in seq_len(nrow(d) - 1)) {
        g <- ly[i] - ly[i + 1]
        if (g < dy) {
          push <- (dy - g) / 2
          ly[i] <- ly[i] + push; ly[i + 1] <- ly[i + 1] - push; moved <- TRUE
        }
      }
      if (ly[1] > ylim[2]) ly <- ly - (ly[1] - ylim[2])
      if (ly[nrow(d)] < ylim[1]) ly <- ly + (ylim[1] - ly[nrow(d)])
      if (!moved) break
    }
    d$lab_y <- ly
    d$lab_x <- vapply(seq_len(nrow(d)), function(i) {
      sel <- df$UMAP_2 >= ly[i] - dy & df$UMAP_2 <= ly[i] + dy
      edge <- if (!any(sel)) {
        if (s == "R") xr[2] else xr[1]
      } else if (s == "R") max(df$UMAP_1[sel]) else min(df$UMAP_1[sel])
      if (s == "R") edge + gap_data else edge - gap_data
    }, numeric(1))
    d$hjust <- if (s == "R") 0 else 1
    d$seg_x <- if (s == "R") d$lab_x - gap_data * 0.5 else d$lab_x + gap_data * 0.5
    d
  }

  lab <- bind_rows(place("L"), place("R"))

  # Tighten the panel to the text that is actually drawn. The margin solve above
  # reserves room at the DATA edge, but each label then steps in to the cloud's
  # outline, so that reservation leaves dead space at both sides and squeezes
  # the embedding. Re-solve the limits from the real label positions (a few
  # passes: the glyph width itself depends on the limits).
  for (it in seq_len(4)) {
    w_data <- nchar(lab$label) * char_in * diff(xlim) / Uw
    left  <- min(c(xr[1], lab$lab_x[lab$side == "L"] - w_data[lab$side == "L"]))
    right <- max(c(xr[2], lab$lab_x[lab$side == "R"] + w_data[lab$side == "R"]))
    xlim <- c(left - gap_data, right + gap_data)
  }

  p +
    geom_segment(data = lab, aes(x = seg_x, y = lab_y, xend = x, yend = y),
                 colour = "black", linewidth = line.width, inherit.aes = FALSE) +
    geom_point(data = lab, aes(x = x, y = y, fill = label), shape = 21,
               colour = "black", stroke = 0.25, size = dot.size,
               inherit.aes = FALSE, show.legend = FALSE) +
    scale_fill_manual(values = colours, guide = "none") +
    geom_text(data = lab, aes(x = lab_x, y = lab_y, label = label, hjust = hjust),
              colour = "black", size = label.size, family = label.family,
              vjust = 0.5, inherit.aes = FALSE) +
    coord_cartesian(xlim = xlim, ylim = ylim, clip = "off")
}

# How the cell-type panels (05c, 05e) name their types.
#   "leader"  each type is named OUTSIDE the cloud in a left/right column and
#             joined to its own cells by a thin black line ending on a
#             colour-filled dot -- the layout settled on for 04c.
#   "legend"  the previous behaviour: no text on the clusters, a right-hand key.
#   "on-plot" the original look: label = TRUE + NoLegend().
# 05d is NOT affected: the cluster-ID review panel keeps its numbers on the
# clusters, because reading an id off the plot is the whole point of it.
# ---- Cells drawn per cell type (FIGURES ONLY) -------------------------------
# The panels equalise visual weight by sampling at most CELLS_PER_TYPE cells
# from each plotted identity, as 04_umap_epithelial_LauLabel and
# 09_umap05e_ADSER_polyp_celltype do.
#
# READ THIS BEFORE CHANGING IT: the cap applies to the PLOTS and nothing else.
# Every table, the dominance test, the annotation assignment and both saved RDS
# objects above are computed on ALL cells -- a composition percentage or a
# dominant-type call derived from a 4,000-cell cap would be a different
# quantity, and 'Review' status would flip for clusters near the threshold.
# The sampling is seeded and written out per panel, so a figure can always be
# tied back to the cells behind it.
#
# Set CELLS_PER_TYPE=all (or 0 / none) to draw every cell.
cap_env <- tolower(trimws(Sys.getenv("CELLS_PER_TYPE", "4000")))
CELLS_PER_TYPE <- if (cap_env %in% c("all", "none", "0")) {
  NA_integer_
} else {
  as.integer(cap_env)
}
SAMPLE_SEED <- 20260923L

# 05d is NOT capped: it is the cluster-ID review panel, and a per-cell-type cap
# would distort the cluster sizes it exists to show.
sample_cells_by <- function(group, cap = CELLS_PER_TYPE, seed = SAMPLE_SEED) {
  ids <- colnames(obj)
  if (is.na(cap)) return(ids)
  g <- as.character(obj@meta.data[[group]])
  set.seed(seed)
  # sample() is only reached when the group is larger than the cap, so no
  # per-group min() is needed.
  unlist(lapply(split(ids, g), function(x) {
    if (length(x) <= cap) x else sort(sample(x, cap))
  }), use.names = FALSE)
}

CELLTYPE_LABEL_MODE <- "leader"
if (!CELLTYPE_LABEL_MODE %in% c("leader", "legend", "on-plot")) {
  stop("CELLTYPE_LABEL_MODE must be leader, legend or on-plot; got: ",
       CELLTYPE_LABEL_MODE)
}

# "leader" needs horizontal room for the label columns. The margin is solved
# from the real text width, so these two numbers are passed to the labeller as
# well as to save_plot() and must agree.
CT_PANEL_W <- if (identical(CELLTYPE_LABEL_MODE, "leader")) 11 else 8
CT_PANEL_H <- if (identical(CELLTYPE_LABEL_MODE, "leader")) 7 else 6

AXIS_X <- "UMAP 1"
AXIS_Y <- "UMAP 2"

# ---- Typeface ---------------------------------------------------------------
# One size for everything on the panel: axis titles, tick labels, legend and
# the leader labels.
BASE_PT <- 11

# Naming a family is not the same as getting it. fontconfig SUBSTITUTES
# silently, so a Linux run asking for Arial can be rendered in DejaVu Sans with
# no warning and slightly different metrics -- which matters here because the
# label margin is solved from character counts. So the family is verified
# before use, exactly as the four W:/VUMC_projects/scTWAS/Viz Manhattan scripts
# do: the resolved family must be the one requested, its file must exist, it
# must be scalable, and glyphs must actually be present. The first preference
# that passes wins; if none do, "" is used, which is the device default.
#   PLOT_FONT_OVERRIDE <- ""       force the device default
#   PLOT_FONT_OVERRIDE <- "Arial"  force Arial even if it does not verify
PLOT_FONT_PREFS    <- c("Arial", "Helvetica", "Liberation Sans", "DejaVu Sans")
PLOT_FONT_OVERRIDE <- NA_character_

font_is_usable <- function(fam) {
  if (!nzchar(fam)) return(FALSE)
  if (!requireNamespace("systemfonts", quietly = TRUE)) return(FALSE)
  ok <- tryCatch({
    fi <- systemfonts::font_info(family = fam)
    gi <- systemfonts::glyph_info(c("A", "g", "1"), family = fam)
    identical(as.character(fi$family)[1], fam) &&
      file.exists(as.character(fi$path)[1]) &&
      isTRUE(as.logical(fi$scalable)[1]) &&
      all(gi$index > 0)
  }, error = function(e) FALSE)
  isTRUE(ok)
}

resolve_plot_font <- function(prefs = PLOT_FONT_PREFS,
                              override = PLOT_FONT_OVERRIDE) {
  if (!is.na(override)) {
    message("Font: forced to '", override, "' by PLOT_FONT_OVERRIDE",
            if (!nzchar(override)) " (device default)" else "")
    return(override)
  }
  if (!requireNamespace("systemfonts", quietly = TRUE)) {
    message("Font: systemfonts is not installed, so '", prefs[1], "' cannot be ",
            "verified; using the device default. install.packages('systemfonts')",
            " to get the requested family.")
    return("")
  }
  for (f in prefs) {
    if (font_is_usable(f)) {
      message("Font: using '", f, "'",
              if (!identical(f, prefs[1])) {
                paste0(" ('", prefs[1], "' did not verify on this machine)")
              } else "")
      return(f)
    }
  }
  message("Font: none of ", paste(prefs, collapse = ", "),
          " verified; using the device default")
  ""
}

PLOT_FONT <- resolve_plot_font()

# geom_text() sizes in MILLIMETRES, so a point size has to be divided by .pt.
LABEL_SIZE_MM <- BASE_PT / .pt

font_theme <- function() {
  theme(
    text         = element_text(family = PLOT_FONT, size = BASE_PT),
    axis.title   = element_text(family = PLOT_FONT, size = BASE_PT),
    axis.text    = element_text(family = PLOT_FONT, size = BASE_PT),
    legend.text  = element_text(family = PLOT_FONT, size = BASE_PT),
    legend.title = element_text(family = PLOT_FONT, size = BASE_PT)
  )
}

save_plot <- function(p, name, w = 8, h = 6) {
  # DimPlot titles a panel with its group.by column, so these came out headed
  # "celltype_annot", "cluster" and "celltype_annot_reviewed". Dropped here,
  # together with the axis titles that were previously repeated at each call
  # site, so one place governs every panel and none can drift.
  p <- p + labs(x = AXIS_X, y = AXIS_Y, title = NULL) + font_theme()
  ggsave(
    filename = file.path(FIG_DIR, paste0(name, ".pdf")),
    plot = p,
    width = w,
    height = h,
    device = cairo_pdf
  )
  # PNG twin at 300 dpi, as the other figure scripts in this project deliver.
  ggsave(
    filename = file.path(FIG_DIR, paste0(name, ".png")),
    plot = p,
    width = w,
    height = h,
    dpi = 300
  )
  message("Wrote: ", file.path(FIG_DIR, name), ".{pdf,png}")
}

# One builder for both cell-type panels, so 05c and 05e cannot diverge.
celltype_panel <- function(group, name) {

  cells_use <- sample_cells_by(group)

  # What was available and what was drawn, per identity, beside the tables.
  g_all <- as.character(obj@meta.data[[group]])
  g_use <- g_all[match(cells_use, colnames(obj))]
  sampling <- merge(
    as.data.frame(table(identity = g_all), responseName = "available"),
    as.data.frame(table(identity = g_use), responseName = "plotted"),
    by = "identity", all = TRUE
  )
  write.csv(sampling, file.path(RESULTS_DIR, paste0(name, "_sampling.csv")),
            row.names = FALSE)
  message("\n", name, ": plotted ", length(cells_use), " of ", ncol(obj),
          " cells (cap ",
          if (is.na(CELLS_PER_TYPE)) "none" else CELLS_PER_TYPE,
          " per identity, seed ", SAMPLE_SEED, ")")
  print(sampling, row.names = FALSE)

  p <- DimPlot(
    obj,
    reduction = "umap_harmony",
    group.by = group,
    cells = cells_use,
    label = identical(CELLTYPE_LABEL_MODE, "on-plot"),
    raster = FALSE
  )

  if (identical(CELLTYPE_LABEL_MODE, "legend")) {
    p <- p +
      labs(colour = "Cell type") +
      guides(colour = guide_legend(override.aes = list(size = 3))) +
      # sizes come from font_theme() in save_plot(); override.aes stays because
      # legend keys would otherwise inherit the tiny UMAP point size
      theme(legend.position = "right")
  } else {
    p <- p + NoLegend()
  }

  if (identical(CELLTYPE_LABEL_MODE, "leader")) {
    lv <- sort(unique(as.character(obj@meta.data[[group]])))
    # This script sets no manual palette, so the anchor dots have to reproduce
    # the colours DimPlot picked: scales::hue_pal() walked along the identity
    # levels in SORTED order, which is how Seurat factors them.
    cols <- setNames(scales::hue_pal()(length(lv)), lv)
    # Anchors from the PLOTTED cells only: a label placed off the full cell set
    # could point at a spot where this panel draws nothing.
    e <- as.data.frame(Embeddings(obj, "umap_harmony"))[cells_use, 1:2]
    d <- data.frame(
      UMAP_1 = e[[1]],
      UMAP_2 = e[[2]],
      label = g_use,
      stringsAsFactors = FALSE
    )
    p <- side_leader_labels(p, d, "label", cols, label.size = LABEL_SIZE_MM,
                            label.family = PLOT_FONT,
                            plot_w = CT_PANEL_W, plot_h = CT_PANEL_H)
    message(name, ": ", length(lv), " cell types named on leader lines")
  }

  save_plot(p, name, w = CT_PANEL_W, h = CT_PANEL_H)
}

celltype_panel("celltype_annot",
               "05c_umap_dominant_celltype_res0.8.rm23_sub")

# 05d is the cluster-ID review panel: the numbers stay ON the clusters, because
# reading an id off the plot is what makes it usable. Axis titles now come from
# save_plot(), which also drops DimPlot's "cluster" panel title.
save_plot(
  DimPlot(
    obj,
    reduction = "umap_harmony",
    group.by = "cluster",
    label = TRUE,
    raster = FALSE
  ) + NoLegend(),
  "05d_umap_cluster_for_annotation_review_res0.8.rm23_sub"
)

# =============================================================================
# 8. Reviewed annotation: force selected Review clusters to dominant cell type
# =============================================================================
# After inspecting cluster composition / markers, clusters 8 and 10 are accepted
# as their dominant original cell type even if dominant_percent < DOMINANT_MIN_PCT.
REVIEW_USE_DOMINANT <- c("8", "10")

dominant_map_raw <- setNames(
  dominant_celltype$dominant_celltype,
  dominant_celltype$cluster
)

missing_review_clusters <- setdiff(
  REVIEW_USE_DOMINANT,
  names(dominant_map_raw)
)
if (length(missing_review_clusters) > 0) {
  stop(
    "No dominant cell type available for reviewed cluster(s): ",
    paste(missing_review_clusters, collapse = ", ")
  )
}

# Keep the automatic annotation intact and make a separate reviewed column.
obj$celltype_annot_reviewed <- as.character(obj$celltype_annot)

idx_review <- as.character(obj$cluster) %in% REVIEW_USE_DOMINANT
obj$celltype_annot_reviewed[idx_review] <- unname(
  dominant_map_raw[as.character(obj$cluster[idx_review])]
)

message("\nReviewed cluster -> final annotation:")
print(
  table(
    obj$cluster,
    obj$celltype_annot_reviewed,
    useNA = "ifany"
  )
)

reviewed_annotation_summary <- dominant_celltype %>%
  transmute(
    cluster = as.character(cluster),
    dominant_celltype,
    dominant_percent,
    automatic_annotation = proposed_annotation,
    manual_override = cluster %in% REVIEW_USE_DOMINANT,
    reviewed_annotation = ifelse(
      cluster %in% REVIEW_USE_DOMINANT,
      dominant_celltype,
      proposed_annotation
    )
  ) %>%
  arrange(as.numeric(cluster))

write.csv(
  reviewed_annotation_summary,
  file.path(
    RESULTS_DIR,
    "cluster_reviewed_annotation_res0.8.rm23_sub.csv"
  ),
  row.names = FALSE
)

# 05e: reviewed annotation UMAP. Clusters 8 and 10 now use dominant cell type.
celltype_panel("celltype_annot_reviewed",
               "05e_umap_reviewed_celltype_res0.8.rm23_sub")

# Save an additional reviewed object so the automatic final2 object is preserved.
saveRDS(
  obj,
  file.path(
    RESULTS_DIR,
    "integrated_seurat_final2_reviewed.res0.8.rm23_sub.rds"
  )
)

write.csv(
  obj@meta.data,
  file.path(
    RESULTS_DIR,
    "integrated_seurat_final2_reviewed.res0.8.rm23_sub_metadata.csv"
  ),
  row.names = FALSE
)

message("Done. Annotation outputs written to: ", RESULTS_DIR)
