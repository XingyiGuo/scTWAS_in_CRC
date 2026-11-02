#!/usr/bin/env Rscript
# Wellcome Sanger healthy epithelium: BPCells -> Seurat -> three labelled UMAPs.
# Run twice if desired: first run creates fine_annotation_mapping.csv; edit its
# simplified_annotation column, then rerun to use your reviewed names.

options(bitmapType = "cairo", Seurat.object.assay.version = "v5")
set.seed(20260923)

# ---- Settings ---------------------------------------------------------------
SEED <- 20260923L
RES <- as.numeric(Sys.getenv("WS_CLUSTER_RES", "0.8"))
CAP <- Sys.getenv("WS_CELLS_PER_TYPE", "4000") # "all" for no sampling cap
CAP <- if (tolower(CAP) %in% c("all", "none", "0")) NA_integer_ else as.integer(CAP)
MIN_DIST <- as.numeric(Sys.getenv("WS_MIN_DIST", "0.5"))
N_PCS <- 50L; K <- 15L; N_HVG <- 2000L
FINE_MIN_CELLS <- as.integer(Sys.getenv("WS_FINE_MIN_CELLS", "300")) # 04c only
if (is.na(FINE_MIN_CELLS) || FINE_MIN_CELLS < 1L) stop("WS_FINE_MIN_CELLS must be a positive integer")

resolve <- function(paths, kind = "file") {
  paths <- paths[!is.na(paths) & nzchar(paths)]
  ok <- if (kind == "dir") dir.exists(paths) else file.exists(paths)
  if (!any(ok)) stop("Not found: ", paste(paths, collapse = " | "))
  paths[which(ok)[1]]
}
subdir <- "Healthy_TI_then_rectum_278samples_outputs_LauLabel"
OUT <- resolve(c(Sys.getenv("WS_LAULABEL_OUT_DIR"),
                 file.path(Sys.getenv("WS_EXPRESSION_DIR"), subdir),
                 file.path("/data/l2_bioinfo1/liq17/scTWAS/data/WellcomeSanger/expression", subdir),
                 file.path("/v5000/data/l2_bioinfo1/liq17/scTWAS/data/WellcomeSanger/expression", subdir),
                 file.path("/v5000/data/l2_bioinfo1/scRNA-seq/WellcomeSanger/processed/E-MTAB-16986", subdir),
                 file.path("L:/scRNA-seq/WellcomeSanger", subdir)), "dir")
activate <- resolve(c(Sys.getenv("WS_RENV_ACTIVATE"),
                      "/data/l2_bioinfo1/liq17/scTWAS/src/renv/activate.R",
                      "/data/l2_bioinfo1/liq17/scTWAS/data/VUMC_HTAN_GSE266546/src/renv/activate.R"))
source(activate)
suppressPackageStartupMessages({
  library(BPCells); library(Seurat); library(ggplot2); library(dplyr)
})
FIG <- file.path(OUT, "umap"); dir.create(FIG, recursive = TRUE, showWarnings = FALSE)
# ---- Which ORIGINAL -> SIMPLIFIED mapping to read ---------------------------
# The mapping is what turns a raw fine annotation into the "LauLabel:
# simplified_annotation" text panel 04c prints, so it is always used. A curated
# copy can be supplied as the FIRST ARGUMENT, or via WS_FINE_MAP:
#
#   Rscript 04_umap_epithelial_LauLabel_04c_filtered.R \
#           /path/to/fine_annotation_mapping_fixed.csv
#
# A file named that way is READ ONLY. The resolved mapping (collisions settled,
# display_label built) is always written to <FIG>/fine_annotation_mapping.csv,
# so hand-edited names are never overwritten. Passing a path that does not
# exist is an error rather than a silent fall-back to the default.
cli_args <- commandArgs(trailingOnly = TRUE)
MAP_ARG <- if (length(cli_args) >= 1L && nzchar(cli_args[1])) {
  cli_args[1]
} else {
  Sys.getenv("WS_FINE_MAP", "")
}
MAP_OUT <- file.path(FIG, "fine_annotation_mapping.csv")
if (nzchar(MAP_ARG)) {
  if (!file.exists(MAP_ARG)) stop("Mapping CSV not found: ", MAP_ARG)
  MAP_IN <- MAP_ARG
  message("Fine-annotation mapping IN : ", MAP_IN)
  if (!identical(normalizePath(MAP_IN, mustWork = FALSE),
                 normalizePath(MAP_OUT, mustWork = FALSE))) {
    message("  read-only; resolved mapping will be written to ", MAP_OUT)
  }
} else {
  MAP_IN <- MAP_OUT
  message("Fine-annotation mapping IN : ", MAP_IN,
          if (file.exists(MAP_IN)) " (existing)" else " (none yet; will be created)")
}

# Same 05-style palette used by 09.R's active match05 setting.
COLORS <- c(STM = "#00B6EB", TAC = "#A58AFF", ABS = "#F8766D",
            CT = "#C49A00", GOB = "#00C094", EE = "#53B400",
            TUF = "#FB61D7", ASC = "#2455C7", SSC = "#AD255B")
write.csv(data.frame(cell_type = names(COLORS), hex = unname(COLORS)),
          file.path(FIG, "cell_type_colors.csv"), row.names = FALSE)

# ---- Read and sample exactly as in original ---------------------------------
mat <- BPCells::open_matrix_dir(file.path(OUT, "Epithelial_expression_BPCELLS"))
meta <- read.csv(file.path(OUT, "05_Epithelial_cell_metadata_LauLabel.csv.gz"),
                 row.names = 1, check.names = FALSE)
stopifnot(all(colnames(mat) %in% rownames(meta)))
missing <- setdiff(rownames(meta), colnames(mat))
if (length(missing)) warning(length(missing), " metadata cells missing from BPCells matrix")
obs <- meta[colnames(mat), , drop = FALSE]
if (!"LauLabel" %in% names(obs)) stop("Metadata lacks LauLabel")
fine_col <- if ("finelabel" %in% names(obs)) "finelabel" else if ("Final_labels" %in% names(obs)) "Final_labels" else stop("No finelabel or Final_labels column")
if (anyNA(obs$LauLabel) || anyNA(obs[[fine_col]])) stop("Missing LauLabel or fine annotation")
if (!all(unique(obs$LauLabel) %in% names(COLORS))) stop("Unknown LauLabel: ", paste(setdiff(unique(obs$LauLabel), names(COLORS)), collapse = ", "))
set.seed(SEED)
keep <- unlist(lapply(split(rownames(obs), obs$LauLabel), function(x)
  if (is.na(CAP) || length(x) <= CAP) x else sort(sample(x, CAP))), use.names = FALSE)
md <- obs[keep, , drop = FALSE]
obj <- CreateSeuratObject(mat[, keep, drop = FALSE], meta.data = md,
                          project = "WellcomeSanger_epithelial")
obj <- NormalizeData(obj, verbose = FALSE)
obj <- FindVariableFeatures(obj, nfeatures = N_HVG, verbose = FALSE)
obj <- ScaleData(obj, features = VariableFeatures(obj), scale.max = 10, verbose = FALSE)
obj <- RunPCA(obj, features = VariableFeatures(obj), npcs = N_PCS, seed.use = 42, verbose = FALSE)
obj <- FindNeighbors(obj, dims = 1:N_PCS, k.param = K, verbose = FALSE)
obj <- FindClusters(obj, resolution = RES, algorithm = 1, random.seed = SEED, verbose = FALSE)
obj <- RunUMAP(obj, dims = 1:N_PCS, n.neighbors = K, min.dist = MIN_DIST,
               umap.method = "uwot", metric = "euclidean", seed.use = 0, verbose = FALSE)

# ---- Cluster composition, majority identity and ARI -------------------------
xy <- as.data.frame(Embeddings(obj, "umap")); names(xy) <- c("UMAP_1", "UMAP_2")
df <- data.frame(cell_id = rownames(xy), xy,
                 original_cluster = as.character(obj$seurat_clusters),
                 LauLabel = as.character(obj$LauLabel),
                 original_annotation = as.character(obj@meta.data[[fine_col]]),
                 stringsAsFactors = FALSE)
size_order <- names(sort(table(df$original_cluster), decreasing = TRUE))
df$cluster <- match(df$original_cluster, size_order)
stopifnot(identical(df$cell_id, colnames(obj)))
comp <- df %>% count(cluster, LauLabel, name = "n_cells") %>%
  group_by(cluster) %>% mutate(cluster_size = sum(n_cells), fraction = n_cells / cluster_size) %>%
  arrange(cluster, desc(n_cells)) %>% ungroup()
call <- comp %>% group_by(cluster) %>% slice_max(n_cells, n = 1, with_ties = FALSE) %>%
  transmute(cluster, cluster_size, majority_LauLabel = LauLabel, purity = round(fraction, 3)) %>% ungroup()
df$majority_LauLabel <- call$majority_LauLabel[match(df$cluster, call$cluster)]
write.csv(comp, file.path(FIG, "cluster_LauLabel_composition.csv"), row.names = FALSE)
write.csv(call, file.path(FIG, "cluster_majority_LauLabel.csv"), row.names = FALSE)
comb2 <- function(x) sum(x * (x - 1) / 2)
tab <- table(df$cluster, df$LauLabel); n <- sum(tab)
a <- comb2(tab); r <- comb2(rowSums(tab)); c <- comb2(colSums(tab))
ari <- (a - r*c/choose(n, 2)) / ((r+c)/2 - r*c/choose(n, 2))

# ---- Reviewable ORIGINAL -> SIMPLIFIED fine-annotation mapping ---------------
# This table is made from the actual dataset, not guessed from cell-type names.
# Preserve original_annotation; edit ONLY simplified_annotation between runs.
# Keep one unique simplified label within each LauLabel so no fine groups merge.
make_short <- function(x) {
  x <- trimws(x)
  x <- gsub("[ _]+", " ", x)
  x <- gsub("\\b(cells?|cell type|epithelial|epithelium)\\b", "", x, ignore.case = TRUE)
  x <- gsub("[[:space:]]+", " ", x)
  # Preserve biologically meaningful terminal +/++ markers.
  x <- trimws(x)
  ifelse(nzchar(x), x, "Unspecified")
}
map <- df %>% count(LauLabel, original_annotation, name = "n_cells") %>%
  arrange(factor(LauLabel, levels = names(COLORS)), desc(n_cells), original_annotation) %>%
  mutate(simplified_annotation = make_short(original_annotation))
if (file.exists(MAP_IN)) {
  old <- read.csv(MAP_IN, check.names = FALSE, stringsAsFactors = FALSE)
  required <- c("LauLabel", "original_annotation", "simplified_annotation")
  if (!all(required %in% names(old))) stop("Mapping CSV missing required columns: ", MAP_IN)
  key <- paste(map$LauLabel, map$original_annotation, sep = "\r")
  oldkey <- paste(old$LauLabel, old$original_annotation, sep = "\r")
  if (anyDuplicated(oldkey)) stop("Duplicate original annotation keys in mapping CSV")
  ix <- match(key, oldkey)
  edited <- !is.na(ix) & !is.na(old$simplified_annotation[ix]) & nzchar(trimws(old$simplified_annotation[ix]))
  map$simplified_annotation[edited] <- trimws(old$simplified_annotation[ix[edited]])
}
# Write the mapping before drawing panels 04b/04c, retaining any edits.
map$display_label <- paste0(map$LauLabel, ": ", map$simplified_annotation)
write.csv(map, MAP_OUT, row.names = FALSE)
# Resolve collisions from an existing edited CSV using the original annotation.
# This preserves all distinct fine groups instead of aborting or merging them.
dup <- duplicated(map$display_label) | duplicated(map$display_label, fromLast = TRUE)
if (any(dup)) {
  warning("Resolving ", sum(dup), " colliding fine labels using their original annotations")
  map$simplified_annotation[dup] <- make_short(map$original_annotation[dup])
  map$display_label <- paste0(map$LauLabel, ": ", map$simplified_annotation)
}
if (anyDuplicated(map$display_label)) {
  # Last-resort stable disambiguation for user-provided edited labels.
  map$display_label <- make.unique(map$display_label, sep = "_")
}
write.csv(map, MAP_OUT, row.names = FALSE)
ix <- match(paste(df$LauLabel, df$original_annotation, sep = "\r"),
            paste(map$LauLabel, map$original_annotation, sep = "\r"))
stopifnot(!anyNA(ix))
df$simplified_annotation <- map$simplified_annotation[ix]
df$display_label <- map$display_label[ix]
write.csv(df, file.path(FIG, "umap_WellcomeSanger_epithelial_clusters_LauLabel_sourcedata.csv"), row.names = FALSE)

# ---- Palette and UMAPs -------------------------------------------------------
# A = numbered clusters, B = LauLabel, C = "GOB: simplified annotation".
lev <- intersect(names(COLORS), unique(df$LauLabel))
obj$cluster <- factor(as.character(df$cluster), levels = as.character(sort(unique(df$cluster))))
obj$LauLabel <- factor(df$LauLabel, levels = lev)
# Order fine identities by LauLabel, then abundance. Use family shades.
shade <- function(hex, f) {
  v <- grDevices::col2rgb(hex)[, 1] / 255
  v <- v + (1 - v) * f
  grDevices::rgb(v[1], v[2], v[3])
}
map <- map %>% group_by(LauLabel) %>%
  mutate(shade_fraction = if (n() == 1) 0 else seq(0, 0.50, length.out = n()),
         hex = mapply(shade, COLORS[LauLabel], shade_fraction)) %>% ungroup()
fine_colors <- setNames(map$hex, map$display_label)
obj$fine_display <- factor(df$display_label, levels = map$display_label)
cluster_colors <- setNames(COLORS[call$majority_LauLabel[order(call$cluster)]],
                           as.character(sort(call$cluster)))

# ---- 04c label geometry: side columns hugging the cloud ---------------------
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
              colour = "black", size = label.size, vjust = 0.5,
              inherit.aes = FALSE) +
    coord_cartesian(xlim = xlim, ylim = ylim, clip = "off")
}

save_plot <- function(p, stem, w = 8, h = 6) {
  p <- p + labs(x = "UMAP 1", y = "UMAP 2")
  ggsave(file.path(FIG, paste0(stem, ".pdf")), p, width = w, height = h,
         device = cairo_pdf, limitsize = FALSE)
  ggsave(file.path(FIG, paste0(stem, ".png")), p, width = w, height = h,
         dpi = 300, limitsize = FALSE)
}
panel <- function(group, colors, text_size = 3, point_size = 0.3) {
  DimPlot(obj, reduction = "umap", group.by = group, cols = colors,
          label = TRUE, repel = TRUE, label.size = text_size,
          pt.size = point_size, raster = FALSE) + NoLegend()
}
tag <- paste0("res", format(RES, trim = TRUE))
p_a <- panel("cluster", cluster_colors, 4)
p_b <- panel("LauLabel", COLORS[lev], 4)
# ---- 04c only: remove fine-annotation groups with <300 sampled cells ------
# Counts are from the exact cells plotted, not from the full un-sampled dataset.
# Filter the display, NOT the embedding: 04a/04b and all UMAP coordinates remain intact.
fine_counts <- df %>% count(LauLabel, original_annotation, display_label, name = "n_cells") %>%
  arrange(desc(n_cells), display_label) %>%
  mutate(keep_04c = n_cells >= FINE_MIN_CELLS)
write.csv(fine_counts, file.path(FIG, "04c_fine_annotation_filter_counts.csv"), row.names = FALSE)
keep_labels <- fine_counts$display_label[fine_counts$keep_04c]
if (!length(keep_labels)) stop("No fine groups reach the 04c threshold: ", FINE_MIN_CELLS)
fine_df <- df[df$display_label %in% keep_labels, , drop = FALSE]
fine_df$display_label <- factor(fine_df$display_label, levels = keep_labels)

# One anchor per retained fine annotation, at the median position of its cells.
anchors <- fine_df %>% group_by(display_label) %>%
  summarise(UMAP_1 = median(UMAP_1), UMAP_2 = median(UMAP_2), .groups = "drop")

# Panel c: points only, then colour dot -> black leader -> text outside.
# plot_w / plot_h MUST match the save_plot() call below: the horizontal margin
# and the vertical label spacing are both solved from the figure's real size.
C_W <- 11; C_H <- 7
p_c <- side_leader_labels(
  ggplot(fine_df, aes(UMAP_1, UMAP_2, colour = display_label)) +
    geom_point(size = 0.3, alpha = 0.85, show.legend = FALSE) +
    scale_colour_manual(values = fine_colors[keep_labels], drop = FALSE) +
    theme_classic() +
    theme(legend.position = "none", plot.margin = margin(10, 10, 10, 10)),
  fine_df, "display_label", fine_colors[keep_labels],
  label.size = 2.6, plot_w = C_W, plot_h = C_H
)
message("04c: retaining ", length(keep_labels), "/", nrow(fine_counts),
        " fine groups (>=", FINE_MIN_CELLS, " cells); plotted ", nrow(fine_df),
        "/", nrow(df), " cells. See 04c_fine_annotation_filter_counts.csv")
if (length(keep_labels) > 26) {
  message("04c: ", length(keep_labels), " labels; the side columns hold about 26 at ",
          "this text size. Raise WS_FINE_MIN_CELLS or increase C_H if they look tight")
}
save_plot(p_a, paste0("04a_umap_clusters_", tag))
save_plot(p_b, paste0("04b_umap_LauLabel_", tag))
save_plot(p_c, paste0("04c_umap_finelabel_", tag), w = C_W, h = C_H)

# The three-panel composite (umap_WellcomeSanger_epithelial_clusters_LauLabel)
# is deliberately NOT drawn: panel c's outside labels need the full 11 in of
# width to stay legible, and the composite squeezed it to a third of that.
# Only the single panels 04a / 04b / 04c are written.

writeLines(capture.output(sessionInfo()), file.path(FIG, "04_umap_R_sessionInfo.txt"))
message("Wrote mapping: ", MAP_OUT, " | cell colors: ", file.path(FIG, "cell_type_colors.csv"))
message("Cells: ", nrow(df), " | clusters: ", nrow(call), " | ARI: ", round(ari, 3))
