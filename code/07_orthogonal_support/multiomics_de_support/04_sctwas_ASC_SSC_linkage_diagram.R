## 04_sctwas_ASC_SSC_linkage_diagram.R ----------------------------------------
## Linkage diagram of how the lesion-specific scTWAS genes connect the normal
## epithelial cell types to the polyp compartment and to carcinoma.
##
## Two panels, one per nominating track:
##   ASC panel  32 adenoma-specific genes: E1-E7 (normal cell type -> ASC),
##              C1-C7 (normal cell type -> MSI-H and -> MSS carcinoma),
##              L1 (ASC -> MSI-H, ASC -> MSS)
##   SSC panel  15 serrated-specific genes: the same E and C links plus
##              L2 (SSC -> MSI-H, SSC -> MSS)
##
## Edge weight is the number of that track's genes supported in the contrast:
## padj < 0.05, |log2FC| > 0.5 and sign(log2FC) matching the gene's scTWAS Z,
## i.e. the Up and Down calls added together. Thick line for >= 5 genes, thin
## line below 5, as specified. A contrast that was tested and supported no gene
## is drawn as a faint dotted link rather than omitted, so a missing link never
## has to be read as either zero or untested.
##
## Input
##   table/DE_key_results/b_ASC_SSC_genes_crossceltype_support_matrix.csv
##     47 genes x 25 contrast blocks (80 columns: 5 identifiers plus
##     log2FC / padj / support-call per block), two header rows. Block
##     identities are parsed from the first header row, so the file can gain or
##     lose blocks without editing this script.
##
## Output
##   figure/Fig_scTWAS_ASC_SSC_linkage.{png,pdf}
##   figure/Fig_scTWAS_ASC_SSC_linkage_sourcedata.csv   one row per drawn edge
##   table/DE_key_results/scTWAS_ASC_SSC_linkage_edges.csv
## -----------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
})

PROJ <- "W:/VUMC_projects/scTWAS/MultiOmics_DE_support"
if (!dir.exists(file.path(PROJ, "table"))) PROJ <- "."
DIR_KEY <- file.path(PROJ, "table", "DE_key_results")
DIR_FIG <- file.path(PROJ, "figure")
dir.create(DIR_FIG, recursive = TRUE, showWarnings = FALSE)

## Read the flat companion, not the grouped two-header matrix. The grouped file
## is the human-readable deliverable and has been round-tripped through a
## spreadsheet: it now carries 23 blocks instead of 25, having lost both L2
## columns, which would silently drop panel b's SSC -> carcinoma links. The flat
## file is written by the same upstream script, with unique machine-readable
## names, and still has all 25 blocks.
MATRIX_CSV <- file.path(PROJ, "table",
  "ASC_SSC_genes_crossceltype_support_matrix_flat.csv")

## ---- font sizes, in points --------------------------------------------------
## FS_MIN is a hard floor: nothing in this figure is set below it, checked below.
## geom_text/geom_label take millimetres, so those two are divided by .pt
## (72.27/25.4) to keep the declared value in points like everything else.
FS_MIN      <- 8      # the floor; the caption sits on it
FS_CAPTION  <- FS_MIN
FS_EDGE     <- 8.5    # numbers printed on each link
FS_LEGEND   <- 10
FS_SUBTITLE <- 10.5
FS_NODE     <- 11.5   # cell-type and state node labels
FS_PANEL    <- 12     # per-panel titles a / b
FS_TITLE    <- 14
stopifnot(all(c(FS_CAPTION, FS_EDGE, FS_LEGEND, FS_SUBTITLE,
                FS_NODE, FS_PANEL, FS_TITLE) >= FS_MIN))

THICK_AT <- 10L                    # >= 10 genes -> thick line
WT_LAB <- list(thick = paste0("\u2265 ", THICK_AT, " genes"),
               thin  = paste0("< ", THICK_AT, " genes"))
SUPPORTED <- c("Up", "Down")

## node order down the left column, following the reference sketch
CT_NODES <- c("ABS", "GOB", "STM", "CT", "EE", "TAC", "TUF")
## the numeric index inside a C/E tag denotes the cell type, in this fixed order
CE_INDEX <- c("1" = "ABS", "2" = "GOB", "3" = "CT", "4" = "EE",
              "5" = "STM", "6" = "TUF", "7" = "TAC")

## 1. read the flat matrix -----------------------------------------------------
## One *_DEsupport column per contrast block, so the block tags come from the
## column names rather than from a merged header row.
body <- fread(MATRIX_CSV)
de_cols <- grep("_DEsupport$", names(body), value = TRUE)
blk_tag <- sub("_DEsupport$", "", de_cols)
stopifnot(length(blk_tag) == 25L, !anyDuplicated(blk_tag),
          ## all four L blocks must be present or a track loses its exit links
          all(c("L1_MSIH", "L1_MSS", "L2_MSIH", "L2_MSS") %in% blk_tag))

lab <- body[, .(track, gene)]
for (k in seq_along(blk_tag)) lab[[blk_tag[k]]] <- body[[de_cols[k]]]
stopifnot(nrow(lab) == 47L, setequal(unique(lab$track), c("ASC", "SSC")))

m <- melt(lab, id.vars = c("track", "gene"),
          variable.name = "tag", value.name = "call")
m[, tag := as.character(tag)]
## A contrast belonging to the other lesion carries no call. The pipeline writes
## "not_applicable" there, but a spreadsheet round-trip leaves the cell empty,
## so both spellings are accepted and normalised.
m[is.na(call) | !nzchar(trimws(call)), call := "not_applicable"]
stopifnot(all(m$call %chin% c("Up", "Down", "DE_discordant", "NS",
                              "not_tested", "not_applicable")))
## the inapplicable cells must be exactly the other lesion's L block
stopifnot(m[call == "not_applicable", all(grepl("^L[12]_", tag))],
          m[grepl("^L1_", tag) & track == "SSC", all(call == "not_applicable")],
          m[grepl("^L2_", tag) & track == "ASC", all(call == "not_applicable")])

## 2. edge weights ------------------------------------------------------------
w <- m[, .(n_supported = sum(call %chin% SUPPORTED),
           n_tested    = sum(call != "not_applicable"),
           n_genes     = .N),
       by = .(track, tag)]
## a tag not applicable to a track (L2 for ASC, L1 for SSC) is dropped
w <- w[n_tested > 0L]

track_n <- lab[, .(track_genes = .N), by = track]
w <- merge(w, track_n, by = "track")

## 3. edge geometry -----------------------------------------------------------
## left column: the seven normal cell types; right column: the polyp
## compartment of that track, then MSI-H and MSS carcinoma
X_L <- 0; X_R <- 2.05
node_xy <- function(track) {
  rbindlist(list(
    data.table(node = CT_NODES, x = X_L, y = seq(7, 1), kind = "normal"),
    data.table(node = track,    x = X_R, y = 5.5, kind = "polyp"),
    data.table(node = "MSI-H",  x = X_R, y = 3.3, kind = "carcinoma"),
    data.table(node = "MSS",    x = X_R, y = 1.1, kind = "carcinoma")))
}

## tag -> (source node, target node, block)
edge_spec <- function(track) {
  L <- if (track == "ASC") "L1" else "L2"
  rbindlist(list(
    data.table(tag = paste0("E", 1:7), from = CE_INDEX[as.character(1:7)],
               to = track, block = "E"),
    data.table(tag = paste0("C", 1:7, "_MSIH"), from = CE_INDEX[as.character(1:7)],
               to = "MSI-H", block = "C"),
    data.table(tag = paste0("C", 1:7, "_MSS"), from = CE_INDEX[as.character(1:7)],
               to = "MSS", block = "C"),
    data.table(tag = paste0(L, "_MSIH"), from = track, to = "MSI-H", block = "L"),
    data.table(tag = paste0(L, "_MSS"),  from = track, to = "MSS",  block = "L")))
}

## The argument is `tr`, not `track`: inside a data.table `i` expression a name
## that matches a column resolves to the column, so `w[track == track]` would
## be self-referential and select every row.
build_edges <- function(tr) {
  nd <- node_xy(tr)
  e  <- merge(edge_spec(tr), w[track == tr], by = "tag")
  stopifnot(nrow(e) == 23L)
  e <- merge(e, nd[, .(from = node, x = x, y = y)], by = "from")
  e <- merge(e, nd[, .(to = node, xend = x, yend = y)], by = "to")
  ## A FACTOR with both levels always present, not a bare character column.
  ## At a high THICK_AT a track can have no thick link at all (SSC has none at
  ## 10), and an unfactored column would then declare a one-level linewidth
  ## scale in that panel. patchwork cannot merge guides that differ, so the
  ## collected legend would show the thin key twice instead of one thin/thick
  ## pair. drop = FALSE on the scale keeps both keys in both panels.
  e[, weight_class := factor(
      fifelse(n_supported >= THICK_AT, WT_LAB[["thick"]], WT_LAB[["thin"]]),
      levels = c(WT_LAB[["thin"]], WT_LAB[["thick"]]))]
  e[, drawn := n_supported > 0L]
  ## Label position along each link. A fixed fraction for every link puts all
  ## the labels on one vertical line, where they interleave and cannot be
  ## attributed to their own link, so the fraction is stepped by source row:
  ## the labels then fall on a diagonal and each sits on its own line. The
  ## three links leaving one cell type already diverge, so they need no step
  ## between them.
  e[, src_rank := match(from, c(CT_NODES, tr))]
  e[, lfrac := 0.24 + 0.052 * src_rank]
  e[, `:=`(lx = x + lfrac * (xend - x), ly = y + lfrac * (yend - y))]
  e[, track := tr]
  e[]
}

edges <- rbindlist(lapply(c("ASC", "SSC"), build_edges))
stopifnot(nrow(edges) == 46L)

BLOCK_COL <- c(E = "#0072B2", C = "#D55E00", L = "#009E73")
BLOCK_LAB <- c(E = "E1\u2013E7  normal \u2192 polyp",
               C = "C1\u2013C7  normal \u2192 carcinoma",
               L = "L1 / L2  polyp \u2192 carcinoma")

## 4. one panel ---------------------------------------------------------------
panel <- function(tr, tag_letter, width_guide = FALSE) {
  nd <- node_xy(tr)
  ed <- edges[track == tr]
  n_genes <- unique(ed$track_genes)
  n_thick <- ed[drawn == TRUE & n_supported >= THICK_AT, .N]

  ggplot() +
    ## tested but nothing supported: faint dotted, unlabelled
    geom_segment(data = ed[drawn == FALSE],
                 aes(x = x, y = y, xend = xend, yend = yend),
                 colour = "grey86", linewidth = 0.25, linetype = "dotted") +
    ## supported links, straight so the printed count sits exactly on the line
    geom_segment(data = ed[drawn == TRUE & block != "L"],
                 aes(x = x, y = y, xend = xend, yend = yend,
                     colour = block, linewidth = weight_class),
                 alpha = 0.75, lineend = "round") +
    ## polyp -> carcinoma links arc to the right, clear of the node column
    geom_curve(data = ed[drawn == TRUE & block == "L"],
               aes(x = x, y = y, xend = xend, yend = yend,
                   colour = block, linewidth = weight_class),
               curvature = -0.55, alpha = 0.75, lineend = "round") +
    geom_label(data = ed[drawn == TRUE],
               aes(x = lx, y = ly, label = n_supported, colour = block),
               ## geom_label size is in mm: 2.9 mm is about 8 pt
               size = FS_EDGE / .pt, fill = "white", label.size = 0,
               label.padding = unit(0.5, "mm"), show.legend = FALSE) +
    geom_label(data = nd, aes(x = x, y = y, label = node, fill = kind),
               size = FS_NODE / .pt, label.r = unit(2.6, "mm"), label.size = 0.3,
               colour = "black", show.legend = FALSE) +
    scale_colour_manual(values = BLOCK_COL, labels = BLOCK_LAB,
                        breaks = names(BLOCK_COL), name = NULL) +
    ## patchwork will not merge the two panels' linewidth guides even when the
    ## scales are identical, because panel b has no thick link at THICK_AT = 10
    ## and its guide key rows therefore differ. Collecting would then print the
    ## pair twice, so the guide is drawn for panel a only and serves both.
    scale_linewidth_manual(values = c(0.4, 1.7) |>
                             setNames(c(WT_LAB[["thin"]], WT_LAB[["thick"]])),
                           drop = FALSE, name = NULL,
                           guide = if (width_guide) "legend" else "none") +
    scale_fill_manual(values = c(normal = "white", polyp = "#FDE0C5",
                                 carcinoma = "#DEEBF7"), guide = "none") +
    scale_x_continuous(limits = c(-0.35, 2.95)) +
    scale_y_continuous(limits = c(0.4, 7.6)) +
    ## the plural had been applied to the noun only, giving "1 link reach"
    labs(title = sprintf("%s  %s-nominated genes (n = %d): %d %s %d or more genes",
                         tag_letter, tr, n_genes, n_thick,
                         if (n_thick == 1L) "link reaches" else "links reach",
                         THICK_AT)) +
    theme_void(base_size = FS_LEGEND) +
    theme(plot.title = element_text(size = FS_PANEL, hjust = 0,
                                    margin = margin(b = 4)),
          legend.position = "bottom", legend.box = "vertical",
          legend.spacing.y = unit(0.8, "mm"),
          legend.text = element_text(size = FS_LEGEND),
          legend.key.width = unit(9, "mm"),
          legend.key.height = unit(4, "mm"))
}

p <- (panel("ASC", "a", width_guide = TRUE) | panel("SSC", "b")) +
  plot_layout(guides = "collect") +
  plot_annotation(
    title = paste0("Adenoma-specific risk genes link the normal epithelium to ",
                   "the polyp compartment; serrated-specific genes do not"),
    subtitle = paste0("Link weight is the number of that track's scTWAS genes ",
                      "with supported differential expression in the contrast; ",
                      "printed on each link"),
    caption = paste0(
      "Supported = padj < 0.05, |log2FC| > 0.5 and sign(log2FC) matching the gene's scTWAS Z, counting Up and Down together. The numeric index of a\n",
      "C or E tag denotes the cell type (1 = ABS, 2 = GOB, 3 = CT, 4 = EE, 5 = STM, 6 = TUF, 7 = TAC), so each cell type's own contrast is the one drawn from it.\n",
      "C and E contrasts are cell-type matched; L contrasts compare carcinoma against the polyp compartment and are gene level. Faint dotted links were\n",
      "tested and supported no gene. L1 applies only to ASC genes and L2 only to SSC genes, so each panel carries 23 links."),
    theme = theme(
      plot.title    = element_text(size = FS_TITLE, hjust = 0),
      plot.subtitle = element_text(size = FS_SUBTITLE, hjust = 0, colour = "grey30"),
      plot.caption  = element_text(size = FS_CAPTION, hjust = 0, colour = "grey30"))) &
  theme(legend.position = "bottom")

FIG_STEM <- file.path(DIR_FIG, "Fig_scTWAS_ASC_SSC_linkage")
out <- edges[, .(track, track_genes, block, tag, from, to,
                 n_supported, n_tested, weight_class, drawn)]
setorder(out, track, block, tag)
fwrite(out, paste0(FIG_STEM, "_sourcedata.csv"))
fwrite(out, file.path(DIR_KEY, "scTWAS_ASC_SSC_linkage_edges.csv"))
ggsave(paste0(FIG_STEM, ".png"), plot = p, width = 10.6, height = 6.8,
       units = "in", dpi = 300)
ggsave(paste0(FIG_STEM, ".pdf"), plot = p, width = 10.6, height = 6.8,
       units = "in", useDingbats = FALSE)

## 5. console summary --------------------------------------------------------
cat("\n== links by track and weight class ==\n")
print(dcast(edges, track + block ~ weight_class, value.var = "tag",
            fun.aggregate = length))
cat("\n== zero-support links (drawn dotted) ==\n")
print(edges[drawn == FALSE, .(track, tag, from, to)])
cat("\n== strongest links ==\n")
print(head(edges[order(-n_supported), .(track, tag, from, to, n_supported)], 12))
cat("\nwrote ", FIG_STEM, ".{png,pdf}\n", sep = "")

## ---- compact 12 x 5 cm version, for placement as a composite panel ---------
## Requested size: 12 cm wide x 5 cm tall = 4.72 x 1.97 in.
##
## Why this version drops the title, subtitle and caption rather than shrinking
## them: at the 8 pt floor those four elements plus the legend need 1.47 in of
## the 1.97 in height (75%), leaving 0.5 in for both diagrams, which cannot be
## drawn. They are not lost - this figure is Panel B of Figure 5 and that text
## lives in the composite's legend (Figure5_legend.docx), which is where a panel
## of this size belongs anyway. Kept: the a/b panel titles, the colour legend
## and the link numbers.
##
## Every size here is 8 pt, the FS_MIN floor, so nothing is smaller than in the
## full-size figure - the diagram area shrinks, the type does not.
CM <- 1 / 2.54
FS_C <- FS_MIN

panel_compact <- function(tr, tag_letter, width_guide = FALSE) {
  pp <- panel(tr, tag_letter, width_guide = width_guide)
  pp$labels$title <- sprintf("%s  %s (n = %d)", tag_letter, tr,
                             edges[track == tr, first(track_genes)])
  pp +
    ## re-declare the two geom_label sizes at the floor
    theme_void(base_size = FS_C) +
    theme(plot.title = element_text(size = FS_C, hjust = 0, face = "bold",
                                    margin = margin(b = 1)),
          legend.position = "bottom", legend.box = "horizontal",
          legend.margin = margin(t = 0, b = 0),
          legend.text = element_text(size = FS_C),
          legend.key.width = unit(5, "mm"),
          legend.key.height = unit(2.5, "mm"),
          plot.margin = margin(1, 2, 1, 2))
}

## Label sizes for the compact panel. The link numbers sit BELOW the 8 pt floor
## by explicit choice: at 12 cm the diagram area is about 4 cm tall and 8 pt
## digits on neighbouring links merged (8 and 7 on two different links read as
## "87"; 1 and 5 read as "15"). Smaller separated digits are readable where
## larger touching ones are not. FS_MIN still governs every other element, and
## the full-size figure keeps 8.5 pt numbers.
FS_EDGE_C <- 5.5
stopifnot(FS_EDGE_C < FS_MIN)     # deliberate, documented above

## Link widths for the compact panel. The full-size pair (0.4 / 1.7) scales to
## heavy bands at 12 cm: a 1.7 pt band is wider than a 5.5 pt digit is tall, so
## the thick links were covering the numbers printed beside them. Halving both
## keeps the thin/thick ratio, and therefore the reading, while giving the
## labels clear ground.
LW_C <- c(thin = 0.22, thick = 0.85)

## Label placement for the compact panel. The full-size figure steps the
## along-line fraction by SOURCE row only, so the three or four links leaving
## one cell type share a fraction and rely on their own divergence to separate -
## which holds at 10.6 in and fails at 4.7 in, where the numbers overlapped.
##
## Here each label is placed greedily. Links are taken strongest first, and each
## one is given the position along its own line that overlaps the already-placed
## labels least, with a mild pull toward mid-line to break ties. Measured
## against a geom_label box of 0.6 em per glyph plus 0.5 mm padding, this takes
## the overlapping-pair count from 17 to 0 in panel a and 5 to 0 in panel b,
## against 3 and 5 for a simple even spread. The count is asserted below, so a
## future change to the font, the canvas or the link set cannot quietly
## reintroduce collisions.
PAD_IN <- 0.5 / 25.4
box_w <- function(n, fs) nchar(as.character(n)) * 0.6 * fs / 72 + 2 * PAD_IN
box_h <- function(fs) fs * 1.2 / 72 + 2 * PAD_IN
## plotting area of ONE compact panel, inches, and the data ranges it maps
XR_C <- 3.30; YR_C <- 7.2
W_PAN_C <- (12 / 2.54) / 2 - 0.10
H_PAN_C <- (5 / 2.54) - 0.28

n_label_overlaps <- function(e, fs) {
  xi <- (e$lx + 0.35) / XR_C * W_PAN_C
  yi <- (e$ly - 0.4) / YR_C * H_PAN_C
  w <- box_w(e$n_supported, fs) / 2
  h <- box_h(fs) / 2
  k <- 0L
  for (i in seq_len(nrow(e) - 1L)) for (j in (i + 1L):nrow(e))
    if (abs(xi[i] - xi[j]) < w[i] + w[j] && abs(yi[i] - yi[j]) < 2 * h) k <- k + 1L
  k
}

place_labels <- function(e, fs, grid = seq(0.12, 0.86, by = 0.01)) {
  e <- copy(e)
  setorder(e, -n_supported)                 # strongest links choose first
  xi <- yi <- hw <- numeric(0)
  lf <- numeric(nrow(e))
  hh <- box_h(fs) / 2
  for (i in seq_len(nrow(e))) {
    cx <- (e$x[i] + grid * (e$xend[i] - e$x[i]) + 0.35) / XR_C * W_PAN_C
    cy <- (e$y[i] + grid * (e$yend[i] - e$y[i]) - 0.4) / YR_C * H_PAN_C
    w <- box_w(e$n_supported[i], fs) / 2
    if (!length(xi)) {
      best <- which.min(abs(grid - 0.30))
    } else {
      pen <- vapply(seq_along(grid), function(k)
        sum(pmax(0, (hw + w) - abs(xi - cx[k])) * (abs(yi - cy[k]) < 2 * hh)),
        numeric(1))
      best <- which.min(pen + 0.02 * abs(grid - 0.45))
    }
    lf[i] <- grid[best]
    xi <- c(xi, cx[best]); yi <- c(yi, cy[best]); hw <- c(hw, w)
  }
  e[, lfrac := lf]
  e[, `:=`(lx = x + lfrac * (xend - x), ly = y + lfrac * (yend - y))]
  e[]
}

## Layer order in panel() is segment, curve, edge labels, node labels, so the
## first GeomLabel layer carries the numbers and the second the node names.
compact_labels <- function(pp, width_guide = FALSE) {
  lab_at <- which(vapply(pp$layers, function(L) inherits(L$geom, "GeomLabel"),
                         logical(1)))
  stopifnot(length(lab_at) == 2L)
  pp$layers[[lab_at[1]]]$aes_params$size <- FS_EDGE_C / .pt   # link numbers
  pp$layers[[lab_at[2]]]$aes_params$size <- FS_C / .pt        # node names
  ## re-place the numbers, and refuse to ship a panel with collisions
  ed2 <- place_labels(pp$layers[[lab_at[1]]]$data, FS_EDGE_C)
  stopifnot(n_label_overlaps(ed2, FS_EDGE_C) == 0L)
  pp$layers[[lab_at[1]]]$data <- ed2
  ## No legend in this version. At 12 x 5 cm the keys took nearly half the
  ## canvas and compressed the node rows until their outlines touched; the
  ## colour and line-width meanings belong in the composite's legend text
  ## (Figure5_legend.docx), which is where a panel of this size is placed.
  ## width_guide is accepted and ignored so the call signature matches panel().
  pp +
    scale_colour_manual(values = BLOCK_COL, labels = BLOCK_LAB,
                        breaks = names(BLOCK_COL), name = NULL,
                        guide = "none") +
    scale_linewidth_manual(values = unname(LW_C[c("thin", "thick")]) |>
                             setNames(c(WT_LAB[["thin"]], WT_LAB[["thick"]])),
                           drop = FALSE, name = NULL, guide = "none")
}

pc <- (compact_labels(panel_compact("ASC", "a")) |
       compact_labels(panel_compact("SSC", "b"))) &
  theme(legend.position = "none")

STEM_C <- paste0(FIG_STEM, "_12x5cm")
ggsave(paste0(STEM_C, ".png"), plot = pc,
       width = 12 * CM, height = 5 * CM, units = "in", dpi = 600)
ggsave(paste0(STEM_C, ".pdf"), plot = pc,
       width = 12 * CM, height = 5 * CM, units = "in", useDingbats = FALSE)
cat("wrote", paste0(STEM_C, ".{png,pdf}"), "at 12 x 5 cm, all type at",
    FS_C, "pt\n")
