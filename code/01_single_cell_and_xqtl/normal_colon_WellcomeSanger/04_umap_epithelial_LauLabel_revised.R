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
MAP_FILE <- file.path(FIG, "fine_annotation_mapping.csv")

# Same 05-style palette used by 09.R's active match05 setting.
COLORS <- c(STM = "#00B6EB", TAC = "#A58AFF", ABS = "#F8766D",
            CT = "#C49A00", GOB = "#00C094", EE = "#53B400",
            TUF = "#FB61D7", ASC = "#1F3B73", SSC = "#4F5311")
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
               metric = "euclidean", seed.use = 0, verbose = FALSE)

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
  x <- trimws(gsub("^[[:punct:][:space:]]+|[[:punct:][:space:]]+$", "", x))
  ifelse(nzchar(x), x, "Unspecified")
}
map <- df %>% count(LauLabel, original_annotation, name = "n_cells") %>%
  arrange(factor(LauLabel, levels = names(COLORS)), desc(n_cells), original_annotation) %>%
  mutate(simplified_annotation = make_short(original_annotation))
if (file.exists(MAP_FILE)) {
  old <- read.csv(MAP_FILE, check.names = FALSE, stringsAsFactors = FALSE)
  required <- c("LauLabel", "original_annotation", "simplified_annotation")
  if (!all(required %in% names(old))) stop("Mapping CSV missing required columns: ", MAP_FILE)
  key <- paste(map$LauLabel, map$original_annotation, sep = "\r")
  oldkey <- paste(old$LauLabel, old$original_annotation, sep = "\r")
  if (anyDuplicated(oldkey)) stop("Duplicate original annotation keys in mapping CSV")
  ix <- match(key, oldkey)
  edited <- !is.na(ix) & !is.na(old$simplified_annotation[ix]) & nzchar(trimws(old$simplified_annotation[ix]))
  map$simplified_annotation[edited] <- trimws(old$simplified_annotation[ix[edited]])
}
# Write the mapping before drawing panels 04b/04c, retaining any edits.
map$display_label <- paste0(map$LauLabel, ": ", map$simplified_annotation)
write.csv(map, MAP_FILE, row.names = FALSE)
if (anyDuplicated(map$display_label)) stop("Duplicate display_label in mapping. Edit simplified_annotation in: ", MAP_FILE)
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
p_c <- panel("fine_display", fine_colors, 2.5)
save_plot(p_a, paste0("04a_umap_clusters_", tag))
save_plot(p_b, paste0("04b_umap_LauLabel_", tag))
save_plot(p_c, paste0("04c_umap_finelabel_", tag))

# Same three-panel composite; 04c identities are now actual on-cluster labels.
if (requireNamespace("patchwork", quietly = TRUE)) {
  small <- function(p, title) p + labs(title = title) +
    theme(plot.title = element_text(size = 8), axis.title = element_text(size = 7),
          axis.text = element_text(size = 6))
  pp <- small(p_a, "a  Louvain clusters") |
        small(p_b, "b  LauLabel") |
        small(p_c, "c  Simplified fine annotation")
  save_plot(pp, "umap_WellcomeSanger_epithelial_clusters_LauLabel", 12, 4.5)
} else warning("Install patchwork to obtain composite; all single panels were saved")
writeLines(capture.output(sessionInfo()), file.path(FIG, "04_umap_R_sessionInfo.txt"))
message("Wrote mapping: ", MAP_FILE, " | cell colors: ", file.path(FIG, "cell_type_colors.csv"))
message("Cells: ", nrow(df), " | clusters: ", nrow(call), " | ARI: ", round(ari, 3))
