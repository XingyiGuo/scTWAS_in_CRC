#!/usr/bin/env Rscript
# =============================================================================
# UMAP of the cells that built the HTAN ADSER scTWAS models (9 cell types)
# -----------------------------------------------------------------------------
# Runs on the cluster:
#   script  /data/l2_bioinfo1/liq17/scTWAS/data/VUMC_HTAN/VUMC_HTAN_DIS_VAL_ADSER/src
#   data    /data/l2_bioinfo1/liq17/scTWAS/data/VUMC_HTAN/VUMC_HTAN_DIS_VAL_ADSER/Results
# and unchanged on the Windows mirror (L:/scRNA-seq/VUMC_HTAN/...); ROOT is a
# first-existing-wins candidate list, so no edit is needed to move between them.
# Submit with 08a_umap_ADSER_model_cells_sbatch.sh.
#
# Input object : objH_DISVAL_PCA83_EPI_afterCell200_Gene3QC_beforeConditionFilter.rds
#   Merged DIS + VAL + PCA83 epithelium, saved at line 715 of
#   HTAN_ADSER_CellOracleRsave_mergeQC_PCA83_HPtoSER_ASCSSCfilter.R, i.e. AFTER
#   per-cell QC (nFeature >= 200, applied per h5ad before the merge) and gene QC
#   (detected in >= 3 cells, applied after the merge) and BEFORE any disease
#   condition / ASC-SSC lineage filtering. Counts only: no normalisation, no
#   variable features, no PCA, no UMAP. It is the only object on L: that still
#   contains the three PCA83 add-on specimens (HTA11_6407_2000001021,
#   HTA11_8890_2000001011, HTA11_8930_2000001011 = 3,175 model cells); the
#   12:25 CellOracle outputs predate them.
#
# Cell set    : the 81,929 cells listed in HTAN_ADSER_model_cells_81929_keys.csv,
#   which is the cell-level TWAS match table (9CT_expression_id_matched_cell_
#   metadata.csv) keyed by merged_cell_id = paste0(htan_file, "_", h5ad_barcode)
#   -- the same identifier the merged object uses for its columns.
#
# Outputs     : UMAP coordinate CSV (cell_id + UMAPHARMONY_1/2 + annotation),
#   per-cell-type count table, PDF + PNG figure, optional processed RDS.
#
# Resources   : ~82k cells x 2,000 scaled HVGs is ~1.3 GB dense plus the sparse
#   counts; budget 32 GB and 8 cores. Do NOT run this on a 15 GB laptop with the
#   full object loaded -- subset immediately after readRDS (this script does).
# =============================================================================

# ---- 0. Thread hygiene (must precede library loads) -------------------------
# A threaded/dense allocation of ~323 GiB killed the first full-panel run
# ("libgomp: Out of memory allocating 346711037184 bytes"). The request scales
# with the OpenMP/BLAS thread count, so the compiled layers are pinned to one
# thread here and parallelism is left to the explicit future plan below. Also
# export OMP_NUM_THREADS=1 in the sbatch wrapper -- libraries that read it at
# dlopen time will not see a value set after they are loaded.
for (v in c("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
            "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS")) {
  if (!nzchar(Sys.getenv(v))) do.call(Sys.setenv, setNames(list("1"), v))
}

suppressPackageStartupMessages({
  library(Seurat)
  library(harmony)
  library(ggplot2)
  library(dplyr)
  library(readr)
})

# ---- 1. Config --------------------------------------------------------------
# First existing wins: cluster path, then the Windows mirror.
ROOT_CANDIDATES <- c(
  "/data/l2_bioinfo1/liq17/scTWAS/data/VUMC_HTAN/VUMC_HTAN_DIS_VAL_ADSER",
  "L:/scRNA-seq/VUMC_HTAN/VUMC_HTAN_DIS_VAL_ADSER"
)
ROOT <- ROOT_CANDIDATES[which(dir.exists(ROOT_CANDIDATES))[1]]
if (is.na(ROOT)) stop("No ROOT candidate exists: ", paste(ROOT_CANDIDATES, collapse = " | "))
RES_DIR <- file.path(ROOT, "Results")
message("ROOT: ", ROOT)

OBJ_RDS  <- file.path(RES_DIR, "objH_DISVAL_PCA83_EPI_afterCell200_Gene3QC_beforeConditionFilter.rds")
OUT_DIR  <- file.path(RES_DIR, "UMAP_ADSER_model_cells")

# Model-cell key set, in order of preference:
#   1. the prepared key file (merged_cell_id already built), if it was copied over
#   2. the cell-level TWAS match table, which is always present on both sides --
#      merged_cell_id is then rebuilt here exactly as the merge script built it.
KEYS_CANDIDATES <- c(
  file.path(RES_DIR, "TWAS_scRNA_meta", "HTAN_ADSER_model_cells_81929_keys.csv"),
  file.path(RES_DIR, "TWAS_scRNA_meta", "9CT_expression_id_matched_cell_metadata.csv")
)
N_MODEL_CELLS_EXPECTED <- 81929

SEED           <- 20260923
N_HVG          <- 2000
N_PC           <- 30
HARMONY_GROUP  <- "SampleID"   # donor-level correction; "batch" is the coarser alternative
UMAP_MIN_DIST  <- 0.30
UMAP_NN        <- 30
SAVE_OBJECT    <- TRUE         # processed Seurat with PCA/harmony/UMAP
ALLOW_PARTIAL  <- FALSE        # see the missing-cell guard in section 2

# Cells drawn per cell type, sampled without replacement within each type after
# the model-cell subset. 4000 caps the five large types and takes every cell of
# CT (1,317), EE (456), TUF (778); NA embeds all 81,929. The cap is what makes
# the run fit in memory, and it also stops ABS (34,144 cells, 42% of the panel)
# from dominating the HVGs and the PCs.
CELLS_PER_TYPE <- 4000
REFILTER_GENES <- TRUE         # drop genes detected in <3 of the retained cells
USE_FUTURE     <- FALSE        # TRUE only if a run is CPU-bound, not memory-bound

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
set.seed(SEED)

n_cpus <- suppressWarnings(as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "1")))
if (is.na(n_cpus) || n_cpus < 1) n_cpus <- 1
if (USE_FUTURE && n_cpus > 1 && requireNamespace("future", quietly = TRUE)) {
  future::plan("multicore", workers = min(n_cpus, 4))
  options(future.globals.maxSize = 8 * 1024^3)
  message("future: multicore x ", min(n_cpus, 4))
} else {
  if (requireNamespace("future", quietly = TRUE)) future::plan("sequential")
  message("future: sequential (memory-safe default)")
}

# Step timer: prints elapsed seconds and current R heap so the next failure is
# localised instead of appearing as a bare allocation error.
.t0 <- Sys.time()
step <- function(label) {
  gc_now <- gc(verbose = FALSE)
  message(sprintf("[%6.1fs] %-34s heap %5.1f GB",
                  as.numeric(difftime(Sys.time(), .t0, units = "secs")),
                  label, sum(gc_now[, "used"] * c(8, 8)) / 1024^2 / 1024))
}

# Canonical project palette (CT_SIG_COLORS, NPG-derived, darkened to >=4.5:1 on
# white) so this embedding threads with the Manhattan panels in W:/.../Viz.
CT_COLORS <- c(
  STM = "#4433A8",  # indigo
  TAC = "#17695E",  # teal
  ABS = "#B5382A",  # red
  CT  = "#A76A00",  # amber
  GOB = "#17557F",  # blue
  EE  = "#6E4B3A",  # brown
  TUF = "#6E6E22",  # olive
  ASC = "#8A5015",  # ochre
  SSC = "#882255"   # wine
)
# NOTE: the usual allowance that CT/EE/TUF may reuse the ASC/SSC hues does NOT
# apply here -- this panel draws all nine types at once, so all nine hues must
# stay distinct. The block above is already hue-separated; do not collapse it.

CT_NAMES <- c(
  STM = "Stem Cell Type",
  TAC = "Transit-Amplifying Cell Type",
  ABS = "Absorptive Cell Type",
  CT  = "Crypt-Top Cell Type",
  GOB = "Goblet Cell Type",
  EE  = "Enteroendocrine Cell Type",
  TUF = "Tuft Cell Type",
  ASC = "Adenoma-Specific Cell Type",
  SSC = "Serrated-Polyp-Specific Cell Type"
)
CT_ORDER <- names(CT_COLORS)   # crypt order, lesion-specific types last

# ---- 2. Load and subset to the model cells ----------------------------------
message("Reading merged object: ", OBJ_RDS)
obj  <- readRDS(OBJ_RDS)
message("  full object: ", nrow(obj), " genes x ", ncol(obj), " cells")

KEYS_CSV <- KEYS_CANDIDATES[which(file.exists(KEYS_CANDIDATES))[1]]
if (is.na(KEYS_CSV)) stop("No key file found: ", paste(KEYS_CANDIDATES, collapse = " | "))
message("Key file: ", KEYS_CSV)
keys <- readr::read_csv(KEYS_CSV, show_col_types = FALSE)

# Rebuild merged_cell_id when reading the raw match table. This is the identifier
# the merge script assigned: merged_cell_id <- paste0(nm, "_", obs$h5ad_barcode)
# where nm is the h5ad file stem, recorded per cell in htan_file.
if (!"merged_cell_id" %in% names(keys)) {
  stopifnot(all(c("htan_file", "h5ad_barcode") %in% names(keys)))
  keys$merged_cell_id <- paste0(keys$htan_file, "_", keys$h5ad_barcode)
}
stopifnot(!any(duplicated(keys$merged_cell_id)))
if (nrow(keys) != N_MODEL_CELLS_EXPECTED) {
  warning("key file has ", nrow(keys), " rows, expected ", N_MODEL_CELLS_EXPECTED,
          " -- check that the match table has not been regenerated")
}
message("  key file: ", nrow(keys), " model cells")

present <- keys$merged_cell_id %in% colnames(obj)
if (!all(present)) {
  miss <- keys[!present, c("SampleID", "expression_celltype")]
  readr::write_csv(miss, file.path(OUT_DIR, "QC_model_cells_absent_from_object.csv"))
  msg <- paste0(sum(!present), " of ", nrow(keys), " model cells are not in this ",
                "object (see QC_model_cells_absent_from_object.csv). Expected cause: ",
                "a filtered CellOracle object was used -- AllCT_ASCSSCfiltered and the ",
                "nine by_celltype RDS were written before the three PCA83 add-on ",
                "specimens (6407/8890/8930, 3,175 cells) were merged in, so they hold ",
                "only 78,754 of the 81,929. Use the beforeConditionFilter object.")
  if (!ALLOW_PARTIAL) stop(msg)
  warning(msg, " Continuing on the present cells because ALLOW_PARTIAL = TRUE.")
  keys <- keys[present, ]
}

# Stratified sample BEFORE touching the matrix, so the object is subset once.
keys_use <- keys
if (!is.na(CELLS_PER_TYPE)) {
  set.seed(SEED)
  # slice_sample(n = ) returns the whole group when it holds fewer than n rows,
  # so no per-group min() is needed (and `n =` must be a constant anyway).
  keys_use <- keys |>
    group_by(expression_celltype) |>
    slice_sample(n = CELLS_PER_TYPE, replace = FALSE) |>
    ungroup()
  samp_tab <- keys |>
    count(expression_celltype, name = "available") |>
    left_join(count(keys_use, expression_celltype, name = "sampled"),
              by = "expression_celltype")
  readr::write_csv(samp_tab, file.path(OUT_DIR, "sampling_by_celltype.csv"))
  message("Stratified sample: ", nrow(keys_use), " of ", nrow(keys),
          " model cells, <= ", CELLS_PER_TYPE, " per cell type")
  print(as.data.frame(samp_tab), row.names = FALSE)
}

obj <- subset(obj, cells = keys_use$merged_cell_id)
stopifnot(ncol(obj) == nrow(keys_use))
keys <- keys_use
invisible(gc(verbose = FALSE))
step("subset to model cells")

# The merged object's gene set was filtered across all 196,672 cells; in this
# subset thousands of genes are all-zero and only cost memory downstream.
if (REFILTER_GENES) {
  det <- Matrix::rowSums(GetAssayData(obj, layer = "counts") > 0)
  keep_genes <- names(det)[det >= 3]
  message("Gene refilter: ", length(keep_genes), " of ", nrow(obj),
          " genes detected in >=3 retained cells")
  obj <- subset(obj, features = keep_genes)
  invisible(gc(verbose = FALSE))
  step("gene refilter")
}

# Annotation carried over from the TWAS match table. CellType_cb agrees with
# expression_celltype for all 81,929 cells; Cell_Type is NA for every PCA83
# cell, so it is never used for colouring.
idx <- match(colnames(obj), keys$merged_cell_id)
obj$CellType_model    <- factor(keys$expression_celltype[idx], levels = CT_ORDER)
obj$SampleID          <- keys$SampleID[idx]
obj$PatientID         <- keys$PatientID[idx]
obj$cohort_model      <- keys$cohort[idx]
obj$batch_model       <- keys$batch[idx]
obj$condition_model   <- keys$sample_condition[idx]
stopifnot(!any(is.na(obj$CellType_model)))

readr::write_csv(
  as.data.frame(table(CellType = obj$CellType_model, Cohort = obj$cohort_model)),
  file.path(OUT_DIR, "cells_by_celltype_cohort.csv")
)

# ---- 3. Embedding ------------------------------------------------------------
obj <- NormalizeData(obj, normalization.method = "LogNormalize", scale.factor = 1e4, verbose = FALSE)
step("NormalizeData")
obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = N_HVG, verbose = FALSE)
step("FindVariableFeatures")
# Scale and PCA on the HVGs only -- never `features = rownames(obj)`, which is
# the dense-allocation trap on a panel this size.
obj <- ScaleData(obj, features = VariableFeatures(obj), verbose = FALSE)
step("ScaleData (HVGs only)")
obj <- RunPCA(obj, features = VariableFeatures(obj), npcs = N_PC, seed.use = SEED, verbose = FALSE)
step("RunPCA")

# Donor-level batch correction. harmony >= 1.0 uses reduction.use/reduction.save;
# older wrappers used reduction=. Try the current names, fall back once.
obj <- tryCatch(
  harmony::RunHarmony(obj, group.by.vars = HARMONY_GROUP, reduction.use = "pca",
                      dims.use = 1:N_PC, reduction.save = "harmony", verbose = FALSE),
  error = function(e) harmony::RunHarmony(obj, group.by.vars = HARMONY_GROUP,
                                          reduction = "pca", dims.use = 1:N_PC, verbose = FALSE)
)

step("RunHarmony")

obj <- RunUMAP(obj, reduction = "harmony", dims = 1:N_PC,
               min.dist = UMAP_MIN_DIST, n.neighbors = UMAP_NN,
               seed.use = SEED, reduction.name = "umapharmony",
               reduction.key = "UMAPHARMONY_", verbose = FALSE)
step("RunUMAP")

emb <- as.data.frame(Embeddings(obj, "umapharmony"))
names(emb) <- c("UMAPHARMONY_1", "UMAPHARMONY_2")
coords <- data.frame(
  cell_id        = rownames(emb),
  emb,
  CellType       = as.character(obj$CellType_model),
  SampleID       = obj$SampleID,
  PatientID      = obj$PatientID,
  cohort         = obj$cohort_model,
  sample_condition = obj$condition_model,
  row.names = NULL
)
COORD_CSV <- file.path(OUT_DIR, "UMAP_ADSER_model_cells_coords.csv")
readr::write_csv(coords, COORD_CSV)
message("Wrote coordinates: ", COORD_CSV)
# Column names (cell_id, UMAPHARMONY_1/2) match the PANEL_B_CSV candidate list
# in Final_Viz_4panel.R, so this drops into the composite without an edit.

if (SAVE_OBJECT) {
  saveRDS(obj, file.path(OUT_DIR, "objH_ADSER_model_cells_umap.rds"), compress = TRUE)
}

# ---- 4. Figure ---------------------------------------------------------------
n_cells <- nrow(coords)
pt_size <- if (n_cells >= 2e5) 0.035 else if (n_cells >= 5e4) 0.08 else 0.22

set.seed(SEED)
plot_df <- coords[sample(nrow(coords)), ]            # shuffled draw order
plot_df$CellType <- factor(plot_df$CellType, levels = CT_ORDER)

cent <- plot_df |>
  group_by(CellType) |>
  summarise(x = median(UMAPHARMONY_1), y = median(UMAPHARMONY_2), .groups = "drop")

rng_x <- range(plot_df$UMAPHARMONY_1); rng_y <- range(plot_df$UMAPHARMONY_2)
pad   <- 0.04 * diff(rng_x)
ax0x  <- rng_x[1] - pad; ax0y <- rng_y[1] - pad
axlen <- 0.14 * diff(rng_x)

label_layer <- if (requireNamespace("shadowtext", quietly = TRUE)) {
  shadowtext::geom_shadowtext(data = cent, aes(x, y, label = CellType, colour = CellType),
                              bg.colour = "white", bg.r = 0.18, fontface = "bold",
                              size = 3.0, show.legend = FALSE)
} else {
  list(
    geom_text(data = cent, aes(x, y, label = CellType), colour = "white",
              fontface = "bold", size = 3.4, show.legend = FALSE),
    geom_text(data = cent, aes(x, y, label = CellType, colour = CellType),
              fontface = "bold", size = 3.0, show.legend = FALSE)
  )
}

p <- ggplot(plot_df, aes(UMAPHARMONY_1, UMAPHARMONY_2, colour = CellType)) +
  geom_point(size = pt_size, stroke = 0, alpha = 0.85) +
  label_layer +
  scale_colour_manual(
    values = CT_COLORS,
    breaks = CT_ORDER,
    labels = paste0(CT_ORDER, " \u2014 ", CT_NAMES[CT_ORDER]),
    name   = NULL
  ) +
  annotate("segment", x = ax0x, xend = ax0x + axlen, y = ax0y, yend = ax0y,
           arrow = arrow(length = unit(0.06, "in")), linewidth = 0.3, colour = "grey30") +
  annotate("segment", x = ax0x, xend = ax0x, y = ax0y, yend = ax0y + axlen,
           arrow = arrow(length = unit(0.06, "in")), linewidth = 0.3, colour = "grey30") +
  annotate("text", x = ax0x + axlen / 2, y = ax0y - 0.35 * pad, label = "UMAP 1",
           size = 2.1, colour = "grey30", vjust = 1) +
  annotate("text", x = ax0x - 0.35 * pad, y = ax0y + axlen / 2, label = "UMAP 2",
           size = 2.1, colour = "grey30", angle = 90, vjust = 0) +
  coord_fixed(ratio = 1, clip = "off") +
  guides(colour = guide_legend(override.aes = list(size = 2.2, alpha = 1))) +
  labs(caption = sprintf(
    "%s cells%s from %d specimens / %d participants; Harmony on %s, %d PCs",
    format(n_cells, big.mark = ","),
    if (is.na(CELLS_PER_TYPE)) "" else sprintf(
      " sampled from %s model cells at \u2264%s per cell type",
      format(N_MODEL_CELLS_EXPECTED, big.mark = ","),
      format(CELLS_PER_TYPE, big.mark = ",")),
    dplyr::n_distinct(plot_df$SampleID), dplyr::n_distinct(plot_df$PatientID),
    HARMONY_GROUP, N_PC)) +
  theme_void(base_size = 8) +
  theme(
    plot.caption      = element_text(size = 5.6, colour = "grey35", hjust = 0),
    legend.position   = "right",
    legend.key        = element_blank(),
    legend.background = element_blank(),
    legend.text       = element_text(size = 6.5),
    legend.key.height = unit(0.34, "cm"),
    plot.margin       = margin(6, 6, 10, 10)
  )

FIG_BASE <- file.path(OUT_DIR, "UMAP_ADSER_model_cells_celltype")
# cairo_pdf keeps the em-dash in the legend; fall back to the base pdf device on
# an R build without cairo.
pdf_dev <- if (isTRUE(capabilities("cairo"))) grDevices::cairo_pdf else grDevices::pdf
ggsave(paste0(FIG_BASE, ".pdf"), plot = p, width = 7.4, height = 4.8, device = pdf_dev)
ggsave(paste0(FIG_BASE, ".png"), plot = p, width = 7.4, height = 4.8, dpi = 300)
message("Wrote figure: ", FIG_BASE, ".{pdf,png}")

# Companion panel: cohort, to show the three cohorts are mixed after Harmony.
p_coh <- ggplot(plot_df[sample(nrow(plot_df)), ],
                aes(UMAPHARMONY_1, UMAPHARMONY_2, colour = cohort)) +
  geom_point(size = pt_size, stroke = 0, alpha = 0.8) +
  scale_colour_manual(values = c(DIS = "#1B7837", VAL = "#762A83", PCA83 = "#C0C0C0"), name = NULL) +
  coord_fixed(ratio = 1) +
  guides(colour = guide_legend(override.aes = list(size = 2.2, alpha = 1))) +
  theme_void(base_size = 8) +
  theme(legend.position = "right", legend.key = element_blank())
ggsave(file.path(OUT_DIR, "UMAP_ADSER_model_cells_cohort.png"),
       plot = p_coh, width = 6.4, height = 4.4, dpi = 300)

writeLines(capture.output(sessionInfo()), file.path(OUT_DIR, "sessionInfo.txt"))
message("Done. n cells = ", n_cells, " | point size = ", pt_size)
