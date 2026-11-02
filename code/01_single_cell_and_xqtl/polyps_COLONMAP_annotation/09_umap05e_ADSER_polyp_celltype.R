#!/usr/bin/env Rscript
# =============================================================================
# 09_umap05e_ADSER_polyp_celltype.R
# -----------------------------------------------------------------------------
# 05e-style UMAP for the VUMC HTAN AD/SER polyp panel (DIS + VAL + PCA83), drawn
# over the cells that built the 9 scTWAS cell-type models.
#
# The figure call is the one used for
#   VUMC_HTAN_GSE266546/src/05_annotation_UMAP.R -> 05e_umap_reviewed_celltype*
# i.e.  DimPlot(reduction = "umap_harmony", group.by = <annotation>,
#               label = TRUE, raster = FALSE) + NoLegend()
#       ggsave(..., width = 8, height = 6, device = cairo_pdf)
# The only difference is the annotation column: GSE266546 labels clusters
# (celltype_annot_reviewed); here the label is the cell type that the TWAS
# expression matrices were built from (expression_celltype -> CellType_model),
# so the panel shows exactly the nine model cell types and nothing else.
#
# WHERE THE CELL TYPES COME FROM (see the two scripts named below for detail)
#   HTAN_ADSER_CellOracleRsave_mergeQC_PCA83_HPtoSER_ASCSSCfilter.R
#     - cell QC  : >= 200 detected genes, applied per h5ad before the merge
#     - gene QC  : detected in >= 3 cells, applied after DIS+VAL+PCA83 merge
#     - DIS/VAL  : Sample_Classification canonicalised to AD / SER, and
#                  CellType_cb restricted to ABS GOB CT EE TUF TAC STM ASC SSC
#     - PCA83    : TISSUE == "Polyp" only, then ADTYPE Tubular/Tubulovillous
#                  -> AD, ADTYPE Sessile/Traditional Serrated -> SER, and only
#                  where ADTYPE is missing or "Not stated" HPTYPE
#                  Microvesicular/Goblet-cell rich -> SER. Unclassified polyps
#                  are dropped. All non-missing PCA83 CellType_cb labels kept.
#     - lineage  : ASC retained only in AD tissue, SSC only in SER tissue
#   The TWAS panel is then the genotype-matched subset of that object, listed
#   cell by cell in TWAS_scRNA_meta/HTAN_ADSER_model_cells_81929_keys.csv
#   (81,929 cells / 92 specimens / 84 participants).
#
# INPUTS (first existing wins, cluster path then the Windows mirror)
#   Results/UMAP_ADSER_model_cells/objH_ADSER_model_cells_umap.rds
#       -> reuse the embedding written by 08_umap_ADSER_model_cells.R
#   Results/objH_DISVAL_PCA83_EPI_afterCell200_Gene3QC_beforeConditionFilter.rds
#       -> otherwise compute the embedding here with 08's parameters
#   Results/TWAS_scRNA_meta/HTAN_ADSER_model_cells_81929_keys.csv
#       (or 9CT_expression_id_matched_cell_metadata.csv, from which
#        merged_cell_id is rebuilt as paste0(htan_file, "_", h5ad_barcode))
#
# OUTPUTS -> Results/UMAP_ADSER_model_cells/figures/
#   05e_umap_ADSER_polyp_celltype.pdf / .png     <- the cell-type panel
#   05d_umap_ADSER_polyp_cluster_res0.8.pdf      <- cluster-id review panel
#   05e_umap_ADSER_polyp_cohort.pdf              <- Harmony mixing check
#   05e_umap_ADSER_polyp_condition.pdf           <- AD vs SER
#   ../UMAP_ADSER_polyp_coords.csv               <- cell_id + UMAPHARMONY_1/2
#                                                   (+ cluster_res0.8)
#   ../UMAP_ADSER_polyp_counts_by_celltype.csv
#   ../UMAP_ADSER_polyp_ASCSSC_lineage_check.csv
#   ../UMAP_ADSER_polyp_cluster_res0.8_celltype_composition.csv
#
# The cell-type panel identifies types in a LEGEND, not with text drawn on the
# clusters (CELLTYPE_LABELS_ON_PLOT). The cluster-id panel keeps its numbers on
# the plot, as 05d does. Clustering is additive: it never relabels a cell.
#
# NOTE on reuse: clustering is recomputed every run because it is cheap, but it
# is only persisted into PROCESSED_RDS when this script computes the embedding.
# A run that reuses a saved embedding re-derives the partition from it, which
# is deterministic at a fixed SEED, N_PC and CLUSTER_RES.
#
# SUBMIT   sbatch 09a_umap05e_ADSER_polyp_celltype_sbatch.sh
#
# READING THE LOGS: every progress line here comes from message(), which R
# writes to STDERR, while print() of the count tables goes to STDOUT. So the
# .error file carries the whole narrative (timings, paths, the lineage check)
# and a NON-EMPTY .error is normal, not a failure. Judge success by the
# "exit=0" line the wrapper echoes and by the three "Wrote:" lines.
#
# Note on PROCESSED_RDS: this script writes the SAME filename 08 writes,
# objH_ADSER_model_cells_umap.rds, because at equal CELLS_PER_TYPE the two
# embeddings are the same recipe over the same cells. Whichever runs second
# overwrites it, and the cell-identity check in section 3 is what keeps a
# mismatched object from being reused rather than the filename.
# =============================================================================

# ---- 0. Thread hygiene (must precede library loads) -------------------------
# A threaded dense allocation of ~323 GiB killed the first full-panel run of 08
# ("libgomp: Out of memory allocating 346711037184 bytes"). The request scales
# with the OpenMP/BLAS thread count, so the compiled layers are pinned to one
# thread here; the sbatch wrapper exports the same variables because some
# builds read them only at dlopen time.
THREAD_VARS <- c("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
                 "OMP_THREAD_LIMIT", "VECLIB_MAXIMUM_THREADS",
                 "NUMEXPR_NUM_THREADS", "RCPP_PARALLEL_NUM_THREADS")
omp_env_before <- Sys.getenv("OMP_NUM_THREADS", unset = "")
for (v in THREAD_VARS) {
  if (!nzchar(Sys.getenv(v))) do.call(Sys.setenv, setNames(list("1"), v))
}

options(bitmapType = "cairo")     # as in 05_annotation_UMAP.R
set.seed(123)

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(readr)
})

# ---- 0b. Clamp threads at RUNTIME, not just in the environment --------------
# Setting OMP_NUM_THREADS above is necessary but NOT sufficient: libgomp reads
# its internal control variables in its own library constructor, i.e. when it
# is dlopen'd, which happens before this script's Sys.setenv can run if the
# wrapper did not already export the variable. A run launched without the
# exports therefore still asks libgomp for one thread per core on the node, the
# cgroup refuses, and the job dies at RunHarmony with
#   libgomp: Thread creation failed: Resource temporarily unavailable
# (EAGAIN from pthread_create -- a thread/process LIMIT, not a memory problem;
# it is a different failure from the ~323 GiB dense allocation that killed the
# first full-panel run of 08). RhpcBLASctl calls omp_set_num_threads() and the
# OpenBLAS setter through their C APIs, which do take effect after load, so it
# is the only in-process fix. Preflight rather than discover this after a
# 904 MB readRDS and a PCA.
# Escape hatch: set TRUE only if you know the node's thread limits are fine.
ALLOW_UNPINNED_THREADS <- FALSE

clamped_at_runtime <- FALSE
if (requireNamespace("RhpcBLASctl", quietly = TRUE)) {
  RhpcBLASctl::blas_set_num_threads(1)
  RhpcBLASctl::omp_set_num_threads(1)
  clamped_at_runtime <- TRUE
}
message("Thread control:")
message("  OMP_NUM_THREADS in the environment at startup: ",
        if (nzchar(omp_env_before)) omp_env_before else "<unset>")
message("  RhpcBLASctl runtime clamp: ",
        if (clamped_at_runtime) "applied (BLAS 1, OMP 1)" else "NOT AVAILABLE")
message("  cores visible: ", parallel::detectCores(),
        " | SLURM_CPUS_PER_TASK: ",
        Sys.getenv("SLURM_CPUS_PER_TASK", unset = "<unset>"))
if (clamped_at_runtime) {
  # blas_get_num_procs() probes OpenBLAS directly and returns an uninitialised
  # value on a FlexiBLAS-dispatched R build -- ACCRE's r/4.5.0 loads
  # flexiblas/3.3.1 over openblas + blis, and the readback came back as
  # 1596085528. The setter still took effect; only this readback is unreliable,
  # so an implausible number is reported as such instead of being printed as if
  # it meant something.
  nb <- tryCatch(RhpcBLASctl::blas_get_num_procs(), error = function(e) NA_integer_)
  message("  effective: BLAS ",
          if (is.na(nb) || nb < 1 || nb > 1024) {
            paste0("readback unreliable on this BLAS (got ", nb,
                   "); FlexiBLAS/BLIS dispatch does not answer the ",
                   "OpenBLAS query -- the OMP figure below is the one to trust")
          } else nb,
          " | OMP ", RhpcBLASctl::omp_get_max_threads())
}
if (!identical(omp_env_before, "1") && !clamped_at_runtime &&
    !ALLOW_UNPINNED_THREADS) {
  stop("OMP_NUM_THREADS was not exported before R started (it is '",
       if (nzchar(omp_env_before)) omp_env_before else "<unset>",
       "') and RhpcBLASctl is not installed, so the compiled layers cannot be ",
       "pinned to one thread from inside R. This run would very likely die at ",
       "RunHarmony with 'libgomp: Thread creation failed: Resource ",
       "temporarily unavailable'. Fix either way:\n",
       "  (a) submit through 09a_umap05e_ADSER_polyp_celltype_sbatch.sh, which ",
       "exports OMP_NUM_THREADS=1 before Rscript; or run\n",
       "      OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 Rscript --vanilla ",
       "09_umap05e_ADSER_polyp_celltype.R\n",
       "  (b) install RhpcBLASctl once: ",
       "install.packages('RhpcBLASctl')\n",
       "Set ALLOW_UNPINNED_THREADS <- TRUE above to override this check.")
}

# The thread limit bites at whichever OpenMP region the run reaches first, so
# it is NOT a Harmony-specific fault: one run died at RunHarmony after PCA, the
# next died before PCA in the normalise/scale block. Every compiled stage is
# therefore wrapped by with_thread_diag(), which names the stage and explains
# the failure class instead of letting a bare libgomp line be the whole story.
stop_if_thread_failure <- function(e, where) {
  msg <- conditionMessage(e)
  if (!grepl("Thread creation failed|Resource temporarily unavailable|libgomp|pthread",
             msg, ignore.case = TRUE)) {
    return(invisible(NULL))
  }
  stop("Thread/process limit hit during ", where, ": ", msg, "\n",
       "This is EAGAIN from pthread_create -- the job cannot create another ",
       "thread. It is a LIMIT, not a shortage of RAM, so raising --mem will ",
       "not help (this run's heap was under 2 GB). Checklist:\n",
       "  1. Submit via 09a_umap05e_ADSER_polyp_celltype_sbatch.sh so ",
       "OMP_NUM_THREADS=1 is exported BEFORE R starts -- libgomp reads it in ",
       "its constructor and ignores a later Sys.setenv. The wrapper also ",
       "echoes 'ulimit -u' and the cgroup pids.max; read those from the .out.\n",
       "  2. Install RhpcBLASctl so threads can be clamped after load (this ",
       "run: ", if (clamped_at_runtime) "clamp applied" else "NOT AVAILABLE",
       ", OMP_NUM_THREADS at startup '",
       if (nzchar(omp_env_before)) omp_env_before else "<unset>", "').\n",
       "  3. If the clamp WAS applied and it still fails, the job's thread ",
       "allowance is genuinely exhausted: raise --cpus-per-task, or ask ACCRE ",
       "support about ulimit -u / pids.max for the partition.\n",
       "  4. R-level parallelism is not used here (no future/parallel plan), ",
       "so any remaining threads come from OpenBLAS/OpenMP inside Seurat, ",
       "RcppArmadillo or harmony.",
       call. = FALSE)
}
with_thread_diag <- function(where, expr) {
  tryCatch(expr, error = function(e) {
    stop_if_thread_failure(e, where)
    stop(e)
  })
}

# ---- 1. Config --------------------------------------------------------------
ROOT_CANDIDATES <- c(
  "/data/l2_bioinfo1/liq17/scTWAS/data/VUMC_HTAN/VUMC_HTAN_DIS_VAL_ADSER",
  "L:/scRNA-seq/VUMC_HTAN/VUMC_HTAN_DIS_VAL_ADSER"
)
ROOT <- ROOT_CANDIDATES[which(dir.exists(ROOT_CANDIDATES))[1]]
if (is.na(ROOT)) {
  stop("No ROOT candidate exists: ", paste(ROOT_CANDIDATES, collapse = " | "))
}
RES_DIR <- file.path(ROOT, "Results")
OUT_DIR <- file.path(RES_DIR, "UMAP_ADSER_model_cells")
FIG_DIR <- file.path(OUT_DIR, "figures")
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)
message("ROOT   : ", ROOT)
message("FIG_DIR: ", FIG_DIR)

PROCESSED_RDS <- file.path(OUT_DIR, "objH_ADSER_model_cells_umap.rds")
MERGED_RDS <- file.path(
  RES_DIR, "objH_DISVAL_PCA83_EPI_afterCell200_Gene3QC_beforeConditionFilter.rds"
)
KEYS_CANDIDATES <- c(
  file.path(RES_DIR, "TWAS_scRNA_meta", "HTAN_ADSER_model_cells_81929_keys.csv"),
  file.path(RES_DIR, "TWAS_scRNA_meta", "9CT_expression_id_matched_cell_metadata.csv")
)
N_MODEL_CELLS_EXPECTED <- 81929

# Embedding parameters, identical to 08_umap_ADSER_model_cells.R so the two
# panels are the same embedding when this script has to compute it.
SEED          <- 20260923
N_HVG         <- 2000
N_PC          <- 30
HARMONY_GROUP <- "SampleID"     # donor level; "batch" is the coarser option
UMAP_MIN_DIST <- 0.30
UMAP_NN       <- 30
REFILTER_GENES <- TRUE

# Cells drawn per cell type, sampled without replacement within each type after
# the model-cell subset. 4000 is enough to draw the panel and is the same cap
# 08_umap_ADSER_model_cells.R settled on: it takes 26,551 of the 81,929 cells
# (ABS/ASC/GOB/SSC/STM/TAC capped at 4,000; CT 1,317, TUF 778, EE 456 taken
# whole) while still covering all 92 specimens and 84 participants. The cap is
# also what stops ABS -- 34,144 cells, 42% of the panel -- from dominating the
# HVGs and the PCs. Set NA to embed every model cell.
CELLS_PER_TYPE <- 4000L

# Reuse an embedding already on disk instead of recomputing it. Set FALSE to
# force a fresh embedding (it will overwrite PROCESSED_RDS).
REUSE_EMBEDDING <- TRUE
SAVE_OBJECT     <- TRUE

# Take the cell set -- and therefore the specimens and participants -- from
# 08_umap_ADSER_model_cells.R's own output rather than re-drawing it, so 08 and
# 09 are the same cells rather than two draws that merely share a seed. See
# section 2b. Set FALSE to draw independently at SEED / CELLS_PER_TYPE.
CELLS_FROM_08 <- TRUE

# Cells with no AD/SER label (the three PCA83 add-on specimens) are kept by
# default because the TWAS models were fitted on all 81,929 cells. See 4b.
DROP_UNLABELLED_CONDITION <- FALSE

MAKE_EXTRA_PANELS <- TRUE       # cohort + condition companion UMAPs

# How the cell-type panel names its types.
#   "leader"  each type is named OUTSIDE the cloud in a left/right column and
#             joined to its own cells by a thin black line ending on a
#             colour-filled dot -- the layout settled on for 04c in
#             04_umap_epithelial_LauLabel_04c_filtered.R. Nothing is typed over
#             the points and no key has to be cross-referenced.
#   "legend"  the previous behaviour: no text on the clusters, colours read off
#             an inside or right-hand key.
#   "on-plot" the original 05e look: label = TRUE + NoLegend().
CELLTYPE_LABEL_MODE <- "leader"

# Kept for backward compatibility: TRUE still forces the original 05e look.
CELLTYPE_LABELS_ON_PLOT <- FALSE
if (CELLTYPE_LABELS_ON_PLOT) CELLTYPE_LABEL_MODE <- "on-plot"
if (!CELLTYPE_LABEL_MODE %in% c("leader", "legend", "on-plot")) {
  stop("CELLTYPE_LABEL_MODE must be leader, legend or on-plot; got: ",
       CELLTYPE_LABEL_MODE)
}

# Figure size for the cell-type panel. "leader" needs horizontal room for the
# label columns -- the margin is solved from the real text width, so these two
# numbers are passed to the labeller as well as to save_plot() and must agree.
CT_PANEL_W <- if (identical(CELLTYPE_LABEL_MODE, "leader")) 11 else 8
CT_PANEL_H <- if (identical(CELLTYPE_LABEL_MODE, "leader")) 7 else 6
# FALSE = bare codes (ABS, CT, ...), as 05_annotation_UMAP.R shows them. TRUE
# expands each entry to "ABS - Absorptive Cell Type".
LEGEND_FULL_NAMES <- FALSE

# Legend placement. "topleft" puts it INSIDE the panel at the upper-left
# corner. Note the attached 05_annotation_UMAP.R actually uses "right"; set
# this to "right" for literal parity with that script.
LEGEND_PLACEMENT <- "topleft"   # "topleft" | "right"
LEGEND_TITLE     <- "Cell type" # as labs(colour = ) in 05; NULL for no title

# Axis titles. DimPlot derives them from the reduction KEY, so this embedding
# would otherwise be labelled "UMAPHARMONY_1" / "UMAPHARMONY_2". These are the
# titles Final_Viz_4panel.R uses for its UMAP panels, so the panels read the
# same way side by side. Applied in save_plot(), which every panel goes
# through, so no panel can drift from the others.
AXIS_X <- "UMAP 1"
AXIS_Y <- "UMAP 2"

# De-novo clustering of these same cells, at the resolution 04_sub_cluster.R
# uses for the GSE266546 object. Additive only: it never relabels anything --
# the panel is always coloured by the TWAS cell-type label. See section 4c.
DO_CLUSTER  <- TRUE
CLUSTER_RES <- 0.8
# The cluster-ID review panel keeps numbers drawn on the clusters, because the
# numbers are what make that panel usable -- the legend convention above is for
# the cell-type panels only, exactly as 05d does.
CLUSTER_PANEL_LABELS <- TRUE

CT_ORDER <- c("STM", "TAC", "ABS", "CT", "GOB", "EE", "TUF", "ASC", "SSC")
CT_NAMES <- c(                  # canonical display names for figure panels
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

# ---- Palette ----------------------------------------------------------------
#   "match05"  each cell type gets the colour 05_annotation_UMAP.R gives it
#   "project"  CT_SIG_COLORS, as used by the Manhattan/coloc figures in
#              W:/VUMC_projects/scTWAS/Viz
#   "seurat"   no manual scale; DimPlot's own hue_pal() over these 9 levels.
#              NOTE this does NOT reproduce 05's colours -- see below.
PALETTE <- "match05"

# Why "seurat" is not the same as matching 05. DimPlot colours identities with
# scales::hue_pal()(n) walked along the factor levels, so every colour depends
# on n. 05 has SEVEN levels (ABS CT EE GOB STM TAC TUF) and this panel has NINE
# -- the same seven plus ASC and SSC -- so the ring is resampled and only ABS
# happens to land on the same hue:
#   05 (n=7):  ABS #F8766D  CT #C49A00  EE #53B400  GOB #00C094  STM #00B6EB
#              TAC #A58AFF  TUF #FB61D7
#   here (n=9) if left to DimPlot: CT #93AA00, EE #00BA38, GOB #00C19F,
#              STM #619CFF, TAC #DB72FB, TUF #FF61C3 -- six of seven differ.
# So "match05" recomputes hue_pal() at 05's OWN level count and assigns by
# NAME. It is computed rather than hardcoded so it tracks whatever scales
# version the cluster has, the same way 05 does.
PROJECT_COLORS <- c(            # NPG-derived, darkened to >= 4.5:1 on white
  STM = "#4433A8", TAC = "#17695E", ABS = "#B5382A", CT  = "#A76A00",
  GOB = "#17557F", EE  = "#6E4B3A", TUF = "#6E6E22", ASC = "#8A5015",
  SSC = "#882255"
)

# ASC and SSC under PALETTE = "match05". These are NOT the project ochre
# (#8A5015) / wine (#882255) pair: those two both sit in the dark red-brown
# family and read as the same colour at UMAP point size. Distances in CIE Lab,
# where the bar to clear is 42.4 -- the minimum separation WITHIN 05's own
# seven colours:
#   ochre + wine : 55.5 apart, 40.6 from the nearest of 05's seven
#   navy + olive : 74.5 apart, 48.9 from the nearest of 05's seven
# Chosen from a 20-colour dark-candidate search maximising the worse of those
# two numbers. Opposite hue families and both much darker than 05's light
# pastels, so they separate by hue AND by lightness. Swap the two values if
# adenoma reads better as the green one.
LESION_COLORS <- c(ASC = "#1F3B73",   # dark blue
                   SSC = "#4F5311")   # dark olive

# 05's identity levels. Read from its own output when reachable so this is data
# rather than a constant: if that label set ever changes, n changes and ALL
# seven colours shift, which is exactly the failure this block exists to catch.
CT05_LEVELS_RECORDED <- c("ABS", "CT", "EE", "GOB", "STM", "TAC", "TUF")
CT05_CSV_CANDIDATES <- c(
  "/data/l2_bioinfo1/liq17/scTWAS/data/VUMC_HTAN_GSE266546/raw/rds2/cluster_reviewed_annotation_res0.8.rm23_sub.csv",
  "L:/scRNA-seq/VUMC_HTAN_GSE266546/raw/rds2/cluster_reviewed_annotation_res0.8.rm23_sub.csv"
)

build_match05_palette <- function() {
  hit <- CT05_CSV_CANDIDATES[which(file.exists(CT05_CSV_CANDIDATES))[1]]
  if (!is.na(hit)) {
    rev_annot <- readr::read_csv(hit, show_col_types = FALSE)
    lv <- sort(unique(as.character(rev_annot$reviewed_annotation)))
    message("05 identity levels read from ", basename(hit), ": ",
            paste(lv, collapse = ", "))
  } else {
    lv <- CT05_LEVELS_RECORDED
    message("05's annotation table not found; using the recorded level set: ",
            paste(lv, collapse = ", "))
  }
  if (!setequal(lv, CT05_LEVELS_RECORDED)) {
    stop("05's identity levels are now {", paste(lv, collapse = ", "),
         "} but this palette was matched against {",
         paste(CT05_LEVELS_RECORDED, collapse = ", "), "}. Because hue_pal() ",
         "depends on the level COUNT, every shared colour has shifted and the ",
         "two panels no longer agree. Re-derive the mapping (and update ",
         "CT05_LEVELS_RECORDED) before trusting this figure.")
  }
  # Seurat factors identities in sorted order, so this reproduces 05 exactly.
  base <- setNames(scales::hue_pal()(length(lv)), lv)
  # ASC and SSC do not exist in 05, so they need hues of their own -- see
  # LESION_COLORS above for why they are NOT the project ochre/wine pair.
  extra <- LESION_COLORS[c("ASC", "SSC")]
  out <- c(base, extra)[CT_ORDER]
  stopifnot(!any(is.na(out)), !any(duplicated(out)))
  out
}

CT_COLORS <- switch(
  PALETTE,
  match05 = build_match05_palette(),
  project = PROJECT_COLORS[CT_ORDER],
  seurat  = NULL,
  stop("PALETTE must be one of: match05, project, seurat")
)
if (!is.null(CT_COLORS)) {
  message("Palette (", PALETTE, "): ",
          paste(names(CT_COLORS), unname(CT_COLORS), sep = "=", collapse = " "))
}

.t0 <- Sys.time()
step <- function(label) {
  g <- gc(verbose = FALSE)
  message(sprintf("[%6.1fs] %-34s heap %5.1f GB",
                  as.numeric(difftime(Sys.time(), .t0, units = "secs")),
                  label, sum(g[, "used"] * c(8, 8)) / 1024^2 / 1024))
}

# ---- 2. Model-cell key table ------------------------------------------------
KEYS_CSV <- KEYS_CANDIDATES[which(file.exists(KEYS_CANDIDATES))[1]]
if (is.na(KEYS_CSV)) {
  stop("No key file found: ", paste(KEYS_CANDIDATES, collapse = " | "))
}
message("Key file: ", KEYS_CSV)
keys <- readr::read_csv(KEYS_CSV, show_col_types = FALSE)
if (!"merged_cell_id" %in% names(keys)) {
  stopifnot(all(c("htan_file", "h5ad_barcode") %in% names(keys)))
  keys$merged_cell_id <- paste0(keys$htan_file, "_", keys$h5ad_barcode)
}
stopifnot(!any(duplicated(keys$merged_cell_id)))
if (nrow(keys) != N_MODEL_CELLS_EXPECTED) {
  warning("key file has ", nrow(keys), " rows, expected ",
          N_MODEL_CELLS_EXPECTED, " -- has the match table been regenerated?")
}
message("Model cells in key table: ", nrow(keys))
print(table(keys$expression_celltype, keys$cohort, useNA = "ifany"))

# ---- 2b. Cell set: the same individuals and cells as 08 ---------------------
# 08_umap_ADSER_model_cells.R writes one row per embedded cell to
# UMAP_ADSER_model_cells_coords.csv, so its exact draw is recoverable. Reading
# that cell_id column is what makes 08 and 09 the same panel: it pins the
# specimens and participants to 08's, and does not depend on dplyr's RNG stream
# being identical across R/dplyr versions the way re-running slice_sample at a
# shared seed would. The seeded draw below is only the fallback for a clean
# tree where 08 has not run yet.
COORDS_08 <- file.path(OUT_DIR, "UMAP_ADSER_model_cells_coords.csv")

if (CELLS_FROM_08 && file.exists(COORDS_08)) {
  c08 <- readr::read_csv(COORDS_08, show_col_types = FALSE)
  if (!"cell_id" %in% names(c08)) {
    stop("08's coordinate file has no cell_id column: ", COORDS_08)
  }
  ids08 <- unique(as.character(c08$cell_id))
  absent <- setdiff(ids08, keys$merged_cell_id)
  if (length(absent)) {
    stop(length(absent), " of 08's ", length(ids08), " cells are absent from ",
         "the key table. 08 and the key table describe different cell sets -- ",
         "the match table was probably regenerated after 08 ran.")
  }
  keys_use <- keys[match(ids08, keys$merged_cell_id), ]
  CELL_SET_SOURCE <- paste0("08 coords (", basename(COORDS_08), ")")
} else if (!is.na(CELLS_PER_TYPE)) {
  if (CELLS_FROM_08) {
    message("CELLS_FROM_08 is TRUE but ", COORDS_08, " does not exist; ",
            "falling back to the seeded stratified draw.")
  }
  set.seed(SEED)
  # slice_sample(n =) returns the whole group when it holds fewer than n rows,
  # so no per-group min() is needed.
  keys_use <- keys |>
    group_by(expression_celltype) |>
    slice_sample(n = CELLS_PER_TYPE, replace = FALSE) |>
    ungroup()
  readr::write_csv(
    keys |>
      count(expression_celltype, name = "available") |>
      left_join(count(keys_use, expression_celltype, name = "sampled"),
                by = "expression_celltype"),
    file.path(OUT_DIR, "sampling_by_celltype_09.csv")
  )
  CELL_SET_SOURCE <- paste0("seeded stratified draw, <=", CELLS_PER_TYPE,
                            " per cell type, seed ", SEED)
} else {
  keys_use <- keys
  CELL_SET_SOURCE <- "all model cells"
}

message("\nCell set: ", CELL_SET_SOURCE)
message("  cells: ", nrow(keys_use),
        " | specimens: ", dplyr::n_distinct(keys_use$SampleID),
        " | participants: ", dplyr::n_distinct(keys_use$PatientID))
print(table(keys_use$expression_celltype, keys_use$cohort, useNA = "ifany"))

# ---- 3. Get an embedded object ----------------------------------------------
# 05e's DimPlot names the reduction "umap_harmony"; 08 wrote "umapharmony".
# Accept either and normalise, so this script works on both.
normalise_reduction <- function(obj) {
  have <- Reductions(obj)
  hit <- intersect(c("umap_harmony", "umapharmony", "umap"), have)[1]
  if (is.na(hit)) {
    stop("No UMAP reduction on the object; found: ", paste(have, collapse = ", "))
  }
  if (hit != "umap_harmony") {
    obj[["umap_harmony"]] <- obj[[hit]]
    message("Renamed reduction '", hit, "' -> 'umap_harmony'")
  }
  obj
}

# The cell set resolved in 2b is what any saved embedding is checked against,
# so an object built over different cells -- a full 81,929-cell run, or a run at
# a different CELLS_PER_TYPE -- is not silently reused under this one.
n_expected <- nrow(keys_use)

reuse_ok <- FALSE
if (REUSE_EMBEDDING && file.exists(PROCESSED_RDS)) {
  obj <- readRDS(PROCESSED_RDS)
  if (ncol(obj) == n_expected &&
      !length(setdiff(colnames(obj), keys_use$merged_cell_id))) {
    obj <- normalise_reduction(obj)
    reuse_ok <- TRUE
    message("Reusing embedding: ", PROCESSED_RDS)
    step("readRDS (processed)")
  } else {
    message("Saved embedding at ", PROCESSED_RDS, " holds ", ncol(obj),
            " cells and this run's set has ", n_expected,
            " (", length(setdiff(colnames(obj), keys_use$merged_cell_id)),
            " of its cells are outside that set) -- recomputing instead of ",
            "reusing it.")
    rm(obj)
    invisible(gc(verbose = FALSE))
  }
}

if (!reuse_ok) {
  if (!file.exists(MERGED_RDS)) {
    stop("Neither a processed embedding nor the merged object exists:\n  ",
         PROCESSED_RDS, "\n  ", MERGED_RDS)
  }
  message("Computing the embedding from: ", MERGED_RDS)
  obj <- readRDS(MERGED_RDS)
  message("  full object: ", nrow(obj), " genes x ", ncol(obj), " cells")
  step("readRDS (merged)")

  present <- keys$merged_cell_id %in% colnames(obj)
  if (!all(present)) {
    readr::write_csv(
      keys[!present, intersect(c("SampleID", "expression_celltype"), names(keys))],
      file.path(OUT_DIR, "QC_model_cells_absent_from_object.csv")
    )
    stop(sum(!present), " of ", nrow(keys), " model cells are not in this ",
         "object (see QC_model_cells_absent_from_object.csv). Expected cause: ",
         "a filtered CellOracle object was used -- AllCT_ASCSSCfiltered and the ",
         "nine by_celltype RDS were written before the three PCA83 add-on ",
         "specimens (HTA11_6407/8890/8930, 3,175 cells) were merged in. Use ",
         "the beforeConditionFilter object.")
  }

  # Subset to the cell set resolved in 2b before touching the matrix.
  obj <- subset(obj, cells = keys_use$merged_cell_id)
  stopifnot(ncol(obj) == nrow(keys_use))
  invisible(gc(verbose = FALSE))
  step("subset to model cells")

  if (REFILTER_GENES) {
    # The merged gene set was filtered across all 196,672 cells; in this subset
    # thousands of genes are all-zero and only cost memory downstream.
    det <- Matrix::rowSums(GetAssayData(obj, layer = "counts") > 0)
    keep_genes <- names(det)[det >= 3]
    message("Gene refilter: ", length(keep_genes), " of ", nrow(obj),
            " genes detected in >=3 retained cells")
    obj <- subset(obj, features = keep_genes)
    invisible(gc(verbose = FALSE))
    step("gene refilter")
  }

  # Each stage is timed AND wrapped, so a thread failure names the stage it
  # died in rather than appearing as a bare libgomp line after the last step().
  obj <- with_thread_diag("NormalizeData",
    NormalizeData(obj, normalization.method = "LogNormalize",
                  scale.factor = 1e4, verbose = FALSE))
  step("NormalizeData")
  obj <- with_thread_diag("FindVariableFeatures",
    FindVariableFeatures(obj, selection.method = "vst",
                         nfeatures = N_HVG, verbose = FALSE))
  step("FindVariableFeatures")
  # Scale and PCA on the HVGs only -- never features = rownames(obj), which is
  # the dense-allocation trap on a panel this size.
  obj <- with_thread_diag("ScaleData",
    ScaleData(obj, features = VariableFeatures(obj), verbose = FALSE))
  step("ScaleData (HVGs only)")
  obj <- with_thread_diag("RunPCA",
    RunPCA(obj, features = VariableFeatures(obj), npcs = N_PC,
           seed.use = SEED, verbose = FALSE))
  step("RunPCA")

  # RunHarmony is where the OpenMP layer is first exercised hard, so it is where
  # an unpinned thread count surfaces. Re-clamp immediately before the call --
  # a package loaded in between (harmony pulls in RcppArmadillo) can reset the
  # OpenMP thread count when it initialises.
  if (clamped_at_runtime) {
    RhpcBLASctl::blas_set_num_threads(1)
    RhpcBLASctl::omp_set_num_threads(1)
    message("Re-clamped before Harmony: BLAS ", RhpcBLASctl::blas_get_num_procs(),
            " | OMP ", RhpcBLASctl::omp_get_max_threads())
  }

  # harmony >= 1.0 uses reduction.use/reduction.save; older wrappers used
  # reduction=. Try the current names, fall back once on an ARGUMENT error --
  # but never swallow a thread/resource failure, which the blind fallback in
  # earlier versions of this script turned into a second identical crash.
  run_harmony <- function(use_new_args) {
    if (use_new_args) {
      harmony::RunHarmony(obj, group.by.vars = HARMONY_GROUP,
                          reduction.use = "pca", dims.use = 1:N_PC,
                          reduction.save = "harmony", verbose = FALSE)
    } else {
      harmony::RunHarmony(obj, group.by.vars = HARMONY_GROUP,
                          reduction = "pca", dims.use = 1:N_PC, verbose = FALSE)
    }
  }
  obj <- tryCatch(
    run_harmony(TRUE),
    error = function(e) {
      # A thread failure must never be retried as if it were an argument-name
      # problem: the blind fallback in the first version of this script turned
      # it into a second identical crash.
      stop_if_thread_failure(e, "RunHarmony")
      message("RunHarmony rejected the current argument names (",
              conditionMessage(e), "); retrying with the pre-1.0 names.")
      with_thread_diag("RunHarmony (pre-1.0 argument names)", run_harmony(FALSE))
    }
  )
  step("RunHarmony")

  obj <- with_thread_diag("RunUMAP",
    RunUMAP(obj, reduction = "harmony", dims = 1:N_PC,
            min.dist = UMAP_MIN_DIST, n.neighbors = UMAP_NN,
            seed.use = SEED, reduction.name = "umap_harmony",
            reduction.key = "UMAPHARMONY_", verbose = FALSE))
  step("RunUMAP")

  if (SAVE_OBJECT) {
    saveRDS(obj, PROCESSED_RDS, compress = TRUE)
    message("Saved processed object: ", PROCESSED_RDS)
  }
}

# ---- 4. Annotation column ---------------------------------------------------
# expression_celltype is the label the TWAS expression matrices were built on.
# CellType_cb agrees with it for all 81,929 cells; Cell_Type is NA for every
# PCA83 cell, so it is never used here. Match on barcode, never row position.
idx <- match(colnames(obj), keys$merged_cell_id)
if (any(is.na(idx))) {
  stop(sum(is.na(idx)), " cells on the object are absent from the key table; ",
       "the embedding and the key table describe different cell sets.")
}
obj$CellType_model  <- factor(keys$expression_celltype[idx], levels = CT_ORDER)
obj$SampleID        <- keys$SampleID[idx]
obj$PatientID       <- keys$PatientID[idx]
obj$cohort_model    <- keys$cohort[idx]
obj$condition_model <- keys$sample_condition[idx]
stopifnot(!any(is.na(obj$CellType_model)))

# ---- 4b. ASC / SSC lineage rule --------------------------------------------
# ASC is the adenoma-specific cell type and SSC the serrated-polyp-specific
# one, so ASC must only ever be drawn from AD tissue and SSC only from SER
# tissue. The rule is applied upstream in
# HTAN_ADSER_CellOracleRsave_mergeQC_PCA83_HPtoSER_ASCSSCfilter.R
# (drop_asc_nonad / drop_ssc_nonser, audited in
#  HTAN_DISVAL_PCA83_ADSER_HPtoSER_ASCSSC_lineage_drop_AUDIT.csv: DIS 441 ASC /
#  113 SSC, VAL 314 / 93, PCA83 507 / 826 cells dropped), and the panel here
# inherits it through the key table. It is re-proved rather than assumed,
# because a regenerated match table would break it silently.
ct_chr <- as.character(obj$CellType_model)
cond_chr <- as.character(obj$condition_model)
bad_asc <- which(ct_chr == "ASC" & (is.na(cond_chr) | cond_chr != "AD"))
bad_ssc <- which(ct_chr == "SSC" & (is.na(cond_chr) | cond_chr != "SER"))
lineage_audit <- data.frame(
  CellType        = c("ASC", "SSC"),
  required_tissue = c("AD (adenoma)", "SER (serrated polyp)"),
  cells           = c(sum(ct_chr == "ASC"), sum(ct_chr == "SSC")),
  violations      = c(length(bad_asc), length(bad_ssc)),
  stringsAsFactors = FALSE
)
readr::write_csv(lineage_audit,
                 file.path(OUT_DIR, "UMAP_ADSER_polyp_ASCSSC_lineage_check.csv"))
message("\nASC/SSC lineage check:")
print(lineage_audit, row.names = FALSE)
if (length(bad_asc) || length(bad_ssc)) {
  readr::write_csv(
    data.frame(
      cell_id  = colnames(obj)[c(bad_asc, bad_ssc)],
      CellType = ct_chr[c(bad_asc, bad_ssc)],
      sample_condition = cond_chr[c(bad_asc, bad_ssc)],
      SampleID = obj$SampleID[c(bad_asc, bad_ssc)]
    ),
    file.path(OUT_DIR, "QC_ASCSSC_lineage_violations.csv")
  )
  stop(length(bad_asc), " ASC cell(s) outside AD tissue and ", length(bad_ssc),
       " SSC cell(s) outside SER tissue are present in the key table (see ",
       "QC_ASCSSC_lineage_violations.csv). The ASC/SSC lineage filter has been ",
       "lost upstream -- rerun the merge script before drawing this panel.")
}

# Cells whose AD/SER label is missing entirely. These are the three PCA83
# add-on specimens (HTA11_6407/8890/8930), whose ADTYPE/HPTYPE is NA, Unknown
# or "Not stated", so the pathology mapping never classified them. They carry
# no ASC or SSC cells -- only the shared types -- so they cannot violate the
# rule above. They stay in the panel by default because the TWAS models were
# fitted on all 81,929 cells; set DROP_UNLABELLED_CONDITION to exclude them.
n_unlab <- sum(is.na(cond_chr))
if (n_unlab) {
  message(n_unlab, " cell(s) have no AD/SER label, from specimen(s): ",
          paste(sort(unique(obj$SampleID[is.na(cond_chr)])), collapse = ", "))
  if (DROP_UNLABELLED_CONDITION) {
    obj <- subset(obj, cells = colnames(obj)[!is.na(cond_chr)])
    message("Dropped them: ", ncol(obj), " cells remain")
    ct_chr <- as.character(obj$CellType_model)
    cond_chr <- as.character(obj$condition_model)
  }
}
# DimPlot drops cells whose group.by value is NA, so the condition panel uses a
# column where the missing label is an explicit level.
obj$condition_plot <- ifelse(is.na(cond_chr), "unlabelled", cond_chr)
Idents(obj) <- "CellType_model"

message("\nCells per model cell type:")
print(table(obj$CellType_model))
message("cohort:");    print(table(obj$cohort_model, useNA = "ifany"))
message("condition:"); print(table(obj$condition_model, useNA = "ifany"))
message("specimens: ", dplyr::n_distinct(obj$SampleID),
        " | participants: ", dplyr::n_distinct(obj$PatientID))

# ---- 4c. De-novo clustering at resolution 0.8 -------------------------------
# This panel does NOT need clustering to exist: it is coloured by
# expression_celltype, the label the TWAS models were fitted on, so the groups
# are given rather than discovered. The clustering here is additive -- an
# independent partition of the same cells, at the resolution 04_sub_cluster.R
# uses for the GSE266546 object (FindNeighbors on the Harmony embedding over
# the same dims, then FindClusters at 0.8) -- and its value is the
# cross-tabulation against the model labels: it shows whether the nine TWAS
# cell types are recovered as distinct neighbourhoods, and which of them a
# graph-based clustering would merge or split. It never overrides the labels.
if (DO_CLUSTER) {
  cl_red <- if ("harmony" %in% Reductions(obj)) "harmony" else "pca"
  if (cl_red != "harmony") {
    message("No 'harmony' reduction on this object; clustering on '", cl_red,
            "' instead, so the partition is NOT batch-corrected.")
  }
  obj <- with_thread_diag("FindNeighbors",
    FindNeighbors(obj, reduction = cl_red, dims = 1:N_PC, verbose = FALSE))
  obj <- with_thread_diag("FindClusters",
    FindClusters(obj, resolution = CLUSTER_RES, random.seed = SEED,
                 verbose = FALSE))
  step(paste0("FindClusters (res ", CLUSTER_RES, ")"))

  # FindClusters moves Idents to the new partition; the cell-type label is the
  # identity this script plots, so put it back.
  cl <- as.character(Idents(obj))
  obj$cluster <- factor(cl, levels = as.character(sort(as.integer(unique(cl)))))
  Idents(obj) <- "CellType_model"
  message("Clusters at resolution ", CLUSTER_RES, ": ",
          dplyr::n_distinct(obj$cluster))
  print(table(obj$cluster))

  # Composition of each cluster in terms of the model cell types, with the
  # dominant type flagged -- the same quantity 05_annotation_UMAP.R annotates
  # from, reported here rather than used to relabel anything.
  comp <- data.frame(cluster = obj$cluster,
                     CellType = as.character(obj$CellType_model)) |>
    count(cluster, CellType, name = "cells") |>
    group_by(cluster) |>
    mutate(cluster_total = sum(cells),
           pct_of_cluster = round(100 * cells / cluster_total, 2),
           is_dominant = pct_of_cluster == max(pct_of_cluster)) |>
    ungroup() |>
    arrange(cluster, desc(pct_of_cluster))
  readr::write_csv(comp, file.path(
    OUT_DIR, paste0("UMAP_ADSER_polyp_cluster_res", CLUSTER_RES,
                    "_celltype_composition.csv")))
  dom <- comp |> filter(is_dominant) |>
    select(cluster, dominant_celltype = CellType, cells, cluster_total,
           pct_of_cluster)
  message("\nDominant model cell type per cluster (res ", CLUSTER_RES, "):")
  print(as.data.frame(dom), row.names = FALSE)

  # How many clusters each cell type is spread across, and vice versa: the two
  # numbers that say whether the labels and the partition agree.
  message("\nClusters per cell type:")
  print(comp |> group_by(CellType) |>
          summarise(clusters = n_distinct(cluster),
                    clusters_where_dominant = sum(is_dominant),
                    .groups = "drop") |> as.data.frame(), row.names = FALSE)
}

# ---- 5. Coordinates + count tables -----------------------------------------
emb <- as.data.frame(Embeddings(obj, "umap_harmony"))
names(emb) <- c("UMAPHARMONY_1", "UMAPHARMONY_2")
coords <- data.frame(
  cell_id          = rownames(emb),
  emb,
  CellType         = as.character(obj$CellType_model),
  SampleID         = obj$SampleID,
  PatientID        = obj$PatientID,
  cohort           = obj$cohort_model,
  sample_condition = obj$condition_model,
  row.names        = NULL
)
if (DO_CLUSTER) {
  coords[[paste0("cluster_res", CLUSTER_RES)]] <- as.character(obj$cluster)
}
readr::write_csv(coords, file.path(OUT_DIR, "UMAP_ADSER_polyp_coords.csv"))
# Column names (cell_id, UMAPHARMONY_1/2) match the PANEL_B_CSV candidate list
# in Final_Viz_4panel.R, so this drops into the composite without an edit.

counts_tab <- coords |>
  group_by(CellType) |>
  summarise(cells = n(),
            specimens    = n_distinct(SampleID),
            participants = n_distinct(PatientID),
            AD_cells     = sum(sample_condition == "AD", na.rm = TRUE),
            SER_cells    = sum(sample_condition == "SER", na.rm = TRUE),
            .groups = "drop") |>
  mutate(CellType = factor(CellType, levels = CT_ORDER)) |>
  arrange(CellType)
readr::write_csv(counts_tab,
                 file.path(OUT_DIR, "UMAP_ADSER_polyp_counts_by_celltype.csv"))
message("\nCounts by cell type:")
print(as.data.frame(counts_tab), row.names = FALSE)

# ---- 6. 05e figure ----------------------------------------------------------
# ---- Cell-type label geometry: side columns hugging the cloud ---------------
# Ported VERBATIM from 04_umap_epithelial_LauLabel_04c_filtered.R (which took
# the densest-bin anchor rule from Final_Viz_4panel.R) so the two panels place
# their labels identically. Both functions read UMAP_1/UMAP_2 off the data
# frame handed to them, so this embedding's UMAPHARMONY_1/2 columns are copied
# under those names at the call site rather than the helpers being forked.
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

# Same helper and geometry as save_plot() in 05_annotation_UMAP.R.
save_plot <- function(p, name, w = 8, h = 6) {
  # Every panel in this script is a UMAP of the same embedding, so the axis
  # titles are set here rather than at four call sites -- one place to change,
  # and a panel added later cannot forget them. Overriding DimPlot's key-derived
  # "UMAPHARMONY_1"/"UMAPHARMONY_2".
  # DimPlot titles the panel with its group.by column, so the cell-type panel
  # came out headed "CellType_model" and the cluster panel "cluster". Dropped
  # here, with the axis titles, so no panel can reintroduce it.
  p <- p + labs(x = AXIS_X, y = AXIS_Y, title = NULL)
  ggsave(file.path(FIG_DIR, paste0(name, ".pdf")), plot = p,
         width = w, height = h, device = cairo_pdf)
  ggsave(file.path(FIG_DIR, paste0(name, ".png")), plot = p,
         width = w, height = h, dpi = 300)
  message("Wrote: ", file.path(FIG_DIR, name), ".{pdf,png}")
}

# The text each type is named by, shared between the legend and the leader
# labels so the two modes cannot disagree.
CT_LABEL_TEXT <- if (LEGEND_FULL_NAMES) {
  setNames(paste0(CT_ORDER, " \u2014 ", CT_NAMES[CT_ORDER]), CT_ORDER)
} else {
  setNames(CT_ORDER, CT_ORDER)
}

p05e <- DimPlot(
  obj,
  reduction = "umap_harmony",
  group.by  = "CellType_model",
  label     = identical(CELLTYPE_LABEL_MODE, "on-plot"),
  raster    = FALSE
)
if (identical(CELLTYPE_LABEL_MODE, "on-plot")) {
  p05e <- p05e + NoLegend()          # original 05e: text on clusters, no key
} else if (identical(CELLTYPE_LABEL_MODE, "leader")) {
  p05e <- p05e + NoLegend()          # the labels ARE the key
} else {
  # Legend instead of on-plot text, with 05_annotation_UMAP.R's own legend
  # settings: labs(colour = "Cell type"), override.aes size 3, text 10 /
  # title 11. override.aes is needed because the plotted cells are drawn at
  # DimPlot's point size, which is unreadable as a key glyph.
  leg_theme <- if (identical(LEGEND_PLACEMENT, "topleft")) {
    # Inside the panel, anchored by its own top-left corner so the block grows
    # down and right and cannot overhang the axes.
    theme(legend.position = c(0.01, 0.99),
          legend.justification = c(0, 1),
          legend.background = element_rect(fill = alpha("white", 0.7),
                                           colour = NA),
          legend.key = element_blank(),
          legend.margin = margin(2, 4, 2, 2),
          legend.spacing.y = unit(0.02, "cm"))
  } else {
    theme(legend.position = "right", legend.key = element_blank())
  }
  p05e <- p05e +
    labs(colour = LEGEND_TITLE) +
    guides(colour = guide_legend(override.aes = list(size = 3, alpha = 1))) +
    leg_theme +
    theme(legend.text = element_text(size = 10),
          legend.title = element_text(size = 11))
}
if (!is.null(CT_COLORS)) {
  # Assign by NAME, so each cell type keeps its colour whatever order or subset
  # of types a given panel happens to contain. `breaks` fixes the legend order
  # to crypt order with the two lesion-specific types last, rather than the
  # alphabetical order DimPlot would use.
  p05e <- p05e + scale_colour_manual(
    values = CT_COLORS,
    breaks = CT_ORDER,
    labels = if (LEGEND_FULL_NAMES) {
      paste0(CT_ORDER, " \u2014 ", CT_NAMES[CT_ORDER])
    } else CT_ORDER
  )
} else if (LEGEND_FULL_NAMES) {
  # PALETTE = "seurat": no manual scale to carry labels, so relabel via the
  # discrete scale while leaving DimPlot's own colours untouched.
  p05e <- p05e + scale_colour_discrete(
    breaks = CT_ORDER,
    labels = paste0(CT_ORDER, " \u2014 ", CT_NAMES[CT_ORDER])
  )
}
if (identical(CELLTYPE_LABEL_MODE, "leader")) {
  # The helpers read UMAP_1/UMAP_2, so this embedding's columns are copied under
  # those names -- the helpers themselves stay byte-identical to 04c's.
  ct_lab_df <- data.frame(
    UMAP_1 = emb$UMAPHARMONY_1,
    UMAP_2 = emb$UMAPHARMONY_2,
    label  = unname(CT_LABEL_TEXT[as.character(obj$CellType_model)]),
    stringsAsFactors = FALSE
  )
  # The dot on each anchor is filled with that type's colour, so it needs the
  # palette keyed by the LABEL TEXT rather than by the bare code. Under
  # PALETTE = "seurat" there is no manual palette, so DimPlot's own hue ring is
  # reproduced at this panel's level count.
  lab_hex <- if (!is.null(CT_COLORS)) {
    unname(CT_COLORS[CT_ORDER])
  } else {
    scales::hue_pal()(length(CT_ORDER))
  }
  ct_lab_cols <- setNames(lab_hex, unname(CT_LABEL_TEXT[CT_ORDER]))

  p05e <- side_leader_labels(p05e, ct_lab_df, "label", ct_lab_cols,
                             label.size = 3, plot_w = CT_PANEL_W,
                             plot_h = CT_PANEL_H)
  message("Cell-type panel: ", length(unique(ct_lab_df$label)),
          " types named on leader lines outside the cloud")
}
save_plot(p05e, "05e_umap_ADSER_polyp_celltype", w = CT_PANEL_W, h = CT_PANEL_H)

if (DO_CLUSTER) {
  # 05d-equivalent: the cluster-ID review panel. Numbers stay ON the clusters
  # here (CLUSTER_PANEL_LABELS) because reading an id off the plot is the whole
  # point of it, and the colours are left to DimPlot -- the cell-type palette
  # must not be reused for cluster ids.
  save_plot(
    DimPlot(obj, reduction = "umap_harmony", group.by = "cluster",
            label = CLUSTER_PANEL_LABELS, raster = FALSE) + NoLegend(),
    paste0("05d_umap_ADSER_polyp_cluster_res", CLUSTER_RES)
  )
}

if (MAKE_EXTRA_PANELS) {
  # Cohort: shows that DIS, VAL and PCA83 are mixed after Harmony rather than
  # forming three separate islands.
  save_plot(
    DimPlot(obj, reduction = "umap_harmony", group.by = "cohort_model",
            label = FALSE, raster = FALSE) +
      scale_colour_manual(values = c(DIS = "#1B7837", VAL = "#762A83",
                                     PCA83 = "#8C8C8C"), name = NULL),
    "05e_umap_ADSER_polyp_cohort"
  )
  # Condition: AD vs SER. NA marks the three PCA83 add-on specimens, whose
  # AD/SER label is absent from the genotype-matched key table.
  save_plot(
    DimPlot(obj, reduction = "umap_harmony", group.by = "condition_model",
            label = FALSE, raster = FALSE) +
      scale_colour_manual(values = c(AD = "#B5382A", SER = "#17557F",
                                     unlabelled = "#D9D9D9"), name = NULL),
    "05e_umap_ADSER_polyp_condition"
  )
}

writeLines(capture.output(sessionInfo()), file.path(OUT_DIR, "sessionInfo_09.txt"))
message("\nDone. cells = ", ncol(obj), " | figures in ", FIG_DIR)
