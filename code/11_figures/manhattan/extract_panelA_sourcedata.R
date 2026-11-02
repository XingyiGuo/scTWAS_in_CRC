#!/usr/bin/env Rscript
###############################################################################
# extract_panelA_sourcedata.R
#
# WHAT IT DOES
#   Reads a Seurat .rds and writes the small cell-level table that
#   Final_Viz_4panel.R needs in order to DRAW panel A itself:
#       cell_id, UMAPHARMONY_1, UMAPHARMONY_2, annotation
#   one row per cell, plus any extra metadata columns asked for with --extra.
#   Nothing is plotted here; this script only exports coordinates.
#
# WHY IT HAS TO RUN ON ACCRE
#   Those coordinates exist nowhere else.  Seurat keeps reduction embeddings
#   outside @meta.data, so the metadata CSV that 05_annotation.R writes next to
#   the object has celltype_annot_reviewed but no UMAP columns, and the object
#   itself is 5.76 GB on disk and needs tens of GB of RAM to load.  The figure
#   script then reads only the CSV, so this is a one-off.
#
# HOW TO RUN ON ACCRE
#   sbatch extract_panelA_sourcedata.sh          # wrapper in this folder
#   # or interactively, inside the same container the rest of the pipeline uses:
#   Rscript extract_panelA_sourcedata.R
#   # any default can be overridden, e.g. for a different object/reduction:
#   Rscript extract_panelA_sourcedata.R \
#     --rds=/data/.../integrated_seurat_final2_reviewed.res0.8.rm23_sub.rds \
#     --out=/data/.../umap_HTAN_GSE266546_reviewed_celltype_sourcedata.csv \
#     --reduction=umap_harmony --group=celltype_annot_reviewed \
#     --extra=cluster,sample,study
#   Positional use still works:  Rscript extract_panelA_sourcedata.R <rds> <out>
#
# AFTERWARDS
#   Copy the CSV (a few MB) to the Windows side, into
#     L:/scRNA-seq/VUMC_HTAN_GSE266546/raw/rds2/
#   which is the first path Final_Viz_4panel.R looks in, and re-run the figure.
#
# EXIT CODES
#   0 = table written; non-zero = nothing written, and the reason is the last
#   line of the job log.  The script never writes a partial file: it validates
#   the object first and only then opens the output.
###############################################################################

## --- defaults (ACCRE paths) --------------------------------------------------
PROJECT <- "/data/l2_bioinfo1/liq17/scTWAS/data/VUMC_HTAN_GSE266546/src"
RDS_DIR <- "/data/l2_bioinfo1/liq17/scTWAS/data/VUMC_HTAN_GSE266546/raw/rds2"

DEFAULT_RDS <- file.path(RDS_DIR, "integrated_seurat_final2_reviewed.res0.8.rm23_sub.rds")
DEFAULT_OUT <- file.path(RDS_DIR, "umap_HTAN_GSE266546_reviewed_celltype_sourcedata.csv")
DEFAULT_REDUCTION <- "umap_harmony"            # as used by 05_annotation.R
DEFAULT_GROUP     <- "celltype_annot_reviewed" # the column step 05e colours by
DEFAULT_EXTRA     <- ""                        # e.g. "cluster,sample,study"

## --- renv, exactly as 05_annotation.R does it --------------------------------
# Harmless when the project is not renv-managed or when the script runs inside
# the RStudio container with its own library.
activate <- file.path(PROJECT, "renv", "activate.R")
if (file.exists(activate)) {
  message("Activating renv: ", activate)
  source(activate)
}

suppressPackageStartupMessages({
  has_dt <- requireNamespace("data.table", quietly = TRUE)
})

## --- arguments: --key=value, or positional <rds> <out> -----------------------
args <- commandArgs(trailingOnly = TRUE)
flags <- grepl("^--", args)
opts <- list()
if (any(flags)) {
  kv <- strsplit(sub("^--", "", args[flags]), "=", fixed = TRUE)
  for (p in kv) opts[[p[1]]] <- if (length(p) > 1) paste(p[-1], collapse = "=") else "TRUE"
}
pos <- args[!flags]
opt <- function(name, default, i = NA) {
  if (!is.null(opts[[name]])) return(opts[[name]])
  if (!is.na(i) && length(pos) >= i) return(pos[i])
  default
}

RDS       <- opt("rds", DEFAULT_RDS, 1)
OUT       <- opt("out", DEFAULT_OUT, 2)
REDUCTION <- opt("reduction", DEFAULT_REDUCTION)
GROUP     <- opt("group", DEFAULT_GROUP)
EXTRA     <- setdiff(trimws(strsplit(opt("extra", DEFAULT_EXTRA), ",")[[1]]), "")

message("\n================ extract_panelA_sourcedata ================")
message("host      : ", Sys.info()[["nodename"]])
message("R         : ", R.version.string)
message("rds       : ", RDS)
message("out       : ", OUT)
message("reduction : ", REDUCTION)
message("group     : ", GROUP)
if (length(EXTRA)) message("extra     : ", paste(EXTRA, collapse = ", "))
message("start     : ", format(Sys.time()))

if (!file.exists(RDS)) stop("no such .rds: ", RDS)
message("size      : ", format(round(file.size(RDS) / 1e9, 2)), " GB on disk")

## --- read the object ---------------------------------------------------------
# This is the slow, memory-hungry step: allow ~20x the file size in RAM and
# expect several minutes on a shared filesystem.
t0 <- Sys.time()
obj <- readRDS(RDS)
message("loaded    : ", format(round(difftime(Sys.time(), t0, units = "mins"), 1)),
        " min, class ", paste(class(obj), collapse = "/"))

## --- embeddings: three routes, so this works with or without Seurat attached -
emb <- NULL
route <- NA_character_
if (requireNamespace("SeuratObject", quietly = TRUE)) {
  emb <- try(SeuratObject::Embeddings(obj, reduction = REDUCTION), silent = TRUE)
  if (inherits(emb, "try-error")) emb <- NULL else route <- "SeuratObject::Embeddings"
}
if (is.null(emb)) {
  emb <- try(obj@reductions[[REDUCTION]]@cell.embeddings, silent = TRUE)
  if (inherits(emb, "try-error")) emb <- NULL else route <- "@reductions slot"
}
if (is.null(emb)) {
  # last resort: reach the slots through attributes, which does not need the
  # Seurat S4 class definitions to be loaded at all
  red <- try(attr(obj, "reductions")[[REDUCTION]], silent = TRUE)
  if (!inherits(red, "try-error") && !is.null(red)) {
    emb <- attr(red, "cell.embeddings")
    route <- "attributes"
  }
}
if (is.null(emb)) {
  available <- try(names(attr(obj, "reductions")), silent = TRUE)
  stop("reduction '", REDUCTION, "' could not be read. Reductions present: ",
       if (inherits(available, "try-error")) "unknown" else paste(available, collapse = ", "))
}
emb <- as.matrix(emb)
if (ncol(emb) < 2) stop("reduction '", REDUCTION, "' has fewer than 2 dimensions")
message("embedding : ", nrow(emb), " cells x ", ncol(emb), " dims, via ", route)

## --- metadata, aligned to the embedding BY BARCODE ---------------------------
md <- try(obj@meta.data, silent = TRUE)
if (inherits(md, "try-error") || is.null(md)) md <- attr(obj, "meta.data")
if (is.null(md)) stop("could not reach @meta.data")

cells <- rownames(emb)
if (!is.null(cells) && !is.null(rownames(md)) && all(cells %in% rownames(md))) {
  md <- md[cells, , drop = FALSE]      # never assume the two share a row order
  message("alignment : metadata matched to the embedding by barcode")
} else {
  if (nrow(md) != nrow(emb))
    stop("metadata has ", nrow(md), " rows but the embedding has ", nrow(emb),
         ", and barcodes do not match - refusing to guess the pairing")
  message("alignment : no usable barcodes; falling back to row order")
}

want <- c(GROUP, EXTRA)
missing <- setdiff(want, names(md))
if (length(missing))
  stop("column(s) not in @meta.data: ", paste(missing, collapse = ", "),
       "\nAvailable: ", paste(names(md), collapse = ", "))

## --- assemble ----------------------------------------------------------------
out <- data.frame(
  cell_id       = if (is.null(cells)) rownames(md) else cells,
  UMAPHARMONY_1 = as.numeric(emb[, 1]),
  UMAPHARMONY_2 = as.numeric(emb[, 2]),
  annotation    = as.character(md[[GROUP]]),
  stringsAsFactors = FALSE
)
for (cl in EXTRA) out[[cl]] <- as.character(md[[cl]])

# Free the object before writing, so the write has headroom on a tight job.
rm(obj, emb, md); invisible(gc())

## --- checks worth failing on -------------------------------------------------
n_bad_xy <- sum(!is.finite(out$UMAPHARMONY_1) | !is.finite(out$UMAPHARMONY_2))
n_bad_an <- sum(is.na(out$annotation) | out$annotation == "")
if (n_bad_xy) message("WARNING   : ", n_bad_xy, " cell(s) have a non-finite coordinate")
if (n_bad_an) message("WARNING   : ", n_bad_an, " cell(s) have no annotation")
if (nrow(out) == 0) stop("nothing to write")

# The figure keeps only the seven crypt codes, so flag anything else HERE,
# where the object is still in front of you, rather than in the figure script.
CRYPT <- c("STM", "TAC", "ABS", "CT", "GOB", "EE", "TUF")
extra_lv <- setdiff(unique(out$annotation), c(CRYPT, NA_character_, ""))
if (length(extra_lv))
  message("NOTE      : annotation values outside the seven crypt codes, which ",
          "Final_Viz_4panel.R will drop: ", paste(extra_lv, collapse = ", "))

## --- write -------------------------------------------------------------------
dir.create(dirname(OUT), recursive = TRUE, showWarnings = FALSE)
if (file.exists(OUT)) message("NOTE      : overwriting ", OUT)
if (has_dt) data.table::fwrite(out, OUT) else write.csv(out, OUT, row.names = FALSE)
if (!file.exists(OUT)) stop("write failed: ", OUT)

message("\nwrote     : ", OUT, "  (", format(round(file.size(OUT) / 1e6, 1)), " MB, ",
        format(nrow(out), big.mark = ","), " cells)")
message("cells per annotation:")
print(sort(table(out$annotation), decreasing = TRUE))
message("coordinate ranges:")
print(signif(sapply(out[, c("UMAPHARMONY_1", "UMAPHARMONY_2")], range), 4))
message("finished  : ", format(Sys.time()))
message("\n--- sessionInfo ---")
print(sessionInfo())
