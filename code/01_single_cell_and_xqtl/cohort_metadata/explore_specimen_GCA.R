#!/usr/bin/env Rscript

# explore_specimen_GCA.R
# ------------------------------------------------------------------------------
# Interrogate one normal-colon specimen in the GSE266546 + HTAN integrated
# object: how many cells it contributes per cell type, how deep those cells are,
# and - the part that usually explains a surprising count - what the SAME
# PATIENT's other specimens look like.
#
#   object   integrated_seurat_final2_reviewed.res0.8.rm23_sub.rds   (~5.7 GB)
#   keys     meta.data$SampleID x meta.data$celltype_annot_reviewed
#            (the two columns 07_sample_level_pseudobulk_counts.R aggregates on)
#
# Usage:
#   Rscript explore_specimen_GCA.R                        # default specimen
#   Rscript explore_specimen_GCA.R HTA11_4255_2000002011
#
# Memory: readRDS on the 5.7 GB object needs roughly 20-40 GB - this is an
# ACCRE job, not a laptop one. Set READ_METADATA_ONLY <- TRUE to answer the
# same questions from the exported cell-level metadata CSV (50 MB, identical
# SampleID x celltype_annot_reviewed counts, no Seurat dependency). Only the
# per-cell QC summaries below need the object itself.
# ------------------------------------------------------------------------------

SPECIMEN <- {
  a <- commandArgs(trailingOnly = TRUE)
  if (length(a)) a[1] else "HTA11_4255_2000002011"
}

READ_METADATA_ONLY <- TRUE

# ---- paths -------------------------------------------------------------------

GCA_ROOT <- Filter(dir.exists, c(
  "/data/l2_bioinfo1/liq17/scTWAS/data/VUMC_HTAN_GSE266546/raw/rds2",
  file.path("L:", "scRNA-seq", "VUMC_HTAN_GSE266546", "raw", "rds2")
))[1]
if (is.na(GCA_ROOT)) stop("GSE266546 rds2 directory not found")

OBJ_RDS  <- file.path(GCA_ROOT, "integrated_seurat_final2_reviewed.res0.8.rm23_sub.rds")
OBJ_META <- file.path(GCA_ROOT, "integrated_seurat_final2_reviewed.res0.8.rm23_sub_metadata.csv")

CT7 <- c("ABS", "GOB", "CT", "EE", "STM", "TUF", "TAC")

# Columns worth carrying. QC columns exist in the metadata CSV too, but are NA
# for the HTAN batches, so per-cell depth is only reliable from the object.
KEEP <- c("SampleID", "PatientID", "celltype_annot_reviewed", "batch", "study",
          "sample_condition", "Tissue", "cluster", "nCount_RNA", "nFeature_RNA")

# ---- load --------------------------------------------------------------------

if (READ_METADATA_ONLY) {

  stopifnot(file.exists(OBJ_META))
  hdr <- names(read.csv(OBJ_META, nrows = 1, check.names = FALSE))
  md  <- read.csv(OBJ_META, colClasses = "character",
                  check.names = FALSE)[, intersect(KEEP, hdr), drop = FALSE]
  src <- paste("metadata CSV:", basename(OBJ_META))

} else {

  stopifnot(file.exists(OBJ_RDS))
  if (!requireNamespace("SeuratObject", quietly = TRUE)) {
    stop("SeuratObject is required to read ", basename(OBJ_RDS),
         "; or set READ_METADATA_ONLY <- TRUE.")
  }
  suppressPackageStartupMessages(library(SeuratObject))

  message("Reading ", OBJ_RDS, " (needs roughly 20-40 GB) ...")
  obj <- readRDS(OBJ_RDS)
  md  <- obj@meta.data[, intersect(KEEP, colnames(obj@meta.data)), drop = FALSE]
  src <- paste("Seurat object:", basename(OBJ_RDS))
}

for (cc in c("SampleID", "PatientID", "celltype_annot_reviewed"))
  md[[cc]] <- as.character(md[[cc]])
for (cc in intersect(c("nCount_RNA", "nFeature_RNA"), names(md)))
  md[[cc]] <- suppressWarnings(as.numeric(md[[cc]]))

hit <- md[md$SampleID == SPECIMEN, , drop = FALSE]
if (!nrow(hit)) stop("SampleID '", SPECIMEN, "' not found in ", src)

patient <- unique(hit$PatientID)[1]

# ==============================================================================
# 1. The specimen itself
# ==============================================================================

cat("\n=== ", SPECIMEN, "  (patient ", patient, ") ===\n", sep = "")
cat("source: ", src, "\n", sep = "")
cat("cells in object: ", nrow(hit), "\n", sep = "")
cat("batch / study / condition: ",
    paste(unique(hit$batch), collapse = ","), " / ",
    paste(unique(hit$study), collapse = ","), " / ",
    paste(unique(hit$sample_condition), collapse = ","), "\n\n", sep = "")

per_ct <- table(factor(hit$celltype_annot_reviewed, levels = CT7))
print(data.frame(CellType = names(per_ct), cells = as.integer(per_ct)),
      row.names = FALSE)

# Which clusters the cells sit in - a 1-cell cell type that sits in the right
# cluster is a real but rare cell; one sitting alone in an odd cluster is worth
# a second look.
if ("cluster" %in% names(md)) {
  cat("\ncells by cell type x cluster:\n")
  print(table(hit$celltype_annot_reviewed, hit$cluster))
}

# Per-cell depth: only populated when the object (or a batch that carries QC
# columns) is used.
if (all(c("nCount_RNA", "nFeature_RNA") %in% names(md)) &&
    any(!is.na(hit$nCount_RNA))) {
  cat("\nper-cell depth by cell type (median nCount / nFeature):\n")
  d <- aggregate(cbind(nCount_RNA, nFeature_RNA) ~ celltype_annot_reviewed,
                 data = hit, FUN = function(x) round(stats::median(x)))
  print(d, row.names = FALSE)
} else {
  cat("\nper-cell depth: not available from this source",
      " (QC columns are NA for the HTAN batches in the CSV)\n", sep = "")
}

# ==============================================================================
# 2. The same patient's other specimens - the usual explanation
# ==============================================================================

sibs <- md[md$PatientID == patient, , drop = FALSE]

cat("\n--- all specimens of patient ", patient, " ---\n", sep = "")
tb <- table(sibs$SampleID, factor(sibs$celltype_annot_reviewed, levels = CT7))
tb <- cbind(as.data.frame.matrix(tb), total = rowSums(tb))
tb$in_scTWAS_panel <- rownames(tb) == SPECIMEN
print(tb)

if (nrow(tb) > 1) {
  other <- tb[rownames(tb) != SPECIMEN, , drop = FALSE]
  deeper <- rownames(other)[other$total > tb[SPECIMEN, "total"]]
  if (length(deeper)) {
    cat("\nNOTE: ", paste(deeper, collapse = ", "),
        " has more cells than the specimen used in the panel.\n",
        "The panel takes one specimen per donor by a fixed rule (the ",
        "_2000002011 segment is preferred over _2000002021), so a deeper ",
        "sibling is skipped by design, not by error.\n", sep = "")
  }
}
