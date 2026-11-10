#!/usr/bin/env Rscript
# Download datasets from:  https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE201336
#
# 01 -- Build ArchR Arrow files from GSE201336 fragment files.
#
# Arrow names are set to the authors' ATAC_Sample IDs (e.g. "A001-C-002-D") so that
# ArchR cell names ("<ATAC_Sample>#<barcode>") match the cell-type annotation files
# distributed with Becker et al. 2022.
#
#   Rscript 01_build_arrows.R --fragments ./fragments --out ./archr --threads 16
#
# Resources: ~1-2 h on 16 cores for ~80 samples; peak RSS ~6 GB/thread. Arrows total ~40 GB.

suppressPackageStartupMessages({
  library(ArchR)
  library(optparse)
})

opt <- parse_args(OptionParser(option_list = list(
  make_option("--fragments", type = "character"),
  make_option("--manifest",  type = "character", default = "GSE201336_atac_manifest.csv"),
  make_option("--out",       type = "character", default = "./archr"),
  make_option("--threads",   type = "integer",   default = 16),
  make_option("--conditions",type = "character", default = "Normal,Unaffected,Adenoma",
              help = "comma-separated subset of Condition column [default %default]"),
  make_option("--minTSS",    type = "double",    default = 6),
  make_option("--minFrags",  type = "integer",   default = 1000)
)))

stopifnot(!is.null(opt$fragments))

# Resolve input paths before changing the working directory. This keeps relative
# --manifest/--fragments paths relative to the directory from which the script was launched.
manifestPath <- normalizePath(opt$manifest, mustWork = TRUE)
fragmentsDir <- normalizePath(opt$fragments, mustWork = TRUE)
dir.create(opt$out, showWarnings = FALSE, recursive = TRUE)
outDir <- normalizePath(opt$out, mustWork = TRUE)

manifest <- read.csv(manifestPath, stringsAsFactors = FALSE)
requiredCols <- c("GSM", "ATAC_Sample", "Condition")
missingCols <- setdiff(requiredCols, colnames(manifest))
if (length(missingCols) > 0L) {
  stop(sprintf("manifest is missing required column(s): %s",
               paste(missingCols, collapse = ", ")))
}

keep     <- trimws(strsplit(opt$conditions, ",")[[1]])
manifest <- manifest[manifest$Condition %in% keep, , drop = FALSE]
if (nrow(manifest) == 0L) {
  stop(sprintf("no manifest rows matched --conditions: %s",
               paste(keep, collapse = ", ")))
}
if (anyNA(manifest$GSM) || any(!nzchar(manifest$GSM))) {
  stop("manifest contains missing/empty GSM values")
}
if (anyNA(manifest$ATAC_Sample) || any(!nzchar(manifest$ATAC_Sample))) {
  stop("manifest contains missing/empty ATAC_Sample values")
}
if (anyDuplicated(manifest$GSM)) {
  dup <- unique(manifest$GSM[duplicated(manifest$GSM)])
  stop(sprintf("GSM values must be unique after filtering; duplicated: %s",
               paste(dup, collapse = ", ")))
}
if (anyDuplicated(manifest$ATAC_Sample)) {
  dup <- unique(manifest$ATAC_Sample[duplicated(manifest$ATAC_Sample)])
  stop(sprintf("ATAC_Sample values must be unique for ArchR sampleNames; duplicated: %s",
               paste(dup, collapse = ", ")))
}

message(sprintf("Selected %d samples across conditions: %s",
                nrow(manifest), paste(keep, collapse = ", ")))

# GEO fragment files are stored flat under --fragments, e.g.
# GSM6058770_A001-C-002-D_20200702_fragments.tsv.gz.
# Match by GSM accession because ATAC_Sample is the normalized Arrow/annotation name
# and can intentionally differ from the sample token embedded in a GEO filename.
# The anchored suffix excludes tabix files such as *.fragments.tsv.gz.tbi.gz.
fragmentFiles <- list.files(
  fragmentsDir,
  pattern = "_fragments\\.tsv\\.gz$",
  full.names = TRUE,
  recursive = FALSE
)
if (length(fragmentFiles) == 0L) {
  stop(sprintf("no *_fragments.tsv.gz files found in %s", fragmentsDir))
}

inputFiles <- vapply(seq_len(nrow(manifest)), function(i) {
  gsm <- manifest$GSM[i]
  sample <- manifest$ATAC_Sample[i]
  hits <- fragmentFiles[startsWith(basename(fragmentFiles), paste0(gsm, "_"))]

  if (length(hits) != 1L) {
    stop(sprintf(
      "expected 1 fragments file for GSM %s (ATAC_Sample %s), found %d%s",
      gsm, sample, length(hits),
      if (length(hits)) paste0(": ", paste(basename(hits), collapse = ", ")) else ""
    ))
  }
  hits
}, character(1))
names(inputFiles) <- manifest$ATAC_Sample

message("Resolved fragment files:")
for (i in seq_along(inputFiles)) {
  message(sprintf("  %s <- %s", names(inputFiles)[i], basename(inputFiles[i])))
}

setwd(outDir)
addArchRThreads(threads = opt$threads)
addArchRGenome("hg38")   # Becker et al. aligned with Cell Ranger ATAC v1.2.0 against hg38

# QC thresholds: the paper filtered on TSS enrichment and unique fragments.
ArrowFiles <- createArrowFiles(
  inputFiles      = inputFiles,
  sampleNames     = names(inputFiles),
  minTSS          = opt$minTSS,
  minFrags        = opt$minFrags,
  addTileMat      = TRUE,
  addGeneScoreMat = TRUE,
  force           = FALSE
)

# Doublet scores. The authors additionally supply an explicit post-filter cell list;
# 02_relabel.R subsets to annotated cells, so the doublet scores are informational.
doubScores <- addDoubletScores(input = ArrowFiles, k = 10, knnMethod = "UMAP", LSIMethod = 1)

proj <- ArchRProject(ArrowFiles = ArrowFiles, outputDirectory = "proj_all_cells", copyArrows = FALSE)
saveArchRProject(proj, outputDirectory = "proj_all_cells", load = FALSE)

message(sprintf("Done: %d arrows, %d cells passing QC.", length(ArrowFiles), nCells(proj)))
