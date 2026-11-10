#!/usr/bin/env Rscript
# 04 -- Export the reproducible peaks to BED, the union peak set, and a pseudobulk
#       accessibility matrix (broken out by Condition x CellType). Differentials are opt-in.
#
#   Rscript 04_export_peaks.R --archr ./archr --out ./peaks
#
# With the project's default (step 03 --byCondition), per_group/ contains one BED per
# Condition x CellType: Normal__{STM,TAC,ABS,CT,GOB,EE} and Adenoma__{...same six}.

suppressPackageStartupMessages({
  library(ArchR)
  library(optparse)
  library(rtracklayer)
  library(SummarizedExperiment)
  library(GenomicRanges)
  library(Matrix)
})

opt <- parse_args(OptionParser(option_list = list(
  make_option("--archr",   type = "character", default = "./archr"),
  make_option("--project", type = "character", default = "proj_epithelial_peaks"),
  make_option("--out",     type = "character", default = "./peaks"),
  make_option("--threads", type = "integer",   default = 16),
  make_option("--diff",    action = "store_true", default = FALSE,
              help = "also compute differential-peak comparisons (off by default)"),
  make_option("--bigwigs", action = "store_true", default = FALSE)
)))

dir.create(opt$out, recursive = TRUE, showWarnings = FALSE)
outdir <- normalizePath(opt$out, mustWork = TRUE)

setwd(normalizePath(opt$archr, mustWork = TRUE))
addArchRThreads(opt$threads)
addArchRGenome("hg38")

dir.create(file.path(outdir, "per_group"), recursive = TRUE, showWarnings = FALSE)

proj <- loadArchRProject(opt$project)

## ------------------------------------------- per-group reproducible peaks -> BED
peakFiles <- list.files(
  file.path(opt$project, "PeakCalls"),
  pattern = "-reproduciblePeaks\\.gr\\.rds$",
  full.names = TRUE,
  recursive = TRUE
)

if (length(peakFiles) == 0) {
  stop("No reproducible peak files found under: ",
       file.path(opt$project, "PeakCalls"))
}
message("Found ", length(peakFiles), " reproducible peak files.")

summary <- do.call(rbind, lapply(peakFiles, function(f) {
  group <- sub("-reproduciblePeaks\\.gr\\.rds$", "", basename(f))
  gr    <- readRDS(f)

  bed <- data.frame(
    chrom       = as.character(seqnames(gr)),
    chromStart  = start(gr) - 1L,          # BED is 0-based half-open
    chromEnd    = end(gr),
    name        = paste0(group, "_peak_", seq_along(gr)),
    score       = round(pmin(mcols(gr)$score, 1000)),
    strand      = ".",
    signalValue = mcols(gr)$Reproducibility,
    negLog10p   = mcols(gr)$score,
    peakType    = mcols(gr)$peakType
  )
  bed <- bed[order(bed$chrom, bed$chromStart), ]
  write.table(bed, file.path(outdir, "per_group", paste0(group, ".narrowPeak.bed")),
              sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

  data.frame(Group = group, nPeaks = length(gr), medianWidth = median(width(gr)),
             pctDistal   = round(100 * mean(mcols(gr)$peakType == "Distal"), 1),
             pctPromoter = round(100 * mean(mcols(gr)$peakType == "Promoter"), 1))
}))

write.csv(summary, file.path(outdir, "peak_summary_by_group.csv"), row.names = FALSE)
print(summary)

## ------------------------------------------- union peak set
union <- getPeakSet(proj)
export.bed(union, file.path(outdir, "union_peakset_epithelial.bed"))
saveRDS(union, file.path(outdir, "union_peakset_epithelial.rds"))

## ------------------------------------------- pseudobulk CPM matrix
# Peaks were called per cell type, but we aggregate accessibility by Condition x CellType
# so each cell type can be compared directly between normal and adenoma on shared features.
se  <- getMatrixFromProject(proj, useMatrix = "PeakMatrix")
ord <- match(colnames(se), proj$cellNames)
grp <- factor(paste(proj$Condition[ord], proj$CellTypeLau[ord], sep = "__"))
agg <- as.matrix(assay(se) %*% Matrix::sparse.model.matrix(~ 0 + grp))
colnames(agg) <- levels(grp)
rownames(agg) <- paste0(seqnames(rowRanges(se)), ":", start(rowRanges(se)), "-", end(rowRanges(se)))
cpm <- t(t(agg) / colSums(agg)) * 1e6
write.csv(round(cpm, 3), gzfile(file.path(outdir, "peak_by_condition_celltype_cpm.csv.gz")))

## ------------------------------------------- differential peaks (opt-in)
# Off by default: the task is to call peaks per cell type per condition, not compare them.
# Pass --diff to also emit Adenoma-vs-Normal (per cell type) and secretory-vs-absorptive
# marker-peak tables.
if (opt$diff) {
  diff <- function(sub, groupBy, fg, bg, fname, cut = "FDR <= 0.05 & abs(Log2FC) >= 1") {
    if (!all(c(fg, bg) %in% getCellColData(sub, groupBy)[[groupBy]])) {
      message("skip ", fname, " (missing group)"); return(invisible(NULL))
    }
    mk <- getMarkerFeatures(sub, useMatrix = "PeakMatrix", groupBy = groupBy,
                            useGroups = fg, bgdGroups = bg,
                            bias = c("TSSEnrichment", "log10(nFrags)"), testMethod = "wilcoxon")
    write.csv(as.data.frame(getMarkers(mk, cutOff = cut)[[1]]),
              file.path(outdir, fname), row.names = FALSE)
  }

  for (cond in c("Normal", "Adenoma")) {
    sub <- proj[proj$Condition == cond, ]
    diff(sub, "CellTypeLau", "GOB", "CT",  sprintf("diff_%s_GOB_vs_CT.csv", cond))
    diff(sub, "CellTypeLau", "GOB", "ABS", sprintf("diff_%s_GOB_vs_ABS.csv", cond))
  }

  shared <- intersect(unique(proj$CellTypeLau[proj$Condition == "Adenoma"]),
                      unique(proj$CellTypeLau[proj$Condition == "Normal"]))
  for (ct in shared) {
    sub <- proj[proj$CellTypeLau == ct & proj$Condition %in% c("Adenoma", "Normal"), ]
    diff(sub, "Condition", "Adenoma", "Normal", sprintf("diff_Adenoma_vs_Normal_%s.csv", ct))
  }
}

## ------------------------------------------- optional bigWigs
if (opt$bigwigs) {
  getGroupBW(proj, groupBy = "PeakCallingGroup", normMethod = "ReadsInTSS",
             tileSize = 100, maxCells = 1000)
}

message("Wrote peaks to ", outdir)
