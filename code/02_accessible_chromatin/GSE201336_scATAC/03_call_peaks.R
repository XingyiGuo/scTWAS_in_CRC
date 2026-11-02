#!/usr/bin/env Rscript
# 03 -- Reproducible peak set for every CELL TYPE, pooling cells across tissue
#       conditions (normal + adenoma/polyp).  ArchR pseudo-bulk replicates -> MACS2
#       -> iterative overlap merging.
#
# DEFAULT for this project is --byCondition: a separate peak set per Condition x CellType,
# i.e. Normal__{STM,TAC,ABS,CT,GOB,EE} and Adenoma__{STM,TAC,ABS,CT,GOB,EE}.
# Without --byCondition it instead builds one consensus set per cell type, pooling cells of
# that type across conditions. Condition is kept as per-cell metadata either way.
#
#   Rscript 03_call_peaks.R --archr ./archr --conditions "Normal,Adenoma" --minCells 200
#
# NOTE ON TERMINOLOGY: in GSE201336, disease state "Polyp" == conventional adenoma. There
# is no separate "adenoma" vs "polyp" axis in this dataset; both map to Condition "Adenoma".
#
# To instead keep peak sets split by condition (the old behaviour), pass --byCondition.
#
# Resources: ~6-10 h on 16 cores, ~64-96 GB RAM, ~150 GB scratch. Set TMPDIR to local SSD.

suppressPackageStartupMessages({
  library(ArchR)
  library(optparse)
  library(BSgenome.Hsapiens.UCSC.hg38)
})

opt <- parse_args(OptionParser(option_list = list(
  make_option("--archr",       type = "character", default = "./archr"),
  make_option("--project",     type = "character", default = "proj_epithelial_lau"),
  make_option("--out",         type = "character", default = "proj_epithelial_peaks"),
  make_option("--threads",     type = "integer",   default = 16),
  make_option("--minCells",    type = "integer",   default = 100),
  make_option("--conditions",  type = "character", default = "Normal,Adenoma",
              help = "conditions to pool for peak calling [default %default]. Options: Normal, Unaffected, Unaffected_dysplastic, Adenoma, CRC"),
  make_option("--byCondition", action = "store_true", default = FALSE,
              help = "call peaks per Condition x CellType instead of per CellType"),
  make_option("--macs2",       type = "character", default = NULL)
)))

setwd(normalizePath(opt$archr))
addArchRThreads(opt$threads)
addArchRGenome("hg38")
`%ni%` <- Negate(`%in%`)

proj <- loadArchRProject(opt$project)

## ------------------------------------------------- select the sample pool
keep <- trimws(strsplit(opt$conditions, ",")[[1]])
avail <- unique(proj$Condition)
missing <- setdiff(keep, avail)
if (length(missing))
  warning("Requested conditions not present in project: ", paste(missing, collapse = ", "),
          "\n(Did step 01 build arrows for them? Available: ", paste(avail, collapse = ", "), ")")
proj <- proj[proj$Condition %in% keep, ]
message(sprintf("Pooling %d cells across conditions: %s",
                nCells(proj), paste(intersect(keep, avail), collapse = ", ")))

## ------------------------------------------------- define peak-calling groups
# Default: one group per cell type, pooled across the selected conditions.
if (opt$byCondition) {
  proj$PeakCallingGroup <- paste(proj$Condition, proj$CellTypeLau, sep = "__")
  groupLabel <- "Condition x CellType"
} else {
  proj$PeakCallingGroup <- proj$CellTypeLau
  groupLabel <- "CellType (pooled across conditions)"
}

counts <- sort(table(proj$PeakCallingGroup))
cat(sprintf("\nCells per %s group:\n", groupLabel)); print(counts)

# Cross-tab so you can see how each cell-type peak set is composed by condition
cat("\nComposition (cells) by Condition x CellType:\n")
print(table(proj$Condition, proj$CellTypeLau))

## ------------------------------------------------- drop under-powered groups
# MACS2 on a handful of cells produces noise. Becker et al. required >=300 cells and
# hand-merged rare groups; here we drop anything under --minCells and record what was lost.
# Pooling across conditions means most groups clear the bar comfortably; EE may still be thin.
tooSmall <- names(counts)[counts < opt$minCells]
write.csv(data.frame(Group = names(counts), nCells = as.integer(counts),
                     Dropped = names(counts) %in% tooSmall),
          "peak_calling_group_sizes.csv", row.names = FALSE)

if (length(tooSmall)) {
  warning(sprintf("Dropping %d group(s) with < %d cells: %s",
                  length(tooSmall), opt$minCells, paste(tooSmall, collapse = ", ")))
  proj <- proj[proj$PeakCallingGroup %ni% tooSmall, ]
}

## ------------------------------------------------- reproducible peaks
proj <- addGroupCoverages(ArchRProj = proj, groupBy = "PeakCallingGroup", force = TRUE)

pathToMacs2 <- if (is.null(opt$macs2)) findMacs2() else opt$macs2
message("MACS2: ", pathToMacs2)

# 500 bp fixed-width peaks, 150k cap per group, >=2 reproducible pseudo-bulk replicates
proj <- addReproduciblePeakSet(
  ArchRProj   = proj,
  groupBy     = "PeakCallingGroup",
  pathToMacs2 = pathToMacs2,
  genomeSize  = 2.7e9,
  force       = TRUE
)

proj <- addPeakMatrix(proj, force = TRUE)
saveArchRProject(proj, outputDirectory = opt$out, load = FALSE)

message(sprintf("\nUnion peak set: %d peaks across %d %s groups.",
                length(getPeakSet(proj)), length(unique(proj$PeakCallingGroup)), groupLabel))
