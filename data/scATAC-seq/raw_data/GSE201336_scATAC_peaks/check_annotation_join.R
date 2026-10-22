#!/usr/bin/env Rscript
# PRE-FLIGHT CHECK -- run this before 01/02. Costs seconds; saves hours.
#
# Confirms that the sample prefixes in the consolidated ATAC cell-type annotation file
# match the ATAC_Sample column of the manifest. If they don't, ArchR will happily build
# arrows and then match zero cells.
#
# Expects the single file you downloaded: epithelial_celltypes_atac.tsv
#   columns: Cell ("<ATAC_Sample>#<barcode>"), Sample, Donor, CellType
#
#   Rscript check_annotation_join.R --labels epithelial_celltypes_atac.tsv

suppressPackageStartupMessages(library(optparse))
opt <- parse_args(OptionParser(option_list = list(
  make_option("--labels",   type = "character", default = "epithelial_celltypes_atac.tsv",
              help = "consolidated ATAC cell-type TSV [default %default]"),
  make_option("--manifest", type = "character", default = "GSE201336_atac_manifest.csv")
)))

manifest <- read.csv(opt$manifest, stringsAsFactors = FALSE)

ann <- read.table(opt$labels, sep = "\t", header = TRUE, stringsAsFactors = FALSE,
                  quote = "\"", check.names = FALSE)
need <- c("Cell", "CellType")
if (!all(need %in% colnames(ann)))
  stop("Annotation file must have columns Cell and CellType. Found: ",
       paste(colnames(ann), collapse = ", "))

# Prefer the file's own Sample column; fall back to parsing it out of Cell.
ann$Sample <- if ("Sample" %in% colnames(ann)) ann$Sample else sub("#.*", "", ann$Cell)
samples <- sort(unique(ann$Sample))

inManifest <- intersect(samples, manifest$ATAC_Sample)
notInMan   <- setdiff(samples, manifest$ATAC_Sample)
notInAnnot <- setdiff(manifest$ATAC_Sample, samples)

cat(sprintf("Annotation file: %d cells, %d distinct samples\n", nrow(ann), length(samples)))
cat(sprintf("  matched to manifest    : %d\n", length(inManifest)))
cat(sprintf("  in annot, NOT manifest : %d %s\n", length(notInMan),
            if (length(notInMan)) paste0("-> ", paste(notInMan, collapse = ", ")) else ""))
cat(sprintf("  in manifest, NOT annot : %d %s\n", length(notInAnnot),
            if (length(notInAnnot)) paste0("-> ", paste(notInAnnot, collapse = ", ")) else ""))

cat("\nEpithelial cell types present:\n")
print(sort(table(ann$CellType), decreasing = TRUE))

cat("\nCells per Condition x CellType (projected from manifest):\n")
ann$Condition <- manifest$Condition[match(ann$Sample, manifest$ATAC_Sample)]
print(table(ann$Condition, ann$CellType))

if (length(notInMan) > 0) {
  stop("Sample-name mismatch. Fix ARROW_NAME_OVERRIDES in build_manifest.py and re-run.")
}
cat("\nOK -- names line up. Safe to run 01_build_arrows.R\n")
