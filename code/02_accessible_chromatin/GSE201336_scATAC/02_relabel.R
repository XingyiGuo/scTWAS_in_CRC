#!/usr/bin/env Rscript
# 02 -- Re-label Becker et al. scATAC epithelial cells into the Chen/Lau (Cell 2021)
#       nomenclature and attach sample metadata. No label transfer, no external reference.
#
# Normal / Unaffected / Adenoma:  ABS, CT, EE, GOB, STM, TAC
#   (Best4+ Enterocytes dropped; TUF and ASC are not available in this dataset -- see README)
#
#   Rscript 02_relabel.R --archr ./archr --labels epithelial_celltypes_atac.tsv
#
# Resources: minutes; light on memory (no LSI, no integration).

suppressPackageStartupMessages({
  library(ArchR)
  library(optparse)
})

opt <- parse_args(OptionParser(option_list = list(
  make_option("--archr",    type = "character", default = "./archr"),
  make_option("--labels",   type = "character", default = "epithelial_celltypes_atac.tsv",
              help = "consolidated ATAC cell-type TSV (Cell, Sample, Donor, CellType)"),
  make_option("--manifest", type = "character", default = "GSE201336_atac_manifest.csv"),
  make_option("--labelMap", type = "character", default = "lau_label_map.csv"),
  make_option("--threads",  type = "integer",   default = 16),
  make_option("--project",  type = "character", default = "proj_epithelial_lau")
)))

archr     <- normalizePath(opt$archr)
labelsCsv <- normalizePath(opt$labels)          # resolve before setwd()
manifest  <- read.csv(normalizePath(opt$manifest), stringsAsFactors = FALSE)
labelMap  <- read.csv(normalizePath(opt$labelMap), stringsAsFactors = FALSE)

setwd(archr)
addArchRThreads(opt$threads)
addArchRGenome("hg38")

## ------------------------------------------------------- 1) Becker labels -> project
proj <- loadArchRProject("proj_all_cells")

# Single consolidated annotation file: one row per epithelial cell, already de-duplicated,
# covering all conditions. Columns: Cell ("<Sample>#<barcode>"), Sample, Donor, CellType.
ann <- read.table(labelsCsv, sep = "\t", header = TRUE, stringsAsFactors = FALSE,
                  quote = "\"", check.names = FALSE)
if (!all(c("Cell", "CellType") %in% colnames(ann)))
  stop("Annotation file must have columns Cell and CellType. Found: ",
       paste(colnames(ann), collapse = ", "))

labels <- ann[!duplicated(ann$Cell), c("Cell", "CellType")]
labels <- labels[labels$Cell %in% proj$cellNames, ]
if (nrow(labels) == 0) stop("No annotation cell names matched. Run check_annotation_join.R.")
message(sprintf("Matched %d / %d annotated epithelial cells to the ArchR project.",
                nrow(labels), nrow(ann)))

proj <- proj[labels$Cell, ]
proj <- addCellColData(proj, data = labels$CellType, cells = labels$Cell,
                       name = "BeckerCellType", force = TRUE)

## ------------------------------------------------------- 2) map to Lau nomenclature
unmapped <- setdiff(unique(proj$BeckerCellType), labelMap$BeckerATAC)
if (length(unmapped)) stop("Unmapped Becker labels: ", paste(unmapped, collapse = ", "))

lau <- labelMap$LauLabel[match(proj$BeckerCellType, labelMap$BeckerATAC)]
proj <- addCellColData(proj, data = lau, cells = proj$cellNames, name = "CellTypeLau", force = TRUE)

# Best4+ Enterocytes are a distinct BEST4/OTOP2/CA7 lineage, not crypt-top colonocytes.
# Per project decision they are excluded rather than merged into CT.
nDrop <- sum(proj$CellTypeLau == "DROP")
message(sprintf("Dropping %d Best4+ enterocytes (not part of the ABS/CT schema).", nDrop))
proj <- proj[proj$CellTypeLau != "DROP", ]

# Not available in this dataset (both documented in README):
#   TUF -- never annotated in Becker's scATAC epithelial reference.
#   ASC -- adenoma cells were projected into the normal manifold, so they carry normal-like
#          labels; deriving ASC needs an external scRNA reference, which we are not using.
message("NOTE: TUF and ASC are not produced -- neither exists in this scATAC annotation.")

## ------------------------------------------------------- 3) sample metadata
m <- manifest[match(proj$Sample, manifest$ATAC_Sample), ]
for (col in c("Donor", "FAP", "DiseaseState", "Condition", "Location")) {
  proj <- addCellColData(proj, data = as.character(m[[col]]),
                         cells = proj$cellNames, name = col, force = TRUE)
}

## ------------------------------------------------------- 4) summary + save
cat("\nFinal cells per Condition x CellTypeLau:\n")
print(table(proj$Condition, proj$CellTypeLau))

write.csv(as.data.frame(table(proj$Condition, proj$CellTypeLau)),
          file.path(archr, "celltype_counts_lau.csv"), row.names = FALSE)

proj$PeakCallingGroup <- paste(proj$Condition, proj$CellTypeLau, sep = "__")
saveArchRProject(proj, outputDirectory = opt$project, load = FALSE)
message("\nSaved ", opt$project)
