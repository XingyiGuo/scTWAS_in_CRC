#!/usr/bin/env Rscript
# 02 -- Re-label Becker et al. scATAC epithelial cells into the Chen/Lau (Cell 2021)
#       nomenclature, and recover the adenoma-specific (ASC) population by label
#       transfer from the Lau-lab polyp scRNA atlas.
#
# Normal / Unaffected:  ABS, CT, EE, GOB, STM, TAC        (Best4+ dropped; TUF unavailable)
# Adenoma:              same six + ASC                    (+ SSC if serrated cells present)
#
#   Rscript 02_relabel_and_transfer.R --archr ./archr --annot ./cell_type_annotations \
#       --refRNA /path/to/chen_lau_polyp_epithelial.rds
#
# The reference must be a Seurat object with a metadata column (--refCellTypeCol,
# default "CellType") containing at least "ASC". Obtain from the HTAN portal
# (https://data.humantumoratlas.org/, Vanderbilt pre-cancer atlas) -- see README.
#
# Resources: ~2-3 h on 16 cores, ~64 GB RAM.

suppressPackageStartupMessages({
  library(ArchR)
  library(Seurat)
  library(optparse)
})

opt <- parse_args(OptionParser(option_list = list(
  make_option("--archr",     type = "character", default = "./archr"),
  make_option("--annot",     type = "character", default = "./cell_type_annotations"),
  make_option("--manifest",  type = "character", default = "GSE201336_atac_manifest.csv"),
  make_option("--labelMap",  type = "character", default = "lau_label_map.csv"),
  make_option("--refRNA",    type = "character", default = NULL,
              help = "Seurat .rds of the Chen/Lau polyp epithelial atlas"),
  make_option("--refCellTypeCol", type = "character", default = "CellType"),
  make_option("--minScore",  type = "double",    default = 0.5,
              help = "min predictedScore to accept an ASC/SSC transfer [default %default]"),
  make_option("--threads",   type = "integer",   default = 16),
  make_option("--project",   type = "character", default = "proj_epithelial_lau")
)))

archr    <- normalizePath(opt$archr)
annot    <- normalizePath(opt$annot)
manifest <- read.csv(normalizePath(opt$manifest), stringsAsFactors = FALSE)
labelMap <- read.csv(normalizePath(opt$labelMap), stringsAsFactors = FALSE)

setwd(archr)
addArchRThreads(opt$threads)
addArchRGenome("hg38")
`%ni%` <- Negate(`%in%`)

## ------------------------------------------------------- 1) Becker labels -> project
proj <- loadArchRProject("proj_all_cells")

projections <- read.csv(file.path(annot, "all_samples_cell_types_projections.csv"),
                        row.names = 1, stringsAsFactors = FALSE)
hubmap <- read.table(file.path(annot, "HuBMAPCellAnnotations_hubmap_colon_epithelial.tsv"),
                     sep = "\t", header = TRUE, row.names = 1, stringsAsFactors = FALSE)
colnames(hubmap) <- c("Cell", "CellType")

# HuBMAP normals were the projection reference; keep their native labels, drop the copies
normalSamples      <- manifest$ATAC_Sample[manifest$Condition == "Normal"]
projections$Sample <- sub("#.*", "", projections$Cell)
projections        <- projections[projections$Sample %ni% normalSamples, c("Cell", "CellType")]

labels <- rbind(projections, hubmap[, c("Cell", "CellType")])
labels <- labels[!duplicated(labels$Cell), ]
labels <- labels[labels$Cell %in% proj$cellNames, ]
if (nrow(labels) == 0) stop("No annotation cell names matched. Run check_annotation_join.R.")

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

# TUF is absent by construction: Becker et al. annotated 11 epithelial types in scATAC and
# Tuft was not among them (it appears only on the snRNA axis of their confusion matrix).
message("NOTE: TUF has no scATAC counterpart in this dataset and is not recoverable here.")

## ------------------------------------------------------- 3) sample metadata
m <- manifest[match(proj$Sample, manifest$ATAC_Sample), ]
for (col in c("Donor", "FAP", "DiseaseState", "Condition", "Location")) {
  proj <- addCellColData(proj, data = as.character(m[[col]]),
                         cells = proj$cellNames, name = col, force = TRUE)
}

## ------------------------------------------------------- 4) ASC by label transfer
# Becker projected polyp cells into the *normal* manifold, so every adenoma cell currently
# carries a normal-like label. ASC does not exist in their annotation and must be
# transferred in from a reference that defines it.
if (is.null(opt$refRNA)) {
  warning("--refRNA not supplied: skipping ASC transfer. Adenoma cells keep normal-like labels.")
} else {
  adenomaCells <- proj$cellNames[proj$Condition == "Adenoma"]
  if (length(adenomaCells) < 100) stop("Too few adenoma cells for integration.")

  projAd <- proj[adenomaCells, ]

  # Integration needs a reduced dim computed on the adenoma cells themselves
  projAd <- addIterativeLSI(projAd, useMatrix = "TileMatrix", name = "LSI_adenoma",
                            iterations = 4, varFeatures = 15000, dimsToUse = 1:25,
                            clusterParams = list(resolution = c(0.1, 0.1, 0.2),
                                                 sampleCells = 10000, n.start = 10),
                            sampleCellsPre = NULL, force = TRUE)

  seRNA <- readRDS(opt$refRNA)
  ctCol <- opt$refCellTypeCol
  if (ctCol %ni% colnames(seRNA@meta.data))
    stop("Reference lacks metadata column '", ctCol, "'. Columns: ",
         paste(colnames(seRNA@meta.data), collapse = ", "))
  if ("ASC" %ni% unique(seRNA@meta.data[[ctCol]]))
    stop("Reference '", ctCol, "' has no 'ASC' level. Levels: ",
         paste(unique(seRNA@meta.data[[ctCol]]), collapse = ", "))

  projAd <- addGeneIntegrationMatrix(
    ArchRProj    = projAd,
    useMatrix    = "GeneScoreMatrix",
    matrixName   = "GeneIntegrationMatrix",
    reducedDims  = "LSI_adenoma",
    seRNA        = seRNA,
    addToArrow   = FALSE,
    groupRNA     = ctCol,
    nameCell     = "predictedCell_Lau",
    nameGroup    = "predictedGroup_Lau",
    nameScore    = "predictedScore_Lau",
    force        = TRUE
  )

  transferred <- data.frame(
    Cell  = projAd$cellNames,
    Group = projAd$predictedGroup_Lau,
    Score = projAd$predictedScore_Lau,
    stringsAsFactors = FALSE
  )
  write.csv(transferred, file.path(archr, "adenoma_label_transfer_Lau.csv"), row.names = FALSE)

  cat("\nTransferred label distribution (adenoma cells):\n")
  print(table(transferred$Group))
  cat("\nMedian predictedScore by group:\n")
  print(tapply(transferred$Score, transferred$Group, median))

  # Conservative: only OVERRIDE with the tumour-specific labels (ASC/SSC), and only when the
  # transfer is confident. Every other adenoma cell keeps its Becker-derived ATAC label,
  # which is more trustworthy than a cross-modality transfer for the common lineages.
  tumourLabels <- intersect(c("ASC", "SSC"), unique(transferred$Group))
  accept <- transferred$Group %in% tumourLabels & transferred$Score >= opt$minScore

  message(sprintf("Accepting %d / %d adenoma cells as %s (predictedScore >= %.2f).",
                  sum(accept), nrow(transferred),
                  paste(tumourLabels, collapse = "/"), opt$minScore))

  # NB: addCellColData() sets every cell NOT in `cells` to NA. So rebuild the full vector
  # rather than passing a subset -- otherwise all normal cells silently lose their labels.
  if (sum(accept) > 0) {
    newLabels <- proj$CellTypeLau
    names(newLabels) <- proj$cellNames
    newLabels[transferred$Cell[accept]] <- transferred$Group[accept]
    proj <- addCellColData(proj, data = unname(newLabels), cells = proj$cellNames,
                           name = "CellTypeLau", force = TRUE)
  }
}

## ------------------------------------------------------- 5) summary + save
cat("\nFinal cells per Condition x CellTypeLau:\n")
print(table(proj$Condition, proj$CellTypeLau))

write.csv(as.data.frame(table(proj$Condition, proj$CellTypeLau)),
          file.path(archr, "celltype_counts_lau.csv"), row.names = FALSE)

proj$PeakCallingGroup <- paste(proj$Condition, proj$CellTypeLau, sep = "__")
saveArchRProject(proj, outputDirectory = opt$project, load = FALSE)
message("\nSaved ", opt$project)
