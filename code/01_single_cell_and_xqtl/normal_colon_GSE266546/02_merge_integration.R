suppressPackageStartupMessages({
  library(Seurat)        # >= 4.3 (works with v5 too)
  library(Matrix)
  library(dplyr)
  library(stringr)
  library(harmony)
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(R.utils)
  library(future)
})


set.seed(123)

# ---- paths ------------------------------------------------------------------
# Locate script directory robustly (works under Rscript and source())
.find_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  hit  <- grep("--file=", args, value = TRUE)
  if (length(hit) > 0) return(dirname(normalizePath(sub("--file=", "", hit[1]))))
  if (!is.null(sys.frames()) && length(sys.frames()) > 0) {
    f <- tryCatch(sys.frame(1)$ofile, error = function(e) NULL)
    if (!is.null(f)) return(dirname(normalizePath(f)))
  }
  getwd()
}
SCRIPT_DIR <- .find_script_dir()

ROOT <- "/data/l2_bioinfo1/lyul1/single_cell_data/IBD/"
RESULTS_DIR <- file.path(ROOT, "Results")
dir.create(RESULTS_DIR, showWarnings = FALSE, recursive = TRUE)


obj266<- readRDS(file.path(RESULTS_DIR, "obj266_preqc_pref_AC.rds")) # obj266_qc_pref_AC.rds
objH<- readRDS(file.path(RESULTS_DIR, "objH_preqc.rds"))

message(
  "    ", nrow(obj266), " genes x ",
  ncol(obj266), " cells"
)
message(
  "    ", length(unique(obj266$SampleID)), " samples, ",
  length(unique(obj266$PatientID)), " patients"
)

message(
  "    ", nrow(objH), " genes x ",
  ncol(objH), " cells"
)
message(
  "    ", length(unique(objH$SampleID)), " samples, ",
  length(unique(objH$PatientID)), " patients"
)
# ============================================================================
# 4) Keep only epithelial cells
# ============================================================================

table(obj266$Celltypes, useNA = "ifany")
table(objH$Cell_Type, useNA = "ifany")
table(objH$Cell.Type, useNA = "ifany")
table(objH$CellType_cb, useNA = "ifany")


immune_labels <- c(
  "Plasma",
  "Plasma Proliferating",
  "NK",
  "T",
  "T Proliferating",
  "B",
  "B Proliferating",
  "Resident Macro.",
  "Recruited Macro.",
  "Neutrophils",
  "Mast",
  "20",
  "Fibro.",
  "ABCA8+ Fibro.",
  "PDGFRA+ Fibro.",
  "PDPN+ Fibro.",
  "Endo.",
  "LND",
  "SSC",
  "PLA", "MYE", 
  "END",
  "FIB", "MAS", "ASC"
)

obj266 <- subset(
  obj266,
  subset = !(Celltypes %in% immune_labels)
)


message("GSE after immune exclusion:")
print(table(obj266$Celltypes, useNA = "ifany"))
message("  ", ncol(obj266), " cells, ",
        length(unique(obj266$PatientID)), " patients")
message(
  "    ", nrow(obj266), " genes x ",
  ncol(obj266), " cells"
)
message(
  "    ", length(unique(obj266$SampleID)), " samples, ",
  length(unique(obj266$PatientID)), " patients"
)

saveRDS(obj266, file.path(RESULTS_DIR, "obj266_pref_AC_immu.rds"))

write.csv(obj266@meta.data, file.path(RESULTS_DIR, "obj266_metadata_pref_AC_immu.csv"), row.names = F)

patients_before <- unique(objH$PatientID)

#test<- objH@meta.data %>%
#  filter(PatientID == "HTA11_5212") 
  
keep_cells <- rep(TRUE, ncol(objH))

if ("Cell_Type" %in% colnames(objH@meta.data)) {
  keep_cells <- keep_cells & !(objH$Cell_Type %in% immune_labels)
}

if ("Cell.Type" %in% colnames(objH@meta.data)) {
  keep_cells <- keep_cells & !(objH$Cell.Type %in% immune_labels)
}

objH <- subset(objH, cells = colnames(objH)[keep_cells])


saveRDS(objH, file.path(RESULTS_DIR, "objH_immu.rds"))

write.csv(objH@meta.data, file.path(RESULTS_DIR, "objH_metadata_immu.csv"), row.names = F)

message("HTAN after immune exclusion:")
print(table(objH$Cell_Type, useNA = "ifany"))
print(table(objH$Cell.Type, useNA = "ifany"))
print(table(objH$CellType_cb, useNA = "ifany"))

message(
  "    ", nrow(objH), " genes x ",
  ncol(objH), " cells"
)
message(
  "    ", length(unique(objH$SampleID)), " samples, ",
  length(unique(objH$PatientID)), " patients"
)

patients_after  <- unique(objH$PatientID)

lost_patients <- setdiff(patients_before, patients_after)
lost_patients

# QC plots before integration immu
qc_pdf <- file.path(RESULTS_DIR, "00_qc_violins_ex_imm.pdf")
pdf(qc_pdf, width = 14, height = 8)
print(VlnPlot(obj266, c("nFeature_RNA", "nCount_RNA", "pct_mito"),
              pt.size = 0, ncol = 3) + plot_annotation("GSE266546 (post-QC)"))
print(VlnPlot(objH,   c("nFeature_RNA", "nCount_RNA", "pct_mito"),
              pt.size = 0, ncol = 3) + plot_annotation("HTAN Islam 2024 (post-QC)"))
dev.off()



# ============================================================================
# 5) Merge, normalize (SCT), PCA, Harmony
# ============================================================================
message("\n[4/5] Merge + SCTransform + PCA + Harmony...")


plan(sequential)
options(future.globals.maxSize = 10 * 1024^3)

merged <- merge(obj266, y = objH, add.cell.ids = c("V266", "HTAN"),
                project = "IBD")

# Ensure pct_mito exists before SCTransform
merged[["pct_mito"]] <- PercentageFeatureSet(merged, pattern = "^MT-")

saveRDS(merged, file.path(RESULTS_DIR, "merged.rds"))
write.csv(merged@meta.data, file.path(RESULTS_DIR, "merged_metadata.csv"), row.names = F)

message(
  "    ", nrow(merged), " genes x ",
  ncol(merged), " cells"
)
message(
  "    ", length(unique(merged$SampleID)), " samples, ",
  length(unique(merged$PatientID)), " patients"
)


rm(obj266, objH); gc()

# Split-SCT-merge keeps the normalization study-aware (preferred)
merged_list <- SplitObject(merged, split.by = "study")

print(sapply(merged_list, function(o) "pct_mito" %in% colnames(o@meta.data)))

merged_list <- lapply(merged_list, function(o)
  SCTransform(o, verbose = FALSE, vars.to.regress = "pct_mito"))

# pick integration features then re-merge SCT slots
features <- SelectIntegrationFeatures(merged_list, nfeatures = 3000)
merged   <- merge(merged_list[[1]], merged_list[-1])
VariableFeatures(merged) <- features
rm(merged_list); gc()

merged <- RunPCA(merged, npcs = 50, verbose = FALSE)

message(
  "    ", nrow(merged), " genes x ",
  ncol(merged), " cells"
)
message(
  "    ", length(unique(merged$SampleID)), " samples, ",
  length(unique(merged$PatientID)), " patients"
)

length(rownames(merged[["RNA"]]))
length(rownames(merged[["SCT"]]))

# UMAP BEFORE batch correction (for the comparison panel)
merged <- RunUMAP(merged, dims = 1:30, reduction = "pca",
                  reduction.name = "umap_pca", verbose = FALSE)

table(merged$batch)
table(merged$study, merged$batch)

message(
  "    ", nrow(merged), " genes x ",
  ncol(merged), " cells"
)
message(
  "    ", length(unique(merged$SampleID)), " samples, ",
  length(unique(merged$PatientID)), " patients"
)

# Harmony on the "study" covariate (the batch we care about). If you also
# want to correct sample-to-sample drift, add "sample".
# merged <- RunHarmony(merged, group.by.vars = "batch",
#                      assay.use = "SCT", reduction = "pca",
#                      dims.use = 1:30, plot_convergence = FALSE)

DefaultAssay(merged) <- "SCT"

merged <- RunHarmony(
  object = merged,
  group.by.vars = "batch",
  reduction.use = "pca",
  dims.use = 1:30,
  reduction.save = "harmony",
  plot_convergence = FALSE
)

# UMAP after Harmony
merged <- RunUMAP(merged, dims = 1:30, reduction = "harmony",
                  reduction.name = "umap_harmony", verbose = FALSE)

# Default reduction = harmony UMAP for downstream (Seurat-v5-safe)
merged[["umap"]] <- merged[["umap_harmony"]]

message(
  "    ", nrow(merged), " genes x ",
  ncol(merged), " cells"
)
message(
  "    ", length(unique(merged$SampleID)), " samples, ",
  length(unique(merged$PatientID)), " patients"
)

length(rownames(merged[["RNA"]]))
length(rownames(merged[["SCT"]]))

# ============================================================================
# 5) Save intermediate object
# ============================================================================
message("\n[5/5] Saving integrated (pre-clustering) object...")
saveRDS(merged, file.path(RESULTS_DIR, "integrated_seurat_postharmony.rds"))
message("  Wrote ", file.path(RESULTS_DIR, "integrated_seurat_postharmony.rds"))
message("Done.")



