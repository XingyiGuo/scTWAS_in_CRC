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
message("Seurat version: ", packageVersion("Seurat"))
message("Harmony version: ", packageVersion("harmony"))
stopifnot(exists("SplitObject"))
set.seed(123)

source("/data/l2_bioinfo1/liq17/scTWAS/data/VUMC_HTAN_GSE266546/src/renv/activate.R")

# ---- paths ------------------------------------------------------------------
# ROOT <- "/data/l2_bioinfo1/lyul1/single_cell_data/IBD/"
# RESULTS_DIR <- file.path(ROOT, "Results")
# dir.create(RESULTS_DIR, showWarnings = FALSE, recursive = TRUE)
# FIG_DIR     <- file.path(RESULTS_DIR, "figures"); dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)
####
#obj<- readRDS(file.path(RESULTS_DIR, "integrated_seurat_postharmony_sub.rds"))
RESULTS_DIR<-"/data/l2_bioinfo1/liq17/scTWAS/data/VUMC_HTAN_GSE266546/raw/rds2"
FIG_DIR     <- file.path(RESULTS_DIR, "figures"); dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)
merged<- readRDS(file.path(RESULTS_DIR, "merged.rds")) #
prev_obj<- readRDS(file.path(RESULTS_DIR, "integrated_seurat_final.res0.8.rm23.rds"))

#keep_clusters <- c(0, 1, 3, 4, 6, 8, 9, 11, 14,  21,  22, 24, 12, 15, 13) #circle to distinguish ABS/STM/TAC
# ABS = c(0, 1, 4, 6, 7, 8, 9, 11, 13, 14, 15, 16, 17, 19, 21, 22, 24, 25), - (7, 9, 15, 16, 17, 19, 25)
# GOB = c(2, 5),
# STM = c(12),
# TAC = c(3),
# TUF = c(18),
# CT  = c(10),
# EE  = c(20)
keep_clusters <- c(0, 1, 4, 6, 8, 11, 13, 14, 21, 22, 24, 2, 5, 12, 3, 18, 10, 20)
keep_clusters <- as.character(unique(keep_clusters))

Idents(prev_obj) <- "cluster"
cells_keep <- WhichCells(
  prev_obj,
  idents = keep_clusters
)

# keep only cells that are also present in merge_obj
cells_keep <- intersect(cells_keep, colnames(merged))

# subset from the original merged object
merged <- subset(
  merged,
  cells = cells_keep
)

message("Batch table before Harmony:")
print(table(merged$batch, useNA = "ifany"))

message("Study x batch table before Harmony:")
print(table(merged$study, merged$batch, useNA = "ifany"))

saveRDS(merged, file.path(RESULTS_DIR, "merged_sub.rds"))
write.csv(merged@meta.data, file.path(RESULTS_DIR, "merged_sub_metadata.csv"), row.names = F)

message(
  "    ", nrow(merged), " genes x ",
  ncol(merged), " cells"
)
message(
  "    ", length(unique(merged$SampleID)), " samples, ",
  length(unique(merged$PatientID)), " patients"
)
rm(prev_obj)

        # GSE266546 HTAN_PCA83samples         HTAN_VUMC 
            # 45513             17552             57826 
# Study x batch table before Harmony:
                
                 # GSE266546 HTAN_PCA83samples HTAN_VUMC
  # GSE266546          45513                 0         0
  # HTAN_islam2024         0             17552     57826
    # 32261 genes x 120891 cells
    # 157 samples, 117 patients

# ============================================================================
# 5) SCT -> PCA -> Harmony -> UMAP
# merged
# ├── RNA
# ├── SCT
# ├── pca
# ├── umap_pca       # before Harmony
# ├── harmony
# └── umap_harmony   # after Harmony
#
# merged
   # │
   # ├── GSE266546
   # │      ↓
   # │   SCTransform
   # │   regress pct_mito
   # │
   # └── HTAN_islam2024
          # ↓
       # SCTransform
       # regress pct_mito
          # ↓
# SelectIntegrationFeatures(3000)
          # ↓
# merge
          # ↓
# PCA
RESULTS_DIR <- "/data/l2_bioinfo1/liq17/scTWAS/data/VUMC_HTAN_GSE266546/raw/rds2"
INPUT_RDS  <- file.path(RESULTS_DIR, "merged_sub.rds")
PCA_RDS    <- file.path(RESULTS_DIR, "merged_postSCT_PCA_sub.rds")
FINAL_RDS  <- file.path(RESULTS_DIR, "integrated_seurat_postharmony_sub.rds")
N_CPU <- 6L
# ----------------------------------------------------------------------------
# Prevent nested BLAS/OpenMP threading
# ----------------------------------------------------------------------------
Sys.setenv(
  OMP_NUM_THREADS = "1",
  OPENBLAS_NUM_THREADS = "1",
  MKL_NUM_THREADS = "1",
  BLIS_NUM_THREADS = "1"
)
options(future.globals.maxSize = 10 * 1024^3)
future::plan(future::sequential)

# ============================================================================
# 5A) Load
# ============================================================================
message("\n[4/5] Loading object...")
merged <- readRDS(INPUT_RDS)
message("    ", nrow(merged), " genes x ", ncol(merged), " cells")
message("    ", length(unique(merged$SampleID)), " samples, ", length(unique(merged$PatientID)), " patients")
print(table(merged$study, merged$batch))

# ============================================================================
# 5B) Split + SCTransform
# ============================================================================
message("\n[4/5] SCTransform by study...")
# Prevent uncontrolled OpenMP / BLAS threads
Sys.setenv(
  OMP_NUM_THREADS = "1",
  OPENBLAS_NUM_THREADS = "1",
  MKL_NUM_THREADS = "1",
  BLIS_NUM_THREADS = "1"
)

# Large Seurat objects need a higher future export limit
options(future.globals.maxSize = 10 * 1024^3)

# Avoid nested parallelism during SCT
future::plan(future::sequential)
merged_list <- Seurat::SplitObject(
  merged,
  split.by = "study"
)
print(sapply(merged_list, ncol))

# Check pct_mito
stopifnot(
  all(vapply(
    merged_list,
    function(x) "pct_mito" %in% colnames(x@meta.data),
    logical(1)
  ))
)
# Run studies ONE AT A TIME
for (i in seq_along(merged_list)) {
  message("    SCTransform: ",names(merged_list)[i]," (",ncol(merged_list[[i]])," cells)")
  merged_list[[i]] <- Seurat::SCTransform(
    merged_list[[i]],
    vars.to.regress = "pct_mito",
    verbose = TRUE
  )
  gc()
}
message("    SCTransform complete.")

# ============================================================================
# 5C) Select features + merge
# ============================================================================
message("\n[4/5] Selecting integration features and merging...")
features <- Seurat::SelectIntegrationFeatures(object.list = merged_list,nfeatures = 3000)
merged <- merge(merged_list[[1]],y = merged_list[-1])
DefaultAssay(merged) <- "SCT"
VariableFeatures(merged) <- features
rm(merged_list)
gc()
message("    ", nrow(merged), " genes x ", ncol(merged), " cells")
message("    RNA genes: ", nrow(merged[["RNA"]]), "; SCT genes: ", nrow(merged[["SCT"]]))

# ============================================================================
# 5D) PCA
# ============================================================================
message("\n[4/5] PCA...")
merged <- Seurat::RunPCA(
  merged,
  assay = "SCT",
  features = features,
  npcs = 50,
  verbose = FALSE
)
message(
  "    PCA: ",
  nrow(Embeddings(merged, "pca")),
  " cells x ",
  ncol(Embeddings(merged, "pca")),
  " PCs"
)

# ----------------------------------------------------------------------------
# IMPORTANT checkpoint
# ----------------------------------------------------------------------------
message("    Saving PCA checkpoint...")
saveRDS(merged,PCA_RDS,compress = FALSE)
gc()

message("\n[4/5] Loading object after PCA...")
merged <- readRDS(PCA_RDS)
# ============================================================================
# 5E) UMAP BEFORE Harmony
# ============================================================================
message("\n[4/5] UMAP before Harmony...")
set.seed(123)
merged <- Seurat::RunUMAP(
  merged,
  reduction = "pca",
  dims = 1:30,
  reduction.name = "umap_pca",
  reduction.key = "UMAPPCA_",
  umap.method = "uwot",
  metric = "cosine",
  n_threads = N_CPU,
  n_sgd_threads = 1,
  seed.use = 123,
  verbose = FALSE
)
gc()

# ============================================================================
# 5F) Harmony batch correction
# ============================================================================
message("\n[4/5] Harmony batch correction...")
print(table(merged$study, merged$batch))
merged <- harmony::RunHarmony(
  object = merged,
  group.by.vars = "batch",
  reduction.use = "pca",
  dims.use = 1:30,
  reduction.save = "harmony",
  ncores = 1,
  project.dim = FALSE,
  plot_convergence = FALSE,
  verbose = TRUE
)
message(
  "    Harmony: ",
  nrow(Embeddings(merged, "harmony")),
  " cells x ",
  ncol(Embeddings(merged, "harmony")),
  " dimensions"
)
gc()


# ============================================================================
# 5G) UMAP AFTER Harmony
# ============================================================================
message("\n[4/5] UMAP after Harmony...")
set.seed(123)
merged <- Seurat::RunUMAP(
  merged,
  reduction = "harmony",
  dims = 1:30,
  reduction.name = "umap_harmony",
  reduction.key = "UMAPHARMONY_",
  umap.method = "uwot",
  metric = "cosine",
  n_threads = N_CPU,
  n_sgd_threads = 1,
  seed.use = 123,
  verbose = FALSE
)
gc()

# ============================================================================
# 5H) Standard UMAP = Harmony UMAP
# ============================================================================
merged[["umap"]] <- merged[["umap_harmony"]]
Key(merged[["umap"]]) <- "UMAP_"

# ============================================================================
# 5I) Final QC
# ============================================================================
message("\nFinal reductions:")
print(Reductions(merged))
message("    ", nrow(merged), " genes x ",ncol(merged), " cells")
message("    ", length(unique(merged$SampleID)), " samples, ",length(unique(merged$PatientID)), " patients")
message(
  "    RNA genes: ",
  nrow(merged[["RNA"]]),
  "; SCT genes: ",
  nrow(merged[["SCT"]])
)
print(table(merged$study, merged$batch))

# ============================================================================
# 5J) Save
# ============================================================================
message("\n[5/5] Saving integrated object...")
saveRDS(merged,FINAL_RDS)
message("    Wrote: ", FINAL_RDS)
message("Done.")

# ============================================================================
# 6. Clustering
# ============================================================================
message("FindNeighbors on harmony embedding...")
obj<- merged
obj <- FindNeighbors(obj, reduction = "harmony", dims = 1:30, verbose = FALSE)
message("FindClusters at multiple resolutions...")
for (res in c(0.4, 0.8, 1.2)) {obj <- FindClusters(obj, resolution = res, verbose = FALSE)}
Idents(obj) <- "SCT_snn_res.0.8" # SCT_snn_res.0.4 0.8
obj$cluster <- Idents(obj)

# ============================================================================
# Per-cluster marker genes
# ============================================================================
message("FindAllMarkers...")
Idents(obj) <- "cluster"
DefaultAssay(obj) <- "SCT"
obj <- PrepSCTFindMarkers(obj)
markers <- FindAllMarkers(obj,
                          only.pos      = TRUE,
                          min.pct       = 0.25,
                          logfc.threshold = 0.25, recorrect_umi = FALSE,
                          verbose       = FALSE)
write.csv(markers, file.path(RESULTS_DIR, "cluster_markers_res0.8.rm23_sub.csv"),
          row.names = FALSE)


# Top-5 marker heatmap
### saving figures
save_plot <- function(p, name, w = 8, h = 6) {
  ggsave(
    filename = file.path(FIG_DIR, paste0(name, ".pdf")),
    plot = p,
    width = w,
    height = h,
    device = cairo_pdf
  )
}

top5 <- markers %>% group_by(cluster) %>% slice_max(avg_log2FC, n = 5)
if (nrow(top5) > 0) {
  p_hm <- DoHeatmap(obj, features = unique(top5$gene)) + NoLegend()
  save_plot(p_hm, "07_top5_markers_heatmap_res0.8.rm23_sub", w = 12, h = 14)
}


# Cell-count tables
tab_study   <- table(obj$cluster, obj$study)
tab_cond    <- table(obj$cluster, obj$sample_condition)
write.csv(tab_study, file.path(RESULTS_DIR, "cells_per_cluster_by_study.res0.8.rm23_sub.csv"))
write.csv(tab_cond,  file.path(RESULTS_DIR, "cells_per_cluster_by_condition.res0.8.rm23_sub.csv"))


# save final object
message(
  "    ", nrow(obj), " genes x ",
  ncol(obj), " cells"
)
message(
  "    ", length(unique(obj$SampleID)), " samples, ",
  length(unique(obj$PatientID)), " patients"
)

message(
  "    ", length(unique(obj$CellType_cb)), " original clusters, ")

saveRDS(obj, file.path(RESULTS_DIR, "integrated_seurat_final.res0.8.rm23_sub.rds"))

write.csv(obj@meta.data, file.path(RESULTS_DIR, "integrated_seurat_final.res0.8.rm23_metadata_sub.csv"), row.names = F)


## clean cols and save final objs
new_cell_names <- colnames(obj)

new_cell_names <- gsub("_h5ad_", "_", new_cell_names)
new_cell_names <- sub("^h5ad_", "", new_cell_names)

new_cell_names <- make.unique(new_cell_names)

obj <- RenameCells(obj, new.names = new_cell_names)


saveRDS(obj, file.path(RESULTS_DIR, "integrated_seurat_final1.res0.8.rm23_sub.rds"))
write.csv(obj@meta.data, file.path(RESULTS_DIR, "integrated_seurat_final1.res0.8.rm23_sub_metadata.csv"), row.names = F)



# Standard UMAPs
save_plot(DimPlot(obj, reduction = "umap_harmony", group.by = "study",     raster = FALSE),
           "02_umap_by_study_res0.8.rm23_sub")
save_plot(DimPlot(obj, reduction = "umap_harmony", group.by = "batch",    raster = FALSE) + NoLegend(),
          "03_umap_by_sample_res0.8.rm23_sub")
save_plot(DimPlot(obj, reduction = "umap_harmony", group.by = "sample_condition", raster = FALSE),
          "04_umap_by_condition_res0.8.rm23_sub")
save_plot(DimPlot(obj, reduction = "umap_harmony", group.by = "cluster",
                  label = TRUE, raster = FALSE) + NoLegend(),
         "05_umap_by_cluster_res0.8.rm23_sub")

for (res in c(0.8)) {
  col <- paste0("SCT_snn_res.", res)
  save_plot(DimPlot(obj, reduction = "umap_harmony", group.by = "cluster",
                    label = TRUE, raster = FALSE) + NoLegend() +
               ggtitle(paste0("Clusters @ res ", res)),
            paste0("05b_umap_clusters_res", res, ".rm23_sub"))
 }




message("Done.  Outputs in: ", RESULTS_DIR)

### write to count matrix
DefaultAssay(obj) <- "RNA"

#Layers(obj[["RNA"]])
obj <- JoinLayers(obj, assay = "RNA")

counts_mat <- GetAssayData(obj, assay = "RNA", layer = "counts")
saveRDS(counts_mat, file.path(RESULTS_DIR, "counts_matrix_res0.8.rm23_sub.rds"))
