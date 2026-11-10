#!/usr/bin/env Rscript
# 03_cluster_annot.R
# ----------------------------------------------------------------------------
# Clusters the Harmony-integrated Seurat object, exports marker and QC
# summaries, and saves the annotated checkpoint for subclustering.
# ----------------------------------------------------------------------------
options(bitmapType = "cairo")
#install.packages("pheatmap", repos = "https://cloud.r-project.org")
suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(R.utils)
  library(pheatmap)
  
})

set.seed(123)

# ---- paths ------------------------------------------------------------------
args <- commandArgs(trailingOnly = FALSE)

ROOT <- "/data/l2_bioinfo1/lyul1/single_cell_data/IBD"
RESULTS_DIR <- file.path(ROOT, "Results")

FIG_DIR     <- file.path(RESULTS_DIR, "figures"); dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)

rds_in <- file.path(RESULTS_DIR, "integrated_seurat_postharmony.rds") # integrated_seurat_postharmony.rds integrated_seurat_final1.res0.8.rm23.rds integrated_seurat_postharmony_pref_AC.rds
stopifnot(file.exists(rds_in))
message("Loading ", rds_in)
obj <- readRDS(rds_in)

# ============================================================================
# Clustering
# ============================================================================
message("FindNeighbors on harmony embedding...")
obj <- FindNeighbors(obj, reduction = "harmony", dims = 1:30, verbose = FALSE)

message("FindClusters at multiple resolutions...")
for (res in c(0.4, 0.8, 1.2)) {
  obj <- FindClusters(obj, resolution = res, verbose = FALSE)
}

Idents(obj) <- "SCT_snn_res.0.8" # SCT_snn_res.0.4 0.8
obj$cluster <- Idents(obj)

# Remove all cells in cluster 23 before downstream summaries/figures/exports
obj$cluster <- as.character(obj$cluster)

message("Clusters before removing cluster 23:")
print(table(obj$cluster, useNA = "ifany"))

obj <- subset(
  obj,
  subset = cluster != "23"
)
stopifnot(!"23" %in% obj$cluster)
Idents(obj) <- "cluster"

message("Clusters after removing cluster 23:")
print(table(obj$cluster, useNA = "ifany"))

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
write.csv(markers, file.path(RESULTS_DIR, "cluster_markers_res0.8.rm23.csv"),
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
  save_plot(p_hm, "07_top5_markers_heatmap_res0.8.rm23", w = 12, h = 14)
}

# Cell-count tables
tab_study   <- table(obj$cluster, obj$study)
tab_cond    <- table(obj$cluster, obj$sample_condition)
write.csv(tab_study, file.path(RESULTS_DIR, "cells_per_cluster_by_study.res0.8.rm23.csv"))
write.csv(tab_cond,  file.path(RESULTS_DIR, "cells_per_cluster_by_condition.res0.8.rm23.csv"))

## annotation

obj$CellType_cb2 <- as.character(obj$CellType_cb)

obj$CellType_cb2[grepl("Goblet", obj$CellType_cb2, ignore.case = TRUE)] <- "GOB"
obj$CellType_cb2[grepl("Stem", obj$CellType_cb2, ignore.case = TRUE)] <- "STM"
obj$CellType_cb2[grepl("Tuft", obj$CellType_cb2, ignore.case = TRUE)] <- "TUF"
obj$CellType_cb2[grepl("EEC|Enteroendocrine", obj$CellType_cb2, ignore.case = TRUE)] <- "EE"

#obj$CellType_cb_simple[grepl("BEST4|OTOP2", obj$CellType_cb_simple, ignore.case = TRUE)] <- "CT"

table(obj$CellType_cb, obj$cluster)
a1<- table(obj$CellType_cb2, obj$cluster) %>% as.data.frame() %>% group_by(Var2) %>% arrange(desc(Freq), .by_group = TRUE)
# write.csv(a1, file.path(RESULTS_DIR, "cluster_by_original_CellType_cb_res0.8.rm23.csv"), row.names = FALSE)


obj$CellType_cb3 <- as.character(obj$CellType_cb2)
obj$CellType_cb3[grepl("Colonocyte|Enterocyte", obj$CellType_cb3, ignore.case = TRUE)] <- "ABS"
obj$CellType_cb3[grepl("Cycling TA", obj$CellType_cb3, ignore.case = TRUE)] <- "TAC"
a2<- table(obj$CellType_cb3, obj$cluster) %>% as.data.frame() %>% group_by(Var2) %>% arrange(desc(Freq), .by_group = TRUE)
# write.csv(a2, file.path(RESULTS_DIR, "cluster_by_original_CellType_cb3_res0.8.rm23.csv"), row.names = FALSE)

obj$celltype_annot <- NA_character_

abs_clusters <- c("0", "1", "11", "13", "14", "15", "4", 
                  "16", "17", "19", "21", "22", "24", "25", "6", "7", "8", "9")

obj$celltype_annot[as.character(obj$cluster) %in% abs_clusters] <- "ABS"
obj$celltype_annot[as.character(obj$cluster) == "18"] <- "TUF"
obj$celltype_annot[as.character(obj$cluster) == "20"]  <- "EE"
obj$celltype_annot[as.character(obj$cluster) %in% c("3")]  <- "TAC" # , "23"
obj$celltype_annot[as.character(obj$cluster) %in% c("2", "5")] <- "GOB"
obj$celltype_annot[as.character(obj$cluster) %in% c("10")]  <- "CT"
obj$celltype_annot[as.character(obj$cluster) == "12"]  <- "STM"

a3<- table(obj$CellType_cb3, obj$celltype_annot) %>% as.data.frame() %>% group_by(Var2) %>% arrange(desc(Freq), .by_group = TRUE)
# write.csv(a3, file.path(RESULTS_DIR, "cluster_by_original_CellType_cb_annot_res0.8.rm23.csv"), row.names = FALSE)

table(obj$cluster, obj$celltype_annot, useNA = "ifany")

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
  "    ", length(unique(obj$celltype_annot)), " annotated clusters, ",
  "    ", length(unique(obj$CellType_cb)), " original clusters, ")

# Step 04 matches this checkpoint against merged.rds using the original cell IDs.
saveRDS(obj, file.path(RESULTS_DIR, "integrated_seurat_final.res0.8.rm23.rds"))

message("Done.  Outputs in: ", RESULTS_DIR)
