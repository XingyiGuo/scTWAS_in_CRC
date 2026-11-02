#!/usr/bin/env Rscript
# 03_cluster_plots.R
# ----------------------------------------------------------------------------
# Loads the Harmony-integrated Seurat object, runs clustering, makes UMAP
# panels (study / sample / condition / cluster), computes per-cluster
# markers, and writes everything to ../results/.
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

#obj$celltype_annot <- NA_character_
#abs_clusters <- c( "1", "10", "11", "13", "14", "4", 
#                   "17", "18", "19", "20", "21", "22", "24", "25", "27", "28", "30", "33", "4", "5", "6", "7", "8", "9")
#obj$celltype_annot[as.character(obj$cluster) %in% abs_clusters] <- "ABS"
#obj$celltype_annot[as.character(obj$cluster) == "23"] <- "TUF"
#obj$celltype_annot[as.character(obj$cluster) == "26"]  <- "EE"
#obj$celltype_annot[as.character(obj$cluster) %in% c("0", "31", "32")]  <- "TAC"
#obj$celltype_annot[as.character(obj$cluster) %in% c("16", "2", "29", "3")] <- "GOB"
#obj$celltype_annot[as.character(obj$cluster) %in% c("12")]  <- "CT"
#obj$celltype_annot[as.character(obj$cluster) == "15"]  <- "STM" # STM 1604 TAC 1247

#a3<- table(obj$CellType_cb3, obj$celltype_annot) %>% as.data.frame() %>% group_by(Var2) %>% arrange(desc(Freq), .by_group = TRUE)
#write.csv(a3, file.path(RESULTS_DIR, "cluster_by_original_CellType_cb_annot_res1.2.csv"), row.names = FALSE)

a4 <- obj@meta.data %>%
  mutate(
    cluster = as.character(cluster),
    CellType_cb = as.character(CellType_cb),
    CellType_cb3 = as.character(CellType_cb3),
    celltype_annot = as.character(celltype_annot)
  ) %>%
  group_by(cluster, CellType_cb, CellType_cb3, celltype_annot) %>%
  summarise(n_cells = n(), .groups = "drop") %>%
  group_by(cluster) %>%
  mutate(
    cluster_total = sum(n_cells),
    percent_within_cluster = round(100 * n_cells / cluster_total, 2)
  ) %>%
  arrange(as.numeric(cluster), desc(percent_within_cluster)) %>%
  ungroup()

# write.csv(a4, file.path(RESULTS_DIR, "cluster_by_original_CellType_cb_annot_pct_res0.8.rm23.csv"), row.names = FALSE)


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

# saveRDS(obj, file.path(RESULTS_DIR, "integrated_seurat_final.res0.8.rm23.rds"))

# write.csv(obj@meta.data, file.path(RESULTS_DIR, "integrated_seurat_final.res0.8.rm23_metadata.csv"), row.names = F)

barcode_both <- data.frame(
  barcode = rownames(obj@meta.data),
  CellType_cb = obj$CellType_cb,
  celltype_annot = obj$celltype_annot,
  stringsAsFactors = FALSE
)

# write.csv(
#   barcode_both,
#   file.path(RESULTS_DIR, "barcode_original_annotation_res0.8.rm23.csv"),
#   row.names = FALSE
# )


# Standard UMAPs
# save_plot(DimPlot(obj, reduction = "umap_harmony", group.by = "study",     raster = FALSE),
#           "02_umap_by_study_res0.8.rm23")
# save_plot(DimPlot(obj, reduction = "umap_harmony", group.by = "batch",    raster = FALSE) + NoLegend(),
#           "03_umap_by_sample_res0.8.rm23")
# save_plot(DimPlot(obj, reduction = "umap_harmony", group.by = "sample_condition", raster = FALSE),
#           "04_umap_by_condition_res0.8.rm23")
# save_plot(DimPlot(obj, reduction = "umap_harmony", group.by = "cluster",
#                   label = TRUE, raster = FALSE) + NoLegend(),
#           "05_umap_by_cluster_res0.8.rm23")

# for (res in c(0.8)) {
#   col <- paste0("SCT_snn_res.", res)
#   save_plot(DimPlot(obj, reduction = "umap_harmony", group.by = "cluster",
#                     label = TRUE, raster = FALSE) + NoLegend() +
#               ggtitle(paste0("Clusters @ res ", res)),
#             paste0("05b_umap_clusters_res", res, ".rm23"))
# }
# save_plot(
#   DimPlot(obj, reduction = "umap_harmony", group.by = "celltype_annot",
#           label = TRUE, raster = FALSE) + NoLegend(),
#   "05c_umap_manual_celltype_res0.8.rm23"
# )





## clean cols and save final objs
new_cell_names <- colnames(obj)

new_cell_names <- gsub("_h5ad_", "_", new_cell_names)
new_cell_names <- sub("^h5ad_", "", new_cell_names)

new_cell_names <- make.unique(new_cell_names)

obj <- RenameCells(obj, new.names = new_cell_names)

obj$CellType_cb2 <- NULL
obj$CellType_cb3 <- NULL

# saveRDS(obj, file.path(RESULTS_DIR, "integrated_seurat_final1.res0.8.rm23.rds"))
# write.csv(obj@meta.data, file.path(RESULTS_DIR, "integrated_seurat_final1.res0.8.rm23_metadata.csv"), row.names = F)


barcode_both <- data.frame(
  barcode = rownames(obj@meta.data),
  CellType_cb = obj$CellType_cb,
  celltype_annot = obj$celltype_annot,
  stringsAsFactors = FALSE
)

# write.csv(
#   barcode_both,  file.path(RESULTS_DIR, "barcode_original_annotation1_res0.8.rm23.csv"),  row.names = FALSE)


message("Done.  Outputs in: ", RESULTS_DIR)

### write to count matrix
DefaultAssay(obj) <- "RNA"

#Layers(obj[["RNA"]])
obj <- JoinLayers(obj, assay = "RNA")

counts_mat <- GetAssayData(obj, assay = "RNA", layer = "counts")
# saveRDS(counts_mat, file.path(RESULTS_DIR, "counts_matrix_res0.8.rm23.rds"))


# Save MatrixMarket sparse matrix
# Matrix::writeMM( counts_mat, file.path(RESULTS_DIR, "counts_matrix_sparse_res0.8.rm23.mtx"))

# Save gene names
# write.table(
#   data.frame(gene = rownames(counts_mat)),
#   file = file.path(RESULTS_DIR, "counts_matrix_sparse_res0.8.rm23.gene.tsv"),
#   sep = "\t",
#   quote = FALSE,
#   row.names = FALSE,
#   col.names = FALSE
# )

# Save cell/barcode names
# write.table(
#   data.frame(cell = colnames(counts_mat)),
#   file = file.path(RESULTS_DIR, "counts_matrix_sparse_res0.8.rm23_barcodes.tsv"),
#   sep = "\t",
#   quote = FALSE,
#   row.names = FALSE,
#   col.names = FALSE
# )

#counts_df <- as.data.frame(as.matrix(counts_mat))
#counts_df <- cbind(gene = rownames(counts_mat), counts_df)

#write.table(
#  counts_df,
#  file = file.path(RESULTS_DIR, "counts_matrix_gse266HTAN_res0.8.tsv"),
#  sep = "\t",
#  quote = FALSE,
#  row.names = FALSE
#)


