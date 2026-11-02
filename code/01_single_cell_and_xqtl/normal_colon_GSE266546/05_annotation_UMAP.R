#!/usr/bin/env Rscript
# 05_annotation.R
# ----------------------------------------------------------------------------
# Annotate clusters using the dominant ORIGINAL cell type in each cluster.
# Workflow:
#   1) Load the final clustered object from 04_sub_cluster.R
#   2) Harmonize original CellType_cb labels
#   3) Calculate original cell-type composition within each cluster
#   4) Identify the dominant cell type and its percentage in each cluster
#   5) Automatically annotate only clusters with sufficient dominance
#      (default >= 60%); otherwise label as "Review"
#   6) Save composition tables, annotated object, metadata, and UMAP
# ----------------------------------------------------------------------------

options(bitmapType = "cairo")
set.seed(123)

# ---- project / renv ---------------------------------------------------------
PROJECT <- "/data/l2_bioinfo1/liq17/scTWAS/data/VUMC_HTAN_GSE266546/src"
source(file.path(PROJECT, "renv/activate.R"))

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
})

# ---- paths ------------------------------------------------------------------
RESULTS_DIR <- "/data/l2_bioinfo1/liq17/scTWAS/data/VUMC_HTAN_GSE266546/raw/rds2"
FIG_DIR <- file.path(RESULTS_DIR, "figures")
dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)

rds_in <- file.path(
  RESULTS_DIR,
  "integrated_seurat_final1.res0.8.rm23_sub.rds"
)

stopifnot(file.exists(rds_in))
message("Loading: ", rds_in)
obj <- readRDS(rds_in)

# ---- input checks ------------------------------------------------------------
required_meta <- c("cluster", "CellType_cb")
missing_meta <- setdiff(required_meta, colnames(obj@meta.data))
if (length(missing_meta) > 0) {
  stop("Missing metadata columns: ", paste(missing_meta, collapse = ", "))
}

if (!"umap_harmony" %in% Reductions(obj)) {
  stop("Reduction 'umap_harmony' is missing from the object.")
}

message("Clusters from 04:")
print(table(obj$cluster, useNA = "ifany"))

# =============================================================================
# 1. Harmonize the ORIGINAL annotations
# =============================================================================
obj$CellType_cb2 <- as.character(obj$CellType_cb)

obj$CellType_cb2[
  grepl("Goblet", obj$CellType_cb2, ignore.case = TRUE)
] <- "GOB"

obj$CellType_cb2[
  grepl("Stem", obj$CellType_cb2, ignore.case = TRUE)
] <- "STM"

obj$CellType_cb2[
  grepl("Tuft", obj$CellType_cb2, ignore.case = TRUE)
] <- "TUF"

obj$CellType_cb2[
  grepl("EEC|Enteroendocrine", obj$CellType_cb2, ignore.case = TRUE)
] <- "EE"

obj$CellType_cb3 <- as.character(obj$CellType_cb2)

obj$CellType_cb3[
  grepl("Colonocyte|Enterocyte", obj$CellType_cb3, ignore.case = TRUE)
] <- "ABS"

obj$CellType_cb3[
  grepl("Cycling TA", obj$CellType_cb3, ignore.case = TRUE)
] <- "TAC"

# Save original-label cross-tabs for inspection
original_by_cluster <- as.data.frame(
  table(obj$CellType_cb, obj$cluster, useNA = "ifany")
) %>%
  rename(CellType_cb = Var1, cluster = Var2, n = Freq) %>%
  filter(n > 0) %>%
  group_by(cluster) %>%
  arrange(desc(n), .by_group = TRUE) %>%
  ungroup()

write.csv(
  original_by_cluster,
  file.path(
    RESULTS_DIR,
    "cluster_by_original_CellType_cb_res0.8.rm23_sub.csv"
  ),
  row.names = FALSE
)

# =============================================================================
# 2. Cell-type composition within each cluster
# =============================================================================
# Cluster total includes all cells.
cluster_totals <- obj@meta.data %>%
  transmute(cluster = as.character(cluster)) %>%
  count(cluster, name = "cluster_total")

# Composition is calculated from cells that have a usable original annotation.
cluster_comp <- obj@meta.data %>%
  transmute(
    cluster = as.character(cluster),
    celltype = as.character(CellType_cb3)
  ) %>%
  filter(!is.na(celltype), celltype != "") %>%
  count(cluster, celltype, name = "n_cells") %>%
  group_by(cluster) %>%
  mutate(
    annotated_total = sum(n_cells),
    percent_of_annotated = round(100 * n_cells / annotated_total, 2)
  ) %>%
  ungroup() %>%
  left_join(cluster_totals, by = "cluster") %>%
  mutate(
    percent_of_cluster = round(100 * n_cells / cluster_total, 2)
  ) %>%
  arrange(as.numeric(cluster), desc(percent_of_annotated))

write.csv(
  cluster_comp,
  file.path(
    RESULTS_DIR,
    "cluster_celltype_composition_res0.8.rm23_sub.csv"
  ),
  row.names = FALSE
)

message("\nCell-type composition by cluster:")
print(cluster_comp, n = Inf)

# =============================================================================
# 3. Dominant original cell type for each cluster
# =============================================================================
dominant_celltype <- cluster_comp %>%
  group_by(cluster) %>%
  slice_max(
    order_by = percent_of_annotated,
    n = 1,
    with_ties = FALSE
  ) %>%
  ungroup() %>%
  arrange(as.numeric(cluster)) %>%
  rename(
    dominant_celltype = celltype,
    dominant_n = n_cells,
    dominant_percent = percent_of_annotated
  )

# ---- confidence threshold ----------------------------------------------------
# A cluster is automatically assigned to its dominant original cell type only
# if that type accounts for at least this percentage of annotated cells.
# Change this value if you want a stricter or more permissive criterion.
DOMINANT_MIN_PCT <- 60

dominant_celltype <- dominant_celltype %>%
  mutate(
    annotation_status = ifelse(
      dominant_percent >= DOMINANT_MIN_PCT,
      "auto",
      "Review"
    ),
    proposed_annotation = ifelse(
      dominant_percent >= DOMINANT_MIN_PCT,
      dominant_celltype,
      "Review"
    )
  )

write.csv(
  dominant_celltype,
  file.path(
    RESULTS_DIR,
    "cluster_dominant_celltype_res0.8.rm23_sub.csv"
  ),
  row.names = FALSE
)

message("\nDominant cell type for each cluster:")
print(
  dominant_celltype %>%
    select(
      cluster,
      dominant_celltype,
      dominant_n,
      cluster_total,
      dominant_percent,
      annotation_status,
      proposed_annotation
    ),
  n = Inf
)

# =============================================================================
# 4. Assign cluster-level annotation
# =============================================================================
annotation_map <- setNames(
  dominant_celltype$proposed_annotation,
  dominant_celltype$cluster
)

obj$celltype_annot <- unname(
  annotation_map[as.character(obj$cluster)]
)

# Any cluster absent from the composition table should be reviewed manually.
obj$celltype_annot[is.na(obj$celltype_annot)] <- "Review"

message("\nCluster -> final annotation:")
print(table(obj$cluster, obj$celltype_annot, useNA = "ifany"))

# =============================================================================
# 5. Detailed comparison table
# =============================================================================
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

write.csv(
  a4,
  file.path(
    RESULTS_DIR,
    "cluster_by_original_CellType_cb_annot_pct_res0.8.rm23_sub.csv"
  ),
  row.names = FALSE
)

# =============================================================================
# 6. Save annotated object and metadata
# =============================================================================
message(
  "    ", nrow(obj), " genes x ",
  ncol(obj), " cells"
)

message(
  "    ", length(unique(obj$SampleID)), " samples, ",
  length(unique(obj$PatientID)), " patients"
)

message(
  "    ", length(unique(obj$cluster)), " clusters; ",
  length(unique(obj$celltype_annot)), " annotation labels"
)

saveRDS(
  obj,
  file.path(
    RESULTS_DIR,
    "integrated_seurat_final2.res0.8.rm23_sub.rds"
  )
)

# write.csv(
  # obj@meta.data,
  # file.path(
    # RESULTS_DIR,
    # "integrated_seurat_final2.res0.8.rm23_sub_metadata.csv"
  # ),
  # row.names = FALSE
# )

# barcode_both <- data.frame(
  # barcode = rownames(obj@meta.data),
  # CellType_cb = obj$CellType_cb,
  # CellType_cb3 = obj$CellType_cb3,
  # cluster = as.character(obj$cluster),
  # celltype_annot = obj$celltype_annot,
  # stringsAsFactors = FALSE
# )

# write.csv(
  # barcode_both,
  # file.path(
    # RESULTS_DIR,
    # "barcode_original_annotation2_res0.8.rm23_sub.csv"
  # ),
  # row.names = FALSE
# )

# =============================================================================
# 7. UMAP plots
# =============================================================================
save_plot <- function(p, name, w = 8, h = 6) {
  ggsave(
    filename = file.path(FIG_DIR, paste0(name, ".pdf")),
    plot = p,
    width = w,
    height = h,
    device = cairo_pdf
  )
}

save_plot(
  DimPlot(
    obj,
    reduction = "umap_harmony",
    group.by = "celltype_annot",
    label = FALSE,
    raster = FALSE
  ) +
    labs(colour = "Cell type") +
    guides(colour = guide_legend(override.aes = list(size = 3))) +
    theme(
      legend.position = "right",
      legend.text = element_text(size = 10),
      legend.title = element_text(size = 11)
    ),
  "05c_umap_dominant_celltype_res0.8.rm23_sub"
)

save_plot(
  DimPlot(
    obj,
    reduction = "umap_harmony",
    group.by = "cluster",
    label = TRUE,
    raster = FALSE
  ) + NoLegend(),
  "05d_umap_cluster_for_annotation_review_res0.8.rm23_sub"
)

# =============================================================================
# 8. Reviewed annotation: force selected Review clusters to dominant cell type
# =============================================================================
# After inspecting cluster composition / markers, clusters 8 and 10 are accepted
# as their dominant original cell type even if dominant_percent < DOMINANT_MIN_PCT.
REVIEW_USE_DOMINANT <- c("8", "10")

dominant_map_raw <- setNames(
  dominant_celltype$dominant_celltype,
  dominant_celltype$cluster
)

missing_review_clusters <- setdiff(
  REVIEW_USE_DOMINANT,
  names(dominant_map_raw)
)
if (length(missing_review_clusters) > 0) {
  stop(
    "No dominant cell type available for reviewed cluster(s): ",
    paste(missing_review_clusters, collapse = ", ")
  )
}

# Keep the automatic annotation intact and make a separate reviewed column.
obj$celltype_annot_reviewed <- as.character(obj$celltype_annot)

idx_review <- as.character(obj$cluster) %in% REVIEW_USE_DOMINANT
obj$celltype_annot_reviewed[idx_review] <- unname(
  dominant_map_raw[as.character(obj$cluster[idx_review])]
)

message("\nReviewed cluster -> final annotation:")
print(
  table(
    obj$cluster,
    obj$celltype_annot_reviewed,
    useNA = "ifany"
  )
)

reviewed_annotation_summary <- dominant_celltype %>%
  transmute(
    cluster = as.character(cluster),
    dominant_celltype,
    dominant_percent,
    automatic_annotation = proposed_annotation,
    manual_override = cluster %in% REVIEW_USE_DOMINANT,
    reviewed_annotation = ifelse(
      cluster %in% REVIEW_USE_DOMINANT,
      dominant_celltype,
      proposed_annotation
    )
  ) %>%
  arrange(as.numeric(cluster))

write.csv(
  reviewed_annotation_summary,
  file.path(
    RESULTS_DIR,
    "cluster_reviewed_annotation_res0.8.rm23_sub.csv"
  ),
  row.names = FALSE
)

# 05e: reviewed annotation UMAP. Clusters 8 and 10 now use dominant cell type.
save_plot(
  DimPlot(
    obj,
    reduction = "umap_harmony",
    group.by = "celltype_annot_reviewed",
    label = FALSE,
    raster = FALSE
  ) +
    labs(colour = "Cell type") +
    guides(colour = guide_legend(override.aes = list(size = 3))) +
    theme(
      legend.position = "right",
      legend.text = element_text(size = 10),
      legend.title = element_text(size = 11)
    ),
  "05e_umap_reviewed_celltype_res0.8.rm23_sub"
)

# Save an additional reviewed object so the automatic final2 object is preserved.
saveRDS(
  obj,
  file.path(
    RESULTS_DIR,
    "integrated_seurat_final2_reviewed.res0.8.rm23_sub.rds"
  )
)

write.csv(
  obj@meta.data,
  file.path(
    RESULTS_DIR,
    "integrated_seurat_final2_reviewed.res0.8.rm23_sub_metadata.csv"
  ),
  row.names = FALSE
)

message("Done. Annotation outputs written to: ", RESULTS_DIR)
