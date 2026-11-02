#!/usr/bin/env Rscript

# 06_cell_type_marker_genes_rank.R
# ------------------------------------------------------------------------------
# Rank cluster marker genes from:
#   cluster_markers_res0.8.rm23_sub.csv
#
# Ranking rules:
#   1) smaller P value first
#   2) if tied, larger avg_log2FC first
#
# Two ranks are generated:
#   - rank_p:    p_val ascending, avg_log2FC descending
#   - rank_padj: p_val_adj ascending, avg_log2FC descending
#
# Known colorectal epithelial cell-type markers are then extracted for every
# cluster. Markers not present in the FindAllMarkers output are retained with
# status = "not_found".
#
# In addition, the reviewed cluster annotation from 05 is validated by
# checking whether the appropriate canonical markers occur among each
# cluster's Top 20 DE genes.
# ------------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
})

# ---- paths -------------------------------------------------------------------
RESULTS_DIR <- "/data/l2_bioinfo1/liq17/scTWAS/data/VUMC_HTAN_GSE266546/raw/rds2"

MARKER_FILE <- file.path(
  RESULTS_DIR,
  "cluster_markers_res0.8.rm23_sub.csv"
)

stopifnot(file.exists(MARKER_FILE))

# ---- known cell-type markers -------------------------------------------------
known_markers <- tibble::tribble(
  ~cell_type, ~gene,
  "ABS", "KRT20",
  "ABS", "GUCA2A",
  "ABS", "ALDOB",
  "GOB", "ATOH1",
  "GOB", "MUC2",
  "GOB", "TFF3",
  "STM", "LGR5",
  "STM", "OLFM4",
  "STM", "ASCL2",
  "EE",  "CHGA",
  "EE",  "NEUROD1",
  "CT",  "OTOP2",
  "CT",  "MEIS1",
  "TAC", "PCNA",
  "TAC", "MKI67",
  "TUF", "POU2F3",
  "TUF", "SOX9"
)

# ---- load marker table -------------------------------------------------------
markers <- read.csv(
  MARKER_FILE,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

required_cols <- c(
  "gene", "cluster", "avg_log2FC",
  "p_val", "p_val_adj", "pct.1", "pct.2"
)

missing_cols <- setdiff(required_cols, colnames(markers))

if (length(missing_cols) > 0) {
  stop(
    "Missing required columns in marker file: ",
    paste(missing_cols, collapse = ", ")
  )
}

markers <- markers %>%
  mutate(
    gene = toupper(as.character(gene)),
    cluster = as.character(cluster),
    avg_log2FC = as.numeric(avg_log2FC),
    p_val = as.numeric(p_val),
    p_val_adj = as.numeric(p_val_adj)
  )

message("Loaded marker table:")
message("  ", nrow(markers), " marker rows")
message("  ", length(unique(markers$cluster)), " clusters")

# ---- rank all genes within each cluster --------------------------------------
#
# rank_p:
#   primary = raw p_val ascending
#   secondary = avg_log2FC descending
#
# rank_padj:
#   primary = adjusted p_val ascending
#   secondary = avg_log2FC descending
#
# gene is used as a final deterministic tie-breaker.

ranked_markers <- markers %>%
  group_by(cluster) %>%
  arrange(
    p_val,
    desc(avg_log2FC),
    gene,
    .by_group = TRUE
  ) %>%
  mutate(
    rank_p = row_number()
  ) %>%
  arrange(
    p_val_adj,
    desc(avg_log2FC),
    gene,
    .by_group = TRUE
  ) %>%
  mutate(
    rank_padj = row_number()
  ) %>%
  ungroup()

# ---- save complete ranked marker table ---------------------------------------
full_rank_file <- file.path(
  RESULTS_DIR,
  "cluster_markers_res0.8.rm23_sub_ranked.csv"
)

write.csv(
  ranked_markers %>%
    arrange(
      as.numeric(cluster),
      rank_padj,
      rank_p
    ),
  full_rank_file,
  row.names = FALSE
)

message("Wrote: ", full_rank_file)

# ---- extract known markers for every cluster ---------------------------------
clusters <- sort(
  unique(ranked_markers$cluster),
  na.last = TRUE
)

expected_grid <- tidyr::crossing(
  cluster = clusters,
  known_markers
)

known_marker_ranks <- expected_grid %>%
  left_join(
    ranked_markers %>%
      select(
        cluster,
        gene,
        p_val,
        p_val_adj,
        avg_log2FC,
        pct.1,
        pct.2,
        rank_p,
        rank_padj
      ),
    by = c("cluster", "gene")
  ) %>%
  mutate(
    status = if_else(
      is.na(rank_p) & is.na(rank_padj),
      "not_found",
      "found"
    )
  ) %>%
  arrange(
    as.numeric(cluster),
    cell_type,
    rank_padj,
    rank_p
  )

known_rank_file <- file.path(
  RESULTS_DIR,
  "known_cell_type_marker_gene_ranks_res0.8.rm23_sub.csv"
)

write.csv(
  known_marker_ranks,
  known_rank_file,
  row.names = FALSE
)

message("Wrote: ", known_rank_file)

# ---- best cluster for each known marker --------------------------------------
#
# This table helps answer questions such as:
#   "In which cluster is OLFM4 ranked most strongly?"
#
best_cluster_per_marker <- known_marker_ranks %>%
  filter(status == "found") %>%
  group_by(cell_type, gene) %>%
  arrange(
    rank_padj,
    rank_p,
    p_val_adj,
    p_val,
    desc(avg_log2FC),
    .by_group = TRUE
  ) %>%
  mutate(best_cluster_rank = row_number()) %>%
  ungroup()

best_cluster_file <- file.path(
  RESULTS_DIR,
  "known_cell_type_markers_best_clusters_res0.8.rm23_sub.csv"
)

write.csv(
  best_cluster_per_marker,
  best_cluster_file,
  row.names = FALSE
)

message("Wrote: ", best_cluster_file)

# ---- concise top result printed to log ---------------------------------------
message("\nTop cluster hits for each known marker:")
print(
  best_cluster_per_marker %>%
    filter(best_cluster_rank <= 3) %>%
    select(
      cell_type,
      gene,
      cluster,
      best_cluster_rank,
      rank_padj,
      rank_p,
      avg_log2FC,
      p_val_adj,
      p_val,
      pct.1,
      pct.2
    ) %>%
    arrange(
      cell_type,
      gene,
      best_cluster_rank
    ),
  n = Inf
)

# ---- optional summary: top known markers within each cluster -----------------
top_known_by_cluster <- known_marker_ranks %>%
  filter(status == "found") %>%
  group_by(cluster) %>%
  arrange(
    rank_padj,
    rank_p,
    desc(avg_log2FC),
    .by_group = TRUE
  ) %>%
  mutate(known_marker_rank_within_cluster = row_number()) %>%
  ungroup()

top_known_file <- file.path(
  RESULTS_DIR,
  "known_markers_ranked_within_each_cluster_res0.8.rm23_sub.csv"
)

write.csv(
  top_known_by_cluster,
  top_known_file,
  row.names = FALSE
)

message("Wrote: ", top_known_file)


# ==============================================================================
# Top-20 DE validation of reviewed cluster annotations
# ==============================================================================
# Biological question:
#   For each cluster's assigned/reviewed cell type, are its canonical marker
#   genes among that cluster's Top 20 DE genes?
#
# Top-20 DE ranking:
#   1) p_val ascending (smallest P first)
#   2) avg_log2FC descending when P values are tied
#   3) gene name as a deterministic final tie-breaker
#
# Outputs:
#   - top20_DE_genes_per_cluster_res0.8.rm23_sub.csv
#   - cluster_annotation_marker_top20_check_res0.8.rm23_sub.csv
#   - cluster_annotation_marker_top20_summary_res0.8.rm23_sub.csv
# ------------------------------------------------------------------------------

ANNOTATION_FILE <- file.path(
  RESULTS_DIR,
  "cluster_reviewed_annotation_res0.8.rm23_sub.csv"
)

if (!file.exists(ANNOTATION_FILE)) {
  stop(
    "Reviewed annotation file not found: ",
    ANNOTATION_FILE,
    "\nRun 05_annotation_reviewed.R first."
  )
}

# ---- Top 20 DE genes per cluster ---------------------------------------------
top20_de <- markers %>%
  group_by(cluster) %>%
  arrange(
    p_val,
    desc(avg_log2FC),
    gene,
    .by_group = TRUE
  ) %>%
  mutate(DE_rank = row_number()) %>%
  filter(DE_rank <= 20) %>%
  ungroup() %>%
  arrange(as.numeric(cluster), DE_rank)

top20_file <- file.path(
  RESULTS_DIR,
  "top20_DE_genes_per_cluster_res0.8.rm23_sub.csv"
)

write.csv(
  top20_de,
  top20_file,
  row.names = FALSE
)

message("Wrote: ", top20_file)

# ---- Load reviewed annotation from 05 ----------------------------------------
cluster_annotation <- read.csv(
  ANNOTATION_FILE,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

required_annotation_cols <- c("cluster", "reviewed_annotation")
missing_annotation_cols <- setdiff(
  required_annotation_cols,
  colnames(cluster_annotation)
)

if (length(missing_annotation_cols) > 0) {
  stop(
    "Missing required columns in reviewed annotation file: ",
    paste(missing_annotation_cols, collapse = ", ")
  )
}

cluster_annotation <- cluster_annotation %>%
  transmute(
    cluster = as.character(cluster),
    annotation = as.character(reviewed_annotation)
  ) %>%
  distinct() %>%
  arrange(as.numeric(cluster))

# Safety check: one reviewed annotation per cluster.
annotation_check <- cluster_annotation %>%
  count(cluster, name = "n_annotations")

bad_annotation <- annotation_check %>%
  filter(n_annotations != 1)

if (nrow(bad_annotation) > 0) {
  stop(
    "Some clusters have more than one reviewed annotation: ",
    paste(bad_annotation$cluster, collapse = ", ")
  )
}

message("\nReviewed cluster annotations used for Top-20 validation:")
print(cluster_annotation, row.names = FALSE)

# ---- Expected canonical markers for each assigned cell type ------------------
#
# Reuse known_markers defined above:
#   cell_type, gene
#
# Clusters still annotated as Review (or another type without known markers)
# are not included in this canonical-marker validation table.
expected_markers <- cluster_annotation %>%
  inner_join(
    known_markers,
    by = c("annotation" = "cell_type")
  ) %>%
  rename(marker_gene = gene)

# ---- Check whether each canonical marker is in the cluster's Top 20 DE -------
marker_top20_check <- expected_markers %>%
  left_join(
    top20_de %>%
      select(
        cluster,
        gene,
        DE_rank,
        p_val,
        p_val_adj,
        avg_log2FC,
        pct.1,
        pct.2
      ),
    by = c(
      "cluster" = "cluster",
      "marker_gene" = "gene"
    )
  ) %>%
  mutate(
    in_top20 = if_else(is.na(DE_rank), "No", "Yes")
  ) %>%
  select(
    cluster,
    annotation,
    marker_gene,
    in_top20,
    DE_rank,
    p_val,
    p_val_adj,
    avg_log2FC,
    pct.1,
    pct.2
  ) %>%
  arrange(
    as.numeric(cluster),
    annotation,
    desc(in_top20),
    DE_rank,
    marker_gene
  )

marker_top20_file <- file.path(
  RESULTS_DIR,
  "cluster_annotation_marker_top20_check_res0.8.rm23_sub.csv"
)

write.csv(
  marker_top20_check,
  marker_top20_file,
  row.names = FALSE
)

message("Wrote: ", marker_top20_file)

# ---- Concise cluster-level summary -------------------------------------------
marker_top20_summary <- marker_top20_check %>%
  group_by(cluster, annotation) %>%
  summarise(
    n_expected_markers = n(),
    n_markers_in_top20 = sum(in_top20 == "Yes"),
    markers_in_top20 = {
      keep <- in_top20 == "Yes"

      if (!any(keep)) {
        "None"
      } else {
        paste0(
          marker_gene[keep],
          " (rank ",
          DE_rank[keep],
          ", P=",
          format(p_val[keep], scientific = TRUE, digits = 3),
          ", log2FC=",
          round(avg_log2FC[keep], 3),
          ")",
          collapse = "; "
        )
      }
    },
    .groups = "drop"
  ) %>%
  arrange(as.numeric(cluster))

marker_top20_summary_file <- file.path(
  RESULTS_DIR,
  "cluster_annotation_marker_top20_summary_res0.8.rm23_sub.csv"
)

write.csv(
  marker_top20_summary,
  marker_top20_summary_file,
  row.names = FALSE
)

message("Wrote: ", marker_top20_summary_file)

# ---- Print exact validation result to log ------------------------------------
message("\nCanonical markers in each cluster's Top 20 DE genes:")
print(
  marker_top20_check,
  row.names = FALSE
)

message("\nTop-20 canonical-marker summary by cluster:")
print(
  marker_top20_summary,
  row.names = FALSE
)


message("\nDone.")
