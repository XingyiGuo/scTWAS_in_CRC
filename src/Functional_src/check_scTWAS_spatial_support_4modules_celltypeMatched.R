#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
})

# ==============================================================================
# Check scTWAS genes against four Visium spatial-expression result tables
#
# IMPORTANT:
#   * scTWAS files store ENSG IDs in column "gene".
#   * Spatial files store HGNC gene symbols in column "gene".
#   * Therefore matching is performed using scTWAS gene_name_* -> spatial gene.
#
# Four independent modules:
#   M1: de_senescence_low_grade_vs_normal.csv
#   M2: de_advanced_vs_nonadvanced.csv
#   M3: de_stemness_high_vs_low.csv
#   M4: cluster_markers.csv
#       - additionally tests whether a scTWAS gene is a marker of a
#         biologically matching Visium cluster for its scTWAS cell type.
#
# For each module the script writes:
#   1) *_all_scTWAS_annotated.csv
#      - all ADSER + NL scTWAS rows
#      - spatial_present = Yes/No
#      - spatial_FDR05_any = Yes/No
#   2) *_overlap_detail.csv
#      - only matched scTWAS rows, joined to the spatial result
#      - cluster_markers can have >1 row per gene because a gene may mark >1 cluster
#   3) *_summary_by_dataset.csv
#   4) *_summary_by_celltype.csv
#
# It also writes:
#   all_4modules_summary_by_dataset.csv
#   all_4modules_summary_by_celltype.csv
# ==============================================================================

# ------------------------------------------------------------------------------
# 0. USER CONFIGURATION
# ------------------------------------------------------------------------------

ADSER_FILE <- "/data/l2_bioinfo1/liq17/scTWAS/TWAS_HTAN_GSE266546_WellcomeSanger_res/SPrediXcan_peakSNP_association/VUMC_HTAN_ADSER_EUR_Multi_sig/VUMC_HTAN_ADSER_9CTs_peakSNP.EUR_Multi.sig.csv"

NL_FILE <- "/data/l2_bioinfo1/liq17/scTWAS/TWAS_HTAN_GSE266546_WellcomeSanger_res/SPrediXcan_peakSNP_association/Cauchy_Meta_NL_WellcomeSanger_NL_EUR_Multi/HTAN_NL_WellcomeSanger_NL_7celltypes_peakSNP.EUR_Multi.CauchyMeta.sig.csv"

# Directory from:
# https://github.com/grady-lab/Spatial-multi-omics/tree/main/Visium/results/tables
SPATIAL_DIR <- "/data/l2_bioinfo1/liq17/scTWAS/data/SpatialExpr/CRC_SpatilTrans_PMID42620287/Spatial-multi-omics/Visium/results/tables"

OUT_DIR <- "/data/l2_bioinfo1/liq17/scTWAS/Functional_res/scTWAS_spatial_support_4modules"

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)


# ------------------------------------------------------------------------------
# 0B. scTWAS cell type -> matching Visium cluster map for MODULE 4
# ------------------------------------------------------------------------------
#
# scTWAS cell types:
#   ABS = absorptive cells
#   CT  = crypt-top colonocytes
#   EE  = enteroendocrine cells
#   GOB = goblet cells
#   STM = stem cells
#   TAC = transit-amplifying cells
#   TUF = tuft cells
#   ASC = adenoma-specific cells
#   SSC = serrated-specific cells
#
# IMPORTANT:
#   The Visium clusters are SPOT-LEVEL spatial states, not pure single-cell
#   populations. Therefore this is a biological correspondence map rather than
#   a one-to-one cell-type annotation.
#
# Conservative mapping used here:
#
#   ABS -> C3
#     C3 has the strongest absorptive-colonocyte program in the manuscript.
#
#   CT -> C3, C6
#     C3/C6 are differentiated normal epithelium; CT correspondence is based on
#     this manuscript classification plus the marker profiles.
#
#   GOB -> C0
#     Marker-informed: C0 is enriched for RETNLB, CLCA1, TFF3, BEST2, ATOH1.
#     The manuscript does not explicitly name C0 "goblet".
#
#   STM -> C1, C8 in NL
#     C1 is manuscript-described crypt-base/stem-like epithelium.
#     C8 is marker-informed from LGR5/OLFM4/SMOC2 enrichment.
#
#   STM -> C1, C2, C8 in ADSER
#     C2 is added because it carries a dysplastic stemness + proliferation
#     program in adenoma.
#
#   EE -> C6
#     Marker-informed: C6 contains PYY, CHGA, GCG, CHGB.
#
#   TAC -> C2
#     C2 is strongly proliferative (MKI67/TOP2A/PCNA).
#
#   TUF -> not assessable
#     No dedicated tuft-like Visium cluster is defined in this manuscript.
#
#   ASC -> C2, C4
#     C2/C4 are the manuscript's principal dysplastic epithelial clusters.
#
#   SSC -> not assessable
#     The Visium cohort contains conventional adenomas, not serrated lesions,
#     so there is no direct serrated-specific matching cluster.
#
# NOTE:
#   This mapping is intentionally explicit and editable. If you later obtain
#   cell2location-based cluster-by-cell-type enrichment from the study, replace
#   this correspondence map with the empirical deconvolution mapping.
# ------------------------------------------------------------------------------

CELLTYPE_CLUSTER_MAP <- data.table(
  scTWAS_dataset = c(
    rep("NL", 7),
    rep("ADSER", 9)
  ),
  CellType = c(
    "ABS", "GOB", "CT", "STM", "EE", "TUF", "TAC",
    "ABS", "GOB", "CT", "STM", "EE", "TUF", "TAC", "ASC", "SSC"
  ),
  matching_clusters = c(
    # NL
    "3",       # ABS
    "0",       # GOB
    "3,6",     # CT
    "1,8",     # STM
    "6",       # EE
    NA,        # TUF: no dedicated Visium cluster
    "2",       # TAC

    # ADSER
    "3",       # ABS
    "0",       # GOB
    "3,6",     # CT
    "1,2,8",   # STM
    "6",       # EE
    NA,        # TUF: no dedicated Visium cluster
    "2",       # TAC
    "2,4",     # ASC: principal dysplastic clusters
    NA         # SSC: no serrated lesions in Visium cohort
  ),
  mapping_basis = c(
    # NL
    "Manuscript: C3 has strongest absorptive-colonocyte markers",
    "Marker-informed: C0 RETNLB/CLCA1/TFF3/BEST2/ATOH1",
    "Manuscript: C3/C6 differentiated normal epithelium; CT correspondence inferred",
    "C1 manuscript crypt-base/stem-like; C8 marker-informed LGR5/OLFM4/SMOC2",
    "Marker-informed: C6 PYY/CHGA/GCG/CHGB",
    "Not assessable: no dedicated tuft-like Visium cluster",
    "Manuscript: C2 proliferative MKI67/TOP2A/PCNA",

    # ADSER
    "Manuscript: C3 has strongest absorptive-colonocyte markers",
    "Marker-informed: C0 RETNLB/CLCA1/TFF3/BEST2/ATOH1",
    "Manuscript: C3/C6 differentiated normal epithelium; CT correspondence inferred",
    "C1/C8 stem-like plus C2 dysplastic stemness/proliferation",
    "Marker-informed: C6 PYY/CHGA/GCG/CHGB",
    "Not assessable: no dedicated tuft-like Visium cluster",
    "Manuscript: C2 proliferative MKI67/TOP2A/PCNA",
    "Manuscript: C2/C4 are principal dysplastic epithelial clusters",
    "Not assessable: Visium cohort does not contain serrated lesions"
  ),
  mapping_confidence = c(
    "High",
    "Moderate",
    "Moderate",
    "Moderate",
    "Moderate",
    "Not_assessable",
    "High",

    "High",
    "Moderate",
    "Moderate",
    "Moderate_to_high",
    "Moderate",
    "Not_assessable",
    "High",
    "High",
    "Not_assessable"
  )
)


# ------------------------------------------------------------------------------
# 1. Helper functions
# ------------------------------------------------------------------------------

clean_char <- function(x) {
  x <- trimws(as.character(x))
  x[x %in% c("", "NA", "NaN", "NULL", ".")] <- NA_character_
  x
}


first_nonmissing_column <- function(dt, candidates) {

  candidates <- candidates[candidates %in% names(dt)]

  if (length(candidates) == 0L) {
    stop(
      "None of the expected gene-symbol columns were found. Expected one of:\n",
      paste(candidates, collapse = ", ")
    )
  }

  xlist <- lapply(candidates, function(cc) clean_char(dt[[cc]]))

  out <- xlist[[1]]

  if (length(xlist) > 1L) {
    for (ii in 2:length(xlist)) {
      idx <- is.na(out) & !is.na(xlist[[ii]])
      out[idx] <- xlist[[ii]][idx]
    }
  }

  out
}


check_symbol_consistency <- function(dt, candidates, label) {

  candidates <- candidates[candidates %in% names(dt)]
  if (length(candidates) < 2L) return(invisible(NULL))

  mat <- lapply(candidates, function(cc) toupper(clean_char(dt[[cc]])))

  conflict <- vapply(
    seq_len(nrow(dt)),
    function(i) {
      vals <- unique(na.omit(vapply(mat, `[`, character(1), i)))
      length(vals) > 1L
    },
    logical(1)
  )

  n_conflict <- sum(conflict)

  if (n_conflict > 0L) {
    warning(
      label, ": ", n_conflict,
      " rows have inconsistent gene symbols across gene_name_* columns. ",
      "The first non-missing symbol in the configured priority order will be used."
    )
  }
}


prepare_sctwas <- function(file, dataset) {

  message("Reading scTWAS: ", file)

  dt <- fread(file, data.table = TRUE, check.names = FALSE)

  required <- c("gene", "CellType")
  missing_required <- setdiff(required, names(dt))

  if (length(missing_required) > 0L) {
    stop(
      dataset, " scTWAS file is missing required columns: ",
      paste(missing_required, collapse = ", ")
    )
  }

  # Preserve ENSG identifier separately.
  setnames(dt, "gene", "scTWAS_ENSG")

  if (dataset == "ADSER") {

    symbol_candidates <- c(
      "gene_name_EUR",
      "gene_name_Multi"
    )

  } else if (dataset == "NL") {

    # Priority:
    # HTAN NL EUR -> Wellcome Sanger EUR -> HTAN NL Multi -> Wellcome Sanger Multi
    symbol_candidates <- c(
      "gene_name_HTAN_NL_EUR",
      "gene_name_WellcomeSanger_NL_EUR",
      "gene_name_HTAN_NL_Multi",
      "gene_name_WellcomeSanger_NL_Multi"
    )

  } else {
    stop("Unknown scTWAS dataset: ", dataset)
  }

  check_symbol_consistency(dt, symbol_candidates, dataset)

  # Create one complete overall HGNC-symbol column across gene_name_* fields.
  dt[, gene_name_overall := first_nonmissing_column(.SD, symbol_candidates),
     .SDcols = symbol_candidates[symbol_candidates %in% names(dt)]]

  dt[, gene_name_overall := toupper(clean_char(gene_name_overall))]

  # Keep gene_symbol as the standardized join key used by downstream modules.
  dt[, gene_symbol := gene_name_overall]

  dt[, CellType := toupper(clean_char(CellType))]
  dt[, scTWAS_dataset := dataset]

  if (anyNA(dt$gene_symbol)) {
    warning(
      dataset, ": ", sum(is.na(dt$gene_symbol)),
      " rows have no usable HGNC symbol and cannot be matched to spatial tables."
    )
  }

  # Stable key for summary counting. One scTWAS result = dataset + cell type + ENSG.
  dt[, scTWAS_key := paste(scTWAS_dataset, CellType, scTWAS_ENSG, sep = "|")]

  dt[]
}


prepare_spatial <- function(file) {

  message("Reading spatial table: ", file)

  sp <- fread(file, data.table = TRUE, check.names = FALSE)

  if (!"gene" %in% names(sp)) {
    stop("Spatial file has no 'gene' column: ", file)
  }

  # Avoid confusing spatial gene symbol with scTWAS ENSG.
  setnames(sp, "gene", "spatial_gene")

  sp[, spatial_gene := clean_char(spatial_gene)]
  sp[, spatial_gene_key := toupper(spatial_gene)]

  sp[]
}


summarize_annotation <- function(annot, by_cols) {

  annot[
    ,
    .(
      n_scTWAS_rows = .N,
      n_scTWAS_unique_genes = uniqueN(gene_symbol[!is.na(gene_symbol)]),

      n_present_rows = sum(spatial_present == "Yes", na.rm = TRUE),
      n_present_unique_genes = uniqueN(
        gene_symbol[spatial_present == "Yes" & !is.na(gene_symbol)]
      ),

      n_FDR05_rows = sum(spatial_FDR05_any == "Yes", na.rm = TRUE),
      n_FDR05_unique_genes = uniqueN(
        gene_symbol[spatial_FDR05_any == "Yes" & !is.na(gene_symbol)]
      )
    ),
    by = by_cols
  ]
}


run_module <- function(
    module_id,
    module_name,
    comparison,
    spatial_file,
    sctwas_all,
    out_dir
) {

  message("")
  message("======================================================================")
  message(module_id, ": ", module_name)
  message("Comparison: ", comparison)
  message("======================================================================")

  sp <- prepare_spatial(spatial_file)

  # Gene-level presence, irrespective of p-value.
  present_genes <- unique(
    sp[!is.na(spatial_gene_key), spatial_gene_key]
  )

  # Gene-level FDR<0.05 support.
  if ("p_val_adj" %in% names(sp)) {
    fdr05_genes <- unique(
      sp[
        !is.na(spatial_gene_key) &
          !is.na(p_val_adj) &
          p_val_adj < 0.05,
        spatial_gene_key
      ]
    )
  } else {
    fdr05_genes <- character(0)
    warning(module_id, ": p_val_adj is absent; spatial_FDR05_any will be 'NA'.")
  }

  # ---------------------------------------------------------------------------
  # A. Annotate every scTWAS row with spatial presence.
  # ---------------------------------------------------------------------------

  annot <- copy(sctwas_all)

  annot[, spatial_module := module_id]
  annot[, spatial_analysis := module_name]
  annot[, spatial_comparison := comparison]

  annot[
    ,
    spatial_present := fifelse(
      !is.na(gene_symbol) & gene_symbol %chin% present_genes,
      "Yes",
      "No"
    )
  ]

  if ("p_val_adj" %in% names(sp)) {
    annot[
      ,
      spatial_FDR05_any := fifelse(
        !is.na(gene_symbol) & gene_symbol %chin% fdr05_genes,
        "Yes",
        "No"
      )
    ]
  } else {
    annot[, spatial_FDR05_any := NA_character_]
  }

  # ---------------------------------------------------------------------------
  # B. Detailed overlap.
  #
  # allow.cartesian=TRUE is intentional:
  # cluster_markers.csv can contain the same gene in several spatial clusters.
  # ---------------------------------------------------------------------------

  overlap <- merge(
    annot[spatial_present == "Yes"],
    sp,
    by.x = "gene_symbol",
    by.y = "spatial_gene_key",
    all = FALSE,
    allow.cartesian = TRUE,
    sort = FALSE
  )

  # Bring important identifiers to the front.
  front_cols <- intersect(
    c(
      "spatial_module",
      "spatial_analysis",
      "spatial_comparison",
      "scTWAS_dataset",
      "CellType",
      "scTWAS_ENSG",
      "gene_symbol",
      "spatial_gene",
      "spatial_present",
      "spatial_FDR05_any",
      "p_val",
      "avg_log2FC",
      "pct.1",
      "pct.2",
      "p_val_adj",
      "cluster",
      "delta_pct"
    ),
    names(overlap)
  )

  setcolorder(overlap, c(front_cols, setdiff(names(overlap), front_cols)))

  # ---------------------------------------------------------------------------
  # C. Summaries.
  # Counts are based on the annotated scTWAS rows, NOT the cartesian overlap,
  # so genes appearing in multiple clusters do not inflate the summary.
  # ---------------------------------------------------------------------------

  summary_dataset <- summarize_annotation(
    annot,
    by_cols = c(
      "spatial_module",
      "spatial_analysis",
      "spatial_comparison",
      "scTWAS_dataset"
    )
  )

  summary_celltype <- summarize_annotation(
    annot,
    by_cols = c(
      "spatial_module",
      "spatial_analysis",
      "spatial_comparison",
      "scTWAS_dataset",
      "CellType"
    )
  )

  prefix <- file.path(
    out_dir,
    paste0(module_id, "_", module_name)
  )

  fwrite(
    annot,
    paste0(prefix, "_all_scTWAS_annotated.csv"),
    na = ""
  )

  fwrite(
    overlap,
    paste0(prefix, "_overlap_detail.csv"),
    na = ""
  )

  fwrite(
    summary_dataset,
    paste0(prefix, "_summary_by_dataset.csv"),
    na = ""
  )

  fwrite(
    summary_celltype,
    paste0(prefix, "_summary_by_celltype.csv"),
    na = ""
  )

  message("")
  message("Summary by scTWAS dataset:")
  print(summary_dataset)

  invisible(
    list(
      annotation = annot,
      overlap = overlap,
      summary_dataset = summary_dataset,
      summary_celltype = summary_celltype
    )
  )
}



# ==============================================================================
# Specialized MODULE 4:
# Does each cell-type-specific scTWAS gene mark its BIOLOGICALLY MATCHING
# Visium cluster(s)?
# ==============================================================================

run_cluster_marker_module <- function(
    module_id,
    module_name,
    comparison,
    spatial_file,
    sctwas_all,
    out_dir,
    celltype_cluster_map
) {

  message("")
  message("======================================================================")
  message(module_id, ": ", module_name, " -- CELL-TYPE MATCHED CLUSTERS")
  message("======================================================================")

  sp <- prepare_spatial(spatial_file)

  if (!"cluster" %in% names(sp)) {
    stop("cluster_markers.csv must contain a 'cluster' column.")
  }

  sp[, cluster_chr := as.character(cluster)]

  # ---------------------------------------------------------------------------
  # A. Add expected matching cluster(s) to every scTWAS row
  # ---------------------------------------------------------------------------

  annot <- merge(
    copy(sctwas_all),
    celltype_cluster_map,
    by = c("scTWAS_dataset", "CellType"),
    all.x = TRUE,
    sort = FALSE
  )

  # Cell type missing from map is a configuration error.
  missing_map <- annot[
    is.na(mapping_basis),
    unique(paste(scTWAS_dataset, CellType, sep = ":"))
  ]

  if (length(missing_map) > 0L) {
    stop(
      "Missing cell-type -> Visium cluster mapping for: ",
      paste(missing_map, collapse = ", ")
    )
  }

  annot[, matching_cluster_assessable :=
          !is.na(matching_clusters) & nzchar(matching_clusters)]

  # ---------------------------------------------------------------------------
  # B. Build marker-cluster lookup for each spatial gene
  # ---------------------------------------------------------------------------

  gene_any_clusters <- sp[
    !is.na(spatial_gene_key),
    .(
      marker_clusters_any = paste(
        sort(unique(as.integer(cluster_chr))),
        collapse = ","
      )
    ),
    by = spatial_gene_key
  ]

  if ("p_val_adj" %in% names(sp)) {

    gene_fdr_clusters <- sp[
      !is.na(spatial_gene_key) &
        !is.na(p_val_adj) &
        p_val_adj < 0.05,
      .(
        marker_clusters_FDR05 = paste(
          sort(unique(as.integer(cluster_chr))),
          collapse = ","
        )
      ),
      by = spatial_gene_key
    ]

  } else {

    gene_fdr_clusters <- data.table(
      spatial_gene_key = character(),
      marker_clusters_FDR05 = character()
    )
  }

  setnames(gene_any_clusters, "spatial_gene_key", "gene_symbol")
  setnames(gene_fdr_clusters, "spatial_gene_key", "gene_symbol")

  annot <- merge(
    annot,
    gene_any_clusters,
    by = "gene_symbol",
    all.x = TRUE,
    sort = FALSE
  )

  annot <- merge(
    annot,
    gene_fdr_clusters,
    by = "gene_symbol",
    all.x = TRUE,
    sort = FALSE
  )

  # ---------------------------------------------------------------------------
  # C. Helper to intersect "expected cluster string" with observed clusters
  # ---------------------------------------------------------------------------

  intersect_cluster_strings <- function(expected, observed) {

    if (is.na(expected) || !nzchar(expected)) {
      return(NA_character_)
    }

    if (is.na(observed) || !nzchar(observed)) {
      return("")
    }

    exp_vec <- trimws(strsplit(expected, ",", fixed = TRUE)[[1]])
    obs_vec <- trimws(strsplit(observed, ",", fixed = TRUE)[[1]])

    hit <- intersect(exp_vec, obs_vec)

    if (length(hit) == 0L) {
      return("")
    }

    paste(sort(as.integer(hit)), collapse = ",")
  }


  annot[, matched_marker_clusters :=
          mapply(
            intersect_cluster_strings,
            matching_clusters,
            marker_clusters_any,
            USE.NAMES = FALSE
          )]

  annot[, matched_marker_clusters_FDR05 :=
          mapply(
            intersect_cluster_strings,
            matching_clusters,
            marker_clusters_FDR05,
            USE.NAMES = FALSE
          )]


  # ---------------------------------------------------------------------------
  # D. Cell-type-matched marker flags
  #
  # Not_assessable is used instead of "No" for TUF and SSC because the study
  # does not define a direct matching Visium cluster for those cell types.
  # ---------------------------------------------------------------------------

  annot[
    ,
    celltype_matched_cluster_marker := fifelse(
      !matching_cluster_assessable,
      "Not_assessable",
      fifelse(
        !is.na(matched_marker_clusters) & nzchar(matched_marker_clusters),
        "Yes",
        "No"
      )
    )
  ]

  if ("p_val_adj" %in% names(sp)) {

    annot[
      ,
      celltype_matched_cluster_marker_FDR05 := fifelse(
        !matching_cluster_assessable,
        "Not_assessable",
        fifelse(
          !is.na(matched_marker_clusters_FDR05) &
            nzchar(matched_marker_clusters_FDR05),
          "Yes",
          "No"
        )
      )
    ]

  } else {

    annot[, celltype_matched_cluster_marker_FDR05 := NA_character_]
  }


  # Keep the original broad "gene appears in ANY spatial cluster" information
  # as a secondary comparison.
  annot[
    ,
    spatial_present_any_cluster := fifelse(
      !is.na(marker_clusters_any) & nzchar(marker_clusters_any),
      "Yes",
      "No"
    )
  ]

  annot[
    ,
    spatial_FDR05_any_cluster := fifelse(
      !is.na(marker_clusters_FDR05) & nzchar(marker_clusters_FDR05),
      "Yes",
      "No"
    )
  ]

  annot[, spatial_module := module_id]
  annot[, spatial_analysis := module_name]
  annot[, spatial_comparison := comparison]


  # ---------------------------------------------------------------------------
  # E. Detailed rows ONLY from the biologically matching clusters
  # ---------------------------------------------------------------------------

  assessable <- annot[
    matching_cluster_assessable == TRUE
  ]

  detail <- merge(
    assessable,
    sp,
    by.x = "gene_symbol",
    by.y = "spatial_gene_key",
    all = FALSE,
    allow.cartesian = TRUE,
    sort = FALSE
  )

  # Retain spatial marker rows only when the observed cluster belongs to that
  # scTWAS row's expected/matching cluster set.
  detail[
    ,
    keep_matching_cluster := mapply(
      function(expected, observed_cluster) {
        expected_vec <- trimws(
          strsplit(expected, ",", fixed = TRUE)[[1]]
        )
        as.character(observed_cluster) %in% expected_vec
      },
      matching_clusters,
      cluster_chr,
      USE.NAMES = FALSE
    )
  ]

  detail <- detail[keep_matching_cluster == TRUE]

  detail[, keep_matching_cluster := NULL]

  front_cols <- intersect(
    c(
      "spatial_module",
      "spatial_analysis",
      "scTWAS_dataset",
      "CellType",
      "scTWAS_ENSG",
      "gene_symbol",
      "matching_clusters",
      "mapping_confidence",
      "mapping_basis",
      "celltype_matched_cluster_marker",
      "celltype_matched_cluster_marker_FDR05",
      "matched_marker_clusters",
      "matched_marker_clusters_FDR05",
      "cluster",
      "spatial_gene",
      "p_val",
      "avg_log2FC",
      "pct.1",
      "pct.2",
      "p_val_adj",
      "delta_pct"
    ),
    names(detail)
  )

  setcolorder(detail, c(front_cols, setdiff(names(detail), front_cols)))


  # ---------------------------------------------------------------------------
  # F. Summary
  #
  # Counts are computed from scTWAS rows, NOT from detail, so genes that are
  # markers of more than one matching cluster do not inflate counts.
  # ---------------------------------------------------------------------------

  summarize_matched_clusters <- function(x, by_cols) {

    x[
      ,
      .(
        n_scTWAS_rows = .N,
        n_scTWAS_unique_genes =
          uniqueN(gene_symbol[!is.na(gene_symbol)]),

        n_assessable_rows =
          sum(matching_cluster_assessable, na.rm = TRUE),

        n_assessable_unique_genes =
          uniqueN(
            gene_symbol[
              matching_cluster_assessable == TRUE &
                !is.na(gene_symbol)
            ]
          ),

        n_matched_cluster_marker_rows =
          sum(
            celltype_matched_cluster_marker == "Yes",
            na.rm = TRUE
          ),

        n_matched_cluster_marker_unique_genes =
          uniqueN(
            gene_symbol[
              celltype_matched_cluster_marker == "Yes" &
                !is.na(gene_symbol)
            ]
          ),

        n_matched_cluster_marker_FDR05_rows =
          sum(
            celltype_matched_cluster_marker_FDR05 == "Yes",
            na.rm = TRUE
          ),

        n_matched_cluster_marker_FDR05_unique_genes =
          uniqueN(
            gene_symbol[
              celltype_matched_cluster_marker_FDR05 == "Yes" &
                !is.na(gene_symbol)
            ]
          ),

        n_not_assessable_rows =
          sum(!matching_cluster_assessable, na.rm = TRUE),

        # Secondary broad comparison:
        n_marker_any_cluster_unique_genes =
          uniqueN(
            gene_symbol[
              spatial_present_any_cluster == "Yes" &
                !is.na(gene_symbol)
            ]
          )
      ),
      by = by_cols
    ]
  }


  summary_dataset <- summarize_matched_clusters(
    annot,
    by_cols = c(
      "spatial_module",
      "spatial_analysis",
      "scTWAS_dataset"
    )
  )

  summary_celltype <- summarize_matched_clusters(
    annot,
    by_cols = c(
      "spatial_module",
      "spatial_analysis",
      "scTWAS_dataset",
      "CellType",
      "matching_clusters",
      "mapping_confidence",
      "mapping_basis"
    )
  )


  # ---------------------------------------------------------------------------
  # G. Save
  # ---------------------------------------------------------------------------

  prefix <- file.path(
    out_dir,
    paste0(module_id, "_", module_name)
  )

  fwrite(
    annot,
    paste0(prefix, "_all_scTWAS_celltype_matched_annotated.csv"),
    na = ""
  )

  fwrite(
    detail,
    paste0(prefix, "_celltype_matched_cluster_marker_detail.csv"),
    na = ""
  )

  fwrite(
    summary_dataset,
    paste0(prefix, "_celltype_matched_summary_by_dataset.csv"),
    na = ""
  )

  fwrite(
    summary_celltype,
    paste0(prefix, "_celltype_matched_summary_by_celltype.csv"),
    na = ""
  )

  fwrite(
    celltype_cluster_map,
    paste0(prefix, "_celltype_to_cluster_mapping_used.csv"),
    na = ""
  )

  message("")
  message("Cell-type-matched cluster-marker summary by scTWAS dataset:")
  print(summary_dataset)

  invisible(
    list(
      annotation = annot,
      overlap = detail,
      summary_dataset = summary_dataset,
      summary_celltype = summary_celltype
    )
  )
}


# ------------------------------------------------------------------------------
# 2. Read and harmonize the two scTWAS result sets once
# ------------------------------------------------------------------------------

ADSER <- prepare_sctwas(
  ADSER_FILE,
  dataset = "ADSER"
)

NL <- prepare_sctwas(
  NL_FILE,
  dataset = "NL"
)

scTWAS_all <- rbindlist(
  list(ADSER, NL),
  use.names = TRUE,
  fill = TRUE
)

message("")
message("scTWAS input summary:")
print(
  scTWAS_all[
    ,
    .(
      rows = .N,
      unique_ENSG = uniqueN(scTWAS_ENSG),
      unique_gene_symbols = uniqueN(gene_symbol)
    ),
    by = scTWAS_dataset
  ]
)


# ==============================================================================
# MODULE 1
# Senescence-associated low-grade adenoma vs normal
# ==============================================================================

M1 <- run_module(
  module_id = "M1",
  module_name = "senescence_low_grade_vs_normal",
  comparison = paste0(
    "Low-grade adenoma vs normal among senescence-high, ",
    "non-stromal Visium spots"
  ),
  spatial_file = file.path(
    SPATIAL_DIR,
    "de_senescence_low_grade_vs_normal.csv"
  ),
  sctwas_all = scTWAS_all,
  out_dir = OUT_DIR
)


# ==============================================================================
# MODULE 2
# Advanced vs nonadvanced adenoma
# ==============================================================================

M2 <- run_module(
  module_id = "M2",
  module_name = "advanced_vs_nonadvanced",
  comparison = paste0(
    "Large vs small adenoma within dysplastic spatial clusters C2/C4"
  ),
  spatial_file = file.path(
    SPATIAL_DIR,
    "de_advanced_vs_nonadvanced.csv"
  ),
  sctwas_all = scTWAS_all,
  out_dir = OUT_DIR
)


# ==============================================================================
# MODULE 3
# Stemness-high vs stemness-low dysplastic regions
# ==============================================================================

M3 <- run_module(
  module_id = "M3",
  module_name = "stemness_high_vs_low",
  comparison = paste0(
    "Stemness-high vs stemness-low within dysplastic spatial clusters C2/C4"
  ),
  spatial_file = file.path(
    SPATIAL_DIR,
    "de_stemness_high_vs_low.csv"
  ),
  sctwas_all = scTWAS_all,
  out_dir = OUT_DIR
)


# ==============================================================================
# MODULE 4
# Spatial cluster markers -- cell-type-matched analysis
# ==============================================================================

M4 <- run_cluster_marker_module(
  module_id = "M4",
  module_name = "cluster_markers_celltype_matched",
  comparison = paste0(
    "scTWAS gene is a marker of a biologically matching Visium cluster ",
    "for its scTWAS cell type"
  ),
  spatial_file = file.path(
    SPATIAL_DIR,
    "cluster_markers.csv"
  ),
  sctwas_all = scTWAS_all,
  out_dir = OUT_DIR,
  celltype_cluster_map = CELLTYPE_CLUSTER_MAP
)


# ------------------------------------------------------------------------------
# 3. Combined summary across all four modules
# ------------------------------------------------------------------------------

all_summary_dataset <- rbindlist(
  list(
    M1$summary_dataset,
    M2$summary_dataset,
    M3$summary_dataset,
    M4$summary_dataset
  ),
  use.names = TRUE,
  fill = TRUE
)

all_summary_celltype <- rbindlist(
  list(
    M1$summary_celltype,
    M2$summary_celltype,
    M3$summary_celltype,
    M4$summary_celltype
  ),
  use.names = TRUE,
  fill = TRUE
)

fwrite(
  all_summary_dataset,
  file.path(OUT_DIR, "all_4modules_summary_by_dataset.csv"),
  na = ""
)

fwrite(
  all_summary_celltype,
  file.path(OUT_DIR, "all_4modules_summary_by_celltype.csv"),
  na = ""
)

message("")
message("======================================================================")
message("All four modules finished.")
message("Output directory: ", OUT_DIR)
message("======================================================================")
print(all_summary_dataset)
