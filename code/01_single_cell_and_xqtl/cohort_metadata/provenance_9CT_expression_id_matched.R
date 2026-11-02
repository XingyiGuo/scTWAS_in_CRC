#!/usr/bin/env Rscript

# provenance_9CT_expression_id_matched.R
# ------------------------------------------------------------------------------
# Subject:
#   Results/TWAS_scRNA_meta/9CT_expression_id_matched_sample_metadata.csv
#   Results/TWAS_scRNA_meta/9CT_expression_id_matched_cell_metadata.csv
#   Results/TWAS_scRNA_meta/9CT_expression_id_match_summary.csv
#   (all three written together, 2026-09-07 21:24)
#
# WHERE THE FILE COMES FROM
# -------------------------
# The script that WROTE these three files is not on any mounted drive and not in
# the session history. Searching L:\scRNA-seq\** , W:\VUMC_projects\scTWAS\**
# (*.R, *.py, *.sh, *.Rmd) for the string "9CT_expression_id_matched" returns
# only CONSUMERS - 08_umap_ADSER_model_cells.R, 09_umap05e_ADSER_polyp_celltype.R
# and cells_per_sample_by_celltype.R. The match step itself ran on ACCRE and was
# never synced back. What follows is therefore a RECONSTRUCTION of its logic,
# verified numerically against the shipped file, not a copy of the original.
#
# The step is an inner join of two ingredients, done per cell type:
#
#   A. CELL UNIVERSE  (on disk, reproducible)
#      Results/HTAN_DISVAL_PCA83_ADSER_HPtoSER_AllCT_ASCSSCfiltered_
#        CellOracleR_metadata.addHTA11_6407_8890_8930.csv
#      100,992 cells x 107 SampleIDs, cell type in CellType_cb. This is the
#      ASC/SSC-filtered, HP->SER-recoded ADSER cell table AFTER the three late
#      specimens HTA11_6407_2000001021 / HTA11_8890_2000001011 /
#      HTA11_8930_2000001011 were appended by
#      src/combine_PCA83_generate_metadata_add3.sh (hence the ".add..." suffix;
#      the pre-add file, ..._CellOracleR_metadata.csv, is 3 specimens short).
#      It also carries a spurious CellType_cb level "20" (414 cells) that the
#      9CT match drops.
#
#   B. EXPRESSION SAMPLE LIST  (NOT on disk)
#      For each cell type, the ordered sample list of that cell type's
#      pseudobulk expression matrix - the columns the PredictDB models were
#      trained on. It survives only inside the shipped file, as the pair
#      (expression_celltype, expression_sample_order); the matrices themselves
#      live on ACCRE. Per-cell-type lengths are
#        ABS 84, GOB 84, TAC 84, STM 79, TUF 78, CT 70, EE 66, ASC 49, SSC 34
#      and expression_sample_order runs 1..N within each cell type.
#
#   JOIN: expression sample  ->  cells of that SampleID and that cell type.
#         metadata_match = TRUE for all 628 pairs, Missing_samples = 0 in the
#         match summary, i.e. every expression column found its cells.
#         n_cells = number of joined cells = 81,929 over 628 pairs / 92 samples.
#
# WHAT THIS SCRIPT DOES
# ---------------------
#   verify_n_cells()   recomputes n_cells from ingredient A and compares, so the
#                      shipped column is checked rather than trusted.
#   rebuild_match()    the reconstructed join, for re-running if the expression
#                      sample lists are brought back from ACCRE.
#   universe_counts()  the alternative count matrix - every cell of the donor in
#                      the universe, whether or not that cell type entered the
#                      expression panel. Use this when the question is "how many
#                      cells were sequenced", not "how many cells are behind the
#                      model".
#
# CAVEAT THE TWO COUNTS DISAGREE ON
# ---------------------------------
#   A 0 in the panel-based matrix means "this cell type has no expression column
#   for this donor". For 6 donor x cell-type pairs that is NOT the same as zero
#   cells - all six belong to the three late-added specimens, which entered the
#   panel for the other cell types but never for ASC/SSC:
#       ASC HTA11_6407_2000001021   90 cells      SSC HTA11_8930_2000001011  295
#       ASC HTA11_8930_2000001011   40 cells      SSC HTA11_6407_2000001021   88
#       ASC HTA11_8890_2000001011   25 cells      SSC HTA11_8890_2000001011    2
#   Everywhere else a panel 0 really is 0 cells.
# ------------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(readxl)      # only for the donor list used by the gap report
})

# ---- paths -------------------------------------------------------------------

SUPP_XLSX <- file.path(
  "F:", "Dropbox", "GuoLab", "COREAD_ToDrGuo", "scTWAS", "ms-guo",
  "SupplementaryTables.xlsx"
)

ADSER_ROOT <- file.path(
  "L:", "scRNA-seq", "VUMC_HTAN", "VUMC_HTAN_DIS_VAL_ADSER", "Results"
)

# Ingredient A: the cell universe (with the three late specimens appended).
UNIVERSE_CSV <- file.path(
  ADSER_ROOT,
  paste0(
    "HTAN_DISVAL_PCA83_ADSER_HPtoSER_AllCT_ASCSSCfiltered_",
    "CellOracleR_metadata.addHTA11_6407_8890_8930.csv"
  )
)

MATCH_DIR     <- file.path(ADSER_ROOT, "TWAS_scRNA_meta")
SAMPLE_META   <- file.path(MATCH_DIR, "9CT_expression_id_matched_sample_metadata.csv")
MATCH_SUMMARY <- file.path(MATCH_DIR, "9CT_expression_id_match_summary.csv")

CT9 <- c("ABS", "GOB", "CT", "EE", "STM", "TUF", "TAC", "ASC", "SSC")

stopifnot(file.exists(UNIVERSE_CSV), file.exists(SAMPLE_META))

# ---- ingredient A: cells per SampleID x cell type ----------------------------

# Only SampleID and CellType_cb are needed. The universe CSV is ~40 MB over ~52
# columns; reading it whole as character costs several hundred MB of R memory
# for no benefit, so take the two columns and nothing else.
load_universe <- function(path = UNIVERSE_CSV) {

  u <- if (requireNamespace("data.table", quietly = TRUE)) {
    as.data.frame(
      data.table::fread(
        path,
        select     = c("SampleID", "CellType_cb"),
        colClasses = "character",
        showProgress = FALSE
      )
    )
  } else {
    read.csv(
      path,
      colClasses  = "character",
      check.names = FALSE
    )[, c("SampleID", "CellType_cb")]
  }

  message(
    "Universe: ", nrow(u), " cells, ",
    dplyr::n_distinct(u$SampleID), " SampleIDs, cell types: ",
    paste(sort(unique(u$CellType_cb)), collapse = ", ")
  )

  u %>%
    filter(CellType_cb %in% CT9) %>%     # drops the spurious "20" level
    count(
      celltype      = CellType_cb,
      Expression_id = SampleID,
      name          = "n_cells_universe"
    )
}

# ---- (1) verify the shipped n_cells column -----------------------------------

verify_n_cells <- function(universe = load_universe()) {

  shipped <- read.csv(SAMPLE_META, check.names = FALSE) %>%
    transmute(
      celltype      = as.character(expression_celltype),
      Expression_id = as.character(SampleID),
      n_cells       = as.integer(n_cells)
    )

  chk <- shipped %>%
    left_join(universe, by = c("celltype", "Expression_id"))

  n_missing <- sum(is.na(chk$n_cells_universe))
  n_differ  <- sum(chk$n_cells != chk$n_cells_universe, na.rm = TRUE)

  message(
    "Shipped pairs: ", nrow(chk),
    " | absent from universe: ", n_missing,
    " | n_cells disagreeing: ", n_differ,
    " | total cells: ", sum(chk$n_cells)
  )

  stopifnot(n_missing == 0, n_differ == 0)

  # Reproduce 9CT_expression_id_match_summary.csv from the shipped file.
  smry <- shipped %>%
    group_by(CellType = celltype) %>%
    summarise(
      Expression_samples = n(),
      Matched_samples    = n(),            # Missing_samples is 0 throughout
      Missing_samples    = 0L,
      Matched_cell_rows  = sum(n_cells),
      .groups = "drop"
    )

  if (file.exists(MATCH_SUMMARY)) {
    ship_smry <- read.csv(MATCH_SUMMARY, check.names = FALSE)
    message(
      "Match summary reproduced: ",
      isTRUE(all.equal(
        as.data.frame(smry[order(smry$CellType), ]),
        as.data.frame(ship_smry[order(ship_smry$CellType), ]),
        check.attributes = FALSE
      ))
    )
  }

  invisible(chk)
}

# ---- (2) the reconstructed join ----------------------------------------------
# expr_samples: data.frame(celltype, Expression_id, expression_sample_order)
#               i.e. ingredient B, recovered from ACCRE or from the shipped file.

rebuild_match <- function(expr_samples, universe = load_universe()) {

  expr_samples %>%
    left_join(universe, by = c("celltype", "Expression_id")) %>%
    mutate(
      metadata_match = !is.na(n_cells_universe),
      n_cells        = tidyr::replace_na(n_cells_universe, 0L)
    ) %>%
    select(
      expression_celltype     = celltype,
      expression_sample_order,
      SampleID                = Expression_id,
      metadata_match,
      n_cells
    ) %>%
    arrange(
      match(expression_celltype, CT9),
      expression_sample_order
    )
}

# Ingredient B as preserved inside the shipped file, for re-running the join.
expr_samples_from_shipped <- function() {
  read.csv(SAMPLE_META, check.names = FALSE) %>%
    transmute(
      celltype                = as.character(expression_celltype),
      Expression_id           = as.character(SampleID),
      expression_sample_order = as.integer(expression_sample_order)
    )
}

# ---- (3) universe-based counts, as a samples x cell-types matrix -------------

universe_counts <- function(ids, universe = load_universe()) {

  universe %>%
    filter(Expression_id %in% ids) %>%
    pivot_wider(
      names_from  = celltype,
      values_from = n_cells_universe,
      values_fill = 0
    ) %>%
    right_join(tibble(Expression_id = ids), by = "Expression_id") %>%
    mutate(across(-Expression_id, ~ tidyr::replace_na(.x, 0L))) %>%
    select(all_of(c("Expression_id", CT9))) %>%
    arrange(match(Expression_id, ids))
}

# Donor x cell-type pairs where the panel says 0 but the donor does have cells.
# `ids` restricts the comparison to a donor set - pass the 84 Expression_ids of
# Supplementary Table S2 to see only the pairs that affect the delivered matrix;
# leave NULL to see all 107 universe samples (most of which are simply not
# TWAS donors at all and are expected to be absent from the panel).
panel_vs_universe_gaps <- function(universe = load_universe(), ids = NULL) {

  shipped <- read.csv(SAMPLE_META, check.names = FALSE) %>%
    transmute(
      celltype      = as.character(expression_celltype),
      Expression_id = as.character(SampleID)
    ) %>%
    mutate(in_panel = TRUE)

  out <- universe %>%
    left_join(shipped, by = c("celltype", "Expression_id")) %>%
    filter(is.na(in_panel), n_cells_universe > 0)

  if (!is.null(ids)) {
    out <- out %>% filter(Expression_id %in% ids)
  }

  out %>%
    select(celltype, Expression_id, n_cells_universe) %>%
    arrange(celltype, desc(n_cells_universe))
}

# ---- run ---------------------------------------------------------------------

if (sys.nframe() == 0L) {

  uni <- load_universe()

  verify_n_cells(uni)

  rebuilt <- rebuild_match(expr_samples_from_shipped(), uni)

  shipped <- read.csv(SAMPLE_META, check.names = FALSE)[
    , c("expression_celltype", "expression_sample_order",
        "SampleID", "metadata_match", "n_cells")
  ]
  # metadata_match is stored as the Python literal "True" - further evidence the
  # original match step was Python, not R. Coerce before comparing.
  shipped$metadata_match <- shipped$metadata_match %in% c("True", "TRUE", "true")

  norm <- function(d) {
    d <- as.data.frame(d[order(d$expression_celltype,
                               d$expression_sample_order), ])
    rownames(d) <- NULL
    d
  }

  message(
    "Rebuilt join identical to shipped file: ",
    isTRUE(all.equal(norm(rebuilt), norm(shipped), check.attributes = FALSE))
  )

  # Pairs where the panel-based matrix shows 0 but cells exist. Restrict to the
  # 84 Supplementary Table S2 polyp donors - those are the rows of the delivered
  # matrix. Without the restriction the same call returns 57 pairs, the extra
  # ones belonging to the 8 panel SampleIDs that have no matched genotype and so
  # never reach the delivered matrix at all.
  donor_ids <- readxl::read_excel(
    SUPP_XLSX,
    sheet = "scRNA-seq patients",
    skip  = 1
  ) %>%
    filter(Tissue == "Colon polyps") %>%
    pull(Expression_id) %>%
    as.character()

  gaps <- panel_vs_universe_gaps(uni, ids = donor_ids)
  message("Panel shows 0 but cells exist, among the 84 TWAS polyp donors:")
  print(as.data.frame(gaps))
}
