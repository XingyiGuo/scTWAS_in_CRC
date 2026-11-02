#!/usr/bin/env Rscript

# cells_per_sample_by_celltype.R
# ------------------------------------------------------------------------------
# Goal:
#   Build a samples x cell-types matrix of CELL COUNTS for every donor in the
#   scTWAS donor list (Supplementary Table S2), covering both tissues:
#
#     rows    = Expression_id (one donor / sample)
#     columns = ABS, GOB, CT, EE, STM, TUF, TAC, ASC, SSC
#     values  = number of single cells of that cell type from that sample
#
# Donor list (the row universe):
#   SupplementaryTables.xlsx, sheet "scRNA-seq patients".
#   The sheet title sits on Excel row 1 and the real header on Excel row 2,
#   so it is read with skip = 1. 465 donors:
#     Colon polyps / VUMC_HTAN        n =  84   (9 cell types, ADSER models)
#     Normal colon / VUMC_HTAN_GCA    n = 103   (7 cell types)
#     Normal colon / WellcomeSanger   n = 278   (7 cell types)
#
# Where the counts come from (one source per cohort; none of them is a
# ready-made combined table):
#
#   (1) COLON POLYPS - VUMC_HTAN (ADSER)
#       9CT_expression_id_matched_sample_metadata.csv
#       Already one row per expression_celltype x SampleID with an n_cells
#       column (628 rows, 92 SampleIDs). Verified: n_cells equals the number of
#       rows for the same pair in the companion cell-level file
#       9CT_expression_id_matched_cell_metadata.csv (81,929 cells) for all 628
#       pairs, so n_cells is a plain cell count and is used as-is.
#       Set COUNT_ADSER_FROM_CELLS = TRUE below to recount from the cell-level
#       file instead of trusting n_cells.
#
#       SPECIMEN vs PARTICIPANT. The ADSER expression panel has 92 SPECIMENS
#       from 84 PARTICIPANTS - 8 participants contributed a second specimen,
#       each with its own expression column and its own PredictDB-trained
#       sample. Supplementary Table S2 lists one specimen per participant, so
#       joining the counts onto the donor list alone silently discards those 8
#       second specimens (1,490 cells, 1,404 of them SSC). The 8 pairs are
#         HTA11_10466_2000001021 / ..._2000001011   HTA11_3997_2000001021 / ..._2000001011
#         HTA11_11156_2000001011 / ..._2000001021   HTA11_6436_2000001011 / ..._2000001021
#         HTA11_12480_2000001031 / ..._2000001011   HTA11_6801_2000001021 / ..._2000001011
#         HTA11_12910_2000001031 / ..._2000001011   HTA11_8755_2000001021 / ..._2000001011
#       (left = the extra specimen, right = the one named in S2). The polyp
#       sheet therefore carries all 92 specimens, flagged by In_SupplTable_S2,
#       and a second sheet aggregates them to the 84 participants. Both total
#       81,929 cells, matching "Model cells in key table: 81929" and the
#       cell-type x cohort table printed by 09_umap05e_ADSER_polyp_celltype.R.
#       Participant is parsed from the specimen id as HTA11_<participant>_<spec>.
#
#   (2) NORMAL COLON - VUMC_HTAN_GCA (GSE266546)
#       No precomputed count table exists. Cells are counted directly from
#       integrated_seurat_final2_reviewed.res0.8.rm23_sub_metadata.csv, grouped
#       by SampleID x celltype_annot_reviewed - the same two columns that
#       07_sample_level_pseudobulk_counts.R aggregates on to build the
#       pseudobulk matrices, so the counts describe exactly the cells behind
#       the expression models. That object carries 157 samples / 120,891 cells;
#       restricting to the 103 genotype-matched GCA donors gives 69,177 cells.
#       SampleID in this file is the same identifier as Expression_id in the
#       supplementary sheet (e.g. "7408-YX-8").
#       No specimen/participant ambiguity here: the 103 selected SampleIDs are
#       103 distinct PatientIDs and are exactly the 103 column names of
#       pseudobulk_counts_ABS.tpmlog2.qc.qn.rmlowcounts.inv.N103.matched.
#       expression_id.geno_order.csv, so rows are 1:1 with the expression panel.
#       (30 further specimens from those same patients exist in the object but
#       were not selected into the panel; excluding them is correct.)
#       Per-cell-type totals over all 157 samples reproduce the "05e available"
#       figures printed by 05_annotation_UMAP.R exactly: ABS 74,320, CT 4,539,
#       EE 1,207, GOB 23,036, STM 4,046, TAC 11,807, TUF 1,936. The "Review"
#       exclusion below drops 0 cells at this stage - clusters 8 and 10 are
#       Review only in the earlier dominant-cell-type pass (05c) and are
#       resolved to CT and STM in celltype_annot_reviewed.
#       13 of the 103 donors are TERMINAL ILEUM, not colon (3,540 cells); the
#       donor sheet records this in Location, and Tissue = "Normal colon" is a
#       cohort label. Filter on Location if a colon-only table is wanted.
#
#   (3) NORMAL COLON - WellcomeSanger
#       06_subject_x_LauLabel_ncells.csv - one row per individual_id x LauLabel
#       with n_cells (1942 rows = 278 subjects x 7 cell types, minus 4 absent
#       combinations). Verified: n_cells equals the row count for the same pair
#       in 05_Epithelial_cell_metadata_LauLabel.csv.gz (358,207 cells) for all
#       1942 pairs. This file is keyed by individual_id, which is the
#       Genotype_id of the supplementary sheet (IBDverse-1-...), NOT the
#       Expression_id (OTARscRNA...), so it is remapped via the donor list.
#
# Conventions in the output:
#   - A cell type that was profiled but contributed no cells is 0.
#   - ASC and SSC are left NA for normal-colon rows: those two populations exist
#     only in the polyp annotation, so a 0 there would wrongly read as
#     "profiled, none found".
#
# Outputs (written to OUT_DIR):
#   scTWAS_cells_per_sample_by_celltype.xlsx   5 sheets (Summary, combined,
#                                              one per cohort)
#   scTWAS_cells_per_sample_by_celltype.csv    flat 465-row combined table
# ------------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(openxlsx)
})

# ---- options -----------------------------------------------------------------

# TRUE  = recount the ADSER cells from 9CT_expression_id_matched_cell_metadata.csv
# FALSE = use the precomputed n_cells column (default; verified identical)
COUNT_ADSER_FROM_CELLS <- FALSE

# Which rows the delivered table carries.
#   "sheet_only" (default, as requested): exactly the Expression_ids listed in
#       the "scRNA-seq patients" sheet. The 8 ADSER second specimens are NOT
#       rows and their cells are NOT added to the listed specimen of the same
#       participant, so the polyp block totals 80,439 rather than the panel's
#       81,929. The omitted specimens are written to their own sheet so the
#       1,490-cell difference is visible rather than silent.
#   "all_panel": every specimen in the expression panel (92 polyp rows), plus a
#       participant-level sheet; polyp totals then match the cell-type x cohort
#       table printed by 09_umap05e_ADSER_polyp_celltype.R.
ROW_SET <- "sheet_only"
stopifnot(ROW_SET %in% c("sheet_only", "all_panel"))

# TRUE: a sample x cell type that has NO column in that cell type's final
# PEER-residual matrix is written as NA rather than 0, so "not modelled for
# this cell type" is visually distinct from "modelled, zero cells". Membership
# is read from the column headers of
# expression_for_scTWAS/<cohort>/*PEER_residuals*.csv.gz - the matrices the
# scTWAS models were actually fitted on.
MASK_TO_FINAL_PANEL <- TRUE

# ---- paths -------------------------------------------------------------------

SUPP_XLSX <- file.path(
  "F:", "Dropbox", "GuoLab", "COREAD_ToDrGuo", "scTWAS", "ms-guo",
  "SupplementaryTables.xlsx"
)

# The 9CT match files moved out of VUMC_HTAN_DIS_VAL_ADSER\Results\TWAS_scRNA_meta
# into the shared TWAS_scRNA_samples_cells_meta root; try the current location
# first and fall back to the old one.
ADSER_DIR <- Filter(
  dir.exists,
  c(
    file.path("L:", "scRNA-seq", "TWAS_scRNA_samples_cells_meta"),
    file.path("L:", "scRNA-seq", "VUMC_HTAN", "VUMC_HTAN_DIS_VAL_ADSER",
              "Results", "TWAS_scRNA_meta")
  )
)[1]
if (is.na(ADSER_DIR)) stop("9CT match metadata directory not found")
ADSER_SAMPLE_META <- file.path(
  ADSER_DIR, "9CT_expression_id_matched_sample_metadata.csv"
)
ADSER_CELL_META <- file.path(
  ADSER_DIR, "9CT_expression_id_matched_cell_metadata.csv"
)

GCA_CELL_META <- file.path(
  "L:", "scRNA-seq", "VUMC_HTAN_GSE266546", "raw", "rds2",
  "integrated_seurat_final2_reviewed.res0.8.rm23_sub_metadata.csv"
)

SANGER_NCELLS <- file.path(
  "L:", "scRNA-seq", "WellcomeSanger",
  "Healthy_TI_then_rectum_278samples_outputs_LauLabel",
  "06_subject_x_LauLabel_ncells.csv"
)

# Final per-cell-type PEER-residual matrices - the expression actually used for
# scTWAS. One subdirectory per cohort; columns are samples, keyed by
# Expression_id for ADSER and VUMC_HTAN_GCA and by Genotype_id for
# WellcomeSanger. The column sets define which sample x cell type combinations
# have a model at all.
PANEL_DIR <- file.path(
  "L:", "scRNA-seq", "TWAS_scRNA_samples_cells_meta", "expression_for_scTWAS"
)
PANEL_SUBDIR <- c(
  "VUMC_HTAN"      = "ADSER",
  "VUMC_HTAN_GCA"  = "VUMC_HTAN_GCA",
  "WellcomeSanger" = "WellcomeSanger"
)

OUT_DIR <- ADSER_DIR

stopifnot(
  file.exists(SUPP_XLSX),
  file.exists(ADSER_SAMPLE_META),
  file.exists(GCA_CELL_META),
  file.exists(SANGER_NCELLS)
)

# ---- cell-type column order --------------------------------------------------

CT9 <- c("ABS", "GOB", "CT", "EE", "STM", "TUF", "TAC", "ASC", "SSC")
CT7 <- c("ABS", "GOB", "CT", "EE", "STM", "TUF", "TAC")

# Donor annotation carried into every output. Set after the sheet is read, so
# the delivered table always carries EVERY column of "scRNA-seq patients"
# (currently Tissue, Project, Genotype_id, Age, Sex, Location, PC1-PC5) with
# Expression_id first. Adding a column to the sheet propagates automatically.
ANNOT_COLS <- NULL   # assigned in section 0

# ==============================================================================
# 0. Donor list - defines which samples appear and in which order
# ==============================================================================

donors <- read_excel(
  SUPP_XLSX,
  sheet = "scRNA-seq patients",
  skip  = 1                       # title paragraph occupies Excel row 1
) %>%
  mutate(
    across(c(Expression_id, Genotype_id, Tissue, Project), as.character)
  )

ANNOT_COLS <- c("Expression_id", setdiff(names(donors), "Expression_id"))

message(
  "Donor-sheet columns carried: ", paste(ANNOT_COLS, collapse = ", ")
)

message(
  "Donors: ", nrow(donors), " (",
  paste(
    sprintf("%s=%d", names(table(donors$Project)), as.integer(table(donors$Project))),
    collapse = ", "
  ),
  ")"
)

polyp_ids <- donors %>% filter(Tissue  == "Colon polyps")   %>% pull(Expression_id)
gca_ids   <- donors %>% filter(Project == "VUMC_HTAN_GCA")  %>% pull(Expression_id)
sanger    <- donors %>% filter(Project == "WellcomeSanger") %>%
  select(Expression_id, Genotype_id)

# Helper: long (sample, celltype, n) -> wide matrix on a fixed row/column order,
# absent combinations becoming 0.
to_wide <- function(long_df, ids, cts) {
  long_df %>%
    filter(Expression_id %in% ids) %>%
    group_by(Expression_id, celltype) %>%
    summarise(n_cells = sum(n_cells), .groups = "drop") %>%
    pivot_wider(
      names_from  = celltype,
      values_from = n_cells,
      values_fill = 0
    ) %>%
    right_join(tibble(Expression_id = ids), by = "Expression_id") %>%
    mutate(across(-Expression_id, ~ tidyr::replace_na(.x, 0))) %>%
    # A cell type with no rows at all in this subset yields no column; add it
    # as 0 so the output always carries the full cell-type set.
    bind_cols(
      setNames(
        as.data.frame(matrix(0L, nrow = length(ids),
                             ncol = length(setdiff(cts, names(.))))),
        setdiff(cts, names(.))
      )
    ) %>%
    select(all_of(c("Expression_id", cts))) %>%
    arrange(match(Expression_id, ids))
}

# ---- final-panel membership --------------------------------------------------
# Column headers of the PEER-residual matrices, one set per cohort x cell type.
# Reading only the first line keeps this cheap (the Sanger files are ~30 MB each).

panel_membership <- function(subdir) {

  d <- file.path(PANEL_DIR, subdir)
  if (!dir.exists(d)) stop("Panel directory not found: ", d)

  fs <- list.files(d, pattern = "PEER_residuals.*\\.csv\\.gz$", full.names = TRUE)
  if (!length(fs)) stop("No PEER_residuals files in ", d)

  out <- lapply(fs, function(f) {
    # Cell type from either naming scheme:
    #   pseudobulk_counts_<CT>.tpmlog2...PEER_residuals.expression_id.csv.gz
    #   <CT>.TPM.log2.QN.INV.PEER_residuals.csv.gz
    b  <- basename(f)
    ct <- if (grepl("^pseudobulk_counts_", b)) {
      sub("^pseudobulk_counts_([A-Za-z0-9]+)\\..*$", "\\1", b)
    } else {
      sub("^([A-Za-z0-9]+)\\..*$", "\\1", b)
    }
    con <- gzfile(f, "rt"); on.exit(close(con))
    hdr <- scan(con, what = "", sep = ",", nlines = 1, quiet = TRUE)
    tibble(celltype = ct, key = gsub('^"|"$', "", hdr[-1]))   # drop the gene column
  })

  bind_rows(out)
}

PANEL <- lapply(PANEL_SUBDIR, panel_membership)

message(
  "Final panels: ",
  paste(sprintf("%s %d cell types / %d sample-columns",
                names(PANEL),
                vapply(PANEL, function(x) dplyr::n_distinct(x$celltype), 1L),
                vapply(PANEL, nrow, 1L)),
        collapse = "; ")
)

# Set a cell to NA where `keys` has no column in that cell type's residual
# matrix. `keys` is the sample identifier used by that cohort's files, aligned
# row-for-row with wide_df.
mask_to_panel <- function(wide_df, subdir_name, keys) {

  if (!MASK_TO_FINAL_PANEL) return(wide_df)

  memb <- PANEL[[subdir_name]]
  cts  <- setdiff(names(wide_df), "Expression_id")

  for (ct in cts) {
    present <- keys %in% memb$key[memb$celltype == ct]
    wide_df[[ct]][!present] <- NA_integer_
  }
  wide_df
}

# ==============================================================================
# 1. Colon polyps - VUMC_HTAN (ADSER), 9 cell types
# ==============================================================================

if (COUNT_ADSER_FROM_CELLS) {

  adser_long <- read.csv(
    ADSER_CELL_META,
    colClasses = "character",
    check.names = FALSE
  ) %>%
    count(
      Expression_id = SampleID,
      celltype      = expression_celltype,
      name          = "n_cells"
    )

} else {

  adser_long <- read.csv(ADSER_SAMPLE_META, check.names = FALSE) %>%
    transmute(
      Expression_id = as.character(SampleID),
      celltype      = as.character(expression_celltype),
      n_cells       = as.integer(n_cells)
    )
}

# Participant id, parsed from the specimen id: HTA11_<participant>_<specimen>.
participant_of <- function(x) sub("^(HTA11_[^_]+)_.*$", "\\1", x)

# Rows = every specimen in the expression panel (92), S2-listed first so the
# delivered order is unchanged, the 8 extra specimens appended.
panel_ids <- unique(adser_long$Expression_id)
extra_ids <- sort(setdiff(panel_ids, polyp_ids))

polyp_rows <- if (ROW_SET == "all_panel") c(polyp_ids, extra_ids) else polyp_ids

message(
  "ADSER source rows: ", nrow(adser_long),
  " across ", length(panel_ids), " panel specimens from ",
  dplyr::n_distinct(participant_of(panel_ids)), " participants; ",
  length(extra_ids), " are second specimens of a listed participant (",
  if (ROW_SET == "all_panel") "kept as their own rows, flagged" else
    "EXCLUDED - ROW_SET is 'sheet_only'", ")"
)
stopifnot(length(setdiff(participant_of(extra_ids),
                         participant_of(polyp_ids))) == 0)

polyp_wide <- to_wide(adser_long, polyp_rows, CT9) %>%
  { mask_to_panel(., "VUMC_HTAN", .$Expression_id) }

# Participant-level view: sum a participant's specimens (all_panel only).
polyp_participant_wide <- adser_long %>%
  mutate(Expression_id = participant_of(Expression_id)) %>%
  to_wide(participant_of(polyp_ids), CT9) %>%
  rename(Participant = Expression_id)

# The specimens left out under "sheet_only", so the omission is auditable.
excluded_wide <- to_wide(adser_long, extra_ids, CT9) %>%
  mutate(
    Participant            = participant_of(Expression_id),
    Listed_specimen_in_S2  = polyp_ids[match(Participant, participant_of(polyp_ids))],
    Cells_not_counted      = rowSums(across(all_of(CT9)))
  ) %>%
  select(all_of(c("Expression_id", "Participant", "Listed_specimen_in_S2",
                  CT9, "Cells_not_counted")))

# ==============================================================================
# 2. Normal colon - VUMC_HTAN_GCA (GSE266546), 7 cell types
#    Counted from cell-level metadata on the same keys the pseudobulk script
#    aggregates on: SampleID x celltype_annot_reviewed.
# ==============================================================================

gca_cells <- read.csv(
  GCA_CELL_META,
  colClasses  = "character",
  check.names = FALSE
)

message(
  "GSE266546 cell metadata: ", nrow(gca_cells), " cells across ",
  dplyr::n_distinct(gca_cells$SampleID), " samples"
)

gca_long <- gca_cells %>%
  filter(
    !is.na(SampleID), SampleID != "",
    !is.na(celltype_annot_reviewed),
    celltype_annot_reviewed != "",
    celltype_annot_reviewed != "Review"     # same exclusions as script 07
  ) %>%
  count(
    Expression_id = SampleID,
    celltype      = celltype_annot_reviewed,
    name          = "n_cells"
  )

gca_wide <- to_wide(gca_long, gca_ids, CT7) %>%
  { mask_to_panel(., "VUMC_HTAN_GCA", .$Expression_id) }

# ==============================================================================
# 3. Normal colon - WellcomeSanger, 7 cell types
#    Keyed by individual_id == Genotype_id; remap to Expression_id.
# ==============================================================================

sanger_long <- read.csv(SANGER_NCELLS, check.names = FALSE) %>%
  transmute(
    Genotype_id = as.character(individual_id),
    celltype    = as.character(LauLabel),
    n_cells     = as.integer(n_cells)
  ) %>%
  inner_join(sanger, by = "Genotype_id") %>%
  select(Expression_id, celltype, n_cells)

message(
  "WellcomeSanger rows mapped to donor list: ", nrow(sanger_long),
  " / ", dplyr::n_distinct(sanger_long$Expression_id), " donors"
)

# Sanger residual matrices are keyed by Genotype_id, not Expression_id.
sanger_wide <- to_wide(sanger_long, sanger$Expression_id, CT7) %>%
  { mask_to_panel(., "WellcomeSanger",
                  sanger$Genotype_id[match(.$Expression_id, sanger$Expression_id)]) }

# ==============================================================================
# 4. Attach donor annotation, totals; assemble sheets
# ==============================================================================

# Total_cells sums only the modelled cell types; N_celltypes_modelled counts the
# non-NA entries and N_celltypes_present the non-NA entries that are > 0.
decorate <- function(wide_df) {
  cts <- setdiff(names(wide_df), "Expression_id")
  donors %>%
    select(all_of(ANNOT_COLS)) %>%
    inner_join(wide_df, by = "Expression_id") %>%
    arrange(match(Expression_id, wide_df$Expression_id)) %>%
    mutate(
      Total_cells          = rowSums(across(all_of(cts)), na.rm = TRUE),
      N_celltypes_modelled = rowSums(across(all_of(cts), ~ !is.na(.x))),
      N_celltypes_present  = rowSums(across(all_of(cts),
                                            ~ !is.na(.x) & .x > 0))
    )
}

gca_tab    <- decorate(gca_wide)
sanger_tab <- decorate(sanger_wide)

# Polyps: donor annotation is joined on PARTICIPANT, so a second specimen
# inherits its participant's Genotype_id / Age / Sex / Location from S2.
polyp_annot <- donors %>%
  filter(Tissue == "Colon polyps") %>%
  select(all_of(ANNOT_COLS)) %>%
  mutate(Participant = participant_of(Expression_id)) %>%
  select(-Expression_id)

polyp_tab <- polyp_wide %>%
  mutate(Participant = participant_of(Expression_id)) %>%
  left_join(polyp_annot, by = "Participant") %>%
  mutate(In_SupplTable_S2 = Expression_id %in% polyp_ids) %>%
  arrange(match(Expression_id, polyp_wide$Expression_id)) %>%
  mutate(
    Total_cells          = rowSums(across(all_of(CT9)), na.rm = TRUE),
    N_celltypes_modelled = rowSums(across(all_of(CT9), ~ !is.na(.x))),
    N_celltypes_present  = rowSums(across(all_of(CT9), ~ !is.na(.x) & .x > 0))
  ) %>%
  select(all_of(c(ANNOT_COLS[1], "Participant", "In_SupplTable_S2",
                  ANNOT_COLS[-1], CT9, "Total_cells",
                  "N_celltypes_modelled", "N_celltypes_present")))

polyp_participant_tab <- polyp_participant_wide %>%
  left_join(polyp_annot, by = "Participant") %>%
  mutate(
    N_specimens          = as.integer(table(participant_of(polyp_wide$Expression_id))[Participant]),
    Total_cells          = rowSums(across(all_of(CT9)), na.rm = TRUE),
    N_celltypes_modelled = rowSums(across(all_of(CT9), ~ !is.na(.x))),
    N_celltypes_present  = rowSums(across(all_of(CT9), ~ !is.na(.x) & .x > 0))
  ) %>%
  select(all_of(c("Participant", "N_specimens", ANNOT_COLS[-1], CT9,
                  "Total_cells", "N_celltypes_modelled", "N_celltypes_present")))

# Masking turns 0 into NA; it must never discard a counted cell.
stopifnot(
  sum(polyp_tab$Total_cells) ==
    sum(adser_long$n_cells) -
    (if (ROW_SET == "sheet_only") sum(excluded_wide$Cells_not_counted) else 0)
)

# Combined 9-column view. ASC/SSC stay NA for normal colon: not profiled there.
normal_tab <- bind_rows(gca_tab, sanger_tab) %>%
  mutate(
    Participant      = Expression_id,
    In_SupplTable_S2 = TRUE,
    ASC = NA_integer_, SSC = NA_integer_
  )

ALL_COLS <- c(ANNOT_COLS[1], "Participant", "In_SupplTable_S2", ANNOT_COLS[-1],
              CT9, "Total_cells", "N_celltypes_modelled", "N_celltypes_present")

combined <- bind_rows(polyp_tab, normal_tab) %>%
  select(all_of(ALL_COLS))

csum <- function(w, cts) {
  v <- vapply(cts, function(ct) sum(w[[ct]], na.rm = TRUE), numeric(1))
  as.list(v)
}

summary_tab <- bind_rows(
  tibble(Cohort = "Colon polyps (VUMC_HTAN)",      !!!csum(polyp_wide,  CT9)),
  tibble(Cohort = "Normal colon (VUMC_HTAN_GCA)",  !!!csum(gca_wide,    CT7)),
  tibble(Cohort = "Normal colon (WellcomeSanger)", !!!csum(sanger_wide, CT7))
) %>%
  mutate(
    N_samples = c(nrow(polyp_wide), nrow(gca_wide), nrow(sanger_wide)),
    N_donors  = c(nrow(polyp_participant_wide), nrow(gca_wide), nrow(sanger_wide)),
    Total_cells = rowSums(across(all_of(CT9)), na.rm = TRUE)
  ) %>%
  select(Cohort, N_samples, N_donors, Total_cells, all_of(CT9))

# ---- panel depth -------------------------------------------------------------
# For each cohort x cell type: how many sheet samples have a model column, and
# how shallow the shallowest of them are. Masking removes only sample x cell
# type combinations with no column at all; a sample that contributed a single
# cell AND has a column is still a real model input and is NOT masked, so this
# sheet is where very thin inputs become visible.

depth_of <- function(tab, cohort, cts) {
  bind_rows(lapply(cts, function(ct) {
    v  <- tab[[ct]]
    kept <- v[!is.na(v)]
    tibble(
      Cohort            = cohort,
      CellType          = ct,
      Sheet_samples     = length(v),
      Modelled_samples  = length(kept),
      Masked_to_NA      = sum(is.na(v)),
      Cells_total       = sum(kept),
      Min_cells         = if (length(kept)) min(kept) else NA_integer_,
      Median_cells      = if (length(kept)) as.numeric(stats::median(kept)) else NA_real_,
      N_with_1_cell     = sum(kept == 1),
      N_with_le5_cells  = sum(kept <= 5),
      N_with_le10_cells = sum(kept <= 10)
    )
  }))
}

panel_depth <- bind_rows(
  depth_of(polyp_tab,  "Colon polyps (VUMC_HTAN)",      CT9),
  depth_of(gca_tab,    "Normal colon (VUMC_HTAN_GCA)",  CT7),
  depth_of(sanger_tab, "Normal colon (WellcomeSanger)", CT7)
)

print(as.data.frame(summary_tab))
print(as.data.frame(panel_depth))

# ==============================================================================
# 5. Write
# ==============================================================================

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

out_xlsx <- file.path(OUT_DIR, "scTWAS_cells_per_sample_by_celltype.xlsx")
out_csv  <- file.path(OUT_DIR, "scTWAS_cells_per_sample_by_celltype.csv")

sheets <- list(
  "Summary"                     = summary_tab,
  "Panel depth"                 = panel_depth,
  "All samples (9 cell types)"  = combined,
  "Colon polyps ADSER"          = polyp_tab,
  "Normal colon VUMC_HTAN_GCA"  = gca_tab,
  "Normal colon WellcomeSanger" = sanger_tab
)
sheets[[if (ROW_SET == "all_panel") "Colon polyps by participant"
        else "Specimens not in S2 (excluded)"]] <-
  if (ROW_SET == "all_panel") polyp_participant_tab else excluded_wide

openxlsx::write.xlsx(sheets, file = out_xlsx, overwrite = TRUE)

write.csv(combined, out_csv, row.names = FALSE, na = "")

message("Wrote: ", out_xlsx)
message("Wrote: ", out_csv)
