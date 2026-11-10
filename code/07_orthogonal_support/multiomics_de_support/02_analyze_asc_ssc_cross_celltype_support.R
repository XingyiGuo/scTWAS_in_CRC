## 02_analyze_asc_ssc_cross_celltype_support.R ------------------------------------
## Cross-cell-type scRNA-seq DE support for the scTWAS genes nominated in the
## lesion-specific ADSER models, ASC and SSC.
##
## Why this exists
##   ASC (adenoma-specific) and SSC (serrated-specific) cells do not exist in
##   normal colon, so these genes have no matched normal cell-type contrast.
##   In 01_annotate_scrna_de_support.R the C and E blocks are joined on
##   gene symbol AND cell type, which forces all 47 ASC/SSC pairs to
##   not_tested in those blocks: the only evidence they can receive there is
##   the gene-level L block. That is a limitation of the cell-type join, not an
##   absence of data - every one of these genes is measured in all seven
##   normal cell-type contrasts.
##
##   This script drops the cell-type join for the ASC/SSC genes only and scans
##   each of them against the full panel, so the question becomes "in which
##   epithelial lineage is this lesion-specific risk gene differentially
##   expressed" rather than "is it testable at all".
##
## Contrast panel, 23 contrasts per gene
##   C1-C7 x {MSI-H, MSS}   14   carcinoma vs each normal cell type
##                               (ABS GOB CT EE STM TUF TAC = C1..C7)
##   E1-E7, matching arm     7   polyp cell type vs each normal cell type;
##                               the ASC arm for ASC genes, SSC for SSC genes
##   L, matching lesion      2   carcinoma vs the polyp compartment itself:
##                               L1_MSIH/L1_MSS for ASC, L2_MSIH/L2_MSS for SSC
##   The L block is the only cell-type-matched evidence available for these
##   genes and is carried through from 01 unchanged.
##
## Direction
##   Each of the 47 genes has exactly one ASC or SSC scTWAS pair (asserted
##   below), so the pair's own Z score fixes the direction with no tie to break.
##
## Support rule, identical to 01 and 08 so the labels mean the same thing
##   padj < 0.05 AND |log2FC| > 0.5 AND sign(log2FC) == sign(Z)
##     Up / Down      supported, in the direction the risk Z predicts
##     DE_discordant  passes both cuts but opposes the Z
##     NS             tested, fails a cut
##     not_tested     gene absent from that contrast's DE table
##
## Input
##   data/scTWAS_ADSER_NL_collected_results.csv
##   data/de_C_E_L_contrasts_ASC75_SSC46.csv (or .parquet)
##
## Output
##   table/ASC_SSC_genes_crossceltype_support_long.csv      47 x 23 = 1081 rows
##   table/ASC_SSC_genes_crossceltype_support_matrix.csv    47 rows, tag columns
##   table/ASC_SSC_genes_crossceltype_support_summary.csv   per-gene roll-up
##   table/ASC_SSC_genes_crossceltype_support_counts.csv    per-contrast counts
##   table/ASC_SSC_contrast_legend.csv                      tag -> meaning
## -----------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(data.table)
})

## paths -----------------------------------------------------------------------
PROJ <- "W:/VUMC_projects/scTWAS/VUMC_HTAN_scTWAS_genes_DE_support"
if (!dir.exists(file.path(PROJ, "data"))) PROJ <- "."

DIR_DATA <- file.path(PROJ, "data")
DIR_TAB  <- file.path(PROJ, "table")
dir.create(DIR_TAB, recursive = TRUE, showWarnings = FALSE)

TWAS_CSV <- file.path(DIR_DATA, "scTWAS_ADSER_NL_collected_results.csv")
DE_CSV   <- file.path(DIR_DATA, "de_C_E_L_contrasts_ASC75_SSC46.csv")
DE_PARQ  <- file.path(DIR_DATA, "de_C_E_L_contrasts_ASC75_SSC46.parquet")

Q_CUT   <- 0.05
LFC_CUT <- 0.5

## C/E index -> normal cell type, in the order used by the contrast tags
CE_INDEX <- c(`1` = "ABS", `2` = "GOB", `3` = "CT", `4` = "EE",
              `5` = "STM", `6` = "TUF", `7` = "TAC")

label_support <- function(log2FC, padj, z) {
  out  <- rep("NS", length(log2FC))
  sig  <- !is.na(padj) & padj < Q_CUT & !is.na(log2FC) & abs(log2FC) > LFC_CUT
  conc <- sign(log2FC) == sign(z)
  out[sig &  conc & log2FC > 0] <- "Up"
  out[sig &  conc & log2FC < 0] <- "Down"
  out[sig & !conc]              <- "DE_discordant"
  out[is.na(log2FC)]            <- "not_tested"
  out
}

## 1. the ASC and SSC scTWAS genes --------------------------------------------
twas <- fread(TWAS_CSV)
lesion <- twas[CellType %chin% c("ASC", "SSC")]
setnames(lesion, "CellType", "track")
## one ASC/SSC pair per gene, so the direction needs no tie-break
stopifnot(nrow(lesion) == 47L, uniqueN(lesion$gene) == 47L,
          !anyDuplicated(lesion, by = "gene"),
          lesion[, .N, by = track][order(track)]$N == c(32L, 15L))
message(sprintf("lesion-specific scTWAS genes: %d ASC, %d SSC",
                lesion[track == "ASC", .N], lesion[track == "SSC", .N]))

## 2. DE results ---------------------------------------------------------------
de_cols <- c("gene", "cell_type", "contrast_tag", "stage", "group1",
             "log2FC", "padj")
if (requireNamespace("arrow", quietly = TRUE) && file.exists(DE_PARQ)) {
  de <- as.data.table(arrow::read_parquet(DE_PARQ,
                                          col_select = dplyr::all_of(de_cols)))
} else {
  de <- fread(DE_CSV, select = de_cols)
}
de[, `:=`(log2FC = as.numeric(log2FC), padj = as.numeric(padj))]

## 3. contrast panel per track ------------------------------------------------
## C is shared by both tracks; E and L take the arm matching the lesion.
panel <- function(track) {
  c(paste0("C", 1:7, "_MSIH"), paste0("C", 1:7, "_MSS"),
    paste0("E", 1:7, "_", track),
    if (track == "ASC") c("L1_MSIH", "L1_MSS") else c("L2_MSIH", "L2_MSS"))
}
stopifnot(length(panel("ASC")) == 23L, length(panel("SSC")) == 23L)
stopifnot(all(unlist(lapply(c("ASC", "SSC"), panel)) %chin% de$contrast_tag))

## legend, so the tag columns are readable without the manuscript at hand
legend <- unique(de[contrast_tag %chin% unlist(lapply(c("ASC","SSC"), panel)),
                    .(contrast_tag, stage, group1, contrast_cell_type = cell_type)])
legend[, block := substr(contrast_tag, 1, 1)]
legend[, description := fifelse(
  block == "C", paste0(group1, " carcinoma vs normal ", contrast_cell_type),
  fifelse(block == "E", paste0(group1, " polyp cells vs normal ", contrast_cell_type),
          paste0(group1, " carcinoma vs ", contrast_cell_type, " polyp cells")))]
setorder(legend, block, contrast_tag)
fwrite(legend, file.path(DIR_TAB, "ASC_SSC_contrast_legend.csv"))

## 4. long form: every gene against every contrast in its panel ---------------
de_key <- de[, .(contrast_tag, gene, cell_type, stage, group1, log2FC, padj)]
setkey(de_key, contrast_tag, gene)

long <- rbindlist(lapply(c("ASC", "SSC"), function(tr) {
  g <- lesion[track == tr]
  CJ_gc <- CJ(contrast_tag = panel(tr), gene_name = g$gene_name, unique = TRUE)
  CJ_gc <- merge(CJ_gc, g[, .(gene_name, gene, Dataset, track,
                              Sig_resource, Z_score, P_value)],
                 by = "gene_name", sort = FALSE)
  m <- de_key[.(CJ_gc$contrast_tag, CJ_gc$gene_name),
              on = .(contrast_tag, gene)]
  ## Contrast metadata comes from the tag, not from the matched DE row: a gene
  ## absent from a contrast has no row to read it from, which would leave
  ## contrast_cell_type NA and split that contrast in the count tables.
  CJ_gc[, `:=`(log2FC = round(m$log2FC, 4), padj = m$padj)]
  CJ_gc <- merge(CJ_gc, legend[, .(contrast_tag, contrast_cell_type, stage, group1)],
                 by = "contrast_tag", sort = FALSE)
  stopifnot(!anyNA(CJ_gc$contrast_cell_type))
  CJ_gc[, DEsupport := label_support(m$log2FC, m$padj, Z_score)]
  CJ_gc[]
}))
long[, block := substr(contrast_tag, 1, 1)]
long[, arm := sub("^[CEL][0-9]?_", "", contrast_tag)]
setcolorder(long, c("track", "gene", "gene_name", "Dataset", "Sig_resource",
                    "Z_score", "P_value", "block", "contrast_tag",
                    "contrast_cell_type", "stage", "group1", "arm",
                    "log2FC", "padj", "DEsupport"))
setorder(long, track, -Z_score, block, contrast_tag)
stopifnot(nrow(long) == 47L * 23L)
fwrite(long, file.path(DIR_TAB, "ASC_SSC_genes_crossceltype_support_long.csv"),
       na = "")

## 5. wide matrix, laid out as ASC_SSC_genes_crossceltype_support_matrix_example
##
## Two header rows. Row 1 names each contrast block and spans its three
## columns; row 2 names the columns within the block: log2FC, padj,
## DE_discordant. "DE_discordant" is the example's own name for the
## support-call column and is kept verbatim, although the values it carries are
## Up / Down / DE_discordant / NS / not_tested / not_applicable.
##
## Column order follows the example: C1_MSIH, C1_MSS, C2_MSIH, ... C7_MSS, the
## Any C Support roll-up, then E1..E7, then the L blocks.
##
## E has seven blocks rather than fourteen because each gene is scored in its
## own arm only (E*_ASC for an ASC gene, E*_SSC for an SSC gene), which is what
## the example's "matched cell type with model" wording denotes. The normal
## cell type is added in parentheses so E1..E7 are distinguishable, which the
## single-block example had no need to do.
##
## The L section carries all four tags, so both prefixes from the example are
## present, but a gene only ever receives values from its own lesion (L1 for
## ASC, L2 for SSC). The other lesion's three cells read not_applicable rather
## than a cross-lesion contrast - the same restriction applied throughout this
## script, and the reason its matched-L count is 6 where the per-pair table in
## 01 reports 11.
##
## Two files are written from one assembly: the grouped-header CSV for reading
## and pasting into the supplement, and a *_flat.csv with unique
## machine-readable names for downstream code.

## value lookup: one row per gene x contrast actually scored
val <- long[, .(gene, contrast_tag, log2FC, padj, DEsupport)]
setkey(val, gene, contrast_tag)

id <- unique(long[, .(track, gene, gene_name, Z_score, P_value)])
stopifnot(nrow(id) == 47L)
id[, absZ := abs(Z_score)]
setorder(id, track, -absZ)
id[, absZ := NULL]

## block spec: a group header, a flat prefix, and the contrast tag to pull for
## each track (NA = not applicable to that track)
mk_block <- function(header, flat, tag_asc, tag_ssc)
  list(header = header, flat = flat, tag_asc = tag_asc, tag_ssc = tag_ssc)

blocks <- list()
for (i in 1:7) {
  ct <- CE_INDEX[[as.character(i)]]
  for (arm in c("MSIH", "MSS")) {
    tg <- sprintf("C%d_%s", i, arm)
    blocks[[length(blocks) + 1L]] <- mk_block(
      sprintf("C%d_%s_vs_NL(%s)", i, arm, ct), tg, tg, tg)
  }
}
ANY_C_AT <- length(blocks)          # the roll-up column sits after C7_MSS
for (i in 1:7) {
  ct <- CE_INDEX[[as.character(i)]]
  blocks[[length(blocks) + 1L]] <- mk_block(
    sprintf("E%d_Polyps (ASC/SSC, matched cell type with model)_vs_NL(%s)", i, ct),
    sprintf("E%d", i), sprintf("E%d_ASC", i), sprintf("E%d_SSC", i))
}
for (spec in list(c("1", "MSIH"), c("1", "MSS"), c("2", "MSIH"), c("2", "MSS"))) {
  n <- spec[1]; arm <- spec[2]; tg <- paste0("L", n, "_", arm)
  blocks[[length(blocks) + 1L]] <- mk_block(
    sprintf("L%s_%s_vs_Polyps (ASC/SSC, matched cell type with model)", n, arm),
    tg,
    if (n == "1") tg else NA_character_,
    if (n == "2") tg else NA_character_)
}

## pull the three columns of one block, in the row order of `id`
fetch_block <- function(blk) {
  tg <- ifelse(id$track == "ASC", blk$tag_asc, blk$tag_ssc)
  m  <- val[.(id$gene, tg), on = .(gene, contrast_tag)]
  na_tag <- is.na(tg)
  list(log2FC = ifelse(na_tag, NA_real_, m$log2FC),
       padj   = ifelse(na_tag, NA_real_, m$padj),
       lab    = ifelse(na_tag, "not_applicable", m$DEsupport))
}

flat <- data.table(track = id$track, gene = id$gene, gene_name = id$gene_name,
                   Z_score = id$Z_score, P_value = id$P_value)
hdr1 <- rep("", 5L)
hdr2 <- c("track", "gene", "gene_name", "Z_score", "P_value")

for (k in seq_along(blocks)) {
  blk <- blocks[[k]]
  v   <- fetch_block(blk)
  flat[[paste0(blk$flat, "_log2FC")]]    <- v$log2FC
  flat[[paste0(blk$flat, "_padj")]]      <- v$padj
  flat[[paste0(blk$flat, "_DEsupport")]] <- v$lab
  hdr1 <- c(hdr1, blk$header, "", "")
  hdr2 <- c(hdr2, "log2FC", "padj", "DE_discordant")

  if (k == ANY_C_AT) {                # Any C Support, a single column
    c_lab <- grep("^C[1-7]_MS(IH|S)_DEsupport$", names(flat), value = TRUE)
    flat[["Any_C_Support"]] <-
      rowSums(as.matrix(flat[, ..c_lab]) == "Up" |
              as.matrix(flat[, ..c_lab]) == "Down")
    hdr1 <- c(hdr1, "Any C Support")
    hdr2 <- c(hdr2, "")
  }
}
stopifnot(length(hdr1) == ncol(flat), length(hdr2) == ncol(flat))
## 5 ids + 14 C blocks + 1 roll-up + 7 E blocks + 4 L blocks, three columns each
stopifnot(ncol(flat) == 5L + 3L * 14L + 1L + 3L * 7L + 3L * 4L)

fwrite(flat, file.path(DIR_TAB,
                       "ASC_SSC_genes_crossceltype_support_matrix_flat.csv"),
       na = "")

## grouped-header version: the two header rows, then the body with no names
csv_escape <- function(x) {
  x <- as.character(x)
  needs <- grepl('[",\n]', x)
  x[needs] <- paste0('"', gsub('"', '""', x[needs]), '"')
  x
}
MATRIX_CSV <- file.path(DIR_TAB, "ASC_SSC_genes_crossceltype_support_matrix.csv")
writeLines(c(paste(csv_escape(hdr1), collapse = ","),
             paste(csv_escape(hdr2), collapse = ",")), MATRIX_CSV)
fwrite(flat, MATRIX_CSV, append = TRUE, col.names = FALSE, na = "")

## the grouped-header file must be the two headers plus one row per gene
stopifnot(length(readLines(MATRIX_CSV)) == 2L + nrow(flat))

## 6. per-gene summary --------------------------------------------------------
sup_lv <- c("Up", "Down")
gsum <- long[, {
  Cs <- .SD[block == "C"]; Es <- .SD[block == "E"]; Ls <- .SD[block == "L"]
  list(
    Z_score = Z_score[1], P_value = P_value[1],
    direction = ifelse(Z_score[1] > 0, "risk-increasing", "risk-decreasing"),

    C_tested          = sum(Cs$DEsupport != "not_tested"),
    C_supported       = sum(Cs$DEsupport %chin% sup_lv),
    C_discordant      = sum(Cs$DEsupport == "DE_discordant"),
    C_celltypes       = paste(sort(unique(Cs[DEsupport %chin% sup_lv,
                                             contrast_cell_type])), collapse = ","),
    C_arms            = paste(sort(unique(Cs[DEsupport %chin% sup_lv, arm])),
                              collapse = ","),

    E_tested          = sum(Es$DEsupport != "not_tested"),
    E_supported       = sum(Es$DEsupport %chin% sup_lv),
    E_discordant      = sum(Es$DEsupport == "DE_discordant"),
    E_celltypes       = paste(sort(unique(Es[DEsupport %chin% sup_lv,
                                             contrast_cell_type])), collapse = ","),

    L_tested          = sum(Ls$DEsupport != "not_tested"),
    L_supported       = sum(Ls$DEsupport %chin% sup_lv),
    L_arms            = paste(sort(unique(Ls[DEsupport %chin% sup_lv, arm])),
                              collapse = ","),

    n_supported       = sum(DEsupport %chin% sup_lv),
    n_tested          = sum(DEsupport != "not_tested")
  )
}, by = .(track, gene, gene_name)]

gsum[, blocks_supported := paste(c("C", "E", "L")[c(C_supported, E_supported,
                                                    L_supported) > 0],
                                 collapse = ","),
     by = .(track, gene)]
gsum[, any_support := n_supported > 0]
## setorder takes column names, not expressions, so |Z| needs its own column
gsum[, absZ := abs(Z_score)]
setorder(gsum, track, -n_supported, -absZ)
gsum[, absZ := NULL]
fwrite(gsum, file.path(DIR_TAB, "ASC_SSC_genes_crossceltype_support_summary.csv"),
       na = "")

## 7. per-contrast counts -----------------------------------------------------
lv <- c("Up", "Down", "DE_discordant", "NS", "not_tested")
cnt <- long[, {
  t <- table(factor(DEsupport, levels = lv))
  c(as.list(t), list(n_genes = .N, supported = t[["Up"]] + t[["Down"]]))
}, by = .(track, block, contrast_tag, contrast_cell_type, arm)]
setorder(cnt, track, block, contrast_tag)
fwrite(cnt, file.path(DIR_TAB, "ASC_SSC_genes_crossceltype_support_counts.csv"),
       na = "")

## 8. console summary ---------------------------------------------------------
cat("\n== support per contrast ==\n")
print(cnt[, .(track, contrast_tag, cell_type = contrast_cell_type, arm,
              Up, Down, DE_discordant, NS, not_tested, supported)])

cat("\n== genes with any cross-cell-type support ==\n")
print(gsum[, .(n_genes = .N, with_support = sum(any_support),
               C_only = sum(C_supported > 0 & E_supported == 0 & L_supported == 0),
               E_any  = sum(E_supported > 0),
               L_any  = sum(L_supported > 0)), by = track])

cat("\n== per-gene detail, genes with support ==\n")
print(gsum[any_support == TRUE,
           .(track, gene_name, Z = round(Z_score, 2), direction,
             C_supported, C_celltypes, E_supported, E_celltypes,
             L_supported, L_arms, blocks_supported)])

cat("\n== genes with no support anywhere in the panel ==\n")
print(gsum[any_support == FALSE, .(track, gene_name, Z = round(Z_score, 2),
                                   n_tested, C_discordant, E_discordant)])

cat("\nwrote 5 tables to ", DIR_TAB, "\n", sep = "")
