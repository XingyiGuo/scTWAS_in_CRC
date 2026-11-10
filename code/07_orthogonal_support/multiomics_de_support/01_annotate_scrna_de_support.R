## 01_annotate_scrna_de_support.R ---------------------------------------
## scRNA-seq differential-expression support for the scTWAS gene x cell-type
## pairs, using the C / E / L contrast set (ASC n=75, SSC n=46).
##
## Input
##   data/scTWAS_ADSER_NL_collected_results.csv      546 gene x cell-type pairs
##   data/de_C_E_L_contrasts_ASC75_SSC46.csv         pseudobulk DESeq2 results
##   data/resutls_example.xlsx                       target column layout
##
## Output
##   table/scTWAS_ADSER_NL_collected_results_scRNAseq_DEsupport.csv   flat names
##   table/scTWAS_ADSER_NL_collected_results_scRNAseq_DEsupport.xlsx  grouped header
##
## Eight DE blocks are attached to every scTWAS row, each contributing
## Log2FC / padj / DEsupport:
##   MSI-H vs matched normal cell type      C1-C7, MSI-H arm
##   MSS   vs matched normal cell type      C1-C7, MSS arm
##   ASC   vs matched normal cell type      E1-E7, ASC arm
##   SSC   vs matched normal cell type      E1-E7, SSC arm
##   MSI-H vs ASC / MSS vs ASC              L1_MSIH, L1_MSS      (gene level)
##   MSI-H vs SSC / MSS vs SSC              L2_MSIH, L2_MSS      (gene level)
##
## The C and E blocks are cell-type matched: a scTWAS row for cell type X is
## joined to the contrast whose cell_type is X (C1/E1 = ABS, C2/E2 = GOB,
## C3/E3 = CT, C4/E4 = EE, C5/E5 = STM, C6/E6 = TUF, C7/E7 = TAC).  The ADSER
## rows for cell types ASC and SSC have no matched normal-epithelium contrast,
## so their C and E blocks are necessarily not_tested.  The L blocks compare
## carcinoma against the polyp compartment itself and are therefore gene level:
## every scTWAS row carrying the same gene receives the same L values.
##
## DEsupport labelling (house rule, reproduces data/resutls_example.xlsx exactly)
##   not_tested     gene absent from that contrast's DE table (no log2FC)
##   NS             padj is NA, or padj >= 0.05, or |log2FC| <= 0.5
##   Up             padj < 0.05 & |log2FC| > 0.5 & log2FC > 0 & sign matches Z
##   Down           padj < 0.05 & |log2FC| > 0.5 & log2FC < 0 & sign matches Z
##   DE_discordant  padj < 0.05 & |log2FC| > 0.5 & sign opposes the scTWAS Z
## -----------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(data.table)
})

## paths -----------------------------------------------------------------------
PROJ <- "W:/VUMC_projects/scTWAS/VUMC_HTAN_scTWAS_genes_DE_support"
if (!dir.exists(file.path(PROJ, "data"))) PROJ <- "."   # local copies fallback

DIR_DATA <- file.path(PROJ, "data")
DIR_TAB  <- file.path(PROJ, "table")
dir.create(DIR_TAB, recursive = TRUE, showWarnings = FALSE)

TWAS_CSV <- file.path(DIR_DATA, "scTWAS_ADSER_NL_collected_results.csv")
DE_CSV   <- file.path(DIR_DATA, "de_C_E_L_contrasts_ASC75_SSC46.csv")
DE_PARQ  <- file.path(DIR_DATA, "de_C_E_L_contrasts_ASC75_SSC46.parquet")

OUT_CSV  <- file.path(DIR_TAB, "scTWAS_ADSER_NL_collected_results_scRNAseq_DEsupport.csv")
OUT_XLSX <- file.path(DIR_TAB, "scTWAS_ADSER_NL_collected_results_scRNAseq_DEsupport.xlsx")

## thresholds ------------------------------------------------------------------
Q_CUT   <- 0.05       # BH q-value ceiling
LFC_CUT <- 0.5        # |log2 fold change| floor, strict

## 1. scTWAS pairs -------------------------------------------------------------
twas <- fread(TWAS_CSV)
stopifnot(all(c("Dataset", "CellType", "gene", "gene_name",
                "Sig_resource", "Z_score", "P_value") %in% names(twas)))
twas[, row_id := .I]
message(sprintf("scTWAS: %d pairs, %d unique symbols, %d cell types",
                nrow(twas), uniqueN(twas$gene_name), uniqueN(twas$CellType)))

## 2. scRNA-seq DE results -----------------------------------------------------
## Only the five columns used downstream are read; the full CSV is ~300 MB.
de_cols <- c("gene", "cell_type", "contrast_tag", "stage", "group1",
             "log2FC", "padj")
if (requireNamespace("arrow", quietly = TRUE) && file.exists(DE_PARQ)) {
  de <- as.data.table(arrow::read_parquet(DE_PARQ,
                                          col_select = dplyr::all_of(de_cols)))
} else {
  de <- fread(DE_CSV, select = de_cols)
}
de[, log2FC := as.numeric(log2FC)]
de[, padj   := as.numeric(padj)]
message(sprintf("DE: %d rows, %d contrast tags, %d cell types",
                nrow(de), uniqueN(de$contrast_tag), uniqueN(de$cell_type)))

## 3. DEsupport label ----------------------------------------------------------
## Vectorised; `z` is the scTWAS Z score of the same row.
label_support <- function(log2FC, padj, z) {
  out <- rep("NS", length(log2FC))
  sig <- !is.na(padj) & padj < Q_CUT & !is.na(log2FC) & abs(log2FC) > LFC_CUT
  conc <- sign(log2FC) == sign(z)
  out[sig &  conc & log2FC > 0] <- "Up"
  out[sig &  conc & log2FC < 0] <- "Down"
  out[sig & !conc]              <- "DE_discordant"
  out[is.na(log2FC)]            <- "not_tested"   # gene missing from contrast
  out
}

## 4. block definitions --------------------------------------------------------
## `match_ct = TRUE`  -> join on gene symbol AND cell type (C and E blocks)
## `match_ct = FALSE` -> join on gene symbol only          (L blocks)
blocks <- list(
  list(key = "C_MSIH",  match_ct = TRUE,
       filter = quote(stage == "normal_vs_carcinoma"       & group1 == "MSI-H"),
       group  = "MSI-H vs normal cell type*"),
  list(key = "C_MSS",   match_ct = TRUE,
       filter = quote(stage == "normal_vs_carcinoma"       & group1 == "MSS"),
       group  = "MSS vs normal cell type*"),
  list(key = "E_ASC",   match_ct = TRUE,
       filter = quote(stage == "normal_vs_polyp_celltype"  & group1 == "ASC"),
       group  = "ASC vs normal cell type*"),
  list(key = "E_SSC",   match_ct = TRUE,
       filter = quote(stage == "normal_vs_polyp_celltype"  & group1 == "SSC"),
       group  = "SSC vs normal cell type*"),
  list(key = "L1_MSIH", match_ct = FALSE,
       filter = quote(contrast_tag == "L1_MSIH"), group = "MSI-H vs ASC"),
  list(key = "L1_MSS",  match_ct = FALSE,
       filter = quote(contrast_tag == "L1_MSS"),  group = "MSS vs ASC"),
  list(key = "L2_MSIH", match_ct = FALSE,
       filter = quote(contrast_tag == "L2_MSIH"), group = "MSI-H vs SSC"),
  list(key = "L2_MSS",  match_ct = FALSE,
       filter = quote(contrast_tag == "L2_MSS"),  group = "MSS vs SSC")
)

res <- copy(twas)
for (b in blocks) {
  sub <- de[eval(b$filter)]
  if (nrow(sub) == 0L) stop("no DE rows selected for block ", b$key)
  if (b$match_ct) {
    stopifnot(!anyDuplicated(sub, by = c("gene", "cell_type")))
    sub <- sub[, .(gene, cell_type, .lfc = log2FC, .padj = padj)]
    m <- sub[.(res$gene_name, res$CellType), on = .(gene, cell_type)]
  } else {
    stopifnot(!anyDuplicated(sub, by = "gene"))
    sub <- sub[, .(gene, .lfc = log2FC, .padj = padj)]
    m <- sub[.(res$gene_name), on = .(gene)]
  }
  res[[paste0(b$key, "_log2FC")]] <- round(m$.lfc, 4)
  res[[paste0(b$key, "_padj")]]   <- m$.padj
  res[[paste0(b$key, "_DEsupport")]] <-
    label_support(m$.lfc, m$.padj, res$Z_score)
  tags <- sort(unique(de[eval(b$filter)]$contrast_tag))
  message(sprintf("  %-8s %2d tag(s) [%s] -> %d of %d rows matched", b$key,
                  length(tags), paste(tags, collapse = ", "),
                  sum(!is.na(m$.lfc)), nrow(res)))
}

## 5. roll-up columns ----------------------------------------------------------
## TRUE  : at least one arm of the pair is supported (Up or Down)
## FALSE : both arms tested, neither supported
## NA    : neither arm testable for this gene x cell type
either <- function(a, b) {
  sa <- res[[paste0(a, "_DEsupport")]]; sb <- res[[paste0(b, "_DEsupport")]]
  ifelse(sa == "not_tested" & sb == "not_tested", NA,
         sa %in% c("Up", "Down") | sb %in% c("Up", "Down"))
}
res[, either_C   := either("C_MSIH",  "C_MSS")]
res[, either_E   := either("E_ASC",   "E_SSC")]
res[, either_L1  := either("L1_MSIH", "L1_MSS")]
res[, either_L2  := either("L2_MSIH", "L2_MSS")]

sup_cols <- paste0(vapply(blocks, `[[`, "", "key"), "_DEsupport")
res[, any_support := Reduce(`|`, lapply(sup_cols,
      function(cc) res[[cc]] %in% c("Up", "Down")))]
res[, any_support := ifelse(
      Reduce(`&`, lapply(sup_cols, function(cc) res[[cc]] == "not_tested")),
      NA, any_support)]

## 6. assemble in the resutls_example.xlsx column order ------------------------
id_cols <- c("Dataset", "CellType", "gene", "gene_name",
             "Sig_resource", "Z_score", "P_value")
ord <- id_cols
for (k in c("C_MSIH", "C_MSS")) ord <- c(ord, paste0(k, c("_log2FC", "_padj", "_DEsupport")))
ord <- c(ord, "either_C")
for (k in c("E_ASC", "E_SSC"))  ord <- c(ord, paste0(k, c("_log2FC", "_padj", "_DEsupport")))
ord <- c(ord, "either_E")
for (k in c("L1_MSIH", "L1_MSS")) ord <- c(ord, paste0(k, c("_log2FC", "_padj", "_DEsupport")))
ord <- c(ord, "either_L1")
for (k in c("L2_MSIH", "L2_MSS")) ord <- c(ord, paste0(k, c("_log2FC", "_padj", "_DEsupport")))
ord <- c(ord, "either_L2", "any_support")
out <- res[, ..ord]
stopifnot(nrow(out) == nrow(twas), ncol(out) == 36L)

fwrite(out, OUT_CSV)
message("wrote ", OUT_CSV)

## 7. same table with the two-row grouped header ------------------------------
## Row 1 blank, row 2 block labels, row 3 column names, data from row 4 -
## the layout of data/resutls_example.xlsx.
if (requireNamespace("openxlsx", quietly = TRUE)) {
  grp <- rep("", 36L)
  grp[8]  <- "MSI-H vs normal cell type*"
  grp[11] <- "MSS vs normal cell type*"
  grp[14] <- "Either MSI-H/MSS vs normal support"
  grp[15] <- "ASC vs normal cell type*"
  grp[18] <- "SSC vs normal cell type*"
  grp[21] <- "Either ASC/SSC vs normal support"
  grp[22] <- "MSI-H vs ASC"
  grp[25] <- "MSS vs ASC"
  grp[28] <- "Either MSI-H/MSS vs ASC support"
  grp[29] <- "MSI-H vs SSC"
  grp[32] <- "MSS vs SSC"
  grp[35] <- "Either MSI-H/MSS vs SSC support"
  grp[36] <- "Any support"
  nm <- c(id_cols, rep(c("Log2FC", "padj", "DEsupport"), 2), "True/False",
          rep(c("Log2FC", "padj", "DEsupport"), 2), "True/False",
          rep(c("Log2FC", "padj", "DEsupport"), 2), "True/False",
          rep(c("Log2FC", "padj", "DEsupport"), 2), "True/False", "")

  wb <- openxlsx::createWorkbook()
  openxlsx::addWorksheet(wb, "Sheet1")
  openxlsx::writeData(wb, 1, as.data.table(as.list(grp)), startRow = 2,
                      colNames = FALSE)
  openxlsx::writeData(wb, 1, as.data.table(as.list(nm)),  startRow = 3,
                      colNames = FALSE)
  openxlsx::writeData(wb, 1, out, startRow = 4, colNames = FALSE)
  hdr <- openxlsx::createStyle(textDecoration = "bold", halign = "center",
                               valign = "center", wrapText = TRUE)
  openxlsx::addStyle(wb, 1, hdr, rows = 2:3, cols = 1:36, gridExpand = TRUE)
  for (st in c(8, 11, 15, 18, 22, 25, 29, 32))
    openxlsx::mergeCells(wb, 1, cols = st:(st + 2), rows = 2)
  openxlsx::setColWidths(wb, 1, cols = 1:36, widths = "auto")
  openxlsx::freezePane(wb, 1, firstActiveRow = 4, firstActiveCol = 5)
  openxlsx::saveWorkbook(wb, OUT_XLSX, overwrite = TRUE)
  message("wrote ", OUT_XLSX)
} else {
  message("openxlsx not installed - CSV written, xlsx skipped")
}

## 8. console summary ----------------------------------------------------------
lv <- c("Up", "Down", "DE_discordant", "NS", "not_tested")
tab <- sapply(sup_cols, function(cc) table(factor(out[[cc]], levels = lv)))
print(t(tab))
cat("\nsupported (Up/Down) pairs per block:\n")
print(colSums(tab[c("Up", "Down"), , drop = FALSE]))
cat("\ngenes with support in at least one block:",
    uniqueN(out[any_support == TRUE, gene_name]), "of",
    uniqueN(out$gene_name), "\n")
cat("pairs with support in at least one block:",
    sum(out$any_support, na.rm = TRUE), "of", nrow(out), "\n")
cat("pairs never testable (all eight blocks not_tested):",
    sum(is.na(out$any_support)), "\n")
