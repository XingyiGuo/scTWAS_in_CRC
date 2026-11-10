## 03_summarize_scrna_support_counts.R ---------------------------------------
## Summary of how many scTWAS genes each cell-type model has supported by
## scRNA-seq differential expression, split Up / Down, for the eight C / E / L
## contrasts.
##
## Input
##   table/DE_key_results/a_scTWAS_ADSER_NL_collected_results_scRNAseq_DEsupport.csv
##   or table/scTWAS_ADSER_NL_collected_results_scRNAseq_DEsupport.csv from step 01.
##   Both contain the matched-cell-type pairs and eight DEsupport blocks.
##
## Output
##   table/DE_key_results/scTWAS_DEsupport_summary_UpDown.csv
##     the requested layout: two header rows, Models x Cell types down the side,
##     eight contrasts across the top, each split Up / Down
##   table/DE_key_results/scTWAS_DEsupport_summary_UpDown_long.csv
##     one row per model x contrast x direction, with the model's denominator
##   The grouped bars are rendered separately by 04_de_support_grouped_bars.R.
##
## What a cell counts
##   The number of scTWAS gene x cell-type pairs of that model whose DEsupport
##   label for that contrast is Up (or Down). DEsupport was assigned upstream in
##   01_annotate_scrna_de_support.R under the project rule padj < 0.05 AND
##   |log2FC| > 0.5 AND sign(log2FC) == sign(Z); Up and Down therefore both mean
##   "supported", and record the direction of the expression change rather than
##   agreement or disagreement with the risk Z. DE_discordant, NS and not_tested
##   pairs are not counted in either column, so a row does not sum to its
##   denominator.
##
## Denominators
##   Pairs per model, not genes: a gene nominated in several cell-type models is
##   counted once in each. The C and E blocks are cell-type matched upstream, so
##   for a model of cell type X the "Normal (Matched cell types)" contrasts are
##   that cell type's own contrast (ABS = C1/E1 ... TAC = C7/E7). The L blocks
##   compare carcinoma against the polyp compartment and are gene level, so they
##   apply to every model regardless of its cell type.
##
## Scope note
##   The ADSER ASC and SSC models are absent from the input by design - they
##   have no matched normal-epithelium contrast. Their cross-cell-type treatment
##   is in 02_analyze_asc_ssc_cross_celltype_support.R.
## -----------------------------------------------------------------------------

suppressPackageStartupMessages(library(data.table))

default_root <- "W:/VUMC_projects/scTWAS/VUMC_HTAN_scTWAS_genes_DE_support"
if (!dir.exists(file.path(default_root, "table"))) default_root <- "."
PROJ <- Sys.getenv("SC_TWAS_SUPPORT_ROOT", default_root)
DIR_KEY <- file.path(PROJ, "table", "DE_key_results")

KEY_CSV <- file.path(DIR_KEY,
  "a_scTWAS_ADSER_NL_collected_results_scRNAseq_DEsupport.csv")
if (!file.exists(KEY_CSV)) {
  KEY_CSV <- file.path(PROJ, "table", "scTWAS_ADSER_NL_collected_results_scRNAseq_DEsupport.csv")
}
dir.create(DIR_KEY, recursive = TRUE, showWarnings = FALSE)

## contrast blocks, in the order of the requested layout -----------------------
## block       = the DEsupport column prefix in the input
## header      = group header for the output table, verbatim as specified
## short       = compact label for the figure axis
BLOCKS <- data.table(
  block  = c("C_MSIH", "C_MSS", "L1_MSIH", "L1_MSS",
             "L2_MSIH", "L2_MSS", "E_ASC", "E_SSC"),
  header = c("Normal (Matched cell types) -> MISH",
             "Normal(Matched cell types) -> MSS",
             "ASC -> MISH", "ASC -> MSS",
             "SSC -> MISH", "SSC -> MSS",
             "Normal (Matched cell types) -> ASC",
             "Normal (Matched cell types)-> SSC"),
  short  = c("Normal \u2192 MSI-H", "Normal \u2192 MSS",
             "ASC \u2192 MSI-H", "ASC \u2192 MSS",
             "SSC \u2192 MSI-H", "SSC \u2192 MSS",
             "Normal \u2192 ASC", "Normal \u2192 SSC"))

CT_ORDER    <- c("ABS", "GOB", "CT", "EE", "STM", "TUF", "TAC")
MODEL_LABEL <- c(NL = "NL models", ADSER = "Polyps models")

## 1. input -------------------------------------------------------------------
d <- fread(KEY_CSV)
supcols <- paste0(BLOCKS$block, "_DEsupport")
stopifnot(all(c("Dataset", "CellType", "gene") %in% names(d)),
          all(supcols %in% names(d)))
stopifnot(all(unique(unlist(d[, ..supcols])) %chin%
                c("Up", "Down", "DE_discordant", "NS", "not_tested")))
stopifnot(setequal(unique(d$Dataset), c("NL", "ADSER")),
          setequal(unique(d$CellType), CT_ORDER))
message(sprintf("input: %d pairs, %d models", nrow(d),
                uniqueN(d[, .(Dataset, CellType)])))

## 2. counts ------------------------------------------------------------------
m <- melt(d[, c("Dataset", "CellType", "gene", supcols), with = FALSE],
          id.vars = c("Dataset", "CellType", "gene"),
          variable.name = "block", value.name = "DEsupport")
m[, block := sub("_DEsupport$", "", as.character(block))]

## complete grid first, so a model x contrast with no supported pair is a
## printed zero rather than a missing row
grid <- CJ(Dataset = c("NL", "ADSER"), CellType = CT_ORDER,
           block = BLOCKS$block, unique = TRUE)
obs  <- m[, .(Up = sum(DEsupport == "Up"), Down = sum(DEsupport == "Down")),
          by = .(Dataset, CellType, block)]
cnt  <- merge(grid, obs, by = c("Dataset", "CellType", "block"), all.x = TRUE)
cnt[is.na(Up), Up := 0L][is.na(Down), Down := 0L]

## denominator: pairs tested in that model
den <- d[, .(n_pairs = .N, n_genes = uniqueN(gene)), by = .(Dataset, CellType)]
cnt <- merge(cnt, den, by = c("Dataset", "CellType"))

cnt[, Dataset  := factor(Dataset, levels = c("NL", "ADSER"))]
cnt[, CellType := factor(CellType, levels = CT_ORDER)]
cnt[, block    := factor(block, levels = BLOCKS$block)]
setorder(cnt, Dataset, CellType, block)

## the two counted labels are a subset of the model's pairs
stopifnot(cnt[, all(Up + Down <= n_pairs)])
## and every pair of every model is accounted for by some label
stopifnot(nrow(m) == nrow(d) * nrow(BLOCKS))

## 3. long output -------------------------------------------------------------
lng <- melt(cnt, id.vars = c("Dataset", "CellType", "block", "n_pairs", "n_genes"),
            measure.vars = c("Up", "Down"),
            variable.name = "direction", value.name = "n_supported")
lng <- merge(lng, BLOCKS, by = "block", sort = FALSE)
lng[, Models := MODEL_LABEL[as.character(Dataset)]]
lng <- lng[, .(Models, Dataset, CellType, n_pairs, n_genes,
               block, contrast = header, contrast_short = short,
               direction, n_supported)]
setorder(lng, Dataset, CellType, block, direction)
fwrite(lng, file.path(DIR_KEY, "scTWAS_DEsupport_summary_UpDown_long.csv"))

## 4. requested layout -------------------------------------------------------
## Row 1 names each contrast and spans its two columns; row 2 carries Up/Down.
## The Models label is repeated on every row rather than merged, because a CSV
## cannot express a merged cell; merge the cells after pasting if wanted.
wide <- dcast(cnt, Dataset + CellType + n_pairs ~ block,
              value.var = c("Up", "Down"))
ordcols <- unlist(lapply(BLOCKS$block, function(b) paste0(c("Up_", "Down_"), b)))
out <- data.table(Models = MODEL_LABEL[as.character(wide$Dataset)],
                  `Cell types` = as.character(wide$CellType))
for (cc in ordcols) out[[cc]] <- wide[[cc]]

hdr1 <- c("Models", "Cell types",
          as.vector(rbind(BLOCKS$header, rep("", nrow(BLOCKS)))))
hdr2 <- c("", "", rep(c("Up", "Down"), nrow(BLOCKS)))
stopifnot(length(hdr1) == ncol(out), length(hdr2) == ncol(out),
          ncol(out) == 2L + 2L * nrow(BLOCKS))

csv_escape <- function(x) {
  x <- as.character(x)
  needs <- grepl('[",\n]', x)
  x[needs] <- paste0('"', gsub('"', '""', x[needs]), '"')
  x
}
SUM_CSV <- file.path(DIR_KEY, "scTWAS_DEsupport_summary_UpDown.csv")
writeLines(c(paste(csv_escape(hdr1), collapse = ","),
             paste(csv_escape(hdr2), collapse = ",")), SUM_CSV)
fwrite(out, SUM_CSV, append = TRUE, col.names = FALSE)
stopifnot(length(readLines(SUM_CSV)) == 2L + nrow(out))

## 5. console summary --------------------------------------------------------
tot <- cnt[, .(Up = sum(Up), Down = sum(Down)), by = block]
tot[, total := Up + Down][, pct_up := round(100 * Up / total, 1)]
cat("\n== support per contrast, both model sets ==\n"); print(tot)
cat("\n== MSS vs MSI-H within each pair ==\n")
print(data.table(
  pair  = c("Normal->carcinoma", "ASC->carcinoma", "SSC->carcinoma"),
  MSIH  = tot[block %chin% c("C_MSIH", "L1_MSIH", "L2_MSIH"), total],
  MSS   = tot[block %chin% c("C_MSS", "L1_MSS", "L2_MSS"), total]
)[, fold := round(MSS / MSIH, 2)][])
cat("\n== totals by model set ==\n")
print(cnt[, .(Up = sum(Up), Down = sum(Down)), by = Dataset])
cat("\nwrote ", SUM_CSV, " and the long count table\n", sep = "")
