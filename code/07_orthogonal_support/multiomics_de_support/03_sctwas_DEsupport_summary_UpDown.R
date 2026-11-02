## 03_sctwas_DEsupport_summary_UpDown.R ---------------------------------------
## Summary of how many scTWAS genes each cell-type model has supported by
## scRNA-seq differential expression, split Up / Down, for the eight C / E / L
## contrasts.
##
## Input
##   table/DE_key_results/a_scTWAS_ADSER_NL_collected_results_scRNAseq_DEsupport.csv
##     499 gene x cell-type pairs (the 546-pair table with the 47 ASC/SSC-model
##     pairs removed), 36 columns, eight DEsupport blocks. Byte-identical to
##     table/scTWAS_ADSER_NL_collected_results_scRNAseq_DEsupport.csv.
##
## Output
##   table/DE_key_results/scTWAS_DEsupport_summary_UpDown.csv
##     the requested layout: two header rows, Models x Cell types down the side,
##     eight contrasts across the top, each split Up / Down
##   table/DE_key_results/scTWAS_DEsupport_summary_UpDown_long.csv
##     one row per model x contrast x direction, with the model's denominator
##   figure/Fig_scTWAS_DEsupport_summary_UpDown.{png,pdf}
##   figure/Fig_scTWAS_DEsupport_summary_UpDown_sourcedata.csv
##
## What a cell counts
##   The number of scTWAS gene x cell-type pairs of that model whose DEsupport
##   label for that contrast is Up (or Down). DEsupport was assigned upstream in
##   01_sctwas_scRNAseq_DEsupport_C_E_L.R under the project rule padj < 0.05 AND
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
##   is in 02_sctwas_ASC_SSC_crossceltype_support.R.
## -----------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

PROJ <- "W:/VUMC_projects/scTWAS/VUMC_HTAN_scTWAS_genes_DE_support"
if (!dir.exists(file.path(PROJ, "table"))) PROJ <- "."
DIR_KEY <- file.path(PROJ, "table", "DE_key_results")
DIR_FIG <- file.path(PROJ, "figure")
dir.create(DIR_FIG, recursive = TRUE, showWarnings = FALSE)

KEY_CSV <- file.path(DIR_KEY,
  "a_scTWAS_ADSER_NL_collected_results_scRNAseq_DEsupport.csv")

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
obs  <- dcast(m[DEsupport %chin% c("Up", "Down")],
              Dataset + CellType + block ~ DEsupport,
              value.var = "gene", fun.aggregate = length)
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

## 5. figure ------------------------------------------------------------------
## Up and Down are opposing directions, so each gets its own sequential ramp in
## the project's palette - red for up, blue for down - rather than sharing one
## neutral ramp. The two ramps carry identical limits, so a 7 in the Up half and
## a 7 in the Down half are the same depth of colour and the halves stay
## comparable. Each ramp keeps its own colourbar beneath the half it applies to.
## Every cell prints its value, so the figure is exact rather than indicative.
##
## Two fill scales cannot coexist in one ggplot, so the halves are built
## separately and joined with patchwork.
suppressPackageStartupMessages(library(patchwork))

UP_COL   <- "#B2182B"   # up, project convention
DOWN_COL <- "#2166AC"   # down, project convention
UP_LOW   <- "#FFF5F0"
DOWN_LOW <- "#F7FBFF"

fg <- copy(lng)
fg[, Models := factor(Models, levels = MODEL_LABEL)]
fg[, contrast_short := factor(contrast_short, levels = BLOCKS$short)]
fg[, ct_label := sprintf("%s (%d)", CellType, n_pairs)]
lv <- fg[, .(ct_label = unique(ct_label)), by = .(Models, CellType)]
setorder(lv, Models, CellType)
fg[, ct_label := factor(ct_label, levels = rev(unique(lv$ct_label)))]

## shared limits across both halves, and one text-contrast threshold that
## therefore applies to both
CNT_MAX <- max(fg$n_supported)
TXT_CUT <- ceiling(0.55 * CNT_MAX)
fg[, txt_col := fifelse(n_supported >= TXT_CUT, "white", "black")]

half_plot <- function(dir_label, low, high, show_y) {
  ggplot(fg[direction == dir_label],
         aes(x = contrast_short, y = ct_label, fill = n_supported)) +
    geom_tile(colour = "white", linewidth = 0.7) +
    geom_text(aes(label = n_supported, colour = txt_col), size = 2.5) +
    facet_grid(Models ~ ., scales = "free_y", space = "free_y") +
    scale_fill_gradient(low = low, high = high, limits = c(0, CNT_MAX),
                        breaks = c(0, 5, 10, 15, 20),
                        name = sprintf("%s: scTWAS pairs supported", dir_label),
                        guide = guide_colourbar(title.position = "top",
                                                barheight = unit(2.4, "mm"),
                                                barwidth  = unit(28, "mm"))) +
    scale_colour_identity() +
    scale_x_discrete(position = "top") +
    labs(x = NULL,
         y = if (show_y) "scTWAS cell-type model  (pairs tested)" else NULL) +
    theme_classic(base_size = 8) +
    theme(axis.text.x  = element_text(size = 6, angle = 40, hjust = 0,
                                      colour = "black"),
          axis.text.y  = if (show_y) element_text(size = 6, colour = "black")
                         else element_blank(),
          axis.title.y = element_text(size = 8),
          axis.ticks   = element_blank(),
          axis.line    = element_blank(),
          strip.background = element_blank(),
          strip.text.y = element_text(size = 8, angle = -90),
          legend.position  = "bottom",
          legend.direction = "horizontal",
          legend.title = element_text(size = 7),
          legend.text  = element_text(size = 6),
          panel.spacing = unit(2, "mm"),
          plot.margin = margin(2, 4, 2, 2))
}

p_up   <- half_plot("Up",   UP_LOW,   UP_COL,   show_y = TRUE)
p_down <- half_plot("Down", DOWN_LOW, DOWN_COL, show_y = FALSE)

p <- (p_up | p_down) +
  plot_annotation(
    title = paste0("Every MSS contrast supports more scTWAS genes than its ",
                   "MSI-H counterpart, by 1.6 to 2.5 fold"),
    subtitle = paste0("scRNA-seq DE support for 499 scTWAS gene \u00d7 cell-type ",
                      "pairs across eight contrasts; red = up in the later ",
                      "state, blue = down, and cells print the pair count"),
    caption = paste0(
      "Supported = padj < 0.05, |log2FC| > 0.5 and sign(log2FC) matching the scTWAS Z. Up and Down are both supported and record the direction of\n",
      "the expression change, not agreement with the risk Z; discordant, non-significant and untested pairs appear in neither half, so rows do not sum\n",
      "to the denominator in brackets. Both colour ramps share the same 0-", CNT_MAX,
      " limits, so equal counts read as equal depth across the two halves. 'Normal (Matched cell\n",
      "types)' contrasts are cell-type matched (ABS = C1/E1 ... TAC = C7/E7); the ASC and SSC contrasts are gene level and apply to every model. The\n",
      "ADSER ASC and SSC models have no matched normal contrast and are analysed separately."),
    theme = theme(
      plot.title    = element_text(size = 8, hjust = 0),
      plot.subtitle = element_text(size = 7, hjust = 0, colour = "grey30"),
      plot.caption  = element_text(size = 5.5, hjust = 0, colour = "grey30")))

FIG_STEM <- file.path(DIR_FIG, "Fig_scTWAS_DEsupport_summary_UpDown")
fwrite(fg[, .(Models, Dataset, CellType, n_pairs, n_genes, block,
              contrast, contrast_short, direction, n_supported)],
       paste0(FIG_STEM, "_sourcedata.csv"))
ggsave(paste0(FIG_STEM, ".png"), plot = p, width = 9.0, height = 5.8,
       units = "in", dpi = 300)
ggsave(paste0(FIG_STEM, ".pdf"), plot = p, width = 9.0, height = 5.8,
       units = "in", useDingbats = FALSE)

## 6. console summary --------------------------------------------------------
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
cat("\nwrote ", SUM_CSV, "\n  and ", FIG_STEM, ".{png,pdf}\n", sep = "")
