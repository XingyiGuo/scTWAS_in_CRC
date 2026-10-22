# conda activate twasguo

library(stringr)
library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)

# ==============================================================================
# User settings
# ==============================================================================
root <- "/data/l2_bioinfo1/liq17/scTWAS/TWAS_HTAN_GSE266546_WellcomeSanger_res/SPrediXcan_peakSNP_association/"
OUTDIR <- "/data/l2_bioinfo1/liq17/scTWAS/Functional_res/"

test <- "NL"   # choose "ADSER" or "NL"
dis  <- "Colorectal Adenocarcinoma"

MEDIAN_CUTOFF <- -0.5
N_EXTRA <- 2

# Requested cell-type order in the combined figure.
# Any unexpected cell types are appended alphabetically at the end.
CELLTYPE_ORDER <- c("ASC", "SSC", "ABS", "GOB", "CT", "EE", "STM", "TAC", "TUF")

dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

# ==============================================================================
# 1. Read significant TWAS genes and retain gene_name + CellType
# ==============================================================================
if (test == "ADSER") {
  ctMWAS_genes <- fread(
    paste0(
      root,
      "/VUMC_HTAN_ADSER_EUR_Multi_sig/",
      "VUMC_HTAN_ADSER_9CTs_peakSNP.EUR_Multi.sig.csv"
    )
  )
} else if (test == "NL") {
  ctMWAS_genes <- fread(
    paste0(
      root,
      "/Cauchy_Meta_NL_WellcomeSanger_NL_EUR_Multi/",
      "HTAN_NL_WellcomeSanger_NL_7celltypes_peakSNP.EUR_Multi.CauchyMeta.sig.csv"
    )
  )
} else {
  stop("test must be either 'ADSER' or 'NL'.")
}

required_cols <- c("gene_name", "CellType")
missing_cols <- setdiff(required_cols, names(ctMWAS_genes))
if (length(missing_cols) > 0) {
  stop("Missing required TWAS columns: ", paste(missing_cols, collapse = ", "))
}

sig_df <- ctMWAS_genes %>%
  transmute(
    gene_name = str_trim(as.character(gene_name)),
    CellType  = str_trim(as.character(CellType))
  ) %>%
  filter(!is.na(gene_name), gene_name != "", !is.na(CellType), CellType != "") %>%
  distinct(gene_name, CellType)

print(dim(sig_df))
print(n_distinct(sig_df$gene_name))
print(
  sig_df %>%
    group_by(CellType) %>%
    summarise(n_genes = n_distinct(gene_name), .groups = "drop")
)

genes <- sort(unique(sig_df$gene_name))

# ==============================================================================
# 2. Query DepMap CRC dependency scores
# ==============================================================================
source("/data/l2_bioinfo1/liq17/LQ_references/GeneFunctioalSupports/FunctionalBasic.R")

dep_scores <- get_depmap_ceres_scores(genes, dis)
dep_scores <- as.data.frame(dep_scores, check.names = FALSE)

# ==============================================================================
# Helper: identify the gene-symbol column returned by get_depmap_ceres_scores()
# ==============================================================================
find_gene_column <- function(df) {
  # get_depmap_ceres_scores() currently returns the gene-symbol column as
  # "Genes" for the wide gene x CRC-cell-line matrix.
  candidates <- c(
    "gene_name",
    "Genes", "genes",
    "gene", "Gene", "GENE",
    "gene_symbol", "GeneSymbol",
    "symbol", "Symbol",
    "gene_symbol.x"
  )

  hit <- candidates[candidates %in% names(df)]
  if (length(hit) > 0) return(hit[1])

  # Case-insensitive fallback for minor naming differences.
  nm_lower <- tolower(names(df))
  candidate_lower <- tolower(candidates)
  idx <- match(candidate_lower, nm_lower, nomatch = 0)
  idx <- idx[idx > 0]

  if (length(idx) > 0) {
    return(names(df)[idx[1]])
  }

  stop(
    "Could not identify the gene column in dep_scores. Columns are: ",
    paste(names(df), collapse = ", ")
  )
}

gene_col <- find_gene_column(dep_scores)
if (gene_col != "gene_name") {
  names(dep_scores)[names(dep_scores) == gene_col] <- "gene_name"
}

dep_scores$gene_name <- str_trim(as.character(dep_scores$gene_name))
cat("Detected DepMap gene column and standardized it to: gene_name\n")
cat("DepMap dimensions:", nrow(dep_scores), "genes x", ncol(dep_scores) - 1, "CRC cell-line score columns\n")

# ==============================================================================
# 3. ADD CellType to the DepMap output
#
# Important: if a TWAS gene belongs to >1 cell type, its DepMap record is
# intentionally repeated once for each matching gene_name--CellType pair.
# ==============================================================================
if ("CellType" %in% names(dep_scores)) {
  dep_scores$CellType <- NULL
}

dep_scores_with_ct <- dep_scores %>%
  inner_join(sig_df, by = "gene_name") %>%
  relocate(CellType, .after = gene_name)

output_csv <- file.path(OUTDIR, paste0(test, "_DepMap_CRC_with_CellType.csv"))
fwrite(dep_scores_with_ct, output_csv)
cat("Saved CellType-annotated DepMap file:\n", output_csv, "\n")

# ==============================================================================
# 4. Convert DepMap results to long format for boxplots
#
# The helper supports two common return formats:
#   A) long: gene + cell-line column + score column
#   B) wide: one row per gene and one numeric column per CRC cell line
# ==============================================================================
first_existing <- function(x, candidates) {
  hit <- candidates[candidates %in% x]
  if (length(hit) == 0) return(NULL)
  hit[1]
}

make_depmap_long <- function(df) {
  cell_line_col <- first_existing(
    names(df),
    c(
      "cell_line", "CellLine", "cellline", "cell_line_name",
      "cell_line_display_name", "model_name", "ModelName", "model"
    )
  )

  score_col <- first_existing(
    names(df),
    c(
      "dependency_score", "dep_score", "DepScore", "CERES", "ceres",
      "gene_effect", "GeneEffect", "score", "Score"
    )
  )

  # Long format already present.
  if (!is.null(cell_line_col) && !is.null(score_col)) {
    out <- df %>%
      transmute(
        gene_name = as.character(gene_name),
        CellType = as.character(CellType),
        CellLine = as.character(.data[[cell_line_col]]),
        DepScore = suppressWarnings(as.numeric(.data[[score_col]]))
      )
    return(out)
  }

  # Otherwise assume a wide gene x cell-line matrix.
  # Numeric columns are interpreted as cell-line dependency scores.
  numeric_cols <- names(df)[vapply(df, is.numeric, logical(1))]

  # Exclude obvious numeric metadata columns if present.
  metadata_numeric <- grep(
    "(^|_)(id|entrez|start|end|chr|chrom|n|num|count)$|entrez|chromosome",
    numeric_cols,
    ignore.case = TRUE,
    value = TRUE
  )
  score_cols <- setdiff(numeric_cols, metadata_numeric)

  if (length(score_cols) == 0) {
    stop(
      "No numeric DepMap score columns were detected. Columns are: ",
      paste(names(df), collapse = ", ")
    )
  }

  df %>%
    select(gene_name, CellType, all_of(score_cols)) %>%
    pivot_longer(
      cols = all_of(score_cols),
      names_to = "CellLine",
      values_to = "DepScore"
    ) %>%
    mutate(
      gene_name = as.character(gene_name),
      CellType = as.character(CellType),
      CellLine = as.character(CellLine),
      DepScore = suppressWarnings(as.numeric(DepScore))
    )
}

dep_long <- make_depmap_long(dep_scores_with_ct) %>%
  filter(!is.na(DepScore), is.finite(DepScore))

if (nrow(dep_long) == 0) {
  stop("No usable DepMap dependency scores were available after reshaping.")
}

# ==============================================================================
# 5. Calculate median score across CRC cell lines for each TWAS gene/cell type
# ==============================================================================
median_df <- dep_long %>%
  group_by(CellType, gene_name) %>%
  summarise(
    median_score = median(DepScore, na.rm = TRUE),
    n_cell_lines = sum(!is.na(DepScore)),
    .groups = "drop"
  ) %>%
  filter(is.finite(median_score))

# Preserve requested cell-type order; append unexpected types at the end.
observed_ct <- unique(median_df$CellType)
ct_levels <- c(
  CELLTYPE_ORDER[CELLTYPE_ORDER %in% observed_ct],
  sort(setdiff(observed_ct, CELLTYPE_ORDER))
)

median_df <- median_df %>%
  mutate(CellType = factor(CellType, levels = ct_levels))

# ==============================================================================
# 6. Select genes for plotting
#
# For EACH cell type:
#   - keep every gene with median DepMap score < -0.5
#   - add the two genes with the smallest medians among genes > -0.5
# ==============================================================================
strong_genes <- median_df %>%
  filter(median_score < MEDIAN_CUTOFF) %>%
  mutate(Selection = paste0("Median < ", MEDIAN_CUTOFF))

extra_genes <- median_df %>%
  filter(median_score > MEDIAN_CUTOFF) %>%
  group_by(CellType) %>%
  arrange(median_score, gene_name, .by_group = TRUE) %>%
  slice_head(n = N_EXTRA) %>%
  ungroup() %>%
  mutate(Selection = paste0("Median > ", MEDIAN_CUTOFF))

selected_df <- bind_rows(strong_genes, extra_genes) %>%
  distinct(CellType, gene_name, .keep_all = TRUE) %>%
  arrange(CellType, median_score, gene_name)

if (nrow(selected_df) == 0) {
  stop("No genes were selected for plotting.")
}

# Composite ID keeps the same gene separate if it is significant in >1 cell type.
selected_df <- selected_df %>%
  mutate(
    CellType_chr = as.character(CellType),
    GeneCell = paste(CellType_chr, gene_name, sep = "___")
  )

gene_levels <- selected_df$GeneCell

dep_plot <- dep_long %>%
  mutate(CellType = factor(CellType, levels = ct_levels)) %>%
  inner_join(
    selected_df %>%
      select(CellType, gene_name, median_score, n_cell_lines, Selection, GeneCell),
    by = c("CellType", "gene_name")
  ) %>%
  mutate(GeneCell = factor(GeneCell, levels = gene_levels))

# Save median/selection summary used to build the figure.
summary_csv <- file.path(
  OUTDIR,
  paste0(test, "_DepMap_CRC_boxplot_gene_medians_selected.csv")
)
fwrite(as.data.frame(selected_df), summary_csv)
cat("Saved selected-gene median summary:\n", summary_csv, "\n")

# Save the actual long-format data used in the boxplot.
plotdata_csv <- file.path(
  OUTDIR,
  paste0(test, "_DepMap_CRC_boxplot_selected_long.csv")
)
fwrite(as.data.frame(dep_plot), plotdata_csv)
cat("Saved boxplot source data:\n", plotdata_csv, "\n")

# ==============================================================================
# 7. One combined figure
#
# Cell-type panels appear from left to right as:
# ASC -> SSC -> ABS -> GOB -> CT -> EE -> STM -> TAC -> TUF
# (only cell types present in the selected dataset are shown).
# Within each cell type, genes run from smallest median to largest median.
# ==============================================================================
p <- ggplot(
  dep_plot,
  aes(x = GeneCell, y = DepScore, fill = Selection)
) +
  geom_hline(
    yintercept = MEDIAN_CUTOFF,
    linetype = "dashed",
    linewidth = 0.45
  ) +
  geom_boxplot(
    width = 0.68,
    outlier.size = 0.8,
    outlier.alpha = 0.65,
    linewidth = 0.45
  ) +
  facet_grid(
    . ~ CellType,
    scales = "free_x",
    space = "free_x",
    drop = TRUE
  ) +
  scale_x_discrete(
    labels = function(x) sub("^.*___", "", x)
  ) +
  scale_fill_manual(
    values = setNames(
      c("#D55E00", "#0072B2"),
      c(
        paste0("Median < ", MEDIAN_CUTOFF),
        paste0("Median > ", MEDIAN_CUTOFF)
      )
    )
  ) +
  labs(
    x = NULL,
    y = "DepMap gene dependency score",
    fill = NULL,
    title = paste0(test, " TWAS genes: CRC cell-line dependency"),
    subtitle = paste0(
      "All genes with median score < ", MEDIAN_CUTOFF,
      " plus the ", N_EXTRA,
      " lowest-median genes with median score > ", MEDIAN_CUTOFF,
      " per cell type; genes ordered by median score"
    )
  ) +
  theme_bw(base_size = 11) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "grey95", color = "grey60"),
    strip.text = element_text(face = "bold", size = 11),
    axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1, size = 8),
    legend.position = "top",
    plot.title = element_text(face = "bold")
  )

# Dynamic width helps keep gene labels readable in the combined figure.
n_plot_genes <- nrow(selected_df)
fig_width <- max(10, min(30, 0.42 * n_plot_genes + 4))
fig_height <- 6.5

pdf_file <- file.path(
  OUTDIR,
  paste0(test, "_DepMap_CRC_boxplot_by_TWAS_CellType.pdf")
)
png_file <- file.path(
  OUTDIR,
  paste0(test, "_DepMap_CRC_boxplot_by_TWAS_CellType_300dpi.png")
)

ggsave(pdf_file, p, width = fig_width, height = fig_height, units = "in")
ggsave(
  png_file, p,
  width = fig_width,
  height = fig_height,
  units = "in",
  dpi = 300
)

cat("Saved combined boxplot PDF:\n", pdf_file, "\n")
cat("Saved combined boxplot PNG (300 dpi):\n", png_file, "\n")

# Console summary of selected genes in final plotting order.
print(
  selected_df %>%
    select(CellType, gene_name, median_score, n_cell_lines, Selection) %>%
    arrange(CellType, median_score)
)
