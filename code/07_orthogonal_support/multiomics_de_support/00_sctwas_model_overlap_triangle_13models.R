## 09_sctwas_model_overlap_triangle_13models.R
## Gene-set overlap between 13 scTWAS cell-type models:
## 7 normal-colon (NL) models + 6 Polyps models.
## Polyps CT, EE, and TUF are excluded.
## A "model" = one Dataset x CellType combination; overlap is counted on distinct
## gene symbols. Output: lower-triangle count heatmap (PNG + PDF + source data)
## and the 13x13 matrix / 78-pair tables.

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
})

root    <- "W:/VUMC_projects/scTWAS/VUMC_HTAN_scTWAS_genes_DE_support"
in_csv  <- file.path(root, "data", "scTWAS_ADSER_NL_collected_results.csv")
fig_dir <- file.path(root, "figure")
tab_dir <- file.path(root, "table")
dir.create(fig_dir, showWarnings = FALSE)
dir.create(tab_dir, showWarnings = FALSE)

fig_stem <- "Fig_scTWAS_model_overlap_triangle"

## ---- 1. models -------------------------------------------------------------
res <- read.csv(in_csv, stringsAsFactors = FALSE)
stopifnot(
  nrow(res) == 546,
  all(c("Dataset", "CellType", "gene_name") %in% names(res))
)

## Keep all 7 NL models, but remove CT, EE, and TUF from the Polyps dataset.
res <- res %>%
  filter(!(Dataset != "NL" & CellType %in% c("CT", "EE", "TUF")))

## Use "Polyps" consistently instead of AD / ADSER in model labels and plot.
res$model <- paste(ifelse(res$Dataset == "NL", "NL", "Polyps"), res$CellType)

model_levels <- c(
  paste("NL", c("ABS", "GOB", "STM", "CT", "TAC", "TUF", "EE")),
  paste("Polyps", c("ABS", "ASC", "GOB", "TAC", "SSC", "STM"))
)

## Diagnostic information if model names do not match the expected 13 models.
missing_models <- setdiff(model_levels, unique(res$model))
extra_models   <- setdiff(unique(res$model), model_levels)

if (length(missing_models) > 0) {
  stop("Expected models missing from data: ", paste(missing_models, collapse = ", "))
}
if (length(extra_models) > 0) {
  stop("Unexpected models found in data: ", paste(extra_models, collapse = ", "))
}

stopifnot(setequal(model_levels, unique(res$model)))

gene_sets <- split(res$gene_name, factor(res$model, levels = model_levels))
gene_sets <- lapply(gene_sets, unique)
n_models  <- length(model_levels)
stopifnot(n_models == 13)

## ---- 2. 13 x 13 overlap matrix --------------------------------------------
M <- outer(
  seq_len(n_models), seq_len(n_models),
  Vectorize(function(i, j) length(intersect(gene_sets[[i]], gene_sets[[j]])))
)

J <- outer(
  seq_len(n_models), seq_len(n_models),
  Vectorize(function(i, j) {
    length(intersect(gene_sets[[i]], gene_sets[[j]])) /
      length(union(gene_sets[[i]], gene_sets[[j]]))
  })
)

dimnames(M) <- dimnames(J) <- list(model_levels, model_levels)
stopifnot(diag(M) == vapply(gene_sets, length, integer(1)))

write.csv(
  data.frame(model = model_levels, M, check.names = FALSE),
  file.path(tab_dir, "scTWAS_model_gene_overlap_matrix.csv"),
  row.names = FALSE
)

## Long-form pair table: choose(13, 2) = 78 unordered pairs.
pairs <- which(lower.tri(M), arr.ind = TRUE)

pair_tab <- data.frame(
  model_A   = model_levels[pairs[, "col"]],
  model_B   = model_levels[pairs[, "row"]],
  n_A       = lengths(gene_sets)[pairs[, "col"]],
  n_B       = lengths(gene_sets)[pairs[, "row"]],
  n_shared  = M[pairs],
  jaccard   = round(J[pairs], 4),
  row.names = NULL,
  stringsAsFactors = FALSE
)

pair_tab$comparison <- with(
  pair_tab,
  ifelse(
    startsWith(model_A, "NL") & startsWith(model_B, "NL"),
    "within_NL",
    ifelse(
      startsWith(model_A, "Polyps") & startsWith(model_B, "Polyps"),
      "within_Polyps",
      "NL_vs_Polyps"
    )
  )
)

pair_tab$shared_genes <- mapply(
  function(a, b) {
    paste(sort(intersect(gene_sets[[a]], gene_sets[[b]])), collapse = ";")
  },
  pair_tab$model_A,
  pair_tab$model_B,
  USE.NAMES = FALSE
)

pair_tab <- pair_tab[order(-pair_tab$n_shared), ]

write.csv(
  pair_tab,
  file.path(tab_dir, "scTWAS_model_gene_overlap_pairs.csv"),
  row.names = FALSE
)

## ---- 3. plotted values (source data) --------------------------------------
plot_df <- as.data.frame(as.table(M), stringsAsFactors = FALSE)
names(plot_df) <- c("row_model", "col_model", "n_shared")
plot_df$i <- match(plot_df$row_model, model_levels)
plot_df$j <- match(plot_df$col_model, model_levels)
plot_df <- plot_df[plot_df$i >= plot_df$j, ]  # lower triangle + diagonal
plot_df$cell <- ifelse(
  plot_df$i == plot_df$j,
  "diagonal (model size)",
  "pair overlap"
)
plot_df$x <- factor(plot_df$col_model, levels = model_levels)
plot_df$y <- factor(plot_df$row_model, levels = rev(model_levels))

write.csv(
  plot_df[order(plot_df$i, plot_df$j),
          c("row_model", "col_model", "n_shared", "cell")],
  file.path(fig_dir, paste0(fig_stem, "_sourcedata.csv")),
  row.names = FALSE
)

diag_df <- plot_df[plot_df$i == plot_df$j, ]
off_df  <- plot_df[plot_df$i != plot_df$j, ]
vmax    <- max(off_df$n_shared)

## ---- 4. figure -------------------------------------------------------------
## x-axis: NL models occupy positions 1-7, Polyps positions 8-13.
## y-axis is reversed, so the NL | Polyps boundary is at n_models - 7 + 0.5.
blk <- n_models - 7 + 0.5

p <- ggplot() +
  geom_tile(
    data = off_df,
    aes(x, y, fill = n_shared),
    colour = "white",
    linewidth = 0.4
  ) +
  geom_tile(
    data = diag_df,
    aes(x, y),
    fill = "#d9d9d9",
    colour = "white",
    linewidth = 0.4
  ) +
  geom_text(
    data = off_df,
    aes(x, y, label = n_shared, colour = n_shared > 0.6 * vmax),
    size = 2.5
  ) +
  geom_text(
    data = diag_df,
    aes(x, y, label = n_shared),
    size = 2.5,
    fontface = "bold",
    colour = "#2b2b2b"
  ) +
  scale_colour_manual(
    values = c("FALSE" = "#2b2b2b", "TRUE" = "white"),
    guide = "none"
  ) +
  scale_fill_gradient(
    low = "#f7fbff",
    high = "#08519c",
    limits = c(0, vmax),
    name = "shared\ngenes"
  ) +
  geom_vline(xintercept = 7.5, colour = "#444444", linewidth = 0.4) +
  geom_hline(yintercept = blk, colour = "#444444", linewidth = 0.4) +
  annotate(
    "text",
    x = 4,
    y = n_models + 0.9,
    label = "normal colon (NL)",
    size = 2.8,
    colour = "#444444"
  ) +
  annotate(
    "text",
    x = 10.5,
    y = blk - 0.6,
    label = "Polyps",
    size = 2.8,
    colour = "#444444"
  ) +
  scale_x_discrete(position = "bottom", expand = expansion(add = 0.5)) +
  scale_y_discrete(expand = expansion(add = 0.5)) +
  coord_cartesian(clip = "off") +
  theme_minimal(base_size = 9) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
    axis.text.y = element_text(size = 7),
    axis.ticks = element_blank(),
    plot.title = element_text(hjust = 0, size = 9, margin = margin(b = 12)),
    plot.caption = element_text(hjust = 0, size = 6.5, colour = "grey40"),
    plot.margin = margin(10, 10, 6, 6),
    legend.title = element_text(size = 7),
    legend.text = element_text(size = 6.5),
    legend.key.width = grid::unit(0.28, "cm"),
    legend.key.height = grid::unit(0.5, "cm")
  )

ggsave(
  file.path(fig_dir, paste0(fig_stem, ".png")),
  plot = p,
  width = 8.4,
  height = 7.6,
  dpi = 300
)

ggsave(
  file.path(fig_dir, paste0(fig_stem, ".pdf")),
  plot = p,
  width = 8.4,
  height = 7.6,
  device = cairo_pdf
)

## ---- 5. summary ------------------------------------------------------------
cat("Models retained:", n_models, "(7 NL + 6 Polyps)\n")
cat("Pairwise comparisons:", nrow(pair_tab), "\n")

cat(
  "within_NL mean/max:",
  mean(pair_tab$n_shared[pair_tab$comparison == "within_NL"]),
  max(pair_tab$n_shared[pair_tab$comparison == "within_NL"]),
  "\n"
)

cat(
  "within_Polyps mean/max:",
  mean(pair_tab$n_shared[pair_tab$comparison == "within_Polyps"]),
  max(pair_tab$n_shared[pair_tab$comparison == "within_Polyps"]),
  "\n"
)

cat(
  "NL_vs_Polyps mean/max:",
  mean(pair_tab$n_shared[pair_tab$comparison == "NL_vs_Polyps"]),
  max(pair_tab$n_shared[pair_tab$comparison == "NL_vs_Polyps"]),
  "\n"
)
