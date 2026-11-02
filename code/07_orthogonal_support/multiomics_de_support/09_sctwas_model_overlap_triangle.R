## 09_sctwas_model_overlap_triangle.R
## Gene-set overlap between the 16 scTWAS cell-type models (7 normal colon + 9 polyp).
## A "model" = one Dataset x CellType combination; overlap is counted on distinct
## gene symbols. Output: lower-triangle count heatmap (PNG + PDF + source data)
## and the 16x16 matrix / 120-pair tables.
##
## Inputs (first reachable one is used):
##   RAW_CSV    data/scTWAS_ADSER_NL_collected_results.csv  (546 gene x cell-type rows)
##   MATRIX_CSV table/scTWAS_model_gene_overlap_matrix.csv   (precomputed 16x16 counts)
## The RAW path also regenerates the two tables; the MATRIX path only redraws the figure.

suppressPackageStartupMessages(library(ggplot2))

root       <- "W:/VUMC_projects/scTWAS/MultiOmics_DE_support"
RAW_CSV    <- file.path(root, "data", "scTWAS_ADSER_NL_collected_results.csv")
MATRIX_CSV <- file.path(root, "table", "scTWAS_model_gene_overlap_matrix.csv")
OUT_DIR    <- root                      # figure/ and table/ are created under here
fig_stem   <- "Fig_scTWAS_model_overlap_triangle"

fig_dir <- file.path(OUT_DIR, "figure"); dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
tab_dir <- file.path(OUT_DIR, "table");  dir.create(tab_dir, showWarnings = FALSE, recursive = TRUE)

## ---- model order (as requested: matched cell types first, then polyp-only ASC/SSC) ----
nl_types  <- c("ABS", "GOB", "STM", "CT", "TAC", "TUF", "EE")
pol_types <- c("ABS", "GOB", "STM", "CT", "TAC", "TUF", "EE", "ASC", "SSC")
model_levels <- c(paste("NL", nl_types), paste("Polyps", pol_types))
n_models <- length(model_levels)
n_nl     <- length(nl_types)

## ---- 1. overlap matrix -----------------------------------------------------
if (file.exists(RAW_CSV)) {
  res <- read.csv(RAW_CSV, stringsAsFactors = FALSE)
  stopifnot(nrow(res) == 546, all(c("Dataset", "CellType", "gene_name") %in% names(res)))
  ## ADSER (adenoma / serrated) is reported as "colorectal polyps" in the figure
  res$model <- paste(ifelse(res$Dataset == "NL", "NL", "Polyps"), res$CellType)
  stopifnot(setequal(model_levels, unique(res$model)))

  gene_sets <- lapply(split(res$gene_name, factor(res$model, levels = model_levels)), unique)
  M <- outer(seq_len(n_models), seq_len(n_models),
             Vectorize(function(i, j) length(intersect(gene_sets[[i]], gene_sets[[j]]))))
  J <- outer(seq_len(n_models), seq_len(n_models),
             Vectorize(function(i, j) length(intersect(gene_sets[[i]], gene_sets[[j]])) /
                                      length(union(gene_sets[[i]], gene_sets[[j]]))))
  dimnames(M) <- dimnames(J) <- list(model_levels, model_levels)
  stopifnot(diag(M) == vapply(gene_sets, length, integer(1)))

  write.csv(data.frame(model = model_levels, M, check.names = FALSE),
            file.path(tab_dir, "scTWAS_model_gene_overlap_matrix.csv"), row.names = FALSE)

  ## 120 unordered pairs with Jaccard + shared gene lists
  idx <- which(lower.tri(M), arr.ind = TRUE)
  pair_tab <- data.frame(model_A = model_levels[idx[, "col"]],
                         model_B = model_levels[idx[, "row"]],
                         n_A = lengths(gene_sets)[idx[, "col"]],
                         n_B = lengths(gene_sets)[idx[, "row"]],
                         n_shared = M[idx], jaccard = round(J[idx], 4),
                         row.names = NULL, stringsAsFactors = FALSE)
  pair_tab$comparison <- with(pair_tab, ifelse(
    startsWith(model_A, "NL") & startsWith(model_B, "NL"), "within_NL",
    ifelse(startsWith(model_A, "Polyps") & startsWith(model_B, "Polyps"),
           "within_polyps", "NL_vs_polyps")))
  pair_tab$shared_genes <- mapply(function(a, b)
    paste(sort(intersect(gene_sets[[a]], gene_sets[[b]])), collapse = ";"),
    pair_tab$model_A, pair_tab$model_B, USE.NAMES = FALSE)
  write.csv(pair_tab[order(-pair_tab$n_shared), ],
            file.path(tab_dir, "scTWAS_model_gene_overlap_pairs.csv"), row.names = FALSE)
} else {
  ## figure-only redraw from the precomputed matrix (row names may use the old AD prefix)
  mt <- read.csv(MATRIX_CSV, stringsAsFactors = FALSE, check.names = FALSE)
  rn <- sub("^AD ", "Polyps ", mt[[1]])
  M  <- as.matrix(mt[, -1]); rownames(M) <- rn
  colnames(M) <- sub("^AD ", "Polyps ", colnames(M))
  M <- M[model_levels, model_levels]
}

## ---- 2. plotted values (source data) --------------------------------------
plot_df <- as.data.frame(as.table(M), stringsAsFactors = FALSE)
names(plot_df) <- c("row_model", "col_model", "n_shared")
plot_df$i <- match(plot_df$row_model, model_levels)
plot_df$j <- match(plot_df$col_model, model_levels)
plot_df <- plot_df[plot_df$i >= plot_df$j, ]          # lower triangle + diagonal
plot_df$cell <- ifelse(plot_df$i == plot_df$j, "diagonal (model size)", "pair overlap")
plot_df$x <- factor(plot_df$col_model, levels = model_levels)
plot_df$y <- factor(plot_df$row_model, levels = rev(model_levels))
write.csv(plot_df[order(plot_df$i, plot_df$j), c("row_model", "col_model", "n_shared", "cell")],
          file.path(fig_dir, paste0(fig_stem, "_sourcedata.csv")), row.names = FALSE)

diag_df <- plot_df[plot_df$i == plot_df$j, ]
off_df  <- plot_df[plot_df$i != plot_df$j, ]
vmax    <- max(off_df$n_shared)

## ---- 3. figure -------------------------------------------------------------
blk <- n_models - n_nl + 0.5            # normal | polyp boundary on the reversed y scale
p <- ggplot() +
  geom_tile(data = off_df, aes(x, y, fill = n_shared), colour = "white", linewidth = 0.4) +
  geom_tile(data = diag_df, aes(x, y), fill = "#d9d9d9", colour = "white", linewidth = 0.4) +
  geom_text(data = off_df, aes(x, y, label = n_shared,
                               colour = n_shared > 0.6 * vmax), size = 2.5) +
  geom_text(data = diag_df, aes(x, y, label = n_shared),
            size = 2.5, fontface = "bold", colour = "#2b2b2b") +
  scale_colour_manual(values = c("FALSE" = "#2b2b2b", "TRUE" = "white"), guide = "none") +
  scale_fill_gradient(low = "#f7fbff", high = "#08519c", limits = c(0, vmax),
                      name = "shared\ngenes") +
  geom_vline(xintercept = n_nl + 0.5, colour = "#444444", linewidth = 0.4) +
  geom_hline(yintercept = blk, colour = "#444444", linewidth = 0.4) +
  annotate("text", x = 3.5, y = n_models + 0.9, label = "normal colon (NL)",
           size = 2.8, colour = "#444444") +
  annotate("text", x = 12.5, y = blk - 0.6, label = "colorectal polyps",
           size = 2.8, colour = "#444444") +
  scale_x_discrete(expand = expansion(add = 0.5)) +
  scale_y_discrete(expand = expansion(add = 0.5)) +
  coord_cartesian(clip = "off") +
  labs(x = NULL, y = NULL,
       title = "scTWAS gene sets overlap within normal colon, but barely within polyps") +
  theme_minimal(base_size = 9) +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
        axis.text.y = element_text(size = 7),
        axis.ticks = element_blank(),
        plot.title = element_text(hjust = 0, size = 9, margin = margin(b = 12)),
        plot.margin = margin(10, 10, 6, 6),
        legend.title = element_text(size = 7), legend.text = element_text(size = 6.5),
        legend.key.width = unit(0.28, "cm"), legend.key.height = unit(0.5, "cm"))

ggsave(file.path(fig_dir, paste0(fig_stem, ".png")), plot = p,
       width = 8.4, height = 7.6, dpi = 300)
ggsave(file.path(fig_dir, paste0(fig_stem, ".pdf")), plot = p,
       width = 8.4, height = 7.6, device = cairo_pdf)

cat("diagonal:", paste(model_levels, diag(M), sep = "="), "\n")
cat("max off-diagonal:", vmax, "\n")
