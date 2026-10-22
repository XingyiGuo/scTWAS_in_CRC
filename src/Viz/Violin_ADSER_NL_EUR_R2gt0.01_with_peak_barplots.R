#!/usr/bin/env Rscript
set.seed(123)

# ==============================================================================
# Violin plots of predictable genes in ADSER and Normal-colon (NL)
# cell-type-specific EUR TWAS
#
# Predictable/tested gene definition (identical to Bonferroni denominator):
#   1) valid EUR TWAS pvalue (0 <= pvalue <= 1), AND
#   2) R2_global > 0.01
#
# This single script processes BOTH:
#   - ADSER: 9 cell types
#   - NL:    7 cell types
#
# For each dataset it writes:
#   - one PNG with all cell types
#   - one PDF with all cell types
#   - one summary TSV
#   - one gene-level TSV
#   - one peak-region count bar plot (PNG + PDF)
#   - one peak-region count TSV
# ==============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(ggplot2)
  library(scales)
})

# ------------------------------------------------------------------------------
# User settings
# ------------------------------------------------------------------------------

ASSOC_ROOT <- "/v5000/data/l2_bioinfo1/liq17/scTWAS/TWAS_HTAN_GSE266546_WellcomeSanger_res/SPrediXcan_peakSNP_association"

R2_THRESHOLD <- 0.01
ALPHA <- 0.05

# Same cell-type palette used in the Manhattan plots.
CT_COLORS <- c(
  "ASC" = "#A6611A",  # ochre-brown
  "SSC" = "#882255",  # wine/magenta
  "ABS" = "#CC79A7",  # pink-purple
  "GOB" = "#0072B2",  # blue
  "CT"  = "#E69F00",  # orange
  "EE"  = "#8C564B",  # brown
  "STM" = "#332288",  # deep indigo
  "TUF" = "#999933",  # olive
  "TAC" = "#56B4E9"   # sky blue
)

# Dataset-specific configuration.
DATASETS <- list(

  ADSER = list(
    input_pattern = file.path(
      ASSOC_ROOT,
      "VUMC_HTAN_ADSER_%s_peakSNP.EUR.TWAS.csv"
    ),
    cell_types = c("ABS", "GOB", "CT", "EE", "STM", "TUF", "TAC", "ASC", "SSC"),
    outdir = file.path(
      ASSOC_ROOT,
      "VUMC_HTAN_ADSER_EUR_Multi_sig/Violin_ADSER_EUR_R2gt0.01"
    ),
    prefix = "ADSER_9CT_EUR",
    title = "Prediction performance of ADSER cell-type-specific EUR TWAS models"
  ),

  NL = list(
    input_pattern = file.path(
      ASSOC_ROOT,
      "VUMC_HTAN_GSE266546_NL_%s_peakSNP.EUR.TWAS.csv"
    ),
    cell_types = c("ABS", "GOB", "CT", "STM", "EE", "TUF", "TAC"),
    outdir = file.path(
      ASSOC_ROOT,
      "Cauchy_Meta_NL_WellcomeSanger_NL_EUR/Violin_VUMC_HTAN_GSE266546_NL_EUR_R2gt0.01"
    ),
    prefix = "NL_7CT_EUR",
    title = "Prediction performance of normal-colon cell-type-specific EUR TWAS models"
  )
)


# ------------------------------------------------------------------------------
# Peak-region counts for bar plots
# ------------------------------------------------------------------------------

PEAK_COUNTS <- tibble::tribble(
  ~Group,  ~CellType, ~Num_regions_after_chain_merge, ~Note,
  "NL",    "ABS",     150086,                         "VUMC_GSE201342",
  "NL",    "CT",      152781,                         "VUMC_GSE201343",
  "NL",    "EE",       68521,                         "VUMC_GSE201344",
  "NL",    "GOB",     152548,                         "VUMC_GSE201345",
  "NL",    "STM",     145960,                         "VUMC_GSE201346",
  "NL",    "TAC",     149672,                         "VUMC_GSE201347",
  "NL",    "TUF",      55454,                         "VUMC",
  "ADSER", "ABS",     188833,                         "VUMC_GSE201336",
  "ADSER", "CT",      135528,                         "VUMC_GSE201337",
  "ADSER", "EE",      109919,                         "VUMC_GSE201338",
  "ADSER", "GOB",     154237,                         "VUMC_GSE201339",
  "ADSER", "STM",     147288,                         "VUMC_GSE201340",
  "ADSER", "TAC",     165047,                         "VUMC_GSE201341",
  "ADSER", "TUF",      16018,                         "VUMC",
  "ADSER", "ASC",      28652,                         "VUMC",
  "ADSER", "SSC",      13061,                         "VUMC"
)

# ------------------------------------------------------------------------------
# Peak-region count bar plot
# ------------------------------------------------------------------------------

make_peak_barplot <- function(dataset_name, cfg) {

  CTS <- cfg$cell_types

  peak_df <- PEAK_COUNTS %>%
    dplyr::filter(Group == dataset_name) %>%
    dplyr::mutate(
      CellType = factor(CellType, levels = CTS)
    ) %>%
    dplyr::arrange(CellType)

  missing_ct <- setdiff(CTS, as.character(peak_df$CellType))
  if (length(missing_ct) > 0) {
    stop(
      dataset_name,
      ": missing peak-region count(s) for cell type(s): ",
      paste(missing_ct, collapse = ", ")
    )
  }

  # Save the hard-coded peak-count table in the same violin output folder.
  PEAK_TSV <- file.path(
    cfg$outdir,
    paste0(cfg$prefix, "_peak_region_counts.tsv")
  )

  data.table::fwrite(
    as.data.frame(peak_df),
    PEAK_TSV,
    sep = "\t"
  )

  y_max <- max(peak_df$Num_regions_after_chain_merge, na.rm = TRUE)

  p_bar <- ggplot2::ggplot(
    peak_df,
    ggplot2::aes(
      x = CellType,
      y = Num_regions_after_chain_merge,
      fill = CellType
    )
  ) +

    ggplot2::geom_col(
      width = 0.72,
      color = "black",
      linewidth = 0.25
    ) +

    ggplot2::geom_text(
      ggplot2::aes(
        label = scales::comma(Num_regions_after_chain_merge)
      ),
      vjust = -0.35,
      size = 3.1
    ) +

    ggplot2::scale_fill_manual(
      values = CT_COLORS,
      guide = "none"
    ) +

    ggplot2::scale_y_continuous(
      labels = scales::comma,
      limits = c(0, y_max * 1.14),
      expand = ggplot2::expansion(mult = c(0, 0))
    ) +

    ggplot2::labs(
      title = ifelse(
        dataset_name == "NL",
        "Normal colon cell-type-specific scATAC peak regions",
        "Adenoma/serrated cell-type-specific scATAC peak regions"
      ),
      x = "Cell type",
      y = "Number of merged peak regions"
    ) +

    ggplot2::theme_classic(base_size = 11) +

    ggplot2::theme(
      legend.position = "none",
      plot.title = ggplot2::element_text(
        size = 13,
        face = "bold",
        hjust = 0.5
      ),
      axis.title = ggplot2::element_text(size = 10),
      axis.text.x = ggplot2::element_text(
        size = 9,
        face = "bold",
        color = "black"
      ),
      axis.text.y = ggplot2::element_text(
        size = 9,
        color = "black"
      ),
      axis.ticks = ggplot2::element_line(linewidth = 0.3),
      plot.margin = ggplot2::margin(8, 12, 8, 8)
    )

  PNG_OUT <- file.path(
    cfg$outdir,
    paste0(cfg$prefix, "_peak_region_counts_barplot.png")
  )

  PDF_OUT <- file.path(
    cfg$outdir,
    paste0(cfg$prefix, "_peak_region_counts_barplot.pdf")
  )

  ggplot2::ggsave(
    PNG_OUT,
    plot = p_bar,
    width = 8.5,
    height = 5.2,
    dpi = 300,
    bg = "white"
  )

  ggplot2::ggsave(
    PDF_OUT,
    plot = p_bar,
    width = 8.5,
    height = 5.2,
    device = grDevices::cairo_pdf,
    bg = "white"
  )

  message("  Peak bar PNG: ", PNG_OUT)
  message("  Peak bar PDF: ", PDF_OUT)
  message("  Peak counts  : ", PEAK_TSV)

  invisible(
    list(
      data = peak_df,
      plot = p_bar,
      png = PNG_OUT,
      pdf = PDF_OUT,
      tsv = PEAK_TSV
    )
  )
}

# ------------------------------------------------------------------------------
# Read/filter one cell type
# ------------------------------------------------------------------------------

read_one_ct <- function(ct, input_pattern, dataset_name) {

  infile <- sprintf(input_pattern, ct)

  if (!file.exists(infile)) {
    stop(
      dataset_name, " / ", ct,
      ": input file not found:\n  ", infile
    )
  }

  message("Reading ", dataset_name, " / ", ct, ": ", infile)

  df <- data.table::fread(infile, data.table = FALSE)

  required <- c("gene", "gene_name", "pvalue", "R2_global")
  missing_cols <- setdiff(required, names(df))

  if (length(missing_cols) > 0) {
    stop(
      dataset_name, " / ", ct,
      ": missing required column(s): ",
      paste(missing_cols, collapse = ", "),
      "\nFile: ", infile
    )
  }

  p_all  <- suppressWarnings(as.numeric(df$pvalue))
  r2_all <- suppressWarnings(as.numeric(df$R2_global))

  # Exact EUR TWAS testing universe used for Bonferroni correction.
  valid_p <-
    !is.na(p_all) &
    is.finite(p_all) &
    p_all >= 0 &
    p_all <= 1

  eligible_test <-
    valid_p &
    !is.na(r2_all) &
    is.finite(r2_all) &
    r2_all > R2_THRESHOLD

  n_total_rows <- nrow(df)
  n_valid_p    <- sum(valid_p)
  n_tests      <- sum(eligible_test)

  if (n_tests == 0) {
    stop(
      dataset_name, " / ", ct,
      ": no eligible EUR TWAS tests with R2_global > ",
      R2_THRESHOLD
    )
  }

  bonf_cutoff <- ALPHA / n_tests

  # Diagnostic: denominator is eligible TWAS rows/tests.
  tested_gene_ids <- sub("\\..*$", "", as.character(df$gene[eligible_test]))
  n_unique_genes  <- dplyr::n_distinct(tested_gene_ids)

  if (n_unique_genes != n_tests) {
    warning(
      dataset_name, " / ", ct, ": ",
      n_tests, " eligible EUR TWAS rows but ",
      n_unique_genes, " unique gene IDs. ",
      "Bonferroni denominator remains eligible test rows."
    )
  }

  out <- df %>%
    dplyr::transmute(
      Dataset = dataset_name,
      CellType = ct,
      gene = sub("\\..*$", "", as.character(gene)),
      gene_name = as.character(gene_name),
      pvalue = suppressWarnings(as.numeric(pvalue)),
      R2_global = suppressWarnings(as.numeric(R2_global))
    ) %>%
    dplyr::filter(
      !is.na(pvalue),
      is.finite(pvalue),
      pvalue >= 0,
      pvalue <= 1,
      !is.na(R2_global),
      is.finite(R2_global),
      R2_global > R2_THRESHOLD
    ) %>%
    dplyr::mutate(
      N_total_rows = n_total_rows,
      N_valid_EUR_pvalues = n_valid_p,
      N_TWAS_tests = n_tests,
      Bonferroni_cutoff = bonf_cutoff
    )

  message(
    "  Total rows = ", n_total_rows,
    " | valid EUR P values = ", n_valid_p,
    " | predictable/tested genes = ", n_tests,
    " | Bonferroni cutoff = ",
    format(bonf_cutoff, scientific = TRUE, digits = 4)
  )

  out
}

# ------------------------------------------------------------------------------
# Process one dataset and make violin plot
# ------------------------------------------------------------------------------

run_dataset <- function(dataset_name, cfg) {

  CTS <- cfg$cell_types

  stopifnot(all(CTS %in% names(CT_COLORS)))

  dir.create(
    cfg$outdir,
    recursive = TRUE,
    showWarnings = FALSE
  )

  message("")
  message("==============================================================")
  message("PROCESSING DATASET: ", dataset_name)
  message("==============================================================")

  ALL <- lapply(
    CTS,
    function(ct) {
      read_one_ct(
        ct = ct,
        input_pattern = cfg$input_pattern,
        dataset_name = dataset_name
      )
    }
  ) %>%
    dplyr::bind_rows() %>%
    dplyr::mutate(
      CellType = factor(CellType, levels = CTS)
    )

  # --------------------------------------------------------------------------
  # Summary table
  # --------------------------------------------------------------------------

  SUMMARY <- ALL %>%
    dplyr::group_by(Dataset, CellType) %>%
    dplyr::summarise(
      N_total_rows = dplyr::first(N_total_rows),
      N_valid_EUR_pvalues = dplyr::first(N_valid_EUR_pvalues),
      N_predictable_genes_R2gt0.01 = dplyr::n(),
      Bonferroni_cutoff = dplyr::first(Bonferroni_cutoff),
      R2_min = min(R2_global, na.rm = TRUE),
      R2_Q1 = as.numeric(stats::quantile(R2_global, 0.25, na.rm = TRUE)),
      R2_median = stats::median(R2_global, na.rm = TRUE),
      R2_mean = mean(R2_global, na.rm = TRUE),
      R2_Q3 = as.numeric(stats::quantile(R2_global, 0.75, na.rm = TRUE)),
      R2_max = max(R2_global, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      CellType = factor(CellType, levels = CTS)
    ) %>%
    dplyr::arrange(CellType)

  SUMMARY_OUT <- file.path(
    cfg$outdir,
    paste0(cfg$prefix, "_predictable_genes_R2gt0.01_summary.tsv")
  )

  GENE_OUT <- file.path(
    cfg$outdir,
    paste0(cfg$prefix, "_predictable_genes_R2gt0.01_gene_level.tsv")
  )

  data.table::fwrite(
    as.data.frame(SUMMARY),
    SUMMARY_OUT,
    sep = "\t"
  )

  data.table::fwrite(
    as.data.frame(ALL),
    GENE_OUT,
    sep = "\t"
  )

  print(SUMMARY)

  # --------------------------------------------------------------------------
  # Violin plot
  # --------------------------------------------------------------------------

  y_data_max <- max(ALL$R2_global, na.rm = TRUE)
  y_label    <- y_data_max * 1.075
  y_top      <- y_data_max * 1.16

  LABEL_DF <- SUMMARY %>%
    dplyr::mutate(
      CellType = factor(CellType, levels = CTS),
      y = y_label,
      label = paste0(
        "N = ",
        scales::comma(N_predictable_genes_R2gt0.01)
      )
    )

  p <- ggplot2::ggplot(
    ALL,
    ggplot2::aes(
      x = CellType,
      y = R2_global,
      fill = CellType,
      color = CellType
    )
  ) +

    ggplot2::geom_violin(
      scale = "width",
      trim = TRUE,
      linewidth = 0.45,
      alpha = 0.72
    ) +

    ggplot2::geom_boxplot(
      width = 0.13,
      outlier.shape = NA,
      fill = "white",
      color = "black",
      linewidth = 0.35
    ) +

    ggplot2::stat_summary(
      fun = median,
      geom = "point",
      shape = 21,
      size = 1.8,
      stroke = 0.35,
      fill = "white",
      color = "black"
    ) +

    ggplot2::geom_hline(
      yintercept = R2_THRESHOLD,
      linetype = "dashed",
      linewidth = 0.45,
      color = "grey35"
    ) +

    ggplot2::geom_text(
      data = LABEL_DF,
      ggplot2::aes(
        x = CellType,
        y = y,
        label = label
      ),
      inherit.aes = FALSE,
      size = 3.1,
      fontface = "bold"
    ) +

    ggplot2::scale_fill_manual(
      values = CT_COLORS,
      guide = "none"
    ) +

    ggplot2::scale_color_manual(
      values = CT_COLORS,
      guide = "none"
    ) +

    ggplot2::scale_y_continuous(
      limits = c(0, y_top),
      breaks = scales::pretty_breaks(n = 6),
      expand = ggplot2::expansion(mult = c(0, 0))
    ) +

    ggplot2::labs(
      title = cfg$title,
      subtitle = paste0(
        "Predictable/tested genes: valid EUR TWAS P value and ",
        "R2_global > ", R2_THRESHOLD,
        "; N = Bonferroni denominator"
      ),
      x = "Cell type",
      y = expression(R^2[global])
    ) +

    ggplot2::theme_classic(base_size = 11) +

    ggplot2::theme(
      legend.position = "none",
      plot.title = ggplot2::element_text(
        size = 13,
        face = "bold",
        hjust = 0.5
      ),
      plot.subtitle = ggplot2::element_text(
        size = 9,
        hjust = 0.5
      ),
      axis.title = ggplot2::element_text(size = 10),
      axis.text.x = ggplot2::element_text(
        size = 9,
        face = "bold",
        color = "black"
      ),
      axis.text.y = ggplot2::element_text(
        size = 9,
        color = "black"
      ),
      axis.ticks = ggplot2::element_line(linewidth = 0.3),
      plot.margin = ggplot2::margin(8, 12, 8, 8)
    )

  print(p)

  # --------------------------------------------------------------------------
  # Save
  # --------------------------------------------------------------------------

  PNG_OUT <- file.path(
    cfg$outdir,
    paste0(cfg$prefix, "_predictable_genes_R2global_violin.png")
  )

  PDF_OUT <- file.path(
    cfg$outdir,
    paste0(cfg$prefix, "_predictable_genes_R2global_violin.pdf")
  )

  ggplot2::ggsave(
    PNG_OUT,
    plot = p,
    width = 9.2,
    height = 5.8,
    dpi = 300,
    bg = "white"
  )

  ggplot2::ggsave(
    PDF_OUT,
    plot = p,
    width = 9.2,
    height = 5.8,
    device = grDevices::cairo_pdf,
    bg = "white"
  )

  message("")
  message("Done: ", dataset_name)
  message("  PNG     : ", PNG_OUT)
  message("  PDF     : ", PDF_OUT)
  message("  Summary : ", SUMMARY_OUT)
  message("  Gene TSV: ", GENE_OUT)

  # --------------------------------------------------------------------------
  # Peak-region count bar plot
  # --------------------------------------------------------------------------

  PEAK_BAR <- make_peak_barplot(
    dataset_name = dataset_name,
    cfg = cfg
  )

  invisible(
    list(
      data = ALL,
      summary = SUMMARY,
      plot = p,
      png = PNG_OUT,
      pdf = PDF_OUT,
      peak_bar = PEAK_BAR
    )
  )
}

# ------------------------------------------------------------------------------
# Run ADSER and NL from the same script
# ------------------------------------------------------------------------------

RESULTS <- lapply(
  names(DATASETS),
  function(dataset_name) {
    run_dataset(
      dataset_name = dataset_name,
      cfg = DATASETS[[dataset_name]]
    )
  }
)

names(RESULTS) <- names(DATASETS)

message("")
message("==============================================================")
message("ALL DATASETS COMPLETE")
message("==============================================================")
message("Processed: ", paste(names(DATASETS), collapse = ", "))
