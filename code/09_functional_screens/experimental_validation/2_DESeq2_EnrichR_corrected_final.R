#!/usr/bin/env Rscript

# ============================================================
# DESeq2 significant genes -> EnrichR pathway enrichment
#
# Input:
#   HCT116_SW480_7KD_DESeq2_significant_results.csv
#
# Enrichment is performed separately for:
#   CellType x KD_Gene x Direction
#
# Databases:
#   1. Reactome Pathways 2024
#   2. KEGG 2021 Human
#   3. MSigDB Hallmark 2020
#
# Significant pathway:
#   EnrichR Adjusted.P.value < 0.10
#
# Outputs:
#   1. Individual enrichment results for each
#      CellType x KD_Gene x Direction
#   2. One combined enrichment TSV
# ============================================================


# ============================================================
# 0. Packages
# ============================================================

suppressPackageStartupMessages({
  library(data.table)
  library(enrichR)
})


# ============================================================
# 1. Paths and parameters
# ============================================================

INPUT_FILE <- paste0(
  "/data/l2_bioinfo1/liq17/scTWAS/Experiments_res/RNA_seq_DE/",
  "HCT116_SW480_7KD_DESeq2_significant_results.csv"
)

ENRICHR_OUTDIR <- "/data/l2_bioinfo1/liq17/scTWAS/Experiments_res/RNA_seq_DE/EnrichR/"

COMBINED_OUTDIR <- "/data/l2_bioinfo1/liq17/scTWAS/Experiments_res/RNA_seq_DE/"

COMBINED_OUTFILE <- file.path(
  COMBINED_OUTDIR,
  "DESeq2_sig_EnrichR_pathways.tsv"
)

# EnrichR FDR threshold
ENRICHR_FDR <- 0.10

# EnrichR databases
dbs_lib <- c(
  "Reactome_Pathways_2024",
  "KEGG_2021_Human",
  "MSigDB_Hallmark_2020"
)


# ============================================================
# 2. Create output directories
# ============================================================

dir.create(
  ENRICHR_OUTDIR,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  COMBINED_OUTDIR,
  recursive = TRUE,
  showWarnings = FALSE
)


# ============================================================
# 3. Configure EnrichR
# ============================================================

message("Configuring EnrichR...")

setEnrichrSite("Enrichr")


# ============================================================
# 4. Function: format and filter EnrichR output
# ============================================================

format_output <- function(
    enriched,
    output_prefix,
    write_file = TRUE,
    fdr_cutoff = 0.10
) {

  # Standard EnrichR columns that we want to retain
  desired_columns <- c(
    "Term",
    "Overlap",
    "P.value",
    "Adjusted.P.value",
    "Old.P.value",
    "Old.Adjusted.P.value",
    "Odds.Ratio",
    "Combined.Score",
    "Genes"
  )

  filtered_list <- list()


  # ----------------------------------------------------------
  # Process each enrichment database
  # ----------------------------------------------------------

  for (source_name in names(enriched)) {

    df <- as.data.frame(enriched[[source_name]])

    if (nrow(df) == 0) {
      next
    }


    # Make sure adjusted P-value exists
    if (!"Adjusted.P.value" %in% colnames(df)) {

      warning(
        "Adjusted.P.value not found for enrichment database: ",
        source_name
      )

      next
    }


    # --------------------------------------------------------
    # Filter significant pathways
    # --------------------------------------------------------

    df_sig <- df[
      !is.na(df$Adjusted.P.value) &
        df$Adjusted.P.value < fdr_cutoff,
      ,
      drop = FALSE
    ]

    if (nrow(df_sig) == 0) {
      next
    }


    # Add enrichment source/database
    df_sig$Source <- source_name


    # --------------------------------------------------------
    # Make output robust to differences between EnrichR
    # versions by creating missing columns as NA
    # --------------------------------------------------------

    for (col in desired_columns) {

      if (!col %in% colnames(df_sig)) {
        df_sig[[col]] <- NA
      }
    }


    # Keep desired columns
    df_sig <- df_sig[
      ,
      c(desired_columns, "Source"),
      drop = FALSE
    ]

    filtered_list[[source_name]] <- df_sig
  }


  # ----------------------------------------------------------
  # No significant pathway
  # ----------------------------------------------------------

  if (length(filtered_list) == 0) {

    message(
      "No significant pathways for ",
      output_prefix,
      " (Adjusted.P.value < ",
      fdr_cutoff,
      ")."
    )

    return(NULL)
  }


  # ----------------------------------------------------------
  # Combine databases
  # ----------------------------------------------------------

  combined_df <- rbindlist(
    filtered_list,
    use.names = TRUE,
    fill = TRUE
  )


  # ----------------------------------------------------------
  # Write individual result
  # ----------------------------------------------------------

  if (write_file) {

    outfile <- file.path(
      ENRICHR_OUTDIR,
      paste0(
        output_prefix,
        "_Pathway_enriched_genes.tsv"
      )
    )

    fwrite(
      combined_df,
      file = outfile,
      sep = "\t",
      quote = FALSE,
      na = ""
    )

    message(
      "Written: ",
      outfile
    )
  }


  return(combined_df)
}


# ============================================================
# 5. Read DESeq2 significant results
# ============================================================

message("Reading DESeq2 significant results...")

genes_df <- fread(
  INPUT_FILE,
  data.table = TRUE
)

message(
  "Input rows: ",
  nrow(genes_df)
)


# ============================================================
# 6. Check required columns
# ============================================================

required_columns <- c(
  "Cell",
  "KD_Gene",
  "Direction",
  "gene_name"
)

missing_columns <- setdiff(
  required_columns,
  colnames(genes_df)
)

if (length(missing_columns) > 0) {

  stop(
    "The following required columns are missing from the input file: ",
    paste(missing_columns, collapse = ", ")
  )
}


# ============================================================
# 7. Clean identifiers
# ============================================================

genes_df[
  ,
  Cell := trimws(as.character(Cell))
]

genes_df[
  ,
  KD_Gene := trimws(as.character(KD_Gene))
]

genes_df[
  ,
  Direction := trimws(as.character(Direction))
]

genes_df[
  ,
  gene_name := trimws(as.character(gene_name))
]


# Convert blank strings to NA
genes_df[
  Cell == "",
  Cell := NA_character_
]

genes_df[
  KD_Gene == "",
  KD_Gene := NA_character_
]

genes_df[
  Direction == "",
  Direction := NA_character_
]

genes_df[
  gene_name == "",
  gene_name := NA_character_
]


# ============================================================
# 8. Check input structure
# ============================================================

message("\nCell types:")
print(
  sort(unique(na.omit(genes_df$Cell)))
)

message("\nKD genes:")
print(
  sort(unique(na.omit(genes_df$KD_Gene)))
)

message("\nDirections:")
print(
  sort(unique(na.omit(genes_df$Direction)))
)


# Check expected directions
unexpected_directions <- setdiff(
  unique(na.omit(genes_df$Direction)),
  c("Up", "Down")
)

if (length(unexpected_directions) > 0) {

  warning(
    "Unexpected values found in Direction: ",
    paste(unexpected_directions, collapse = ", "),
    ". Only 'Up' and 'Down' will be analyzed."
  )
}


# ============================================================
# 9. Summarize number of significant DE genes
# ============================================================

gene_count_summary <- genes_df[
  !is.na(Cell) &
    !is.na(KD_Gene) &
    !is.na(Direction),
  .(
    N_rows = .N,
    N_unique_genes = uniqueN(
      gene_name[!is.na(gene_name)]
    )
  ),
  by = .(
    Cell,
    KD_Gene,
    Direction
  )
]

setorder(
  gene_count_summary,
  Cell,
  KD_Gene,
  Direction
)

message("\nDE gene counts:")
print(gene_count_summary)


# Write gene-count summary
GENE_COUNT_OUTFILE <- file.path(
  COMBINED_OUTDIR,
  "DESeq2_sig_gene_counts_by_CellType_KD_Direction.tsv"
)

fwrite(
  gene_count_summary,
  file = GENE_COUNT_OUTFILE,
  sep = "\t",
  quote = FALSE,
  na = ""
)


# ============================================================
# 10. Run EnrichR
#
# IMPORTANT:
# Enrichment is separated by:
#
#     CellType x KD_Gene x Direction
#
# Therefore HCT116 and SW480 are NOT pooled.
# ============================================================

celltype_enrichr_list <- list()


cell_types <- sort(
  unique(
    na.omit(genes_df$Cell)
  )
)


for (cell_type in cell_types) {

  message(
    "\n",
    paste(rep("=", 70), collapse = "")
  )

  message(
    "Cell type: ",
    cell_type
  )

  message(
    paste(rep("=", 70), collapse = "")
  )


  # ----------------------------------------------------------
  # Subset cell type
  # ----------------------------------------------------------

  cell_df <- genes_df[
    Cell == cell_type
  ]


  kd_genes <- sort(
    unique(
      na.omit(cell_df$KD_Gene)
    )
  )


  # ----------------------------------------------------------
  # Loop over KD genes
  # ----------------------------------------------------------

  for (kd_gene in kd_genes) {

    model_df_raw <- cell_df[
      KD_Gene == kd_gene
    ]


    # --------------------------------------------------------
    # Separate upregulated and downregulated genes
    # --------------------------------------------------------

    for (dic in c("Up", "Down")) {

      model_df <- model_df_raw[
        Direction == dic
      ]


      # ------------------------------------------------------
      # Prepare unique valid gene symbols
      # ------------------------------------------------------

      input <- unique(
        na.omit(
          model_df$gene_name
        )
      )

      input <- trimws(input)

      input <- input[
        input != ""
      ]

      input <- sort(
        unique(input)
      )


      message(
        "\nRunning EnrichR: ",
        "Cell=", cell_type,
        "; KD_Gene=", kd_gene,
        "; Direction=", dic,
        "; N_genes=", length(input)
      )


      # ------------------------------------------------------
      # Skip empty gene list
      # ------------------------------------------------------

      if (length(input) == 0) {

        message(
          "Skipping because no valid genes were found."
        )

        next
      }


      # ------------------------------------------------------
      # Run EnrichR
      # ------------------------------------------------------

      enriched <- tryCatch(

        {
          enrichr(
            input,
            dbs_lib
          )
        },

        error = function(e) {

          warning(
            "EnrichR failed for ",
            cell_type,
            " / ",
            kd_gene,
            " / ",
            dic,
            ": ",
            conditionMessage(e)
          )

          return(NULL)
        }
      )


      if (is.null(enriched)) {
        next
      }


      # ------------------------------------------------------
      # Output prefix
      # ------------------------------------------------------

      output_prefix <- paste(
        cell_type,
        kd_gene,
        dic,
        sep = "_"
      )


      # ------------------------------------------------------
      # Filter and write significant enrichment results
      # ------------------------------------------------------

      sig_df <- format_output(
        enriched = enriched,
        output_prefix = output_prefix,
        write_file = TRUE,
        fdr_cutoff = ENRICHR_FDR
      )


      # ------------------------------------------------------
      # Collect significant pathways
      # ------------------------------------------------------

      if (
        !is.null(sig_df) &&
        nrow(sig_df) > 0
      ) {

        sig_df <- as.data.table(sig_df)


        # Add analysis identifiers
        sig_df[
          ,
          Cell := cell_type
        ]

        sig_df[
          ,
          KD_Gene := kd_gene
        ]

        sig_df[
          ,
          Direction := dic
        ]

        sig_df[
          ,
          N_DE_genes := length(input)
        ]


        # ----------------------------------------------------
        # Preserve Dataset if it exists in input
        # ----------------------------------------------------

        if ("Dataset" %in% colnames(model_df)) {

          dataset_values <- unique(
            na.omit(
              as.character(model_df$Dataset)
            )
          )

          if (length(dataset_values) == 1) {

            sig_df[
              ,
              Dataset := dataset_values[1]
            ]

          } else if (length(dataset_values) > 1) {

            sig_df[
              ,
              Dataset := paste(
                sort(dataset_values),
                collapse = ";"
              )
            ]
          }
        }


        # ----------------------------------------------------
        # Store
        # ----------------------------------------------------

        key <- paste(
          cell_type,
          kd_gene,
          dic,
          sep = "__"
        )

        celltype_enrichr_list[[key]] <- sig_df
      }
    }
  }
}


# ============================================================
# 11. Combine all significant enrichment results
# ============================================================

if (length(celltype_enrichr_list) > 0) {

  celltype_enrichr_df <- rbindlist(
    celltype_enrichr_list,
    use.names = TRUE,
    fill = TRUE
  )


  # ==========================================================
  # Put identifiers first
  #
  # Using intersect() prevents setcolorder() from failing
  # when Dataset or version-specific EnrichR columns are absent.
  # ==========================================================

  desired_order <- c(
    "Dataset",
    "Cell",
    "KD_Gene",
    "Direction",
    "N_DE_genes",
    "Source",
    "Term",
    "Overlap",
    "P.value",
    "Adjusted.P.value",
    "Old.P.value",
    "Old.Adjusted.P.value",
    "Odds.Ratio",
    "Combined.Score",
    "Genes"
  )


  existing_order <- intersect(
    desired_order,
    colnames(celltype_enrichr_df)
  )


  remaining_columns <- setdiff(
    colnames(celltype_enrichr_df),
    existing_order
  )


  setcolorder(
    celltype_enrichr_df,
    c(
      existing_order,
      remaining_columns
    )
  )


  # ==========================================================
  # Sort output
  # ==========================================================

  sort_columns <- intersect(
    c(
      "Cell",
      "KD_Gene",
      "Direction",
      "Source",
      "Adjusted.P.value"
    ),
    colnames(celltype_enrichr_df)
  )

  setorderv(
    celltype_enrichr_df,
    sort_columns
  )


  # ==========================================================
  # Write combined result
  # ==========================================================

  fwrite(
    celltype_enrichr_df,
    file = COMBINED_OUTFILE,
    sep = "\t",
    quote = FALSE,
    na = ""
  )


  message(
    "\n",
    paste(rep("=", 70), collapse = "")
  )

  message(
    "Combined significant EnrichR pathways written to:"
  )

  message(
    COMBINED_OUTFILE
  )

  message(
    "Total significant pathway rows: ",
    nrow(celltype_enrichr_df)
  )

  message(
    paste(rep("=", 70), collapse = "")
  )


} else {

  message(
    "\nNo significant EnrichR pathways were found ",
    "at Adjusted.P.value < ",
    ENRICHR_FDR,
    "."
  )
}


# ============================================================
# 12. Finish
# ============================================================

message("\nEnrichR analysis completed.")
