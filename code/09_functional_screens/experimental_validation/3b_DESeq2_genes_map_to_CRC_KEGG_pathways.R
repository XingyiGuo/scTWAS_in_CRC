#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
})

# ============================================================
# User settings
# ============================================================

criteria_keep <- c("High", "Moderate")

REF_FILE <- "/v5000/data/l2_bioinfo1/liq17/LQ_references/pmid_reported_genes/ms_supplementary_tables/CRC_pathway_family_reference_KEGG_hsa_202609.csv"
KEGG_FILE <- "/v5000/data/l2_bioinfo1/liq17/LQ_references/pmid_reported_genes/ms_supplementary_tables/KEGG_Genes_2026.csv"
TARGET_FILE <- "/v5000/data/l2_bioinfo1/liq17/scTWAS/Experiments_res/RNA_seq_DE/HCT116_SW480_7KD_DESeq2_significant_results.csv"

OUT_DIR <- "/v5000/data/l2_bioinfo1/liq17/scTWAS/Experiments_res/RNA_seq_DE/annotate_CRC_pathways"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

OUT_REF_SELECTED <- file.path(
  OUT_DIR, "CRC_pathway_reference_selected.tsv"
)
OUT_GENE_PATHWAY_LONG <- file.path(
  OUT_DIR, "CRC_pathway_KEGG_gene_annotation_long.tsv"
)
OUT_GENE_PATHWAY_3COL <- file.path(
  OUT_DIR, "CRC_pathway_KEGG_gene_annotation_3col.tsv"
)
OUT_GENE_PATHWAY_COLLAPSED <- file.path(
  OUT_DIR, "CRC_pathway_KEGG_gene_annotation_by_gene.tsv"
)
OUT_TARGET_ANNOTATED <- file.path(
  OUT_DIR, "HCT116_SW480_7KD_DESeq2_significant_results_CRCpathway_annotated.csv"
)

# ============================================================
# Step 1. Load CRC pathway reference table and select
#         High + Moderate CRC-relevant pathway families
# ============================================================

# The reference CSV is Windows-1252 encoded (e.g. byte 0x96 for an en dash),
# so reading it directly as UTF-8 can trigger:
#   "input string 1 is invalid UTF-8"
# Read with the correct source encoding and normalize to UTF-8 first.
read_windows1252_fread <- function(file) {
  con <- file(file, open = "r", encoding = "windows-1252")
  on.exit(close(con), add = TRUE)

  txt <- readLines(con, warn = FALSE)
  txt <- enc2utf8(txt)

  fread(text = paste(txt, collapse = "\n"))
}

crc_pathway_reference <- read_windows1252_fread(REF_FILE)

# Extra safeguard: explicitly normalize all character columns to UTF-8.
char_cols <- names(crc_pathway_reference)[
  vapply(crc_pathway_reference, is.character, logical(1))
]

for (j in char_cols) {
  set(
    crc_pathway_reference,
    j = j,
    value = enc2utf8(crc_pathway_reference[[j]])
  )
}

required_ref_cols <- c(
  "CRC_relevance",
  "CRC_pathway_family",
  "Representative_KEGG_hsa_ID"
)

missing_ref_cols <- setdiff(required_ref_cols, names(crc_pathway_reference))

if (length(missing_ref_cols) > 0) {
  stop(
    paste0(
      "Reference table is missing required column(s): ",
      paste(missing_ref_cols, collapse = ", "),
      "\nPlease use CRC_pathway_family_reference_table_with_KEGG_hsa.tsv ",
      "rather than the version without KEGG hsa IDs."
    )
  )
}

# Keep only High and Moderate CRC-relevance categories.
# Prefix matching is used so this remains robust to the descriptive text
# after "High" / "Moderate" and to different dash characters.
crc_pathway_reference_select <- crc_pathway_reference[
  grepl(
    paste0("^(", paste(criteria_keep, collapse = "|"), ")\\b"),
    trimws(CRC_relevance)
  ) &
    !is.na(Representative_KEGG_hsa_ID) &
    nzchar(trimws(Representative_KEGG_hsa_ID))
]

if (nrow(crc_pathway_reference_select) == 0) {
  stop("No CRC pathway rows passed the selected CRC_relevance criteria.")
}

fwrite(
  crc_pathway_reference_select,
  OUT_REF_SELECTED,
  sep = "\t",
  na = ""
)

# ------------------------------------------------------------
# Expand multiple hsa IDs:
# e.g. "hsa04012; hsa04014; hsa04010"
# -> one PathwayID per row
# ------------------------------------------------------------

crc_pathway_reference_hsa <- crc_pathway_reference_select[
  ,
  .(
    PathwayID = trimws(
      unlist(strsplit(
        Representative_KEGG_hsa_ID,
        split = "[;,]"
      ))
    )
  ),
  by = .(
    CRC_relevance,
    CRC_pathway_family
  )
]

crc_pathway_reference_hsa <- unique(
  crc_pathway_reference_hsa[
    !is.na(PathwayID) &
      nzchar(PathwayID)
  ]
)

# ============================================================
# Step 2. Retrieve genes belonging to selected KEGG pathways
# ============================================================

crc_pathway_KEGG_genes <- fread(KEGG_FILE)

required_kegg_cols <- c(
  "PathwayID",
  "PathwayName",
  "All_Genes"
)

missing_kegg_cols <- setdiff(
  required_kegg_cols,
  names(crc_pathway_KEGG_genes)
)

if (length(missing_kegg_cols) > 0) {
  stop(
    paste0(
      "KEGG gene table is missing required column(s): ",
      paste(missing_kegg_cols, collapse = ", ")
    )
  )
}

# Join CRC pathway family/relevance to the KEGG gene table.
crc_kegg_selected <- merge(
  crc_pathway_reference_hsa,
  crc_pathway_KEGG_genes[
    ,
    .(
      PathwayID,
      PathwayName,
      All_Genes
    )
  ],
  by = "PathwayID",
  all.x = TRUE,
  sort = FALSE
)

# Report any reference hsa IDs absent from KEGG_Genes.csv.
missing_hsa <- unique(
  crc_kegg_selected[
    is.na(PathwayName) | !nzchar(trimws(PathwayName)),
    PathwayID
  ]
)

if (length(missing_hsa) > 0) {
  warning(
    paste0(
      "The following hsa IDs were not found in KEGG_Genes.csv: ",
      paste(missing_hsa, collapse = ", ")
    )
  )
}

crc_kegg_selected <- crc_kegg_selected[
  !is.na(All_Genes) &
    nzchar(trimws(All_Genes))
]

# ------------------------------------------------------------
# Expand All_Genes:
# "APC;CTNNB1;AXIN1;..." -> one gene per row
# ------------------------------------------------------------

crc_pathway_KEGG_genes_select_df_annot_long <- crc_kegg_selected[
  ,
  .(
    gene_name = trimws(
      unlist(strsplit(
        All_Genes,
        split = ";",
        fixed = TRUE
      ))
    )
  ),
  by = .(
    PathwayID,
    KEGG_PathwayName = PathwayName,
    pathway = CRC_pathway_family,
    criteria = CRC_relevance
  )
]

crc_pathway_KEGG_genes_select_df_annot_long <-
  crc_pathway_KEGG_genes_select_df_annot_long[
    !is.na(gene_name) &
      nzchar(gene_name)
  ]

# Human gene symbols are normally uppercase; normalize matching.
crc_pathway_KEGG_genes_select_df_annot_long[
  ,
  gene_name := toupper(gene_name)
]

crc_pathway_KEGG_genes_select_df_annot_long <- unique(
  crc_pathway_KEGG_genes_select_df_annot_long
)

# Put High before Moderate for readable output.
crc_pathway_KEGG_genes_select_df_annot_long[
  ,
  criteria_order := fifelse(
    grepl("^High\\b", criteria), 1L,
    fifelse(grepl("^Moderate\\b", criteria), 2L, 3L)
  )
]

setorder(
  crc_pathway_KEGG_genes_select_df_annot_long,
  gene_name,
  criteria_order,
  pathway,
  PathwayID
)

crc_pathway_KEGG_genes_select_df_annot_long[
  ,
  criteria_order := NULL
]

# Full long annotation with hsa ID and KEGG pathway name.
fwrite(
  crc_pathway_KEGG_genes_select_df_annot_long,
  OUT_GENE_PATHWAY_LONG,
  sep = "\t",
  na = ""
)

# ------------------------------------------------------------
# Requested 3-column table:
# [gene_name, pathway, criteria]
#
# A gene can have multiple rows because it can participate in
# multiple CRC-relevant pathways.
# ------------------------------------------------------------

crc_pathway_KEGG_genes_select_df_annot <- unique(
  crc_pathway_KEGG_genes_select_df_annot_long[
    ,
    .(
      gene_name,
      pathway,
      criteria
    )
  ]
)

fwrite(
  crc_pathway_KEGG_genes_select_df_annot,
  OUT_GENE_PATHWAY_3COL,
  sep = "\t",
  na = ""
)

# ============================================================
# Step 2b. Collapse pathway annotation to one row per gene
#          BEFORE joining to DESeq2 results.
#
# This prevents duplicated DESeq2 rows when one gene belongs
# to several pathways.
# ============================================================

crc_pathway_gene_annotation_collapsed <-
  crc_pathway_KEGG_genes_select_df_annot_long[
    ,
    .(
      `CRC known pathways` =
        paste(unique(pathway), collapse = "; "),

      `CRC pathway relevance` =
        paste(unique(criteria), collapse = "; "),

      `CRC KEGG hsa IDs` =
        paste(unique(PathwayID), collapse = "; "),

      `CRC KEGG pathway names` =
        paste(unique(KEGG_PathwayName), collapse = "; "),

      `CRC pathway annotation` =
        paste(
          unique(
            paste0(
              pathway,
              " [",
              criteria,
              "]"
            )
          ),
          collapse = "; "
        )
    ),
    by = gene_name
  ]

fwrite(
  crc_pathway_gene_annotation_collapsed,
  OUT_GENE_PATHWAY_COLLAPSED,
  sep = "\t",
  na = ""
)

# ============================================================
# Step 3. Annotate target DESeq2 genes
# ============================================================

target_genes_df <- fread(TARGET_FILE)

if (!"gene_name" %chin% names(target_genes_df)) {
  stop("TARGET_FILE must contain a column named 'gene_name'.")
}

# Preserve original gene symbols but use an uppercase matching key.
target_genes_df[
  ,
  gene_name_match := toupper(trimws(gene_name))
]

crc_pathway_gene_annotation_collapsed[
  ,
  gene_name_match := gene_name
]

# Remove old annotation columns, if present, so the new annotation
# replaces them cleanly rather than producing .x/.y columns.
annotation_cols <- c(
  "CRC known pathways",
  "CRC pathway relevance",
  "CRC KEGG hsa IDs",
  "CRC KEGG pathway names",
  "CRC pathway annotation"
)

old_annotation_cols <- intersect(
  annotation_cols,
  names(target_genes_df)
)

if (length(old_annotation_cols) > 0) {
  target_genes_df[
    ,
    (old_annotation_cols) := NULL
  ]
}

# Add a row ID to guarantee that the original DESeq2 row order
# is restored after merge().
target_genes_df[
  ,
  .original_row_order := .I
]

target_genes_annotated <- merge(
  target_genes_df,
  crc_pathway_gene_annotation_collapsed[
    ,
    c(
      "gene_name_match",
      annotation_cols
    ),
    with = FALSE
  ],
  by = "gene_name_match",
  all.x = TRUE,
  sort = FALSE
)

setorder(
  target_genes_annotated,
  .original_row_order
)

target_genes_annotated[
  ,
  c(
    "gene_name_match",
    ".original_row_order"
  ) := NULL
]

# Put pathway annotation immediately after gene_name.
annotation_cols_present <- intersect(
  annotation_cols,
  names(target_genes_annotated)
)

other_cols <- setdiff(
  names(target_genes_annotated),
  annotation_cols_present
)

gene_name_pos <- match(
  "gene_name",
  other_cols
)

new_col_order <- append(
  other_cols,
  annotation_cols_present,
  after = gene_name_pos
)

setcolorder(
  target_genes_annotated,
  new_col_order
)

fwrite(
  target_genes_annotated,
  OUT_TARGET_ANNOTATED,
  na = ""
)

# ============================================================
# Summary
# ============================================================

n_target_rows_annotated <- target_genes_annotated[
  !is.na(`CRC known pathways`) &
    nzchar(`CRC known pathways`),
  .N
]

n_target_genes_annotated <- uniqueN(
  target_genes_annotated[
    !is.na(`CRC known pathways`) &
      nzchar(`CRC known pathways`),
    gene_name
  ]
)

cat("\n================ CRC pathway annotation summary ================\n")
cat(
  "Selected reference pathway families:",
  nrow(crc_pathway_reference_select),
  "\n"
)
cat(
  "Unique selected KEGG hsa IDs:",
  uniqueN(crc_pathway_reference_hsa$PathwayID),
  "\n"
)
cat(
  "KEGG hsa IDs matched in KEGG_Genes.csv:",
  uniqueN(crc_kegg_selected$PathwayID),
  "\n"
)
cat(
  "Unique genes in selected CRC-relevant KEGG pathways:",
  uniqueN(crc_pathway_KEGG_genes_select_df_annot$gene_name),
  "\n"
)
cat(
  "Target DESeq2 rows annotated:",
  n_target_rows_annotated,
  "of",
  nrow(target_genes_annotated),
  "\n"
)
cat(
  "Unique target genes annotated:",
  n_target_genes_annotated,
  "\n"
)

cat("\nOutputs:\n")
cat("  ", OUT_REF_SELECTED, "\n")
cat("  ", OUT_GENE_PATHWAY_LONG, "\n")
cat("  ", OUT_GENE_PATHWAY_3COL, "\n")
cat("  ", OUT_GENE_PATHWAY_COLLAPSED, "\n")
cat("  ", OUT_TARGET_ANNOTATED, "\n")
cat("================================================================\n")
