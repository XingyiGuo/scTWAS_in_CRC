#!/usr/bin/env Rscript

# ============================================================
# Annotate EnrichR pathways for CRC relevance and propagate
# High / Moderate CRC-relevant pathway membership to DE genes
#
# Inputs:
#   1) DESeq2_sig_EnrichR_pathways.tsv
#   2) CRC_pathway_family_reference_KEGG_hsa_202609.csv
#   3) HCT116_SW480_7KD_DESeq2_significant_results.csv
#
# Matching is performed at:
#   Cell x KD_Gene x Direction
#
# Pathway annotation strategy:
#   1) Exact normalized pathway-name match to KEGG reference
#      (if duplicated in reference, highest CRC relevance wins)
#   2) Conservative cross-database family rules derived from
#      the CRC pathway-family reference, for Reactome/Hallmark
#   3) Remaining unmatched pathways are assigned Low /
#      "Other uncategorized", consistent with the reference
#
# Outputs:
#   A) DESeq2_sig_EnrichR_pathways_CRC_annotated.tsv
#   B) DESeq2_CRC_HighModerate_pathway_gene_membership.tsv
#   C) HCT116_SW480_7KD_DESeq2_significant_results_CRC_pathway_annotated.csv
#   D) DESeq2_CRC_pathway_gene_summary_by_comparison.tsv
# ============================================================


suppressPackageStartupMessages({
  library(data.table)
})


# ============================================================
# 1. Paths
# ============================================================

BASE_DIR <- "/data/l2_bioinfo1/liq17/scTWAS/Experiments_res/RNA_seq_DE"
OUT_DIR <- paste0(BASE_DIR, "/annotate_CRC_pathways/")

ENRICHR_FILE <- file.path(
  BASE_DIR,
  "DESeq2_sig_EnrichR_pathways.tsv"
)

CRC_REF_FILE <- file.path(
  "/v5000/data/l2_bioinfo1/liq17/LQ_references/pmid_reported_genes/ms_supplementary_tables/CRC_pathway_family_reference_KEGG_hsa_202609.csv"
)

DE_FILE <- file.path(
  BASE_DIR,
  "HCT116_SW480_7KD_DESeq2_significant_results.csv"
)

PATHWAY_OUTFILE <- file.path(
  OUT_DIR,
  "DESeq2_sig_EnrichR_pathways_CRC_annotated.tsv"
)

MEMBERSHIP_OUTFILE <- file.path(
  OUT_DIR,
  "DESeq2_CRC_HighModerate_pathway_gene_membership.tsv"
)

DE_OUTFILE <- file.path(
  OUT_DIR,
  "HCT116_SW480_7KD_DESeq2_significant_results_CRC_pathway_annotated.csv"
)

SUMMARY_OUTFILE <- file.path(
  OUT_DIR,
  "DESeq2_CRC_pathway_gene_summary_by_comparison.tsv"
)


# ============================================================
# 2. Helper functions
# ============================================================

normalize_pathway <- function(x) {

  x <- as.character(x)
  x[is.na(x)] <- ""

  x <- tolower(trimws(x))

  # Normalize common spelling / punctuation differences
  x <- gsub("signalling", "signaling", x, fixed = TRUE)
  x <- gsub("&", " and ", x, fixed = TRUE)

  # Remove KEGG-style hsa IDs if they appear in names
  x <- gsub("\\bhsa[0-9]{5}\\b", " ", x, perl = TRUE)

  # Keep only letters/numbers as comparison tokens
  x <- gsub("[^a-z0-9]+", " ", x, perl = TRUE)
  x <- gsub("\\s+", " ", x, perl = TRUE)

  trimws(x)
}


simple_relevance <- function(x) {

  x <- as.character(x)

  out <- rep(NA_character_, length(x))

  out[grepl("^High", x, ignore.case = TRUE)] <- "High"
  out[grepl("^Moderate", x, ignore.case = TRUE)] <- "Moderate"
  out[grepl("^Low", x, ignore.case = TRUE)] <- "Low"

  out
}


relevance_rank <- function(x) {

  x <- simple_relevance(x)

  out <- rep(0L, length(x))
  out[x == "Low"] <- 1L
  out[x == "Moderate"] <- 2L
  out[x == "High"] <- 3L

  out
}


# ------------------------------------------------------------
# Cross-database annotation rules
#
# These rules translate Reactome / Hallmark wording into the
# pathway families represented in the supplied CRC reference.
#
# Specific clearly non-CRC contexts are checked before broader
# signaling patterns to avoid false positive "High" assignments.
# ------------------------------------------------------------

classify_by_reference_family <- function(term) {

  x <- tolower(as.character(term))

  hit <- function(pattern) {
    grepl(pattern, x, ignore.case = TRUE, perl = TRUE)
  }


  # ==========================================================
  # LOW: specific non-CRC / indirect contexts
  # ==========================================================

  if (hit(
    "sars|viral|virus|infection|hepatitis|toxoplas|amoeb|chagas|papillomavirus|pathogen"
  )) {
    return(c(
      CRC_relevance_simple = "Low",
      CRC_relevance_detail = "Low indirect-nonspecific",
      CRC_pathway_family = "Infection pathogen-related",
      CRC_mapping_method = "Reference_family_rule",
      CRC_mapping_confidence = "Medium"
    ))
  }


  if (hit(
    "neurodegener|neuronal|neurotransmitter|synap|nmda|kainate|glutamate|dopamin|serotonin|norepinephrine|acetylcholine|neurexin|neuroligin|long[- ]term potentiation|long[- ]term depression|addiction|opioid|darpp|neurite|ncam|netrin|l1cam|nrcam"
  )) {
    return(c(
      CRC_relevance_simple = "Low",
      CRC_relevance_detail = "Low indirect-nonspecific",
      CRC_pathway_family = "Neurologic neurodegenerative pathway",
      CRC_mapping_method = "Reference_family_rule",
      CRC_mapping_confidence = "Medium"
    ))
  }


  if (hit(
    "cardiac muscle|cardiomyocyte|adrenergic signaling in cardiomyocytes"
  )) {
    return(c(
      CRC_relevance_simple = "Low",
      CRC_relevance_detail = "Low indirect-nonspecific",
      CRC_pathway_family = "Cardiovascular developmental pathway",
      CRC_mapping_method = "Reference_family_rule",
      CRC_mapping_confidence = "Medium"
    ))
  }


  if (hit(
    "notch|hedgehog|pluripotency|axon guidance|myogenesis"
  )) {
    return(c(
      CRC_relevance_simple = "Low",
      CRC_relevance_detail = "Low indirect-nonspecific",
      CRC_pathway_family = "Developmental differentiation pathway",
      CRC_mapping_method = "Reference_family_rule",
      CRC_mapping_confidence = "Medium"
    ))
  }


  if (hit(
    "^diseases? of |diseases associated with|generic disease"
  )) {
    return(c(
      CRC_relevance_simple = "Low",
      CRC_relevance_detail = "Low indirect-nonspecific",
      CRC_pathway_family = "Generic disease pathway",
      CRC_mapping_method = "Reference_family_rule",
      CRC_mapping_confidence = "Medium"
    ))
  }


  if (hit(
    "small cell lung cancer|lung cancer|breast cancer|prostate cancer|melanoma|leukemia"
  ) && !hit("colorectal|colon cancer|rectal cancer")) {
    return(c(
      CRC_relevance_simple = "Low",
      CRC_relevance_detail = "Low indirect-nonspecific",
      CRC_pathway_family = "Other cancer-specific pathway",
      CRC_mapping_method = "Reference_family_rule",
      CRC_mapping_confidence = "Medium"
    ))
  }


  # ==========================================================
  # HIGH: canonical CRC pathways
  # ==========================================================

  if (hit(
    "colorectal cancer|colon cancer|rectal cancer"
  )) {
    return(c(
      CRC_relevance_simple = "High",
      CRC_relevance_detail = "High-canonical direct",
      CRC_pathway_family = "CRC-specific composite pathway",
      CRC_mapping_method = "Reference_family_rule",
      CRC_mapping_confidence = "Medium"
    ))
  }


  if (hit(
    "wnt|beta[- ]?catenin|ctnnb1|apc truncation|apc mutant|signaling by apc|axin|destruction complex|gsk3beta mutant|amer1 mutant"
  )) {
    return(c(
      CRC_relevance_simple = "High",
      CRC_relevance_detail = "High-canonical direct",
      CRC_pathway_family = "WNT beta-catenin",
      CRC_mapping_method = "Reference_family_rule",
      CRC_mapping_confidence = "Medium"
    ))
  }


  if (hit(
    "erb[bB]|egfr|ras signaling|ras activation|raf activation|raf map|mapk|map kinase|signaling by met|met activates|met promotes|signaling by ntrk|signaling by pdgf|neurotrophin signaling"
  )) {
    return(c(
      CRC_relevance_simple = "High",
      CRC_relevance_detail = "High-canonical direct",
      CRC_pathway_family = "RTK-RAS-MAPK",
      CRC_mapping_method = "Reference_family_rule",
      CRC_mapping_confidence = "Medium"
    ))
  }


  if (hit(
    "pi3k|akt|mtor|mTORC1"
  )) {
    return(c(
      CRC_relevance_simple = "High",
      CRC_relevance_detail = "High-canonical direct",
      CRC_pathway_family = "PI3K-AKT-mTOR",
      CRC_mapping_method = "Reference_family_rule",
      CRC_mapping_confidence = "Medium"
    ))
  }


  if (hit(
    "tgf[- ]?beta|smad"
  )) {
    return(c(
      CRC_relevance_simple = "High",
      CRC_relevance_detail = "High-canonical direct",
      CRC_pathway_family = "TGF-beta SMAD",
      CRC_mapping_method = "Reference_family_rule",
      CRC_mapping_confidence = "Medium"
    ))
  }


  if (hit(
    "tp53|p53"
  )) {
    return(c(
      CRC_relevance_simple = "High",
      CRC_relevance_detail = "High-canonical direct",
      CRC_pathway_family = "TP53 DNA-damage checkpoint",
      CRC_mapping_method = "Reference_family_rule",
      CRC_mapping_confidence = "Medium"
    ))
  }


  if (hit(
    "mismatch repair|msh2|msh3|msh6|mutsalpha|mutsbeta"
  )) {
    return(c(
      CRC_relevance_simple = "High",
      CRC_relevance_detail = "High-canonical direct",
      CRC_pathway_family = "Mismatch repair MSI",
      CRC_mapping_method = "Reference_family_rule",
      CRC_mapping_confidence = "Medium"
    ))
  }


  if (hit(
    "cell cycle|mitotic|g1 s|g1/s|g2 m|g2/m|cyclin|e2f|spindle checkpoint|sister chromatid|kinetochore|orc complex|pre-replicative|pre-replication|s phase"
  )) {
    return(c(
      CRC_relevance_simple = "High",
      CRC_relevance_detail = "High-canonical direct",
      CRC_pathway_family = "Cell cycle proliferation",
      CRC_mapping_method = "Reference_family_rule",
      CRC_mapping_confidence = "Medium"
    ))
  }


  # ==========================================================
  # MODERATE: established CRC-associated pathways
  # ==========================================================

  if (hit(
    "hippo|yap|taz"
  )) {
    return(c(
      CRC_relevance_simple = "Moderate",
      CRC_relevance_detail = "Moderate-established CRC-associated",
      CRC_pathway_family = "Hippo YAP-TAZ",
      CRC_mapping_method = "Reference_family_rule",
      CRC_mapping_confidence = "Medium"
    ))
  }


  if (hit(
    "interferon|interleukin|cytokine|jak[- ]?stat|nf[- ]?kappa|tnf|tlr|myd88|innate immune|ctla4|cd28|immunoregulatory|immune inflammatory"
  )) {
    return(c(
      CRC_relevance_simple = "Moderate",
      CRC_relevance_detail = "Moderate-established CRC-associated",
      CRC_pathway_family = "Immune inflammatory signaling",
      CRC_mapping_method = "Reference_family_rule",
      CRC_mapping_confidence = "Medium"
    ))
  }


  if (hit(
    "epithelial mesenchymal|emt|extracellular matrix|ecm|focal adhesion|adherens|tight junction|integrin|laminin|collagen|cadherin|cell[- ]cell junction|cell junction|cell-extracellular matrix|anchoring fibril"
  )) {
    return(c(
      CRC_relevance_simple = "Moderate",
      CRC_relevance_detail = "Moderate-established CRC-associated",
      CRC_pathway_family = "EMT / ECM adhesion migration",
      CRC_mapping_method = "Reference_family_rule",
      CRC_mapping_confidence = "Medium"
    ))
  }


  if (hit(
    "vegf|vascular endothelial|angiogenesis|hif[- ]?1"
  )) {
    return(c(
      CRC_relevance_simple = "Moderate",
      CRC_relevance_detail = "Moderate-established CRC-associated",
      CRC_pathway_family = "Angiogenesis VEGF",
      CRC_mapping_method = "Reference_family_rule",
      CRC_mapping_confidence = "Medium"
    ))
  }


  if (hit(
    "dna replication|replication stress|base excision|nucleotide excision|homologous recombination|non-homologous|fanconi|brca1|brca2|palb2|rad51|activation of atr"
  )) {
    return(c(
      CRC_relevance_simple = "Moderate",
      CRC_relevance_detail = "Moderate-established CRC-associated",
      CRC_pathway_family = "DNA replication repair (non-MMR)",
      CRC_mapping_method = "Reference_family_rule",
      CRC_mapping_confidence = "Medium"
    ))
  }


  if (hit(
    "apoptosis|ferroptosis|necroptosis|programmed cell death"
  )) {
    return(c(
      CRC_relevance_simple = "Moderate",
      CRC_relevance_detail = "Moderate-established CRC-associated",
      CRC_pathway_family = "Regulated cell death",
      CRC_mapping_method = "Reference_family_rule",
      CRC_mapping_confidence = "Medium"
    ))
  }


  if (hit(
    "senescence"
  )) {
    return(c(
      CRC_relevance_simple = "Moderate",
      CRC_relevance_detail = "Moderate-established CRC-associated",
      CRC_pathway_family = "Cellular senescence",
      CRC_mapping_method = "Reference_family_rule",
      CRC_mapping_confidence = "Medium"
    ))
  }


  if (hit(
    "oxidative phosphorylation|respiratory electron transport|aerobic respiration|glycolysis|gluconeogenesis|citrate cycle|tca cycle|pentose phosphate|glutathione|central carbon|pyruvate metabolism|glucose metabolism"
  )) {
    return(c(
      CRC_relevance_simple = "Moderate",
      CRC_relevance_detail = "Moderate-established CRC-associated",
      CRC_pathway_family = "Metabolism redox mitochondrial",
      CRC_mapping_method = "Reference_family_rule",
      CRC_mapping_confidence = "Medium"
    ))
  }


  if (hit(
    "protein processing in endoplasmic reticulum|unfolded protein|upr|ubiquitin mediated|proteasome|autophagy|protein folding|chaperonin|proteostasis"
  )) {
    return(c(
      CRC_relevance_simple = "Moderate",
      CRC_relevance_detail = "Moderate-established CRC-associated",
      CRC_pathway_family = "Stress proteostasis",
      CRC_mapping_method = "Reference_family_rule",
      CRC_mapping_confidence = "Medium"
    ))
  }


  if (hit(
    "transcriptional misregulation in cancer|pathways in cancer"
  )) {
    return(c(
      CRC_relevance_simple = "Moderate",
      CRC_relevance_detail = "Moderate-established CRC-associated",
      CRC_pathway_family = "Shared oncogenic transcription",
      CRC_mapping_method = "Reference_family_rule",
      CRC_mapping_confidence = "Medium"
    ))
  }


  if (hit(
    "signaling by bmp|bmp signaling"
  )) {
    return(c(
      CRC_relevance_simple = "Moderate",
      CRC_relevance_detail = "Moderate-established CRC-associated",
      CRC_pathway_family = "EMT stromal program",
      CRC_mapping_method = "Reference_family_rule",
      CRC_mapping_confidence = "Medium"
    ))
  }


  # ==========================================================
  # LOW: indirect / nonspecific pathway families
  # ==========================================================

  if (hit(
    "xenobiotic|cytochrome p450|drug metabolism"
  )) {
    return(c(
      CRC_relevance_simple = "Low",
      CRC_relevance_detail = "Low indirect-nonspecific",
      CRC_pathway_family = "Drug xenobiotic metabolism",
      CRC_mapping_method = "Reference_family_rule",
      CRC_mapping_confidence = "Medium"
    ))
  }


  if (hit(
    "ampk|foxo"
  )) {
    return(c(
      CRC_relevance_simple = "Low",
      CRC_relevance_detail = "Low indirect-nonspecific",
      CRC_pathway_family = "General cellular metabolic regulation",
      CRC_mapping_method = "Reference_family_rule",
      CRC_mapping_confidence = "Medium"
    ))
  }


  if (hit(
    "insulin|glucagon|steroid|cholesterol|vitamin|biotin|glycogen|carbohydrate|peroxisome|heme metabolism|sphingolipid|starch and sucrose"
  )) {
    return(c(
      CRC_relevance_simple = "Low",
      CRC_relevance_detail = "Low indirect-nonspecific",
      CRC_pathway_family = "General metabolism endocrine",
      CRC_mapping_method = "Reference_family_rule",
      CRC_mapping_confidence = "Medium"
    ))
  }


  if (hit(
    "endocytosis|lysosome|gap junction"
  )) {
    return(c(
      CRC_relevance_simple = "Low",
      CRC_relevance_detail = "Low indirect-nonspecific",
      CRC_pathway_family = "General cellular process",
      CRC_mapping_method = "Reference_family_rule",
      CRC_mapping_confidence = "Medium"
    ))
  }


  if (hit(
    "uv response|generic stress"
  )) {
    return(c(
      CRC_relevance_simple = "Low",
      CRC_relevance_detail = "Low indirect-nonspecific",
      CRC_pathway_family = "Generic stress response",
      CRC_mapping_method = "Reference_family_rule",
      CRC_mapping_confidence = "Medium"
    ))
  }


  # ==========================================================
  # Fallback
  #
  # The supplied CRC reference contains a Low category called
  # "Other uncategorized". Therefore terms not supported by a
  # High/Moderate/reference-family rule are retained as Low,
  # but flagged with Low mapping confidence.
  # ==========================================================

  c(
    CRC_relevance_simple = "Low",
    CRC_relevance_detail = "Low indirect-nonspecific",
    CRC_pathway_family = "Other uncategorized",
    CRC_mapping_method = "Fallback_other_uncategorized",
    CRC_mapping_confidence = "Low"
  )
}


collapse_unique <- function(x, sep = " | ") {

  x <- unique(as.character(x))
  x <- x[!is.na(x) & trimws(x) != ""]

  if (length(x) == 0) {
    return("")
  }

  paste(sort(x), collapse = sep)
}


# ============================================================
# 3. Read input files
# ============================================================

message("Reading input files...")

enrich_df <- fread(
  ENRICHR_FILE,
  sep = "\t",
  data.table = TRUE
)

crc_ref <- fread(
  CRC_REF_FILE,
  data.table = TRUE,
  encoding = "UTF-8"
)

de_df <- fread(
  DE_FILE,
  data.table = TRUE
)


# ============================================================
# 4. Check required columns
# ============================================================

required_enrich <- c(
  "Cell",
  "KD_Gene",
  "Direction",
  "Source",
  "Term",
  "Genes"
)

required_ref <- c(
  "CRC_relevance",
  "CRC_pathway_family",
  "Representative_KEGG_pathway_name"
)

required_de <- c(
  "Cell",
  "KD_Gene",
  "Direction",
  "gene_name"
)


missing_enrich <- setdiff(
  required_enrich,
  names(enrich_df)
)

missing_ref <- setdiff(
  required_ref,
  names(crc_ref)
)

missing_de <- setdiff(
  required_de,
  names(de_df)
)


if (length(missing_enrich) > 0) {
  stop(
    "Missing columns in EnrichR file: ",
    paste(missing_enrich, collapse = ", ")
  )
}

if (length(missing_ref) > 0) {
  stop(
    "Missing columns in CRC reference file: ",
    paste(missing_ref, collapse = ", ")
  )
}

if (length(missing_de) > 0) {
  stop(
    "Missing columns in DESeq2 file: ",
    paste(missing_de, collapse = ", ")
  )
}


# ============================================================
# 5. Build exact pathway-name mapping from CRC reference
# ============================================================

message("Building exact pathway-name reference map...")


# Add missing optional reference fields if necessary
optional_ref_cols <- c(
  "Reference_key",
  "Annotation_rationale",
  "Representative_KEGG_hsa_ID"
)

for (cc in optional_ref_cols) {
  if (!cc %in% names(crc_ref)) {
    crc_ref[, (cc) := NA_character_]
  }
}


# Split pathway names separated by "|"
ref_long <- crc_ref[
  !is.na(Representative_KEGG_pathway_name) &
    trimws(Representative_KEGG_pathway_name) != "",
  {

    paths <- unlist(
      strsplit(
        Representative_KEGG_pathway_name,
        "\\s*\\|\\s*",
        perl = TRUE
      )
    )

    paths <- trimws(paths)
    paths <- paths[paths != ""]

    .(
      Reference_pathway_name = paths
    )
  },
  by = .(
    CRC_relevance,
    CRC_pathway_family,
    Representative_KEGG_hsa_ID,
    Annotation_rationale,
    Reference_key
  )
]


ref_long[
  ,
  pathway_norm := normalize_pathway(
    Reference_pathway_name
  )
]

ref_long[
  ,
  CRC_relevance_simple := simple_relevance(
    CRC_relevance
  )
]

ref_long[
  ,
  relevance_rank := relevance_rank(
    CRC_relevance
  )
]


# If the same KEGG pathway occurs in multiple reference families
# with different relevance levels, use the highest level:
# High > Moderate > Low.
ref_exact <- ref_long[
  pathway_norm != "",
  {

    max_rank <- max(
      relevance_rank,
      na.rm = TRUE
    )

    z <- .SD[
      relevance_rank == max_rank
    ]

    .(
      CRC_relevance_simple_ref = z$CRC_relevance_simple[1],
      CRC_relevance_detail_ref = z$CRC_relevance[1],
      CRC_pathway_family_ref = collapse_unique(
        z$CRC_pathway_family
      ),
      CRC_reference_pathway_ref = collapse_unique(
        z$Reference_pathway_name
      ),
      CRC_reference_key_ref = collapse_unique(
        z$Reference_key
      ),
      CRC_annotation_rationale_ref = collapse_unique(
        z$Annotation_rationale
      )
    )
  },
  by = pathway_norm
]


message(
  "Unique exact pathway names in CRC reference: ",
  nrow(ref_exact)
)


# ============================================================
# 6. Annotate EnrichR pathways
# ============================================================

message("Annotating EnrichR pathways...")


enrich_df[
  ,
  original_row_id := .I
]

enrich_df[
  ,
  pathway_norm := normalize_pathway(Term)
]


# Exact normalized-name match
enrich_annot <- merge(
  enrich_df,
  ref_exact,
  by = "pathway_norm",
  all.x = TRUE,
  sort = FALSE
)


# Initialize final annotation columns
enrich_annot[
  ,
  `:=`(
    CRC_relevance_simple = NA_character_,
    CRC_relevance_detail = NA_character_,
    CRC_pathway_family = NA_character_,
    CRC_mapping_method = NA_character_,
    CRC_mapping_confidence = NA_character_,
    CRC_reference_pathway = NA_character_,
    CRC_reference_key = NA_character_,
    CRC_annotation_rationale = NA_character_
  )
]


# Use exact reference-name matches first
exact_idx <- which(
  !is.na(enrich_annot$CRC_relevance_simple_ref)
)


if (length(exact_idx) > 0) {

  enrich_annot[
    exact_idx,
    `:=`(
      CRC_relevance_simple = CRC_relevance_simple_ref,
      CRC_relevance_detail = CRC_relevance_detail_ref,
      CRC_pathway_family = CRC_pathway_family_ref,
      CRC_mapping_method = "Exact_reference_name",
      CRC_mapping_confidence = "High",
      CRC_reference_pathway = CRC_reference_pathway_ref,
      CRC_reference_key = CRC_reference_key_ref,
      CRC_annotation_rationale = CRC_annotation_rationale_ref
    )
  ]
}


# Apply cross-database family rules only to pathways that did not
# have an exact reference-name match.
rule_idx <- which(
  is.na(enrich_annot$CRC_relevance_simple)
)


if (length(rule_idx) > 0) {

  rule_list <- lapply(
    enrich_annot$Term[rule_idx],
    classify_by_reference_family
  )

  rule_mat <- do.call(
    rbind,
    rule_list
  )

  rule_dt <- as.data.table(
    rule_mat
  )

  enrich_annot[
    rule_idx,
    `:=`(
      CRC_relevance_simple =
        rule_dt$CRC_relevance_simple,

      CRC_relevance_detail =
        rule_dt$CRC_relevance_detail,

      CRC_pathway_family =
        rule_dt$CRC_pathway_family,

      CRC_mapping_method =
        rule_dt$CRC_mapping_method,

      CRC_mapping_confidence =
        rule_dt$CRC_mapping_confidence
    )
  ]
}


# Restore original enrichment-table row order
setorder(
  enrich_annot,
  original_row_id
)


# Remove helper/reference intermediate columns
drop_cols <- intersect(
  c(
    "pathway_norm",
    "original_row_id",
    "CRC_relevance_simple_ref",
    "CRC_relevance_detail_ref",
    "CRC_pathway_family_ref",
    "CRC_reference_pathway_ref",
    "CRC_reference_key_ref",
    "CRC_annotation_rationale_ref"
  ),
  names(enrich_annot)
)

enrich_annot[
  ,
  (drop_cols) := NULL
]


# Put CRC annotation immediately after pathway term
front_cols <- c(
  "Cell",
  "KD_Gene",
  "Direction",
  "N_DE_genes",
  "Source",
  "Term",
  "CRC_relevance_simple",
  "CRC_relevance_detail",
  "CRC_pathway_family",
  "CRC_mapping_method",
  "CRC_mapping_confidence",
  "CRC_reference_pathway",
  "CRC_reference_key",
  "CRC_annotation_rationale"
)

front_cols <- intersect(
  front_cols,
  names(enrich_annot)
)

setcolorder(
  enrich_annot,
  c(
    front_cols,
    setdiff(names(enrich_annot), front_cols)
  )
)


fwrite(
  enrich_annot,
  PATHWAY_OUTFILE,
  sep = "\t",
  quote = FALSE,
  na = ""
)


message(
  "Annotated pathway table written to: ",
  PATHWAY_OUTFILE
)


message("\nCRC relevance among enriched pathway rows:")
print(
  enrich_annot[
    ,
    .N,
    by = CRC_relevance_simple
  ][
    order(
      factor(
        CRC_relevance_simple,
        levels = c("High", "Moderate", "Low")
      )
    )
  ]
)


message("\nMapping methods:")
print(
  enrich_annot[
    ,
    .N,
    by = .(
      CRC_mapping_method,
      CRC_mapping_confidence
    )
  ][
    order(-N)
  ]
)


# ============================================================
# 7. Expand genes from High / Moderate CRC-relevant pathways
# ============================================================

message(
  "\nExpanding genes from High/Moderate CRC-relevant pathways..."
)


hm_pathways <- enrich_annot[
  CRC_relevance_simple %in% c(
    "High",
    "Moderate"
  ) &
    !is.na(Genes) &
    trimws(Genes) != ""
]


if (nrow(hm_pathways) == 0) {

  warning(
    "No High or Moderate CRC-relevant enriched pathways were found."
  )

  pathway_gene_long <- data.table(
    Cell = character(),
    KD_Gene = character(),
    Direction = character(),
    gene_name_norm = character()
  )

} else {

  pathway_gene_long <- hm_pathways[
    ,
    {

      gene_vec <- unlist(
        strsplit(
          as.character(Genes),
          "[;,]",
          perl = TRUE
        )
      )

      gene_vec <- trimws(gene_vec)
      gene_vec <- gene_vec[
        !is.na(gene_vec) &
          gene_vec != ""
      ]

      gene_vec <- unique(gene_vec)

      .(
        gene_name = gene_vec,
        gene_name_norm = toupper(gene_vec)
      )
    },
    by = .(
      Cell,
      KD_Gene,
      Direction,
      Source,
      Term,
      CRC_relevance_simple,
      CRC_relevance_detail,
      CRC_pathway_family,
      CRC_mapping_method,
      CRC_mapping_confidence
    )
  ]


  setorder(
    pathway_gene_long,
    Cell,
    KD_Gene,
    Direction,
    gene_name_norm,
    CRC_relevance_simple,
    Source,
    Term
  )
}


fwrite(
  pathway_gene_long,
  MEMBERSHIP_OUTFILE,
  sep = "\t",
  quote = FALSE,
  na = ""
)


message(
  "High/Moderate pathway-gene membership table written to: ",
  MEMBERSHIP_OUTFILE
)


# ============================================================
# 8. Summarize CRC pathway membership per DE gene
# ============================================================

if (nrow(pathway_gene_long) > 0) {

  gene_crc_summary <- pathway_gene_long[
    ,
    .(
      CRC_High_pathway_enriched = ifelse(
        any(CRC_relevance_simple == "High"),
        "Yes",
        "No"
      ),

      CRC_Moderate_pathway_enriched = ifelse(
        any(CRC_relevance_simple == "Moderate"),
        "Yes",
        "No"
      ),

      N_CRC_High_pathways = uniqueN(
        Term[
          CRC_relevance_simple == "High"
        ]
      ),

      N_CRC_Moderate_pathways = uniqueN(
        Term[
          CRC_relevance_simple == "Moderate"
        ]
      ),

      CRC_High_pathways = collapse_unique(
        paste0(
          Source[
            CRC_relevance_simple == "High"
          ],
          ": ",
          Term[
            CRC_relevance_simple == "High"
          ]
        )
      ),

      CRC_Moderate_pathways = collapse_unique(
        paste0(
          Source[
            CRC_relevance_simple == "Moderate"
          ],
          ": ",
          Term[
            CRC_relevance_simple == "Moderate"
          ]
        )
      ),

      CRC_High_pathway_families = collapse_unique(
        CRC_pathway_family[
          CRC_relevance_simple == "High"
        ]
      ),

      CRC_Moderate_pathway_families = collapse_unique(
        CRC_pathway_family[
          CRC_relevance_simple == "Moderate"
        ]
      )
    ),
    by = .(
      Cell,
      KD_Gene,
      Direction,
      gene_name_norm
    )
  ]


  gene_crc_summary[
    ,
    CRC_High_or_Moderate_pathway_enriched := ifelse(
      CRC_High_pathway_enriched == "Yes" |
        CRC_Moderate_pathway_enriched == "Yes",
      "Yes",
      "No"
    )
  ]


  gene_crc_summary[
    ,
    CRC_max_relevance := fifelse(
      CRC_High_pathway_enriched == "Yes",
      "High",
      fifelse(
        CRC_Moderate_pathway_enriched == "Yes",
        "Moderate",
        "None"
      )
    )
  ]

} else {

  gene_crc_summary <- data.table(
    Cell = character(),
    KD_Gene = character(),
    Direction = character(),
    gene_name_norm = character(),
    CRC_High_pathway_enriched = character(),
    CRC_Moderate_pathway_enriched = character(),
    N_CRC_High_pathways = integer(),
    N_CRC_Moderate_pathways = integer(),
    CRC_High_pathways = character(),
    CRC_Moderate_pathways = character(),
    CRC_High_pathway_families = character(),
    CRC_Moderate_pathway_families = character(),
    CRC_High_or_Moderate_pathway_enriched = character(),
    CRC_max_relevance = character()
  )
}


# ============================================================
# 9. Annotate original DESeq2 gene table
#
# Critical:
# Join on Cell + KD_Gene + Direction + gene_name
# so pathway membership is comparison-specific.
# ============================================================

message(
  "\nAnnotating original DESeq2 significant-gene table..."
)


de_df[
  ,
  original_row_id := .I
]

de_df[
  ,
  gene_name_norm := toupper(
    trimws(
      as.character(gene_name)
    )
  )
]


de_annot <- merge(
  de_df,
  gene_crc_summary,
  by = c(
    "Cell",
    "KD_Gene",
    "Direction",
    "gene_name_norm"
  ),
  all.x = TRUE,
  sort = FALSE
)


# Fill genes not present in any High/Moderate enriched pathway
flag_cols <- c(
  "CRC_High_pathway_enriched",
  "CRC_Moderate_pathway_enriched",
  "CRC_High_or_Moderate_pathway_enriched"
)

for (cc in flag_cols) {
  de_annot[
    is.na(get(cc)),
    (cc) := "No"
  ]
}


count_cols <- c(
  "N_CRC_High_pathways",
  "N_CRC_Moderate_pathways"
)

for (cc in count_cols) {
  de_annot[
    is.na(get(cc)),
    (cc) := 0L
  ]
}


text_cols <- c(
  "CRC_High_pathways",
  "CRC_Moderate_pathways",
  "CRC_High_pathway_families",
  "CRC_Moderate_pathway_families"
)

for (cc in text_cols) {
  de_annot[
    is.na(get(cc)),
    (cc) := ""
  ]
}


de_annot[
  is.na(CRC_max_relevance),
  CRC_max_relevance := "None"
]


# Restore original DESeq2 row order
setorder(
  de_annot,
  original_row_id
)

de_annot[
  ,
  c(
    "original_row_id",
    "gene_name_norm"
  ) := NULL
]


# Put original columns first, CRC annotations after
original_de_cols <- names(de_df)

original_de_cols <- setdiff(
  original_de_cols,
  c(
    "original_row_id",
    "gene_name_norm"
  )
)

crc_gene_cols <- c(
  "CRC_High_or_Moderate_pathway_enriched",
  "CRC_max_relevance",
  "CRC_High_pathway_enriched",
  "CRC_Moderate_pathway_enriched",
  "N_CRC_High_pathways",
  "N_CRC_Moderate_pathways",
  "CRC_High_pathway_families",
  "CRC_Moderate_pathway_families",
  "CRC_High_pathways",
  "CRC_Moderate_pathways"
)

setcolorder(
  de_annot,
  c(
    intersect(
      original_de_cols,
      names(de_annot)
    ),
    intersect(
      crc_gene_cols,
      names(de_annot)
    ),
    setdiff(
      names(de_annot),
      c(
        original_de_cols,
        crc_gene_cols
      )
    )
  )
)


fwrite(
  de_annot,
  DE_OUTFILE,
  sep = ",",
  quote = TRUE,
  na = ""
)


message(
  "Annotated DESeq2 gene table written to: ",
  DE_OUTFILE
)


# ============================================================
# 10. Summary by Cell x KD_Gene x Direction
# ============================================================

comparison_summary <- de_annot[
  ,
  .(
    N_DE_genes = uniqueN(
      gene_name[
        !is.na(gene_name) &
          trimws(gene_name) != ""
      ]
    ),

    N_genes_in_CRC_High_pathways = uniqueN(
      gene_name[
        CRC_High_pathway_enriched == "Yes"
      ]
    ),

    N_genes_in_CRC_Moderate_pathways = uniqueN(
      gene_name[
        CRC_Moderate_pathway_enriched == "Yes"
      ]
    ),

    N_genes_in_CRC_High_or_Moderate_pathways = uniqueN(
      gene_name[
        CRC_High_or_Moderate_pathway_enriched == "Yes"
      ]
    )
  ),
  by = .(
    Cell,
    KD_Gene,
    Direction
  )
]


comparison_summary[
  ,
  Pct_genes_in_CRC_High_or_Moderate_pathways :=
    ifelse(
      N_DE_genes > 0,
      100 *
        N_genes_in_CRC_High_or_Moderate_pathways /
        N_DE_genes,
      NA_real_
    )
]


setorder(
  comparison_summary,
  Cell,
  KD_Gene,
  Direction
)


fwrite(
  comparison_summary,
  SUMMARY_OUTFILE,
  sep = "\t",
  quote = FALSE,
  na = ""
)


message(
  "Comparison-level summary written to: ",
  SUMMARY_OUTFILE
)


# ============================================================
# 11. Console summary
# ============================================================

message("\n============================================================")
message("FINAL SUMMARY")
message("============================================================")

message(
  "EnrichR pathway rows: ",
  nrow(enrich_annot)
)

message(
  "High CRC-relevance pathway rows: ",
  enrich_annot[
    CRC_relevance_simple == "High",
    .N
  ]
)

message(
  "Moderate CRC-relevance pathway rows: ",
  enrich_annot[
    CRC_relevance_simple == "Moderate",
    .N
  ]
)

message(
  "Low CRC-relevance pathway rows: ",
  enrich_annot[
    CRC_relevance_simple == "Low",
    .N
  ]
)

message(
  "DESeq2 rows annotated as belonging to >=1 High/Moderate CRC pathway: ",
  de_annot[
    CRC_High_or_Moderate_pathway_enriched == "Yes",
    .N
  ]
)

message("\nDone.")
