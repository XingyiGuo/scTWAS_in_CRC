library(dplyr)
library(readr)

# ============================================================
# Input files
# ============================================================

NL_FILE <- paste0(
  "/data/l2_bioinfo1/liq17/scTWAS/",
  "TWAS_HTAN_GSE266546_WellcomeSanger_res/",
  "SPrediXcan_peakSNP_association/",
  "Cauchy_Meta_NL_WellcomeSanger_NL_EUR_Multi/",
  "HTAN_NL_WellcomeSanger_NL_7celltypes_peakSNP.EUR_Multi.CauchyMeta.sig.csv"
)

ADSER_FILE <- paste0(
  "/data/l2_bioinfo1/liq17/scTWAS/",
  "TWAS_HTAN_GSE266546_WellcomeSanger_res/",
  "SPrediXcan_peakSNP_association/",
  "VUMC_HTAN_ADSER_EUR_Multi_sig/",
  "VUMC_HTAN_ADSER_9CTs_peakSNP.EUR_Multi.sig.csv"
)


# ============================================================
# Output directory
# ============================================================

OUT_DIR <- paste0(
  "/data/l2_bioinfo1/liq17/scTWAS/",
  "TWAS_HTAN_GSE266546_WellcomeSanger_res/",
  "SPrediXcan_peakSNP_association/"
)

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)


# ============================================================
# 1. ADSER
# ============================================================

adser <- read_csv(
  ADSER_FILE,
  show_col_types = FALSE
)

adser_out <- adser %>%
  mutate(
    Z_score = case_when(
      Sig_resource == "EUR_sig" ~ zscore_EUR,
      Sig_resource == "Multi_sig_EUR_notsig" ~ zscore_Multi,
      TRUE ~ NA_real_
    ),

    P_value = case_when(
      Sig_resource == "EUR_sig" ~ pvalue_EUR,
      Sig_resource == "Multi_sig_EUR_notsig" ~ pvalue_Multi,
      TRUE ~ NA_real_
    )
  ) %>%
  transmute(
    Dataset = "ADSER",
    CellType,
    gene,
    gene_name,
    Sig_resource,
    Z_score,
    P_value
  )


# ============================================================
# 2. Normal colon (NL)
# ============================================================

nl <- read_csv(
  NL_FILE,
  show_col_types = FALSE
)

nl_out <- nl %>%
  mutate(

    # --------------------------------------------------------
    # Select Z score according to ancestry + P-driven dataset
    # --------------------------------------------------------
    Z_score = case_when(

      # EUR-significant
      Sig_resource == "EUR_sig" &
        `P driven_EUR` == "VUMC_HTAN_GSE266546" ~
        zscore_HTAN_NL_EUR,

      Sig_resource == "EUR_sig" &
        `P driven_EUR` == "WellcomeSanger" ~
        zscore_WellcomeSanger_NL_EUR,

      # Multi-significant but EUR non-significant
      Sig_resource == "Multi_sig_EUR_notsig" &
        `P driven_Multi` == "VUMC_HTAN_GSE266546" ~
        zscore_HTAN_NL_Multi,

      Sig_resource == "Multi_sig_EUR_notsig" &
        `P driven_Multi` == "WellcomeSanger" ~
        zscore_WellcomeSanger_NL_Multi,

      TRUE ~ NA_real_
    ),

    # --------------------------------------------------------
    # Meta-analysis P value
    # --------------------------------------------------------
    P_value = case_when(
      Sig_resource == "EUR_sig" ~ META_P_EUR,
      Sig_resource == "Multi_sig_EUR_notsig" ~ META_P_Multi,
      TRUE ~ NA_real_
    )
  ) %>%
  transmute(
    Dataset = "NL",
    CellType,
    gene,
    gene_name,
    Sig_resource,
    Z_score,
    P_value
  )


# ============================================================
# 3. Combine ADSER and NL
# ============================================================

TWAS_results <- bind_rows(
  adser_out,
  nl_out
) %>%
  arrange(
    Dataset,
    CellType,
    P_value
  )


# ============================================================
# 4. QC
# ============================================================

cat("========================================\n")
cat("ADSER rows:", nrow(adser_out), "\n")
cat("NL rows:", nrow(nl_out), "\n")
cat("Combined rows:", nrow(TWAS_results), "\n")
cat("========================================\n")

cat(
  "Missing Z scores:",
  sum(is.na(TWAS_results$Z_score)),
  "\n"
)

cat(
  "Missing P values:",
  sum(is.na(TWAS_results$P_value)),
  "\n"
)


# Check Sig_resource values
cat("\nSig_resource counts:\n")
print(table(TWAS_results$Dataset, TWAS_results$Sig_resource))


# Check NL P-driven datasets
cat("\nNL EUR P-driven resources:\n")
print(table(nl$`P driven_EUR`, useNA = "ifany"))

cat("\nNL Multi P-driven resources:\n")
print(table(nl$`P driven_Multi`, useNA = "ifany"))


# ============================================================
# 5. Save output
# ============================================================

write_csv(
  adser_out,
  file.path(
    OUT_DIR,
    "TWAS_ADSER_collected_results.csv"
  )
)

write_csv(
  nl_out,
  file.path(
    OUT_DIR,
    "TWAS_NL_collected_results.csv"
  )
)

write_csv(
  TWAS_results,
  file.path(
    OUT_DIR,
    "TWAS_ADSER_NL_collected_results.csv"
  )
)


cat("\nFiles written to:\n")

cat(
  file.path(
    OUT_DIR,
    "TWAS_ADSER_collected_results.csv"
  ),
  "\n"
)

cat(
  file.path(
    OUT_DIR,
    "TWAS_NL_collected_results.csv"
  ),
  "\n"
)

cat(
  file.path(
    OUT_DIR,
    "TWAS_ADSER_NL_collected_results.csv"
  ),
  "\n"
)