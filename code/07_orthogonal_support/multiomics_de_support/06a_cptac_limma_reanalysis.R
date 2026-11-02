## 06_cptac_limma_reanalysis.R -------------------------------------------------
## Alternative differential-abundance analysis for the CPTAC-COAD protein arm,
## replacing the per-protein t-test in CPTAC_COAD_TWAS_annotation_streamfix_
## updated.py with the field-standard moderated-t approach.
##
## Why replace the t-test
##   The Python script runs an independent paired t-test per protein and then
##   BH-adjusts across the proteome. That is valid but leaves power on the
##   table: with ~100 pairs and thousands of proteins, the per-protein variance
##   estimate is noisy, and empirical-Bayes shrinkage of those variances
##   towards a common prior (limma) recovers a meaningful amount of power at no
##   cost in type-I error. For TMT data such as CPTAC, DEqMS extends limma by
##   letting the prior depend on how many PSMs or peptides quantified each
##   protein, which is the documented recommendation for this data type.
##
## What this will and will not fix
##   Diagnosed on the existing summary table (see
##   table/DE_key_results/CPTAC_threshold_diagnosis.csv): of the 186 scTWAS
##   genes with protein data, 121 are ALREADY significant at FDR < 0.05. Only
##   27 of those clear |log2FC| > 0.5, and 15 of the 27 also agree in sign with
##   the scTWAS Z. The bottleneck is therefore the effect-size floor, not the
##   test - better statistics will move 121 a little, not 15 a lot. The two
##   changes that actually move the count are the floor (see TREAT below) and
##   protein coverage (186 of 347 genes measured).
##
## Effect-size floor
##   A post-hoc |log2FC| filter on top of an FDR cut has no inferential status:
##   it discards genes without testing anything about them. limma::treat tests
##   the interval null H0: |log2FC| <= lfc directly, so the floor becomes part
##   of the inference and the resulting FDR is the FDR of the thresholded
##   hypothesis. Use TREAT if the floor is kept. Note that protein log-ratios
##   are compressed relative to RNA in this data - median |log2FC| 0.25 against
##   0.49 for TCGA - so an RNA-derived floor of 0.5 is not the same stringency
##   on protein. The quantile-matched protein equivalent of the RNA 0.5 floor
##   is 0.25.
##
## Input: a gene (or protein) by sample matrix and a sample table. This script
##   does NOT download CPTAC - point MATRIX_RDS / SAMPLES_CSV at the matrix the
##   Python script already builds in memory, e.g. by adding
##       prot.to_csv(OUT_DIR / "CPTAC_COAD_protein_matrix.csv")
##   right after the gene-level collapse in that script.
##
## Output
##   table/DE_key_results/CPTAC_COAD_protein_DE_limma.csv
##     gene_name, CPTAC_log2FC, CPTAC_P, CPTAC_FDR  (same columns the
##     downstream annotation step expects), plus the TREAT columns and the
##     DEqMS columns when PSM counts are supplied.
## -----------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(data.table)
  library(limma)
})

PROJ <- "W:/VUMC_projects/scTWAS/VUMC_HTAN_scTWAS_genes_DE_support"
DIR_KEY <- file.path(PROJ, "table", "DE_key_results")
dir.create(DIR_KEY, recursive = TRUE, showWarnings = FALSE)

## --- configuration -----------------------------------------------------------
## MATRIX: rows = genes/proteins, columns = samples, values = log2 abundance
##         (CPTAC umich proteomics is already log2, so no transform is applied)
## SAMPLES: one row per column of MATRIX, with columns
##         sample, patient, group  where group is "tumor" or "normal"
MATRIX_CSV  <- file.path(PROJ, "data", "CPTAC_COAD_protein_matrix.csv")
SAMPLES_CSV <- file.path(PROJ, "data", "CPTAC_COAD_samples.csv")
## optional: gene_name, psm_count  (or peptide count) for DEqMS
PSM_CSV     <- file.path(PROJ, "data", "CPTAC_COAD_psm_counts.csv")

LFC_FLOOR  <- 0.25     # protein-scale equivalent of the RNA 0.5 floor
FDR_CUT    <- 0.05
MIN_PRESENT <- 0.5     # keep a protein quantified in >= this fraction per group

## =============================================================================
## 1. run the analysis
## =============================================================================

#' Moderated-t differential abundance for a paired tumour / normal design.
#'
#' Pairing is handled with duplicateCorrelation rather than by putting patient
#' in the design matrix: CPTAC has unmatched tumours as well as matched pairs,
#' and a patient factor would discard every unpaired sample. duplicateCorrelation
#' estimates one within-patient correlation and uses it as a block structure, so
#' paired and unpaired samples both contribute.
#'
#' @param mat  numeric matrix, genes x samples, log2 scale
#' @param samples data.table with sample / patient / group
#' @param psm optional data.table with gene_name / psm_count -> switches to DEqMS
#' @return data.table with one row per gene
cptac_limma_de <- function(mat, samples, psm = NULL,
                           lfc_floor = LFC_FLOOR, min_present = MIN_PRESENT) {

  stopifnot(all(c("sample", "patient", "group") %in% names(samples)))
  samples <- as.data.table(samples)[match(colnames(mat), sample)]
  stopifnot(!anyNA(samples$sample), all(samples$group %chin% c("tumor", "normal")))

  ## presence filter, applied per group so a protein absent from one arm is
  ## dropped rather than contributing a one-sided comparison
  keep <- apply(mat, 1, function(v) {
    ok <- !is.na(v)
    mean(ok[samples$group == "tumor"])  >= min_present &&
      mean(ok[samples$group == "normal"]) >= min_present
  })
  message(sprintf("proteins: %d input, %d retained after presence filter",
                  nrow(mat), sum(keep)))
  mat <- mat[keep, , drop = FALSE]

  grp <- factor(samples$group, levels = c("normal", "tumor"))
  design <- model.matrix(~ grp)
  colnames(design) <- c("Intercept", "tumor_vs_normal")

  ## within-patient correlation, so matched pairs are not treated as independent
  corfit <- duplicateCorrelation(mat, design, block = samples$patient)
  message(sprintf("consensus within-patient correlation: %.3f",
                  corfit$consensus.correlation))

  fit <- lmFit(mat, design, block = samples$patient,
               correlation = corfit$consensus.correlation)

  ## DEqMS when PSM counts are available: the variance prior is allowed to
  ## depend on the number of spectra quantifying each protein
  used_deqms <- FALSE
  if (!is.null(psm) && requireNamespace("DEqMS", quietly = TRUE)) {
    psm <- as.data.table(psm)
    fit$count <- psm[match(rownames(mat), gene_name), psm_count]
    if (!anyNA(fit$count)) {
      fit <- DEqMS::spectraCounteBayes(eBayes(fit))
      used_deqms <- TRUE
      message("variance model: DEqMS spectraCounteBayes")
    } else {
      message("PSM counts incomplete; falling back to limma eBayes")
    }
  }
  if (!used_deqms) {
    fit <- eBayes(fit, trend = TRUE, robust = TRUE)
    message("variance model: limma eBayes (trend, robust)")
  }

  res <- as.data.table(topTable(fit, coef = "tumor_vs_normal",
                                number = Inf, sort.by = "none"),
                       keep.rownames = "gene_name")
  setnames(res, c("logFC", "P.Value", "adj.P.Val"),
           c("CPTAC_log2FC", "CPTAC_P", "CPTAC_FDR"), skip_absent = TRUE)

  if (used_deqms) {
    dq <- as.data.table(DEqMS::outputResult(fit, coef_col = "tumor_vs_normal"),
                        keep.rownames = "gene_name")
    res <- merge(res, dq[, .(gene_name, sca.P.Value, sca.adj.pval, count)],
                 by = "gene_name", all.x = TRUE)
    setnames(res, c("sca.P.Value", "sca.adj.pval"),
             c("CPTAC_P_deqms", "CPTAC_FDR_deqms"))
  }

  ## TREAT: the effect-size floor as a tested hypothesis rather than a filter
  tfit <- treat(fit, lfc = lfc_floor, trend = TRUE, robust = TRUE)
  tt <- as.data.table(topTreat(tfit, coef = "tumor_vs_normal",
                               number = Inf, sort.by = "none"),
                      keep.rownames = "gene_name")
  res <- merge(res,
               tt[, .(gene_name, CPTAC_P_treat = P.Value,
                      CPTAC_FDR_treat = adj.P.Val)],
               by = "gene_name", all.x = TRUE)
  res[, lfc_floor_used := lfc_floor]
  res[]
}

## =============================================================================
## 2. drive it
## =============================================================================

if (file.exists(MATRIX_CSV) && file.exists(SAMPLES_CSV)) {

  mt <- fread(MATRIX_CSV)
  gcol <- names(mt)[1]
  mat <- as.matrix(mt[, -1, with = FALSE])
  rownames(mat) <- mt[[gcol]]

  samples <- fread(SAMPLES_CSV)
  psm <- if (file.exists(PSM_CSV)) fread(PSM_CSV) else NULL

  res <- cptac_limma_de(mat, samples, psm)
  fwrite(res, file.path(DIR_KEY, "CPTAC_COAD_protein_DE_limma.csv"))

  cat("\n== proteins passing each rule ==\n")
  print(res[, .(
    FDR_only          = sum(CPTAC_FDR < FDR_CUT, na.rm = TRUE),
    FDR_and_lfc_post  = sum(CPTAC_FDR < FDR_CUT &
                              abs(CPTAC_log2FC) > LFC_FLOOR, na.rm = TRUE),
    TREAT             = sum(CPTAC_FDR_treat < FDR_CUT, na.rm = TRUE))])
  cat("\nwrote ", file.path(DIR_KEY, "CPTAC_COAD_protein_DE_limma.csv"), "\n", sep = "")

} else {

  ## Smoke test on a synthetic matrix shaped like CPTAC-COAD, so the code path
  ## above is known to run before it is pointed at the real data.
  message("Protein matrix not found; running the smoke test instead.\n  expected: ",
          MATRIX_CSV)
  set.seed(1)
  n_gene <- 800; n_pair <- 90; n_extra_tumor <- 20
  pat <- sprintf("C3L-%04d", seq_len(n_pair + n_extra_tumor))
  samples <- rbind(
    data.table(sample = paste0(pat[seq_len(n_pair)], ".N"),
               patient = pat[seq_len(n_pair)], group = "normal"),
    data.table(sample = pat, patient = pat, group = "tumor"))
  base <- matrix(rnorm(n_gene * nrow(samples), 0, 0.45), n_gene)
  rownames(base) <- sprintf("GENE%04d", seq_len(n_gene))
  colnames(base) <- samples$sample
  ## 120 genes genuinely shifted in tumour, by a protein-scale amount
  shifted <- sample(n_gene, 120)
  base[shifted, samples$group == "tumor"] <-
    base[shifted, samples$group == "tumor"] + 0.4
  base[sample(length(base), 0.12 * length(base))] <- NA   # TMT missingness

  res <- cptac_limma_de(base, samples,
                        psm = data.table(gene_name = rownames(base),
                                         psm_count = sample(2:40, n_gene, TRUE)))
  cat("\n== smoke test ==\n")
  cat("genes returned:", nrow(res), "\n")
  print(res[, .(FDR_sig = sum(CPTAC_FDR < FDR_CUT, na.rm = TRUE),
                TREAT_sig = sum(CPTAC_FDR_treat < FDR_CUT, na.rm = TRUE),
                recovered_of_120 = sum(CPTAC_FDR < FDR_CUT &
                                         gene_name %chin% rownames(base)[shifted],
                                       na.rm = TRUE))])
  cat("\nTREAT is stricter than the FDR-only rule, as expected:",
      res[, sum(CPTAC_FDR_treat < FDR_CUT, na.rm = TRUE) <=
            sum(CPTAC_FDR < FDR_CUT, na.rm = TRUE)], "\n")
}
