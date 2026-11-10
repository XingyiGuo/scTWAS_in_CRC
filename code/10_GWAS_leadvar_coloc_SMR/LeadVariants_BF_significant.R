# conda activate tftwas_coloc
# ------------------------------------------------------------------------------
# Bonferroni-significant lead-variant x gene associations across the 7 colon
# epithelial cell types, cis and trans counted and corrected SEPARATELY.
#
#   cutoff_cis   = 0.05 / (number of cis   tests in that cell type)
#   cutoff_trans = 0.05 / (number of trans tests in that cell type)
#
# Correction is PER CELL TYPE: each cell type's own test counts set its cutoffs.
# A study-wide correction would divide by the pooled totals instead.
#
# Each results file is several GB, so it is streamed in one pass: every block is
# counted and only rows with p <= PREFILTER are held as significance candidates.
# PREFILTER is checked against both cutoffs and the file re-scanned if it is too
# tight.
#
# Outputs (in BF_OUT_DIR, default the results directory):
#   <CT>_cis_BFsig.csv / <CT>_trans_BFsig.csv     per cell type
#   AllCellTypes_cis_BFsig.csv / AllCellTypes_trans_BFsig.csv
#   AllCellTypes_BF_summary.csv                   counts + cutoffs per cell type
#
# Every output carries a CellTypes column: the semicolon-separated list of cell
# types in which that lead variant x gene pair is significant in the same class
# (CellType stays the per-row cell type; n_CellTypes is the count).
#
# Usage: Rscript LeadVariants_BF_significant.R VUMC_HTAN_GSE266546          # all cell types of that analysis
#        BF_CELL_TYPES=ABS,GOB Rscript LeadVariants_BF_significant.R VUMC_HTAN_GSE266546   # a subset
# ------------------------------------------------------------------------------

suppressPackageStartupMessages(library(data.table))

args <- if (exists("BF_ARGS")) BF_ARGS else commandArgs(trailingOnly = TRUE)

# ------------------------------------------------------------------------------
# Analysis type: WellcomeSanger (normal colon, N=278), VUMC_HTAN_GSE266546
# (normal colon, N=103) or VUMC_Polyps (ADSER). This one switch fixes the
# results directory, the cell-type list, the expression tree and the
# phenotype-BED filename together, so regression results are never paired
# with another cohort's BED.
#
#   Rscript LeadVariants_BF_significant.R [WellcomeSanger|VUMC_HTAN_GSE266546|VUMC_Polyps]
#   BF_ANALYSIS=VUMC_Polyps Rscript LeadVariants_BF_significant.R
# ------------------------------------------------------------------------------

ANALYSIS <- if (length(args) >= 1 && nzchar(args[1]))
    args[1] else Sys.getenv("BF_ANALYSIS", unset = "WellcomeSanger")

if (!ANALYSIS %in% c("WellcomeSanger", "VUMC_HTAN_GSE266546", "VUMC_Polyps"))
    stop("BF_ANALYSIS must be 'WellcomeSanger', 'VUMC_HTAN_GSE266546' or 'VUMC_Polyps', got: ", ANALYSIS)

## Only the analysis type is positional. Cell types follow from it; the old
## "CTs analysis" two-argument form is rejected rather than guessed at.
if (length(args) >= 2)
    stop("usage: Rscript LeadVariants_BF_significant.R <analysis>\n",
         "  cell types are set by the analysis; to subset use BF_CELL_TYPES=ABS,GOB")

first_existing <- function(...) {
    for (p in c(...)) if (!is.na(p) && nzchar(p) && file.exists(p)) return(p)
    NA_character_
}

LOCAL_ROOT <- "W:/VUMC_projects/scTWAS/GWAS_leadvar_coloc_SMR"

if (ANALYSIS == "WellcomeSanger") {

    CTs_all      <- c("ABS", "CT", "EE", "GOB", "STM", "TAC", "TUF")
    RES_DEFAULT  <- "/v5000/data/l2_bioinfo1/liq17/scTWAS/GWAS_leadvar_coloc_SMR_res/WellcomeSanger"
    CLUSTER_ROOT <- "/data/l2_bioinfo1/liq17/scTWAS/data/WellcomeSanger"
    EXPR_DEFAULT <- file.path(CLUSTER_ROOT,
                              "expression/Healthy_TI_then_rectum_278samples_outputs_LauLabel",
                              "pseudobulk_Epithelial_LauLabel_x_individual/processed_expression")
    bed_name     <- function(ct) paste0(ct, ".TPM.log2.QN.INV.PEER_residuals.individual_id.tensorqtl.bed.gz")
    bed_subdirs  <- function(ct) c(paste0("PEER_", ct), ".")

} else if (ANALYSIS == "VUMC_Polyps") {

    CTs_all      <- c("ABS", "CT", "EE", "GOB", "STM", "TAC", "TUF", "ASC", "SSC")
    RES_DEFAULT  <- "/v5000/data/l2_bioinfo1/liq17/scTWAS/GWAS_leadvar_coloc_SMR_res/Polyps"
    CLUSTER_ROOT <- "/data/l2_bioinfo1/liq17/scTWAS/data/VUMC_HTAN"
    EXPR_DEFAULT <- file.path(CLUSTER_ROOT, "expression/ADSER_pseudobulk_expression_for_peer")
    bed_name     <- function(ct) paste0("pseudobulk_counts_", ct,
                                        ".tpmlog2.qc.qn.rmlowcounts.inv.PEER_residuals.genotype_id.tensorqtl.bed.gz")
    bed_subdirs  <- function(ct) c(file.path(ct, paste0("PEER_", ct)), ".")

} else {

    ## VUMC_HTAN_GSE266546: normal colon (NL), N=103, same 7 cell types as
    ## WellcomeSanger. Regression results come from LeadVariants_Profile_lm_complete_VUMC.R
    ## (OUT_DIR = .../GWAS_leadvar_coloc_SMR_res/VUMC_HTAN_GSE266546/); the BED
    ## is the NL tensorQTL phenotype file, same filename convention as that script.
    CTs_all      <- c("ABS", "CT", "EE", "GOB", "STM", "TAC", "TUF")
    RES_DEFAULT  <- "/v5000/data/l2_bioinfo1/liq17/scTWAS/GWAS_leadvar_coloc_SMR_res/VUMC_HTAN_GSE266546"
    CLUSTER_ROOT <- "/data/l2_bioinfo1/liq17/scTWAS/data/VUMC_HTAN_GSE266546"
    EXPR_DEFAULT <- file.path(CLUSTER_ROOT, "expression/pseudobulk_expression_for_peer_NL")
    bed_name     <- function(ct) paste0("pseudobulk_counts_", ct,
                                        ".tpmlog2.qc.qn.rmlowcounts.inv.PEER_residuals.genotype_id.tensorqtl.bed.gz")
    bed_subdirs  <- function(ct) c(file.path(ct, paste0("PEER_", ct)), ".")

}

RES_DIR  <- Sys.getenv("BF_RES_DIR", unset = RES_DEFAULT)
EXPR_DIR <- Sys.getenv("LEADVAR_EXPR_DIR", unset = first_existing(EXPR_DEFAULT, LOCAL_ROOT))

## Optional subset: BF_CELL_TYPES=ABS,GOB (default and "all": every cell type of the analysis)
CTs_env <- Sys.getenv("BF_CELL_TYPES", unset = "")
CTs <- if (nzchar(CTs_env) && !identical(CTs_env, "all"))
    trimws(strsplit(CTs_env, ",")[[1]]) else CTs_all
bad_ct <- setdiff(CTs, CTs_all)
if (length(bad_ct))
    stop("BF_CELL_TYPES not valid for ", ANALYSIS, ": ", paste(bad_ct, collapse = ", "),
         " (expected subset of ", paste(CTs_all, collapse = ", "), ")")

OUT_DIR   <- Sys.getenv("BF_OUT_DIR",   unset = RES_DIR)
CHUNK     <- as.integer(Sys.getenv("BF_CHUNK",     unset = "500000"))
PREFILTER <- as.numeric(Sys.getenv("BF_PREFILTER", unset = "1e-4"))
ALPHA     <- 0.05

## Phenotype BEDs: for each cell type the significant cis and trans genes are
## extracted from the whole-gene BED into <...>.sigcis.tensorqtl.bed.gz /
## <...>.sigtrans.tensorqtl.bed.gz, so the tensorQTL trans run has a small
## enough phenotype set to fit in memory. The BED name and location come from
## ANALYSIS above; LEADVAR_EXPR_DIR overrides the directory and BF_BED_TEMPLATE
## ({CT} substituted) the filename. Subset BEDs go to BF_BED_OUT_DIR, which
## defaults to the results directory, not the shared expression tree.

BED_DIR      <- Sys.getenv("BF_BED_DIR",      unset = RES_DIR)
BED_TEMPLATE <- Sys.getenv("BF_BED_TEMPLATE", unset = "")
BED_OUT_DIR  <- Sys.getenv("BF_BED_OUT_DIR",  unset = OUT_DIR)
WRITE_BEDS   <- Sys.getenv("BF_WRITE_BEDS",   unset = "1")

## Fraction of significant genes that must be present in the whole-gene BED
## before the subset is trusted. A low value almost always means the BED comes
## from a different cohort (or a differently filtered gene universe) than the
## regression results - e.g. WellcomeSanger results against an ADSER BED.
MIN_KEPT_FRAC <- as.numeric(Sys.getenv("BF_MIN_KEPT_FRAC", unset = "0.5"))

if (WRITE_BEDS == "1" && is.na(EXPR_DIR))
    message("NOTE: no expression directory found for ", ANALYSIS,
            "; set LEADVAR_EXPR_DIR if the subset BEDs are wanted.")

expr_file <- function(ct) {
    f <- if (nzchar(BED_TEMPLATE)) sub("{CT}", ct, BED_TEMPLATE, fixed = TRUE) else bed_name(ct)
    cand <- character(0)
    if (!is.na(EXPR_DIR))
        cand <- c(cand, file.path(EXPR_DIR, bed_subdirs(ct), f))
    cand <- c(cand, file.path(BED_DIR, f))
    first_existing(cand)
}

sig_bed_name <- function(bed_in, cls) {
    b <- basename(bed_in)
    if (grepl("whole_gene", b, fixed = TRUE))
        return(sub("whole_gene", paste0("sig", cls), b, fixed = TRUE))
    if (grepl("\\.tensorqtl\\.bed\\.gz$", b))
        return(sub("\\.tensorqtl\\.bed\\.gz$", paste0(".sig", cls, ".tensorqtl.bed.gz"), b))
    sub("\\.bed\\.gz$", paste0(".sig", cls, ".bed.gz"), b)
}

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
if (WRITE_BEDS == "1") dir.create(BED_OUT_DIR, recursive = TRUE, showWarnings = FALSE)

CIS_LAB   <- "leadvar_cisExpr"
TRANS_LAB <- "leadvar_transExpr"

res_file <- function(ct) file.path(RES_DIR, paste0(ct, "_leadvar_expr_lm.csv"))

message("Analysis   : ", ANALYSIS)
message("Results dir: ", RES_DIR)
message("Expression : ", ifelse(is.na(EXPR_DIR), "<not found>", EXPR_DIR))
message("Cell types requested: ", paste(CTs, collapse = ", "))


# ------------------------------------------------------------------------------
# One streaming pass over a results file: count tests per class, keep candidates
# ------------------------------------------------------------------------------

scan_file <- function(path, keep_thresh) {
    con <- file(path, "r")
    on.exit(close(con))

    hdr <- readLines(con, n = 1L)
    n_cis <- 0L; n_trans <- 0L; n_other <- 0L; n_rows <- 0L
    cand <- list(); b <- 0L

    repeat {
        lines <- readLines(con, n = CHUNK)
        if (length(lines) == 0L) break

        dt <- fread(text = c(hdr, lines), showProgress = FALSE)
        n_rows  <- n_rows  + nrow(dt)
        n_cis   <- n_cis   + sum(dt$Cis_Trans == CIS_LAB)
        n_trans <- n_trans + sum(dt$Cis_Trans == TRANS_LAB)
        n_other <- n_other + sum(!dt$Cis_Trans %in% c(CIS_LAB, TRANS_LAB))

        hit <- dt[!is.na(p) & p <= keep_thresh]
        if (nrow(hit)) { b <- b + 1L; cand[[b]] <- hit }
    }

    list(n_rows = n_rows, n_cis = n_cis, n_trans = n_trans, n_other = n_other,
         cand = if (b) rbindlist(cand) else NULL)
}


# ------------------------------------------------------------------------------
# Extract a gene subset from a tensorQTL phenotype BED.
# Lines are passed through untouched (no parsing of the sample columns) and the
# match is on the bare ENSG, so GENCODE version suffixes never cause a miss.
# ------------------------------------------------------------------------------

subset_bed <- function(bed_in, gene_ids, bed_out) {
    con_in <- gzfile(bed_in, "rt")
    hdr    <- readLines(con_in, n = 1L)
    body   <- readLines(con_in)
    close(con_in)

    gid  <- sub("^(?:[^\t]*\t){3}([^\t]*)\t.*$", "\\1", body, perl = TRUE)
    base <- unique(sub("\\..*$", "", gene_ids))
    keep <- sub("\\..*$", "", gid) %in% base

    con_out <- gzfile(bed_out, "wt")
    writeLines(c(hdr, body[keep]), con_out)
    close(con_out)

    absent <- setdiff(base, sub("\\..*$", "", gid))
    list(n_in = length(body), n_kept = sum(keep), n_requested = length(base),
         n_absent = length(absent), absent = absent)
}


# ------------------------------------------------------------------------------
# Per cell type
# ------------------------------------------------------------------------------

cis_list <- list(); trans_list <- list(); summ <- list()
missing_ct <- character(0)

for (ct in CTs) {

    f <- res_file(ct)

    if (!file.exists(f)) {
        message(ct, ": results file not found, skipped (", f, ")")
        missing_ct <- c(missing_ct, ct)
        next
    }

    message(ct, ": scanning ", round(file.size(f) / 2^30, 2), " GB ...")
    s <- scan_file(f, PREFILTER)

    cutoff_cis   <- if (s$n_cis   > 0) ALPHA / s$n_cis   else NA_real_
    cutoff_trans <- if (s$n_trans > 0) ALPHA / s$n_trans else NA_real_

    need_rescan <- (!is.na(cutoff_cis)   && cutoff_cis   > PREFILTER) ||
                   (!is.na(cutoff_trans) && cutoff_trans > PREFILTER)
    if (need_rescan) {
        loose <- max(cutoff_cis, cutoff_trans, na.rm = TRUE)
        message("  prefilter ", PREFILTER, " tighter than a cutoff; re-scanning at ", loose)
        s$cand <- scan_file(f, loose)$cand
    }

    cand <- if (is.null(s$cand)) data.table() else s$cand

    cis_sig   <- if (nrow(cand) && !is.na(cutoff_cis))   cand[Cis_Trans == CIS_LAB   & p < cutoff_cis]   else cand[0]
    trans_sig <- if (nrow(cand) && !is.na(cutoff_trans)) cand[Cis_Trans == TRANS_LAB & p < cutoff_trans] else cand[0]

    ## carry the cutoff that produced each call
    if (nrow(cis_sig))   cis_sig[,   bonferroni_cutoff := cutoff_cis]
    if (nrow(trans_sig)) trans_sig[, bonferroni_cutoff := cutoff_trans]

    ## the CellType column comes from the regression output; fall back to ct
    if (nrow(cis_sig)   && !"CellType" %in% names(cis_sig))   cis_sig[,   CellType := ct]
    if (nrow(trans_sig) && !"CellType" %in% names(trans_sig)) trans_sig[, CellType := ct]

    cis_list[[ct]]   <- cis_sig
    trans_list[[ct]] <- trans_sig

    ## significant-gene phenotype BEDs for this cell type
    if (WRITE_BEDS == "1") {
        bed_in <- expr_file(ct)

        if (is.na(bed_in)) {
            message("  whole-gene BED not found under ", EXPR_DIR, " or ", BED_DIR,
                    "; BEDs not written for ", ct)
        } else {
            for (cls in c("cis", "trans")) {
                sig <- if (cls == "cis") cis_sig else trans_sig

                out_name <- sig_bed_name(bed_in, cls)
                bed_out  <- file.path(BED_OUT_DIR, out_name)

                if (!nrow(sig)) {
                    message("  ", cls, ": no significant genes; ", out_name, " not written")
                    next
                }

                b <- subset_bed(bed_in, unique(sig$gene_id), bed_out)
                frac <- if (b$n_requested > 0) b$n_kept / b$n_requested else 0

                message("  ", cls, " BED: kept ", b$n_kept, " of ", b$n_requested,
                        " genes (", b$n_absent, " absent from the BED) -> ", out_name)

                if (b$n_absent > 0) {
                    absent_file <- file.path(BED_OUT_DIR,
                                             sub("\\.bed\\.gz$", ".absent_genes.txt", out_name))
                    writeLines(b$absent, absent_file)
                    message("    absent gene IDs -> ", basename(absent_file))
                }

                if (frac < MIN_KEPT_FRAC) {
                    message("    WARNING: only ", round(100 * frac, 1), "% of the significant ",
                            cls, " genes are in this BED.")
                    message("    WARNING: the BED is probably from a different cohort or gene ",
                            "universe than the regression results.")
                    message("    WARNING: BED = ", bed_in)
                    message("    WARNING: set LEADVAR_EXPR_DIR to the expression tree of the ",
                            "cohort that produced ", basename(RES_DIR), ".")
                }

                if (b$n_kept == 0)
                    message("    NOTE: ", out_name, " contains the header only.")
            }
        }
    }

    summ[[ct]] <- data.table(
        CellType = ct,
        n_tests_total = s$n_rows, n_tests_cis = s$n_cis, n_tests_trans = s$n_trans,
        n_tests_unclassified = s$n_other,
        bonferroni_cutoff_cis = cutoff_cis, bonferroni_cutoff_trans = cutoff_trans,
        n_significant_cis = nrow(cis_sig), n_significant_trans = nrow(trans_sig),
        n_unique_genes_cis   = if (nrow(cis_sig))   uniqueN(cis_sig$gene_id)   else 0L,
        n_unique_genes_trans = if (nrow(trans_sig)) uniqueN(trans_sig$gene_id) else 0L)

    message("  cis: ", nrow(cis_sig), " of ", s$n_cis, " tests (cutoff ",
            signif(cutoff_cis, 5), ") | trans: ", nrow(trans_sig), " of ",
            s$n_trans, " tests (cutoff ", signif(cutoff_trans, 5), ")")
}

if (length(summ) == 0)
    stop("no results files found under ", RES_DIR)


# ------------------------------------------------------------------------------
# Combine, add CellTypes, write
# ------------------------------------------------------------------------------

## CellTypes = all cell types where this lead variant x gene pair is significant
## in the same class; n_CellTypes is how many.
add_celltypes <- function(x) {
    if (!nrow(x)) return(x)
    key_rs <- if ("GWAS_RSID" %in% names(x)) "GWAS_RSID" else "geno_SNP"
    x[, CellTypes   := paste(sort(unique(CellType)), collapse = ";"),
        by = c(key_rs, "gene_id")]
    x[, n_CellTypes := uniqueN(CellType), by = c(key_rs, "gene_id")]
    setorder(x, p)
    x[]
}

cis_all   <- add_celltypes(rbindlist(cis_list,   use.names = TRUE, fill = TRUE))
trans_all <- add_celltypes(rbindlist(trans_list, use.names = TRUE, fill = TRUE))

## per-cell-type files (now carrying CellTypes) and the combined tables
for (ct in names(cis_list)) {
    fwrite(if (nrow(cis_all))   cis_all[CellType == ct]   else cis_all,
           file.path(OUT_DIR, paste0(ct, "_cis_BFsig.csv")))
    fwrite(if (nrow(trans_all)) trans_all[CellType == ct] else trans_all,
           file.path(OUT_DIR, paste0(ct, "_trans_BFsig.csv")))
}

fwrite(cis_all,   file.path(OUT_DIR, "AllCellTypes_cis_BFsig.csv"))
fwrite(trans_all, file.path(OUT_DIR, "AllCellTypes_trans_BFsig.csv"))

summary_dt <- rbindlist(summ)
summary_dt <- rbind(summary_dt, data.table(
    CellType = "TOTAL",
    n_tests_total = sum(summary_dt$n_tests_total),
    n_tests_cis   = sum(summary_dt$n_tests_cis),
    n_tests_trans = sum(summary_dt$n_tests_trans),
    n_tests_unclassified = sum(summary_dt$n_tests_unclassified),
    bonferroni_cutoff_cis = NA_real_, bonferroni_cutoff_trans = NA_real_,
    n_significant_cis   = nrow(cis_all),
    n_significant_trans = nrow(trans_all),
    n_unique_genes_cis   = if (nrow(cis_all))   uniqueN(cis_all$gene_id)   else 0L,
    n_unique_genes_trans = if (nrow(trans_all)) uniqueN(trans_all$gene_id) else 0L))

fwrite(summary_dt, file.path(OUT_DIR, "AllCellTypes_BF_summary.csv"))

message("\n", paste(capture.output(print(summary_dt)), collapse = "\n"))
if (length(missing_ct))
    message("\nNOT RUN (no results file): ", paste(missing_ct, collapse = ", "))
message("\nWrote: AllCellTypes_cis_BFsig.csv (", nrow(cis_all), " rows), ",
        "AllCellTypes_trans_BFsig.csv (", nrow(trans_all), " rows), ",
        "AllCellTypes_BF_summary.csv")
