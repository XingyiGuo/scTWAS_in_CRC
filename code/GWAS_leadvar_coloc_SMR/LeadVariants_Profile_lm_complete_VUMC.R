# conda activate tftwas_coloc
# ------------------------------------------------------------------------------
# Expression ~ GWAS lead-variant genotype, per cell type.
#
# For every CRC GWAS lead variant and every gene in the PEER-residual expression
# file: lm(expression ~ genotype), one row per lead variant x gene, annotated
# leadvar_cisExpr / leadvar_transExpr by whether the gene lies within +/-500 kb
# of the lead variant on the same chromosome.
#
# Usage:  Rscript LeadVariants_Profile_lm_complete.R <CT index|CT name> [chroms]
#         Rscript LeadVariants_Profile_lm_complete.R 1            # ABS, chr1-22
#         Rscript LeadVariants_Profile_lm_complete.R ABS 22       # ABS, chr22
#
# Path overrides (environment): LEADVAR_EXPR_DIR, LEADVAR_TRAW_DIR,
#                               LEADVAR_CSV, LEADVAR_OUT_DIR
# ------------------------------------------------------------------------------

suppressPackageStartupMessages({
    library(data.table)
    library(R.utils)          # fread() needs this to read .gz directly
})

CTs_list <- c("ABS", "CT", "EE", "GOB", "STM", "TAC", "TUF")

args <- if (exists("LEADVAR_ARGS")) LEADVAR_ARGS else commandArgs(trailingOnly = TRUE)

if (length(args) < 1)
    stop("give a cell-type index (1-7) or name, e.g. Rscript LeadVariants_Profile_lm_complete.R 1")

CT <- if (grepl("^[0-9]+$", args[1])) CTs_list[as.integer(args[1])] else args[1]

if (is.na(CT) || !CT %in% CTs_list)
    stop("unknown cell type: ", args[1], " (expected one of ", paste(CTs_list, collapse = ", "), ")")

CHROMS <- if (length(args) >= 2) as.integer(strsplit(args[2], ",")[[1]]) else 1:22

FLANK <- 500000   # cis window around the lead variant, bp


# ------------------------------------------------------------------------------
# Paths: first existing wins (cluster, then the local example folder)
# ------------------------------------------------------------------------------

first_existing <- function(...) {
    for (p in c(...)) if (file.exists(p)) return(p)
    NA_character_
}

CLUSTER_ROOT <- "/data/l2_bioinfo1/liq17/scTWAS/data/VUMC_HTAN_GSE266546/"
LOCAL_ROOT   <- "W:/VUMC_projects/scTWAS/GWAS_leadvar_coloc_SMR"

EXPR_DIR <- Sys.getenv("LEADVAR_EXPR_DIR", unset = first_existing(
    file.path(CLUSTER_ROOT, "/expression/pseudobulk_expression_for_peer_NL"),
    LOCAL_ROOT))

TRAW_DIR <- Sys.getenv("LEADVAR_TRAW_DIR", unset = first_existing(
    file.path(CLUSTER_ROOT, "genotype/VUMC_HTAN_GSE266546/R03/06_traw"),
    LOCAL_ROOT))

LEADVAR_CSV <- Sys.getenv("LEADVAR_CSV", unset = first_existing(
    paste0("/data/l2_bioinfo1/liq17/Cancer_GWAS_SS/CRC_GWAS_SS/Metal_202608/GC_Final/processed/",
           "CRC_GWAS_Clumping/EUR_plus_MULTI_loci_20261025/CRC_238signals_plus_additional_EUR_MULTI_loci.csv"),
    file.path(LOCAL_ROOT, "CRC_238signals_plus_additional_EUR_MULTI_loci.csv")))

OUT_DIR <- "/data/l2_bioinfo1/liq17/scTWAS/GWAS_leadvar_coloc_SMR_res/VUMC_HTAN_GSE266546/"

stopifnot(!is.na(EXPR_DIR), !is.na(TRAW_DIR), !is.na(LEADVAR_CSV), !is.na(OUT_DIR))

expr_file <- function(ct) {
    f <- paste0("pseudobulk_counts_",ct,".tpmlog2.qc.qn.rmlowcounts.inv.PEER_residuals.genotype_id.tensorqtl.bed.gz")  #
    first_existing(file.path(EXPR_DIR, paste0(ct, "/PEER_", ct), f), file.path(EXPR_DIR, f))
}

traw_file <- function(chrom)
    file.path(TRAW_DIR, sprintf("chr%d.HTAN_GSE266546.qc.rsid.traw", chrom)) #

OUT_CSV <- file.path(OUT_DIR, sprintf("%s_leadvar_expr_lm.csv", CT))
LOG_TXT <- file.path(OUT_DIR, sprintf("%s_leadvar_expr_lm.log", CT))

message("Cell type   : ", CT)
message("Expression  : ", expr_file(CT))
message("Genotype dir: ", TRAW_DIR)
message("Lead variants: ", LEADVAR_CSV)
message("Output      : ", OUT_CSV)


# ------------------------------------------------------------------------------
# Expression (tensorqtl BED: #chr, start, end, gene_id, then one column per sample)
# ------------------------------------------------------------------------------

EXPR_FILE <- expr_file(CT)
stopifnot(!is.na(EXPR_FILE))

expr <- fread(EXPR_FILE)
setnames(expr, 1:4, c("gene_chr", "gene_start", "gene_end", "gene_id"))
expr[, gene_chr := sub("^chr", "", gene_chr)]

gene_info <- expr[, .(gene_id, gene_chr, gene_start, gene_end)]
expr_samples <- setdiff(names(expr), c("gene_chr", "gene_start", "gene_end", "gene_id"))

message("Genes: ", nrow(gene_info), " | expression samples: ", length(expr_samples))


# ------------------------------------------------------------------------------
# Lead variants
# ------------------------------------------------------------------------------

lead <- fread(LEADVAR_CSV)
stopifnot(all(c("RSID", "CHR", "BP") %in% names(lead)))
gwas_cols <- names(lead)

message("Lead variants: ", nrow(lead), " on ", uniqueN(lead$CHR), " chromosomes")


# ------------------------------------------------------------------------------
# Vectorised simple regression: expression ~ genotype, one genotype, all genes.
# Identical to lm(y ~ g) per gene (intercept + slope), computed in closed form.
# ------------------------------------------------------------------------------

lm_all_genes <- function(Y, g) {
    ok <- !is.na(g)
    n  <- sum(ok)
    if (n < 10 || length(unique(g[ok])) < 2)
        return(NULL)

    gc  <- g[ok] - mean(g[ok])
    Sgg <- sum(gc^2)
    if (Sgg <= 0)
        return(NULL)

    Yo  <- if (all(ok)) Y else Y[, ok, drop = FALSE]
    Yc  <- Yo - rowMeans(Yo)

    num  <- as.vector(Yc %*% gc)
    beta <- num / Sgg
    rss  <- pmax(rowSums(Yc^2) - num^2 / Sgg, 0)
    dfr  <- n - 2
    se   <- sqrt((rss / dfr) / Sgg)
    tval <- beta / se

    list(n = n, beta = beta, se = se, t = tval,
         p = 2 * pt(-abs(tval), dfr))
}


# ------------------------------------------------------------------------------
# Per chromosome: pull each lead variant's genotype row out of the .traw and
# regress every gene on it. Results are appended to one CSV as they are produced.
# ------------------------------------------------------------------------------

if (file.exists(OUT_CSV)) file.remove(OUT_CSV)

n_written  <- 0L
n_tested   <- 0L
n_missing  <- character(0)
first_write <- TRUE

for (chrom in CHROMS) {

    lead_chr <- lead[CHR == chrom]
    if (nrow(lead_chr) == 0) next

    TRAW <- traw_file(chrom)
    if (!file.exists(TRAW)) {
        warning("missing traw for chr", chrom, ": ", TRAW)
        next
    }

    # Header once, then only the key columns, so the 278-sample matrix is read
    # for the handful of lead-variant rows instead of the whole chromosome.
    traw_hdr <- names(fread(TRAW, nrows = 0))
    traw_key <- fread(TRAW, select = c("CHR", "SNP", "POS"))

    traw_samples <- sub("^[^_]*_", "", setdiff(traw_hdr, c("CHR", "SNP", "(C)M", "POS", "COUNTED", "ALT")))
    geno_cols    <- setdiff(traw_hdr, c("CHR", "SNP", "(C)M", "POS", "COUNTED", "ALT"))

    # Samples shared by expression and genotype, in expression order
    use_samples <- intersect(expr_samples, traw_samples)
    if (length(use_samples) < 10)
        stop("only ", length(use_samples), " samples shared between expression and chr", chrom, " genotype")

    Y <- as.matrix(expr[, ..use_samples])
    geno_col_for_sample <- geno_cols[match(use_samples, traw_samples)]

    message("chr", chrom, ": ", nrow(lead_chr), " lead variant(s), ",
            length(use_samples), " shared samples")

    for (k in seq_len(nrow(lead_chr))) {

        rs <- lead_chr$RSID[k]
        bp <- lead_chr$BP[k]

        # rsID first, then chr:pos for variants the genotype calls differently
        hit <- which(traw_key$SNP == rs)
        matched_by <- "rsID"
        if (length(hit) == 0) {
            hit <- which(traw_key$POS == bp)
            matched_by <- "chr:pos"
        }
        if (length(hit) == 0) {
            n_missing <- c(n_missing, rs)
            message("  ", rs, " (chr", chrom, ":", bp, ") not in genotype after QC; skipped")
            next
        }
        hit <- hit[1]

        # Read just that genotype row (skip = data row index, header is line 1)
        gt_row <- fread(TRAW, skip = hit, nrows = 1, header = FALSE, col.names = traw_hdr)
        g <- as.numeric(unlist(gt_row[1, ..geno_col_for_sample]))

        fit <- lm_all_genes(Y, g)
        if (is.null(fit)) {
            message("  ", rs, ": genotype not usable (monomorphic or too few calls); skipped")
            next
        }

        cis <- gene_info$gene_chr == as.character(chrom) &
               gene_info$gene_end   >= (bp - FLANK) &
               gene_info$gene_start <= (bp + FLANK)

        out <- data.table(
            CellType   = CT,
            lead_chr   = chrom,
            matched_by = matched_by,
            geno_SNP   = gt_row$SNP[1],
            geno_POS   = gt_row$POS[1],
            COUNTED    = gt_row$COUNTED[1],
            ALT        = gt_row$ALT[1],
            gene_id    = gene_info$gene_id,
            gene_chr   = gene_info$gene_chr,
            gene_start = gene_info$gene_start,
            gene_end   = gene_info$gene_end,
            Cis_Trans  = ifelse(cis, "leadvar_cisExpr", "leadvar_transExpr"),
            n          = fit$n,
            beta       = fit$beta,
            se         = fit$se,
            t          = fit$t,
            p          = fit$p
        )

        # every column of the GWAS table, carried through unchanged
        for (cc in gwas_cols)
            set(out, j = paste0("GWAS_", cc), value = lead_chr[[cc]][k])

        setcolorder(out, c("CellType", paste0("GWAS_", gwas_cols), "matched_by",
                           "geno_SNP", "geno_POS", "COUNTED", "ALT",
                           "gene_id", "gene_chr", "gene_start", "gene_end",
                           "Cis_Trans", "n", "beta", "se", "t", "p"))
        out[, lead_chr := NULL]

        fwrite(out, OUT_CSV, append = !first_write, col.names = first_write)
        first_write <- FALSE

        n_written <- n_written + nrow(out)
        n_tested  <- n_tested + 1L

        message("  ", rs, " (", matched_by, "): ", nrow(out), " genes, ",
                sum(cis), " cis, min p = ", signif(min(fit$p), 3))
    }

    rm(Y, traw_key); gc(verbose = FALSE)
}


# ------------------------------------------------------------------------------
# Run log
# ------------------------------------------------------------------------------

writeLines(c(
    paste0("cell_type=", CT),
    paste0("expression_file=", EXPR_FILE),
    paste0("genotype_dir=", TRAW_DIR),
    paste0("leadvar_table=", LEADVAR_CSV),
    paste0("chromosomes=", paste(CHROMS, collapse = ",")),
    paste0("flank_bp=", FLANK),
    paste0("n_genes=", nrow(gene_info)),
    paste0("n_lead_variants_in_table=", nrow(lead)),
    paste0("n_lead_variants_tested=", n_tested),
    paste0("n_lead_variants_absent_from_genotype=", length(n_missing)),
    paste0("lead_variants_absent=", paste(n_missing, collapse = ",")),
    paste0("n_rows_written=", n_written),
    paste0("output_csv=", OUT_CSV)
), LOG_TXT)

message("Done: ", n_tested, " lead variant(s) x ", nrow(gene_info), " genes = ",
        n_written, " rows -> ", OUT_CSV)
