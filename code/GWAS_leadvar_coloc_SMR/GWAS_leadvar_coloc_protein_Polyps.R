# Optimized trans reader: vectorized variant parsing + disk cache.
# Run: Rscript GWAS_leadvar_coloc_protein_fast.R transprotein
# Cache: <COLOC_RES_DIR>/trans_coloc_cache (delete or set COLOC_TRANS_CACHE_DIR="" to disable).
# The first uncached run still scans each large trans parquet once.
# Statistical coloc settings and input/output formats are unchanged.

# conda activate tftwas_coloc
# ==============================================================================
# GWAS_leadvar_coloc_protein.R
#
# Colocalisation between the CRC GWAS and cell-type eQTLs, one test per row of
# the lead-variant significance table. One script, two modes -- the mode picks
# the significance file, the QTL file layout and the QTL reader:
#
#   cisprotein    AllCellTypes_cis_BFsig.csv    -> AllCellTypes_cis_BFsig.coloc.csv
#                 <PREFIX>.cis_qtl_pairs.chr<N>.parquet   (one per chromosome)
#                 columns slope / slope_se / pval_nominal
#
#   transprotein  AllCellTypes_trans_BFsig.csv  -> AllCellTypes_trans_BFsig.coloc.csv
#                 <PREFIX>.GWAS_leadvar_trans.<CT>.trans_qtl_pairs.parquet
#                                                         (one per cell type)
#                 columns b / b_se / pval
#
# Everything else -- GWAS reader, harmonisation, coloc call, output assembly --
# is shared, so cis and trans results are directly comparable.
#
# Conventions follow 01_coloc_gwas_xqtl_sigSource_peakSNP.R:
#   GWAS      : RSID/A1/A2/FRQ/BETA/SE/P, harmonised to the QTL ALT allele
#   QTL       : tensorQTL parquet, variant_id = rsID;chrom_pos_REF_ALT,
#               effect size ALT-oriented
#   coloc     : coloc::coloc.abf(dataset1 = cc GWAS, dataset2 = quant QTL, MAF)
#               with package default priors (p1 = p2 = 1e-4, p12 = 1e-5)
#   filtering : biallelic ACGT, SE > 0, 0 < P <= 1, strand-ambiguous dropped,
#               >= min_overlap_snps harmonised SNPs required
#
# The coloc region is centred on the GWAS LEAD VARIANT (GWAS_BP +/- flank), not
# on the gene TSS, because the task list is lead-variant driven; in trans mode
# the gene usually sits on another chromosome, so the lead variant is the only
# sensible centre. Both sides are restricted to that same window.
#
# TWO THINGS TO KNOW ABOUT THE QTL INPUT
#   1. The input is a tensorQTL parquet, NOT <CT>...sigcis/sigtrans.tensorqtl.bed.gz
#      -- those BEDs are the PHENOTYPE matrices fed to tensorQTL (expression
#      residuals); they carry no betas or standard errors, so coloc cannot be
#      run from them.
#   2. tensorQTL --mode trans always drops pairs within +/-5 Mb of the gene
#      (trans.filter_cis, hardcoded) and writes only pairs below
#      --pval_threshold. The trans run must therefore use
#      TRANS_PVAL_THRESHOLD=1 (the default in xQTL_per_chrom_Wellcome.sh);
#      sparse output leaves too few SNPs per window. Trans rows whose whole
#      window falls inside that +/-5 Mb exclusion are reported as
#      inside_tensorqtl_cis_window_5Mb -- those pairs belong to cis mode.
#
# GWAS cohort is routed by GWAS_Source:
#   PMID38670944_238signal -> Multi,  MULTI_additional -> Multi,
#   EUR_additional         -> EUR
#
# Usage:
#   Rscript GWAS_leadvar_coloc_protein.R cisprotein              # all cell types
#   Rscript GWAS_leadvar_coloc_protein.R transprotein            # all cell types
#   Rscript GWAS_leadvar_coloc_protein.R cisprotein   ABS        # one cell type
#   Rscript GWAS_leadvar_coloc_protein.R transprotein ABS,GOB
#
# Environment overrides (all optional):
#   COLOC_MODE, COLOC_RES_DIR, COLOC_GWAS_ROOT,
#   COLOC_QTL_PREFIX, COLOC_QTL_TEMPLATE   (apply to the selected mode)
# ==============================================================================

suppressPackageStartupMessages({
    library(data.table)
    library(arrow)
    library(coloc)
})

args <- if (exists("COLOC_ARGS")) COLOC_ARGS else commandArgs(trailingOnly = TRUE)

MODE_NAME <- if (length(args) >= 1 && nzchar(args[1]))
    tolower(trimws(args[1])) else tolower(Sys.getenv("COLOC_MODE", unset = ""))

if (!MODE_NAME %in% c("cisprotein", "transprotein"))
    stop("First argument must be 'cisprotein' or 'transprotein' (got: '",
         MODE_NAME, "').\n",
         "  Rscript GWAS_leadvar_coloc_protein.R cisprotein [ABS,GOB]\n",
         "  Rscript GWAS_leadvar_coloc_protein.R transprotein [ABS,GOB]")


# ------------------------------------------------------------------------------
# Configuration
# ------------------------------------------------------------------------------

RES_DIR <- Sys.getenv("COLOC_RES_DIR", unset =
    "/data/l2_bioinfo1/liq17/scTWAS/GWAS_leadvar_coloc_SMR_res/Polyps")

GWAS_ROOT <- Sys.getenv("COLOC_GWAS_ROOT", unset =
    "/data/l2_bioinfo1/liq17/Cancer_GWAS_SS/CRC_GWAS_SS/Metal_202608/GC_Final/processed")

## CRC GWAS sample sizes (cases / controls / total), as supplied:
##   Multi : 132,261 / 1,752,398 / 1,884,659
##   EUR   :  99,260 / 1,119,512 / 1,218,772
## coloc needs N (total) and s (case fraction) for a case/control trait; both
## are derived here from the counts so the three numbers cannot drift apart.
gwas_cc <- function(n_cases, n_controls, n_total, file_template) {
    stopifnot(n_cases + n_controls == n_total)
    list(file_template = file_template,
         n_cases = n_cases, n_controls = n_controls, N = n_total,
         case_fraction = n_cases / n_total)
}

CFG <- list(

    cell_types = c("ABS", "CT", "EE", "GOB", "STM", "TAC", "TUF", "ASC", "SSC"),

    # GWAS_Source -> GWAS cohort
    source_to_cohort = c(
        PMID38670944_238signal = "Multi",
        MULTI_additional       = "Multi",
        EUR_additional         = "EUR"
    ),

    gwas = list(
        EUR = gwas_cc(99260, 1119512, 1218772, file.path(GWAS_ROOT,
            "EUR/rsid_annotated_CHRPOS",
            "final_gc_maf0.01_CRC_EUR_META_1.chr{chrom}.chrpos.rsid.assoc.dosage.gz")),
        Multi = gwas_cc(132261, 1752398, 1884659, file.path(GWAS_ROOT,
            "Multi/rsid_annotated_CHRPOS",
            "final_gc_maf0.01_CRC_Multi_META_1.chr{chrom}.chrpos.rsid.assoc.dosage.gz"))
    ),

    qtl_N                 = 278,          # WellcomeSanger Healthy_TI_Rectum donors
    qtl_batch_size        = as.integer(Sys.getenv("COLOC_TRANS_BATCH_SIZE", "250000")),
                                  # rows per Arrow batch in trans mode
    flank_length          = 500000,
    min_overlap_snps      = 50,
    drop_strand_ambiguous = TRUE,
    ## hardcoded in tensorqtl.trans.filter_cis(); pairs closer than this are
    ## absent from the trans parquet by construction.
    tensorqtl_cis_window  = 5000000
)

## Mode table. Only the selected mode reads COLOC_QTL_PREFIX / COLOC_QTL_TEMPLATE,
## so one override pair serves both modes.
make_mode <- function(mode_name) {
    if (identical(mode_name, "cisprotein")) {
        list(
            name      = "cisprotein",
            label     = "cis",
            sig_file  = file.path(RES_DIR, "AllCellTypes_cis_BFsig.csv"),
            out_file  = file.path(RES_DIR, "AllCellTypes_cis_BFsig.coloc.csv"),
            prefix_template = Sys.getenv("COLOC_QTL_PREFIX",
                unset = "VUMC_HTAN_ADSER_{CT}"),
            file_template = Sys.getenv("COLOC_QTL_TEMPLATE", unset = paste0(
                "/data/l2_bioinfo1/liq17/scTWAS/data/VUMC_HTAN/xQTL_VUMC_HTAN/",
                "{PREFIX}/{PREFIX}.cis_qtl_pairs.chr{chrom}.parquet")),
            ## tensorQTL cis spelling
            col_beta = "slope", col_se = "slope_se", col_p = "pval_nominal",
            col_optional = c("tss_distance", "af", "ma_samples", "ma_count"),
            one_file_per = "chromosome",
            apply_cis_window_check = FALSE,
            reader = function(CT, chrom, want) read_qtl_cis(CT, chrom)
        )
    } else {
		print("No sig trans protein, skip analyses")
		q()
    }
}

MODE <- make_mode(MODE_NAME)

CTs <- if (length(args) >= 2 && nzchar(args[2]) && !identical(args[2], "all"))
    trimws(strsplit(args[2], ",")[[1]]) else CFG$cell_types

stopifnot(file.exists(MODE$sig_file))


# ------------------------------------------------------------------------------
# Helpers (same semantics as the reference pipeline)
# ------------------------------------------------------------------------------

fill_template <- function(template, values) {
    out <- template
    for (nm in names(values))
        out <- gsub(paste0("{", nm, "}"), as.character(values[[nm]]), out, fixed = TRUE)
    out
}

safe_numeric <- function(x) suppressWarnings(as.numeric(x))

normalize_chr <- function(x) sub("^chr", "", trimws(as.character(x)), ignore.case = TRUE)

canonical_gene_id <- function(x) sub("\\..*$", "", trimws(as.character(x)))

complement_allele <- function(x) chartr("ACGT", "TGCA", toupper(x))

is_ambiguous_pair <- function(a1, a2) {
    (a1 == "A" & a2 == "T") | (a1 == "T" & a2 == "A") |
    (a1 == "C" & a2 == "G") | (a1 == "G" & a2 == "C")
}

## tensorQTL variant_id "rs123;22_15270698_T_C" -> rsID + chrom/pos/REF/ALT
parse_tensorqtl_variant <- function(variant_id) {
    x    <- as.character(variant_id)
    has  <- grepl(";", x, fixed = TRUE)
    rsid <- ifelse(has, sub(";.*$", "", x), NA_character_)
    key  <- ifelse(has, sub("^[^;]+;", "", x), x)
    p    <- strsplit(key, "_", fixed = TRUE)
    pick <- function(i) vapply(p, function(z) if (length(z) >= 4L) z[i] else NA_character_,
                               character(1))
    data.frame(
        SNP         = tolower(trimws(rsid)),
        variant_chr = normalize_chr(pick(1)),
        variant_pos = safe_numeric(pick(2)),
        REF         = toupper(trimws(pick(3))),
        ALT         = toupper(trimws(pick(4))),
        stringsAsFactors = FALSE
    )
}

qtl_path <- function(CT, chrom) {
    prefix <- fill_template(MODE$prefix_template, list(CT = CT))
    fill_template(MODE$file_template, list(PREFIX = prefix, CT = CT, chrom = chrom))
}

## shared QC on a QTL table that already carries BETA/SE/P/SNP/REF/ALT
clean_qtl <- function(eq) {
    keep <- grepl("^rs[0-9]+$", eq$SNP, ignore.case = TRUE) &
        eq$REF %in% c("A", "C", "G", "T") &
        eq$ALT %in% c("A", "C", "G", "T") &
        eq$REF != eq$ALT &
        is.finite(eq$BETA) &
        is.finite(eq$SE) & eq$SE > 0 &
        is.finite(eq$P) & eq$P > 0 & eq$P <= 1

    eq <- eq[keep, , drop = FALSE]
    eq <- eq[order(eq$P), , drop = FALSE]
    eq[!duplicated(paste(eq$phenotype_key, eq$SNP, sep = "|")), , drop = FALSE]
}


# ------------------------------------------------------------------------------
# Readers. One GWAS chromosome is held at a time; the QTL cache holds one
# chromosome (cis) or one cell type (trans), because the task list is processed
# in (cell type, chromosome) order.
# ------------------------------------------------------------------------------

GWAS_CACHE <- new.env(parent = emptyenv())
QTL_CACHE  <- new.env(parent = emptyenv())

read_gwas_chr <- function(cohort, chrom) {
    key <- paste(cohort, chrom, sep = "|")
    if (exists(key, envir = GWAS_CACHE, inherits = FALSE))
        return(get(key, envir = GWAS_CACHE, inherits = FALSE))
    rm(list = ls(envir = GWAS_CACHE), envir = GWAS_CACHE)   # single slot

    path <- fill_template(CFG$gwas[[cohort]]$file_template, list(chrom = chrom))
    if (!file.exists(path)) stop("GWAS file not found: ", path)

    g <- fread(path, data.table = FALSE, check.names = FALSE)

    req <- c("RSID", "A1", "A2", "BETA", "SE", "P")
    miss <- setdiff(req, names(g))
    if (length(miss))
        stop("GWAS file missing column(s): ", paste(miss, collapse = ", "), "\nFile: ", path)

    ## SNP_CHR_POS is "chr:pos"; the position drives the lead-variant window.
    pos <- if ("SNP_CHR_POS" %in% names(g))
        safe_numeric(sub("^[^:]+:", "", as.character(g$SNP_CHR_POS))) else NA_real_

    out <- data.frame(
        SNP  = trimws(as.character(g$RSID)),
        POS  = pos,
        A1   = toupper(trimws(as.character(g$A1))),
        A2   = toupper(trimws(as.character(g$A2))),
        BETA = safe_numeric(g$BETA),
        SE   = safe_numeric(g$SE),
        P    = safe_numeric(g$P),
        ## FRQ carries the literal "not_in_1KG" for variants absent from 1000G,
        ## which becomes NA here; the QTL af then supplies the MAF.
        FRQ  = if ("FRQ" %in% names(g)) safe_numeric(g$FRQ) else NA_real_,
        stringsAsFactors = FALSE
    )

    keep <- grepl("^rs[0-9]+$", out$SNP, ignore.case = TRUE) &
        out$A1 %in% c("A", "C", "G", "T") &
        out$A2 %in% c("A", "C", "G", "T") &
        out$A1 != out$A2 &
        is.finite(out$BETA) &
        is.finite(out$SE) & out$SE > 0 &
        is.finite(out$P) & out$P > 0 & out$P <= 1

    out <- out[keep, , drop = FALSE]
    out <- out[order(out$P), , drop = FALSE]
    out <- out[!duplicated(tolower(out$SNP)), , drop = FALSE]
    out$SNP <- tolower(out$SNP)

    assign(key, out, envir = GWAS_CACHE)
    message("  ", cohort, " chr", chrom, " GWAS variants: ", nrow(out))
    out
}

## MODE = cisprotein. One parquet per chromosome, small enough to read whole.
read_qtl_cis <- function(CT, chrom) {
    key <- paste(CT, chrom, sep = "|")
    if (exists(key, envir = QTL_CACHE, inherits = FALSE))
        return(get(key, envir = QTL_CACHE, inherits = FALSE))
    rm(list = ls(envir = QTL_CACHE), envir = QTL_CACHE)     # single slot

    path <- qtl_path(CT, chrom)
    if (!file.exists(path)) stop("cis-QTL parquet not found: ", path)

    schema_names <- as.character(arrow::open_dataset(path, format = "parquet")$schema$names)
    required <- c("phenotype_id", "variant_id", MODE$col_beta, MODE$col_se, MODE$col_p)
    miss <- setdiff(required, schema_names)
    if (length(miss))
        stop("cis-QTL parquet missing column(s): ", paste(miss, collapse = ", "),
             "\nFile: ", path, "\nAvailable: ", paste(schema_names, collapse = ", "))

    optional <- intersect(MODE$col_optional, schema_names)
    ## all_of() keeps arrow's tidyselect-based col_select from warning about an
    ## external character vector (tidyselect is a hard dependency of arrow).
    eq <- as.data.frame(arrow::read_parquet(
        path, col_select = tidyselect::all_of(c(required, optional))))

    eq$phenotype_key <- canonical_gene_id(eq$phenotype_id)
    eq <- cbind(eq, parse_tensorqtl_variant(eq$variant_id))

    eq$BETA <- safe_numeric(eq[[MODE$col_beta]])
    eq$SE   <- safe_numeric(eq[[MODE$col_se]])
    eq$P    <- safe_numeric(eq[[MODE$col_p]])
    eq$AF   <- if ("af" %in% names(eq)) safe_numeric(eq$af) else NA_real_

    eq <- clean_qtl(eq)

    assign(key, eq, envir = QTL_CACHE)
    message("  ", CT, " chr", chrom, " cis-QTL rows: ", nrow(eq))
    eq
}

## MODE = transprotein. Streaming Arrow read with fast vectorized ID parsing
## and a persistent cache containing only candidate gene x GWAS-window pairs.
## The very first scan of a large non-partitioned parquet can still take time.
## Subsequent runs with the same source file and candidate windows reuse the cache.
parse_variant_fast <- function(variant_id) {
    v <- as.character(variant_id)
    n <- length(v)
    if (!n) return(data.frame(SNP = character(0), variant_chr = character(0),
                              variant_pos = numeric(0), REF = character(0),
                              ALT = character(0), stringsAsFactors = FALSE))
    ## tstrsplit uses vectorized splitting, not one R function call per SNP.
    semi <- data.table::tstrsplit(v, ";", fixed = TRUE, keep = 1:2)
    if (length(semi) < 2L)
        stop("Unexpected tensorQTL variant_id format: expected rsID;chr_pos_REF_ALT")
    part <- data.table::tstrsplit(semi[[2]], "_", fixed = TRUE, keep = 1:4)
    if (length(part) < 4L)
        stop("Unexpected tensorQTL variant_id format: expected rsID;chr_pos_REF_ALT")
    data.frame(SNP = tolower(trimws(semi[[1]])),
               variant_chr = normalize_chr(part[[1]]),
               variant_pos = safe_numeric(part[[2]]),
               REF = toupper(trimws(part[[3]])),
               ALT = toupper(trimws(part[[4]])),
               stringsAsFactors = FALSE)
}

read_qtl_trans <- function(CT, want) {
    if (exists(CT, envir = QTL_CACHE, inherits = FALSE))
        return(get(CT, envir = QTL_CACHE, inherits = FALSE))
    rm(list = ls(envir = QTL_CACHE), envir = QTL_CACHE)

    path <- qtl_path(CT, NA)
    if (!file.exists(path)) stop("trans-QTL parquet not found: ", path)

    ## Invalidate a persisted cache whenever source file or candidate windows change.
    info <- file.info(path)
    source_sig <- c(normalizePath(path, mustWork = TRUE),
                    as.character(info$size), as.character(as.numeric(info$mtime)))
    want_sig <- paste(sort(paste(want$key, want$lo, want$hi, sep = ":")),
                      collapse = "\n")
    cache_dir <- Sys.getenv("COLOC_TRANS_CACHE_DIR",
                            unset = file.path(RES_DIR, "trans_coloc_cache"))
    use_cache <- nzchar(cache_dir)
    if (use_cache) {
        dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
        cache_file <- file.path(cache_dir, paste0("trans_coloc_", CT, ".rds"))
        if (file.exists(cache_file)) {
            cached <- tryCatch(readRDS(cache_file), error = function(e) NULL)
            if (is.list(cached) && identical(cached$source_sig, source_sig) &&
                identical(cached$want_sig, want_sig)) {
                message("  ", CT, " loaded filtered trans-QTL cache: ",
                        nrow(cached$eq), " rows")
                assign(CT, cached$eq, envir = QTL_CACHE)
                return(cached$eq)
            }
        }
    }

    t_start <- proc.time()[[3]]
    message("  ", CT, " first pass through trans-QTL parquet: ", path)
    ds <- arrow::open_dataset(path, format = "parquet")
    schema_names <- as.character(ds$schema$names)
    required <- c("phenotype_id", "variant_id", MODE$col_beta, MODE$col_se, MODE$col_p)
    miss <- setdiff(required, schema_names)
    if (length(miss)) stop("trans-QTL parquet missing: ", paste(miss, collapse = ", "))
    optional <- intersect(MODE$col_optional, schema_names)

    genes <- unique(c(as.character(want$gene_id), canonical_gene_id(want$gene_id)))
    expr <- arrow::Expression$create(
        "is_in", arrow::Expression$field_ref("phenotype_id"),
        options = list(value_set = arrow::Array$create(genes), skip_nulls = TRUE))
    sc <- arrow::Scanner$create(ds, filter = expr,
                                projection = c(required, optional),
                                batch_size = CFG$qtl_batch_size)
    rbr <- sc$ToRecordBatchReader()

    acc <- list(); n_seen <- 0; n_kept <- 0; n_batch <- 0L
    repeat {
        batch <- rbr$read_next_batch()
        if (is.null(batch)) break
        n_batch <- n_batch + 1L
        eq <- as.data.frame(batch)
        n_seen <- n_seen + nrow(eq)
        if (!nrow(eq)) next

        parsed <- parse_variant_fast(eq$variant_id)
        phenotype_key <- canonical_gene_id(eq$phenotype_id)
        m <- match(paste(phenotype_key, parsed$variant_chr, sep = "|"), want$key)
        ok <- !is.na(m) & is.finite(parsed$variant_pos) &
              parsed$variant_pos >= want$lo[m] & parsed$variant_pos <= want$hi[m]
        if (any(ok)) {
            kept <- cbind(eq[ok, , drop = FALSE],
                          phenotype_key = phenotype_key[ok],
                          parsed[ok, , drop = FALSE])
            acc[[length(acc) + 1L]] <- kept
            n_kept <- n_kept + nrow(kept)
        }
        if (n_batch == 1L || n_batch %% 10L == 0L)
            message("  ", CT, " scan batch ", n_batch,
                    ": filtered-gene rows ", n_seen, ", regional rows ", n_kept,
                    ", elapsed ", round(proc.time()[[3]] - t_start), "s")
    }

    eq <- if (length(acc)) do.call(rbind, acc) else
        data.frame(phenotype_id = character(0), phenotype_key = character(0),
                   variant_id = character(0), SNP = character(0),
                   variant_chr = character(0), variant_pos = numeric(0),
                   REF = character(0), ALT = character(0),
                   stringsAsFactors = FALSE)
    eq$BETA <- if (nrow(eq)) safe_numeric(eq[[MODE$col_beta]]) else numeric(0)
    eq$SE   <- if (nrow(eq)) safe_numeric(eq[[MODE$col_se]]) else numeric(0)
    eq$P    <- if (nrow(eq)) safe_numeric(eq[[MODE$col_p]]) else numeric(0)
    eq$AF   <- if (nrow(eq) && "af" %in% names(eq)) safe_numeric(eq$af) else
        rep(NA_real_, nrow(eq))
    eq <- clean_qtl(eq)

    message("  ", CT, " scan complete: ", n_seen, " rows after gene filter; ",
            nrow(eq), " after windows/QC; elapsed ",
            round(proc.time()[[3]] - t_start), "s")
    if (use_cache) {
        saveRDS(list(source_sig = source_sig, want_sig = want_sig, eq = eq),
                cache_file, compress = FALSE)
        message("  ", CT, " saved filtered trans-QTL cache: ", cache_file)
    }
    assign(CT, eq, envir = QTL_CACHE)
    eq
}


# ------------------------------------------------------------------------------
# Harmonisation: merge on rsID, align the QTL beta to the GWAS A1
# ------------------------------------------------------------------------------

harmonize_gwas_qtl <- function(gwas, qtl) {
    x <- merge(gwas, qtl, by = "SNP", suffixes = c("_g", "_x"))
    if (!nrow(x)) return(x)

    a1g <- toupper(x$A1); a2g <- toupper(x$A2)
    a1x <- toupper(x$ALT); a2x <- toupper(x$REF)

    ambiguous <- is_ambiguous_pair(a1g, a2g)
    same      <- a1x == a1g & a2x == a2g
    flip      <- a1x == a2g & a2x == a1g
    comp_same <- a1x == complement_allele(a1g) & a2x == complement_allele(a2g)
    comp_flip <- a1x == complement_allele(a2g) & a2x == complement_allele(a1g)

    keep <- same | flip | comp_same | comp_flip
    if (isTRUE(CFG$drop_strand_ambiguous)) keep <- keep & !ambiguous

    x <- x[keep, , drop = FALSE]
    if (!nrow(x)) return(x)

    a1g <- toupper(x$A1); a2g <- toupper(x$A2)
    a1x <- toupper(x$ALT); a2x <- toupper(x$REF)
    flip      <- a1x == a2g & a2x == a1g
    comp_flip <- a1x == complement_allele(a2g) & a2x == complement_allele(a1g)

    x$BETA_x_aligned <- x$BETA_x
    x$BETA_x_aligned[flip | comp_flip] <- -x$BETA_x_aligned[flip | comp_flip]

    maf_x <- pmin(x$AF, 1 - x$AF)
    maf_g <- pmin(x$FRQ, 1 - x$FRQ)
    x$maf <- ifelse(is.finite(maf_x) & maf_x > 0 & maf_x < 0.5, maf_x, maf_g)

    x <- x[is.finite(x$maf) & x$maf > 0 & x$maf <= 0.5 &
           is.finite(x$BETA_g) & is.finite(x$SE_g) & is.finite(x$P_g) &
           is.finite(x$BETA_x_aligned) & is.finite(x$SE_x) & is.finite(x$P_x),
           , drop = FALSE]
    if (!nrow(x)) return(x)

    x <- x[order(x$P_g, x$P_x), , drop = FALSE]
    x[!duplicated(x$SNP), , drop = FALSE]
}


# ------------------------------------------------------------------------------
# One coloc test for one row of the significance table
# ------------------------------------------------------------------------------

NA_RESULT <- function(status, cohort = NA_character_, n_gwas = NA_integer_,
                      n_qtl = NA_integer_, n_overlap = NA_integer_,
                      gwas_cfg = NULL) {
    data.frame(
        coloc_nsnps     = NA_real_,
        coloc_PP.H0.abf = NA_real_, coloc_PP.H1.abf = NA_real_,
        coloc_PP.H2.abf = NA_real_, coloc_PP.H3.abf = NA_real_,
        coloc_PP.H4.abf = NA_real_,
        coloc_status = status, coloc_cohort = cohort,
        coloc_n_gwas = n_gwas, coloc_n_qtl = n_qtl, coloc_n_overlap = n_overlap,
        coloc_gwas_ncase    = if (is.null(gwas_cfg)) NA_integer_ else gwas_cfg$n_cases,
        coloc_gwas_ncontrol = if (is.null(gwas_cfg)) NA_integer_ else gwas_cfg$n_controls,
        coloc_gwas_N        = if (is.null(gwas_cfg)) NA_integer_ else gwas_cfg$N,
        coloc_gwas_s        = if (is.null(gwas_cfg)) NA_real_    else gwas_cfg$case_fraction,
        stringsAsFactors = FALSE, check.names = FALSE
    )
}

run_one_coloc <- function(row, want, read_qtl = MODE$reader) {

    cohort <- unname(CFG$source_to_cohort[as.character(row$GWAS_Source)])
    if (is.na(cohort)) return(NA_RESULT(paste0("unmapped_GWAS_Source:", row$GWAS_Source)))
    gwas_cfg <- CFG$gwas[[cohort]]

    chrom    <- normalize_chr(row$GWAS_CHR)
    lead_bp  <- safe_numeric(row$GWAS_BP)
    gene_key <- canonical_gene_id(row$gene_id)
    lo <- lead_bp - CFG$flank_length
    hi <- lead_bp + CFG$flank_length

    ## trans mode only: tensorQTL has already removed every pair within +/-5 Mb
    ## of the gene, so a window entirely inside that exclusion can never be
    ## colocalised here -- those pairs belong to cis mode.
    if (isTRUE(MODE$apply_cis_window_check) &&
        identical(normalize_chr(row$gene_chr), chrom)) {
        gl <- safe_numeric(row$gene_start) - CFG$tensorqtl_cis_window
        gh <- safe_numeric(row$gene_end)   + CFG$tensorqtl_cis_window
        if (is.finite(gl) && is.finite(gh) && lo >= gl && hi <= gh)
            return(NA_RESULT("inside_tensorqtl_cis_window_5Mb", cohort,
                             gwas_cfg = gwas_cfg))
    }

    gw <- read_gwas_chr(cohort, chrom)
    gw <- gw[is.finite(gw$POS) & gw$POS >= lo & gw$POS <= hi, , drop = FALSE]

    eq <- read_qtl(as.character(row$CellType), chrom, want)
    eq <- eq[eq$phenotype_key == gene_key, , drop = FALSE]
    if (!nrow(eq))
        return(NA_RESULT("gene_not_found_in_qtl", cohort, nrow(gw), 0L, 0L, gwas_cfg))

    eq <- eq[eq$variant_chr == chrom &
             is.finite(eq$variant_pos) & eq$variant_pos >= lo & eq$variant_pos <= hi,
             , drop = FALSE]
    if (!nrow(eq))
        return(NA_RESULT("no_qtl_variant_in_window", cohort, nrow(gw), 0L, 0L, gwas_cfg))

    matched <- harmonize_gwas_qtl(gw, eq)
    n_overlap <- nrow(matched)

    if (n_overlap < CFG$min_overlap_snps)
        return(NA_RESULT("too_few_overlap_snps", cohort, nrow(gw), nrow(eq),
                         n_overlap, gwas_cfg))

    ds1 <- list(beta = matched$BETA_g, varbeta = matched$SE_g^2,
                pvalues = matched$P_g, snp = matched$SNP,
                type = "cc", s = gwas_cfg$case_fraction,
                N = as.integer(round(gwas_cfg$N)))

    ds2 <- list(beta = matched$BETA_x_aligned, varbeta = matched$SE_x^2,
                pvalues = matched$P_x, snp = matched$SNP,
                type = "quant", N = CFG$qtl_N)

    fit <- coloc::coloc.abf(dataset1 = ds1, dataset2 = ds2, MAF = matched$maf)
    s <- fit$summary

    data.frame(
        coloc_nsnps     = unname(s["nsnps"]),
        coloc_PP.H0.abf = unname(s["PP.H0.abf"]),
        coloc_PP.H1.abf = unname(s["PP.H1.abf"]),
        coloc_PP.H2.abf = unname(s["PP.H2.abf"]),
        coloc_PP.H3.abf = unname(s["PP.H3.abf"]),
        coloc_PP.H4.abf = unname(s["PP.H4.abf"]),
        coloc_status = "coloc_completed", coloc_cohort = cohort,
        coloc_n_gwas = nrow(gw), coloc_n_qtl = nrow(eq), coloc_n_overlap = n_overlap,
        coloc_gwas_ncase    = gwas_cfg$n_cases,
        coloc_gwas_ncontrol = gwas_cfg$n_controls,
        coloc_gwas_N        = gwas_cfg$N,
        coloc_gwas_s        = gwas_cfg$case_fraction,
        stringsAsFactors = FALSE, check.names = FALSE
    )
}


# ------------------------------------------------------------------------------
# Run over the significance table
# ------------------------------------------------------------------------------

sig <- fread(MODE$sig_file, data.table = FALSE, check.names = FALSE)
req_cols <- c("CellType", "GWAS_CHR", "GWAS_BP", "GWAS_Source", "gene_id")
if (isTRUE(MODE$apply_cis_window_check))
    req_cols <- c(req_cols, "gene_chr", "gene_start", "gene_end")
stopifnot(all(req_cols %in% names(sig)))

sig$.row <- seq_len(nrow(sig))
todo <- sig[as.character(sig$CellType) %in% CTs, , drop = FALSE]

## Optional third CLI arg for testing one or several lead SNPs.
## Example: Rscript GWAS_leadvar_coloc_protein_fast.R transprotein ABS rs112191583
## A targeted run writes a *separate* file; it cannot overwrite the full output.
if (length(args) >= 3 && nzchar(args[3]) && !identical(args[3], "all")) {
    chosen_rsids <- trimws(strsplit(args[3], ",", fixed = TRUE)[[1]])
    todo <- todo[tolower(as.character(todo$GWAS_RSID)) %in% tolower(chosen_rsids),
                 , drop = FALSE]
    if (!nrow(todo)) stop("No tasks match selected GWAS_RSID: ", args[3])
    safe_label <- gsub("[^A-Za-z0-9_-]", "_", args[3])
    MODE$out_file <- sub("\\.csv$", paste0(".", safe_label, ".csv"), MODE$out_file)
}

message("Mode                 : ", MODE$name, " (one QTL parquet per ", MODE$one_file_per, ")")
message("Significant ", MODE$label, " rows : ", nrow(sig))
message("Rows to test         : ", nrow(todo), " (", paste(CTs, collapse = ", "), ")")
message("Flank                : +/-", CFG$flank_length, " bp around the lead variant")
message("GWAS                 : Multi N=", CFG$gwas$Multi$N,
        " (s=", signif(CFG$gwas$Multi$case_fraction, 4), ") | EUR N=", CFG$gwas$EUR$N,
        " (s=", signif(CFG$gwas$EUR$case_fraction, 4), ")")
message("Output               : ", MODE$out_file)

## (cell type, chromosome) order keeps one GWAS file and one QTL file in memory
todo <- todo[order(as.character(todo$CellType),
                   as.integer(normalize_chr(todo$GWAS_CHR))), , drop = FALSE]

## widest window needed per (cell type, gene, chromosome); used in trans mode to
## cut each Arrow batch down before anything is kept.
want_all <- data.frame(
    CellType = as.character(todo$CellType),
    gene_id  = as.character(todo$gene_id),
    chrom    = normalize_chr(todo$GWAS_CHR),
    lo       = safe_numeric(todo$GWAS_BP) - CFG$flank_length,
    hi       = safe_numeric(todo$GWAS_BP) + CFG$flank_length,
    stringsAsFactors = FALSE
)
want_all$key <- paste(canonical_gene_id(want_all$gene_id), want_all$chrom, sep = "|")
want_all <- do.call(rbind, lapply(
    split(want_all, paste(want_all$CellType, want_all$key)), function(d)
        data.frame(CellType = d$CellType[1], gene_id = d$gene_id[1], chrom = d$chrom[1],
                   lo = min(d$lo), hi = max(d$hi), key = d$key[1],
                   stringsAsFactors = FALSE)))

res <- vector("list", nrow(todo))

for (i in seq_len(nrow(todo))) {
    row  <- todo[i, , drop = FALSE]
    want <- want_all[want_all$CellType == as.character(row$CellType), , drop = FALSE]
    message("[", i, "/", nrow(todo), "] ", row$CellType, " ", row$gene_id,
            " x ", row$GWAS_RSID, " (chr", row$GWAS_CHR, ":", row$GWAS_BP, ")")
    r <- tryCatch(run_one_coloc(row, want),
                  error = function(e) NA_RESULT(paste0("ERROR: ", conditionMessage(e))))
    r$.row <- row$.row
    res[[i]] <- r
    message("    ", r$coloc_status,
            if (identical(r$coloc_status, "coloc_completed"))
                paste0(": nsnps=", r$coloc_nsnps,
                       " PP.H4=", signif(r$coloc_PP.H4.abf, 3)) else "")
}

res <- do.call(rbind, res)

## append the coloc columns to the input table, original row order preserved
out <- merge(sig, res, by = ".row", all.x = TRUE, sort = FALSE)
out <- out[order(out$.row), , drop = FALSE]
out$.row <- NULL

fwrite(out, MODE$out_file)

message("\nWrote ", nrow(out), " rows -> ", MODE$out_file)
message("completed: ", sum(out$coloc_status == "coloc_completed", na.rm = TRUE),
        " | PP.H4 > 0.8: ",
        sum(out$coloc_PP.H4.abf > 0.8, na.rm = TRUE),
        " | PP.H4 > 0.5: ",
        sum(out$coloc_PP.H4.abf > 0.5, na.rm = TRUE))
