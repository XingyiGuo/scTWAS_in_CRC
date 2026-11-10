# conda activate tftwas_coloc
# ==============================================================================
# GWAS_leadvar_SMR_HEIDI_inputs.R
#
# SMR + HEIDI inputs for every significant lead-variant x gene pair, in the same
# two modes as GWAS_leadvar_coloc_protein.R:
#   cisprotein    AllCellTypes_cis_BFsig.csv   + <PREFIX>.cis_qtl_pairs.chr<N>.parquet
#   transprotein  AllCellTypes_trans_BFsig.csv + <PREFIX>.GWAS_leadvar_trans.<CT>.trans_qtl_pairs.parquet
#   all           both, one manifest per mode plus a combined manifest
#
# Per task: region = GWAS lead variant +/- flank_length (not the gene TSS); GWAS
# and QTL harmonised to the GWAS A1 exactly as in the coloc script; target SNP =
# the harmonised SNP with the smallest eQTL p in the region, forced later with
# --target-snp. Writes <base>.gwas_500kb.ma, .eqtl_500kb.esd, .eqtl_500kb.flist,
# .harmonized_snps.txt, .target_snp_qc.tsv and one manifest row per task, which
# GWAS_leadvar_run_SMR_HEIDI.sh consumes.
#
# PROBE PLACEMENT: the BESD probe is annotated at the LEAD VARIANT, not the gene
# TSS, so SMR's --cis-wind window coincides with the region and the same code
# path serves trans pairs whose gene is on another chromosome. It affects
# neither the SMR statistic (target SNP only) nor HEIDI (region SNPs and LD),
# but Probe_bp / ProbeChr in the .smr output are the lead variant's coordinates;
# the gene's own position travels as gene_chr / gene_tss / gene_lead_distance.
#
# QTL INPUT: the tensorQTL parquet, NOT <CT>...sig{cis,trans}.tensorqtl.bed.gz,
# which is the phenotype matrix and carries no betas. tensorQTL --mode trans
# also drops pairs within +/-5 Mb of the gene, so such tasks are reported as
# inside_tensorqtl_cis_window_5Mb and belong to cis mode; the trans run must use
# TRANS_PVAL_THRESHOLD=1 or too few region SNPs survive.
#
# Caches (all under RES_DIR unless overridden, all signature-checked against the
# source file and the window set, so a changed table or flank forces a rescan):
#   gwas_window_cache/gwas_<cohort>_chr<N>.rds   window-trimmed GWAS
#   trans_coloc_cache/trans_coloc_<CT>.rds       filtered trans QTL, SAME name
#                                                and signature as the coloc
#                                                script, so one scan serves both
#   SMR_SKIP_EXISTING=1 (default) reuses a task whose five files are on disk and
#   whose QC table matches the same gene, lead variant and region.
#
# Usage:
#   Rscript GWAS_leadvar_SMR_HEIDI_inputs.R cisprotein
#   Rscript GWAS_leadvar_SMR_HEIDI_inputs.R transprotein ABS
#   Rscript GWAS_leadvar_SMR_HEIDI_inputs.R all ABS,GOB <config.txt>
#
# Env overrides: COLOC_RES_DIR, COLOC_GWAS_ROOT, COLOC_QTL_PREFIX,
#   COLOC_QTL_TEMPLATE, SMR_OUTPUT_DIR, SMR_SKIP_EXISTING,
#   COLOC_GWAS_CACHE_DIR, COLOC_TRANS_CACHE_DIR, COLOC_TRANS_BATCH_SIZE
# ==============================================================================

suppressPackageStartupMessages({
    library(data.table)
    library(arrow)
})

args <- if (exists("SMR_ARGS")) SMR_ARGS else commandArgs(trailingOnly = TRUE)

MODE_ARG <- if (length(args) >= 1 && nzchar(args[1])) tolower(trimws(args[1])) else ""
if (!MODE_ARG %in% c("cisprotein", "transprotein", "all", "both"))
    stop("First argument must be 'cisprotein', 'transprotein' or 'all' (got: '",
         MODE_ARG, "').")
MODES_TO_RUN <- if (MODE_ARG %in% c("all", "both"))
    c("cisprotein", "transprotein") else MODE_ARG

# ---- Configuration -------------------------------------------------------------
## Identical sources to GWAS_leadvar_coloc_protein.R; an optional config file (argument 3) may override any element of CFG.

RES_DIR <- Sys.getenv("COLOC_RES_DIR", unset =
    "/data/l2_bioinfo1/liq17/scTWAS/GWAS_leadvar_coloc_SMR_res/Polyps")

GWAS_ROOT <- Sys.getenv("COLOC_GWAS_ROOT", unset =
    "/data/l2_bioinfo1/liq17/Cancer_GWAS_SS/CRC_GWAS_SS/Metal_202608/GC_Final/processed")

## CRC GWAS sample sizes (cases / controls / total), as supplied:
##   Multi : 132,261 / 1,752,398 / 1,884,659
##   EUR   :  99,260 / 1,119,512 / 1,218,772
gwas_cc <- function(n_cases, n_controls, n_total, file_template) {
    stopifnot(n_cases + n_controls == n_total)
    list(file_template = file_template,
         n_cases = n_cases, n_controls = n_controls, N = n_total,
         case_fraction = n_cases / n_total)
}

CFG <- list(

    analysis_group = "Polyps",
    eqtl_source    = "VUMC_HTAN_ADSER",

    ## Same nine ADSER epithelial cell types as LeadVariants_Profile_lm_complete_Polyps.R.
    ## Only ABS, ASC, CT, GOB, STM, TAC carry BF-significant cis rows; the rest
    ## simply yield no tasks.
    cell_types = c("ABS", "CT", "EE", "GOB", "STM", "TAC", "TUF", "ASC", "SSC"),

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

    ## Not used by this script (kept for parity with the coloc config). The lm
    ## sample size is per cell type: 84 (ABS, GOB, TAC), 79 (STM), 70 (CT),
    ## 49 (ASC) in AllCellTypes_cis_BFsig.csv.
    qtl_N          = 84,
    qtl_batch_size = as.integer(Sys.getenv("COLOC_TRANS_BATCH_SIZE", "250000")),

    ## PLINK2 pgen used by the runner to build the per-task HEIDI LD reference.
    ## TODO confirm on the cluster: derived from the .traw prefix used in
    ## LeadVariants_Profile_lm_complete_Polyps.R (genotype/plinkformat/
    ## COLONMAP_95.chr{chrom}.nodup.maf001.traw). Point this at whatever pgen
    ## set the ADSER tensorQTL run used.
    pfile_template = paste0(
        "/data/l2_bioinfo1/liq17/scTWAS/data/VUMC_HTAN/genotype/plinkformat/",
        "COLONMAP_95.chr{chrom}.nodup.maf001"),

    flank_length          = 500000L,
    min_overlap_snps      = 3L,
    drop_strand_ambiguous = TRUE,
    tensorqtl_cis_window  = 5000000,

    ## SMR/HEIDI options. --target-snp forces the instrument, so peqtl_smr does
    ## not select the target; peqtl_heidi still controls HEIDI SNP selection.
    smr = list(
        peqtl_smr      = 5e-8,
        peqtl_heidi    = 1.57e-3,
        heidi_mtd      = 1L,
        heidi_min_m    = 3L,
        heidi_max_m    = 20L,
        ld_lower_limit = 0.05,
        ld_upper_limit = 0.9,
        maf            = 0.01,
        diff_freq      = 0.20,
        diff_freq_prop = 0.05,
        thread_num     = 8L,
        ## Kb, and must cover the whole region: the probe is annotated at the
        ## lead variant, so flank_length/1000 is exactly right.
        cis_wind       = 500L
    ),

    output_dir = Sys.getenv("SMR_OUTPUT_DIR", unset = file.path(RES_DIR, "SMR_HEIDI_leadvar_500kb"))
)

CT_ARG <- if (length(args) >= 2 && nzchar(args[2]) && !identical(tolower(args[2]), "all"))
    trimws(strsplit(args[2], ",")[[1]]) else NULL

CONFIG_FILE <- if (length(args) >= 3 && nzchar(args[3])) args[3] else ""
if (nzchar(CONFIG_FILE)) {
    stopifnot(file.exists(CONFIG_FILE))
    message("Config file: ", CONFIG_FILE)
    source(CONFIG_FILE, local = TRUE)
    if (exists("SMR_CONFIG", inherits = FALSE)) CFG <- modifyList(CFG, SMR_CONFIG)
}

CTs <- if (!is.null(CT_ARG)) CT_ARG else CFG$cell_types

## Mode table: the mode picks the significance file, the QTL file layout, the
## tensorQTL column spelling and the reader -- same contract as the coloc script.
make_mode <- function(mode_name) {
    if (identical(mode_name, "cisprotein")) {
        list(
            name  = "cisprotein", label = "cis",
            sig_file = file.path(RES_DIR, "AllCellTypes_cis_BFsig.csv"),
            prefix_template = Sys.getenv("COLOC_QTL_PREFIX",
                unset = "VUMC_HTAN_ADSER_{CT}"),
            file_template = Sys.getenv("COLOC_QTL_TEMPLATE", unset = paste0(
                "/data/l2_bioinfo1/liq17/scTWAS/data/VUMC_HTAN/xQTL_VUMC_HTAN/",
                "{PREFIX}/{PREFIX}.cis_qtl_pairs.chr{chrom}.parquet")),
            col_beta = "slope", col_se = "slope_se", col_p = "pval_nominal",
            col_optional = c("tss_distance", "af", "ma_samples", "ma_count"),
            one_file_per = "chromosome",
            apply_cis_window_check = FALSE,
            reader = function(CT, chrom, want) read_qtl_cis(CT, chrom, want)
        )
    } else {
        list(
            name  = "transprotein", label = "trans",
            sig_file = file.path(RES_DIR, "AllCellTypes_trans_BFsig.csv"),
            prefix_template = Sys.getenv("COLOC_QTL_PREFIX",
                unset = "VUMC_HTAN_ADSER_{CT}"),
            file_template = Sys.getenv("COLOC_QTL_TEMPLATE", unset = file.path(
                "/data/l2_bioinfo1/liq17/scTWAS/GWAS_leadvar_coloc_SMR_res/Polyps",
                "{PREFIX}.GWAS_leadvar_trans.{CT}.trans_qtl_pairs.parquet")),
            col_beta = "b", col_se = "b_se", col_p = "pval",
            col_optional = c("af", "r2"),
            one_file_per = "cell type",
            apply_cis_window_check = TRUE,
            reader = function(CT, chrom, want) read_qtl_trans(CT, want)
        )
    }
}

# ---- Helpers (identical to GWAS_leadvar_coloc_protein.R) -----------------------

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
## Fully vectorised: data.table::tstrsplit splits the whole column in C, instead
## of one strsplit list element plus four vapply passes per variant. This is the
## single hottest line in the script on a multi-million-row parquet.
parse_tensorqtl_variant <- function(variant_id) {
    x <- as.character(variant_id)
    if (!length(x))
        return(data.frame(SNP = character(0), variant_chr = character(0),
                          variant_pos = numeric(0), REF = character(0),
                          ALT = character(0), stringsAsFactors = FALSE))

    has <- grepl(";", x, fixed = TRUE)
    if (all(has)) {                       # the usual case: one split, no ifelse
        sp   <- data.table::tstrsplit(x, ";", fixed = TRUE, keep = 1:2)
        rsid <- sp[[1]]; key <- sp[[2]]
    } else {
        rsid <- ifelse(has, sub(";.*$", "", x), NA_character_)
        key  <- ifelse(has, sub("^[^;]+;", "", x), x)
    }

    p <- data.table::tstrsplit(key, "_", fixed = TRUE, keep = 1:4)
    if (length(p) < 4L) p <- c(p, rep(list(rep(NA_character_, length(x))), 4L - length(p)))

    data.frame(
        SNP         = tolower(trimws(rsid)),
        variant_chr = normalize_chr(p[[1]]),
        variant_pos = safe_numeric(p[[2]]),
        REF         = toupper(trimws(p[[3]])),
        ALT         = toupper(trimws(p[[4]])),
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

ensure_dir <- function(x) {
    if (!dir.exists(x)) dir.create(x, recursive = TRUE, showWarnings = FALSE)
    x
}

sanitize_component <- function(x) gsub("[^A-Za-z0-9._-]+", "_", as.character(x))

abs_path <- function(x) {
    y <- suppressWarnings(normalizePath(x, winslash = "/", mustWork = FALSE))
    if (is.na(y) || !nzchar(y)) x else y
}

## keep manifest rows valid TSV records even when a message contains a tab
clean_message <- function(x) {
    if (is.null(x) || all(is.na(x))) return(NA_character_)
    gsub("[\t\r\n]+", " ", as.character(x))
}

GWAS_CACHE <- new.env(parent = emptyenv())
QTL_CACHE  <- new.env(parent = emptyenv())

## Union of the lead-variant windows per (cohort, chromosome), registered before
## the task loop so each GWAS file is read once and only in-window rows kept.
GWAS_WINDOWS <- new.env(parent = emptyenv())

set_gwas_windows <- function(tasks) {
    coh <- unname(CFG$source_to_cohort[as.character(tasks$GWAS_Source)])
    ok  <- !is.na(coh)
    if (!any(ok)) return(invisible(NULL))
    d <- data.frame(
        key   = paste(coh[ok], normalize_chr(tasks$GWAS_CHR)[ok], sep = "|"),
        lo    = safe_numeric(tasks$GWAS_BP)[ok] - CFG$flank_length,
        hi    = safe_numeric(tasks$GWAS_BP)[ok] + CFG$flank_length,
        stringsAsFactors = FALSE)
    for (k in unique(d$key)) {
        w <- d[d$key == k, c("lo", "hi"), drop = FALSE]
        if (exists(k, envir = GWAS_WINDOWS, inherits = FALSE))
            w <- rbind(get(k, envir = GWAS_WINDOWS, inherits = FALSE), w)
        w <- w[order(w$lo), , drop = FALSE]
        ## merge overlaps so the row filter is a single findInterval
        grp <- cumsum(c(TRUE, w$lo[-1] > cummax(w$hi)[-nrow(w)]))
        new_w <- data.frame(lo = as.numeric(tapply(w$lo, grp, min)),
                            hi = as.numeric(tapply(w$hi, grp, max)))
        old_w <- if (exists(k, envir = GWAS_WINDOWS, inherits = FALSE))
            get(k, envir = GWAS_WINDOWS, inherits = FALSE) else NULL
        ## a cached subset built for a narrower window set would be incomplete
        if (!is.null(old_w) && !isTRUE(all.equal(old_w, new_w)) &&
            exists(k, envir = GWAS_CACHE, inherits = FALSE))
            rm(list = k, envir = GWAS_CACHE)
        assign(k, new_w, envir = GWAS_WINDOWS)
    }
    invisible(NULL)
}

## TRUE when every needed interval sits inside one of `stored`, so a cache built
## for a wider window set is reused as is (extra rows are harmless: every task
## filters to its own window).
windows_covered <- function(need, stored) {
    if (is.null(stored) || !nrow(stored)) return(FALSE)
    if (is.null(need) || !nrow(need)) return(TRUE)
    j <- findInterval(need$lo, stored$lo)
    all(j >= 1L & need$hi <= stored$hi[pmax(j, 1L)])
}

in_windows <- function(pos, w) {
    if (is.null(w) || !nrow(w)) return(rep(TRUE, length(pos)))
    i <- findInterval(pos, w$lo)
    ok <- i >= 1L & is.finite(pos)
    ok[ok] <- pos[ok] <= w$hi[i[ok]]
    ok
}

read_gwas_chr <- function(cohort, chrom) {
    key <- paste(cohort, chrom, sep = "|")
    if (exists(key, envir = GWAS_CACHE, inherits = FALSE))
        return(get(key, envir = GWAS_CACHE, inherits = FALSE))

    path <- fill_template(CFG$gwas[[cohort]]$file_template, list(chrom = chrom))
    if (!file.exists(path)) stop("GWAS file not found: ", path)

    w <- if (exists(key, envir = GWAS_WINDOWS, inherits = FALSE))
        get(key, envir = GWAS_WINDOWS, inherits = FALSE) else NULL

    ## Window-trimmed GWAS cached on disk: a chromosome costs seconds to minutes
    ## to parse while the trimmed result is a few thousand rows.
    info <- file.info(path)
    source_sig <- c(normalizePath(path, mustWork = TRUE),
                    as.character(info$size), as.character(as.numeric(info$mtime)))
    gcache_dir <- Sys.getenv("COLOC_GWAS_CACHE_DIR",
                             unset = file.path(RES_DIR, "gwas_window_cache"))
    use_gcache <- nzchar(gcache_dir)
    gcache_file <- NULL
    if (use_gcache) {
        dir.create(gcache_dir, recursive = TRUE, showWarnings = FALSE)
        gcache_file <- file.path(gcache_dir,
                                 paste0("gwas_", cohort, "_chr", chrom, ".rds"))
        if (file.exists(gcache_file)) {
            cached <- tryCatch(readRDS(gcache_file), error = function(e) NULL)
            if (is.list(cached) && identical(cached$source_sig, source_sig) &&
                is.data.frame(cached$out) &&
                all(c("SNP", "POS", "A1", "A2", "BETA", "SE", "P", "FRQ") %in%
                    names(cached$out)) &&
                windows_covered(w, cached$windows)) {
                message("  ", cohort, " chr", chrom,
                        " GWAS loaded from cache: ", nrow(cached$out), " variants")
                assign(key, cached$out, envir = GWAS_CACHE)
                return(cached$out)
            }
        }
    }

    g <- fread(path, data.table = FALSE, check.names = FALSE,
               nThread = getDTthreads())

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

    ## Keep only the lead-variant windows this run will ever ask for. Trimming
    ## AFTER the dedup leaves the retained SNP set identical to the unfiltered
    ## version; the subset is then small enough for every chromosome to stay
    ## cached for the whole run.
    n_all <- nrow(out)
    out <- out[in_windows(out$POS, w), , drop = FALSE]

    if (use_gcache)
        saveRDS(list(source_sig = source_sig, windows = w, out = out),
                gcache_file, compress = FALSE)

    assign(key, out, envir = GWAS_CACHE)
    message("  ", cohort, " chr", chrom, " GWAS variants: ", nrow(out),
            " (of ", n_all, " QC-passing on the chromosome)")
    out
}

## ---- shared by both readers ------------------------------------------------

empty_qtl <- function() data.frame(
    phenotype_id = character(0), phenotype_key = character(0),
    variant_id = character(0), SNP = character(0), variant_chr = character(0),
    variant_pos = numeric(0), REF = character(0), ALT = character(0),
    BETA = numeric(0), SE = numeric(0), P = numeric(0), AF = numeric(0),
    stringsAsFactors = FALSE)

qtl_read_cols <- function(ds, path, label) {
    have <- as.character(ds$schema$names)
    need <- c("phenotype_id", "variant_id", MODE$col_beta, MODE$col_se, MODE$col_p)
    miss <- setdiff(need, have)
    if (length(miss))
        stop(label, " parquet missing column(s): ", paste(miss, collapse = ", "),
             "\nFile: ", path, "\nAvailable: ", paste(have, collapse = ", "),
             if (identical(label, "trans-QTL"))
                 "\n(cis-style slope/slope_se/pval_nominal means this is a cis file)"
             else "")
    c(need, intersect(MODE$col_optional, have))
}

## only the needed phenotypes are decoded
qtl_gene_filter <- function(gene_ids) {
    genes <- unique(c(as.character(gene_ids), canonical_gene_id(gene_ids)))
    arrow::Expression$create("is_in", arrow::Expression$field_ref("phenotype_id"),
        options = list(value_set = arrow::Array$create(genes), skip_nulls = TRUE))
}

## map the mode's column spelling onto BETA/SE/P/AF, then the shared QC
qtl_finalise <- function(eq) {
    eq$BETA <- safe_numeric(eq[[MODE$col_beta]])
    eq$SE   <- safe_numeric(eq[[MODE$col_se]])
    eq$P    <- safe_numeric(eq[[MODE$col_p]])
    eq$AF   <- if ("af" %in% names(eq)) safe_numeric(eq$af) else rep(NA_real_, nrow(eq))
    clean_qtl(eq)
}

## keep only variants inside the widest window needed for their gene
qtl_in_gene_windows <- function(eq, w) {
    if (is.null(w) || !nrow(w) || !nrow(eq)) return(eq)
    m <- match(paste(eq$phenotype_key, eq$variant_chr, sep = "|"), w$key)
    eq[!is.na(m) & is.finite(eq$variant_pos) &
       eq$variant_pos >= w$lo[m] & eq$variant_pos <= w$hi[m], , drop = FALSE]
}

## MODE = cisprotein. One parquet per chromosome; only this cell type's genes on
## this chromosome are read, so a whole-chromosome file never lands in R.
read_qtl_cis <- function(CT, chrom, want) {
    key <- paste(CT, chrom, sep = "|")
    if (exists(key, envir = QTL_CACHE, inherits = FALSE))
        return(get(key, envir = QTL_CACHE, inherits = FALSE))
    rm(list = ls(envir = QTL_CACHE), envir = QTL_CACHE)     # single slot

    path <- qtl_path(CT, chrom)
    if (!file.exists(path)) stop("cis-QTL parquet not found: ", path)

    ds   <- arrow::open_dataset(path, format = "parquet")
    cols <- qtl_read_cols(ds, path, "cis-QTL")
    w    <- want[want$chrom == chrom, , drop = FALSE]

    eq <- as.data.frame(arrow::Scanner$create(
        ds, filter = qtl_gene_filter(w$gene_id), projection = cols)$ToTable())
    eq$phenotype_key <- canonical_gene_id(eq$phenotype_id)
    eq <- cbind(eq, parse_tensorqtl_variant(eq$variant_id))
    eq <- qtl_finalise(qtl_in_gene_windows(eq, w))

    assign(key, eq, envir = QTL_CACHE)
    message("  ", CT, " chr", chrom, " cis-QTL rows in windows: ", nrow(eq))
    eq
}

## MODE = transprotein. One parquet per cell type, far too large to materialise:
## the scan is streamed, and its filtered result is cached on disk under the same
## name and signature the coloc script uses, so whichever script scans first
## serves the other.
read_qtl_trans <- function(CT, want) {
    if (exists(CT, envir = QTL_CACHE, inherits = FALSE))
        return(get(CT, envir = QTL_CACHE, inherits = FALSE))
    rm(list = ls(envir = QTL_CACHE), envir = QTL_CACHE)     # single slot

    path <- qtl_path(CT, NA)
    if (!file.exists(path)) stop("trans-QTL parquet not found: ", path)

    info <- file.info(path)
    source_sig <- c(normalizePath(path, mustWork = TRUE),
                    as.character(info$size), as.character(as.numeric(info$mtime)))
    want_sig <- paste(sort(paste(want$key, want$lo, want$hi, sep = ":")),
                      collapse = "\n")
    cache_dir <- Sys.getenv("COLOC_TRANS_CACHE_DIR",
                            unset = file.path(RES_DIR, "trans_coloc_cache"))
    cache_file <- if (nzchar(cache_dir)) {
        dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
        file.path(cache_dir, paste0("trans_coloc_", CT, ".rds"))
    } else NULL

    if (!is.null(cache_file) && file.exists(cache_file)) {
        cached <- tryCatch(readRDS(cache_file), error = function(e) NULL)
        if (is.list(cached) && identical(cached$source_sig, source_sig) &&
            identical(cached$want_sig, want_sig) && is.data.frame(cached$eq) &&
            all(names(empty_qtl()) %in% names(cached$eq))) {
            message("  ", CT, " loaded filtered trans-QTL cache: ",
                    nrow(cached$eq), " rows")
            assign(CT, cached$eq, envir = QTL_CACHE)
            return(cached$eq)
        }
    }

    t0 <- proc.time()[[3]]
    message("  ", CT, " first pass through trans-QTL parquet: ", path)
    ds   <- arrow::open_dataset(path, format = "parquet")
    cols <- qtl_read_cols(ds, path, "trans-QTL")
    rbr  <- arrow::Scanner$create(ds, filter = qtl_gene_filter(want$gene_id),
                                  projection = cols,
                                  batch_size = CFG$qtl_batch_size)$ToRecordBatchReader()

    acc <- list(); n_seen <- 0; n_kept <- 0; n_batch <- 0L
    repeat {
        batch <- rbr$read_next_batch()
        if (is.null(batch)) break
        n_batch <- n_batch + 1L
        if (!batch$num_rows) next
        n_seen <- n_seen + batch$num_rows

        ## Two-stage parse: chromosome and position for the whole batch, because
        ## that is all the window test needs; rsID, alleles and the data.frame
        ## conversion only for the rows that survive it.
        vid <- batch[["variant_id"]]$as_vector()
        sp  <- data.table::tstrsplit(vid, ";", fixed = TRUE, keep = 2)[[1]]
        nos <- is.na(sp)
        if (any(nos)) sp[nos] <- vid[nos]            # id without an rsID prefix
        p2  <- data.table::tstrsplit(sp, "_", fixed = TRUE, keep = 1:2)
        phk <- canonical_gene_id(batch[["phenotype_id"]]$as_vector())
        pos <- safe_numeric(p2[[2]])

        m  <- match(paste(phk, normalize_chr(p2[[1]]), sep = "|"), want$key)
        ok <- !is.na(m) & is.finite(pos) & pos >= want$lo[m] & pos <= want$hi[m]

        if (any(ok)) {
            idx <- which(ok)
            eqb <- as.data.frame(batch$Take(idx - 1L))      # Take is 0-based
            eqb$phenotype_key <- phk[idx]
            acc[[length(acc) + 1L]] <- cbind(eqb, parse_tensorqtl_variant(eqb$variant_id))
            n_kept <- n_kept + length(idx)
        }
        if (n_batch == 1L || n_batch %% 10L == 0L)
            message("  ", CT, " scan batch ", n_batch, ": gene-filtered rows ",
                    n_seen, ", regional rows ", n_kept, ", elapsed ",
                    round(proc.time()[[3]] - t0), "s")
    }

    eq <- if (length(acc)) qtl_finalise(do.call(rbind, acc)) else empty_qtl()
    message("  ", CT, " scan complete: ", n_seen, " rows after the gene filter; ",
            nrow(eq), " after windows and QC; elapsed ",
            round(proc.time()[[3]] - t0), "s")

    if (!is.null(cache_file)) {
        saveRDS(list(source_sig = source_sig, want_sig = want_sig, eq = eq),
                cache_file, compress = FALSE)
        message("  ", CT, " saved filtered trans-QTL cache: ", cache_file)
    }
    assign(CT, eq, envir = QTL_CACHE)
    eq
}

# ---- Harmonisation: merge on rsID, align the QTL beta AND the QTL frequency to the GWAS A1 ----
## Same allele logic as the coloc script, plus the A1-oriented frequencies and the alignment label that the .ma / .esd files need.

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
    same      <- a1x == a1g & a2x == a2g
    flip      <- a1x == a2g & a2x == a1g
    comp_same <- a1x == complement_allele(a1g) & a2x == complement_allele(a2g)
    comp_flip <- a1x == complement_allele(a2g) & a2x == complement_allele(a1g)
    beta_flip <- flip | comp_flip

    x$alignment <- ifelse(same, "same",
                   ifelse(flip, "flip", ifelse(comp_same, "comp_same", "comp_flip")))

    ## align the QTL effect to the GWAS A1
    x$BETA_x_aligned <- x$BETA_x
    x$BETA_x_aligned[beta_flip] <- -x$BETA_x_aligned[beta_flip]

    ## tensorQTL af is the ALT frequency; after flipping, the aligned effect
    ## allele is REF (= GWAS A1), so the frequency becomes 1-af.
    x$AF_x_aligned <- x$AF
    idx <- beta_flip & is.finite(x$AF_x_aligned)
    x$AF_x_aligned[idx] <- 1 - x$AF_x_aligned[idx]

    ## prefer each study's own frequency, fall back to the other study only when
    ## the required one is absent or out of range (GWAS FRQ is "not_in_1KG" for
    ## variants absent from 1000G and arrives as NA).
    x$GWAS_A1_freq <- x$FRQ
    bad_g <- !is.finite(x$GWAS_A1_freq) | x$GWAS_A1_freq <= 0 | x$GWAS_A1_freq >= 1
    x$GWAS_A1_freq[bad_g] <- x$AF_x_aligned[bad_g]

    x$EQTL_A1_freq <- x$AF_x_aligned
    bad_x <- !is.finite(x$EQTL_A1_freq) | x$EQTL_A1_freq <= 0 | x$EQTL_A1_freq >= 1
    x$EQTL_A1_freq[bad_x] <- x$GWAS_A1_freq[bad_x]

    x <- x[is.finite(x$GWAS_A1_freq) & x$GWAS_A1_freq > 0 & x$GWAS_A1_freq < 1 &
           is.finite(x$EQTL_A1_freq) & x$EQTL_A1_freq > 0 & x$EQTL_A1_freq < 1 &
           is.finite(x$BETA_g) & is.finite(x$SE_g) & x$SE_g > 0 &
           is.finite(x$P_g) & x$P_g > 0 & x$P_g <= 1 &
           is.finite(x$BETA_x_aligned) & is.finite(x$SE_x) & x$SE_x > 0 &
           is.finite(x$P_x) & x$P_x > 0 & x$P_x <= 1,
           , drop = FALSE]
    if (!nrow(x)) return(x)

    ## target selection depends on this order: smallest eQTL p first
    x <- x[order(x$P_x, x$P_g, x$SNP), , drop = FALSE]
    x[!duplicated(x$SNP), , drop = FALSE]
}

# ---- One task: write the SMR inputs for one significant lead-variant x gene pair ----

## the target SNP as it appears in the manifest; the QC table writes the same
## numbers under its own (longer) names
target_fields <- function(target) list(
    target_SNP = target$SNP, target_pos = target$variant_pos,
    target_eQTL_p = target$P_x, target_eQTL_beta = target$BETA_x_aligned,
    target_eQTL_se = target$SE_x, target_eQTL_freq = target$EQTL_A1_freq,
    target_GWAS_p = target$P_g, target_GWAS_beta = target$BETA_g,
    target_GWAS_se = target$SE_g, target_GWAS_freq = target$GWAS_A1_freq,
    target_allele_alignment = target$alignment)

manifest_row <- function(row, status, message_text = NA_character_, extra = list()) {
    cohort <- unname(CFG$source_to_cohort[as.character(row$GWAS_Source)])
    base <- list(
        analysis_group = CFG$analysis_group,
        mode           = MODE$name,
        cell_type      = as.character(row$CellType),
        gene           = as.character(row$gene_id),
        chromosome     = normalize_chr(row$GWAS_CHR),
        cohort         = if (is.na(cohort)) NA_character_ else cohort,
        eqtl_source    = CFG$eqtl_source,
        GWAS_RSID      = as.character(row$GWAS_RSID),
        GWAS_BP        = safe_numeric(row$GWAS_BP),
        GWAS_P         = safe_numeric(row$GWAS_P),
        GWAS_Source    = as.character(row$GWAS_Source),
        lm_beta        = safe_numeric(row$beta),
        lm_p           = safe_numeric(row$p),
        gene_chr       = normalize_chr(row$gene_chr),
        gene_tss       = safe_numeric(row$gene_start),
        gene_lead_distance = NA_real_,
        probe_bp       = NA_real_,
        region_start   = NA_real_, region_end = NA_real_,
        n_gwas_region  = NA_integer_, n_qtl_region = NA_integer_,
        n_overlap_harmonized = NA_integer_,
        target_SNP = NA_character_, target_pos = NA_real_,
        target_eQTL_p = NA_real_, target_eQTL_beta = NA_real_,
        target_eQTL_se = NA_real_, target_eQTL_freq = NA_real_,
        target_GWAS_p = NA_real_, target_GWAS_beta = NA_real_,
        target_GWAS_se = NA_real_, target_GWAS_freq = NA_real_,
        target_allele_alignment = NA_character_,
        gwas_ma = NA_character_, eqtl_esd = NA_character_,
        eqtl_flist = NA_character_, beqtl_prefix = NA_character_,
        snplist = NA_character_, target_qc = NA_character_,
        pfile_prefix = NA_character_, task_out_prefix = NA_character_,
        peqtl_smr = CFG$smr$peqtl_smr, peqtl_heidi = CFG$smr$peqtl_heidi,
        heidi_mtd = CFG$smr$heidi_mtd, heidi_min_m = CFG$smr$heidi_min_m,
        heidi_max_m = CFG$smr$heidi_max_m,
        ld_lower_limit = CFG$smr$ld_lower_limit,
        ld_upper_limit = CFG$smr$ld_upper_limit,
        maf = CFG$smr$maf, diff_freq = CFG$smr$diff_freq,
        diff_freq_prop = CFG$smr$diff_freq_prop,
        thread_num = CFG$smr$thread_num, cis_wind = CFG$smr$cis_wind,
        status = status, message = clean_message(message_text)
    )
    base[names(extra)] <- extra
    as.data.frame(base, stringsAsFactors = FALSE, check.names = FALSE)
}

## Reuse is only safe if the files describe the same pair and region, so the QC
## table is checked, not just the file names.
SKIP_EXISTING <- !identical(Sys.getenv("SMR_SKIP_EXISTING", unset = "1"), "0")

task_paths <- function(row, cohort) {
    CT   <- as.character(row$CellType)
    gene <- as.character(row$gene_id)
    task_dir <- file.path(CFG$output_dir, MODE$name, CT, cohort,
                          sanitize_component(gene))
    base_name <- paste(CFG$analysis_group, MODE$label, CT, cohort,
                       sanitize_component(gene),
                       sanitize_component(as.character(row$GWAS_RSID)), sep = ".")
    list(
        dir             = task_dir,
        base            = base_name,
        gwas_ma         = file.path(task_dir, paste0(base_name, ".gwas_500kb.ma")),
        eqtl_esd        = file.path(task_dir, paste0(base_name, ".eqtl_500kb.esd")),
        eqtl_flist      = file.path(task_dir, paste0(base_name, ".eqtl_500kb.flist")),
        beqtl_prefix    = file.path(task_dir, paste0(base_name, ".eqtl_500kb")),
        snplist         = file.path(task_dir, paste0(base_name, ".harmonized_snps.txt")),
        target_qc       = file.path(task_dir, paste0(base_name, ".target_snp_qc.tsv")),
        task_out_prefix = file.path(task_dir, paste0(base_name, ".SMR_HEIDI"))
    )
}

reuse_prepared_task <- function(row, cohort, geom) {
    tp <- task_paths(row, cohort)
    need <- c(tp$gwas_ma, tp$eqtl_esd, tp$eqtl_flist, tp$snplist, tp$target_qc)
    if (!all(file.exists(need)) || any(file.info(need)$size <= 0)) return(NULL)

    qc <- tryCatch(fread(tp$target_qc, data.table = FALSE, check.names = FALSE),
                   error = function(e) NULL)
    if (is.null(qc) || !nrow(qc)) return(NULL)
    same_region <-
        identical(as.character(qc$gene[1]), as.character(row$gene_id)) &&
        identical(as.character(qc$GWAS_RSID[1]), as.character(row$GWAS_RSID)) &&
        isTRUE(all.equal(safe_numeric(qc$region_start[1]), geom$region_start)) &&
        isTRUE(all.equal(safe_numeric(qc$region_end[1]),   geom$region_end))
    if (!same_region) return(NULL)

    pfile_prefix <- fill_template(CFG$pfile_template,
                                  list(chrom = normalize_chr(row$GWAS_CHR),
                                       CT = as.character(row$CellType)))
    have_pfile <- all(file.exists(paste0(pfile_prefix, c(".pgen", ".pvar", ".psam"))))

    manifest_row(row,
        status = if (have_pfile) "prepared" else "prepared_missing_pfile_prefix",
        message_text = "reused_existing_inputs",
        extra = c(geom, list(
            n_gwas_region = safe_numeric(qc$n_gwas_region[1]),
            n_qtl_region = safe_numeric(qc$n_qtl_region[1]),
            n_overlap_harmonized = safe_numeric(qc$n_overlap_harmonized[1]),
            target_SNP = as.character(qc$target_SNP[1]),
            target_pos = safe_numeric(qc$target_pos[1]),
            target_eQTL_p = safe_numeric(qc$target_eQTL_p[1]),
            target_eQTL_beta = safe_numeric(qc$target_eQTL_beta_aligned_A1[1]),
            target_eQTL_se = safe_numeric(qc$target_eQTL_se[1]),
            target_eQTL_freq = safe_numeric(qc$target_eQTL_A1_freq[1]),
            target_GWAS_p = safe_numeric(qc$target_GWAS_p[1]),
            target_GWAS_beta = safe_numeric(qc$target_GWAS_beta[1]),
            target_GWAS_se = safe_numeric(qc$target_GWAS_se[1]),
            target_GWAS_freq = safe_numeric(qc$target_GWAS_A1_freq[1]),
            target_allele_alignment = as.character(qc$allele_alignment[1]),
            gwas_ma = abs_path(tp$gwas_ma), eqtl_esd = abs_path(tp$eqtl_esd),
            eqtl_flist = abs_path(tp$eqtl_flist),
            beqtl_prefix = abs_path(tp$beqtl_prefix),
            snplist = abs_path(tp$snplist), target_qc = abs_path(tp$target_qc),
            pfile_prefix = pfile_prefix,
            task_out_prefix = abs_path(tp$task_out_prefix))))
}

prepare_one_task <- function(row, want, read_qtl = MODE$reader) {

    CT      <- as.character(row$CellType)
    gene    <- as.character(row$gene_id)
    chrom   <- normalize_chr(row$GWAS_CHR)
    lead_bp <- safe_numeric(row$GWAS_BP)
    gene_key <- canonical_gene_id(gene)

    cohort <- unname(CFG$source_to_cohort[as.character(row$GWAS_Source)])
    if (is.na(cohort))
        return(manifest_row(row, paste0("unmapped_GWAS_Source:", row$GWAS_Source)))

    region_start <- max(1, lead_bp - CFG$flank_length)
    region_end   <- lead_bp + CFG$flank_length
    gene_tss     <- safe_numeric(row$gene_start)
    gene_lead_distance <- if (identical(normalize_chr(row$gene_chr), chrom))
        abs(lead_bp - gene_tss) else NA_real_
    geom <- list(region_start = region_start, region_end = region_end,
                 gene_lead_distance = gene_lead_distance, probe_bp = lead_bp)

    ## trans mode only: tensorQTL has already dropped every pair within +/-5 Mb
    ## of the gene, so such a window holds no trans QTL rows at all.
    if (isTRUE(MODE$apply_cis_window_check) &&
        identical(normalize_chr(row$gene_chr), chrom)) {
        gl <- safe_numeric(row$gene_start) - CFG$tensorqtl_cis_window
        gh <- safe_numeric(row$gene_end)   + CFG$tensorqtl_cis_window
        if (is.finite(gl) && is.finite(gh) && region_start >= gl && region_end <= gh)
            return(manifest_row(row, "inside_tensorqtl_cis_window_5Mb", NA, geom))
    }

    ## Resume: reuse a task already on disk, manifest row rebuilt from its QC table.
    if (SKIP_EXISTING) {
        reused <- reuse_prepared_task(row, cohort, geom)
        if (!is.null(reused)) return(reused)
    }

    gw <- read_gwas_chr(cohort, chrom)
    gw <- gw[is.finite(gw$POS) & gw$POS >= region_start & gw$POS <= region_end,
             , drop = FALSE]
    geom$n_gwas_region <- nrow(gw)

    eq <- read_qtl(CT, chrom, want)
    eq <- eq[eq$phenotype_key == gene_key, , drop = FALSE]
    if (!nrow(eq)) return(manifest_row(row, "gene_not_found_in_qtl", NA, geom))

    eq <- eq[eq$variant_chr == chrom & is.finite(eq$variant_pos) &
             eq$variant_pos >= region_start & eq$variant_pos <= region_end,
             , drop = FALSE]
    geom$n_qtl_region <- nrow(eq)
    if (!nrow(eq)) return(manifest_row(row, "no_qtl_variant_in_region", NA, geom))

    matched <- harmonize_gwas_qtl(gw, eq)
    geom$n_overlap_harmonized <- nrow(matched)
    if (nrow(matched) < CFG$min_overlap_snps)
        return(manifest_row(row, "too_few_harmonized_overlap_snps",
                            paste0(nrow(matched), " < ", CFG$min_overlap_snps), geom))

    ## harmonize_gwas_qtl returns rows ordered by eQTL p-value, so the target is
    ## the first row: the smallest-p eQTL in the lead-variant region.
    target <- matched[1, , drop = FALSE]

    tp <- task_paths(row, cohort)
    ensure_dir(tp$dir)
    gwas_ma    <- tp$gwas_ma
    eqtl_esd   <- tp$eqtl_esd
    eqtl_flist <- tp$eqtl_flist
    beqtl_prefix <- tp$beqtl_prefix
    snplist    <- tp$snplist
    target_qc  <- tp$target_qc
    task_out_prefix <- tp$task_out_prefix

    ## GWAS .ma format required by SMR / GCTA-COJO: SNP A1 A2 freq b se p n
    fwrite(data.frame(
        SNP = matched$SNP, A1 = matched$A1, A2 = matched$A2,
        freq = matched$GWAS_A1_freq, b = matched$BETA_g, se = matched$SE_g,
        p = matched$P_g, n = CFG$gwas[[cohort]]$N, stringsAsFactors = FALSE),
        gwas_ma, sep = "\t", quote = FALSE, na = "NA")

    ## ESD format required by SMR; Beta/Freq are A1-oriented, as in the .ma
    esd_out <- data.frame(
        Chr = as.integer(chrom), SNP = matched$SNP,
        Bp = as.integer(round(matched$variant_pos)),
        A1 = matched$A1, A2 = matched$A2, Freq = matched$EQTL_A1_freq,
        Beta = matched$BETA_x_aligned, se = matched$SE_x, p = matched$P_x,
        stringsAsFactors = FALSE)
    esd_out <- esd_out[order(esd_out$Bp, esd_out$SNP), , drop = FALSE]
    fwrite(esd_out, eqtl_esd, sep = "\t", quote = FALSE, na = "NA")

    ## BESD probe list; ProbeBp is the lead variant (see the header).
    probe_id <- paste(gene, as.character(row$GWAS_RSID), sep = "__")
    fwrite(data.frame(
        Chr = as.integer(chrom), ProbeID = probe_id, GeneticDistance = 0,
        ProbeBp = as.integer(round(lead_bp)), Gene = gene, Orientation = "+",
        PathOfEsd = abs_path(eqtl_esd), stringsAsFactors = FALSE),
        eqtl_flist, sep = "\t", quote = FALSE, na = "NA")

    fwrite(data.frame(SNP = matched$SNP), snplist, sep = "\t", col.names = FALSE)

    fwrite(data.frame(
        analysis_group = CFG$analysis_group, mode = MODE$name, cell_type = CT,
        gene = gene, probe_id = probe_id, chromosome = chrom, cohort = cohort,
        mode_label = MODE$label,
        eqtl_source = CFG$eqtl_source, GWAS_RSID = as.character(row$GWAS_RSID),
        lead_bp = lead_bp, gene_chr = normalize_chr(row$gene_chr),
        gene_tss = gene_tss, gene_lead_distance = gene_lead_distance,
        region_start = region_start, region_end = region_end,
        target_SNP = target$SNP, target_pos = target$variant_pos,
        GWAS_A1 = target$A1, GWAS_A2 = target$A2,
        original_QTL_REF = target$REF, original_QTL_ALT = target$ALT,
        allele_alignment = target$alignment,
        target_eQTL_beta_aligned_A1 = target$BETA_x_aligned,
        target_eQTL_se = target$SE_x, target_eQTL_p = target$P_x,
        target_eQTL_A1_freq = target$EQTL_A1_freq,
        target_GWAS_beta = target$BETA_g, target_GWAS_se = target$SE_g,
        target_GWAS_p = target$P_g, target_GWAS_A1_freq = target$GWAS_A1_freq,
        n_gwas_region = nrow(gw), n_qtl_region = nrow(eq),
        n_overlap_harmonized = nrow(matched), stringsAsFactors = FALSE),
        target_qc, sep = "\t", quote = FALSE, na = "NA")

    pfile_prefix <- fill_template(CFG$pfile_template, list(chrom = chrom, CT = CT))
    have_pfile <- all(file.exists(paste0(pfile_prefix, c(".pgen", ".pvar", ".psam"))))

    manifest_row(row,
        status = if (have_pfile) "prepared" else "prepared_missing_pfile_prefix",
        message_text = if (have_pfile) NA_character_ else
            paste0("LD reference PGEN not found: ", pfile_prefix),
        extra = c(geom, target_fields(target), list(
            gwas_ma = abs_path(gwas_ma), eqtl_esd = abs_path(eqtl_esd),
            eqtl_flist = abs_path(eqtl_flist), beqtl_prefix = abs_path(beqtl_prefix),
            snplist = abs_path(snplist), target_qc = abs_path(target_qc),
            pfile_prefix = pfile_prefix,
            task_out_prefix = abs_path(task_out_prefix))))
}

# ---- Run the requested modes ---------------------------------------------------

ensure_dir(CFG$output_dir)
all_manifests <- list()

## Register every mode's windows before the first GWAS read, so one cached
## subset per (cohort, chromosome) serves cis and trans alike.
for (mode_name in MODES_TO_RUN) {
    m0 <- make_mode(mode_name)
    if (!file.exists(m0$sig_file)) next
    s0 <- fread(m0$sig_file, data.table = FALSE, check.names = FALSE,
                select = c("CellType", "GWAS_CHR", "GWAS_BP", "GWAS_Source"))
    set_gwas_windows(s0[as.character(s0$CellType) %in% CTs, , drop = FALSE])
}

for (mode_name in MODES_TO_RUN) {

    MODE <<- make_mode(mode_name)
    rm(list = ls(envir = QTL_CACHE), envir = QTL_CACHE)

    stopifnot(file.exists(MODE$sig_file))
    sig <- fread(MODE$sig_file, data.table = FALSE, check.names = FALSE)
    req_cols <- c("CellType", "GWAS_RSID", "GWAS_CHR", "GWAS_BP", "GWAS_Source",
                  "gene_id", "gene_chr", "gene_start", "gene_end")
    stopifnot(all(req_cols %in% names(sig)))

    todo <- sig[as.character(sig$CellType) %in% CTs, , drop = FALSE]
    todo <- todo[order(as.character(todo$CellType),
                       as.integer(normalize_chr(todo$GWAS_CHR))), , drop = FALSE]

    message("\n==============================================================")
    message("Mode                 : ", MODE$name,
            " (one QTL parquet per ", MODE$one_file_per, ")")
    message("Significance table   : ", MODE$sig_file)
    message("Significant ", MODE$label, " rows : ", nrow(sig))
    message("Tasks to prepare     : ", nrow(todo), " (", paste(CTs, collapse = ", "), ")")
    message("Region               : GWAS lead variant +/-", CFG$flank_length, " bp")
    message("Output               : ", file.path(CFG$output_dir, MODE$name))

    ## Polyps has no BF-significant trans pair (AllCellTypes_trans_BFsig.csv is
    ## header-only), and an empty task set would reach fwrite() with a NULL
    ## manifest. Skip the mode instead.
    if (nrow(todo) == 0L) {
        message("No significant ", MODE$label, " pair for the requested cell type(s); mode skipped.")
        next
    }

    ## widest region needed per (cell type, gene, chromosome); trans mode uses it
    ## to cut each Arrow batch down before anything is kept.
    want_all <- data.frame(
        CellType = as.character(todo$CellType),
        gene_id  = as.character(todo$gene_id),
        chrom    = normalize_chr(todo$GWAS_CHR),
        lo       = safe_numeric(todo$GWAS_BP) - CFG$flank_length,
        hi       = safe_numeric(todo$GWAS_BP) + CFG$flank_length,
        stringsAsFactors = FALSE)
    want_all$key <- paste(canonical_gene_id(want_all$gene_id), want_all$chrom, sep = "|")
    want_all <- do.call(rbind, lapply(
        split(want_all, paste(want_all$CellType, want_all$key)), function(d)
            data.frame(CellType = d$CellType[1], gene_id = d$gene_id[1],
                       chrom = d$chrom[1], lo = min(d$lo), hi = max(d$hi),
                       key = d$key[1], stringsAsFactors = FALSE)))

    set_gwas_windows(todo)

    rows <- vector("list", nrow(todo))
    want <- NULL; want_ct <- ""
    for (i in seq_len(nrow(todo))) {
        row  <- todo[i, , drop = FALSE]
        if (!identical(as.character(row$CellType), want_ct)) {
            want_ct <- as.character(row$CellType)
            want <- want_all[want_all$CellType == want_ct, , drop = FALSE]
        }
        message("[", i, "/", nrow(todo), "] ", row$CellType, " ", row$gene_id,
                " x ", row$GWAS_RSID, " (chr", row$GWAS_CHR, ":", row$GWAS_BP, ")")
        r <- tryCatch(prepare_one_task(row, want),
                      error = function(e) manifest_row(row, "ERROR",
                                                       conditionMessage(e)))
        rows[[i]] <- r
        message("    ", r$status,
                if (identical(r$status, "prepared"))
                    paste0(": target=", r$target_SNP, " eQTL p=",
                           signif(r$target_eQTL_p, 3), " nSNP=",
                           r$n_overlap_harmonized) else "")
    }

    m <- do.call(rbind, rows)
    mf <- file.path(CFG$output_dir, paste0(
        "SMR_HEIDI.", CFG$analysis_group, ".", MODE$name, ".",
        if (length(CTs) == length(CFG$cell_types)) "all" else
            paste(CTs, collapse = "_"),
        ".manifest.tsv"))
    fwrite(m, mf, sep = "\t", quote = FALSE, na = "NA")
    message("\nManifest: ", mf)
    print(table(m$status, useNA = "ifany"))
    all_manifests[[mode_name]] <- m
}

if (length(all_manifests) > 1L) {
    combined <- do.call(rbind, all_manifests)
    cf <- file.path(CFG$output_dir, paste0(
        "SMR_HEIDI.", CFG$analysis_group, ".cis_and_trans.manifest.tsv"))
    fwrite(combined, cf, sep = "\t", quote = FALSE, na = "NA")
    message("\nCombined manifest: ", cf, "  (", nrow(combined), " tasks)")
}
