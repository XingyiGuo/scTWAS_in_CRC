#!/usr/bin/env Rscript

# Fast extraction of model SNPs for significant TWAS genes
# from either ADSER models or normal-colon (NL) models,
# followed by a one-pass GWAS RSID lookup across chromosomes 1-22.
#
# Set TF_NETWORK_GROUP=NL to extract the normal-colon models.
#
# Major speed improvements:
#   1) Read the significant-gene file once.
#   2) Open each required cell-type/model SQLite DB only once.
#   3) Extract all requested genes from each DB in one SQL query.
#   4) Combine unique model rsIDs across all cell types.
#   5) Scan each GWAS chromosome only once (22 scans total), rather than
#      repeatedly scanning the GWAS for each cell type.

suppressPackageStartupMessages({
  library(DBI)
  library(RSQLite)
  library(data.table)
})

# Respect the CPUs allocated by SLURM. Outside SLURM, auto-detect CPUs.
slurm_txt <- Sys.getenv("SLURM_CPUS_PER_TASK", unset = "")
if (nzchar(slurm_txt)) {
  available_cpus <- suppressWarnings(as.integer(slurm_txt))
} else {
  available_cpus <- suppressWarnings(parallel::detectCores(logical = FALSE))
}
if (is.na(available_cpus) || available_cpus < 1L) available_cpus <- 1L
if (.Platform$OS.type == "windows") available_cpus <- 1L

# GENCODE v26 gene annotation (GRCh38/hg38)
GTF_FILE <- Sys.getenv(
  "TF_NETWORK_GTF",
  "/data/l2_bioinfo1/liq17/LQ_references/gencode/gencode.v26.annotation.gtf"
)


# -------------------------------------------------------------------------
# Extract additional gene-level attributes directly from the GTF.
#
# GENCODE GTF files consistently provide gene_type/gene_name and strand.
# A free-text description is not guaranteed to be present in every GENCODE
# release. The code therefore checks both "gene_description" and
# "description"; if neither attribute is present, gene_description is NA.
# -------------------------------------------------------------------------

extract_gtf_attribute <- function(attribute_text, key) {

  pattern <- paste0(
    "(?:^|;[[:space:]]*)",
    key,
    "[[:space:]]+\"([^\"]+)\""
  )

  m <- regexec(
    pattern,
    attribute_text,
    perl = TRUE
  )

  hits <- regmatches(
    attribute_text,
    m
  )

  vapply(
    hits,
    function(z) {
      if (length(z) >= 2L) {
        z[2]
      } else {
        NA_character_
      }
    },
    character(1)
  )
}


read_gtf_gene_attributes <- function(gtf_file) {

  if (!file.exists(gtf_file)) {
    stop("GTF file not found: ", gtf_file)
  }

  # Read only gene records from the GTF.
  # grep is substantially faster and more memory-efficient than loading
  # every exon/transcript record and filtering afterward.
  cmd <- paste(
    "awk -F '\\t' '$3 == \"gene\" {print}'",
    shQuote(gtf_file)
  )

  gtf_gene <- fread(
    cmd = cmd,
    header = FALSE,
    sep = "\t",
    quote = "",
    fill = TRUE,
    showProgress = FALSE
  )

  if (ncol(gtf_gene) < 9L) {
    stop(
      "Unexpected GTF format: fewer than 9 columns were found in gene records."
    )
  }

  setnames(
    gtf_gene,
    1:9,
    c(
      "chr",
      "source",
      "feature",
      "start",
      "end",
      "score",
      "strand",
      "frame",
      "attributes"
    )
  )

  gtf_gene[, geneid := extract_gtf_attribute(
    as.character(attributes),
    "gene_id"
  )]

  gtf_gene[, GTF_genename := extract_gtf_attribute(
    as.character(attributes),
    "gene_name"
  )]

  gtf_gene[, gene_type := extract_gtf_attribute(
    as.character(attributes),
    "gene_type"
  )]

  # Some GTFs use gene_biotype instead of gene_type.
  gene_biotype <- extract_gtf_attribute(
    as.character(gtf_gene$attributes),
    "gene_biotype"
  )

  gtf_gene[
    is.na(gene_type) | gene_type == "",
    gene_type := gene_biotype[
      is.na(gene_type) | gene_type == ""
    ]
  ]

  gtf_gene[, gene_description := extract_gtf_attribute(
    as.character(attributes),
    "gene_description"
  )]

  # Also check for an attribute named simply "description".
  description_alt <- extract_gtf_attribute(
    as.character(gtf_gene$attributes),
    "description"
  )

  gtf_gene[
    is.na(gene_description) | gene_description == "",
    gene_description := description_alt[
      is.na(gene_description) | gene_description == ""
    ]
  ]

  # Remove ENSG version suffixes.
  gtf_gene[, geneid := sub(
    "\\.[0-9]+$",
    "",
    trimws(as.character(geneid))
  )]

  gtf_gene[
    ,
    .(
      geneid,
      GTF_genename,
      gene_chr = as.character(chr),
      gene_start = as.integer(start),
      gene_end = as.integer(end),
      gene_direction = as.character(strand),
      gene_type = as.character(gene_type),
      gene_description = as.character(gene_description)
    )
  ]
}


gtf_anno <- read_gtf_gene_attributes(GTF_FILE)

# Remove rows without a usable ENSG identifier.
gtf_anno <- gtf_anno[
  !is.na(geneid) &
    geneid != ""
]

# If duplicate ENSG IDs are present, retain the longest gene interval.
gtf_anno[, gene_length := gene_end - gene_start + 1L]

setorder(
  gtf_anno,
  geneid,
  -gene_length
)

gtf_anno <- gtf_anno[
  !duplicated(geneid)
]

gtf_anno[, gene_length := NULL]


# Report GTF annotation availability.
message(
  "GENCODE v26 gene records loaded: ",
  nrow(gtf_anno)
)

message(
  "Genes with gene_type: ",
  sum(
    !is.na(gtf_anno$gene_type) &
      gtf_anno$gene_type != ""
  )
)

message(
  "Genes with a free-text GTF description: ",
  sum(
    !is.na(gtf_anno$gene_description) &
      gtf_anno$gene_description != ""
  )
)

# -------------------------------------------------------------------------
# 1. Fixed settings
# -------------------------------------------------------------------------

test <- toupper(Sys.getenv("TF_NETWORK_GROUP", "ADSER"))

if (test == "ADSER") {

  CELL_TYPES <- c("ABS", "ASC", "CT", "EE", "GOB", "SSC", "STM", "TAC", "TUF")
  DATASET_PREFIX <- "VUMC_HTAN_ADSER"

  SIG_FILE <- paste0(
    "/data/l2_bioinfo1/liq17/scTWAS/",
    "TWAS_HTAN_GSE266546_WellcomeSanger_res/",
    "SPrediXcan_peakSNP_association/",
    "VUMC_HTAN_ADSER_EUR_Multi_sig/",
    "VUMC_HTAN_ADSER_9CTs_peakSNP.EUR_Multi.sig.csv"
  )

  DB_PATTERN <- paste0(
    "/data/l2_bioinfo1/liq17/scTWAS/",
    "TWAS_HTAN_GSE266546_WellcomeSanger_res/",
    "VUMC_HTAN_ADSER_%s_peakSNP/",
    "VUMC_HTAN_ADSER_%s_peakSNP.db"
  )

} else if (test == "NL") {

  CELL_TYPES <- c("ABS", "CT", "EE", "GOB", "STM", "TAC", "TUF")
  DATASET_PREFIX <- "HTAN_NL_WellcomeSanger_NL"

  SIG_FILE <- paste0(
    "/data/l2_bioinfo1/liq17/scTWAS/",
    "TWAS_HTAN_GSE266546_WellcomeSanger_res/",
    "SPrediXcan_peakSNP_association/",
    "Cauchy_Meta_NL_WellcomeSanger_NL_EUR_Multi/",
    "HTAN_NL_WellcomeSanger_NL_7celltypes_peakSNP.EUR_Multi.CauchyMeta.sig.csv"
  )

  VUMC_HTAN_GSE266546_PATTERN <- paste0(
    "/data/l2_bioinfo1/liq17/scTWAS/",
    "TWAS_HTAN_GSE266546_WellcomeSanger_res/",
    "VUMC_HTAN_GSE266546_NL_%s_peakSNP/",
    "VUMC_HTAN_GSE266546_NL_%s_peakSNP.db"
  )

  WellcomeSanger_PATTERN <- paste0(
    "/data/l2_bioinfo1/liq17/scTWAS/",
    "TWAS_HTAN_GSE266546_WellcomeSanger_res/",
    "WellcomeSanger_NL_%s_peakSNP/",
    "WellcomeSanger_NL_%s_peakSNP.db"
  )

} else {
  stop("Unknown test value: ", test, ". Use 'NL' or 'ADSER'.")
}

# DB and GWAS phases run sequentially, so each phase can use all allocated CPUs.
N_DB_CORES <- min(length(CELL_TYPES), available_cpus)
N_GWAS_CORES <- min(22L, available_cpus)

# Keep outputs next to the corresponding significant-gene input.
OUT_DIR <- file.path(dirname(SIG_FILE), "CRC51TFs_occupy_model_snps")

GWAS_PATTERN <- paste0(
  "/data/l2_bioinfo1/GWAS_SS/BIOVU/CRC/metal/GC_Final/processed/EUR/",
  "rsid_annotated_CHRPOS/",
  "final_gc_maf0.01_CRC_EUR_META_1.chr%d.chrpos.rsid.assoc.dosage.gz"
)

GWAS_COLS <- c(
  "SNP", "SNP_CHR_POS", "RSID", "A1", "A2", "FRQ", "BETA", "SE", "P"
)

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
BY_CT_DIR <- file.path(OUT_DIR, "by_celltype")
dir.create(BY_CT_DIR, recursive = TRUE, showWarnings = FALSE)

# -------------------------------------------------------------------------
# 2. Read significant genes once
# -------------------------------------------------------------------------

if (!file.exists(SIG_FILE)) {
  stop("Significant-gene file not found: ", SIG_FILE)
}

sig <- fread(SIG_FILE, showProgress = FALSE)

required_cols <- c("gene", "CellType")
if (test == "NL") {
  required_cols <- c(required_cols, "P driven_Multi")
}

missing_cols <- setdiff(required_cols, names(sig))
if (length(missing_cols) > 0L) {
  stop(
    "Missing required column(s) in SIG_FILE: ",
    paste(missing_cols, collapse = ", ")
  )
}

sig[, CellType := toupper(trimws(as.character(CellType)))]
sig[, gene := trimws(as.character(gene))]

if (test == "ADSER") {

  sig <- unique(sig[
    CellType %chin% CELL_TYPES &
      !is.na(gene) &
      gene != "",
    .(CellType, gene)
  ])

} else if (test == "NL") {

  # P driven_Multi indicates which normal-colon prediction model
  # should supply the SNP weights for each significant gene.
  sig[, Model_source := trimws(as.character(`P driven_Multi`))]

  valid_sources <- c(
    "VUMC_HTAN_GSE266546",
    "WellcomeSanger"
  )

  invalid_sources <- unique(
    sig[
      !is.na(Model_source) &
        Model_source != "" &
        !Model_source %chin% valid_sources,
      Model_source
    ]
  )

  if (length(invalid_sources) > 0L) {
    stop(
      "Unexpected value(s) in 'P driven_Multi': ",
      paste(invalid_sources, collapse = ", ")
    )
  }

  sig <- unique(sig[
    CellType %chin% CELL_TYPES &
      !is.na(gene) &
      gene != "" &
      Model_source %chin% valid_sources,
    .(CellType, gene, Model_source)
  ])
}

message("Significant CellType-gene pairs: ", uniqueN(sig[, .(CellType, gene)]))

if (test == "NL") {
  message("NL genes by model source:")
  print(sig[, .N, by = .(CellType, Model_source)][order(CellType, Model_source)])
}

# -------------------------------------------------------------------------
# 3. Extract model SNPs from one cell-type DB in a single SQL query
# -------------------------------------------------------------------------

extract_one_ct <- function(ct, genes_ct, db_file) {

  if (!file.exists(db_file)) {
    warning("DB file not found for ", ct, ": ", db_file)
    return(list(
      snps = data.table(),
      qc = data.table(
        CellType = ct,
        gene = genes_ct,
        n_snps_extracted = 0L,
        present_in_modelDB = FALSE
      ),
      missing = data.table(
        CellType = ct,
        gene = genes_ct,
        reason = "DB file not found"
      )
    ))
  }

  con <- dbConnect(SQLite(), db_file)
  on.exit(dbDisconnect(con), add = TRUE)

  needed_tables <- c("weights", "extra")
  db_tables <- dbListTables(con)
  if (!all(needed_tables %in% db_tables)) {
    stop(
      "Required table(s) missing in ", db_file, ": ",
      paste(setdiff(needed_tables, db_tables), collapse = ", ")
    )
  }

  wanted <- data.frame(gene = unique(genes_ct), stringsAsFactors = FALSE)
  dbWriteTable(
    con,
    "wanted_genes_tmp",
    wanted,
    temporary = TRUE,
    overwrite = TRUE
  )

  # One SQL query per cell type. This retrieves SNP-level weights and all gene/model
  # information from extra without repeatedly querying individual genes.
  sql <- paste0(
    "SELECT ",
    "w.rsid, w.gene, w.weight, w.ref_allele, w.eff_allele, ",
    "e.genename, e.[pred.perf.r2], e.[pred.perf.pval], e.[pred.perf.qval], ",
    "e.[n.snps.in.model], e.R2_global, e.r_squared_global, e.[r.global] ",
    "FROM weights AS w ",
    "INNER JOIN wanted_genes_tmp AS g ON w.gene = g.gene ",
    "LEFT JOIN extra AS e ON w.gene = e.gene"
  )

  snps <- as.data.table(dbGetQuery(con, sql))

  if (nrow(snps) > 0L) {
    setnames(snps, "rsid", "rs")
    setnames(snps, "n.snps.in.model", "n_snps_in_model")
    snps[, CellType := ct]

    preferred <- c(
      "CellType", "gene", "genename", "n_snps_in_model",
      "rs", "weight", "ref_allele", "eff_allele"
    )
    setcolorder(snps, c(preferred, setdiff(names(snps), preferred)))
  }

  # Gene-level extraction QC.
  counts <- if (nrow(snps) > 0L) {
    snps[, .(
      n_snps_extracted = .N,
      n_snps_in_model = suppressWarnings(as.numeric(first(n_snps_in_model))),
      genename = first(genename)
    ), by = gene]
  } else {
    data.table(
      gene = character(),
      n_snps_extracted = integer(),
      n_snps_in_model = numeric(),
      genename = character()
    )
  }

  qc <- merge(
    data.table(CellType = ct, gene = unique(genes_ct)),
    counts,
    by = "gene",
    all.x = TRUE,
    sort = FALSE
  )

  qc[is.na(n_snps_extracted), n_snps_extracted := 0L]
  qc[, present_in_modelDB := n_snps_extracted > 0L]
  qc[, count_matches_db :=
       !is.na(n_snps_in_model) & n_snps_in_model == n_snps_extracted]

  missing <- qc[present_in_modelDB == FALSE, .(
    CellType,
    gene,
    reason = "Gene has no SNP rows in model DB"
  )]

  list(snps = snps, qc = qc, missing = missing)
}

# -------------------------------------------------------------------------
# 4. Extract all cell types
# -------------------------------------------------------------------------

message("\nAvailable CPUs: ", available_cpus)
message("Parallel DB extraction workers: ", N_DB_CORES)

if (N_DB_CORES == 1L && nzchar(slurm_txt)) {
  message(
    paste0(
      "NOTE: SLURM_CPUS_PER_TASK=1, so DB extraction is serial. ",
      "Request more CPUs (e.g. --cpus-per-task=", length(CELL_TYPES),
      ") to parallelize all ", length(CELL_TYPES), " cell types."
    )
  )
}

extract_results <- parallel::mclapply(
  CELL_TYPES,
  function(ct) {

    sig_ct <- sig[CellType == ct]

    message(
      "Extracting ", ct,
      " | genes = ", uniqueN(sig_ct$gene)
    )

    if (nrow(sig_ct) == 0L) {
      return(list(
        snps = data.table(),
        qc = data.table(),
        missing = data.table()
      ))
    }

    if (test == "ADSER") {

      db_file <- sprintf(DB_PATTERN, ct, ct)
      return(extract_one_ct(ct, sig_ct$gene, db_file))

    } else if (test == "NL") {

      source_patterns <- list(
        VUMC_HTAN_GSE266546 = VUMC_HTAN_GSE266546_PATTERN,
        WellcomeSanger = WellcomeSanger_PATTERN
      )

      source_results <- lapply(
        names(source_patterns),
        function(model_source) {

          genes_source <- sig_ct[
            Model_source == model_source,
            gene
          ]

          if (length(genes_source) == 0L) {
            return(NULL)
          }

          db_file <- sprintf(
            source_patterns[[model_source]],
            ct,
            ct
          )

          message(
            "  ", ct,
            " | model source = ", model_source,
            " | genes = ", length(genes_source)
          )

          res <- extract_one_ct(
            ct = ct,
            genes_ct = genes_source,
            db_file = db_file
          )

          if (nrow(res$snps) > 0L) {
            res$snps[, Model_source := model_source]
          }
          if (nrow(res$qc) > 0L) {
            res$qc[, Model_source := model_source]
          }
          if (nrow(res$missing) > 0L) {
            res$missing[, Model_source := model_source]
          }

          res
        }
      )

      source_results <- Filter(Negate(is.null), source_results)

      if (length(source_results) == 0L) {
        return(list(
          snps = data.table(),
          qc = data.table(),
          missing = data.table()
        ))
      }

      return(list(
        snps = rbindlist(
          lapply(source_results, `[[`, "snps"),
          use.names = TRUE,
          fill = TRUE
        ),
        qc = rbindlist(
          lapply(source_results, `[[`, "qc"),
          use.names = TRUE,
          fill = TRUE
        ),
        missing = rbindlist(
          lapply(source_results, `[[`, "missing"),
          use.names = TRUE,
          fill = TRUE
        )
      ))
    }
  },
  mc.cores = N_DB_CORES,
  mc.preschedule = TRUE
)

names(extract_results) <- CELL_TYPES

model_list <- lapply(extract_results, `[[`, "snps")
qc_list <- lapply(extract_results, `[[`, "qc")
missing_list <- lapply(extract_results, `[[`, "missing")

model_snps <- rbindlist(model_list, use.names = TRUE, fill = TRUE)
model_qc <- rbindlist(qc_list, use.names = TRUE, fill = TRUE)
missing_genes <- rbindlist(missing_list, use.names = TRUE, fill = TRUE)

if (nrow(model_snps) == 0L) {
  stop("No model SNPs were extracted from any of the ", length(CELL_TYPES) ," DB files.")
}

message("\nTotal model SNP rows across ",length(CELL_TYPES) , " CTs: ", nrow(model_snps))
message("Unique rsIDs across ", length(CELL_TYPES) ," CTs: ", uniqueN(model_snps$rs))
message("Unique CellType-gene pairs represented: ",
        uniqueN(model_snps[, .(CellType, gene)]))

# -------------------------------------------------------------------------
# 5. Fast GWAS lookup: scan chromosomes 1-22 only once for all cell types
# -------------------------------------------------------------------------

lookup_gwas_once <- function(model_rs, gwas_pattern, chromosomes = 1:22, n_cores = 1L) {

  model_rs <- unique(trimws(as.character(model_rs)))
  model_rs <- model_rs[!is.na(model_rs) & model_rs != ""]

  message("\nUnique model rsIDs to check in GWAS: ", length(model_rs))

  if (length(model_rs) == 0L) {
    return(data.table())
  }

  if (!nzchar(Sys.which("awk")) || !nzchar(Sys.which("gzip"))) {
    stop("awk and gzip are required for the fast GWAS lookup but were not found in PATH.")
  }

  rs_file <- tempfile(pattern = "allCT_model_rsids_", fileext = ".txt")
  writeLines(model_rs, rs_file)
  on.exit(unlink(rs_file), add = TRUE)

  message("GWAS chromosome scan workers: ", n_cores)

  hits <- parallel::mclapply(
    chromosomes,
    function(chr) {
      gwas_file <- sprintf(gwas_pattern, chr)

      if (!file.exists(gwas_file)) {
        warning("GWAS file missing for chr", chr, ": ", gwas_file)
        return(data.table())
      }

      message("Scanning GWAS chr", chr, " ...")

      # Exact RSID lookup against GWAS column 3.  awk first loads the small
      # model-rsID set into memory, then streams the compressed GWAS file once.
      # This is faster and safer than grepping every GWAS field.
      awk_prog <- "NR==FNR { want[$1]=1; next } ($3 in want)"
      cmd <- sprintf(
        "LC_ALL=C gzip -cd %s 2>/dev/null | awk %s %s -",
        shQuote(gwas_file),
        shQuote(awk_prog),
        shQuote(rs_file)
      )

      x <- tryCatch(
        fread(
          cmd = cmd,
          header = FALSE,
          fill = TRUE,
          showProgress = FALSE
        ),
        error = function(e) {
          warning("Could not parse GWAS chr", chr, ": ", conditionMessage(e))
          data.table()
        }
      )

      if (nrow(x) == 0L) return(data.table())

      if (ncol(x) != length(GWAS_COLS)) {
        stop(
          "Unexpected GWAS column count on chr", chr,
          ": expected ", length(GWAS_COLS), ", found ", ncol(x)
        )
      }

      setnames(x, GWAS_COLS)
      x[, RSID := trimws(as.character(RSID))]
      x <- x[RSID %chin% model_rs]

      if (nrow(x) > 0L) {
        x[, GWAS_chr := as.integer(chr)]
        message("chr", chr, " matched rows: ", nrow(x))
      }
      x
    },
    mc.cores = n_cores,
    mc.preschedule = TRUE
  )

  hits <- rbindlist(hits, use.names = TRUE, fill = TRUE)
  if (nrow(hits) == 0L) return(hits)

  unique(hits)
}

gwas_hits <- lookup_gwas_once(model_snps$rs, GWAS_PATTERN, n_cores = N_GWAS_CORES)

# -------------------------------------------------------------------------
# 6. Join GWAS information back to all model SNP rows
# -------------------------------------------------------------------------

model_snps[, .row_order := .I]

if (nrow(gwas_hits) > 0L) {

  # Count duplicates per RSID; keep one representative GWAS row per RSID so
  # model SNP rows are not duplicated unexpectedly.
  gwas_n <- gwas_hits[, .(GWAS_n_rows = .N), by = RSID]
  gwas_one <- gwas_hits[, .SD[1L], by = RSID]
  gwas_one <- merge(gwas_one, gwas_n, by = "RSID", all.x = TRUE, sort = FALSE)

  setnames(
    gwas_one,
    old = c("SNP", "SNP_CHR_POS", "RSID", "A1", "A2", "FRQ", "BETA", "SE", "P"),
    new = c(
      "GWAS_SNP", "GWAS_SNP_CHR_POS", "rs",
      "GWAS_A1", "GWAS_A2", "GWAS_FRQ",
      "GWAS_BETA", "GWAS_SE", "GWAS_P"
    )
  )

  model_snps <- merge(
    model_snps,
    gwas_one,
    by = "rs",
    all.x = TRUE,
    sort = FALSE
  )

  model_snps[, GWAS_present := !is.na(GWAS_chr)]
  model_snps[is.na(GWAS_n_rows), GWAS_n_rows := 0L]

  # Direct allele-pair consistency check; order does not matter.
  model_snps[, GWAS_alleles_match := fifelse(
    !GWAS_present,
    NA,
    (toupper(ref_allele) == toupper(GWAS_A1) &
       toupper(eff_allele) == toupper(GWAS_A2)) |
      (toupper(ref_allele) == toupper(GWAS_A2) &
         toupper(eff_allele) == toupper(GWAS_A1))
  )]

} else {

  model_snps[, `:=`(
    GWAS_chr = NA_integer_,
    GWAS_n_rows = 0L,
    GWAS_SNP = NA_character_,
    GWAS_SNP_CHR_POS = NA_character_,
    GWAS_A1 = NA_character_,
    GWAS_A2 = NA_character_,
    GWAS_FRQ = NA_character_,
    GWAS_BETA = NA_real_,
    GWAS_SE = NA_real_,
    GWAS_P = NA_real_,
    GWAS_present = FALSE,
    GWAS_alleles_match = NA
  )]
}

setorder(model_snps, .row_order)
model_snps[, .row_order := NULL]

# Put core fields first.
preferred_cols <- c(
  "CellType", "gene", "Model_source", "genename", "n_snps_in_model",
  "rs", "weight", "ref_allele", "eff_allele",
  "GWAS_present", "GWAS_alleles_match", "GWAS_chr", "GWAS_n_rows",
  "GWAS_SNP", "GWAS_SNP_CHR_POS",
  "GWAS_A1", "GWAS_A2", "GWAS_FRQ", "GWAS_BETA", "GWAS_SE", "GWAS_P"
)
preferred_cols <- preferred_cols[preferred_cols %in% names(model_snps)]
setcolorder(model_snps, c(preferred_cols, setdiff(names(model_snps), preferred_cols)))

# -------------------------------------------------------------------------
# 7. Add GWAS coverage to gene-level QC
# -------------------------------------------------------------------------

qc_keys <- c("CellType", "gene")
if ("Model_source" %in% names(model_snps)) {
  qc_keys <- c(qc_keys, "Model_source")
}

gwas_qc <- model_snps[, .(
  n_model_snps_checked = .N,
  n_model_snps_in_GWAS = sum(GWAS_present, na.rm = TRUE),
  pct_model_snps_in_GWAS = round(
    100 * mean(GWAS_present, na.rm = TRUE),
    2
  ),
  all_model_snps_in_GWAS = all(GWAS_present)
), by = qc_keys]

model_qc <- merge(
  model_qc,
  gwas_qc,
  by = qc_keys,
  all.x = TRUE,
  sort = FALSE
)

# -------------------------------------------------------------------------
# 8. Prepare and write outputs: plain CSV, not gzip
# -------------------------------------------------------------------------

COMBINED_PREFIX <- paste0(
  DATASET_PREFIX,
  "_",
  length(CELL_TYPES),
  "CTs"
)

COMBINED_MODEL_FILE <- file.path(
  OUT_DIR,
  paste0(COMBINED_PREFIX, "_sigGenes_modelSNPs.csv")
)

QC_FILE <- file.path(
  OUT_DIR,
  paste0(COMBINED_PREFIX, "_sigGenes_modelSNP_QC.csv")
)

MISSING_FILE <- file.path(
  OUT_DIR,
  paste0(COMBINED_PREFIX, "_sigGenes_missing_from_modelDB.csv")
)

# -------------------------------------------------------------------------
# 8a. Add GENCODE v26 gene annotations
#
# The model 'gene' column contains ENSG IDs and is used as the primary
# annotation key. Existing model-DB gene symbols are retained; if a gene
# symbol is missing, it is filled from GENCODE.
#
# Added GENCODE fields:
#   genename         = gene symbol
#   gene_chr         = chromosome
#   gene_start       = gene start coordinate
#   gene_end         = gene end coordinate
#   gene_direction   = strand (+/-)
#   gene_type        = GENCODE gene type/biotype
#   gene_description = free-text GTF description, when available
# -------------------------------------------------------------------------

annotate_gene_from_gtf <- function(x, gtf_annotation) {

  x <- as.data.table(copy(x))

  # Ensure annotation columns exist even for an empty result table.
  if (nrow(x) == 0L) {

    if (!"genename" %in% names(x)) {
      x[, genename := character()]
    }

    x[, `:=`(
      gene_chr = character(),
      gene_start = integer(),
      gene_end = integer(),
      gene_direction = character(),
      gene_type = character(),
      gene_description = character()
    )]

    return(x)
  }

  if (!"gene" %in% names(x)) {
    stop(
      "Cannot add GTF annotation because the 'gene' column is missing."
    )
  }

  # Preserve the original row order.
  x[, .annotation_row_order := .I]

  # Normalize ENSG IDs by removing optional version suffixes.
  x[, .gene_gtf_key := sub(
    "\\.[0-9]+$",
    "",
    trimws(as.character(gene))
  )]

  # Join GENCODE annotation by ENSG ID.
  x[
    gtf_annotation,
    on = .(.gene_gtf_key = geneid),
    `:=`(
      .gtf_genename = i.GTF_genename,
      gene_chr = i.gene_chr,
      gene_start = i.gene_start,
      gene_end = i.gene_end,
      gene_direction = i.gene_direction,
      gene_type = i.gene_type,
      gene_description = i.gene_description
    )
  ]

  # Preserve the model-DB gene symbol when available;
  # otherwise fill it from the GENCODE annotation.
  if ("genename" %in% names(x)) {

    x[
      is.na(genename) |
        trimws(as.character(genename)) == "",
      genename := .gtf_genename
    ]

  } else {

    x[, genename := .gtf_genename]
  }

  setorder(
    x,
    .annotation_row_order
  )

  x[, c(
    ".annotation_row_order",
    ".gene_gtf_key",
    ".gtf_genename"
  ) := NULL]

  x
}


model_snps <- annotate_gene_from_gtf(
  model_snps,
  gtf_anno
)

model_qc <- annotate_gene_from_gtf(
  model_qc,
  gtf_anno
)

missing_genes <- annotate_gene_from_gtf(
  missing_genes,
  gtf_anno
)


# Place gene annotation fields near the ENSG identifier.
reorder_gene_annotation <- function(x) {

  x <- as.data.table(x)

  first_cols <- c(
    "CellType",
    "gene",
    "genename",
    "gene_chr",
    "gene_start",
    "gene_end",
    "gene_direction",
    "gene_type",
    "gene_description",
    "Model_source"
  )

  first_cols <- first_cols[
    first_cols %in% names(x)
  ]

  setcolorder(
    x,
    c(
      first_cols,
      setdiff(
        names(x),
        first_cols
      )
    )
  )

  x
}


model_snps <- reorder_gene_annotation(
  model_snps
)

model_qc <- reorder_gene_annotation(
  model_qc
)

missing_genes <- reorder_gene_annotation(
  missing_genes
)


# Annotation QC
message(
  "Model-SNP rows with GTF coordinates: ",
  sum(!is.na(model_snps$gene_chr)),
  " / ",
  nrow(model_snps)
)

message(
  "Model-SNP rows with GTF gene type: ",
  sum(
    !is.na(model_snps$gene_type) &
      model_snps$gene_type != ""
  ),
  " / ",
  nrow(model_snps)
)

message(
  "Model-SNP rows with GTF description: ",
  sum(
    !is.na(model_snps$gene_description) &
      model_snps$gene_description != ""
  ),
  " / ",
  nrow(model_snps)
)

message(
  "Gene-level QC rows with GTF coordinates: ",
  sum(!is.na(model_qc$gene_chr)),
  " / ",
  nrow(model_qc)
)

if (nrow(missing_genes) > 0L) {

  message(
    "Missing-model genes with GTF coordinates: ",
    sum(!is.na(missing_genes$gene_chr)),
    " / ",
    nrow(missing_genes)
  )
}



fwrite(model_snps, COMBINED_MODEL_FILE)
fwrite(model_qc, QC_FILE)
fwrite(missing_genes, MISSING_FILE)

# Also write one plain CSV per cell type without repeating the GWAS scans.
for (ct in CELL_TYPES) {

  out_ct <- file.path(
    BY_CT_DIR,
    paste0(
      DATASET_PREFIX,
      "_",
      ct,
      "_sigGenes_modelSNPs.csv"
    )
  )

  fwrite(
    model_snps[CellType == ct],
    out_ct
  )
}

message("\n============================================================")
message("DONE")
message("============================================================")
message("Analysis type: ", test)
message("Combined model SNP file: ", COMBINED_MODEL_FILE)
message("QC file: ", QC_FILE)
message("Missing-gene file: ", MISSING_FILE)
message("Per-cell-type output directory: ", BY_CT_DIR)
message("Total rows: ", nrow(model_snps))
message("Unique rsIDs: ", uniqueN(model_snps$rs))
message(
  "GWAS-present model SNP rows: ",
  sum(model_snps$GWAS_present, na.rm = TRUE)
)
message(
  "GWAS-missing model SNP rows: ",
  sum(!model_snps$GWAS_present, na.rm = TRUE)
)
