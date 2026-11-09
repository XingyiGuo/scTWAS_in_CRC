#!/usr/bin/env bash
set -euo pipefail

# ==============================================================================
# Run SMR + HEIDI from the manifest produced by:
#   GWAS_leadvar_SMR_HEIDI_inputs.R
#
# Adapted from 02b_run_SMR_HEIDI.sh: the region is the GWAS lead variant
# +/-500 kb, the BESD probe is annotated at the lead variant, and --cis-wind is
# taken from the manifest instead of being hardcoded.
#
# Usage:
#   bash GWAS_leadvar_run_SMR_HEIDI.sh \
#     /data/l2_bioinfo1/liq17/scTWAS/TWAS_HTAN_GSE266546_WellcomeSanger_res/SMR_HEIDI_GWAS_eQTL_sigSource/NL_targetTopEQTL_500kb/SMR_HEIDI.NL.ABS.manifest.tsv
#
# Optional executable overrides:
#   SMR_BIN=/path/to/smr PLINK2_BIN=/path/to/plink2 bash ...
#
# Per task this script:
#   1) converts .esd/.flist -> .besd/.esi/.epi
#   2) extracts the harmonized regional SNPs from the source-matched PGEN
#      genotype data and creates a task-specific PLINK bed/bim/fam LD reference
#   3) runs SMR + HEIDI with the PRESELECTED target SNP using --target-snp
#   4) writes a run-status table and a combined SMR result table
# ==============================================================================

if [[ $# -lt 1 ]]; then
    echo "Usage: bash $0 <SMR_HEIDI.manifest.tsv>" >&2
    exit 1
fi

MANIFEST=$1
SMR_BIN=${SMR_BIN:-/data/l2_bioinfo1/liq17/biosoft/smr-1.4.2-linux-x86_64/smr}
PLINK2_BIN=${PLINK2_BIN:-plink2}

if [[ ! -f "$MANIFEST" ]]; then
    echo "ERROR: manifest not found: $MANIFEST" >&2
    exit 1
fi

if ! command -v "$SMR_BIN" >/dev/null 2>&1 && [[ ! -x "$SMR_BIN" ]]; then
    echo "ERROR: SMR executable not found: $SMR_BIN" >&2
    exit 1
fi

if ! command -v "$PLINK2_BIN" >/dev/null 2>&1 && [[ ! -x "$PLINK2_BIN" ]]; then
    echo "ERROR: PLINK2 executable not found: $PLINK2_BIN" >&2
    exit 1
fi

OUTROOT=$(dirname "$MANIFEST")
RUN_STATUS="${MANIFEST%.manifest.tsv}.run_status.tsv"
COMBINED="${MANIFEST%.manifest.tsv}.SMR_HEIDI.results.tsv"

printf 'analysis_group\tcell_type\tgene\tcohort\teqtl_source\ttarget_SNP\trun_status\tn_ld_snps\tsmr_file\tlog_file\tmessage\n' > "$RUN_STATUS"

# Pull only runner-required columns by header name, so the script is robust to
# future additions/reordering of columns in the manifest.
TASK_STREAM=$(mktemp)
trap 'rm -f "$TASK_STREAM"' EXIT

awk -F'\t' 'BEGIN{OFS="\t"}
NR==1 {
  for(i=1;i<=NF;i++) h[$i]=i;
  required[1]="analysis_group"; required[2]="cell_type"; required[3]="gene";
  required[4]="cohort"; required[5]="eqtl_source"; required[6]="target_SNP";
  required[7]="gwas_ma"; required[8]="eqtl_flist"; required[9]="beqtl_prefix";
  required[10]="snplist"; required[11]="pfile_prefix"; required[12]="task_out_prefix";
  required[13]="peqtl_smr"; required[14]="peqtl_heidi"; required[15]="heidi_mtd";
  required[16]="heidi_min_m"; required[17]="heidi_max_m";
  required[18]="ld_lower_limit"; required[19]="ld_upper_limit";
  required[20]="maf"; required[21]="diff_freq"; required[22]="diff_freq_prop";
  required[23]="thread_num"; required[24]="status";
  required[25]="cis_wind";
  for(j=1;j<=25;j++) if(!(required[j] in h)) {
    print "ERROR: manifest missing column " required[j] > "/dev/stderr"; exit 2;
  }
  next
}
{
  print $(h["analysis_group"]), $(h["cell_type"]), $(h["gene"]),
        $(h["cohort"]), $(h["eqtl_source"]), $(h["target_SNP"]),
        $(h["gwas_ma"]), $(h["eqtl_flist"]), $(h["beqtl_prefix"]),
        $(h["snplist"]), $(h["pfile_prefix"]), $(h["task_out_prefix"]),
        $(h["peqtl_smr"]), $(h["peqtl_heidi"]), $(h["heidi_mtd"]),
        $(h["heidi_min_m"]), $(h["heidi_max_m"]),
        $(h["ld_lower_limit"]), $(h["ld_upper_limit"]),
        $(h["maf"]), $(h["diff_freq"]), $(h["diff_freq_prop"]),
        $(h["thread_num"]), $(h["status"]), $(h["cis_wind"])
}' "$MANIFEST" > "$TASK_STREAM"

n_total=$(wc -l < "$TASK_STREAM" | awk '{print $1}')
echo "SMR+HEIDI tasks in manifest: $n_total"
echo "Run status: $RUN_STATUS"

while IFS=$'\t' read -r \
    analysis_group cell_type gene cohort eqtl_source target_snp \
    gwas_ma eqtl_flist beqtl_prefix snplist pfile_prefix task_out_prefix \
    peqtl_smr peqtl_heidi heidi_mtd heidi_min_m heidi_max_m \
    ld_lower_limit ld_upper_limit maf diff_freq diff_freq_prop thread_num prep_status \
    cis_wind

do
    echo
    echo "======================================================================"
    echo "SMR+HEIDI: ${analysis_group} ${cell_type} ${gene} ${cohort} ${eqtl_source}"
    echo "Target SNP: ${target_snp}"
    echo "======================================================================"

    log_file="${task_out_prefix}.run.log"
    smr_file="${task_out_prefix}.smr"
    task_dir=$(dirname "$task_out_prefix")
    mkdir -p "$task_dir"

    # Skip tasks that were not fully prepared.
    if [[ "$prep_status" != "prepared" ]]; then
        msg="prep_status=${prep_status}"
        printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
          "$analysis_group" "$cell_type" "$gene" "$cohort" "$eqtl_source" \
          "$target_snp" "SKIPPED_PREP" "NA" "$smr_file" "$log_file" "$msg" >> "$RUN_STATUS"
        echo "SKIP: $msg"
        continue
    fi

    # Basic input checks.
    missing=0
    for f in "$gwas_ma" "$eqtl_flist" "$snplist"; do
        if [[ ! -f "$f" ]]; then
            echo "Missing input: $f" >&2
            missing=1
        fi
    done
    for ext in pgen pvar psam; do
        if [[ ! -f "${pfile_prefix}.${ext}" ]]; then
            echo "Missing LD PGEN input: ${pfile_prefix}.${ext}" >&2
            missing=1
        fi
    done
    if [[ $missing -ne 0 ]]; then
        printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
          "$analysis_group" "$cell_type" "$gene" "$cohort" "$eqtl_source" \
          "$target_snp" "MISSING_INPUT" "NA" "$smr_file" "$log_file" \
          "one_or_more_required_inputs_missing" >> "$RUN_STATUS"
        continue
    fi

    : > "$log_file"

    # -------------------------------------------------------------------------
    # 1. Convert the single-gene ESD to SMR BESD format.
    # -------------------------------------------------------------------------
    if [[ ! -f "${beqtl_prefix}.besd" || ! -f "${beqtl_prefix}.esi" || ! -f "${beqtl_prefix}.epi" ]]; then
        echo "[make-besd] $beqtl_prefix" | tee -a "$log_file"
        if ! "$SMR_BIN" \
            --eqtl-flist "$eqtl_flist" \
            --make-besd \
            --out "$beqtl_prefix" \
            >> "$log_file" 2>&1; then
            printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
              "$analysis_group" "$cell_type" "$gene" "$cohort" "$eqtl_source" \
              "$target_snp" "MAKE_BESD_FAILED" "NA" "$smr_file" "$log_file" \
              "smr_make_besd_failed" >> "$RUN_STATUS"
            continue
        fi
    else
        echo "[make-besd] existing BESD found; reuse" | tee -a "$log_file"
    fi

    # -------------------------------------------------------------------------
    # 2. Build a local LD reference from the same genotype source used for eQTL.
    #    This also guarantees that only SNPs used by GWAS/eQTL harmonization are
    #    passed to HEIDI.
    # -------------------------------------------------------------------------
    ld_prefix="${task_out_prefix}.LDref"
    if [[ ! -s "${ld_prefix}.bed" || ! -s "${ld_prefix}.bim" || ! -s "${ld_prefix}.fam" ]]; then
        echo "[LD] preparing harmonized SNPs from $pfile_prefix" | tee -a "$log_file"

        # TensorQTL/PGEN IDs can be composite, e.g.
        #   rs75258813;14_...._REF_ALT
        # while the SMR/GWAS/eQTL files and snplist use the rsID only.
        # Build a PLINK2 --update-name map so the LD reference also uses rsIDs.
        # PLINK2 applies --update-name before --extract.
        pvar_file="${pfile_prefix}.pvar"
        rename_map="${ld_prefix}.update_name.tsv"

        awk -F'\t' 'BEGIN{OFS="\t"}
          NR==FNR {
            x=tolower($1); gsub(/^[[:space:]]+|[[:space:]]+$/, "", x);
            if (x != "") wanted[x]=1;
            next
          }
          /^#/ {next}
          {
            id=$3;
            rs=id;
            sub(/;.*/, "", rs);
            rs=tolower(rs);
            if ((rs in wanted) && id != rs) print id, rs;
          }
        ' "$snplist" "$pvar_file" > "$rename_map"

        n_requested=$(awk 'NF{n++} END{print n+0}' "$snplist")
        n_pvar_match=$(awk -F'\t' '
          NR==FNR {
            x=tolower($1); gsub(/^[[:space:]]+|[[:space:]]+$/, "", x);
            if (x != "") wanted[x]=1;
            next
          }
          /^#/ {next}
          {
            rs=$3;
            sub(/;.*/, "", rs);
            rs=tolower(rs);
            if ((rs in wanted) && !seen[rs]++) n++;
          }
          END{print n+0}
        ' "$snplist" "$pvar_file")

        n_dup_rsid=$(awk -F'\t' '
          NR==FNR {
            x=tolower($1); gsub(/^[[:space:]]+|[[:space:]]+$/, "", x);
            if (x != "") wanted[x]=1;
            next
          }
          /^#/ {next}
          {
            rs=$3;
            sub(/;.*/, "", rs);
            rs=tolower(rs);
            if (rs in wanted) count[rs]++;
          }
          END{
            d=0;
            for (rs in count) if (count[rs] > 1) d++;
            print d+0;
          }
        ' "$snplist" "$pvar_file")

        echo "[LD] requested rsIDs: ${n_requested}; found in PVAR after rsID normalization: ${n_pvar_match}" | tee -a "$log_file"
        echo "[LD] PVAR IDs requiring rename to rsID: $(wc -l < "$rename_map" | awk '{print $1}')" | tee -a "$log_file"

        if [[ "$n_pvar_match" -eq 0 ]]; then
            printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
              "$analysis_group" "$cell_type" "$gene" "$cohort" "$eqtl_source" \
              "$target_snp" "LD_PREP_FAILED" "0" "$smr_file" "$log_file" \
              "no_harmonized_rsids_match_pvar_even_after_composite_id_normalization" >> "$RUN_STATUS"
            echo "ERROR: no requested rsIDs match ${pvar_file}, even after stripping the ';...' suffix." | tee -a "$log_file" >&2
            continue
        fi

        if [[ "$n_dup_rsid" -gt 0 ]]; then
            printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
              "$analysis_group" "$cell_type" "$gene" "$cohort" "$eqtl_source" \
              "$target_snp" "LD_PREP_FAILED" "NA" "$smr_file" "$log_file" \
              "duplicate_requested_rsids_in_pvar_after_normalization" >> "$RUN_STATUS"
            echo "ERROR: ${n_dup_rsid} requested rsID(s) map to multiple PVAR records; not safe to rename automatically." | tee -a "$log_file" >&2
            continue
        fi

        plink_args=(
            --pfile "$pfile_prefix"
        )
        if [[ -s "$rename_map" ]]; then
            plink_args+=(--update-name "$rename_map")
        fi
        plink_args+=(
            --extract "$snplist"
            --maf "$maf"
            --make-bed
            --out "$ld_prefix"
        )

        if ! "$PLINK2_BIN" "${plink_args[@]}" >> "$log_file" 2>&1; then
            printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
              "$analysis_group" "$cell_type" "$gene" "$cohort" "$eqtl_source" \
              "$target_snp" "LD_PREP_FAILED" "NA" "$smr_file" "$log_file" \
              "plink2_ld_reference_failed" >> "$RUN_STATUS"
            continue
        fi
    else
        echo "[LD] existing local bed/bim/fam found; reuse" | tee -a "$log_file"
    fi

    n_ld_snps=$(wc -l < "${ld_prefix}.bim" | awk '{print $1}')
    echo "[LD] SNPs retained: $n_ld_snps" | tee -a "$log_file"

    if [[ "$n_ld_snps" -lt "$heidi_min_m" ]]; then
        printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
          "$analysis_group" "$cell_type" "$gene" "$cohort" "$eqtl_source" \
          "$target_snp" "TOO_FEW_LD_SNPS" "$n_ld_snps" "$smr_file" "$log_file" \
          "LD SNP count < heidi_min_m" >> "$RUN_STATUS"
        continue
    fi

    if ! awk -v s="$target_snp" '$2==s{found=1} END{exit(found?0:1)}' "${ld_prefix}.bim"; then
        printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
          "$analysis_group" "$cell_type" "$gene" "$cohort" "$eqtl_source" \
          "$target_snp" "TARGET_NOT_IN_LD" "$n_ld_snps" "$smr_file" "$log_file" \
          "target SNP absent after PGEN->BED LD extraction" >> "$RUN_STATUS"
        continue
    fi

    # -------------------------------------------------------------------------
    # 3. Run SMR + HEIDI. --target-snp forces the user-defined top eQTL target.
    #    --cis-wind is in Kb and comes from the manifest; the probe sits at the
    #    lead variant, so 500 = the lead variant +/-500 kb region.
    # -------------------------------------------------------------------------
    echo "[SMR] running target=${target_snp}" | tee -a "$log_file"

    if "$SMR_BIN" \
        --bfile "$ld_prefix" \
        --gwas-summary "$gwas_ma" \
        --beqtl-summary "$beqtl_prefix" \
        --target-snp "$target_snp" \
        --cis-wind "$cis_wind" \
        --peqtl-smr "$peqtl_smr" \
        --peqtl-heidi "$peqtl_heidi" \
        --heidi-mtd "$heidi_mtd" \
        --heidi-min-m "$heidi_min_m" \
        --heidi-max-m "$heidi_max_m" \
        --ld-lower-limit "$ld_lower_limit" \
        --ld-upper-limit "$ld_upper_limit" \
        --maf "$maf" \
        --diff-freq "$diff_freq" \
        --diff-freq-prop "$diff_freq_prop" \
        --thread-num "$thread_num" \
        --out "$task_out_prefix" \
        >> "$log_file" 2>&1; then

        if [[ -s "$smr_file" ]]; then
            run_status="COMPLETED"
            msg="OK"
        else
            run_status="NO_SMR_OUTPUT"
            msg="SMR command exited 0 but .smr output is absent/empty"
        fi
    else
        run_status="SMR_FAILED"
        msg="SMR command returned non-zero exit status"
    fi

    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
      "$analysis_group" "$cell_type" "$gene" "$cohort" "$eqtl_source" \
      "$target_snp" "$run_status" "$n_ld_snps" "$smr_file" "$log_file" "$msg" >> "$RUN_STATUS"

done < "$TASK_STREAM"

# ==============================================================================
# Collect all successful .smr files and add task metadata.
# ==============================================================================
Rscript --vanilla - "$MANIFEST" "$RUN_STATUS" "$COMBINED" <<'RSCRIPT'
args <- commandArgs(trailingOnly = TRUE)
manifest_file <- args[1]
status_file <- args[2]
out_file <- args[3]

read_tsv <- function(path) {
  read.delim(
    path,
    header = TRUE,
    sep = "\t",
    quote = "",
    comment.char = "",
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
}

write_tsv <- function(x, path) {
  write.table(
    x,
    file = path,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    col.names = TRUE,
    na = "NA"
  )
}

m <- read_tsv(manifest_file)
s <- read_tsv(status_file)

# One task = one .smr path. Keying on it keeps gene x lead-variant pairs
# distinct, which a gene-level key would collapse.
m$smr_file <- paste0(m$task_out_prefix, ".smr")
key <- "smr_file"
meta_cols <- c(key, setdiff(names(m), names(s)))

meta <- m[, meta_cols, drop = FALSE]
meta <- meta[!duplicated(meta[[key]]), , drop = FALSE]

# Preserve runner order while adding manifest metadata.
s$.runner_order <- seq_len(nrow(s))
s2 <- merge(s, meta, by = key, all.x = TRUE, sort = FALSE)
s2 <- s2[order(s2$.runner_order), , drop = FALSE]
s2$.runner_order <- NULL

rows <- list()
k <- 0L

for (i in seq_len(nrow(s2))) {
  if (!identical(as.character(s2$run_status[i]), "COMPLETED")) next

  f <- as.character(s2$smr_file[i])
  if (is.na(f) || !file.exists(f) || file.info(f)$size <= 0) next

  z <- tryCatch(
    read_tsv(f),
    error = function(e) {
      warning("Could not read SMR result: ", f, " :: ", conditionMessage(e))
      NULL
    }
  )
  if (is.null(z) || !nrow(z)) next

  add_cols <- setdiff(names(s2), c("smr_file", "log_file", "message"))
  add <- s2[rep(i, nrow(z)), add_cols, drop = FALSE]

  k <- k + 1L
  rows[[k]] <- cbind(add, z)
}

if (length(rows)) {
  all_cols <- unique(unlist(lapply(rows, names), use.names = FALSE))

  rows <- lapply(rows, function(x) {
    miss <- setdiff(all_cols, names(x))
    if (length(miss)) {
      for (nm in miss) x[[nm]] <- NA
    }
    x[, all_cols, drop = FALSE]
  })

  out <- do.call(rbind, rows)
  write_tsv(out, out_file)
  message("Combined SMR+HEIDI results: ", out_file)
  message("Rows: ", nrow(out))
} else {
  out <- data.frame(
    analysis_group = character(),
    cell_type = character(),
    gene = character(),
    cohort = character(),
    eqtl_source = character(),
    target_SNP = character(),
    run_status = character(),
    stringsAsFactors = FALSE
  )
  write_tsv(out, out_file)
  message("No successful SMR result rows to combine. Empty result file: ", out_file)
}
RSCRIPT

echo
echo "======================================================================"
echo "SMR+HEIDI run finished"
echo "Run status : $RUN_STATUS"
echo "Combined   : $COMBINED"
echo "======================================================================"
