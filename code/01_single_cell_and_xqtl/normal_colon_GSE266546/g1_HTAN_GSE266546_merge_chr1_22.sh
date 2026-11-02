#!/usr/bin/env bash
set -euo pipefail

# ==============================================================================
# 01_merge_HTAN_GSE266546_chr1_22.sh
#
# Merge HTAN + GSE266546 VCFs chromosome-by-chromosome (chr1-22).
#
# Rules:
#   1. Keep ALL samples from both input VCFs.
#   2. Merge the union of variant sites across the two cohorts.
#   3. Do NOT convert cohort-absent genotypes to reference; they remain missing.
#   4. Remove any variant with >20% missing GT across ALL merged subjects:
#          F_MISSING <= 0.20
#   5. Preserve chromosome-specific VCFs (chr1-22).
#
# Parallelization:
#   - 8 chromosomes at a time
#   - 1 bcftools thread per chromosome
#
# Run:
#   bash 01_merge_HTAN_GSE266546_chr1_22.sh
# ==============================================================================


# ==============================================================================
# 0. Software / CPU
# ==============================================================================

BCFTOOLS="/data/l2_bioinfo1/liq17/biosoft/bcftools/bcftools"

N_PARALLEL=8
THREADS_PER_CHR=1
MAX_MISSING=0.20

[[ -x "${BCFTOOLS}" ]] || {
    echo "ERROR: bcftools not executable: ${BCFTOOLS}" >&2
    exit 1
}


# ==============================================================================
# 1. Input / output
# ==============================================================================

BASE="/data/l2_bioinfo1/liq17/scTWAS/data/VUMC_HTAN_GSE266546/genotype"

HTAN_DIR="${BASE}/VUMC_HTAN"

GSE_DIR="${BASE}/VUMC_GSE266546/rsq0.3maf0.01/biallelic_snps"

OUTDIR="${BASE}/VUMC_HTAN_GSE266546/R03"
TMP_DIR="${OUTDIR}/tmp"
QC_DIR="${OUTDIR}/qc"
LOG_DIR="${OUTDIR}/logs"

mkdir -p \
    "${OUTDIR}" \
    "${TMP_DIR}" \
    "${QC_DIR}" \
    "${LOG_DIR}"


# ==============================================================================
# 2. Preflight
# ==============================================================================

echo "============================================================"
echo "PRE-FLIGHT"
echo "============================================================"
echo "HTAN      : ${HTAN_DIR}"
echo "GSE266546 : ${GSE_DIR}"
echo "Output    : ${OUTDIR}"
echo "Filter    : F_MISSING <= ${MAX_MISSING}"
echo "Parallel  : ${N_PARALLEL} chromosomes"
echo


# ------------------------------------------------------------------------------
# 2A. Check chr1-22 input VCFs and indexes
# ------------------------------------------------------------------------------

for i in $(seq 8 22); do

    HTAN_VCF="${HTAN_DIR}/chr${i}.R03.MAF0001.95samples.dose.vcf.gz"

    GSE_VCF="${GSE_DIR}/chr${i}.Rsq0.3.maf0.01.biallelic.snps.vcf.gz"

    for VCF in "${HTAN_VCF}" "${GSE_VCF}"; do

        [[ -s "${VCF}" ]] || {
            echo "ERROR: missing input VCF: ${VCF}" >&2
            exit 1
        }

        if [[ ! -s "${VCF}.tbi" && ! -s "${VCF}.csi" ]]; then
            echo "Index missing; creating TBI:"
            echo "  ${VCF}"
            "${BCFTOOLS}" index -t -f "${VCF}"
        fi

        # Validate that bcftools can read the VCF header.
        "${BCFTOOLS}" view -h "${VCF}" >/dev/null
    done
done

echo "chr1-22 input VCFs + indexes: OK"


# ------------------------------------------------------------------------------
# 2B. Require FORMAT/GT in both cohorts.
#
# F_MISSING is defined from genotype calls, so this workflow should not be run
# on a dosage-only VCF which lacks GT.
# ------------------------------------------------------------------------------

HTAN_TEST="${HTAN_DIR}/chr1.R03.MAF0001.95samples.dose.vcf.gz"
GSE_TEST="${GSE_DIR}/chr1.Rsq0.3.maf0.01.biallelic.snps.vcf.gz"

for VCF in "${HTAN_TEST}" "${GSE_TEST}"; do
    if ! "${BCFTOOLS}" view -h "${VCF}" \
        | grep -q '^##FORMAT=<ID=GT,'; then

        echo "ERROR: FORMAT/GT not found in:" >&2
        echo "  ${VCF}" >&2
        echo "Missingness filter F_MISSING requires GT." >&2
        exit 1
    fi
done

echo "FORMAT/GT present in both cohorts: OK"


# ------------------------------------------------------------------------------
# 2C. Check sample names.
#     Duplicate sample IDs across two VCFs are NOT force-renamed.
# ------------------------------------------------------------------------------

HTAN_SAMPLES="${QC_DIR}/HTAN_samples.txt"
GSE_SAMPLES="${QC_DIR}/GSE266546_samples.txt"
OVERLAP_SAMPLES="${QC_DIR}/HTAN_GSE266546_overlapping_sampleIDs.txt"

"${BCFTOOLS}" query -l "${HTAN_TEST}" \
    | sort -u \
    > "${HTAN_SAMPLES}"

"${BCFTOOLS}" query -l "${GSE_TEST}" \
    | sort -u \
    > "${GSE_SAMPLES}"

comm -12 \
    "${HTAN_SAMPLES}" \
    "${GSE_SAMPLES}" \
    > "${OVERLAP_SAMPLES}"

N_HTAN="$(wc -l < "${HTAN_SAMPLES}")"
N_GSE="$(wc -l < "${GSE_SAMPLES}")"
N_OVERLAP="$(wc -l < "${OVERLAP_SAMPLES}")"
N_EXPECTED=$((N_HTAN + N_GSE))

echo "HTAN samples        : ${N_HTAN}"
echo "GSE266546 samples   : ${N_GSE}"
echo "Overlapping IDs     : ${N_OVERLAP}"

if [[ "${N_OVERLAP}" -gt 0 ]]; then
    echo "ERROR: overlapping sample IDs detected between cohorts." >&2
    echo "These are not automatically renamed because that could alter subject identity." >&2
    cat "${OVERLAP_SAMPLES}" >&2
    exit 1
fi

echo "Expected merged N   : ${N_EXPECTED}"


# ==============================================================================
# 3. Per-chromosome merge function
# ==============================================================================

merge_chr() {
    set -euo pipefail

    local i="$1"

    local HTAN_VCF="${HTAN_DIR}/chr${i}.R03.MAF0001.95samples.dose.vcf.gz"

    local GSE_VCF="${GSE_DIR}/chr${i}.Rsq0.3.maf0.01.biallelic.snps.vcf.gz"

    local TMP_MERGED="${TMP_DIR}/chr${i}.HTAN_GSE266546.unfiltered.vcf.gz"

    local OUT_VCF="${OUTDIR}/chr${i}.HTAN_GSE266546.R03.missing20.vcf.gz"

    local SUMMARY="${QC_DIR}/chr${i}.merge_missing20.summary.txt"

    echo "===== chr${i} START ====="


    # --------------------------------------------------------------------------
    # STEP 1. Merge the two cohorts.
    #
    # bcftools merge uses the UNION of sites.
    #
    # IMPORTANT:
    # We intentionally DO NOT use --missing-to-ref.
    # If a variant is absent from one cohort, those subjects remain missing,
    # which is then correctly reflected in the post-merge F_MISSING filter.
    #
    # -m none prevents creation of new multiallelic records when different ALT
    # alleles occur at the same position.
    # --------------------------------------------------------------------------

    "${BCFTOOLS}" merge \
        -m none \
        --threads "${THREADS_PER_CHR}" \
        -Oz \
        -o "${TMP_MERGED}" \
        "${HTAN_VCF}" \
        "${GSE_VCF}"

    "${BCFTOOLS}" index \
        -t \
        -f \
        "${TMP_MERGED}"

    # Validate merged VCF is readable and indexed.
    "${BCFTOOLS}" view -h "${TMP_MERGED}" >/dev/null
    "${BCFTOOLS}" index -n "${TMP_MERGED}" >/dev/null

    echo "chr${i}: cohort merge OK"


    # --------------------------------------------------------------------------
    # STEP 2. Remove variants missing in >20% of ALL merged subjects.
    #
    # Keep:
    #       F_MISSING <= 0.20
    #
    # This is a SITE filter only. No subjects are removed.
    # --------------------------------------------------------------------------

    "${BCFTOOLS}" view \
        -i "F_MISSING<=${MAX_MISSING}" \
        --threads "${THREADS_PER_CHR}" \
        -Oz \
        -o "${OUT_VCF}" \
        "${TMP_MERGED}"

    "${BCFTOOLS}" index \
        -t \
        -f \
        "${OUT_VCF}"

    # Validate final VCF is readable and indexed.
    "${BCFTOOLS}" view -h "${OUT_VCF}" >/dev/null
    "${BCFTOOLS}" index -n "${OUT_VCF}" >/dev/null

    echo "chr${i}: missingness <=20% filter OK"


    # --------------------------------------------------------------------------
    # STEP 3. QC / verify all samples were retained
    # --------------------------------------------------------------------------

    local N_SAMPLES_MERGED
    local N_SAMPLES_FINAL
    local N_VARIANTS_MERGED
    local N_VARIANTS_FINAL
    local N_REMOVED

    N_SAMPLES_MERGED="$(
        "${BCFTOOLS}" query -l "${TMP_MERGED}" | wc -l
    )"

    N_SAMPLES_FINAL="$(
        "${BCFTOOLS}" query -l "${OUT_VCF}" | wc -l
    )"

    N_VARIANTS_MERGED="$(
        "${BCFTOOLS}" index -n "${TMP_MERGED}"
    )"

    N_VARIANTS_FINAL="$(
        "${BCFTOOLS}" index -n "${OUT_VCF}"
    )"

    N_REMOVED=$((N_VARIANTS_MERGED - N_VARIANTS_FINAL))

    if [[ "${N_SAMPLES_MERGED}" -ne "${N_EXPECTED}" ]]; then
        echo "ERROR chr${i}: merged sample N=${N_SAMPLES_MERGED}, expected ${N_EXPECTED}" >&2
        return 1
    fi

    if [[ "${N_SAMPLES_FINAL}" -ne "${N_EXPECTED}" ]]; then
        echo "ERROR chr${i}: final sample N=${N_SAMPLES_FINAL}, expected ${N_EXPECTED}" >&2
        return 1
    fi

    {
        echo "chromosome=chr${i}"
        echo "HTAN_vcf=${HTAN_VCF}"
        echo "GSE266546_vcf=${GSE_VCF}"
        echo "n_HTAN_samples=${N_HTAN}"
        echo "n_GSE266546_samples=${N_GSE}"
        echo "n_merged_samples=${N_SAMPLES_MERGED}"
        echo "n_final_samples=${N_SAMPLES_FINAL}"
        echo "missingness_rule=F_MISSING<=${MAX_MISSING}"
        echo "n_variants_before_missingness_filter=${N_VARIANTS_MERGED}"
        echo "n_variants_after_missingness_filter=${N_VARIANTS_FINAL}"
        echo "n_variants_removed_missing_gt_gt20pct=${N_REMOVED}"
        echo "output_vcf=${OUT_VCF}"
    } > "${SUMMARY}"

    # Remove only temporary unfiltered merged VCF after successful final QC.
    rm -f \
        "${TMP_MERGED}" \
        "${TMP_MERGED}.tbi"

    echo "chr${i}: samples retained ${N_SAMPLES_FINAL}/${N_EXPECTED}"
    echo "chr${i}: variants ${N_VARIANTS_MERGED} -> ${N_VARIANTS_FINAL}"
    echo "===== chr${i} DONE ====="
}


# ==============================================================================
# 4. Run chr1-22
# ==============================================================================

echo
echo "============================================================"
echo "MERGE HTAN + GSE266546 chr1-22"
echo "============================================================"

# Use native Bash background jobs instead of GNU parallel.
# This avoids bash -c / exported-function quoting issues with `set -u`.
#
# Run chromosomes in batches of N_PARALLEL.  Each chromosome writes its own log.
for ((batch_start=1; batch_start<=22; batch_start+=N_PARALLEL)); do

    batch_end=$((batch_start + N_PARALLEL - 1))
    if (( batch_end > 22 )); then
        batch_end=22
    fi

    echo
    echo "Starting batch: chr${batch_start}-chr${batch_end}"

    pids=()
    chrs=()

    for ((i=batch_start; i<=batch_end; i++)); do
        (
            merge_chr "${i}"
        ) > "${LOG_DIR}/chr${i}.log" 2>&1 &

        pids+=("$!")
        chrs+=("${i}")
    done

    batch_failed=0

    for idx in "${!pids[@]}"; do
        pid="${pids[$idx]}"
        chr="${chrs[$idx]}"

        if wait "${pid}"; then
            echo "chr${chr}: completed"
        else
            echo "ERROR: chr${chr} failed. See ${LOG_DIR}/chr${chr}.log" >&2
            tail -50 "${LOG_DIR}/chr${chr}.log" >&2 || true
            batch_failed=1
        fi
    done

    if (( batch_failed != 0 )); then
        echo "ERROR: stopping because one or more chromosomes failed." >&2
        exit 1
    fi
done


# ==============================================================================
# 5. Final checks
# ==============================================================================

echo
echo "============================================================"
echo "FINAL CHECK"
echo "============================================================"

N_FINAL_VCF="$(
    find "${OUTDIR}" \
        -maxdepth 1 \
        -type f \
        -name 'chr*.HTAN_GSE266546.R03.missing20.vcf.gz' \
        | wc -l
)"

N_FINAL_TBI="$(
    find "${OUTDIR}" \
        -maxdepth 1 \
        -type f \
        -name 'chr*.HTAN_GSE266546.R03.missing20.vcf.gz.tbi' \
        | wc -l
)"

N_SUMMARY="$(
    find "${QC_DIR}" \
        -maxdepth 1 \
        -type f \
        -name 'chr*.merge_missing20.summary.txt' \
        | wc -l
)"

echo "Final VCFs      : ${N_FINAL_VCF}/22"
echo "Final TBI files : ${N_FINAL_TBI}/22"
echo "QC summaries    : ${N_SUMMARY}/22"
echo "Samples retained: ${N_EXPECTED}"

if [[ "${N_FINAL_VCF}" -ne 22 || "${N_FINAL_TBI}" -ne 22 ]]; then
    echo "ERROR: expected 22 final VCFs + indexes." >&2
    exit 1
fi

echo
echo "Output pattern:"
echo "  ${OUTDIR}/chr*.HTAN_GSE266546.R03.missing20.vcf.gz"
echo
echo "Done."
