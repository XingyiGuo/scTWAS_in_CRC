#!/usr/bin/env bash
set -euo pipefail

# ==============================================================================
# Genotype_vcf_to_traw_HTAN_GSE266546_fixed.sh
#
# Input:
#   22 chromosome-specific merged VCFs:
#   /data/l2_bioinfo1/liq17/scTWAS/data/VUMC_HTAN_GSE266546/genotype/
#   VUMC_HTAN_GSE266546/R03/
#   chr1.HTAN_GSE266546.R03.missing20.vcf.gz
#   ...
#   chr22.HTAN_GSE266546.R03.missing20.vcf.gz
#
# Pipeline:
#   STEP 1  rename chromosome label chrN -> N separately for chr1-22
#   STEP 2  TBI-index renamed VCF
#   STEP 3  dbSNP annotation SKIPPED
#   STEP 4  renamed VCF -> PGEN using GT
#   STEP 5  raw genotype statistics
#   STEP 6  remove duplicate IDs + IDs beginning with "chr"
#   STEP 7  MAF / missingness / HWE QC
#   STEP 8  retain standard rsIDs
#   STEP 9  ALT-coded TRAW + TPED/TFAM + allele check
#   STEP 10 LD pruning
#   STEP 11 variant harmonization key
#   STEP 12 per-chromosome summary
#   STEP 13 merge chr1-22 LD-pruned PGEN + genome-wide PCA
#
# Run:
#   bash Genotype_vcf_to_traw_HTAN_GSE266546_fixed.sh
# ==============================================================================


# ==============================================================================
# 0. Software / CPU
# ==============================================================================

BCFTOOLS="/data/l2_bioinfo1/liq17/biosoft/bcftools/bcftools"
PLINK2="/data/l2_bioinfo1/liq17/biosoft/plink2"

N_PARALLEL=8
THREADS_PER_CHR=1

[[ -x "${BCFTOOLS}" ]] || {
    echo "ERROR: bcftools not executable: ${BCFTOOLS}" >&2
    exit 1
}

[[ -x "${PLINK2}" ]] || {
    echo "ERROR: plink2 not executable: ${PLINK2}" >&2
    exit 1
}


# ==============================================================================
# 1. Input / output
# ==============================================================================

BASE="/data/l2_bioinfo1/liq17/scTWAS/data/VUMC_HTAN_GSE266546/genotype/VUMC_HTAN_GSE266546"

# The 22 merged VCFs are already here.
INPUT_VCF_DIR="${BASE}/R03"

# Keep downstream outputs under the same R03 directory.
OUTDIR="${BASE}/R03"

VCF_DIR="${OUTDIR}/01_vcf_renamed"
RAW_PGEN_DIR="${OUTDIR}/02_pgen_raw"
DEDUP_DIR="${OUTDIR}/03_pgen_dedup"
QC_PGEN_DIR="${OUTDIR}/04_pgen_qc"
RSID_DIR="${OUTDIR}/05_pgen_rsid"
TRAW_DIR="${OUTDIR}/06_traw"
LD_DIR="${OUTDIR}/07_ld_prune"
PCA_DIR="${OUTDIR}/08_pgen_ldpruned_for_jointPCA"
JOINT_PCA_DIR="${OUTDIR}/09_jointPCA"
QC_DIR="${OUTDIR}/qc"
LOG_DIR="${OUTDIR}/logs"
TMP_DIR="${OUTDIR}/tmp"

mkdir -p \
    "${VCF_DIR}" \
    "${RAW_PGEN_DIR}" \
    "${DEDUP_DIR}" \
    "${QC_PGEN_DIR}" \
    "${RSID_DIR}" \
    "${TRAW_DIR}" \
    "${LD_DIR}" \
    "${PCA_DIR}" \
    "${JOINT_PCA_DIR}" \
    "${QC_DIR}" \
    "${LOG_DIR}" \
    "${TMP_DIR}"


# ==============================================================================
# 2. Parameters
# ==============================================================================

MAF_CUTOFF=0.01
GENO_CUTOFF=0.01
HWE_CUTOFF=1e-6

LD_WINDOW=200
LD_STEP=50
LD_R2=0.2

N_PCS=20


# ==============================================================================
# 3. Preflight
# ==============================================================================

echo "============================================================"
echo "PRE-FLIGHT"
echo "============================================================"
echo "Input directory       : ${INPUT_VCF_DIR}"
echo "Parallel chromosomes  : ${N_PARALLEL}"
echo "Threads/chromosome    : ${THREADS_PER_CHR}"
echo "PLINK genotype field  : GT"
echo "MAF cutoff            : ${MAF_CUTOFF}"
echo "GENO cutoff           : ${GENO_CUTOFF}"
echo "HWE cutoff            : ${HWE_CUTOFF}"
echo


# ------------------------------------------------------------------------------
# 3A. Check all chr1-22 input VCFs
# ------------------------------------------------------------------------------

for i in $(seq 1 22); do
    VCF="${INPUT_VCF_DIR}/chr${i}.HTAN_GSE266546.R03.missing20.vcf.gz"

    [[ -s "${VCF}" ]] || {
        echo "ERROR: missing input VCF: ${VCF}" >&2
        exit 1
    }

    # Create an index if neither TBI nor CSI exists.
    if [[ ! -s "${VCF}.tbi" && ! -s "${VCF}.csi" ]]; then
        echo "Creating TBI index: ${VCF}"
        "${BCFTOOLS}" index -t -f "${VCF}"
    fi

    # This bcftools installation does not provide `quickcheck`.
    # A readable header is used as the basic input validation.
    "${BCFTOOLS}" view -h "${VCF}" >/dev/null
done

echo "Input VCFs chr1-22: OK"


# ------------------------------------------------------------------------------
# 3B. Chromosome rename map: chr1 -> 1, ..., chr22 -> 22
# ------------------------------------------------------------------------------

CHR_RENAME_FILE="${TMP_DIR}/rename_chr_to_num.txt"
: > "${CHR_RENAME_FILE}"

for i in $(seq 1 22); do
    printf "chr%s\t%s\n" "${i}" "${i}" >> "${CHR_RENAME_FILE}"
done

echo "Chromosome rename map: ${CHR_RENAME_FILE}"


# ==============================================================================
# 4. Per-chromosome function: STEP 1-12
# ==============================================================================

process_chr() {
    set -euo pipefail

    local i="$1"

    local INPUT_VCF="${INPUT_VCF_DIR}/chr${i}.HTAN_GSE266546.R03.missing20.vcf.gz"
    local PREFIX="chr${i}.HTAN_GSE266546"

    # STEP 1 output: the only change is chromosome label chrN -> N.
    local RENAMED_VCF="${VCF_DIR}/chr${i}.HTAN_GSE266546.R03.missing20.renamed.vcf.gz"

    local RAW_PREFIX="${RAW_PGEN_DIR}/${PREFIX}.raw"
    local DEDUP_PREFIX="${DEDUP_DIR}/${PREFIX}.dedup"
    local QC_PREFIX="${QC_PGEN_DIR}/${PREFIX}.qc"
    local RSID_PREFIX="${RSID_DIR}/${PREFIX}.qc.rsid"

    local TRAW_PREFIX="${TRAW_DIR}/${PREFIX}.qc.rsid"
    local TPED_PREFIX="${TRAW_DIR}/${PREFIX}.qc.rsid"

    local LD_PREFIX="${LD_DIR}/${PREFIX}.qc.rsid"
    local PCA_PREFIX="${PCA_DIR}/${PREFIX}.qc.rsid.ldpruned"

    local DUP_FILE="${QC_DIR}/${PREFIX}.duplicate_IDs.tsv"
    local CHR_ID_FILE="${QC_DIR}/${PREFIX}.chr_IDs.txt"
    local EXCLUDE_FILE="${QC_DIR}/${PREFIX}.exclude_duplicate_or_chrID.txt"

    local RSID_KEEP="${QC_DIR}/${PREFIX}.rsid.keep"
    local ALT_ALLELE_FILE="${QC_DIR}/${PREFIX}.ALT_allele.txt"
    local TRAW_MISMATCH="${QC_DIR}/${PREFIX}.traw_ALT_mismatch.tsv"

    local VARIANT_KEY="${QC_DIR}/${PREFIX}.variant_key.tsv"
    local SUMMARY="${QC_DIR}/${PREFIX}.summary.txt"

    local N_DUP_IDS
    local N_CHR_IDS
    local N_EXCLUDE
    local N_RSID_KEEP
    local N_TRAW_MISMATCH
    local N_INPUT
    local N_RAW
    local N_DEDUP
    local N_QC
    local N_RSID
    local N_LD
    local N_SAMPLES

    echo "===== chr${i} START ====="


    # ==========================================================================
    # STEP 1
    # Rename chromosome for this chromosome-specific VCF: chrN -> N
    #
    # No sample subset.
    # No additional R2 filtering.
    # No dbSNP re-annotation.
    # ==========================================================================

    "${BCFTOOLS}" annotate \
        --rename-chrs "${CHR_RENAME_FILE}" \
        -Oz \
        -o "${RENAMED_VCF}" \
        "${INPUT_VCF}"

    echo "chr${i}: STEP1 chr${i} -> ${i} rename OK"


    # ==========================================================================
    # STEP 2
    # Index renamed VCF
    # ==========================================================================

    "${BCFTOOLS}" index \
        -t \
        -f \
        "${RENAMED_VCF}"

    "${BCFTOOLS}" view -h "${RENAMED_VCF}" >/dev/null
    "${BCFTOOLS}" index -n "${RENAMED_VCF}" >/dev/null

    echo "chr${i}: STEP2 renamed VCF TBI index OK"


    # ==========================================================================
    # STEP 3
    # dbSNP annotation SKIPPED
    #
    # Existing IDs from the merged VCF are retained.
    # Missing IDs are filled by PLINK2 in STEP 4.
    # ==========================================================================


    # ==========================================================================
    # STEP 4
    # Renamed VCF -> PGEN using GT
    # ==========================================================================

    "${PLINK2}" \
        --vcf "${RENAMED_VCF}" \
        --set-missing-var-ids '@:#:$r:$a' \
        --make-pgen \
        --threads "${THREADS_PER_CHR}" \
        --out "${RAW_PREFIX}"

    echo "chr${i}: STEP4 renamed VCF -> PGEN (GT) OK"


    # ==========================================================================
    # STEP 5
    # Raw statistics
    # ==========================================================================

    "${PLINK2}" \
        --pfile "${RAW_PREFIX}" \
        --freq \
        --missing \
        --hardy \
        --threads "${THREADS_PER_CHR}" \
        --out "${QC_DIR}/${PREFIX}.raw_stats"

    echo "chr${i}: STEP5 raw statistics OK"


    # ==========================================================================
    # STEP 6
    # Remove all variants belonging to duplicate-ID groups
    # + remove IDs beginning with "chr"
    # ==========================================================================

    awk '
        BEGIN {FS=OFS="\t"}
        !/^#/ {n[$3]++}
        END {
            for (id in n)
                if (n[id] > 1)
                    print id,n[id]
        }
    ' "${RAW_PREFIX}.pvar" \
    | sort \
    > "${DUP_FILE}"

    awk '
        BEGIN {FS=OFS="\t"}
        !/^#/ && $3 ~ /^chr/ {
            print $3
        }
    ' "${RAW_PREFIX}.pvar" \
    | sort -u \
    > "${CHR_ID_FILE}"

    {
        cut -f1 "${DUP_FILE}"
        cat "${CHR_ID_FILE}"
    } \
    | sed '/^$/d' \
    | sort -u \
    > "${EXCLUDE_FILE}"

    N_DUP_IDS="$(wc -l < "${DUP_FILE}")"
    N_CHR_IDS="$(wc -l < "${CHR_ID_FILE}")"
    N_EXCLUDE="$(wc -l < "${EXCLUDE_FILE}")"

    if [[ "${N_EXCLUDE}" -gt 0 ]]; then
        "${PLINK2}" \
            --pfile "${RAW_PREFIX}" \
            --exclude "${EXCLUDE_FILE}" \
            --make-pgen \
            --threads "${THREADS_PER_CHR}" \
            --out "${DEDUP_PREFIX}"
    else
        "${PLINK2}" \
            --pfile "${RAW_PREFIX}" \
            --make-pgen \
            --threads "${THREADS_PER_CHR}" \
            --out "${DEDUP_PREFIX}"
    fi

    echo "chr${i}: STEP6 duplicate-ID groups = ${N_DUP_IDS}"
    echo "chr${i}: STEP6 chr*-ID variants    = ${N_CHR_IDS}"
    echo "chr${i}: STEP6 total excluded      = ${N_EXCLUDE}"


    # ==========================================================================
    # STEP 7
    # Final QC
    # ==========================================================================

    "${PLINK2}" \
        --pfile "${DEDUP_PREFIX}" \
        --maf "${MAF_CUTOFF}" \
        --geno "${GENO_CUTOFF}" \
        --hwe "${HWE_CUTOFF}" \
        --make-pgen \
        --threads "${THREADS_PER_CHR}" \
        --out "${QC_PREFIX}"

    "${PLINK2}" \
        --pfile "${QC_PREFIX}" \
        --freq \
        --missing \
        --hardy \
        --threads "${THREADS_PER_CHR}" \
        --out "${QC_DIR}/${PREFIX}.postQC_stats"

    echo "chr${i}: STEP7 MAF/missingness/HWE QC OK"


    # ==========================================================================
    # STEP 8
    # Keep standard rsIDs only
    # ==========================================================================

    awk '
        !/^#/ && $3 ~ /^rs[0-9]+$/ {
            print $3
        }
    ' "${QC_PREFIX}.pvar" \
    > "${RSID_KEEP}"

    N_RSID_KEEP="$(wc -l < "${RSID_KEEP}")"

    if [[ "${N_RSID_KEEP}" -eq 0 ]]; then
        echo "ERROR chr${i}: no standard rsIDs remain after QC." >&2
        return 1
    fi

    "${PLINK2}" \
        --pfile "${QC_PREFIX}" \
        --extract "${RSID_KEEP}" \
        --make-pgen \
        --threads "${THREADS_PER_CHR}" \
        --out "${RSID_PREFIX}"

    echo "chr${i}: STEP8 rsID variants retained = ${N_RSID_KEEP}"


    # ==========================================================================
    # STEP 9A
    # ID + ALT allele file
    # ==========================================================================

    awk '
        BEGIN {FS=OFS="\t"}
        !/^#/ {
            print $3,$5
        }
    ' "${RSID_PREFIX}.pvar" \
    > "${ALT_ALLELE_FILE}"


    # ==========================================================================
    # STEP 9B
    # ALT-coded TRAW
    #
    # REF/REF = 0
    # REF/ALT = 1
    # ALT/ALT = 2
    # ==========================================================================

    "${PLINK2}" \
        --pfile "${RSID_PREFIX}" \
        --export Av \
        --export-allele "${ALT_ALLELE_FILE}" \
        --threads "${THREADS_PER_CHR}" \
        --out "${TRAW_PREFIX}"

    echo "chr${i}: STEP9B ALT-coded .traw export OK"


    # ==========================================================================
    # STEP 9C
    # TPED/TFAM
    # ==========================================================================

    "${PLINK2}" \
        --pfile "${RSID_PREFIX}" \
        --export tped \
        --threads "${THREADS_PER_CHR}" \
        --out "${TPED_PREFIX}"

    echo "chr${i}: STEP9C TPED/TFAM export OK"


    # ==========================================================================
    # STEP 9D
    # Verify TRAW counted allele == PVAR ALT
    # ==========================================================================

    awk '
        BEGIN {FS=OFS="\t"}

        NR==FNR {
            if ($0 !~ /^#/)
                alt[$3]=$5
            next
        }

        FNR==1 {next}

        {
            id=$2
            counted=$5
            if ((id in alt) && counted != alt[id])
                print id,alt[id],counted
        }
    ' "${RSID_PREFIX}.pvar" "${TRAW_PREFIX}.traw" \
    > "${TRAW_MISMATCH}"

    N_TRAW_MISMATCH="$(wc -l < "${TRAW_MISMATCH}")"

    if [[ "${N_TRAW_MISMATCH}" -ne 0 ]]; then
        echo "ERROR chr${i}: ${N_TRAW_MISMATCH} TRAW ALT-coding mismatches." >&2
        head "${TRAW_MISMATCH}" >&2
        return 1
    fi

    echo "chr${i}: STEP9D ALT coding check PASSED"


    # ==========================================================================
    # STEP 10
    # LD pruning
    # ==========================================================================

    "${PLINK2}" \
        --pfile "${RSID_PREFIX}" \
        --indep-pairwise "${LD_WINDOW}" "${LD_STEP}" "${LD_R2}" \
        --threads "${THREADS_PER_CHR}" \
        --out "${LD_PREFIX}"

    "${PLINK2}" \
        --pfile "${RSID_PREFIX}" \
        --extract "${LD_PREFIX}.prune.in" \
        --make-pgen \
        --threads "${THREADS_PER_CHR}" \
        --out "${PCA_PREFIX}"

    echo "chr${i}: STEP10 LD pruning OK"


    # ==========================================================================
    # STEP 11
    # Variant harmonization key
    # ==========================================================================

    awk '
        BEGIN {
            FS=OFS="\t"
            print "CHR","POS","ID","REF","ALT","variant_key"
        }

        !/^#/ {
            print $1,$2,$3,$4,$5,$1":"$2":"$4":"$5
        }
    ' "${RSID_PREFIX}.pvar" \
    > "${VARIANT_KEY}"

    echo "chr${i}: STEP11 variant key OK"


    # ==========================================================================
    # STEP 12
    # Per-chromosome summary
    # ==========================================================================

    N_INPUT="$("${BCFTOOLS}" index -n "${RENAMED_VCF}")"
    N_SAMPLES="$("${BCFTOOLS}" query -l "${RENAMED_VCF}" | wc -l)"
    N_RAW="$(awk '!/^#/ {n++} END{print n+0}' "${RAW_PREFIX}.pvar")"
    N_DEDUP="$(awk '!/^#/ {n++} END{print n+0}' "${DEDUP_PREFIX}.pvar")"
    N_QC="$(awk '!/^#/ {n++} END{print n+0}' "${QC_PREFIX}.pvar")"
    N_RSID="$(awk '!/^#/ {n++} END{print n+0}' "${RSID_PREFIX}.pvar")"
    N_LD="$(wc -l < "${LD_PREFIX}.prune.in")"

    {
        echo "chromosome=${i}"
        echo "input_vcf=${INPUT_VCF}"
        echo "renamed_vcf=${RENAMED_VCF}"
        echo "genotype_source=GT"
        echo "n_samples=${N_SAMPLES}"
        echo "n_variants_input=${N_INPUT}"
        echo "n_variants_raw_pgen=${N_RAW}"
        echo "duplicate_ID_groups=${N_DUP_IDS}"
        echo "chr_ID_variants_removed=${N_CHR_IDS}"
        echo "total_preQC_excluded=${N_EXCLUDE}"
        echo "n_variants_after_preQC_cleanup=${N_DEDUP}"
        echo "maf_cutoff=${MAF_CUTOFF}"
        echo "geno_cutoff=${GENO_CUTOFF}"
        echo "hwe_cutoff=${HWE_CUTOFF}"
        echo "n_variants_after_MAF_GENO_HWE=${N_QC}"
        echo "n_rsid_variants=${N_RSID}"
        echo "traw_ALT_mismatches=${N_TRAW_MISMATCH}"
        echo "n_ld_pruned_variants=${N_LD}"
        echo "traw=${TRAW_PREFIX}.traw"
        echo "tped=${TPED_PREFIX}.tped"
        echo "tfam=${TPED_PREFIX}.tfam"
        echo "variant_key=${VARIANT_KEY}"
        echo "ld_pruned_pgen_prefix=${PCA_PREFIX}"
    } > "${SUMMARY}"

    echo "===== chr${i} DONE ====="
}


# ==============================================================================
# 5. Run STEP 1-12 for chr1-22
#
# Native Bash batches are used instead of GNU parallel to avoid the prior
# bash -c / $1 / set -u argument-passing problem.
# ==============================================================================

echo
echo "============================================================"
echo "RUN chr1-22"
echo "============================================================"

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
            process_chr "${i}"
        ) > "${LOG_DIR}/chr${i}.pipeline.log" 2>&1 &

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
            echo "ERROR: chr${chr} failed. See ${LOG_DIR}/chr${chr}.pipeline.log" >&2
            tail -60 "${LOG_DIR}/chr${chr}.pipeline.log" >&2 || true
            batch_failed=1
        fi
    done

    if (( batch_failed != 0 )); then
        echo "ERROR: stopping because one or more chromosomes failed." >&2
        exit 1
    fi
done


# ==============================================================================
# STEP 13
# Merge chr1-22 LD-pruned PGEN + genome-wide PCA
# ==============================================================================

echo
echo "============================================================"
echo "STEP 13: MERGE chr1-22 LD-PRUNED PGEN + PCA"
echo "============================================================"

PMERGE_LIST="${JOINT_PCA_DIR}/chr2_22.pmerge_list.txt"
MERGED_PREFIX="${JOINT_PCA_DIR}/HTAN_GSE266546.chr1_22.LDpruned"
PCA_OUT="${JOINT_PCA_DIR}/HTAN_GSE266546.chr1_22.PCA"

for i in $(seq 1 22); do
    P="${PCA_DIR}/chr${i}.HTAN_GSE266546.qc.rsid.ldpruned"

    for ext in pgen pvar psam; do
        if [[ ! -s "${P}.${ext}" ]]; then
            echo "ERROR: missing ${P}.${ext}" >&2
            exit 1
        fi
    done
done

: > "${PMERGE_LIST}"

for i in $(seq 2 22); do
    echo "${PCA_DIR}/chr${i}.HTAN_GSE266546.qc.rsid.ldpruned" \
        >> "${PMERGE_LIST}"
done

"${PLINK2}" \
    --pfile "${PCA_DIR}/chr1.HTAN_GSE266546.qc.rsid.ldpruned" \
    --pmerge-list "${PMERGE_LIST}" \
    --make-pgen \
    --threads "${N_PARALLEL}" \
    --out "${MERGED_PREFIX}"

echo "STEP13: chr1-22 LD-pruned PGEN merge OK"

"${PLINK2}" \
    --pfile "${MERGED_PREFIX}" \
    --pca "${N_PCS}" approx \
    --threads "${N_PARALLEL}" \
    --out "${PCA_OUT}"

echo "STEP13: genome-wide PCA completed"


# ==============================================================================
# Final overview
# ==============================================================================

echo
echo "============================================================"
echo "ALL STEPS COMPLETED"
echo "============================================================"

N_VCF="$(find "${VCF_DIR}" -maxdepth 1 -type f -name 'chr*.renamed.vcf.gz' | wc -l)"
N_TRAW="$(find "${TRAW_DIR}" -maxdepth 1 -type f -name '*.traw' | wc -l)"
N_TPED="$(find "${TRAW_DIR}" -maxdepth 1 -type f -name '*.tped' | wc -l)"
N_TFAM="$(find "${TRAW_DIR}" -maxdepth 1 -type f -name '*.tfam' | wc -l)"
N_SUMMARY="$(find "${QC_DIR}" -maxdepth 1 -type f -name 'chr*.summary.txt' | wc -l)"

echo "Renamed VCFs : ${N_VCF}/22"
echo ".traw files  : ${N_TRAW}/22"
echo ".tped files  : ${N_TPED}/22"
echo ".tfam files  : ${N_TFAM}/22"
echo "summary files: ${N_SUMMARY}/22"

echo
echo "Input VCFs      : ${INPUT_VCF_DIR}/chr*.HTAN_GSE266546.R03.missing20.vcf.gz"
echo "Renamed VCFs    : ${VCF_DIR}"
echo "rsID PGEN       : ${RSID_DIR}"
echo "TRAW/TPED       : ${TRAW_DIR}"
echo "LD prune        : ${LD_DIR}"
echo "PCA chromosome  : ${PCA_DIR}"
echo "PCA merged      : ${JOINT_PCA_DIR}"
echo "Eigenvec        : ${PCA_OUT}.eigenvec"
echo "Eigenval        : ${PCA_OUT}.eigenval"
echo "QC              : ${QC_DIR}"
echo "Chromosome logs : ${LOG_DIR}/chr*.pipeline.log"
