#!/usr/bin/env bash
set -euo pipefail

# ==============================================================================
# genotype_vcf_to_traw_interactive.sh
#
# Full interactive genotype pipeline, chr1-22:
#   - 8 chromosomes in parallel
#   - 1 CPU per chromosome
#   - PLINK2 uses GT (NOT DS)
#
# STEP 1  sample subset + R2 > 0.30 + biallelic SNP + chrN -> N
# STEP 2  TBI index temporary renamed VCF
# STEP 3  chromosome-specific dbSNP rsID annotation + TBI
# STEP 4  annotated VCF -> PGEN (GT)
# STEP 5  raw frequency / missingness / HWE statistics
# STEP 6  remove duplicate IDs + IDs beginning with "chr"
# STEP 7  QC: MAF >= 0.01, geno <= 0.01, HWE >= 1e-6
# STEP 8  retain standard rsIDs
# STEP 9  ALT-coded .traw + PLINK2 .tped/.tfam + allele check
# STEP 10 LD pruning
# STEP 11 variant harmonization key
# STEP 12 per-chromosome summary
# STEP 13 merge chr1-22 LD-pruned PGEN + genome-wide PCA
#
# Run on an interactive node with >=8 CPUs:
#   bash genotype_vcf_to_traw_interactive.sh
# ==============================================================================


# ==============================================================================
# 0. Software / CPU
# ==============================================================================

BCFTOOLS="/data/l2_bioinfo1/liq17/biosoft/bcftools/bcftools"
PLINK2="/data/l2_bioinfo1/liq17/biosoft/plink2"

N_PARALLEL=8
THREADS_PER_CHR=1

command -v parallel >/dev/null 2>&1 || {
    echo "ERROR: GNU parallel not found in PATH." >&2
    exit 1
}

[[ -x "${BCFTOOLS}" ]] || {
    echo "ERROR: bcftools not executable: ${BCFTOOLS}" >&2
    exit 1
}

[[ -x "${PLINK2}" ]] || {
    echo "ERROR: plink2 not executable: ${PLINK2}" >&2
    exit 1
}


# ==============================================================================
# 1. Inputs
# ==============================================================================

BASE="/data/l2_bioinfo1/scRNA-seq/WellcomeSanger/processed/genotype"
VCF="${BASE}/imputed.vcf.gz"

SAMPLE_FILE="/data/l2_bioinfo1/scRNA-seq/WellcomeSanger/processed/E-MTAB-16986/multi_tissue_predicted_summary/one_sample_per_individual_Healthy_TI_then_rectum_individual_metadata_genoID.txt"

DBSNP_DIR="/data/l2_bioinfo1/liq17/LQ_references/NCBI/dbsnp_157/GCF_000001405.40_chr"


# ==============================================================================
# 2. Output directories
# ==============================================================================

OUTDIR="${BASE}/Healthy_TI_Rectum_N278_genotype"

VCF_DIR="${OUTDIR}/01_vcf_annotated"
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
TMP_DIR="${VCF_DIR}/tmp"

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
# 3. Parameters
# ==============================================================================

R2_CUTOFF=0.30
MAF_CUTOFF=0.01
GENO_CUTOFF=0.01
HWE_CUTOFF=1e-6

LD_WINDOW=200
LD_STEP=50
LD_R2=0.2

N_PCS=20


# ==============================================================================
# 4. Preflight
# ==============================================================================

echo "============================================================"
echo "PRE-FLIGHT"
echo "============================================================"
echo "Parallel chromosomes : ${N_PARALLEL}"
echo "Threads/chromosome   : ${THREADS_PER_CHR}"
echo "PLINK genotype field : GT"
echo

[[ -s "${VCF}" ]] || {
    echo "ERROR: missing VCF: ${VCF}" >&2
    exit 1
}

[[ -s "${SAMPLE_FILE}" ]] || {
    echo "ERROR: missing sample file: ${SAMPLE_FILE}" >&2
    exit 1
}

[[ -d "${DBSNP_DIR}" ]] || {
    echo "ERROR: missing dbSNP directory: ${DBSNP_DIR}" >&2
    exit 1
}


# ------------------------------------------------------------------------------
# 4A. Input genotype VCF index
# Existing .csi or .tbi is accepted.
# ------------------------------------------------------------------------------

if [[ ! -s "${VCF}.tbi" && ! -s "${VCF}.csi" ]]; then
    echo "Input VCF index missing; creating CSI index..."
    "${BCFTOOLS}" index -c -f "${VCF}"
fi

"${BCFTOOLS}" quickcheck -v "${VCF}"


# ------------------------------------------------------------------------------
# 4B. Selected sample list
# ------------------------------------------------------------------------------

KEEP_SAMPLES="${TMP_DIR}/keep_samples.txt"

awk '
    NR == 1 && $1 == "individual_id" { next }
    NF > 0 { print $1 }
' "${SAMPLE_FILE}" \
| sort -u \
> "${KEEP_SAMPLES}"

N_REQUESTED="$(wc -l < "${KEEP_SAMPLES}")"

echo "Requested samples: ${N_REQUESTED}"

if [[ "${N_REQUESTED}" -ne 278 ]]; then
    echo "WARNING: expected 278 samples, found ${N_REQUESTED}" >&2
fi


# ------------------------------------------------------------------------------
# 4C. Check selected samples exist in genotype VCF
# ------------------------------------------------------------------------------

VCF_SAMPLES="${TMP_DIR}/vcf_samples.txt"
MISSING_SAMPLES="${QC_DIR}/missing_requested_samples.txt"

"${BCFTOOLS}" query -l "${VCF}" \
| sort -u \
> "${VCF_SAMPLES}"

comm -23 \
    "${KEEP_SAMPLES}" \
    "${VCF_SAMPLES}" \
> "${MISSING_SAMPLES}"

N_MISSING="$(wc -l < "${MISSING_SAMPLES}")"

echo "Requested samples missing from VCF: ${N_MISSING}"

if [[ "${N_MISSING}" -gt 0 ]]; then
    cat "${MISSING_SAMPLES}" >&2
    exit 1
fi


# ------------------------------------------------------------------------------
# 4D. Check chromosome-specific dbSNP VCF + TBI
# ------------------------------------------------------------------------------

echo "Checking chromosome-specific dbSNP files..."

for i in $(seq 1 22); do
    DBSNP="${DBSNP_DIR}/GCF_000001405.40.chr${i}.vcf.gz"

    [[ -s "${DBSNP}" ]] || {
        echo "ERROR: missing dbSNP chr${i}: ${DBSNP}" >&2
        exit 1
    }

    if [[ ! -s "${DBSNP}.tbi" ]]; then
        echo "Creating TBI index for dbSNP chr${i}..."
        "${BCFTOOLS}" index -t -f "${DBSNP}"
    fi
done

echo "dbSNP chr1-22 + TBI indexes: OK"


# ------------------------------------------------------------------------------
# 4E. Chromosome rename map: chr1 -> 1
# ------------------------------------------------------------------------------

CHR_RENAME_FILE="${TMP_DIR}/rename_chr_to_num.txt"

: > "${CHR_RENAME_FILE}"

for i in $(seq 1 22); do
    printf "chr%s\t%s\n" "${i}" "${i}" >> "${CHR_RENAME_FILE}"
done


# ------------------------------------------------------------------------------
# 4F. Detect imputation R2 field
# ------------------------------------------------------------------------------

HEADER="$("${BCFTOOLS}" view -h "${VCF}")"
R2_FIELD=""

for TAG in R2 DR2 INFO_SCORE INFO; do
    if grep -q "^##INFO=<ID=${TAG}," <<< "${HEADER}"; then
        R2_FIELD="${TAG}"
        break
    fi
done

if [[ -z "${R2_FIELD}" ]]; then
    echo "ERROR: cannot detect R2/DR2/INFO_SCORE/INFO field." >&2
    grep '^##INFO=<ID=' <<< "${HEADER}" | head -50 >&2
    exit 1
fi

R2_EXPR="INFO/${R2_FIELD}>${R2_CUTOFF}"

echo "R2 field : ${R2_FIELD}"
echo "Filter   : ${R2_EXPR}"


# ==============================================================================
# 5. Per-chromosome processing function: STEP 1-12
# ==============================================================================

process_chr() {
    set -euo pipefail

    local i="$1"
    local INPUT_CHR="chr${i}"
    local PREFIX="chr${i}.Healthy_TI_rectum"

    local DBSNP="${DBSNP_DIR}/GCF_000001405.40.chr${i}.vcf.gz"

    local TMP_VCF="${TMP_DIR}/${PREFIX}.R03.biallelicSNP.renamed.vcf.gz"
    local ANNOT_VCF="${VCF_DIR}/${PREFIX}.R03.biallelicSNP.rs.vcf.gz"

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
    local N_ANNOT
    local N_RAW
    local N_DEDUP
    local N_QC
    local N_RSID
    local N_LD

    echo "===== chr${i} START ====="


    # ==========================================================================
    # STEP 1
    # Selected samples + R2 + biallelic SNP + chrN -> N
    # ==========================================================================

    "${BCFTOOLS}" view \
        -r "${INPUT_CHR}" \
        -S "${KEEP_SAMPLES}" \
        -i "${R2_EXPR}" \
        -m2 \
        -M2 \
        -v snps \
        -Ou \
        "${VCF}" \
    | "${BCFTOOLS}" annotate \
        --rename-chrs "${CHR_RENAME_FILE}" \
        -Oz \
        -o "${TMP_VCF}"

    echo "chr${i}: STEP1 extraction + R2 + biallelic + rename OK"


    # ==========================================================================
    # STEP 2
    # Temporary VCF TBI index
    # ==========================================================================

    "${BCFTOOLS}" index \
        -t \
        -f \
        "${TMP_VCF}"

    echo "chr${i}: STEP2 temporary TBI index OK"


    # ==========================================================================
    # STEP 3
    # dbSNP rsID annotation + final TBI
    # ==========================================================================

    "${BCFTOOLS}" annotate \
        -a "${DBSNP}" \
        -c ID \
        -Oz \
        -o "${ANNOT_VCF}" \
        "${TMP_VCF}"

    "${BCFTOOLS}" index \
        -t \
        -f \
        "${ANNOT_VCF}"

    "${BCFTOOLS}" quickcheck -v "${ANNOT_VCF}"

    echo "chr${i}: STEP3 dbSNP annotation + TBI OK"


    # ==========================================================================
    # STEP 4
    # Annotated VCF -> PGEN using GT (NOT DS)
    # ==========================================================================

    "${PLINK2}" \
        --vcf "${ANNOT_VCF}" \
        --set-missing-var-ids '@:#:$r:$a' \
        --make-pgen \
        --threads "${THREADS_PER_CHR}" \
        --out "${RAW_PREFIX}"

    echo "chr${i}: STEP4 VCF -> PGEN (GT) OK"


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
    # Remove duplicate IDs + IDs beginning with "chr"
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
    # ALT-coded .traw
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
    # TPED/TFAM directly with PLINK2
    # ==========================================================================

    "${PLINK2}" \
        --pfile "${RSID_PREFIX}" \
        --export tped \
        --threads "${THREADS_PER_CHR}" \
        --out "${TPED_PREFIX}"

    echo "chr${i}: STEP9C TPED/TFAM export OK"


    # ==========================================================================
    # STEP 9D
    # Verify .traw COUNTED allele == .pvar ALT
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
        echo "ERROR chr${i}: ${N_TRAW_MISMATCH} .traw ALT-coding mismatches." >&2
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

    N_ANNOT="$("${BCFTOOLS}" index -n "${ANNOT_VCF}")"
    N_RAW="$(awk '!/^#/ {n++} END{print n+0}' "${RAW_PREFIX}.pvar")"
    N_DEDUP="$(awk '!/^#/ {n++} END{print n+0}' "${DEDUP_PREFIX}.pvar")"
    N_QC="$(awk '!/^#/ {n++} END{print n+0}' "${QC_PREFIX}.pvar")"
    N_RSID="$(awk '!/^#/ {n++} END{print n+0}' "${RSID_PREFIX}.pvar")"
    N_LD="$(wc -l < "${LD_PREFIX}.prune.in")"

    {
        echo "chromosome=${i}"
        echo "genotype_source=GT"
        echo "r2_filter=${R2_EXPR}"
        echo "dbsnp=${DBSNP}"
        echo "annotated_vcf=${ANNOT_VCF}"
        echo "n_variants_annotated=${N_ANNOT}"
        echo "n_variants_raw_pgen=${N_RAW}"
        echo "duplicate_ID_groups=${N_DUP_IDS}"
        echo "chr_ID_variants_removed=${N_CHR_IDS}"
        echo "total_preQC_excluded=${N_EXCLUDE}"
        echo "n_variants_after_preQC_cleanup=${N_DEDUP}"
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

    # Remove temporary renamed VCF only after the chromosome fully succeeds.
    rm -f \
        "${TMP_VCF}" \
        "${TMP_VCF}.tbi"

    echo "===== chr${i} DONE ====="
}


# Export variables + function for GNU parallel's bash children.
export \
    BCFTOOLS \
    PLINK2 \
    VCF \
    DBSNP_DIR \
    KEEP_SAMPLES \
    CHR_RENAME_FILE \
    VCF_DIR \
    RAW_PGEN_DIR \
    DEDUP_DIR \
    QC_PGEN_DIR \
    RSID_DIR \
    TRAW_DIR \
    LD_DIR \
    PCA_DIR \
    QC_DIR \
    TMP_DIR \
    R2_EXPR \
    MAF_CUTOFF \
    GENO_CUTOFF \
    HWE_CUTOFF \
    LD_WINDOW \
    LD_STEP \
    LD_R2 \
    THREADS_PER_CHR

export -f process_chr


# ==============================================================================
# 6. Run STEP 1-12 for chr1-22
# ==============================================================================

echo
echo "============================================================"
echo "RUN chr1-22"
echo "============================================================"

seq 1 22 \
| parallel \
    -j "${N_PARALLEL}" \
    --line-buffer \
    --halt soon,fail=1 \
    --joblog "${LOG_DIR}/parallel.joblog" \
    bash -c 'process_chr "$1"' _ {}


# ==============================================================================
# STEP 13
# Merge chr1-22 LD-pruned PGEN + genome-wide PCA
# ==============================================================================

echo
echo "============================================================"
echo "STEP 13: MERGE chr1-22 LD-PRUNED PGEN + PCA"
echo "============================================================"

PMERGE_LIST="${JOINT_PCA_DIR}/chr2_22.pmerge_list.txt"
MERGED_PREFIX="${JOINT_PCA_DIR}/Healthy_TI_rectum_N278.chr1_22.LDpruned"
PCA_OUT="${JOINT_PCA_DIR}/Healthy_TI_rectum_N278.chr1_22.PCA"

for i in $(seq 1 22); do
    P="${PCA_DIR}/chr${i}.Healthy_TI_rectum.qc.rsid.ldpruned"

    for ext in pgen pvar psam; do
        if [[ ! -s "${P}.${ext}" ]]; then
            echo "ERROR: missing ${P}.${ext}" >&2
            exit 1
        fi
    done
done

: > "${PMERGE_LIST}"

for i in $(seq 2 22); do
    echo "${PCA_DIR}/chr${i}.Healthy_TI_rectum.qc.rsid.ldpruned" \
        >> "${PMERGE_LIST}"
done

"${PLINK2}" \
    --pfile "${PCA_DIR}/chr1.Healthy_TI_rectum.qc.rsid.ldpruned" \
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

N_VCF="$(find "${VCF_DIR}" -maxdepth 1 -type f -name "chr*.rs.vcf.gz" | wc -l)"
N_TRAW="$(find "${TRAW_DIR}" -maxdepth 1 -type f -name "*.traw" | wc -l)"
N_TPED="$(find "${TRAW_DIR}" -maxdepth 1 -type f -name "*.tped" | wc -l)"
N_TFAM="$(find "${TRAW_DIR}" -maxdepth 1 -type f -name "*.tfam" | wc -l)"
N_SUMMARY="$(find "${QC_DIR}" -maxdepth 1 -type f -name "*.summary.txt" | wc -l)"

echo "Annotated VCFs : ${N_VCF}/22"
echo ".traw files    : ${N_TRAW}/22"
echo ".tped files    : ${N_TPED}/22"
echo ".tfam files    : ${N_TFAM}/22"
echo "summary files  : ${N_SUMMARY}/22"

echo
echo "Output root : ${OUTDIR}"
echo "Annotated   : ${VCF_DIR}"
echo "rsID PGEN   : ${RSID_DIR}"
echo "TRAW/TPED   : ${TRAW_DIR}"
echo "LD prune    : ${LD_DIR}"
echo "PCA chr     : ${PCA_DIR}"
echo "PCA merged  : ${JOINT_PCA_DIR}"
echo "Eigenvec    : ${PCA_OUT}.eigenvec"
echo "Eigenval    : ${PCA_OUT}.eigenval"
echo "QC          : ${QC_DIR}"
echo "Job log     : ${LOG_DIR}/parallel.joblog"
