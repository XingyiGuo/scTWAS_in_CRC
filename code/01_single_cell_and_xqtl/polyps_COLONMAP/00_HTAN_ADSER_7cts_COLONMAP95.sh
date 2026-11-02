# ----------------------------------------------------------------------
# STEP 1: match only the remaining 7 CTs to COLONMAP_95 genotype/PCA.
# ----------------------------------------------------------------------

cd /data/l2_bioinfo1/liq17/scTWAS/data/VUMC_HTAN/src
EXPRESSION_ROOT="/v5000/data/l2_bioinfo1/liq17/scTWAS/data/VUMC_HTAN/expression/ADSER_pseudobulk_expression_for_peer"
GENOTYPE_PSAM="/v5000/data/l2_bioinfo1/liq17/scTWAS/data/VUMC_HTAN/genotype/plinkformat/COLONMAP_95.chr1.nodup.maf001.rs.psam"
GENOTYPE_PCA="/data/l2_bioinfo1/liq17/scTWAS/data/VUMC_HTAN_GSE266546/genotype/VUMC_HTAN/plinkformat/genotype_pca/COLONMAP_95.chr1_22.nodup.maf001.ld01.rs.pca.eigenvec"
#for ct in ABS EE GOB CT TUF TAC STM; do
for ct in CT; do
    echo "============================================================"
    echo "Matching ${ct}"
    echo "============================================================"

    python 01_match_HTAN_ASC_SSC_expression_genotype_HTANdirect.py \
        "${ct}" \
        --expression-root "${EXPRESSION_ROOT}" \
        --count-matrix "${EXPRESSION_ROOT}/${ct}/count_matrix_batch12_${ct}.chr1_22.csv" \
        --classifications AD SER \
        --genotype-psam "${GENOTYPE_PSAM}" \
        --pca-eigenvec "${GENOTYPE_PCA}"
done


# ----------------------------------------------------------------------
# STEP 2: expression preprocessing, one ct per process.
# ----------------------------------------------------------------------
cd /data/l2_bioinfo1/liq17/scTWAS/data/VUMC_HTAN/src
SCRIPT_DIR="/data/l2_bioinfo1/liq17/scTWAS/data/VUMC_HTAN/src"
CONFIG="/v5000/data/l2_bioinfo1/liq17/scTWAS/data/VUMC_HTAN/src/HTAN_ADSER_9cts.eqtl.config_COLONMAP95_PCAfixed_pathfix_PEERlocal.txt"
mkdir -p logs_expression
#for ct in ABS EE GOB CT TUF TAC STM; do
for ct in CT; do
    echo "[START] ${ct}"
    (
        python "${SCRIPT_DIR}/02_prepare_pseudobulk_expression_ASC_SSC_bashconfig.py" \
            --config "${CONFIG}" \
            --cell-types "${ct}"
    ) > "logs_expression/${ct}.log" 2>&1 &
done
wait
echo "All 7 expression preprocessing jobs finished."



# ----------------------------------------------------------------------
# STEP 3: Prepare for PEER, one CT per process.
# ----------------------------------------------------------------------
cd /data/l2_bioinfo1/liq17/scTWAS/data/VUMC_HTAN/src

SCRIPT_DIR="/data/l2_bioinfo1/liq17/scTWAS/data/VUMC_HTAN/src"
EXPRESSION_ROOT="/v5000/data/l2_bioinfo1/liq17/scTWAS/data/VUMC_HTAN/expression/ADSER_pseudobulk_expression_for_peer"
mkdir -p logs_peer
for ct in ABS EE GOB CT TUF TAC STM; do
    echo "[START] PEER ${ct}"
    (
        python "${SCRIPT_DIR}/03_Process_peer_HTAN_ASC_SSC_HTANdirect_geneSymbol.py" \
            "${ct}" \
            --expression-root "${EXPRESSION_ROOT}"
    ) > "logs_peer/${ct}.log" 2>&1 &
done
wait
echo "============================================================"
echo "All 7 PEER jobs finished."
echo "============================================================"


# ----------------------------------------------------------------------
# STEP 4: shared tensorQTL genotype.
# ASC/SSC completion usually means this already exists; verify all chr1-22.
# Prepare it only when at least one chromosome is missing.
# ----------------------------------------------------------------------
cd /data/l2_bioinfo1/liq17/scTWAS/data/VUMC_HTAN/src
SCRIPT_DIR="/data/l2_bioinfo1/liq17/scTWAS/data/VUMC_HTAN/src"
CONFIG="/v5000/data/l2_bioinfo1/liq17/scTWAS/data/VUMC_HTAN/src/HTAN_ADSER_9cts.eqtl.config_COLONMAP95_PCAfixed_pathfix_PEERlocal.txt"
LOG_ROOT="${SCRIPT_DIR}/logs_tensorqtl_genotype_$(date +%Y%m%d_%H%M%S)"
mkdir -p "${LOG_ROOT}"
# Load configuration
source "${CONFIG}"
# Expand config templates
expand_template() {
    local x="$1"
    local ct="${2:-}"
    local chrom="${3:-}"
    x="${x//\{work_dir\}/${WORK_DIR}}"
    x="${x//\{dataset_prefix\}/${DATASET_PREFIX}}"
    x="${x//\{tissue\}/${TISSUE}}"

    [[ -z "${ct}" ]] || x="${x//\{CT\}/${ct}}"
    [[ -z "${chrom}" ]] || x="${x//\{chrom\}/${chrom}}"

    printf '%s\n' "${x}"
}

echo
echo "================ 04_tensorqtl_genotype_check ================"
n_missing=0
for CHROM in $(seq 1 22); do
    pfx="$(expand_template \
        "${TENSORQTL_PGEN_PREFIX_TEMPLATE}" \
        "" \
        "${CHROM}")"
    echo "Checking chr${CHROM}: ${pfx}"
    for ext in pgen pvar psam; do
        if [[ ! -s "${pfx}.${ext}" ]]; then
            echo "  [MISSING] ${pfx}.${ext}"
            n_missing=$((n_missing + 1))
        fi
    done
done


echo
echo "Missing tensorQTL genotype files: ${n_missing}"
echo

if (( n_missing == 0 )); then
    echo "[SKIP] tensorQTL COLONMAP_95 genotype already complete for chr1-22."
else
    echo "[RUN] ${n_missing} tensorQTL genotype file(s) missing."
    echo "[RUN] Preparing tensorQTL genotype..."
    bash "${SCRIPT_DIR}/prepare_tensorqtl_genotype.sh" "${CONFIG}" \
        > "${LOG_ROOT}/04_tensorqtl_genotype.log" 2>&1
    status=$?
    if [[ ${status} -eq 0 ]]; then
        echo "[DONE] tensorQTL genotype preparation"
        echo "Log: ${LOG_ROOT}/04_tensorqtl_genotype.log"
    else
        echo "[ERROR] tensorQTL genotype preparation failed"
        echo "See:"
        echo "  ${LOG_ROOT}/04_tensorqtl_genotype.log"
    fi

fi


# ----------------------------------------------------------------------
# STEP 5: PEER residual -> tensorQTL BED for remaining 7 only.
# ----------------------------------------------------------------------
cd /data/l2_bioinfo1/liq17/scTWAS/data/VUMC_HTAN/src
SCRIPT_DIR="/data/l2_bioinfo1/liq17/scTWAS/data/VUMC_HTAN/src"
CONFIG="/v5000/data/l2_bioinfo1/liq17/scTWAS/data/VUMC_HTAN/src/HTAN_ADSER_9cts.eqtl.config_COLONMAP95_PCAfixed_pathfix_PEERlocal.txt"
source "${CONFIG}"
mkdir -p logs_bed

expand_template() {
    local x="$1"
    local ct="${2:-}"
    x="${x//\{work_dir\}/${WORK_DIR}}"
    x="${x//\{dataset_prefix\}/${DATASET_PREFIX}}"
    x="${x//\{tissue\}/${TISSUE}}"
    [[ -z "${ct}" ]] || x="${x//\{CT\}/${ct}}"
    printf '%s\n' "${x}"
}

step5_bed() {
    local CT="$1"
    local expression_file
    local bed_file
    local sample_map
    expression_file="$(expand_template "${EXPRESSION_TEMPLATE}" "${CT}")"
    bed_file="$(expand_template "${PHENOTYPE_BED_TEMPLATE}" "${CT}")"
    sample_map="$(expand_template "${SAMPLE_MAPPING_FILE}" "${CT}")"
    echo "============================================================"
    echo "Cell type  : ${CT}"
    echo "Expression : ${expression_file}"
    echo "Mapping    : ${sample_map}"
    echo "BED output : ${bed_file}"
    echo "============================================================"
    [[ -s "${expression_file}" ]] || {
        echo "ERROR: missing PEER residual expression:"
        echo "${expression_file}"
        return 1
    }
    [[ -s "${sample_map}" ]] || {
        echo "ERROR: missing sample mapping:"
        echo "${sample_map}"
        return 1
    }
    [[ -s "${PEER_TO_BED_SCRIPT}" ]] || {
        echo "ERROR: PEER-to-BED script missing:"
        echo "${PEER_TO_BED_SCRIPT}"
        return 1
    }
    mkdir -p "$(dirname "${bed_file}")"
    "${PYTHON_BIN}" "${PEER_TO_BED_SCRIPT}" \
        --expression "${expression_file}" \
        --mapping "${sample_map}" \
        --gtf "${GENE_ANNOTATION}" \
        --output "${bed_file}" \
        --audit-prefix "${bed_file%.bed.gz}" \
        --expression-id-col "${EXPRESSION_ID_COL}" \
        --genotype-id-col "${GENOTYPE_ID_COL}" \
        --gene-id-style "${BED_GENE_ID_STYLE}"
}
for CT in ABS EE GOB CT TUF TAC STM; do
    echo "[START] BED ${CT}"
    (
        step5_bed "${CT}"
    ) > "logs_bed/${CT}.log" 2>&1 &

done
wait
echo
echo "============================================================"
echo "All 7 tensorQTL BED jobs finished."
echo "============================================================"


# ----------------------------------------------------------------------
# STEP 6: xQTL.
# Seven CT workers run in parallel; each worker runs chr1..22 sequentially.
# xQTL script already skips existing parquet outputs when configured to do so.
# ----------------------------------------------------------------------
cd /data/l2_bioinfo1/liq17/scTWAS/data/VUMC_HTAN/src

SCRIPT_DIR="/data/l2_bioinfo1/liq17/scTWAS/data/VUMC_HTAN/src"

CONFIG="/v5000/data/l2_bioinfo1/liq17/scTWAS/data/VUMC_HTAN/src/HTAN_ADSER_9cts.eqtl.config_COLONMAP95_PCAfixed_pathfix_PEERlocal.txt"

LOG_ROOT="${SCRIPT_DIR}/logs_xqtl_$(date +%Y%m%d_%H%M%S)"
mkdir -p "${LOG_ROOT}"


step6_xqtl_ct() {

    local CT="$1"
    local CHROM

    for CHROM in $(seq 1 22); do

        echo "======================================================================"
        echo "${CT}: chromosome ${CHROM}"
        echo "======================================================================"

        bash "${SCRIPT_DIR}/06_xQTL_per_chrom_ASC_SSC.sh" \
            "${CT}" \
            "${CHROM}" \
            "${CONFIG}"

    done
}


#for CT in ABS EE GOB CT TUF TAC STM; do
for CT in ASC SSC; do
    echo "[START] tensorQTL ${CT}"

    (
        step6_xqtl_ct "${CT}"
    ) > "${LOG_ROOT}/${CT}.log" 2>&1 &

done

wait

echo
echo "======================================================================"
echo "All 7 tensorQTL cell-type workers finished."
echo "Logs: ${LOG_ROOT}"
echo "======================================================================"