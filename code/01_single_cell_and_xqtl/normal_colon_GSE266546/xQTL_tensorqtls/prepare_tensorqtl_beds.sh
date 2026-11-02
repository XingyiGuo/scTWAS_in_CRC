#!/usr/bin/env bash
set -euo pipefail

# ==============================================================================
# prepare_tensorqtl_beds.sh
#
# Generic PEER-residual -> tensorQTL BED preparation.
#
# Dataset-specific paths and sample-ID columns come entirely from the config.
#
# Usage:
#   bash prepare_tensorqtl_beds.sh <eqtl.config.txt>
#
# Examples:
#   bash prepare_tensorqtl_beds.sh HTAN_GSE266546.eqtl.config.txt
#   bash prepare_tensorqtl_beds.sh WellcomeSanger.eqtl.config.txt
# ==============================================================================

[[ $# -ge 1 ]] || {
    echo "Usage: prepare_tensorqtl_beds.sh <eqtl.config.txt>" >&2
    exit 2
}

CONFIG="$1"

[[ -s "${CONFIG}" ]] || {
    echo "ERROR: config not found: ${CONFIG}" >&2
    exit 1
}

# shellcheck disable=SC1090
source "${CONFIG}"

expand_template() {
    local x="$1"
    local ct="$2"

    x="${x//\{work_dir\}/${WORK_DIR}}"
    x="${x//\{dataset_prefix\}/${DATASET_PREFIX}}"
    x="${x//\{tissue\}/${TISSUE}}"
    x="${x//\{CT\}/${ct}}"

    printf '%s\n' "${x}"
}

SAMPLE_MAP="$(expand_template "${SAMPLE_MAPPING_FILE}" NA)"

[[ -n "${PEER_TO_BED_SCRIPT:-}" ]] || {
    echo "ERROR: PEER_TO_BED_SCRIPT is not defined in ${CONFIG}" >&2
    exit 1
}

[[ -s "${PEER_TO_BED_SCRIPT}" ]] || {
    echo "ERROR: PEER-to-BED converter not found: ${PEER_TO_BED_SCRIPT}" >&2
    exit 1
}

[[ -s "${SAMPLE_MAP}" ]] || {
    echo "ERROR: sample mapping file not found: ${SAMPLE_MAP}" >&2
    exit 1
}

[[ -s "${GENE_ANNOTATION}" ]] || {
    echo "ERROR: gene annotation not found: ${GENE_ANNOTATION}" >&2
    exit 1
}

echo "Config             : ${CONFIG}"
echo "Dataset            : ${DATASET_PREFIX}"
echo "Tissue             : ${TISSUE}"
echo "PEER-to-BED script : ${PEER_TO_BED_SCRIPT}"
echo "Sample map         : ${SAMPLE_MAP}"
echo "Expression ID col  : ${EXPRESSION_ID_COL}"
echo "Genotype ID col    : ${GENOTYPE_ID_COL}"

for CT in "${CELL_TYPES[@]}"; do
    EXPRESSION_FILE="$(expand_template "${EXPRESSION_TEMPLATE}" "${CT}")"
    BED_FILE="$(expand_template "${PHENOTYPE_BED_TEMPLATE}" "${CT}")"

    echo
    echo "======================================================================"
    echo "Cell type : ${CT}"
    echo "Expression: ${EXPRESSION_FILE}"
    echo "BED       : ${BED_FILE}"
    echo "======================================================================"

    [[ -s "${EXPRESSION_FILE}" ]] || {
        echo "ERROR: missing expression: ${EXPRESSION_FILE}" >&2
        exit 1
    }

    "${PYTHON_BIN}" "${PEER_TO_BED_SCRIPT}" \
        --expression "${EXPRESSION_FILE}" \
        --mapping "${SAMPLE_MAP}" \
        --gtf "${GENE_ANNOTATION}" \
        --output "${BED_FILE}" \
        --audit-prefix "${BED_FILE%.bed.gz}" \
        --expression-id-col "${EXPRESSION_ID_COL}" \
        --genotype-id-col "${GENOTYPE_ID_COL}" \
        --gene-id-style "${BED_GENE_ID_STYLE}"
done

echo
echo "All tensorQTL BED files completed."
