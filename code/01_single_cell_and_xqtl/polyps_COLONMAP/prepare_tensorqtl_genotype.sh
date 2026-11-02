#!/usr/bin/env bash
set -euo pipefail

# ==============================================================================
# prepare_tensorqtl_genotype.sh
#
# Create tensorQTL-specific PGEN/PVAR/PSAM files under:
#   05_pgen_rsid/geno_for_tensorqtl/
#
# Original variant ID:
#   rs562457791
#
# New variant ID:
#   rs562457791;22_15331114_C_T
#
# where:
#   C = REF
#   T = ALT
#
# tensorQTL's reported "af" is ALT allele frequency, and the nominal slope
# corresponds to the genotype dosage used by tensorQTL. Thus the ALT allele
# encoded at the end of the new ID is the eQTL effect allele.
#
# Usage:
#   bash prepare_tensorqtl_genotype.sh <eqtl.config.txt>
# ==============================================================================

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
[[ $# -ge 1 ]] || {
    echo "Usage: prepare_tensorqtl_genotype.sh <eqtl.config.txt>" >&2
    exit 2
}
CONFIG="$1"

[[ -s "${CONFIG}" ]] || {
    echo "ERROR: config not found: ${CONFIG}" >&2
    exit 1
}

# shellcheck disable=SC1090
source "${CONFIG}"

[[ -x "${PLINK2}" ]] || {
    echo "ERROR: PLINK2 not executable: ${PLINK2}" >&2
    exit 1
}

expand_template() {
    local x="$1"
    local chrom="$2"
    x="${x//\{work_dir\}/${WORK_DIR}}"
    x="${x//\{dataset_prefix\}/${DATASET_PREFIX}}"
    x="${x//\{tissue\}/${TISSUE}}"
    x="${x//\{chrom\}/${chrom}}"
    printf '%s\n' "${x}"
}

OUT_DIR="${TENSORQTL_GENO_DIR//\{work_dir\}/${WORK_DIR}}"
mkdir -p "${OUT_DIR}"

for CHROM in $(seq 1 22); do
    IN_PREFIX="$(expand_template "${ORIGINAL_PGEN_PREFIX_TEMPLATE}" "${CHROM}")"
    OUT_PREFIX="$(expand_template "${TENSORQTL_PGEN_PREFIX_TEMPLATE}" "${CHROM}")"

    for ext in pgen pvar psam; do
        [[ -s "${IN_PREFIX}.${ext}" ]] || {
            echo "ERROR: missing ${IN_PREFIX}.${ext}" >&2
            exit 1
        }
    done

    UPDATE_FILE="${OUT_DIR}/chr${CHROM}.variant_update_name.tsv"
    KEY_FILE="${OUT_DIR}/chr${CHROM}.variant_key.tsv"

    echo
    echo "======================================================================"
    echo "chr${CHROM}"
    echo "Input : ${IN_PREFIX}"
    echo "Output: ${OUT_PREFIX}"
    echo "======================================================================"

    # PLINK2 --update-name default format:
    #   column 1 = old ID
    #   column 2 = new ID
    #
    # New ID:
    #   rsID;CHROM_POS_REF_ALT
    awk '
        BEGIN {FS=OFS="\t"}
        !/^#/ {
            old=$3
            new=old ";" $1 "_" $2 "_" $4 "_" $5
            print old,new
        }
    ' "${IN_PREFIX}.pvar" > "${UPDATE_FILE}"

    {
        printf 'old_id\tnew_id\tchrom\tpos\tref\talt\teffect_allele\n'
        awk '
            BEGIN {FS=OFS="\t"}
            !/^#/ {
                old=$3
                new=old ";" $1 "_" $2 "_" $4 "_" $5
                print old,new,$1,$2,$4,$5,$5
            }
        ' "${IN_PREFIX}.pvar"
    } > "${KEY_FILE}"

    N_INPUT="$(awk '!/^#/ {n++} END{print n+0}' "${IN_PREFIX}.pvar")"
    N_MAP="$(wc -l < "${UPDATE_FILE}")"
    N_DUP_NEW="$(
        cut -f2 "${UPDATE_FILE}" |
        sort |
        uniq -d |
        wc -l
    )"

    [[ "${N_INPUT}" -eq "${N_MAP}" ]] || {
        echo "ERROR: variant mapping count mismatch on chr${CHROM}" >&2
        exit 1
    }

    [[ "${N_DUP_NEW}" -eq 0 ]] || {
        echo "ERROR: duplicated new variant IDs on chr${CHROM}" >&2
        exit 1
    }

    "${PLINK2}" \
        --pfile "${IN_PREFIX}" \
        --update-name "${UPDATE_FILE}" \
        --output-chr chrM \
        --make-pgen \
        --threads "${SLURM_CPUS_PER_TASK:-4}" \
        --out "${OUT_PREFIX}"

    N_OUTPUT="$(awk '!/^#/ {n++} END{print n+0}' "${OUT_PREFIX}.pvar")"

    [[ "${N_OUTPUT}" -eq "${N_INPUT}" ]] || {
        echo "ERROR: output variant count changed on chr${CHROM}: ${N_INPUT} -> ${N_OUTPUT}" >&2
        exit 1
    }

    # Confirm every data-line ID follows:
    #   oldrs;chrom_pos_ref_alt
    BAD="$(
        awk '
            BEGIN {FS="\t"; bad=0}
            !/^#/ {
                if ($3 !~ /;[^_]+_[0-9]+_[^_]+_[^_]+$/) bad++
            }
            END {print bad+0}
        ' "${OUT_PREFIX}.pvar"
    )"

    [[ "${BAD}" -eq 0 ]] || {
        echo "ERROR: ${BAD} malformed tensorQTL variant IDs on chr${CHROM}" >&2
        exit 1
    }

    echo "Variants: ${N_OUTPUT}"
    echo "Example : $(awk '!/^#/ {print $3; exit}' "${OUT_PREFIX}.pvar")"
done

echo
echo "TensorQTL genotype preparation completed:"
echo "  ${OUT_DIR}"
