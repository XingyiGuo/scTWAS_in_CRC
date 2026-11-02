#!/usr/bin/env bash
#SBATCH --job-name=xQTL_cis
#SBATCH --error=%x-%j.error
#SBATCH --output=%x-%j.out
#SBATCH --mem=40G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --time=7-00:00:00
#SBATCH --account=l2_bioinfo1

set -euo pipefail

# ==============================================================================
# conda activate py3_10 
# xQTL_per_chrom.sh
#
# Run tensorQTL cis_nominal for one cell type and one chromosome.
#
# Required preprocessing:
#   bash prepare_tensorqtl_beds.sh
#   bash prepare_tensorqtl_genotype.sh
#
# Usage:
#   sbatch xQTL_per_chrom.sh ABS 1
#   sbatch xQTL_per_chrom.sh 1 1
#   sbatch --array=1-22 xQTL_per_chrom.sh ABS
# ==============================================================================

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

[[ $# -ge 1 ]] || {
    echo "Usage: xQTL_per_chrom.sh <CT_index_or_name> [chrom] [config]" >&2
    exit 2
}

CT_ARG="$1"
CHROM="${2:-${SLURM_ARRAY_TASK_ID:-}}"
#DEFAULT_CONFIG="/data/l2_bioinfo1/liq17/scTWAS/data/VUMC_HTAN_GSE266546/src/xQTL_tensorqtls/WellcomeSanger.eqtl.config.txt"
DEFAULT_CONFIG="/data/l2_bioinfo1/liq17/scTWAS/data/VUMC_HTAN_GSE266546/src/xQTL_tensorqtls/HTAN_GSE266546.eqtl.config.txt"
CONFIG="${3:-${DEFAULT_CONFIG}}"

[[ -n "${CHROM}" ]] || {
    echo "ERROR: chromosome is required." >&2
    exit 2
}

[[ -s "${CONFIG}" ]] || {
    echo "ERROR: config not found: ${CONFIG}" >&2
    exit 1
}

# shellcheck disable=SC1090
source "${CONFIG}"

if ! [[ "${CHROM}" =~ ^([1-9]|1[0-9]|2[0-2])$ ]]; then
    echo "ERROR: chromosome must be 1-22." >&2
    exit 1
fi

# Cell type can be literal or TWAS-style 1-based index.
if [[ "${CT_ARG}" =~ ^[0-9]+$ ]]; then
    IDX="${CT_ARG}"
    (( IDX >= 1 && IDX <= ${#CELL_TYPES[@]} )) || {
        echo "ERROR: CT index out of range." >&2
        exit 1
    }
    CT="${CELL_TYPES[$((IDX - 1))]}"
else
    CT="${CT_ARG}"
    printf '%s\n' "${CELL_TYPES[@]}" | grep -Fxq "${CT}" || {
        echo "ERROR: unknown cell type ${CT}" >&2
        exit 1
    }
fi

expand_template() {
    local x="$1"
    x="${x//\{work_dir\}/${WORK_DIR}}"
    x="${x//\{dataset_prefix\}/${DATASET_PREFIX}}"
    x="${x//\{tissue\}/${TISSUE}}"
    x="${x//\{CT\}/${CT}}"
    x="${x//\{chrom\}/${CHROM}}"
    printf '%s\n' "${x}"
}

BED_ALL="$(expand_template "${PHENOTYPE_BED_TEMPLATE}")"
PGEN_PREFIX="$(expand_template "${TENSORQTL_PGEN_PREFIX_TEMPLATE}")"
OUTPUT_DIR="$(expand_template "${OUTPUT_DIR_TEMPLATE}")"
PREFIX="$(expand_template "${PREFIX_TEMPLATE}")"

CHROM_BED_DIR="${OUTPUT_DIR}/${CHROM_BED_SUBDIR}"
AUDIT_DIR="${OUTPUT_DIR}/${AUDIT_SUBDIR}"
CHROM_BED="${CHROM_BED_DIR}/${PREFIX}.chr${CHROM}.tensorqtl.bed.gz"

mkdir -p "${OUTPUT_DIR}" "${CHROM_BED_DIR}" "${AUDIT_DIR}"

for f in "${BED_ALL}" \
         "${PGEN_PREFIX}.pgen" \
         "${PGEN_PREFIX}.pvar" \
         "${PGEN_PREFIX}.psam"
do
    [[ -s "${f}" ]] || {
        echo "ERROR: required input missing: ${f}" >&2
        echo "Run prepare_tensorqtl_beds.sh and prepare_tensorqtl_genotype.sh first." >&2
        exit 1
    }
done

# Optional environment activation.
if [[ -n "${CONDA_ENV:-}" ]]; then
    source "${HOME}/.bashrc"
    conda activate "${CONDA_ENV}"
fi

# Create a small chromosome-specific BED from the reusable full BED.
if [[ ! -s "${CHROM_BED}" ]]; then
    gzip -cd "${BED_ALL}" |
    awk -v c="chr${CHROM}" '
        NR==1 || $1==c
    ' |
    gzip -c > "${CHROM_BED}"
fi

N_PHENO="$(gzip -cd "${CHROM_BED}" | awk 'NR>1 {n++} END{print n+0}')"
[[ "${N_PHENO}" -gt 0 ]] || {
    echo "ERROR: no phenotypes found on chr${CHROM}" >&2
    exit 1
}

# Verify BED samples are present in PGEN before launching tensorQTL.
"${PYTHON_BIN}" - "${PGEN_PREFIX}" "${CHROM_BED}" "${MIN_COMMON_SAMPLES}" <<'PY'
import sys
import pandas as pd
import tensorqtl

pgen_prefix, bed_file, min_n = sys.argv[1], sys.argv[2], int(sys.argv[3])

# Read phenotype sample IDs using tensorQTL's BED reader.
phenotype_df, _ = tensorqtl.read_phenotype_bed(bed_file)
pheno_samples = list(map(str, phenotype_df.columns))

# Read PLINK2 sample IDs directly from .psam.
# PLINK2 PSAM files normally contain #FID and IID (or #IID).
psam_file = pgen_prefix + ".psam"
psam = pd.read_csv(psam_file, sep=r"\s+", dtype=str)

psam.columns = [str(c).lstrip("#") for c in psam.columns]

if "IID" in psam.columns:
    geno_samples = set(psam["IID"].astype(str))
elif len(psam.columns) >= 1:
    # Defensive fallback for unusual PSAM headers.
    geno_samples = set(psam.iloc[:, 0].astype(str))
else:
    raise RuntimeError(f"Could not read genotype sample IDs from {psam_file}")

common = [x for x in pheno_samples if x in geno_samples]

print(f"BED samples       : {len(pheno_samples)}")
print(f"PGEN/PSAM samples : {len(geno_samples)}")
print(f"Common samples    : {len(common)}")

if len(common) < min_n:
    raise RuntimeError(
        f"Only {len(common)} common BED/PGEN samples; minimum is {min_n}."
    )

if len(common) != len(pheno_samples):
    missing = [x for x in pheno_samples if x not in geno_samples]
    raise RuntimeError(
        "Some BED genotype_id values are absent from PGEN/PSAM IID: "
        + ", ".join(missing[:20])
    )
PY

EXPECTED="${OUTPUT_DIR}/${PREFIX}.cis_qtl_pairs.chr${CHROM}.parquet"

if [[ "${SKIP_IF_OUTPUT_EXISTS}" == "1" && -s "${EXPECTED}" ]]; then
    echo "Already exists: ${EXPECTED}"
    exit 0
fi

echo "======================================================================"
echo "tensorQTL cis_nominal"
echo "======================================================================"
echo "Config      : ${CONFIG}"
echo "Tissue      : ${TISSUE}"
echo "Cell type   : ${CT}"
echo "Chromosome  : ${CHROM}"
echo "BED         : ${CHROM_BED}"
echo "PGEN        : ${PGEN_PREFIX}"
echo "Output      : ${OUTPUT_DIR}"
echo "Prefix      : ${PREFIX}"
echo "Phenotypes  : ${N_PHENO}"
echo "Cis window  : ${CIS_WINDOW}"
echo "MAF         : ${MAF_THRESHOLD}"
echo "======================================================================"

# No --covariates option is passed because the phenotype is already PEER residuals.
#
# Do NOT add --best_only: we want the complete nominal SNP-gene summary
# statistics for later coloc.
"${PYTHON_BIN}" -m tensorqtl \
    "${PGEN_PREFIX}" \
    "${CHROM_BED}" \
    "${PREFIX}" \
    --mode "${TENSORQTL_MODE}" \
    --window "${CIS_WINDOW}" \
    --maf_threshold "${MAF_THRESHOLD}" \
    --output_dir "${OUTPUT_DIR}"

[[ -s "${EXPECTED}" ]] || {
    echo "ERROR: expected tensorQTL output not found: ${EXPECTED}" >&2
    exit 1
}

echo
echo "DONE: ${EXPECTED}"
