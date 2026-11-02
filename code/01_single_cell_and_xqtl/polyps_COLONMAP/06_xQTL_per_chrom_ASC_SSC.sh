#!/usr/bin/env bash
#SBATCH --job-name=xQTL_ASC_SSC
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
# xQTL_per_chrom_ASC_SSC_bashconfig_fixed.sh
#
# tensorQTL cis_nominal for VUMC HTAN GSE266546 ASC / SSC.
#
# IMPORTANT:
#   This version expects a BASH-style *.eqtl.config.txt and loads it with:
#       source "${CONFIG}"
#
# Recommended config:
#   HTAN_GSE266546_ASC_SSC.eqtl.config_COLONMAP95.txt
#
# Usage, one chromosome:
#   sbatch 06_xQTL_per_chrom_ASC_SSC.sh ASC 1 \
#       HTAN_GSE266546_ASC_SSC.eqtl.config_COLONMAP95.txt
#
# Usage, SLURM array:
#   sbatch --array=1-22 xQTL_per_chrom_ASC_SSC_bashconfig_fixed.sh ASC \
#       HTAN_GSE266546_ASC_SSC.eqtl.config_COLONMAP95.txt
#
# In array mode, argument 2 is interpreted as CONFIG when it is non-numeric.
# ==============================================================================

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

[[ $# -ge 1 ]] || {
    echo "Usage:" >&2
    echo "  $0 <ASC|SSC> <chrom> [config]" >&2
    echo "  sbatch --array=1-22 $0 <ASC|SSC> [config]" >&2
    exit 2
}

CT_ARG="$1"

DEFAULT_CONFIG="${SCRIPT_DIR}/HTAN_ADSER_9cts.eqtl.config_COLONMAP95.txt"

# ------------------------------------------------------------------------------
# Parse chromosome/config robustly for both direct and array use.
# ------------------------------------------------------------------------------

if [[ $# -ge 2 && "$2" =~ ^([1-9]|1[0-9]|2[0-2])$ ]]; then
    # Direct mode:
    #   script ASC 1 config.txt
    CHROM="$2"
    CONFIG="${3:-${DEFAULT_CONFIG}}"
else
    # Array mode:
    #   sbatch --array=1-22 script ASC config.txt
    CHROM="${SLURM_ARRAY_TASK_ID:-}"
    CONFIG="${2:-${DEFAULT_CONFIG}}"
fi

[[ -n "${CHROM}" ]] || {
    echo "ERROR: chromosome is required, either as argument 2 or SLURM_ARRAY_TASK_ID." >&2
    exit 2
}

if ! [[ "${CHROM}" =~ ^([1-9]|1[0-9]|2[0-2])$ ]]; then
    echo "ERROR: chromosome must be 1-22; got '${CHROM}'." >&2
    exit 1
fi

[[ -s "${CONFIG}" ]] || {
    echo "ERROR: config not found: ${CONFIG}" >&2
    exit 1
}

# ------------------------------------------------------------------------------
# Load BASH configuration.
# ------------------------------------------------------------------------------

# shellcheck disable=SC1090
source "${CONFIG}"

# Required config variables.
required_vars=(
    WORK_DIR
    DATASET_PREFIX
    TISSUE
    PHENOTYPE_BED_TEMPLATE
    TENSORQTL_PGEN_PREFIX_TEMPLATE
    OUTPUT_DIR_TEMPLATE
    PREFIX_TEMPLATE
    PYTHON_BIN
    CIS_WINDOW
    MAF_THRESHOLD
    MIN_COMMON_SAMPLES
    TENSORQTL_MODE
    CHROM_BED_SUBDIR
    AUDIT_SUBDIR
    SKIP_IF_OUTPUT_EXISTS
)

for v in "${required_vars[@]}"; do
    [[ -n "${!v:-}" ]] || {
        echo "ERROR: required config variable '${v}' is missing/empty in ${CONFIG}" >&2
        exit 1
    }
done

# ------------------------------------------------------------------------------
# Resolve ASC / SSC.
# ------------------------------------------------------------------------------

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
        echo "ERROR: unknown cell type '${CT}'." >&2
        echo "Allowed cell types: ${CELL_TYPES[*]}" >&2
        exit 1
    }
fi

# ------------------------------------------------------------------------------
# Expand config templates.
# ------------------------------------------------------------------------------

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

# ------------------------------------------------------------------------------
# Required-input checks.
# ------------------------------------------------------------------------------

for f in \
    "${BED_ALL}" \
    "${PGEN_PREFIX}.pgen" \
    "${PGEN_PREFIX}.pvar" \
    "${PGEN_PREFIX}.psam"
do
    [[ -s "${f}" ]] || {
        echo "ERROR: required input missing: ${f}" >&2
        echo >&2
        echo "If BED is missing, run:" >&2
        echo "  bash prepare_tensorqtl_beds_ASC_SSC.sh ${CONFIG}" >&2
        echo >&2
        echo "If tensorQTL PGEN is missing, run:" >&2
        echo "  bash prepare_tensorqtl_genotype.sh ${CONFIG}" >&2
        exit 1
    }
done

# ------------------------------------------------------------------------------
# Optional conda activation.
# ------------------------------------------------------------------------------

if [[ -n "${CONDA_ENV:-}" ]]; then
    source "${HOME}/.bashrc"
    conda activate "${CONDA_ENV}"
fi

# ------------------------------------------------------------------------------
# Create chromosome-specific phenotype BED.
# ------------------------------------------------------------------------------

if [[ ! -s "${CHROM_BED}" ]]; then
    gzip -cd "${BED_ALL}" |
    awk -v c="chr${CHROM}" '
        NR == 1 || $1 == c
    ' |
    gzip -c > "${CHROM_BED}"
fi

N_PHENO="$(
    gzip -cd "${CHROM_BED}" |
    awk 'NR > 1 {n++} END {print n+0}'
)"

[[ "${N_PHENO}" -gt 0 ]] || {
    echo "ERROR: no phenotypes found on chr${CHROM}" >&2
    exit 1
}

# ------------------------------------------------------------------------------
# Verify phenotype IDs against PGEN/PSAM.
# ------------------------------------------------------------------------------

"${PYTHON_BIN}" - \
    "${PGEN_PREFIX}" \
    "${CHROM_BED}" \
    "${MIN_COMMON_SAMPLES}" \
    "${AUDIT_DIR}" \
    "${PREFIX}" \
    "${CHROM}" <<'PY'
import os
import sys
import pandas as pd
import tensorqtl

pgen_prefix = sys.argv[1]
bed_file = sys.argv[2]
min_n = int(sys.argv[3])
audit_dir = sys.argv[4]
prefix = sys.argv[5]
chrom = sys.argv[6]

phenotype_df, _ = tensorqtl.read_phenotype_bed(bed_file)
pheno_samples = list(map(str, phenotype_df.columns))

psam_file = pgen_prefix + ".psam"
psam = pd.read_csv(psam_file, sep=r"\s+", dtype=str)
psam.columns = [str(c).lstrip("#") for c in psam.columns]

if "IID" in psam.columns:
    geno_samples = set(psam["IID"].astype(str))
elif len(psam.columns) >= 1:
    geno_samples = set(psam.iloc[:, 0].astype(str))
else:
    raise RuntimeError(f"Could not read genotype sample IDs from {psam_file}")

common = [x for x in pheno_samples if x in geno_samples]
missing = [x for x in pheno_samples if x not in geno_samples]

print(f"BED samples       : {len(pheno_samples)}")
print(f"PGEN/PSAM samples : {len(geno_samples)}")
print(f"Common samples    : {len(common)}")

audit_file = os.path.join(
    audit_dir,
    f"{prefix}.chr{chrom}.sample_match.tsv"
)

pd.DataFrame({
    "phenotype_sample": pheno_samples,
    "in_PGEN": [x in geno_samples for x in pheno_samples],
}).to_csv(audit_file, sep="\t", index=False)

print(f"Sample audit      : {audit_file}")

if len(common) < min_n:
    raise RuntimeError(
        f"Only {len(common)} common BED/PGEN samples; minimum is {min_n}."
    )

if missing:
    raise RuntimeError(
        "Some BED genotype_id values are absent from PGEN/PSAM IID: "
        + ", ".join(missing[:20])
    )
PY

# ------------------------------------------------------------------------------
# tensorQTL cis_nominal
# ------------------------------------------------------------------------------

EXPECTED="${OUTPUT_DIR}/${PREFIX}.cis_qtl_pairs.chr${CHROM}.parquet"

if [[ "${SKIP_IF_OUTPUT_EXISTS}" == "1" && -s "${EXPECTED}" ]]; then
    echo "Already exists: ${EXPECTED}"
    exit 0
fi

echo "======================================================================"
echo "tensorQTL cis_nominal"
echo "======================================================================"
echo "Config      : ${CONFIG}"
echo "Dataset     : ${DATASET_PREFIX}"
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
echo "Min common N: ${MIN_COMMON_SAMPLES}"
echo "======================================================================"

# PEER residual expression is already covariate-adjusted.
# Keep full nominal SNP-gene statistics for coloc / SMR follow-up.
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
