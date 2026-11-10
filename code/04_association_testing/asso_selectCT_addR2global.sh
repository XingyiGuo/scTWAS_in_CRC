#!/bin/bash
#SBATCH --job-name=asso_peakSNP
#SBATCH --error=%x-%j.error
#SBATCH --out=%x-%j.out
#SBATCH --mem=20G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=7-00:00:00
#SBATCH --account=l2_bioinfo1


set -euo pipefail

# ==============================================================================
# SPrediXcan association runner
# Revision: append R2_global and r_squared_global from model DB extra table
#           to each final *.TWAS.csv by gene.
# conda activate spredixcan
#
# Usage:
#   sbatch asso_selectCT_addR2global.sh NL ABS
#   sbatch asso_selectCT_addR2global.sh ADSER ASC
#   sbatch asso_selectCT_addR2global.sh ADSER SSC asso.config.txt
#
# Arguments:
#   1 = tissue/model group: NL or ADSER
#   2 = cell type
#   3 = optional config file; default = config/asso.config.txt
#
# The selected model is tested against all GWAS cohorts listed in COHORTS
# (currently EUR and Multi).
# ==============================================================================

if [[ "$#" -lt 2 || "$#" -gt 3 ]]; then
    echo "Usage: sbatch $(basename "$0") <NL|ADSER> <CELL_TYPE> [asso.config.txt]" >&2
    echo "Examples:" >&2
    echo "  sbatch $(basename "$0") NL ABS" >&2
    echo "  sbatch $(basename "$0") ADSER ASC" >&2
    exit 1
fi

REQUESTED_TISSUE="$1"
REQUESTED_CT="$2"
SCRIPT_DIR=${ASSO_SCRIPT_DIR:-${SLURM_SUBMIT_DIR:-$(dirname "${BASH_SOURCE[0]}")}}
DEFAULT_CONFIG="${SCRIPT_DIR}/config/asso.config.txt"
if [[ -f "${SCRIPT_DIR}/asso.config.txt" ]]; then
    DEFAULT_CONFIG="${SCRIPT_DIR}/asso.config.txt"
fi
CONFIG_FILE="${3:-${DEFAULT_CONFIG}}"

if [[ ! -f "${CONFIG_FILE}" ]]; then
    echo "ERROR: config file not found: ${CONFIG_FILE}. Pass a config or set ASSO_SCRIPT_DIR." >&2
    exit 1
fi

# shellcheck disable=SC1090
source "${CONFIG_FILE}"

mkdir -p "${ASSO_OUTPUT_DIR}"
LOG_DIR="${ASSO_OUTPUT_DIR}/logs"
mkdir -p "${LOG_DIR}"

get_cell_types() {
    local tissue="$1"
    case "${tissue}" in
        NL)
            printf '%s\n' "${NL_CELL_TYPES[@]}"
            ;;
        ADSER)
            printf '%s\n' "${ADSER_CELL_TYPES[@]}"
            ;;
        *)
            echo "ERROR: unsupported tissue/model group: ${tissue}" >&2
            return 1
            ;;
    esac
}

get_gwas_dir() {
    local cohort="$1"
    case "${cohort}" in
        EUR)   printf '%s\n' "${EUR_GWAS_DIR}" ;;
        Multi) printf '%s\n' "${MULTI_GWAS_DIR}" ;;
        *)
            echo "ERROR: unsupported GWAS cohort: ${cohort}" >&2
            return 1
            ;;
    esac
}

get_gwas_pattern() {
    local cohort="$1"
    case "${cohort}" in
        EUR)   printf '%s\n' "${EUR_GWAS_PATTERN}" ;;
        Multi) printf '%s\n' "${MULTI_GWAS_PATTERN}" ;;
        *)
            echo "ERROR: unsupported GWAS cohort: ${cohort}" >&2
            return 1
            ;;
    esac
}


# ------------------------------------------------------------------------------
# Append model-performance metadata from SQLite DB to SPrediXcan output.
#
# SPrediXcan's standard CSV does not automatically expose arbitrary columns
# stored in the model DB "extra" table.  Query the DB directly and merge:
#   R2_global
#   r_squared_global
# by the stable "gene" key.
#
# This uses only Python standard-library modules (sqlite3/csv/os), so it does
# not require pandas.
# ------------------------------------------------------------------------------
append_model_r2() {
    local model_db="$1"
    local twas_csv="$2"

    "${PYTHON_BIN}" - "${model_db}" "${twas_csv}" <<'PY'
import csv
import os
import sqlite3
import sys
import tempfile

db_path, csv_path = sys.argv[1], sys.argv[2]

if not os.path.isfile(db_path):
    raise SystemExit(f"ERROR: model DB not found: {db_path}")
if not os.path.isfile(csv_path):
    raise SystemExit(f"ERROR: TWAS CSV not found: {csv_path}")

required_db_cols = ("gene", "R2_global", "r_squared_global")

con = sqlite3.connect(db_path)
try:
    info = con.execute("PRAGMA table_info(extra)").fetchall()
    db_cols = [row[1] for row in info]

    missing = [c for c in required_db_cols if c not in db_cols]
    if missing:
        raise SystemExit(
            "ERROR: DB extra table is missing required field(s): "
            + ", ".join(missing)
            + f"\nDB: {db_path}"
            + "\nObserved extra fields: " + ", ".join(db_cols)
        )

    rows = con.execute(
        'SELECT gene, R2_global, r_squared_global FROM extra'
    ).fetchall()
finally:
    con.close()

r2_by_gene = {}
duplicate_genes = set()
for gene, r2_global, r_squared_global in rows:
    gene = "" if gene is None else str(gene)
    if gene in r2_by_gene:
        duplicate_genes.add(gene)
    else:
        r2_by_gene[gene] = (r2_global, r_squared_global)

if duplicate_genes:
    examples = ", ".join(sorted(duplicate_genes)[:20])
    raise SystemExit(
        "ERROR: duplicated gene values in DB extra table; refusing ambiguous merge. "
        f"Examples: {examples}"
    )

with open(csv_path, "r", newline="") as fh:
    reader = csv.DictReader(fh)
    if reader.fieldnames is None:
        raise SystemExit(f"ERROR: TWAS CSV has no header: {csv_path}")
    if "gene" not in reader.fieldnames:
        raise SystemExit(
            f"ERROR: TWAS CSV lacks 'gene' column: {csv_path}\n"
            f"Observed columns: {', '.join(reader.fieldnames)}"
        )
    twas_rows = list(reader)
    original_fields = list(reader.fieldnames)

# Keep exact requested field names and avoid duplicate columns on reruns.
new_fields = list(original_fields)
for col in ("R2_global", "r_squared_global"):
    if col not in new_fields:
        new_fields.append(col)

n_matched = 0
n_unmatched = 0

for row in twas_rows:
    gene = str(row.get("gene", ""))
    vals = r2_by_gene.get(gene)

    if vals is None:
        row["R2_global"] = ""
        row["r_squared_global"] = ""
        n_unmatched += 1
    else:
        row["R2_global"] = "" if vals[0] is None else vals[0]
        row["r_squared_global"] = "" if vals[1] is None else vals[1]
        n_matched += 1

directory = os.path.dirname(os.path.abspath(csv_path))
fd, tmp_path = tempfile.mkstemp(
    prefix=os.path.basename(csv_path) + ".",
    suffix=".tmp",
    dir=directory,
    text=True,
)
os.close(fd)

try:
    with open(tmp_path, "w", newline="") as fh:
        writer = csv.DictWriter(
            fh,
            fieldnames=new_fields,
            extrasaction="ignore",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(twas_rows)

    os.replace(tmp_path, csv_path)
except Exception:
    try:
        os.remove(tmp_path)
    except OSError:
        pass
    raise

print(
    f"Appended DB performance fields: matched={n_matched}, "
    f"unmatched={n_unmatched}, output={csv_path}"
)
PY
}

# ------------------------------------------------------------------------------
# Validate requested tissue/cell type
# ------------------------------------------------------------------------------

case "${REQUESTED_TISSUE}" in
    NL|ADSER)
        ;;
    *)
        echo "ERROR: tissue must be NL or ADSER; got: ${REQUESTED_TISSUE}" >&2
        exit 1
        ;;
esac

mapfile -t allowed_cell_types < <(get_cell_types "${REQUESTED_TISSUE}")

valid_ct=false
for x in "${allowed_cell_types[@]}"; do
    if [[ "${x}" == "${REQUESTED_CT}" ]]; then
        valid_ct=true
        break
    fi
done

if [[ "${valid_ct}" != "true" ]]; then
    echo "ERROR: cell type ${REQUESTED_CT} is not configured for ${REQUESTED_TISSUE}." >&2
    echo "Allowed: ${allowed_cell_types[*]}" >&2
    exit 1
fi

# ------------------------------------------------------------------------------
# Preflight
# ------------------------------------------------------------------------------

if [[ ! -f "${SPREDIXCAN}" ]]; then
    echo "ERROR: SPrediXcan.py not found: ${SPREDIXCAN}" >&2
    exit 1
fi

if [[ ! -d "${MODEL_ROOT}" ]]; then
    echo "ERROR: model root not found: ${MODEL_ROOT}" >&2
    exit 1
fi

echo "=============================================================================="
echo "SPrediXcan peakSNP association"
echo "=============================================================================="
echo "Config       : ${CONFIG_FILE}"
echo "SPrediXcan   : ${SPREDIXCAN}"
echo "Model root   : ${MODEL_ROOT}"
echo "Output       : ${ASSO_OUTPUT_DIR}"
echo "Tissue       : ${REQUESTED_TISSUE}"
echo "Cell type    : ${REQUESTED_CT}"
echo "GWAS cohorts : ${COHORTS[*]}"
echo "SNP column   : ${SNP_COLUMN}"
echo "Effect allele: ${EFFECT_ALLELE_COLUMN}"
echo "=============================================================================="

# Validate chr1-22 file discovery for each cohort.
for cohort in "${COHORTS[@]}"; do
    gwas_dir="$(get_gwas_dir "${cohort}")"
    gwas_pattern="$(get_gwas_pattern "${cohort}")"

    if [[ ! -d "${gwas_dir}" ]]; then
        echo "ERROR: GWAS directory not found for ${cohort}: ${gwas_dir}" >&2
        exit 1
    fi

    mapfile -t matched_gwas < <(
        find "${gwas_dir}" -maxdepth 1 -type f -printf '%f\n' \
        | grep -E "${gwas_pattern}" \
        | sort -V || true
    )

    echo "${cohort} GWAS files matched: ${#matched_gwas[@]}"
    if [[ "${#matched_gwas[@]}" -ne 22 ]]; then
        echo "WARNING: expected 22 chromosome files for ${cohort}, found ${#matched_gwas[@]}." >&2
        printf '  %s\n' "${matched_gwas[@]}" >&2
    fi
done

# ------------------------------------------------------------------------------
# Run
# ------------------------------------------------------------------------------

n_run=0
n_skip=0
n_failed=0

tissue="${REQUESTED_TISSUE}"
CT="${REQUESTED_CT}"

model_prefix="${MODEL_DATASET_PREFIX}_${tissue}_${CT}_${MODEL_SUFFIX}"
model_dir="${MODEL_ROOT}/${model_prefix}"

model_db="${model_dir}/${model_prefix}.db"
covariance="${model_dir}/${model_prefix}_cov.txt.gz"

if [[ ! -s "${model_db}" || ! -s "${covariance}" ]]; then
    echo "--------------------------------------------------------------------------"
    echo "MISSING MODEL: ${model_prefix}"
    [[ -s "${model_db}" ]] || echo "  DB : ${model_db}"
    [[ -s "${covariance}" ]] || echo "  COV: ${covariance}"
    exit 1
fi

for cohort in "${COHORTS[@]}"; do
    gwas_dir="$(get_gwas_dir "${cohort}")"
    gwas_pattern="$(get_gwas_pattern "${cohort}")"

    output_file="${ASSO_OUTPUT_DIR}/${model_prefix}.${cohort}.TWAS.csv"
    log_file="${LOG_DIR}/${model_prefix}.${cohort}.SPrediXcan.log"

    if [[ "${SKIP_EXISTING}" == "true" && -s "${output_file}" ]]; then
        # Existing SPrediXcan output may predate the R2_global enhancement.
        # Enrich it in place from the DB without rerunning SPrediXcan.
        header="$(head -n 1 "${output_file}")"
        if [[ "${header}" != *"R2_global"* || "${header}" != *"r_squared_global"* ]]; then
            echo "ENRICH existing TWAS CSV with DB R2 fields: ${output_file}"
            if ! append_model_r2 "${model_db}" "${output_file}"; then
                echo "FAILED to enrich existing result: ${output_file}" >&2
                n_failed=$((n_failed + 1))
                continue
            fi
        fi

        echo "SKIP existing SPrediXcan run: ${output_file}"
        n_skip=$((n_skip + 1))
        continue
    fi

    echo
    echo "=============================================================================="
    echo "Model        : ${model_prefix}"
    echo "Tissue       : ${tissue}"
    echo "Cell type    : ${CT}"
    echo "GWAS cohort  : ${cohort}"
    echo "DB           : ${model_db}"
    echo "Covariance   : ${covariance}"
    echo "GWAS folder  : ${gwas_dir}"
    echo "GWAS pattern : ${gwas_pattern}"
    echo "Output       : ${output_file}"
    echo "=============================================================================="

    # Publish the result only after SPrediXcan and the metadata join succeed.
    run_dir="$(mktemp -d "${ASSO_OUTPUT_DIR}/.${model_prefix}.${cohort}.XXXXXX")"
    tmp_output="${run_dir}/result.csv"
    trap 'rm -f -- "$tmp_output"; rmdir -- "$run_dir"' EXIT

    if "${PYTHON_BIN}" "${SPREDIXCAN}" \
        --model_db_path "${model_db}" \
        --covariance "${covariance}" \
        --gwas_folder "${gwas_dir}" \
        --gwas_file_pattern "${gwas_pattern}" \
        --snp_column "${SNP_COLUMN}" \
        --effect_allele_column "${EFFECT_ALLELE_COLUMN}" \
        --non_effect_allele_column "${NON_EFFECT_ALLELE_COLUMN}" \
        --beta_column "${BETA_COLUMN}" \
        --pvalue_column "${PVALUE_COLUMN}" \
        --output_file "${tmp_output}" \
        --verbosity "${VERBOSITY}" \
        > "${log_file}" 2>&1
    then
        if append_model_r2 "${model_db}" "${tmp_output}" &&
            mv -f -- "${tmp_output}" "${output_file}"; then
            echo "DONE with R2 metadata: ${output_file}"
            n_run=$((n_run + 1))
        else
            status=$?
            echo "FAILED to finalize association (${status}): ${output_file}" >&2
            n_failed=$((n_failed + 1))
        fi
    else
        status=$?
        echo "FAILED (${status}): ${model_prefix} x ${cohort}" >&2
        echo "  Log: ${log_file}" >&2
        n_failed=$((n_failed + 1))
    fi

    rm -f -- "${tmp_output}"
    rmdir -- "${run_dir}"
    trap - EXIT
done

echo
echo "=============================================================================="
echo "Association summary"
echo "=============================================================================="
echo "Completed runs : ${n_run}"
echo "Skipped        : ${n_skip}"
echo "Failed runs    : ${n_failed}"
echo "Output dir     : ${ASSO_OUTPUT_DIR}"
echo "Log dir        : ${LOG_DIR}"
echo "=============================================================================="

if [[ "${n_failed}" -gt 0 ]]; then
    exit 1
fi
