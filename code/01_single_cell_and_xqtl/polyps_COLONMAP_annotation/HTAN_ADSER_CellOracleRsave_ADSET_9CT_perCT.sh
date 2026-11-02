#!/usr/bin/env bash
set -euo pipefail

# ==============================================================================
# HTAN_ADSER_CellOracleRsave_ADSET_9CT_perCT.sh
#
# Environment launcher for:
#   HTAN_ADSER_CellOracleRsave_mergeQC_ADSET_9CT_perCT.R
#
# Usage:
#   bash HTAN_ADSER_CellOracleRsave_ADSET_9CT_perCT.sh
#
# Optional:
#   bash HTAN_ADSER_CellOracleRsave_ADSET_9CT_perCT.sh /data/l2_bioinfo1/liq17/scTWAS/data/VUMC_HTAN/VUMC_HTAN_DIS_VAL_ADSER/src/HTAN_ADSER_CellOracleRsave_mergeQC_PCA83_HPtoSER_ASCSSCfilter.R
# ==============================================================================

# ------------------------------------------------------------------------------
# 1. Load cluster modules
# ------------------------------------------------------------------------------
module load StdEnv/2023
module load r/4.5.0
module load python/3.13.2

# ------------------------------------------------------------------------------
# 2. Activate the Python environment that contains:
#       anndata
#       pandas
#       scipy
#       h5py
# ------------------------------------------------------------------------------
H5AD_VENV="/home/liq17/venvs/r-h5ad"

if [[ ! -x "${H5AD_VENV}/bin/python" ]]; then
    echo "ERROR: Python virtual environment not found:"
    echo "  ${H5AD_VENV}/bin/python"
    exit 1
fi

source "${H5AD_VENV}/bin/activate"

# Explicitly tell the R script which Python interpreter to use.
export HTAN_H5AD_PYTHON="${H5AD_VENV}/bin/python"

# ------------------------------------------------------------------------------
# 3. Check Python dependencies before starting R
# ------------------------------------------------------------------------------
echo "============================================================"
echo "Checking environment"
echo "============================================================"
echo "Rscript:"
which Rscript
Rscript --version

echo
echo "Python:"
echo "${HTAN_H5AD_PYTHON}"
"${HTAN_H5AD_PYTHON}" --version

echo
echo "Checking h5ad Python packages..."
"${HTAN_H5AD_PYTHON}" -c \
'import anndata, pandas, scipy, h5py; print("Python h5ad environment: OK")'

# ------------------------------------------------------------------------------
# 4. Run the R script
# ------------------------------------------------------------------------------
R_SCRIPT="/data/l2_bioinfo1/liq17/scTWAS/data/VUMC_HTAN/VUMC_HTAN_DIS_VAL_ADSER/src/HTAN_ADSER_CellOracleRsave_mergeQC_PCA83_HPtoSER_ASCSSCfilter.R"

if [[ ! -f "${R_SCRIPT}" ]]; then
    echo
    echo "ERROR: R script not found:"
    echo "  ${R_SCRIPT}"
    echo
    echo "Either run this launcher from the directory containing the R script,"
    echo "or provide the full R-script path:"
    echo "  bash run_01a_HTAN_ADSER.sh /full/path/to/HTAN_ADSER_CellOracleRsave_mergeQC_PCA83_HPtoSER_ASCSSCfilter.R"
    exit 1
fi

echo
echo "============================================================"
echo "Running HTAN AD/SER preparation"
echo "============================================================"
echo "R script: ${R_SCRIPT}"
echo "Python : ${HTAN_H5AD_PYTHON}"
echo

Rscript "${R_SCRIPT}"
