#!/usr/bin/env bash
set -euo pipefail

# ==============================================================================
# Create a fresh CPTAC Python environment under scTWAS/src and launch analysis
#
# Environment:
#   /v5000/data/l2_bioinfo1/liq17/scTWAS/src/.venv-cptac
#
# Analysis script:
#   /v5000/data/l2_bioinfo1/liq17/scTWAS/src/Functional_src/
#       CPTAC_COAD_TWAS_annotation.py
#
# This script:
#   1. Finds Python 3.11
#   2. Removes any existing .venv-cptac
#   3. Creates a fresh virtual environment
#   4. Installs/pins required packages
#   5. Tests imports
#   6. Saves exact package versions
#   7. Launches CPTAC_COAD_TWAS_annotation.py
#
# Usage:
#   bash create_and_launch_cptac_env.sh
# ==============================================================================


# ==============================================================================
# 0. Paths
# ==============================================================================

PROJECT_DIR="/v5000/data/l2_bioinfo1/liq17/scTWAS/src"

VENV_DIR="${PROJECT_DIR}/.venv-cptac"

ANALYSIS_SCRIPT="${PROJECT_DIR}/Functional_src/CPTAC_COAD_TWAS_annotation.py"

REQUIREMENTS_FILE="${PROJECT_DIR}/requirements-cptac.txt"


# ==============================================================================
# 1. Enter project directory
# ==============================================================================

echo
echo "============================================================"
echo "CPTAC Python environment setup"
echo "============================================================"
echo "Project directory:"
echo "  ${PROJECT_DIR}"
echo

if [[ ! -d "${PROJECT_DIR}" ]]; then
    echo "ERROR: Project directory does not exist:"
    echo "  ${PROJECT_DIR}"
    exit 1
fi

cd "${PROJECT_DIR}"


# ==============================================================================
# 2. Find Python 3.11
# ==============================================================================

echo "============================================================"
echo "Finding Python 3.11"
echo "============================================================"

PYTHON_BIN=""

# Prefer an explicit python3.11 executable.
if command -v python3.11 >/dev/null 2>&1; then
    PYTHON_BIN="$(command -v python3.11)"
else
    # Fall back to python3 only if it is actually Python 3.11.
    if command -v python3 >/dev/null 2>&1; then
        PY3_VERSION="$(python3 -c 'import sys; print(f"{sys.version_info.major}.{sys.version_info.minor}")')"

        if [[ "${PY3_VERSION}" == "3.11" ]]; then
            PYTHON_BIN="$(command -v python3)"
        fi
    fi
fi

if [[ -z "${PYTHON_BIN}" ]]; then
    echo
    echo "ERROR: Python 3.11 was not found."
    echo
    echo "Current Python executables:"
    command -v python || true
    command -v python3 || true
    python --version 2>/dev/null || true
    python3 --version 2>/dev/null || true
    echo
    echo "Load an available Python 3.11 module first, then rerun this script."
    echo "For example, inspect available modules with:"
    echo
    echo "  module spider python"
    echo
    exit 1
fi

echo "Using:"
echo "  ${PYTHON_BIN}"
"${PYTHON_BIN}" --version


# ==============================================================================
# 3. Remove any existing CPTAC virtual environment
# ==============================================================================

echo
echo "============================================================"
echo "Removing old environment, if present"
echo "============================================================"

if [[ -d "${VENV_DIR}" ]]; then
    echo "Removing:"
    echo "  ${VENV_DIR}"
    rm -rf "${VENV_DIR}"
else
    echo "No existing environment found."
fi


# ==============================================================================
# 4. Create fresh virtual environment
# ==============================================================================

echo
echo "============================================================"
echo "Creating virtual environment"
echo "============================================================"

"${PYTHON_BIN}" -m venv "${VENV_DIR}"

if [[ ! -x "${VENV_DIR}/bin/python" ]]; then
    echo "ERROR: virtual environment creation failed."
    exit 1
fi


# ==============================================================================
# 5. Activate environment
# ==============================================================================

# shellcheck disable=SC1091
source "${VENV_DIR}/bin/activate"

echo
echo "Activated environment:"
echo "  ${VENV_DIR}"

echo
echo "Python executable:"
which python

echo
echo "Python version:"
python --version

echo
echo "pip version:"
python -m pip --version


# ==============================================================================
# 6. Install packaging tools
# ==============================================================================

echo
echo "============================================================"
echo "Installing packaging tools"
echo "============================================================"

# Upgrade pip and wheel.
python -m pip install --upgrade pip wheel

# IMPORTANT:
# cptac -> pyranges -> sorted_nearest imports pkg_resources.
# pkg_resources is unavailable in newer setuptools releases.
# Pin setuptools to a version that still provides pkg_resources.
python -m pip install --force-reinstall "setuptools==80.8.0"


# ==============================================================================
# 7. Install CPTAC dependencies
# ==============================================================================

echo
echo "============================================================"
echo "Installing CPTAC analysis packages"
echo "============================================================"

python -m pip install \
    numpy \
    pandas \
    scipy \
    statsmodels \
    "cptac==1.5.14"

# Re-pin setuptools after dependency installation.
python -m pip install --force-reinstall "setuptools==80.8.0"


# ==============================================================================
# 8. Verify pkg_resources and all imports
# ==============================================================================

echo
echo "============================================================"
echo "Testing installed packages"
echo "============================================================"

python - <<'PY'
import sys

import pkg_resources
import numpy as np
import pandas as pd
import scipy
import statsmodels
import cptac

print("Python       :", sys.version.split()[0])
print("numpy        :", np.__version__)
print("pandas       :", pd.__version__)
print("scipy        :", scipy.__version__)
print("statsmodels  :", statsmodels.__version__)
print("cptac        :", getattr(cptac, "__version__", "version not reported"))
print("pkg_resources: OK")
print()
print("All required imports completed successfully.")
PY


# ==============================================================================
# 9. Save exact package versions
# ==============================================================================

echo
echo "============================================================"
echo "Saving environment versions"
echo "============================================================"

python -m pip freeze > "${REQUIREMENTS_FILE}"

echo "Saved:"
echo "  ${REQUIREMENTS_FILE}"

echo
echo "Core installed packages:"
python -m pip freeze | grep -Ei \
'^(cptac|numpy|pandas|scipy|statsmodels|setuptools|pyranges|sorted-nearest|sorted_nearest)==' \
|| true


# ==============================================================================
# 10. Check analysis script
# ==============================================================================

echo
echo "============================================================"
echo "Checking CPTAC analysis script"
echo "============================================================"

if [[ ! -f "${ANALYSIS_SCRIPT}" ]]; then
    echo "ERROR: analysis script was not found:"
    echo "  ${ANALYSIS_SCRIPT}"
    exit 1
fi

echo "Found:"
echo "  ${ANALYSIS_SCRIPT}"


# ==============================================================================
# 11. Launch CPTAC analysis
# ==============================================================================

echo
echo "============================================================"
echo "Launching CPTAC COAD TWAS annotation"
echo "============================================================"
echo

python "${ANALYSIS_SCRIPT}"


# ==============================================================================
# 12. Finished
# ==============================================================================

echo
echo "============================================================"
echo "CPTAC analysis finished successfully"
echo "============================================================"
echo
echo "Environment:"
echo "  ${VENV_DIR}"
echo
echo "Requirements:"
echo "  ${REQUIREMENTS_FILE}"
echo
echo "To run the analysis again later without reinstalling:"
echo
echo "  cd ${PROJECT_DIR}"
echo "  source .venv-cptac/bin/activate"
echo "  python Functional_src/CPTAC_COAD_TWAS_annotation.py"
echo
