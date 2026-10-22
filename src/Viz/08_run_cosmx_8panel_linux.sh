#!/usr/bin/env bash
set -euo pipefail

# ---- paths ------------------------------------------------------------------
export COSMX_ROOT="${COSMX_ROOT:-/data/l2_bioinfo1/liq17/scTWAS/data/SpatialExpr/CosMx_PMID41986711}"
export COSMX_ARCHIVE_ROOT="${COSMX_ARCHIVE_ROOT:-$COSMX_ROOT}"
export COSMX_FIG_WORKDIR="${COSMX_FIG_WORKDIR:-$PWD/cosmx_8panel_work}"
export COSMX_FIG_OUTDIR="${COSMX_FIG_OUTDIR:-$COSMX_FIG_WORKDIR/figure}"

# Optional HPC module examples -- uncomment/adapt to ACCRE environment:
# module purge
# module load Python/3.11.5
# module load R/4.3.2

# Optional isolated Python environment:
# python3 -m venv .venv-cosmx
# source .venv-cosmx/bin/activate
# python -m pip install --upgrade pip
# python -m pip install numpy pandas

# Optional user R library, if your cluster needs one:
# export R_LIBS_USER="$HOME/R/library"
# mkdir -p "$R_LIBS_USER"

python3 08_prep_cosmx_highlight_8panel_linux.py
Rscript 09_fig_NOR_CORE_highlight_8panel_linux.R
