#!/bin/bash
# =============================================================================
# Submit the ADSER model-cell UMAP.
#   sbatch 08a_umap_ADSER_model_cells_sbatch.sh
# Sizing: the script draws <=4,000 cells per cell type (~26.5k of the 81,929
# model cells), so the working set is small; the floor is set by reading the
# 904 MB / 196,672-cell merged object before the subset. 48 GB is comfortable.
# Threads are pinned to 1 below: the first full-panel run died in the OpenMP
# layer asking for 346,711,037,184 bytes (~323 GiB), and that request grows with
# the thread count. Raise --mem, not the threads, if it recurs.
# =============================================================================
#SBATCH --job-name=umap_ADSER
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=48G
#SBATCH --time=04:00:00
#SBATCH --output=logs_umap/umap_ADSER-%j.out
#SBATCH --error=logs_umap/umap_ADSER-%j.error

set -euo pipefail

PROJ="/data/l2_bioinfo1/liq17/scTWAS/data/VUMC_HTAN/VUMC_HTAN_DIS_VAL_ADSER"
SRC="${PROJ}/src"
mkdir -p "${SRC}/logs_umap"
cd "${SRC}"

# ---- R environment ----------------------------------------------------------
# FILL THIS IN once for this cluster, then leave it alone. Use whichever of the
# three applies and delete the others; the script needs Seurat v5, harmony,
# ggplot2, dplyr, readr (shadowtext optional, only sharpens the on-plot labels).
#   module load R/4.3.1
#   source /data/l2_bioinfo1/liq17/miniconda3/etc/profile.d/conda.sh && conda activate r_seurat
#   export R_LIBS_USER=/data/l2_bioinfo1/liq17/R/library
module list 2>&1 || true

# Pin the compiled layers. Some BLAS/OpenMP builds read these only at load
# time, so exporting here is what actually takes effect -- the same assignment
# inside R is a belt-and-braces fallback.
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OMP_THREAD_LIMIT=1

# --vanilla skips .Rprofile, which also skips renv's autoloader. The last run
# warned "project is out-of-sync"; if the script must use the project library,
# drop --vanilla (or set R_LIBS_USER to the renv library) and run renv::status()
# once to see what drifted.
Rscript --vanilla "${SRC}/08_umap_ADSER_model_cells.R"

echo "exit=$? job=${SLURM_JOB_ID:-NA} $(date)"
