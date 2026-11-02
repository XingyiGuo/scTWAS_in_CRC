#!/bin/bash
# =============================================================================
# Submit the 05e-style ADSER polyp UMAP.
#   sbatch 09a_umap05e_ADSER_polyp_celltype_sbatch.sh
#
# Sizing: 09 embeds all 81,929 model cells by default (CELLS_PER_TYPE = NA),
# and the floor is reading the 904 MB / 196,672-cell merged object before the
# subset, so 64 GB. If it has already run once, the script reuses
# Results/UMAP_ADSER_model_cells/objH_ADSER_model_cells_umap.rds and the job
# needs far less -- drop to --mem=24G for a figure-only rerun.
# Threads are pinned to 1: the first full-panel run of 08 died in the OpenMP
# layer asking for 346,711,037,184 bytes (~323 GiB), and that request grows
# with the thread count. Raise --mem, not the threads, if it recurs.
# =============================================================================
#SBATCH --job-name=umap05e_ADSER
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=06:00:00
#SBATCH --output=logs_umap/umap05e_ADSER-%j.out
#SBATCH --error=logs_umap/umap05e_ADSER-%j.error

set -euo pipefail

PROJ="/data/l2_bioinfo1/liq17/scTWAS/data/VUMC_HTAN/VUMC_HTAN_DIS_VAL_ADSER"
SRC="${PROJ}/src"
mkdir -p "${SRC}/logs_umap"
cd "${SRC}"

# ---- R environment ----------------------------------------------------------
# FILL THIS IN once for this cluster, then leave it alone. Use whichever of the
# three applies and delete the others; the script needs Seurat v5, harmony,
# ggplot2, dplyr, readr, and a cairo-capable R build (cairo_pdf + bitmapType).
#   module load R/4.3.1
#   source /data/l2_bioinfo1/liq17/miniconda3/etc/profile.d/conda.sh && conda activate r_seurat
#   export R_LIBS_USER=/data/l2_bioinfo1/liq17/R/library
module list 2>&1 || true

# Pin the compiled layers. THIS IS THE LOAD-BEARING PART, not a belt-and-braces
# extra: libgomp reads OMP_NUM_THREADS in its own library constructor, so a
# Sys.setenv inside R runs too late. Launch the script any other way -- Rscript
# straight from a login shell, or an srun without these exports -- and libgomp
# asks for one thread per core on the node, the cgroup refuses, and the job
# dies at RunHarmony with
#   libgomp: Thread creation failed: Resource temporarily unavailable
# The script preflights this and stops in seconds rather than after the readRDS.
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OMP_THREAD_LIMIT=1
export OMP_DYNAMIC=FALSE
export OPENBLAS_MAIN_FREE=1
export RCPP_PARALLEL_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
# passive spin-wait: idle OpenMP threads sleep instead of busy-waiting, which
# keeps a stray team from holding slots against the job's thread allowance.
export OMP_WAIT_POLICY=passive

# Raise the per-user process/thread cap if the queue lets us. Harmless when it
# is already at the hard limit; `|| true` keeps `set -e` from killing the job.
ulimit -u unlimited 2>/dev/null || true

# Record the limits the failure actually depends on, so the next .out file
# localises it without a second run.
echo "--- thread/process limits ---"
echo "cpus-per-task = ${SLURM_CPUS_PER_TASK:-NA}  nproc = $(nproc 2>/dev/null || echo NA)"
echo "ulimit -u (max user processes) = $(ulimit -u)"
for f in /sys/fs/cgroup/pids/pids.max /sys/fs/cgroup/pids.max; do
  [ -r "$f" ] && echo "$f = $(cat "$f")"
done
echo "OMP_NUM_THREADS = ${OMP_NUM_THREADS}"
echo "-----------------------------"

# --vanilla skips .Rprofile, which also skips renv's autoloader. If the script
# must use the project library, drop --vanilla (or set R_LIBS_USER to the renv
# library) and run renv::status() once to see what drifted.
Rscript --vanilla "${SRC}/09_umap05e_ADSER_polyp_celltype.R"

echo "exit=$? job=${SLURM_JOB_ID:-NA} $(date)"
