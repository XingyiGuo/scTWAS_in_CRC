#!/bin/bash
#SBATCH --account=l2_bioinfo1
#SBATCH --job-name=umap_ws_allcells
#SBATCH --error=%x-%j.error
#SBATCH --output=%x-%j.out
#SBATCH --mem=160G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=1-0:0:0

# =============================================================================
# run_umap_WS_allcells.sh - regenerate PANEL B of Final_Viz_4panel.R from EVERY
# retained epithelial cell instead of the 4,000-per-type sample.
#
#   sbatch run_umap_WS_allcells.sh
#
# WHY THIS AND NOT AN EXTRACTION SCRIPT
#   Panel B's coordinates are not stored in any .rds, so there is nothing to
#   extract: 04_umap_epithelial_LauLabel_revised.R COMPUTES the embedding on
#   each run, straight from the BPCells matrix
#     <OUT_DIR>/Epithelial_expression_BPCELLS/
#   (LogNormalize 1e4 -> 2,000 vst HVG -> ScaleData -> 50 PCs -> k=15 ->
#    Louvain res 0.5 -> UMAP min.dist 0.5), and it already writes exactly the
#   table the figure reads:
#     <OUT_DIR>/umap/umap_WellcomeSanger_epithelial_clusters_LauLabel_sourcedata.csv
#   The only thing the current table lacks is scope - it holds 28,000 cells
#   (4,000 per LauLabel, seed 20260923) out of the 358,207 retained.
#
#   WS_CELLS_PER_TYPE=all removes the cap.  The script's own header asks for a
#   node with >= 64 GB for that; 160G here matches the rest of this project's
#   jobs, and the 8 CPUs are for the BLAS/ANN threads in PCA, kNN and UMAP.
#
# AFTERWARDS
#   Copy umap_WellcomeSanger_epithelial_clusters_LauLabel_sourcedata.csv to
#     L:/scRNA-seq/WellcomeSanger/Healthy_TI_then_rectum_278samples_outputs_LauLabel/umap/
#   Final_Viz_4panel.R reads that path first, so panel B follows with no edit:
#   it will plot every row and drop its point size to 0.035 automatically.
#   Re-check B_NUDGE in the figure script afterwards - the STM and TAC label
#   offsets were tuned on the capped layout, and an uncapped embedding is a new
#   layout, not a denser version of the same one.
# =============================================================================

set -euo pipefail

# Where 04_umap_epithelial_LauLabel_revised.R lives on ACCRE.
SCRIPT_DIR="${SCRIPT_DIR:-/data/l2_bioinfo1/liq17/scTWAS/src/WellcomeSanger}"
R_SCRIPT="${SCRIPT_DIR}/04_umap_epithelial_LauLabel_revised.R"

# Same module set as run_pipeline.sh in that folder; hdf5 is needed by BPCells.
module load StdEnv/2023
module load gcc/13.3
module load hdf5/1.14.5
module load r/4.5.0

export WS_CELLS_PER_TYPE=all        # <- the whole point of this job
export WS_CLUSTER_RES=0.5           # unchanged, so clusters stay comparable
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-8}
export OPENBLAS_NUM_THREADS=${OMP_NUM_THREADS}

# The script resolves OUT_DIR itself; set WS_LAULABEL_OUT_DIR to override.
# export WS_LAULABEL_OUT_DIR=/data/l2_bioinfo1/liq17/scTWAS/data/WellcomeSanger/expression/Healthy_TI_then_rectum_278samples_outputs_LauLabel

if [[ ! -f "${R_SCRIPT}" ]]; then
    echo "ERROR: R script not found: ${R_SCRIPT}" >&2
    exit 1
fi

echo "============================================================"
echo "WellcomeSanger epithelial UMAP - UNCAPPED (all retained cells)"
echo "start : $(date)"
echo "host  : $(hostname)"
echo "script: ${R_SCRIPT}"
echo "cap   : ${WS_CELLS_PER_TYPE}"
echo "R     : $(which R)"
R --version | head -1
echo "============================================================"

# The script writes the capped and uncapped runs to the SAME filenames, and the
# capped table is what the shipped figure was built from, so keep a copy before
# it is overwritten.  Nothing here deletes anything.
OUT_DIR="${WS_LAULABEL_OUT_DIR:-}"
if [[ -z "${OUT_DIR}" ]]; then
    for c in /data/l2_bioinfo1/liq17/scTWAS/data/WellcomeSanger/expression/Healthy_TI_then_rectum_278samples_outputs_LauLabel \
             /v5000/data/l2_bioinfo1/liq17/scTWAS/data/WellcomeSanger/expression/Healthy_TI_then_rectum_278samples_outputs_LauLabel \
             /v5000/data/l2_bioinfo1/scRNA-seq/WellcomeSanger/processed/E-MTAB-16986/Healthy_TI_then_rectum_278samples_outputs_LauLabel; do
        [[ -d "$c" ]] && OUT_DIR="$c" && break
    done
fi
if [[ -n "${OUT_DIR}" && -d "${OUT_DIR}/umap" ]]; then
    BACKUP="${OUT_DIR}/umap/capped_4000_per_type"
    mkdir -p "${BACKUP}"
    for f in "${OUT_DIR}"/umap/umap_WellcomeSanger_epithelial_clusters_LauLabel*; do
        [[ -f "$f" ]] && cp -n "$f" "${BACKUP}/" || true
    done
    echo "capped run backed up to: ${BACKUP}"
fi

Rscript "${R_SCRIPT}"

echo "============================================================"
echo "done  : $(date)"
echo "outputs in: ${OUT_DIR:-<resolved by the R script>}/umap"
echo "============================================================"
