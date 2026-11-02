#!/usr/bin/env bash
#SBATCH --job-name=step3
#SBATCH --error=step3_%A_%a.error
#SBATCH --output=step3_%A_%a.out
#SBATCH --mem=40G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --time=7-00:00:00
#SBATCH --account=l2_bioinfo1
#SBATCH --array=3

set -euo pipefail

# ============================================================
# Step 3: Run PEER
# One cell type per SLURM array job
# ============================================================

SCRIPT_DIR="/data/l2_bioinfo1/liq17/scTWAS/data/VUMC_HTAN/src"

EXPRESSION_ROOT="/v5000/data/l2_bioinfo1/liq17/scTWAS/data/VUMC_HTAN/expression/ADSER_pseudobulk_expression_for_peer"

CELL_TYPES=(ABS EE GOB CT TUF TAC STM)
#CELL_TYPES=(CT)

CT="${CELL_TYPES[$SLURM_ARRAY_TASK_ID]}"

cd "${SCRIPT_DIR}"

echo "============================================================"
echo "Step 3 PEER"
echo "Job ID       : ${SLURM_JOB_ID}"
echo "Array task   : ${SLURM_ARRAY_TASK_ID}"
echo "Cell type    : ${CT}"
echo "Host         : $(hostname)"
echo "Start time   : $(date)"
echo "============================================================"

python "${SCRIPT_DIR}/03_Process_peer_HTAN_ASC_SSC_HTANdirect_geneSymbol.py" \
    "${CT}" \
    --expression-root "${EXPRESSION_ROOT}"

echo "============================================================"
echo "[DONE] PEER ${CT}"
echo "End time     : $(date)"
echo "============================================================"
