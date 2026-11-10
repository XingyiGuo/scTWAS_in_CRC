#!/bin/bash
#SBATCH --job-name=gse201336_peaks
#SBATCH --cpus-per-task=2
#SBATCH --mem=100G
#SBATCH --time=7-00:00:00
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
#SBATCH	--account=l2_bioinfo1

set -euo pipefail

module load StdEnv/2023
module load r/4.5.0 python/3.9 || true

# Submit from the script directory or set GSE_SCRIPT_DIR to the source directory.
# WORK holds the staged inputs and generated outputs.
SCRIPT_DIR=${GSE_SCRIPT_DIR:-${SLURM_SUBMIT_DIR:-$(dirname "${BASH_SOURCE[0]}")}}
SCRIPT_DIR=$(cd "$SCRIPT_DIR" && pwd)
WORK=$(cd "${WORK:-$PWD}" && pwd)
FRAG=$WORK/fragments
ARCHR=$WORK/archr
LABELS=$WORK/epithelial_celltypes_atac.tsv
LABEL_MAP=$WORK/lau_label_map.csv
PEAKS=$WORK/peaks
THREADS=${SLURM_CPUS_PER_TASK:-16}
MANIFEST=$WORK/GSE201336_atac_manifest.csv

export TMPDIR=${SLURM_TMPDIR:-/tmp}
cd "$WORK"

Rscript "$SCRIPT_DIR/01_build_arrows.R" \
  --fragments "$FRAG" --manifest "$MANIFEST" \
  --out "$ARCHR" --threads "$THREADS" \
  --conditions "Normal,Unaffected,Adenoma"

Rscript "$SCRIPT_DIR/02_relabel.R" \
  --archr "$ARCHR" --labels "$LABELS" --manifest "$MANIFEST" \
  --labelMap "$LABEL_MAP" --threads "$THREADS"

# --- 3) reproducible peak set per Condition x CellType, called SEPARATELY --------
# Normal  -> STM, TAC, ABS, CT, GOB, EE
# Adenoma -> STM, TAC, ABS, CT, GOB, EE
# (--byCondition = each condition's cell types get their own peak set; no pooling.)
macs_args=()
if [[ -n "${MACS2:-}" ]]; then
  macs_args=(--macs2 "$MACS2")
fi
Rscript "$SCRIPT_DIR/03_call_peaks.R" --archr "$ARCHR" --threads "$THREADS" --minCells 200 \
  --byCondition --conditions "Normal,Adenoma" "${macs_args[@]}"

# --- 4) export (no comparison; pass --diff to add differential tables) ---------
Rscript "$SCRIPT_DIR/04_export_peaks.R" --archr "$ARCHR" --out "$PEAKS" --threads "$THREADS"

echo "Done. Peaks in $PEAKS"
