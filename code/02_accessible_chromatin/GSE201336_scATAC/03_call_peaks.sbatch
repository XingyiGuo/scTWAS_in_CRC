#!/bin/bash
#SBATCH --job-name=gse201336_step03
#SBATCH --cpus-per-task=2
#SBATCH --mem=100G
#SBATCH --time=7-00:00:00
#SBATCH --output=/data/l2_bioinfo1/scATAC-seq/scATACseq_GSE201336/GSE201336_scATAC_peaks/logs/%x_%j.out
#SBATCH --error=/data/l2_bioinfo1/scATAC-seq/scATACseq_GSE201336/GSE201336_scATAC_peaks/logs/%x_%j.err
#SBATCH --account=l2_bioinfo1

set -euo pipefail

WORK=/data/l2_bioinfo1/scATAC-seq/scATACseq_GSE201336/GSE201336_scATAC_peaks
ARCHR=$WORK/archr
THREADS=${SLURM_CPUS_PER_TASK:-2}
MACS2=$WORK/macs2_wrapper.sh

mkdir -p "$WORK/logs"
cd "$WORK"

# R / ArchR environment
module purge
module load StdEnv/2023
module load r/4.5.0

# ArchR scratch
export TMPDIR=${SLURM_TMPDIR:-/tmp}

# Fail immediately if MACS2 is not callable.
if [[ ! -x "$MACS2" ]]; then
    echo "ERROR: MACS2 wrapper is missing or not executable: $MACS2" >&2
    exit 1
fi

echo "R: $(which Rscript)"
Rscript --version
echo "MACS2 preflight:"
"$MACS2" --version

# Step 03 only: reproducible peaks per Condition x CellType.
Rscript 03_call_peaks.R \
  --archr "$ARCHR" \
  --threads "$THREADS" \
  --minCells 100 \
  --byCondition \
  --conditions "Normal" \
  --macs2 "$MACS2"

#--conditions "Normal,Adenoma"
echo "Step 03 completed successfully."
