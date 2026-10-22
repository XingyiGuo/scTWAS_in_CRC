#!/bin/bash
#SBATCH --job-name=gse201336_step04
#SBATCH --cpus-per-task=2
#SBATCH --mem=100G
#SBATCH --time=2-00:00:00
#SBATCH --output=/data/l2_bioinfo1/scATAC-seq/scATACseq_GSE201336/GSE201336_scATAC_peaks/logs/%x_%j.out
#SBATCH --error=/data/l2_bioinfo1/scATAC-seq/scATACseq_GSE201336/GSE201336_scATAC_peaks/logs/%x_%j.err
#SBATCH --account=l2_bioinfo1

set -euo pipefail

WORK=/data/l2_bioinfo1/scATAC-seq/scATACseq_GSE201336/GSE201336_scATAC_peaks
ARCHR=$WORK/archr
PEAKS=$WORK/peaks
THREADS=${SLURM_CPUS_PER_TASK:-2}

mkdir -p "$WORK/logs"
cd "$WORK"

# R / ArchR environment
module purge
module load StdEnv/2023
module load r/4.5.0

export TMPDIR=${SLURM_TMPDIR:-/tmp}

echo "R: $(which Rscript)"
Rscript --version

# Step 04 only. This uses proj_epithelial_peaks created by Step 03.
Rscript 04_export_peaks.R \
  --archr "$ARCHR" \
  --out "$PEAKS" \
  --threads "$THREADS"

echo "Step 04 completed successfully. Peaks in: $PEAKS"
