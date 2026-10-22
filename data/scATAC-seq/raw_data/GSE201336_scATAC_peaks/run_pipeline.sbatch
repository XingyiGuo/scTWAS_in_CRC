#!/bin/bash
#SBATCH --job-name=gse201336_peaks
#SBATCH --cpus-per-task=2
#SBATCH --mem=100G
#SBATCH --time=7-00:00:00
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
#SBATCH	--account=l2_bioinfo1

set -euo pipefail
mkdir -p logs

# --- environment -------------------------------------------------------------
# ArchR >= 1.0.2, R >= 4.1, MACS2 (py2.7 or py3), BSgenome.Hsapiens.UCSC.hg38
module load StdEnv/2023
module load r/4.5.0 python/3.9 || true
# pip install --user MACS2==2.2.7.1

# Project root. Override: WORK=/scratch/$USER/gse201336 sbatch run_pipeline.sbatch
# NOTE: fragments + ArchR scratch total ~300 GB. Keep the scripts in Dropbox, but set
# WORK to local/scratch storage -- do not let a sync client touch the working data.
WORK=${WORK:-$PWD}
FRAG=$WORK/fragments
ARCHR=$WORK/archr
LABELS=$WORK/epithelial_celltypes_atac.tsv   # consolidated ATAC cell-type file (downloaded)
PEAKS=$WORK/peaks
THREADS=${SLURM_CPUS_PER_TASK:-16}
MANIFEST=$WORK/GSE201336_atac_manifest.csv

# ArchR writes a lot of temp files; point them at fast local scratch
export TMPDIR=${SLURM_TMPDIR:-/tmp}
cd "$WORK"

# --- 0) manifest + fragments -------------------------------------------------
# Cell-type labels: epithelial_celltypes_atac.tsv is already downloaded next to the scripts.
# python3 build_manifest.py
# bash 00_download_fragments.sh "$FRAG" Normal Unaffected Adenoma

# --- 0b) pre-flight ----------------------------------------------------------
#Rscript check_annotation_join.R --labels "$LABELS"

# --- 1) arrows ---------------------------------------------------------------
Rscript 01_build_arrows.R \
  --fragments "$FRAG" --manifest "$MANIFEST" \
  --out "$ARCHR" --threads "$THREADS" \
  --conditions "Normal,Unaffected,Adenoma"

# --- 2) relabel Becker -> ABS/CT/EE/GOB/STM/TAC (no reference, no ASC) ---------
Rscript 02_relabel.R \
  --archr "$ARCHR" --labels "$LABELS" --threads "$THREADS"

# --- 3) reproducible peak set per Condition x CellType, called SEPARATELY --------
# Normal  -> STM, TAC, ABS, CT, GOB, EE
# Adenoma -> STM, TAC, ABS, CT, GOB, EE
# (--byCondition = each condition's cell types get their own peak set; no pooling.)
Rscript 03_call_peaks.R --archr "$ARCHR" --threads "$THREADS" --minCells 200 \
  --byCondition --conditions "Normal,Adenoma"

# --- 4) export (no comparison; pass --diff to add differential tables) ---------
Rscript 04_export_peaks.R --archr "$ARCHR" --out "$PEAKS" --threads "$THREADS"

echo "Done. Peaks in $PEAKS"
