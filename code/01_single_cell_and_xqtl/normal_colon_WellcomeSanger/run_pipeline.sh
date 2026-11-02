#!/usr/bin/env bash
set -euo pipefail

# ==============================================================================
# Wellcome Sanger E-MTAB-16986
# H5AD -> BPCells -> Seurat
#
# Run 00_install_BPCell.sh once before the first pipeline run.
# ==============================================================================

SCRIPT_DIR="/data/l2_bioinfo1/liq17/scTWAS/src/WellcomeSanger"

R_SCRIPT="${SCRIPT_DIR}/01_h5ad_to_seurat_bpcells.R"

SUMMARY_DIR="/v5000/data/l2_bioinfo1/scRNA-seq/WellcomeSanger/processed/E-MTAB-16986/multi_tissue_predicted_summary"

module load StdEnv/2023
module load gcc/13.3
module load hdf5/1.14.5
module load r/4.5.0

mkdir -p "${SUMMARY_DIR}"

echo "============================================================"
echo "Wellcome Sanger H5AD -> Seurat pipeline"
echo "============================================================"
echo "Start time: $(date)"
echo "Host      : $(hostname)"
echo "R         : $(which R)"
R --version | head -1
echo

# Save the exact module environment for reproducibility.
module list 2> "${SUMMARY_DIR}/modules_h5ad_to_seurat.txt"

if [[ ! -f "${R_SCRIPT}" ]]; then
    echo "ERROR: R script not found:"
    echo "  ${R_SCRIPT}"
    exit 1
fi

Rscript "${R_SCRIPT}"

echo
echo "============================================================"
echo "Pipeline completed successfully"
echo "============================================================"
echo "End time: $(date)"
