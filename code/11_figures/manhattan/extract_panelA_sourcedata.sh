#!/bin/bash
#SBATCH --account=l2_bioinfo1
#SBATCH --job-name=extract_panelA
#SBATCH --error=%x-%j.error
#SBATCH --output=%x-%j.out
#SBATCH --mem=96G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=4:00:00

# =============================================================================
# One-off extraction of panel A's UMAP coordinates for Final_Viz_4panel.R.
#
#   sbatch extract_panelA_sourcedata.sh
#
# Runs extract_panelA_sourcedata.R inside the same RStudio container that
# 04_sub_cluster.sh uses, so the Seurat version that wrote the object is the
# one that reads it.  96G is generous for a 5.76 GB object (readRDS peaks at a
# multiple of the file size); 04_sub_cluster.sh asks for 160G if this ever
# comes back OOM-killed.
#
# The R script writes
#   /data/l2_bioinfo1/liq17/scTWAS/data/VUMC_HTAN_GSE266546/raw/rds2/
#     umap_HTAN_GSE266546_reviewed_celltype_sourcedata.csv
# Copy that CSV to L:/scRNA-seq/VUMC_HTAN_GSE266546/raw/rds2/ on the Windows
# side and re-run Final_Viz_4panel.R, which then draws panel A from it.
# =============================================================================

set -euo pipefail

PROJECT=/data/l2_bioinfo1/liq17/scTWAS/data/VUMC_HTAN_GSE266546/src
SCRIPT_DIR="${SCRIPT_DIR:-$PROJECT}"      # where this .R file lives on ACCRE

cd "$SCRIPT_DIR"

singularity exec \
    --bind /data/l2_bioinfo1:/data/l2_bioinfo1 \
    --bind /panfs/accrepfs.vampire/home/liq17:/panfs/accrepfs.vampire/home/liq17 \
    --bind /cvmfs:/cvmfs \
    /labs/accre_public/singularity/RStudio/4.5.0/R-RStudio-x86-64-v2.sif \
    Rscript "$SCRIPT_DIR/extract_panelA_sourcedata.R"

echo "Done"
