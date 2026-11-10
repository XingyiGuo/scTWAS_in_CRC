#!/bin/bash
#SBATCH --job-name=coloc
#SBATCH --error=%x-%j.error
#SBATCH --out=%x-%j.out
#SBATCH --mem=40G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --time=7-0:0:0
#SBATCH --account=l2_bioinfo1

set -euo pipefail

if [[ $# -ne 2 ]]; then
    echo "Usage: $0 <CT|all> <config>" >&2
    exit 2
fi

SCRIPT_DIR=${COLOC_SMR_SCRIPT_DIR:-${SLURM_SUBMIT_DIR:-$(dirname "${BASH_SOURCE[0]}")}}
exec Rscript "$SCRIPT_DIR/01_coloc_gwas_xqtl_sigSource_peakSNP.R" "$@"
