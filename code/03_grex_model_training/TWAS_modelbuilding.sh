#!/bin/bash
#SBATCH --job-name=TWAS
#SBATCH --error=%x-%j.error
#SBATCH --out=%x-%j.out
#SBATCH --mem=10G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=7-0:0:0
#SBATCH --account=l2_bioinfo1

set -euo pipefail

if [[ $# -lt 2 || $# -gt 3 ]]; then
    echo "Usage: $0 <CT_index> <chrom> [config]" >&2
    exit 2
fi

SCRIPT_DIR=${TWAS_SCRIPT_DIR:-${SLURM_SUBMIT_DIR:-$(dirname "${BASH_SOURCE[0]}")}}
TRAINER="$SCRIPT_DIR/TWAS_modelbuilding_scTWAS.R"
if [[ ! -f "$TRAINER" ]]; then
    echo "Trainer not found at $TRAINER. Submit from the training directory or set TWAS_SCRIPT_DIR." >&2
    exit 1
fi

exec Rscript "$TRAINER" "$@"
