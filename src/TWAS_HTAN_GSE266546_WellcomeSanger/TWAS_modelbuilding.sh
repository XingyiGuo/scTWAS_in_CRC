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

#conda activate tftwas

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
#Rscript TWAS_modelbuilding_scTWAS.R <CT_index> <chrom> [*.config.txt]
Rscript "TWAS_modelbuilding_scTWAS.R" "$1" "$2" "$3"
