#!/bin/bash
#SBATCH --job-name=coloc
#SBATCH --error=%x-%j.error
#SBATCH --out=%x-%j.out
#SBATCH --mem=60G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=7-0:0:0
#SBATCH --account=l2_bioinfo1

#conda activate tftwas_coloc

# usage: sbatch GWAS_leadvar_coloc_protein.sh <Polyps|VUMC|WellcomeSanger> <cisprotein|transprotein>
#
# The cohort is an argument rather than part of this file, because there is one
# wrapper and three cohort scripts. It used to name a combined
# *_Polyps_WellcomeSanger.R that does not exist, so the job failed immediately.

COHORT="${1:?give a cohort: Polyps, VUMC or WellcomeSanger}"
ANALYSIS="${2:?give an analysis: cisprotein or transprotein}"

SCRIPT="GWAS_leadvar_coloc_protein_${COHORT}.R"
[ -f "$SCRIPT" ] || { echo "no such cohort script: $SCRIPT" >&2; exit 1; }

Rscript "$SCRIPT" "$ANALYSIS"
