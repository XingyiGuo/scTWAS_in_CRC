#!/bin/bash
#SBATCH --job-name=merge_vcf
#SBATCH --error=%x-%j.error
#SBATCH --out=%x-%j.out
#SBATCH --mem=160G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=7-0:0:0
#SBATCH --account=l2_bioinfo1

# Submitter for the chr1-22 VCF merge. The job name and the 160G request above
# belong to that merge, not to the sub-clustering lines kept below for reference.
#
# module load StdEnv/2023
# module load r/4.5.0
# Rscript 04_sub_cluster.R

bash g1_HTAN_GSE266546_merge_chr1_22.sh
