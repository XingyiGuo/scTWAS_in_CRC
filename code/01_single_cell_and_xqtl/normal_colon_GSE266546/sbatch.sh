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

# module load StdEnv/2023
# module load r/4.5.0
# Rscript 04_sub_cluster.R
# Rscript 04_sub_cluster_ori.R

bash g1_merge_HTAN_GSE266546_chr1_22_fixed_v2.sh
