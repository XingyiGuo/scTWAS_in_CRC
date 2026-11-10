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

# Submit from this directory to merge the chromosome 1–22 VCFs.
bash g1_HTAN_GSE266546_merge_chr1_22.sh
