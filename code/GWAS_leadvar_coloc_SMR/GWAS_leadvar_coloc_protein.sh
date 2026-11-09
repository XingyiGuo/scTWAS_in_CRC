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

#Rscript GWAS_leadvar_coloc_protein.R $1

Rscript GWAS_leadvar_coloc_protein_Polyps_WellcomeSanger.R $1
