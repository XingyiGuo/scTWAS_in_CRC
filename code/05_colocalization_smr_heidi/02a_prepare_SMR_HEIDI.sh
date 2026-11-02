#!/bin/bash
#SBATCH --job-name=pre_smr
#SBATCH --error=%x-%j.error
#SBATCH --out=%x-%j.out
#SBATCH --mem=40G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --time=7-0:0:0
#SBATCH --account=l2_bioinfo1

#	conda activate tftwas_coloc
# $1 <cellType> ABS
# $2 <config.txt> /data/l2_bioinfo1/liq17/scTWAS/src/coloc_SMRHEIDI_HTAN_GSE266546_WellcomeSanger/SMR_HEIDI_NL_targetTopEQTL_500kb.config.txt
Rscript /data/l2_bioinfo1/liq17/scTWAS/src/coloc_SMRHEIDI_HTAN_GSE266546_WellcomeSanger/02_SMR_HEIDI_inputs_targetTopEQTL_500kb_directTWAS.R $1 $2