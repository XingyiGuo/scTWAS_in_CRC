#!/bin/bash
#SBATCH --job-name=peer_ct_expr
#SBATCH --error=%x-%j.error
#SBATCH --out=%x-%j.out
#SBATCH --mem=40G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=7-0:0:0
#SBATCH --account=l2_bioinfo1

#conda activate polyfun
python Process_peer_HTAN_GSE266546_fixed.py $1
