#!/bin/bash
#SBATCH --account=l2_bioinfo1
#SBATCH --job-name=04_sub 
#SBATCH --error=%x-%A_%a.error 
#SBATCH --output=%x-%A_%a.out
#SBATCH --mem=160G 
#SBATCH --nodes=1 
#SBATCH --ntasks=1 
#SBATCH --cpus-per-task=1 
#SBATCH --time=7-0:0:0

set -euo pipefail

PROJECT=/data/l2_bioinfo1/liq17/scTWAS/data/VUMC_HTAN_GSE266546/src

cd "$PROJECT"
        
singularity exec \
    --bind /data/l2_bioinfo1:/data/l2_bioinfo1 \
    --bind /panfs/accrepfs.vampire/home/liq17:/panfs/accrepfs.vampire/home/liq17 \
    --bind /cvmfs:/cvmfs \
    /labs/accre_public/singularity/RStudio/4.5.0/R-RStudio-x86-64-v2.sif \
    Rscript "$PROJECT/04_sub_cluster.R"
	
echo "Done"