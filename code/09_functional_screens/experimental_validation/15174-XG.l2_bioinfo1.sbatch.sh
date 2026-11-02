#!/bin/bash
#SBATCH --job-name=RNA
#SBATCH --error=%x-%j.error
#SBATCH --output=%x-%j.out
#SBATCH --mem=40G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --time=7-00:00:00
#SBATCH --account=l2_bioinfo1

# Force all SLURM jobs submitted from this Nextflow launcher
# (including nf-* child jobs) to use the l2_bioinfo1 account.
export SBATCH_ACCOUNT=l2_bioinfo1

# Use the classic Nextflow parser for nf-core compatibility.
export NXF_SYNTAX_PARSER=v1

nextflow run nf-core/rnaseq \
    -resume \
    -r 3.26.0 \
    -c /data/l2_bioinfo1/liq17/biosoft/nextflow_tutorials/handson_projects/ref/AA.config \
    -profile singularity \
    --input /data/l2_bioinfo1/RNA-seq/15174-XG/src/15174-XG.samplesheet.csv \
    --outdir /data/l2_bioinfo1/RNA-seq/15174-XG/processed/nf_rnaseq \
    --gtf /data/l2_bioinfo1/liq17/biosoft/nextflow_tutorials/handson_projects/ref/Homo_sapiens.GRCh38.114.gtf.gz \
    --fasta /data/l2_bioinfo1/liq17/biosoft/nextflow_tutorials/handson_projects/ref/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa.gz \
    --save_reference
