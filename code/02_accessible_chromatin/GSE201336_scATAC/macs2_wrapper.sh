#!/bin/bash
set -euo pipefail

# MACS2 was installed in a Python 3.11.2 venv created under StdEnv/2020.
# Keep this environment isolated from the StdEnv/2023 + R/4.5.0 environment
# used by the parent ArchR SLURM job.
module purge
module load StdEnv/2020
module load python/3.11.2

VENV=/data/l2_bioinfo1/scATAC-seq/scATACseq_GSE201336/GSE201336_scATAC_peaks/.venv_macs2
exec "$VENV/bin/macs2" "$@"
