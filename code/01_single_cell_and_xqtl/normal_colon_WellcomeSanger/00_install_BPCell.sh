#!/usr/bin/env bash
set -euo pipefail

# One-time installation of BPCells and required workflow packages
# into the existing renv project.

RENV_PROJECT="/data/l2_bioinfo1/liq17/scTWAS/data/VUMC_HTAN_GSE266546/src"

module load StdEnv/2023
module load gcc/13.3
module load hdf5/1.14.5
module load r/4.5.0

echo "============================================================"
echo "HPC environment"
echo "============================================================"
module list
echo

if [[ ! -f "${RENV_PROJECT}/renv/activate.R" ]]; then
    echo "ERROR: renv activation file not found:"
    echo "  ${RENV_PROJECT}/renv/activate.R"
    exit 1
fi

cd "${RENV_PROJECT}"

Rscript --vanilla - <<'RSCRIPT'
source(
    "/data/l2_bioinfo1/liq17/scTWAS/data/VUMC_HTAN_GSE266546/src/renv/activate.R"
)

# BPCells is installed from R-universe to avoid GitHub API authentication issues.
options(
    repos = c(
        BPCells = "https://bnprks.r-universe.dev",
        SatijaLab = "https://satijalab.r-universe.dev",
        CRAN = "https://cloud.r-project.org"
    )
)

needed <- c("BPCells", "Azimuth", "dplyr")

missing <- needed[
    !vapply(
        needed,
        requireNamespace,
        logical(1),
        quietly = TRUE
    )
]

if (length(missing) > 0) {
    message(
        "Installing missing packages: ",
        paste(missing, collapse = ", ")
    )
    install.packages(
        missing,
        dependencies = NA
    )
} else {
    message("All required workflow packages are already installed.")
}

if (!requireNamespace("Seurat", quietly = TRUE)) {
    stop(
        "Seurat is not installed in the active renv project. ",
        "Restore/install Seurat before running the conversion."
    )
}

cat("\nInstalled package versions:\n")
for (pkg in c("Seurat", "BPCells", "Azimuth", "dplyr")) {
    cat(
        sprintf(
            "%-10s %s\n",
            pkg,
            as.character(packageVersion(pkg))
        )
    )
}

renv::snapshot(prompt = FALSE)

cat("\nInstallation completed successfully.\n")
RSCRIPT
