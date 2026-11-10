#!/usr/bin/env python3

"""
TWAS + CPTAC-COAD protein differential-abundance annotation

Input:
    TWAS_ADSER_NL_collected_results.csv

Output:
    TWAS_ADSER_NL_CPTAC_COAD_DE_support.csv

Positive CPTAC_log2FC = higher protein abundance in tumor
Negative CPTAC_log2FC = lower protein abundance in tumor

DE support (matched to the scRNA-seq support definition):
    support        : FDR < 0.05, |log2FC| > 0.5, and signs agree
    DE_discordant  : FDR < 0.05, |log2FC| > 0.5, but signs disagree
    NS             : tested but FDR >= 0.05 or |log2FC| <= 0.5
    not_tested     : DE result unavailable
"""

from pathlib import Path
import atexit
from datetime import datetime
import importlib.metadata as importlib_metadata
import os
import platform
import sys
import warnings

import numpy as np
import pandas as pd
import requests
from scipy import stats
from statsmodels.stats.multitest import multipletests

try:
    import cptac
except ImportError as exc:
    raise SystemExit(
        "\nPython package 'cptac' is not installed.\n"
        "Use the pinned CPTAC setup in code/07_orthogonal_support/README.md.\n"
    ) from exc


# =============================================================================
# CPTAC / Zenodo download workaround
# =============================================================================
#
# cptac 1.5.14 can fail when its Zenodo downloader cannot determine file_size
# from the record metadata, leading to:
#
#     TypeError: unsupported operand type(s) for //: 'NoneType' and 'int'
#
# The package already has the correct download URL. The failure occurs only
# because its downloader tries to calculate chunk_size from file_size before
# downloading. This replacement streams the same URL directly and therefore
# does not require file_size.
#
# It is installed as a runtime monkey patch; nothing inside site-packages is
# modified on disk.
# =============================================================================

def install_cptac_streaming_download_patch():
    """Patch cptac's downloader to stream files without requiring file size."""
    from cptac.tools import download_tools

    def robust_get_data(url: str, subfolder: str = "", num_threads: int = 4) -> str:
        del num_threads  # retained only for API compatibility

        destination = Path(download_tools.DATA_DIR) / subfolder
        destination.parent.mkdir(parents=True, exist_ok=True)

        print("\nCPTAC streaming-download workaround")
        print("URL :", url)
        print("DEST:", destination)

        headers = dict(getattr(download_tools, "AUTH_HEADER", {}) or {})
        temporary = destination.with_name(destination.name + ".part")

        try:
            with requests.get(
                url,
                headers=headers,
                stream=True,
                timeout=(30, 600),
                allow_redirects=True,
            ) as response:
                response.raise_for_status()

                total = response.headers.get("content-length")
                total = int(total) if total and total.isdigit() else None

                downloaded = 0
                chunk_bytes = 1024 * 1024

                with open(temporary, "wb") as handle:
                    for chunk in response.iter_content(chunk_size=chunk_bytes):
                        if not chunk:
                            continue

                        handle.write(chunk)
                        downloaded += len(chunk)

                        if total:
                            pct = 100.0 * downloaded / total
                            print(
                                f"\rDownloaded {downloaded / 1e6:.1f} MB "
                                f"({pct:.1f}%)",
                                end="",
                                flush=True,
                            )
                        else:
                            print(
                                f"\rDownloaded {downloaded / 1e6:.1f} MB",
                                end="",
                                flush=True,
                            )

            print()
            os.replace(temporary, destination)

        except Exception:
            if temporary.exists():
                temporary.unlink()
            raise

        return subfolder

    download_tools.get_data = robust_get_data

    print("Installed CPTAC streaming-download workaround.")


# Hide pandas fragmentation warnings generated internally by cptac.
warnings.filterwarnings(
    "ignore",
    category=pd.errors.PerformanceWarning,
    module=r"cptac\..*",
)


# =============================================================================
# 1. Paths
# =============================================================================

# IMPORTANT:
# Keep these as pathlib.Path objects.
# Do not replace them with plain strings if using "/" or ".mkdir()".

BASE_DIR = Path(
    "/data/l2_bioinfo1/liq17/scTWAS/"
    "TWAS_HTAN_GSE266546_WellcomeSanger_res/"
    "SPrediXcan_peakSNP_association"
)

TWAS_FILE = BASE_DIR / "TWAS_ADSER_NL_collected_results.csv"

OUT_DIR = Path(
    "/v5000/data/l2_bioinfo1/liq17/scTWAS/"
    "Functional_res/TCGA_CPTAC_COAD_DE"
)

OUT_DIR.mkdir(
    parents=True,
    exist_ok=True
)

CPTAC_DE_CACHE = OUT_DIR / "CPTAC_COAD_protein_DE_allgenes.csv"

FINAL_OUT = OUT_DIR / "TWAS_ADSER_NL_CPTAC_COAD_DE_support.csv"

LOG_FILE = OUT_DIR / "CPTAC_COAD_TWAS_annotation.log"


# =============================================================================
# Logging and software versions
# =============================================================================

class Tee:
    """Write a stream to both the terminal and the analysis log."""

    def __init__(self, *streams):
        self.streams = streams

    def write(self, data):
        for stream in self.streams:
            stream.write(data)
            stream.flush()
        return len(data)

    def flush(self):
        for stream in self.streams:
            stream.flush()

    def isatty(self):
        return any(getattr(stream, "isatty", lambda: False)() for stream in self.streams)


_ORIGINAL_STDOUT = sys.stdout
_ORIGINAL_STDERR = sys.stderr
_LOG_HANDLE = open(LOG_FILE, "w", encoding="utf-8", buffering=1)

sys.stdout = Tee(_ORIGINAL_STDOUT, _LOG_HANDLE)
sys.stderr = Tee(_ORIGINAL_STDERR, _LOG_HANDLE)


def _close_log():
    """Restore terminal streams and close the log file."""
    try:
        sys.stdout.flush()
        sys.stderr.flush()
    finally:
        sys.stdout = _ORIGINAL_STDOUT
        sys.stderr = _ORIGINAL_STDERR
        _LOG_HANDLE.close()


atexit.register(_close_log)


def package_version(package_name):
    """Return installed package version without failing the analysis."""
    try:
        return importlib_metadata.version(package_name)
    except importlib_metadata.PackageNotFoundError:
        return "not installed"


print("============================================================")
print("CPTAC-COAD TWAS DE SUPPORT ANALYSIS")
print("============================================================")
print("Run started:", datetime.now().astimezone().strftime("%Y-%m-%d %H:%M:%S %Z"))
print("Python version:", sys.version.replace("\n", " "))
print("Platform:", platform.platform())
print("\nPackage versions:")
for package_name in [
    "numpy",
    "pandas",
    "scipy",
    "statsmodels",
    "requests",
    "cptac",
]:
    print(f"  {package_name:<14} {package_version(package_name)}")

print("\nDE support thresholds:")
print("  FDR < 0.05")
print("  |log2FC| > 0.5")
print("  sign(log2FC) == sign(scTWAS Z score)")
print("\nSupport labels:")
print("  support        : passes FDR/effect-size thresholds and direction agrees")
print("  DE_discordant  : passes FDR/effect-size thresholds but direction disagrees")
print("  NS             : tested but does not pass FDR/effect-size thresholds")
print("  not_tested     : DE result unavailable")
print("============================================================")

# Install the existing CPTAC streaming workaround after logging is active.
install_cptac_streaming_download_patch()


# -----------------------------------------------------------------------------
# Path QC
# -----------------------------------------------------------------------------

print("\n============================================================")
print("PATH CONFIGURATION")
print("============================================================")
print("BASE_DIR :", BASE_DIR)
print("TWAS_FILE:", TWAS_FILE)
print("OUT_DIR  :", OUT_DIR)
print("Python Path type:", type(OUT_DIR).__name__)

if not TWAS_FILE.exists():
    raise FileNotFoundError(
        f"\nTWAS input file was not found:\n{TWAS_FILE}\n"
    )


# =============================================================================
# 2. Read TWAS results
# =============================================================================

twas = pd.read_csv(
    TWAS_FILE
)

required_cols = [
    "Dataset",
    "CellType",
    "gene",
    "gene_name",
    "Sig_resource",
    "Z_score",
    "P_value",
]

missing_cols = [
    c
    for c in required_cols
    if c not in twas.columns
]

if missing_cols:
    raise ValueError(
        "Missing required TWAS columns: "
        + ", ".join(missing_cols)
    )

twas["gene"] = (
    twas["gene"]
    .astype(str)
    .str.replace(
        r"\..*$",
        "",
        regex=True
    )
)

twas["gene_name"] = twas["gene_name"].astype("string")

twas["Z_score"] = pd.to_numeric(
    twas["Z_score"],
    errors="coerce"
)

twas["P_value"] = pd.to_numeric(
    twas["P_value"],
    errors="coerce"
)

print("\n============================================================")
print("TWAS INPUT")
print("============================================================")
print("Rows:", len(twas))
print("Unique Ensembl genes:", twas["gene"].nunique(dropna=True))
print("Unique gene symbols:", twas["gene_name"].nunique(dropna=True))


# =============================================================================
# 3. DE support
# =============================================================================

def get_de_support(
    z_score,
    log2fc,
    fdr,
    fdr_cutoff=0.05,
    abs_log2fc_cutoff=0.5,
):
    """
    Classify CPTAC tumor-normal differential abundance support for scTWAS.

    A result is considered DE only when:
      FDR < 0.05 AND |log2FC| > 0.5.
    """
    # No usable DE result for this gene/protein.
    if pd.isna(log2fc) or pd.isna(fdr) or pd.isna(z_score):
        return "not_tested"

    # Tested, but does not pass the DE significance/effect-size thresholds.
    if fdr >= fdr_cutoff or abs(log2fc) <= abs_log2fc_cutoff:
        return "NS"

    # scTWAS Z scores are expected to be directional; guard against exact zero.
    if z_score == 0:
        return "NS"

    # Significant DE with sufficient effect size.
    if np.sign(z_score) == np.sign(log2fc):
        return "support"

    return "DE_discordant"


# =============================================================================
# 4. CPTAC differential abundance
# =============================================================================

if CPTAC_DE_CACHE.exists():

    print("\n============================================================")
    print("Loading cached CPTAC-COAD protein DE results")
    print("============================================================")
    print(CPTAC_DE_CACHE)

    cptac_de = pd.read_csv(
        CPTAC_DE_CACHE
    )

else:

    print("\n============================================================")
    print("Loading CPTAC colon adenocarcinoma")
    print("============================================================")

    coad = cptac.Coad()

    # University of Michigan proteomics dataset.
    prot = coad.get_proteomics(
        "umich"
    ).copy()

    print(
        "Original CPTAC proteomics shape:",
        prot.shape
    )

    # -------------------------------------------------------------------------
    # Convert protein columns to gene symbols.
    #
    # CPTAC can use a MultiIndex such as:
    #   Name / Database_ID
    # -------------------------------------------------------------------------

    if isinstance(
        prot.columns,
        pd.MultiIndex
    ):

        if "Name" in prot.columns.names:

            gene_names = (
                prot.columns
                .get_level_values("Name")
                .astype(str)
            )

        else:

            gene_names = (
                prot.columns
                .get_level_values(0)
                .astype(str)
            )

    else:

        gene_names = pd.Index(
            [
                str(x)
                for x in prot.columns
            ]
        )

    prot.columns = gene_names

    # -------------------------------------------------------------------------
    # Collapse multiple protein features mapping to the same gene symbol.
    # -------------------------------------------------------------------------

    prot = (
        prot
        .T
        .groupby(level=0)
        .mean()
        .T
    )

    prot.index = prot.index.astype(str)

    print(
        "Gene-level CPTAC proteomics shape:",
        prot.shape
    )

    # -------------------------------------------------------------------------
    # Identify tumor and adjacent-normal samples.
    #
    # CPTAC adjacent-normal sample IDs commonly end in ".N".
    # -------------------------------------------------------------------------

    normal_mask = prot.index.str.endswith(
        ".N"
    )

    normal = prot.loc[
        normal_mask,
        :
    ].copy()

    tumor = prot.loc[
        ~normal_mask,
        :
    ].copy()

    print(
        "CPTAC tumor samples:",
        tumor.shape[0]
    )

    print(
        "CPTAC normal samples:",
        normal.shape[0]
    )

    print(
        "CPTAC genes:",
        prot.shape[1]
    )

    if normal.shape[0] == 0:
        raise RuntimeError(
            "No CPTAC normal samples were identified using the '.N' suffix."
        )

    if tumor.shape[0] == 0:
        raise RuntimeError(
            "No CPTAC tumor samples were identified."
        )

    # -------------------------------------------------------------------------
    # Identify matched tumor-normal sample pairs.
    #
    # Tumor : C3L-xxxxx
    # Normal: C3L-xxxxx.N
    # -------------------------------------------------------------------------

    tumor_ids = {
        str(x): x
        for x in tumor.index
    }

    normal_ids = {
        str(x).replace(".N", ""): x
        for x in normal.index
    }

    paired_patient_ids = sorted(
        set(tumor_ids.keys()).intersection(
            set(normal_ids.keys())
        )
    )

    print(
        "Potential paired tumor-normal samples:",
        len(paired_patient_ids)
    )

    # -------------------------------------------------------------------------
    # Differential protein abundance
    #
    # Prefer paired t-test when >=5 valid pairs exist for a protein.
    # Otherwise use Welch's two-sample t-test when each group has >=3 values.
    #
    # CPTAC proteomic values are on a log-scale in this workflow, so
    # Tumor mean - Normal mean is used as CPTAC_log2FC.
    # -------------------------------------------------------------------------

    results = []

    for gene in prot.columns:

        tumor_values = pd.to_numeric(
            tumor[gene],
            errors="coerce"
        ).dropna()

        normal_values = pd.to_numeric(
            normal[gene],
            errors="coerce"
        ).dropna()

        paired_tumor = []
        paired_normal = []

        for patient in paired_patient_ids:

            tumor_sample = tumor_ids[
                patient
            ]

            normal_sample = normal_ids[
                patient
            ]

            tumor_value = tumor.at[
                tumor_sample,
                gene
            ]

            normal_value = normal.at[
                normal_sample,
                gene
            ]

            if (
                pd.notna(tumor_value)
                and
                pd.notna(normal_value)
            ):

                paired_tumor.append(
                    float(tumor_value)
                )

                paired_normal.append(
                    float(normal_value)
                )

        paired_tumor = np.asarray(
            paired_tumor,
            dtype=float
        )

        paired_normal = np.asarray(
            paired_normal,
            dtype=float
        )

        n_pairs = len(
            paired_tumor
        )

        # ---------------------------------------------------------------------
        # Paired analysis
        # ---------------------------------------------------------------------

        if n_pairs >= 5:

            log2fc = np.nanmean(
                paired_tumor
                -
                paired_normal
            )

            test_result = stats.ttest_rel(
                paired_tumor,
                paired_normal,
                nan_policy="omit"
            )

            pvalue = test_result.pvalue

        # ---------------------------------------------------------------------
        # Unpaired Welch t-test
        # ---------------------------------------------------------------------

        elif (
            len(tumor_values) >= 3
            and
            len(normal_values) >= 3
        ):

            log2fc = (
                tumor_values.mean()
                -
                normal_values.mean()
            )

            test_result = stats.ttest_ind(
                tumor_values,
                normal_values,
                equal_var=False,
                nan_policy="omit"
            )

            pvalue = test_result.pvalue

        # ---------------------------------------------------------------------
        # Insufficient observations
        # ---------------------------------------------------------------------

        else:

            log2fc = np.nan
            pvalue = np.nan

        results.append(
            {
                "gene_name": str(gene),
                "CPTAC_log2FC": log2fc,
                "CPTAC_P": pvalue,
            }
        )

    cptac_de = pd.DataFrame(
        results
    )

    # -------------------------------------------------------------------------
    # BH FDR across proteins with valid P values.
    # -------------------------------------------------------------------------

    valid = (
        cptac_de["CPTAC_P"].notna()
        &
        np.isfinite(
            cptac_de["CPTAC_P"]
        )
    )

    cptac_de["CPTAC_FDR"] = np.nan

    if valid.sum() > 0:

        adjusted = multipletests(
            cptac_de.loc[
                valid,
                "CPTAC_P"
            ].values,
            method="fdr_bh"
        )[1]

        cptac_de.loc[
            valid,
            "CPTAC_FDR"
        ] = adjusted

    cptac_de.to_csv(
        CPTAC_DE_CACHE,
        index=False
    )

    print(
        "Proteins with valid DE test:",
        int(valid.sum())
    )

    print(
        "Saved CPTAC DE cache:",
        CPTAC_DE_CACHE
    )


# =============================================================================
# 5. Validate CPTAC table
# =============================================================================

required_cptac_cols = [
    "gene_name",
    "CPTAC_log2FC",
    "CPTAC_P",
    "CPTAC_FDR",
]

missing_cptac_cols = [
    c
    for c in required_cptac_cols
    if c not in cptac_de.columns
]

if missing_cptac_cols:
    raise ValueError(
        "CPTAC DE table is missing required columns: "
        + ", ".join(missing_cptac_cols)
    )

cptac_de["gene_name"] = (
    cptac_de["gene_name"]
    .astype("string")
)

for col in [
    "CPTAC_log2FC",
    "CPTAC_P",
    "CPTAC_FDR",
]:
    cptac_de[col] = pd.to_numeric(
        cptac_de[col],
        errors="coerce"
    )

# Make sure one gene symbol occurs only once before many-to-one merge.
if cptac_de["gene_name"].duplicated().any():

    duplicated_symbols = (
        cptac_de.loc[
            cptac_de["gene_name"].duplicated(
                keep=False
            ),
            "gene_name"
        ]
        .dropna()
        .unique()
    )

    print(
        "\nWARNING: duplicated CPTAC gene symbols found:",
        len(duplicated_symbols)
    )

    # Collapse duplicates defensively.
    cptac_de = (
        cptac_de
        .groupby(
            "gene_name",
            as_index=False,
            dropna=False
        )
        .agg(
            CPTAC_log2FC=("CPTAC_log2FC", "mean"),
            CPTAC_P=("CPTAC_P", "min"),
            CPTAC_FDR=("CPTAC_FDR", "min"),
        )
    )

print(
    "\nCPTAC genes with protein data:",
    len(cptac_de)
)


# =============================================================================
# 6. Annotate TWAS
# =============================================================================

final = twas.merge(
    cptac_de[
        [
            "gene_name",
            "CPTAC_log2FC",
            "CPTAC_P",
            "CPTAC_FDR",
        ]
    ],
    how="left",
    on="gene_name",
    validate="many_to_one"
)

final["CPTAC_DE_support"] = [
    get_de_support(
        z,
        fc,
        fdr
    )
    for z, fc, fdr in zip(
        final["Z_score"],
        final["CPTAC_log2FC"],
        final["CPTAC_FDR"]
    )
]

final = final[
    [
        "Dataset",
        "CellType",
        "gene",
        "gene_name",
        "Sig_resource",
        "Z_score",
        "P_value",
        "CPTAC_log2FC",
        "CPTAC_P",
        "CPTAC_FDR",
        "CPTAC_DE_support",
    ]
].sort_values(
    [
        "Dataset",
        "CellType",
        "P_value",
    ],
    kind="stable"
)


# =============================================================================
# 7. QC + save
# =============================================================================

expected_cols = [
    "Dataset",
    "CellType",
    "gene",
    "gene_name",
    "Sig_resource",
    "Z_score",
    "P_value",
    "CPTAC_log2FC",
    "CPTAC_P",
    "CPTAC_FDR",
    "CPTAC_DE_support",
]

if list(final.columns) != expected_cols:
    raise RuntimeError(
        "Final CPTAC column structure does not match expected output."
    )

print("\n============================================================")
print("CPTAC DE SUPPORT SUMMARY")
print("============================================================")

print(
    final["CPTAC_DE_support"]
    .value_counts(
        dropna=False
    )
)

print(
    "\nTWAS rows with CPTAC protein data:",
    final["CPTAC_log2FC"].notna().sum(),
    "/",
    len(final)
)

print(
    "TWAS rows supported by CPTAC:",
    (
        final["CPTAC_DE_support"]
        ==
        "support"
    ).sum()
)


# =============================================================================
# 8. Save final output
# =============================================================================

final.to_csv(
    FINAL_OUT,
    index=False
)

print("\n============================================================")
print("CPTAC ANNOTATED FILE")
print("============================================================")
print(FINAL_OUT)

print("\nCPTAC support labels:")
print(final["CPTAC_DE_support"].value_counts(dropna=False))

print("\nRun completed:", datetime.now().astimezone().strftime("%Y-%m-%d %H:%M:%S %Z"))
print("Log file:", LOG_FILE)
