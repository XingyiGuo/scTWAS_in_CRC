#!/usr/bin/env bash
set -euo pipefail

PYTHON="/home/liq17/venvs/r-h5ad/bin/python"

H5AD="/data/l2_bioinfo1/scRNA-seq/VUMC_HTAN/data/PCA83samples.h5ad"
PCA_META="/data/l2_bioinfo1/scRNA-seq/VUMC_HTAN/data/PCA83samples.metadata.csv"

DEST="/v5000/data/l2_bioinfo1/liq17/scTWAS/data/VUMC_HTAN/VUMC_HTAN_DIS_VAL_ADSER/Results/HTAN_DISVAL_PCA83_ADSER_HPtoSER_AllCT_ASCSSCfiltered_CellOracleR_metadata.csv"

OUT="/v5000/data/l2_bioinfo1/liq17/scTWAS/data/VUMC_HTAN/VUMC_HTAN_DIS_VAL_ADSER/Results/HTAN_DISVAL_PCA83_ADSER_HPtoSER_AllCT_ASCSSCfiltered_CellOracleR_metadata.addHTA11_6407_8890_8930.csv"

echo "============================================================"
echo "1. Generate PCA83samples.metadata.csv"
echo "============================================================"

[[ -s "$H5AD" ]] || {
    echo "ERROR: H5AD not found: $H5AD" >&2
    exit 1
}

"$PYTHON" -u - <<PY
import anndata as ad

H5AD = r"""$H5AD"""
OUT = r"""$PCA_META"""

print("Reading:", H5AD, flush=True)

adata = ad.read_h5ad(H5AD, backed="r")

print("H5AD shape:", adata.shape, flush=True)

obs = adata.obs.copy()

# Preserve cell barcode / obs name
obs.insert(0, "h5ad_barcode", obs.index.astype(str))

obs.to_csv(
    OUT,
    index=False
)

print("Generated:", OUT, flush=True)
print("Metadata rows:", len(obs), flush=True)
print("Metadata columns:", len(obs.columns), flush=True)

adata.file.close()
PY

[[ -s "$PCA_META" ]] || {
    echo "ERROR: metadata CSV was not generated: $PCA_META" >&2
    exit 1
}

echo
echo "============================================================"
echo "2. Extract 3 samples and combine with existing metadata"
echo "============================================================"

[[ -s "$DEST" ]] || {
    echo "ERROR: existing metadata not found: $DEST" >&2
    exit 1
}

"$PYTHON" -u - <<PY
import pandas as pd
import numpy as np
import re
import os

PCA_META = r"""$PCA_META"""
DEST = r"""$DEST"""
OUT = r"""$OUT"""

TARGETS_RAW = [
    "HTA11_6407_200000_1021",
    "HTA11_8890_200000_1011",
    "HTA11_8930_200000_1011",
]

CELL_TYPES = {
    "ABS", "ASC", "CT", "EE", "GOB",
    "SSC", "STM", "TAC", "TUF"
}

def expression_id(x):
    """
    PCA83:
      HTA11_6407_200000_1021

    Expression/combined metadata:
      HTA11_6407_2000001021
    """
    if pd.isna(x):
        return np.nan

    return re.sub(
        r'_(\d{6})_(\d{4})$',
        r'_\1\2',
        str(x).strip()
    )

print("Reading PCA83 metadata:", PCA_META, flush=True)
pca = pd.read_csv(PCA_META, low_memory=False)

required = ["SPECID.htan", "Cell Type"]
for c in required:
    if c not in pca.columns:
        raise RuntimeError(f"Required column missing from PCA83 metadata: {c}")

# ------------------------------------------------------------
# Extract only the 3 target specimens
# ------------------------------------------------------------

new = pca[
    pca["SPECID.htan"]
       .astype(str)
       .isin(TARGETS_RAW)
].copy()

print("\nRaw target counts:", flush=True)

for s in TARGETS_RAW:
    n = (
        new["SPECID.htan"]
        .astype(str)
        .eq(s)
        .sum()
    )
    print(f"  {s}: {n}", flush=True)

missing = [
    s for s in TARGETS_RAW
    if not new["SPECID.htan"].astype(str).eq(s).any()
]

if missing:
    raise RuntimeError(
        "Missing target specimen(s): "
        + ", ".join(missing)
    )

# ------------------------------------------------------------
# Keep the 9 epithelial cell types
# ------------------------------------------------------------

new = new[
    new["Cell Type"]
       .astype(str)
       .isin(CELL_TYPES)
].copy()

print("\nAfter 9-cell-type filtering:", len(new), "rows", flush=True)

print("\nSample x Cell Type:", flush=True)
print(
    pd.crosstab(
        new["SPECID.htan"],
        new["Cell Type"]
    ).to_string(),
    flush=True
)

# ------------------------------------------------------------
# Create combined-metadata fields
# ------------------------------------------------------------

new["orig.ident"] = "PCA83samples"
new["study"] = "HTAN_islam2024"
new["htan_file"] = "PCA83samples"
new["cohort"] = "PCA83"

# Patient ID
if "STUDYID.htan" in new.columns:
    new["PatientID"] = new["STUDYID.htan"].astype(str)
else:
    new["PatientID"] = (
        new["SPECID.htan"]
        .astype(str)
        .str.extract(r'^(HTA11_[^_]+)', expand=False)
    )

# Critical ID conversion for expression files
new["SampleID"] = new["SPECID.htan"].map(expression_id)

# Cell type
new["CellType_cb"] = new["Cell Type"].astype(str)

# Preserve PCA83 pathology fields
if "ADTYPE" in new.columns:
    new["PCA83_ADTYPE"] = new["ADTYPE"]

if "HPTYPE" in new.columns:
    new["PCA83_HPTYPE"] = new["HPTYPE"]

if "TISSUE" in new.columns:
    new["TISSUE_original"] = new["TISSUE"]
    new["sample_condition_original"] = new["TISSUE"]

# ------------------------------------------------------------
# Do NOT force AD/SER for these 3 samples
# ------------------------------------------------------------
# PSHAPE == Sessile alone is not sufficient evidence for SER.

new["sample_condition"] = np.nan
new["condition_mapping_source"] = np.nan
new["Type"] = np.nan

if "sample_condition_source" in new.columns:
    new["sample_condition_source"] = "PCA83_pathology_unclassified"

# Optional Seurat-style QC fields
if "total_counts" in new.columns:
    new["nCount_RNA"] = new["total_counts"]

if "n_genes_by_counts" in new.columns:
    new["nFeature_RNA"] = new["n_genes_by_counts"]

# DIS / VAL-only columns should be empty for PCA83
for c in [
    "HTAN Parent Data File ID",
    "HTAN Specimen ID",
    "Cell_Type",
    "Polyp_Type",
    "Sample_Classification",
]:
    new[c] = np.nan

# ------------------------------------------------------------
# Read existing combined metadata
# ------------------------------------------------------------

print("\nReading existing metadata:", DEST, flush=True)

old = pd.read_csv(
    DEST,
    low_memory=False
)

print("Existing rows:", len(old), flush=True)
print("Existing columns:", len(old.columns), flush=True)

if "SampleID" not in old.columns:
    raise RuntimeError("SampleID missing from existing metadata")

# ------------------------------------------------------------
# Align new rows to exact existing column schema
# ------------------------------------------------------------

for col in old.columns:
    if col not in new.columns:
        new[col] = np.nan

new = new[old.columns].copy()

target_ids = {
    expression_id(x)
    for x in TARGETS_RAW
}

print("\nExpression-format target IDs:", flush=True)
for s in sorted(target_ids):
    print(" ", s, flush=True)

# ------------------------------------------------------------
# Avoid duplicate target rows if script is rerun
# ------------------------------------------------------------

n_existing_targets = (
    old["SampleID"]
    .astype(str)
    .isin(target_ids)
    .sum()
)

print("\nExisting target rows:", n_existing_targets, flush=True)

if n_existing_targets > 0:
    old = old[
        ~old["SampleID"]
          .astype(str)
          .isin(target_ids)
    ].copy()

# ------------------------------------------------------------
# Combine
# ------------------------------------------------------------

combined = pd.concat(
    [old, new],
    ignore_index=True
)

print("\nRows added:", len(new), flush=True)
print("Final rows:", len(combined), flush=True)

# ------------------------------------------------------------
# Validate BEFORE writing
# ------------------------------------------------------------

print("\nValidation:", flush=True)

for s in sorted(target_ids):
    n = (
        combined["SampleID"]
        .astype(str)
        .eq(s)
        .sum()
    )

    print(f"  {s}: {n} rows", flush=True)

    if n == 0:
        raise RuntimeError(
            f"Target sample missing after combine: {s}"
        )

target_rows = combined[
    combined["SampleID"]
    .astype(str)
    .isin(target_ids)
].copy()

print("\nFinal SampleID x CellType_cb:", flush=True)

print(
    pd.crosstab(
        target_rows["SampleID"],
        target_rows["CellType_cb"]
    ).to_string(),
    flush=True
)

check_cols = [
    "SampleID",
    "PatientID",
    "ADTYPE",
    "ATYPIA",
    "PCAD",
    "PCVILLOUS",
    "PSHAPE",
    "PTYPE",
    "SEGMENT",
    "HPTYPE",
    "sample_condition",
]

check_cols = [
    c for c in check_cols
    if c in target_rows.columns
]

print("\nPathology summary:", flush=True)

print(
    target_rows[
        check_cols
    ]
    .drop_duplicates()
    .sort_values("SampleID")
    .to_string(index=False),
    flush=True
)

# ------------------------------------------------------------
# Write NEW combined file
# ------------------------------------------------------------

combined.to_csv(
    OUT,
    index=False
)

print("\nWROTE:", OUT, flush=True)

# ------------------------------------------------------------
# Re-open written file and verify
# ------------------------------------------------------------

verify = pd.read_csv(
    OUT,
    usecols=["SampleID"],
    low_memory=False
)

print("\nWritten-file verification:", flush=True)

for s in sorted(target_ids):
    n = (
        verify["SampleID"]
        .astype(str)
        .eq(s)
        .sum()
    )

    print(f"  {s}: {n} rows", flush=True)

    if n == 0:
        raise RuntimeError(
            f"Target sample missing from written output: {s}"
        )

print("\nSUCCESS", flush=True)
PY

echo
echo "============================================================"
echo "FINAL OUTPUT"
echo "============================================================"

ls -lh "$PCA_META"
ls -lh "$OUT"
