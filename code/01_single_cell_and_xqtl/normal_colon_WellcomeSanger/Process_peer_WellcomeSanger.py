#!/usr/bin/env python3
"""
process_peer_WellcomeSanger.py

Run PEER for one Wellcome Sanger epithelial cell type.

Usage
-----
python process_peer_WellcomeSanger.py ABS
python process_peer_WellcomeSanger.py GOB
python process_peer_WellcomeSanger.py STM

cell_type is used to determine both the expression input file and output names.

Expression input
----------------
<processed_expression>/<cell_type>.TPM.log2.QN.INV.csv.gz

Rows    = genes
Columns = subjects
First column = gene

The expression input is already:
    TPM -> log2(TPM+1) -> QN -> INV

Therefore this script DOES NOT apply QN or INV again.

Covariates
----------
Matched by:
    expression column sample ID
        ==
    covariate individual_id

PEER known covariates:
    sex, age, PC1, PC2, PC3, PC4, PC5
"""

import argparse
import subprocess
from pathlib import Path

import numpy as np
import pandas as pd


# ==============================================================================
# Fixed settings
# ==============================================================================

N_ITERATIONS = 1000

PROCESSED_EXPRESSION_DIR = Path(
    "/data/l2_bioinfo1/liq17/scTWAS/data/WellcomeSanger/expression/"
    "Healthy_TI_then_rectum_278samples_outputs_LauLabel/"
    "pseudobulk_Epithelial_LauLabel_x_individual/"
    "processed_expression"
)

COV_FILE = Path(
    "/data/l2_bioinfo1/liq17/scTWAS/data/WellcomeSanger/genotype/"
    "Healthy_TI_Rectum_N278_age_sex_tissue_geno5PCs.csv"
)

ID_COL = "individual_id"

COVARIATE_COLS = [
    "sex",
    "age",
    "PC1",
    "PC2",
    "PC3",
    "PC4",
    "PC5",
]

VALID_CELL_TYPES = [
    "ABS",
    "CT",
    "EE",
    "GOB",
    "STM",
    "TAC",
    "TUF",
]


# ==============================================================================
# Arguments
# ==============================================================================

parser = argparse.ArgumentParser(
    description="Run PEER for one Wellcome Sanger epithelial cell type."
)

parser.add_argument(
    "cell_type",
    choices=VALID_CELL_TYPES,
    help="Cell type, e.g. ABS, CT, EE, GOB, STM, TAC, TUF",
)

parser.add_argument(
    "--iterations",
    type=int,
    default=N_ITERATIONS,
    help=f"Number of PEER iterations (default: {N_ITERATIONS})",
)

args = parser.parse_args()

cell_type = args.cell_type
n_iterations = args.iterations


# ==============================================================================
# Cell-type-specific input/output
# ==============================================================================

EXPR_FILE = (
    PROCESSED_EXPRESSION_DIR
    / f"{cell_type}.TPM.log2.QN.INV.csv.gz"
)

OUT_DIR = (
    PROCESSED_EXPRESSION_DIR
    / f"PEER_{cell_type}"
)

OUT_DIR.mkdir(parents=True, exist_ok=True)


# ==============================================================================
# Helpers
# ==============================================================================

def peer_factor_count(n_samples: int) -> int:
    if n_samples < 150:
        return 15
    elif n_samples < 250:
        return 30
    elif n_samples < 350:
        return 45
    else:
        return 60


def orient_residuals(raw, genes, samples):
    """Return residuals in genes x samples orientation."""
    g = len(genes)
    n = len(samples)

    if raw.shape == (g, n):
        out = raw.copy()
    elif raw.shape == (n, g):
        out = raw.T
    else:
        raise ValueError(
            f"Unexpected residuals.csv shape {raw.shape}; "
            f"expected {(g, n)} or {(n, g)}."
        )

    out.index = genes
    out.columns = samples
    out.index.name = "gene"
    return out


def orient_x(raw, samples, covariates, n_peer):
    """Orient X as samples x factors when possible."""
    n = len(samples)

    if raw.shape[0] == n:
        x = raw.copy()
    elif raw.shape[1] == n:
        x = raw.T
    else:
        return raw

    x.index = samples
    x.index.name = "individual_id"

    k = x.shape[1]

    if k == n_peer:
        x.columns = [f"PEER{i+1}" for i in range(n_peer)]
    elif k == len(covariates) + n_peer:
        x.columns = list(covariates) + [
            f"PEER{i+1}" for i in range(n_peer)
        ]
    else:
        x.columns = [f"Factor{i+1}" for i in range(k)]

    return x


def orient_w(raw, genes):
    """Orient W as genes x loadings when possible."""
    g = len(genes)

    if raw.shape[0] == g:
        w = raw.copy()
    elif raw.shape[1] == g:
        w = raw.T
    else:
        return raw

    w.index = genes
    w.index.name = "gene"
    w.columns = [f"Loading{i+1}" for i in range(w.shape[1])]
    return w


# ==============================================================================
# 1. Check inputs
# ==============================================================================

if not EXPR_FILE.exists():
    raise FileNotFoundError(
        f"Expression file not found for {cell_type}: {EXPR_FILE}"
    )

if not COV_FILE.exists():
    raise FileNotFoundError(
        f"Covariate file not found: {COV_FILE}"
    )

print("=" * 72)
print(f"PEER: Wellcome Sanger {cell_type}")
print("=" * 72)
print(f"Cell type   : {cell_type}")
print(f"Expression  : {EXPR_FILE}")
print(f"Covariates : {COV_FILE}")
print(f"Output      : {OUT_DIR}")
print()


# ==============================================================================
# 2. Expression: genes x subjects
# ==============================================================================

expr = pd.read_csv(EXPR_FILE)

if "gene" not in expr.columns:
    raise ValueError(
        "Expression file must contain a first column named 'gene'."
    )

expr["gene"] = expr["gene"].astype(str)

if expr["gene"].duplicated().any():
    dup = (
        expr.loc[expr["gene"].duplicated(keep=False), "gene"]
        .drop_duplicates()
    )
    raise ValueError(
        f"Duplicated genes detected. Examples: {dup.head().tolist()}"
    )

expr = expr.set_index("gene")
expr.columns = expr.columns.astype(str)
expr = expr.apply(pd.to_numeric, errors="raise")

if not np.isfinite(expr.to_numpy(dtype=float)).all():
    raise ValueError("Expression matrix contains NaN/Inf.")

print(f"Expression genes    : {expr.shape[0]:,}")
print(f"Expression subjects : {expr.shape[1]:,}")


# ==============================================================================
# 3. Covariates and individual_id matching
# ==============================================================================

cov_all = pd.read_csv(COV_FILE)

required = {ID_COL, *COVARIATE_COLS}
missing = required.difference(cov_all.columns)

if missing:
    raise ValueError(
        "Missing covariate columns: "
        + ", ".join(sorted(missing))
    )

cov_all[ID_COL] = cov_all[ID_COL].astype(str)

if cov_all[ID_COL].duplicated().any():
    raise ValueError(
        f"Duplicated {ID_COL} values detected in covariate file."
    )

cov_all = cov_all.set_index(ID_COL)

expr_samples = expr.columns.tolist()
cov_ids = set(cov_all.index)

# Preserve expression sample order.
common_samples = [
    s for s in expr_samples
    if s in cov_ids
]

expr_only = [
    s for s in expr_samples
    if s not in cov_ids
]

cov_only = [
    s for s in cov_all.index
    if s not in set(expr_samples)
]

print(f"Covariate subjects  : {cov_all.shape[0]:,}")
print(f"Matched subjects    : {len(common_samples):,}")
print(f"Expression only     : {len(expr_only):,}")
print(f"Covariate only      : {len(cov_only):,}")

pd.Series(
    expr_only,
    name=ID_COL,
).to_csv(
    OUT_DIR / f"{cell_type}.expression_only_subjects.csv",
    index=False,
)

pd.Series(
    cov_only,
    name=ID_COL,
).to_csv(
    OUT_DIR / f"{cell_type}.covariate_only_subjects.csv",
    index=False,
)

if len(common_samples) == 0:
    raise RuntimeError(
        "No overlapping subject IDs between expression and covariates."
    )

if expr_only:
    raise RuntimeError(
        f"{len(expr_only)} expression subjects have no covariates."
    )


# ==============================================================================
# 4. Matched PEER inputs
# ==============================================================================

expr_common = expr.loc[:, common_samples]

cov_common = cov_all.loc[
    common_samples,
    COVARIATE_COLS
].copy()

cov_common = cov_common.apply(
    pd.to_numeric,
    errors="raise",
)

if cov_common.isna().any().any():
    bad = cov_common.columns[
        cov_common.isna().any()
    ].tolist()
    raise ValueError(
        f"Missing covariate values detected in: {bad}"
    )

if not np.isfinite(
    cov_common.to_numpy(dtype=float)
).all():
    raise ValueError("Covariates contain NaN/Inf.")


# PEER: samples x genes
expr_NxG = expr_common.T

if expr_NxG.index.tolist() != cov_common.index.tolist():
    raise RuntimeError(
        "Expression/covariate sample order mismatch."
    )

n_samples = expr_NxG.shape[0]
n_genes = expr_NxG.shape[1]
n_peer_factors = peer_factor_count(n_samples)

print()
print(f"PEER samples       : {n_samples}")
print(f"PEER genes         : {n_genes}")
print(f"Known covariates   : {', '.join(COVARIATE_COLS)}")
print(f"PEER factors       : {n_peer_factors}")
print(f"PEER iterations    : {n_iterations}")


# ==============================================================================
# 5. Write PEER inputs
# ==============================================================================

PEER_EXPR_FILE = (
    OUT_DIR
    / f"{cell_type}.TPM.log2.QN.INV.peertool.expression.csv"
)

PEER_COV_FILE = (
    OUT_DIR
    / f"{cell_type}.age_sex_geno5PCs.peertool.covariates.csv"
)

# Numeric-only input for peertool.
expr_NxG.to_csv(
    PEER_EXPR_FILE,
    index=False,
    header=False,
)

cov_common.to_csv(
    PEER_COV_FILE,
    index=False,
    header=False,
)


# Save labeled copies for auditing.
expr_NxG.to_csv(
    OUT_DIR
    / f"{cell_type}.TPM.log2.QN.INV.peertool.expression.labeled.csv"
)

cov_all.loc[
    common_samples
].reset_index().to_csv(
    OUT_DIR
    / f"{cell_type}.age_sex_geno5PCs.matched.labeled.csv",
    index=False,
)

pd.Series(
    common_samples,
    name=ID_COL,
).to_csv(
    OUT_DIR / f"{cell_type}.PEER_samples.csv",
    index=False,
)

pd.Series(
    expr_NxG.columns,
    name="gene",
).to_csv(
    OUT_DIR / f"{cell_type}.PEER_genes.csv",
    index=False,
)


# ==============================================================================
# 6. Run PEER
# ==============================================================================

PEER_OUT_DIR = OUT_DIR / "PEER_outputs"
PEER_OUT_DIR.mkdir(
    parents=True,
    exist_ok=True,
)

peer_cmd = [
    "peertool",
    "-f", str(PEER_EXPR_FILE),
    "-c", str(PEER_COV_FILE),
    "-n", str(n_peer_factors),
    "-i", str(n_iterations),
    "-o", str(PEER_OUT_DIR),
]

print()
print("Running:")
print(" ".join(peer_cmd))

subprocess.run(
    peer_cmd,
    check=True,
)


# ==============================================================================
# 7. Label PEER outputs
# ==============================================================================

gene_names = expr_NxG.columns.tolist()
sample_names = expr_NxG.index.tolist()

residuals_file = PEER_OUT_DIR / "residuals.csv"
x_file = PEER_OUT_DIR / "X.csv"
w_file = PEER_OUT_DIR / "W.csv"


# Residual expression: genes x samples
if residuals_file.exists():
    resid_raw = pd.read_csv(
        residuals_file,
        header=None,
    )

    resid = orient_residuals(
        resid_raw,
        genes=gene_names,
        samples=sample_names,
    )

    resid.to_csv(
        OUT_DIR
        / f"{cell_type}.TPM.log2.QN.INV.PEER_residuals.csv.gz",
        compression="gzip",
    )

    print(
        f"Residuals: {resid.shape} -> genes x samples"
    )


# X
if x_file.exists():
    x_raw = pd.read_csv(
        x_file,
        header=None,
    )

    x = orient_x(
        x_raw,
        samples=sample_names,
        covariates=COVARIATE_COLS,
        n_peer=n_peer_factors,
    )

    x.to_csv(
        OUT_DIR
        / f"{cell_type}.PEER_X.labeled.csv"
    )

    print(f"X: {x.shape}")


# W
if w_file.exists():
    w_raw = pd.read_csv(
        w_file,
        header=None,
    )

    w = orient_w(
        w_raw,
        genes=gene_names,
    )

    w.to_csv(
        OUT_DIR
        / f"{cell_type}.PEER_W.labeled.csv"
    )

    print(f"W: {w.shape}")


# ==============================================================================
# 8. Summary
# ==============================================================================

SUMMARY_FILE = (
    OUT_DIR
    / f"{cell_type}.PEER_summary.txt"
)

with open(SUMMARY_FILE, "w") as fh:
    fh.write(f"cell_type={cell_type}\n")
    fh.write(f"expression_file={EXPR_FILE}\n")
    fh.write(f"covariate_file={COV_FILE}\n")
    fh.write(f"id_column={ID_COL}\n")
    fh.write(f"n_expression_genes={expr.shape[0]}\n")
    fh.write(f"n_expression_samples={expr.shape[1]}\n")
    fh.write(f"n_covariate_samples={cov_all.shape[0]}\n")
    fh.write(f"n_common_samples={len(common_samples)}\n")
    fh.write(f"n_expression_only={len(expr_only)}\n")
    fh.write(f"n_covariate_only={len(cov_only)}\n")
    fh.write(
        f"covariates={','.join(COVARIATE_COLS)}\n"
    )
    fh.write(
        f"n_peer_factors={n_peer_factors}\n"
    )
    fh.write(
        f"n_iterations={n_iterations}\n"
    )
    fh.write(
        "input_normalization=TPM->log2(TPM+1)->QN->INV\n"
    )
    fh.write(
        "additional_INV_before_PEER=False\n"
    )

print()
print("=" * 72)
print("DONE")
print("=" * 72)
print(f"Cell type : {cell_type}")
print(f"Output    : {OUT_DIR}")
print(f"Summary   : {SUMMARY_FILE}")
