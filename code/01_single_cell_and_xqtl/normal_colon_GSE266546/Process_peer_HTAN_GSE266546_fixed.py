#!/usr/bin/env python3
"""
Process_peer_HTAN_GSE266546.py

Run PEER for one epithelial cell type in the HTAN_GSE266546 project.

Usage
-----
python Process_peer_HTAN_GSE266546.py ABS
python Process_peer_HTAN_GSE266546.py GOB
python Process_peer_HTAN_GSE266546.py STM

Expression input
----------------
/data/l2_bioinfo1/liq17/scTWAS/data/VUMC_HTAN_GSE266546/expression/
pseudobulk_expression_for_peer/<cell_type>/
pseudobulk_counts_<cell_type>.tpmlog2.qc.qn.rmlowcounts.inv.csv

Rows    = genes
Columns = expression samples
First column = gene

The expression input is already processed and inverse-normal transformed.
This script DOES NOT apply QN or INV again.

Covariates
----------
/v5000/data/l2_bioinfo1/liq17/scTWAS/data/VUMC_HTAN_GSE266546/genotype/
HTAN_GSE266546_N103_age_sex_tissue_geno5PCs.csv

Samples are matched by:
    expression column sample ID == covariate expression_id

Matched expression samples are reordered to the covariate/genotype row order.
A QTL-ready phenotype copy is written with columns renamed from expression_id
into genotype_id.

PEER known covariates:
    Age
    Sex
    segment_Descending colon
    segment_Sigmoid colon
    segment_Terminal ileum
    PC1-PC5

Note: genotype_id is retained in labeled audit output, but is NOT used as a PEER
covariate.
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

EXPRESSION_ROOT = Path(
    "/data/l2_bioinfo1/liq17/scTWAS/data/VUMC_HTAN_GSE266546/"
    "expression/pseudobulk_expression_for_peer"
)

COV_FILE = Path(
    "/v5000/data/l2_bioinfo1/liq17/scTWAS/data/VUMC_HTAN_GSE266546/"
    "genotype/HTAN_GSE266546_N103_age_sex_tissue_geno5PCs.csv"
)

ID_COL = "expression_id"
GENO_ID_COL = "genotype_id"

COVARIATE_COLS = [
    "Age",
    "Sex",
    "segment_Descending colon",
    "segment_Sigmoid colon",
    "segment_Terminal ileum",
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
    description="Run PEER for one HTAN_GSE266546 epithelial cell type."
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

CELL_DIR = EXPRESSION_ROOT / cell_type

EXPR_FILE = (
    CELL_DIR
    / f"pseudobulk_counts_{cell_type}.tpmlog2.qc.qn.rmlowcounts.inv.csv"
)

OUT_DIR = CELL_DIR / f"PEER_{cell_type}"
OUT_DIR.mkdir(parents=True, exist_ok=True)

PREFIX = f"pseudobulk_counts_{cell_type}.tpmlog2.qc.qn.rmlowcounts.inv"


# ==============================================================================
# Helpers
# ==============================================================================

def peer_factor_count(n_samples: int) -> int:
    """PEER factor count using the same rule as the original script."""
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
    x.index.name = ID_COL

    k = x.shape[1]

    if k == n_peer:
        x.columns = [f"PEER{i + 1}" for i in range(n_peer)]
    elif k == len(covariates) + n_peer:
        x.columns = list(covariates) + [
            f"PEER{i + 1}" for i in range(n_peer)
        ]
    else:
        x.columns = [f"Factor{i + 1}" for i in range(k)]

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
    w.columns = [f"Loading{i + 1}" for i in range(w.shape[1])]
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
print(f"PEER: HTAN_GSE266546 {cell_type}")
print("=" * 72)
print(f"Cell type   : {cell_type}")
print(f"Expression  : {EXPR_FILE}")
print(f"Covariates  : {COV_FILE}")
print(f"Match ID    : {ID_COL}")
print(f"Output      : {OUT_DIR}")
print()


# ==============================================================================
# 2. Expression: genes x samples
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
print(f"Expression samples  : {expr.shape[1]:,}")


# ==============================================================================
# 3. Covariates and expression_id matching
# ==============================================================================

cov_all = pd.read_csv(COV_FILE)

required = {ID_COL, GENO_ID_COL, *COVARIATE_COLS}
missing = required.difference(cov_all.columns)

if missing:
    raise ValueError(
        "Missing covariate columns: "
        + ", ".join(sorted(missing))
    )

cov_all[ID_COL] = cov_all[ID_COL].astype(str)
cov_all[GENO_ID_COL] = cov_all[GENO_ID_COL].astype(str)

if cov_all[ID_COL].duplicated().any():
    dup_ids = (
        cov_all.loc[cov_all[ID_COL].duplicated(keep=False), ID_COL]
        .drop_duplicates()
        .head()
        .tolist()
    )
    raise ValueError(
        f"Duplicated {ID_COL} values detected. Examples: {dup_ids}"
    )

if cov_all[GENO_ID_COL].duplicated().any():
    dup_ids = (
        cov_all.loc[cov_all[GENO_ID_COL].duplicated(keep=False), GENO_ID_COL]
        .drop_duplicates()
        .head()
        .tolist()
    )
    raise ValueError(
        f"Duplicated {GENO_ID_COL} values detected. Examples: {dup_ids}"
    )

# The row order of this file is treated as the genotype sample order.
cov_all = cov_all.set_index(ID_COL)

expr_samples = expr.columns.tolist()
cov_ids = set(cov_all.index)
expr_sample_set = set(expr_samples)

# Preserve COV_FILE row order, which is the genotype sample order.
# This is intentionally NOT the original expression-column order.
common_samples = [
    s for s in cov_all.index
    if s in expr_sample_set
]

expr_only = [
    s for s in expr_samples
    if s not in cov_ids
]

cov_only = [
    s for s in cov_all.index
    if s not in expr_sample_set
]

print(f"Covariate samples   : {cov_all.shape[0]:,}")
print(f"Matched samples     : {len(common_samples):,}")
print(f"Expression only     : {len(expr_only):,}")
print(f"Covariate only      : {len(cov_only):,}")

pd.Series(
    expr_only,
    name=ID_COL,
).to_csv(
    OUT_DIR / f"{cell_type}.expression_only_samples.csv",
    index=False,
)

pd.Series(
    cov_only,
    name=ID_COL,
).to_csv(
    OUT_DIR / f"{cell_type}.covariate_only_samples.csv",
    index=False,
)

if len(common_samples) == 0:
    raise RuntimeError(
        "No overlapping sample IDs between expression columns and "
        f"covariate {ID_COL}."
    )

# Expression-only samples are expected here because the expression object
# contains more subjects than the N103 genotype set. They are excluded rather
# than treated as an error.
if expr_only:
    print(
        f"NOTE: excluding {len(expr_only)} expression samples without genotype/covariate matches."
    )


# ==============================================================================
# 4. Matched phenotype/expression and PEER inputs
# ==============================================================================

# genes x matched samples, ordered exactly as COV_FILE / genotype order.
expr_common = expr.loc[:, common_samples].copy()

genotype_ids = (
    cov_all.loc[common_samples, GENO_ID_COL]
    .astype(str)
    .tolist()
)

if len(set(genotype_ids)) != len(genotype_ids):
    raise RuntimeError("Matched genotype_id values are not unique.")

# Write two audit/phenotype matrices:
#   1) original expression_id column names, but in genotype order
#   2) genotype_id column names, in genotype order -- preferred for QTL phenotype use
MATCHED_EXPR_EXPRID_FILE = (
    CELL_DIR
    / f"{PREFIX}.N{len(common_samples)}.matched.expression_id.geno_order.csv"
)

MATCHED_EXPR_GENOID_FILE = (
    CELL_DIR
    / f"{PREFIX}.N{len(common_samples)}.matched.genotype_id.geno_order.csv"
)

expr_common.to_csv(MATCHED_EXPR_EXPRID_FILE)

expr_common_genoid = expr_common.copy()
expr_common_genoid.columns = genotype_ids
expr_common_genoid.to_csv(MATCHED_EXPR_GENOID_FILE)

sample_map = pd.DataFrame({
    ID_COL: common_samples,
    GENO_ID_COL: genotype_ids,
})
sample_map.to_csv(
    OUT_DIR / f"{cell_type}.matched_expression_genotype_sample_map.csv",
    index=False,
)

print()
print("Matched phenotype files:")
print(f"  expression_id labels : {MATCHED_EXPR_EXPRID_FILE}")
print(f"  genotype_id labels   : {MATCHED_EXPR_GENOID_FILE}")

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

# PEER expects samples x genes.
expr_NxG = expr_common.T

if expr_NxG.index.tolist() != cov_common.index.tolist():
    raise RuntimeError(
        "Expression/covariate sample order mismatch after matching."
    )

# Verify the matched order also corresponds exactly to genotype_id order.
if genotype_ids != cov_all.loc[common_samples, GENO_ID_COL].astype(str).tolist():
    raise RuntimeError("Genotype sample order mismatch after matching.")

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
    / f"{PREFIX}.peertool.expression.csv"
)

PEER_COV_FILE = (
    OUT_DIR
    / f"{cell_type}.age_sex_tissue_geno5PCs.peertool.covariates.csv"
)

# Numeric-only inputs for peertool.
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

# Save labeled copies for auditing/reproducibility.
expr_NxG.to_csv(
    OUT_DIR
    / f"{PREFIX}.peertool.expression.labeled.csv"
)

# Keep genotype_id and any other metadata in this audit copy.
cov_all.loc[
    common_samples
].reset_index().to_csv(
    OUT_DIR
    / f"{cell_type}.age_sex_tissue_geno5PCs.matched.labeled.csv",
    index=False,
)

sample_map.to_csv(
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

# Residual expression: genes x samples.
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
        / f"{PREFIX}.PEER_residuals.expression_id.csv.gz",
        compression="gzip",
    )

    # QTL-ready version: same residuals, columns renamed to genotype_id while
    # retaining the exact genotype sample order.
    resid_genoid = resid.copy()
    resid_genoid.columns = genotype_ids
    resid_genoid.to_csv(
        OUT_DIR
        / f"{PREFIX}.PEER_residuals.genotype_id.csv.gz",
        compression="gzip",
    )

    print(
        f"Residuals: {resid.shape} -> genes x samples"
    )
    print(
        "QTL-ready residuals: "
        f"{OUT_DIR / f'{PREFIX}.PEER_residuals.genotype_id.csv.gz'}"
    )
else:
    print(f"WARNING: PEER residual output not found: {residuals_file}")

# X.
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
else:
    print(f"WARNING: PEER X output not found: {x_file}")

# W.
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
else:
    print(f"WARNING: PEER W output not found: {w_file}")


# ==============================================================================
# 8. Summary
# ==============================================================================

SUMMARY_FILE = (
    OUT_DIR
    / f"{cell_type}.PEER_summary.txt"
)

with open(SUMMARY_FILE, "w") as fh:
    fh.write("project=HTAN_GSE266546\n")
    fh.write(f"cell_type={cell_type}\n")
    fh.write(f"expression_file={EXPR_FILE}\n")
    fh.write(f"covariate_file={COV_FILE}\n")
    fh.write(f"expression_id_column={ID_COL}\n")
    fh.write(f"genotype_id_column={GENO_ID_COL}\n")
    fh.write("matched_sample_order=covariate_file_genotype_order\n")
    fh.write(f"matched_expression_expression_id_file={MATCHED_EXPR_EXPRID_FILE}\n")
    fh.write(f"matched_expression_genotype_id_file={MATCHED_EXPR_GENOID_FILE}\n")
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
        "input_file_suffix=tpmlog2.qc.qn.rmlowcounts.inv\n"
    )
    fh.write(
        "additional_QN_before_PEER=False\n"
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
print(f"Matched phenotype (genotype IDs): {MATCHED_EXPR_GENOID_FILE}")
print(f"Summary   : {SUMMARY_FILE}")
