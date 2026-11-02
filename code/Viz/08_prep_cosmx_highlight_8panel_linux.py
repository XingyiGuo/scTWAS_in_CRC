#!/usr/bin/env python3
"""Prepare CosMx data for an 8-panel Linux figure.

Figure order (2 rows x 4 columns)
  a  Cell composition
  b  CSNK2B
  c  DDX39B
  d  RPS5
  e  LIMA1
  f  CLIC1
  g  CTNNB1
  h  UBE2V1

The script searches the CosMx experiment folders recursively for NumPy .npz
archives with the same schema used by the original plotting pipeline:
    x, y, genes, cellType, and either expr OR expr_q + qscale.

For each gene, the archive containing the *largest number of cells* is used.
For panel a, the largest valid archive containing cellType is used.  Archives
are never merged across coordinate systems because that would superimpose
unrelated slides/field-of-view coordinates.  No cells are subsampled.

Default source path (override with COSMX_ROOT):
  /data/l2_bioinfo1/liq17/scTWAS/data/SpatialExpr/CosMx_PMID41986711
    Nanostring_CosMx_exp1/
    Nanostring_CosMx_exp2/
    Nanostring_GeoMx/        # not needed for the per-cell CosMx maps

Outputs (to COSMX_FIG_WORKDIR, default ./cosmx_8panel_work):
  cosmx_panelA_composition.csv.gz
  cosmx_expr_<GENE>.csv.gz
  cosmx_panel_manifest.csv
  cosmx_archive_inventory.csv

Example Linux run:
  python3 08_prep_cosmx_highlight_8panel_linux.py

Optional environment overrides:
  export COSMX_ROOT=/data/l2_bioinfo1/liq17/scTWAS/data/SpatialExpr/CosMx_PMID41986711
  export COSMX_ARCHIVE_ROOT=$COSMX_ROOT
  export COSMX_FIG_WORKDIR=$PWD/cosmx_8panel_work
"""

from __future__ import annotations

import os
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np
import pandas as pd

DATA_ROOT = Path(
    os.environ.get(
        "COSMX_ROOT",
        "/data/l2_bioinfo1/liq17/scTWAS/data/SpatialExpr/CosMx_PMID41986711",
    )
).expanduser().resolve()

# If the *.npz files made by the extraction pipeline live somewhere else, set
# COSMX_ARCHIVE_ROOT to that directory.  By default search the data root.
ARCHIVE_ROOT = Path(os.environ.get("COSMX_ARCHIVE_ROOT", str(DATA_ROOT))).expanduser().resolve()
WORKDIR = Path(os.environ.get("COSMX_FIG_WORKDIR", "cosmx_8panel_work")).expanduser().resolve()
WORKDIR.mkdir(parents=True, exist_ok=True)

GENES = ["CSNK2B", "DDX39B", "RPS5", "LIMA1", "CLIC1", "CTNNB1", "UBE2V1"]
COSMX_DIRS = ["Nanostring_CosMx_exp1", "Nanostring_CosMx_exp2"]

GROUP = {
    "CRC": "Carcinoma",
    "normal/undefined epithelium": "Normal epithelium",
    "oncofetal": "Oncofetal",
}


def group_celltype(c: str) -> str:
    """Collapse CosMx cell types to the six groups used in the original figure."""
    c = str(c)
    if c in GROUP:
        return GROUP[c]
    cl = c.lower()
    if (
        "fibro" in cl
        or "caf" in cl
        or "trophocyte" in cl
        or "telocyte" in cl
        or c == "stromal.1"
    ):
        return "Fibroblast / stroma"
    if (
        c.startswith("T.")
        or c in ("B cell", "plasma", "macrophage", "mast cells", "MKI67+ lymphocytes")
        or "lymph" in cl
    ):
        return "Immune"
    return "Other"


def candidate_npz_files() -> List[Path]:
    """Find archives, prioritising the two CosMx folders but remaining robust."""
    roots: List[Path] = []
    for d in COSMX_DIRS:
        p = ARCHIVE_ROOT / d
        if p.exists():
            roots.append(p)

    # Extraction outputs sometimes sit directly under ARCHIVE_ROOT rather than
    # inside the public-data folders.  Search ARCHIVE_ROOT as a fallback.
    if not roots:
        roots = [ARCHIVE_ROOT]

    files: List[Path] = []
    seen = set()
    for root in roots:
        for p in root.rglob("*.npz"):
            rp = p.resolve()
            if rp not in seen:
                files.append(rp)
                seen.add(rp)

    # Also include direct-root npz files if folder-specific searching was used.
    if ARCHIVE_ROOT.exists():
        for p in ARCHIVE_ROOT.glob("*.npz"):
            rp = p.resolve()
            if rp not in seen:
                files.append(rp)
                seen.add(rp)

    return sorted(files)


def inspect_archive(path: Path) -> dict | None:
    """Read only metadata needed to decide whether an archive is usable."""
    try:
        with np.load(path, allow_pickle=True) as z:
            keys = set(z.files)
            required = {"x", "y", "genes"}
            if not required.issubset(keys):
                return None
            if not ({"expr", "expr_q"} & keys):
                return None

            genes = [str(g) for g in z["genes"]]
            n_cells = int(z["x"].shape[0])
            n_y = int(z["y"].shape[0])
            if n_cells != n_y:
                return None

            has_celltype = "cellType" in keys and int(z["cellType"].shape[0]) == n_cells
            expr_key = "expr_q" if "expr_q" in keys else "expr"
            expr_shape = tuple(z[expr_key].shape)
            if len(expr_shape) != 2 or expr_shape[0] != n_cells or expr_shape[1] != len(genes):
                return None

            return {
                "path": str(path),
                "n_cells": n_cells,
                "n_genes": len(genes),
                "genes": genes,
                "has_celltype": has_celltype,
                "expr_key": expr_key,
                "qscale": float(z["qscale"]) if "qscale" in keys else np.nan,
            }
    except Exception as e:
        print(f"WARNING: skipping unreadable archive {path}: {e}")
        return None


def load_coordinates_and_celltypes(path: Path) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    with np.load(path, allow_pickle=True) as z:
        x = z["x"].astype(np.float32)
        y = z["y"].astype(np.float32)
        y = y.max() - y  # same orientation as the original figure
        ctype = np.asarray([str(c) for c in z["cellType"]], dtype=object)
    return x, y, ctype


def load_gene(path: Path, gene: str) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Load one gene's values from one archive without altering cell membership."""
    with np.load(path, allow_pickle=True) as z:
        genes = [str(g) for g in z["genes"]]
        if gene not in genes:
            raise KeyError(f"{gene} not in {path}")
        j = genes.index(gene)
        x = z["x"].astype(np.float32)
        y = z["y"].astype(np.float32)
        y = y.max() - y
        if "expr_q" in z:
            qscale = float(z["qscale"])
            v = z["expr_q"][:, j].astype(np.float32) / qscale
        else:
            v = z["expr"][:, j].astype(np.float32)
    return x, y, v


def main() -> None:
    print(f"DATA_ROOT    : {DATA_ROOT}")
    print(f"ARCHIVE_ROOT : {ARCHIVE_ROOT}")
    print(f"WORKDIR      : {WORKDIR}")

    if not ARCHIVE_ROOT.exists():
        raise FileNotFoundError(
            f"Archive root does not exist: {ARCHIVE_ROOT}\n"
            "Set COSMX_ARCHIVE_ROOT to the directory containing the cx_plot_cells*.npz files."
        )

    npz_files = candidate_npz_files()
    if not npz_files:
        raise FileNotFoundError(
            f"No .npz archives found under {ARCHIVE_ROOT}.\n"
            "This plotting workflow expects the per-cell archives produced by the earlier "
            "CosMx extraction step (x, y, genes, cellType, expr/expr_q).\n"
            "Set COSMX_ARCHIVE_ROOT to that extraction-output directory if it is separate "
            "from the public-data root."
        )

    infos = [x for p in npz_files if (x := inspect_archive(p)) is not None]
    if not infos:
        raise RuntimeError("Found .npz files, but none have the expected CosMx plotting schema.")

    inventory_rows = []
    for x in infos:
        inventory_rows.append(
            {
                "path": x["path"],
                "n_cells": x["n_cells"],
                "n_genes": x["n_genes"],
                "has_celltype": x["has_celltype"],
                "expr_key": x["expr_key"],
                "qscale": x["qscale"],
                "highlight_genes_present": ";".join(g for g in GENES if g in x["genes"]),
            }
        )
    pd.DataFrame(inventory_rows).sort_values("n_cells", ascending=False).to_csv(
        WORKDIR / "cosmx_archive_inventory.csv", index=False
    )

    # Panel a: maximum-cell archive that still has cell annotations.
    comp_candidates = [x for x in infos if x["has_celltype"]]
    if not comp_candidates:
        raise RuntimeError("No valid archive contains a cellType vector for panel a.")
    comp_info = max(comp_candidates, key=lambda x: x["n_cells"])
    comp_path = Path(comp_info["path"])
    x, y, ctype = load_coordinates_and_celltypes(comp_path)
    comp = pd.DataFrame(
        {
            "x": np.round(x, 4),
            "y": np.round(y, 4),
            "group": [group_celltype(c) for c in ctype],
        }
    )
    comp.to_csv(WORKDIR / "cosmx_panelA_composition.csv.gz", index=False, compression="gzip")
    print(
        f"panel a: {len(comp):,} cells from {comp_path}; "
        f"groups={sorted(comp['group'].unique())}"
    )

    manifest = [
        {
            "panel": "a",
            "gene": "Cell composition",
            "archive": str(comp_path),
            "n_cells": len(comp),
            "n_expressing": np.nan,
            "expr_vmax": np.nan,
        }
    ]

    panel_letters = dict(zip(GENES, list("bcdefgh")))

    # Each gene uses the archive with the most cells that contains that gene.
    missing = []
    for g in GENES:
        gcand = [x for x in infos if g in x["genes"]]
        if not gcand:
            missing.append(g)
            continue
        best = max(gcand, key=lambda x: x["n_cells"])
        path = Path(best["path"])
        gx, gy, v = load_gene(path, g)
        pos = np.isfinite(v) & (v > 0)
        vmax = max(float(np.quantile(v[pos], 0.98)) if pos.sum() > 50 else 1.0, 0.2)

        out = pd.DataFrame(
            {
                "x": np.round(gx, 4),
                "y": np.round(gy, 4),
                "expr": np.round(v, 4),
                "expressing": pos,
            }
        )
        out_path = WORKDIR / f"cosmx_expr_{g}.csv.gz"
        out.to_csv(out_path, index=False, compression="gzip")

        manifest.append(
            {
                "panel": panel_letters[g],
                "gene": g,
                "archive": str(path),
                "n_cells": len(v),
                "n_expressing": int(pos.sum()),
                "expr_vmax": round(vmax, 4),
            }
        )
        print(
            f"panel {panel_letters[g]} {g}: {len(v):,} cells; "
            f"{int(pos.sum()):,} expressing; q98={vmax:.3f}; source={path}"
        )

    if missing:
        available = sorted({g for x in infos for g in x["genes"]})
        raise KeyError(
            "These requested genes were not found in any valid archive: "
            + ", ".join(missing)
            + "\nInspect cosmx_archive_inventory.csv and/or regenerate the extraction archives "
            "including the requested genes.\n"
            + f"Number of distinct genes available across archives: {len(available):,}"
        )

    mf = pd.DataFrame(manifest)
    mf.to_csv(WORKDIR / "cosmx_panel_manifest.csv", index=False)
    print("\nSelected sources:")
    print(mf.to_string(index=False))
    print(f"\nWrote all plotting inputs to: {WORKDIR}")


if __name__ == "__main__":
    main()
