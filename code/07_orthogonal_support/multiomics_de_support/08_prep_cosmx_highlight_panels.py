#!/usr/bin/env python3
"""08_prep_cosmx_highlight_panels.py

Decode the CosMx per-cell plotting arrays into tidy CSVs that the R figure
script can read, and assemble the statistics block that annotates each panel.

Why a Python step at all
    The per-cell values live in NumPy .npz archives written by the AI-gene-drug
    project (code/29_cosmx_plot_extract.py). R cannot read .npz without
    reticulate, so this converts once and the figure stays pure R.

Panels being assembled (matching the source figures)
    a  cell composition          fig16 panel a   <- cx_plot_cells.npz
    b  GTF2I                     fig18 panel d   <- cx_plot_cells3.npz
    c  TRIM28                    fig16 panel d   <- cx_plot_cells.npz
    d  UBE2V1                    fig18 panel e   <- cx_plot_cells3.npz

    The two archives hold DIFFERENT cell subsets of the same slide (339,739 vs
    260,000 cells, each being every cell expressing one of that archive's genes
    plus a random context sample). Coordinates are therefore not shared between
    panels and each panel keeps its own cell table. Do not merge them.

Conventions copied from the source figures so the panels match
    - y is flipped within each archive: y <- max(y) - y
    - cx_plot_cells3.npz stores expression quantised to uint8 (value x qscale);
      cx_plot_cells.npz stores float16 directly
    - non-expressing cells are drawn first in grey, expressing cells on top,
      ordered by value, with the colour scale capped at the 98th percentile of
      the expressing values (vmax, written into the stats file)
    - cell-type grouping and colours are the GROUP / GCOL maps from
      code/28_validated_figures.py

Outputs (workspace)
    cosmx_panelA_composition.csv.gz   x, y, group
    cosmx_panel_expr.csv.gz           gene, x, y, expr, expressing
    panel_stats.csv                   one row per gene: scTWAS, CosMx, and
                                      every other evidence layer
"""

from __future__ import annotations

import os

import numpy as np
import pandas as pd

OUT = r"W:\VUMC_projects\scTWAS\VUMC_HTAN_scTWAS_genes_DE_support"
GENES = ["GTF2I", "TRIM28", "UBE2V1"]

## which archive supplies which gene; panel a uses ARCH_MAIN
ARCH_MAIN = "cx_plot_cells.npz"      # float16 expr, 339,739 cells
ARCH_ALT = "cx_plot_cells3.npz"      # uint8 expr_q, 260,000 cells
GENE_ARCHIVE = {"TRIM28": ARCH_MAIN, "GTF2I": ARCH_ALT, "UBE2V1": ARCH_ALT}

GROUP = {"CRC": "Carcinoma",
         "normal/undefined epithelium": "Normal epithelium",
         "oncofetal": "Oncofetal"}


def group_celltype(c: str) -> str:
    """Collapse the 24 CosMx cell types into the six plotted groups."""
    if c in GROUP:
        return GROUP[c]
    if "fibro" in c or "CAF" in c or "trophocyte" in c or "telocyte" in c or c == "stromal.1":
        return "Fibroblast / stroma"
    if c.startswith("T.") or c in ("B cell", "plasma", "macrophage", "mast cells",
                                   "MKI67+ lymphocytes"):
        return "Immune"
    return "Other"


def load_archive(path: str):
    """Return (x, y_flipped, expression matrix as float, gene list, cell types)."""
    z = np.load(path, allow_pickle=True)
    x = z["x"].astype(np.float32)
    y = z["y"].astype(np.float32)
    y = y.max() - y                                  # as in the source figures
    if "expr_q" in z:
        e = z["expr_q"].astype(np.float32) / float(z["qscale"])
    else:
        e = z["expr"].astype(np.float32)
    genes = [str(g) for g in z["genes"]]
    ctype = np.array([str(c) for c in z["cellType"]])
    return x, y, e, genes, ctype


def main() -> None:
    ## ---- panel a: composition ------------------------------------------------
    x, y, _, _, ctype = load_archive(ARCH_MAIN)
    comp = pd.DataFrame({"x": np.round(x, 4), "y": np.round(y, 4),
                         "group": [group_celltype(c) for c in ctype]})
    comp.to_csv("cosmx_panelA_composition.csv.gz", index=False)
    print(f"panel a: {len(comp):,} cells, "
          f"{comp.group.nunique()} groups {sorted(comp.group.unique())}")

    ## ---- gene panels ---------------------------------------------------------
    cache = {}
    frames, vmax = [], {}
    for g in GENES:
        arch = GENE_ARCHIVE[g]
        if arch not in cache:
            cache[arch] = load_archive(arch)
        gx, gy, ge, genes, _ = cache[arch]
        if g not in genes:
            raise KeyError(f"{g} not in {arch}: {genes}")
        v = ge[:, genes.index(g)]
        pos = v > 0
        ## colour cap: 98th percentile of the expressing cells, floored at 0.2,
        ## exactly as the source panels do
        vmax[g] = max(float(np.quantile(v[pos], 0.98)) if pos.sum() > 50 else 1.0, 0.2)
        frames.append(pd.DataFrame({"gene": g, "x": np.round(gx, 4), "y": np.round(gy, 4),
                                    "expr": np.round(v, 3), "expressing": pos}))
        print(f"{g}: {arch}, {len(v):,} cells, {int(pos.sum()):,} expressing, "
              f"vmax {vmax[g]:.2f}")
    expr = pd.concat(frames, ignore_index=True)
    expr.to_csv("cosmx_panel_expr.csv.gz", index=False)

    ## ---- statistics block ----------------------------------------------------
    ms = pd.read_csv(os.path.join(OUT, "table",
                                  "scTWAS_347genes_multisource_DE_support.csv"))
    cxr = pd.read_csv("cosmx_geomx_replication.csv")
    ## (source prefix, display label). The first is the GeoMx compartment DE that
    ## the scTWAS tables mislabel as CosMx.
    layers = [("CosMx_NOR_CORE", "GeoMx_PanCKpos_NOR_CORE"),
              ("Visium_M1", "Visium_M1"), ("Visium_M3", "Visium_M3"),
              ("TCGA", "TCGA_bulk_RNA"), ("CPTAC", "CPTAC_bulk_protein")]
    sm = ms.set_index("gene_name")
    rows = []
    for g in GENES:
        r, cx = sm.loc[g], cxr[cxr.gene == g]
        d = {"gene": g, "panel_source": {"GTF2I": "fig18d", "TRIM28": "fig16d",
                                         "UBE2V1": "fig18e"}[g],
             "rep_model": r.rep_model, "rep_Z": r.rep_Z, "rep_P": r.rep_P,
             "n_models": int(r.n_models), "direction_status": r.direction_status,
             "cosmx_delta": float(cx.delta.iloc[0]), "cosmx_p": float(cx.p.iloc[0]),
             "cosmx_fdr": float(cx.fdr.iloc[0]),
             "cosmx_n_patients": int(cx.n_patients.iloc[0]),
             "expr_vmax": round(vmax[g], 3)}
        for pre, lab in layers:
            d[f"{lab}_log2FC"] = r.get(f"{pre}_log2FC")
            d[f"{lab}_padj"] = r.get(f"{pre}_padj")
            d[f"{lab}_label"] = r.get(f"{pre}_DEsupport")
        ## layers counted: GeoMx + single-cell CosMx + every other supported source
        d["n_layers"] = (1
                         + int(d["cosmx_p"] < 0.05 and d["cosmx_delta"] > 0)
                         + sum(d[f"{lab}_label"] in ("Up", "Down")
                               for _, lab in layers[1:]))
        rows.append(d)
    st = pd.DataFrame(rows)
    st.to_csv("panel_stats.csv", index=False)
    print("\n" + st[["gene", "panel_source", "rep_Z", "cosmx_delta", "cosmx_p",
                     "n_layers", "expr_vmax"]].to_string(index=False))


if __name__ == "__main__":
    main()
