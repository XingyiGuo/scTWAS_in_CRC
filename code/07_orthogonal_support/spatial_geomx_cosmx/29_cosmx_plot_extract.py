#!/usr/bin/env python3
"""
29 - Per-cell values for a chosen set of genes, ready to plot (runs on the local machine).

Script 25 wrote a plotting file for the eight strongest validated genes. This does the
same for any gene list, reusing the counts already pulled out of the Seurat object by
script 24, so the 2.1 GB object does not have to be read again.

Cells are kept from one slide: every cell expressing at least one of the listed genes,
plus a random sample of the rest for tissue context. Expression is quantised to one byte
(value x 40, capped at 255) because the output has to cross a slow link and a colour map
cannot resolve more than that anyway.

    python3 29_cosmx_plot_extract.py <genelist.txt> <out.npz> [slide]
"""
import numpy as np, pandas as pd, gzip, sys, os

CX = os.path.expanduser("~/cx")
QSCALE = 40.0

genes_wanted = [ln.strip() for ln in open(sys.argv[1]) if ln.strip()]
out_path = sys.argv[2]
slide_req = int(sys.argv[3]) if len(sys.argv) > 3 else None

z = np.load(f"{CX}/cx_counts.npz", allow_pickle=True)
row, col, val = z["row"].astype(np.int32), z["col"].astype(np.int64), z["val"].astype(np.float32)
gene_of_row = [str(g) for g in z["gene_of_row"]]
rows_sorted = np.array(z["rows"], dtype=np.int32)
ncell = int(z["ncell"])
gene_code = np.searchsorted(rows_sorted, row).astype(np.int32)

have = [g for g in genes_wanted if g in gene_of_row]
missing = [g for g in genes_wanted if g not in gene_of_row]
print(f"{len(have)} of {len(genes_wanted)} genes available" +
      (f"; not extracted: {missing}" if missing else ""))

with gzip.open(f"{CX}/cx_cells.txt.gz", "rt") as fh:
    barcodes = fh.read().split("\n")
meta = pd.read_csv(f"{CX}/meta_sel.tsv.gz", sep="\t", low_memory=False,
                   usecols=["slide_ID_numeric", "x_slide_mm", "y_slide_mm", "cell_ID",
                            "nCount_RNA", "patient_ID", "tumorStage", "cellType"])
meta = meta.set_index("cell_ID").reindex(barcodes)
lib = meta["nCount_RNA"].to_numpy(float)
lib[~np.isfinite(lib) | (lib <= 0)] = np.nan
norm = np.log2(1.0 + 1000.0 * val / lib[col]).astype(np.float32)
norm = np.where(np.isfinite(norm), norm, 0.0)

slide = slide_req if slide_req is not None else int(meta["slide_ID_numeric"].value_counts().index[0])
sel = (meta["slide_ID_numeric"].to_numpy() == slide)
idx = np.flatnonzero(sel)
pos_in_slide = -np.ones(ncell, dtype=np.int64)
pos_in_slide[idx] = np.arange(len(idx))

E = np.zeros((len(idx), len(have)), dtype=np.float32)
for j, g in enumerate(have):
    m = (gene_code == gene_of_row.index(g)) & sel[col]
    E[pos_in_slide[col[m]], j] = norm[m]

rng = np.random.default_rng(29)
expressing = (E > 0).any(1)
context = (~expressing) & (rng.random(len(idx)) < 0.30)
keep = np.flatnonzero(expressing | context)
if len(keep) > 260_000:                       # keep the file small enough to move
    keep = rng.choice(keep, 260_000, replace=False)
    keep.sort()
print(f"slide {slide}: {len(idx):,} cells, keeping {len(keep):,} "
      f"({int(expressing.sum()):,} express at least one listed gene)")

sub = meta.iloc[idx[keep]]
Eq = np.clip(np.round(E[keep] * QSCALE), 0, 255).astype(np.uint8)
np.savez_compressed(
    out_path,
    x=sub["x_slide_mm"].to_numpy(np.float32), y=sub["y_slide_mm"].to_numpy(np.float32),
    expr_q=Eq, qscale=np.float32(QSCALE), genes=np.array(have),
    cellType=np.array(sub["cellType"].astype(str)),
    stage=np.array(sub["tumorStage"].astype(str)),
    patient=np.array(sub["patient_ID"].astype(str)),
    slide=np.int32(slide), n_cells_slide=np.int64(len(idx)))
print(f"wrote {out_path} ({os.path.getsize(out_path)/1e6:.1f} MB)")
