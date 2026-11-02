#!/usr/bin/env python3
"""
25 - Per-cell CosMx expression of the risk genes (runs on the local machine).

Reads the counts pulled out of the Seurat object by 24_rds_reader.py, joins them to the
per-cell metadata, and reduces 1.26 M cells to a small pseudobulk table that can be
carried back for statistics and plotting:

    patient x tumour stage x cell type x gene  ->  n_cells, n_positive, mean expression

Working at the level of patient means is deliberate. Cells from one patient are not
independent, so any test that treats 1.26 M cells as 1.26 M observations will return
absurd p-values; every test downstream is therefore computed across patients on these
means. The per-cell values are kept only for the figures, and only for the genes that
are actually plotted.

Normalisation is library size from the metadata's nCount_RNA (total transcripts detected
in that cell, over the whole 6,175-gene panel), expressed per 1,000 transcripts and
log2-transformed:  log2(1 + 1000 * count / nCount_RNA).

Outputs (all small enough to carry back):
    cx_pseudobulk.csv.gz        patient x stage x cell type x gene
    cx_celltype_cells.csv       cells per patient x stage x cell type
    cx_panel_genes.txt          which risk genes are on the CosMx panel
    cx_plot_cells.npz           per-cell values for the plotted genes, one slide
"""
import numpy as np, pandas as pd, gzip, json, sys, os

CX = os.path.expanduser("~/cx")
VAL = os.path.expanduser("~/mnt/AI-gene-drug/results/validated_risk_genes.csv")
N_PLOT_GENES = 8

# ------------------------------------------------------------------ load
z = np.load(f"{CX}/cx_counts.npz", allow_pickle=True)
row, col, val = z["row"].astype(np.int32), z["col"].astype(np.int64), z["val"].astype(np.float32)
gene_of_row = [str(g) for g in z["gene_of_row"]]
rows_sorted = np.array(z["rows"], dtype=np.int32)
ncell = int(z["ncell"])
gene_code = np.searchsorted(rows_sorted, row).astype(np.int32)   # 0..n_gene-1
n_gene = len(gene_of_row)
print(f"{len(row):,} entries, {n_gene} genes, {ncell:,} cells")

with gzip.open(f"{CX}/cx_cells.txt.gz", "rt") as fh:
    barcodes = fh.read().split("\n")
assert len(barcodes) == ncell, (len(barcodes), ncell)

meta = pd.read_csv(f"{CX}/meta_sel.tsv.gz", sep="\t", low_memory=False,
                   usecols=["fov", "slide_ID_numeric", "x_slide_mm", "y_slide_mm", "cell_ID",
                            "nCount_RNA", "patient_ID", "tumorStage", "cellType",
                            "cellTypeEpithelium", "neighborhood"])
# cell_ID is already the Seurat barcode (c_<slide>_<fov>_<cell>)
meta = meta.set_index("cell_ID").reindex(barcodes)
n_missing = int(meta["patient_ID"].isna().sum())
print(f"metadata joined to matrix columns; unmatched cells: {n_missing:,}")
if n_missing > 0.01 * ncell:
    raise SystemExit("barcode join failed - more than 1 % of columns unmatched")

# ------------------------------------------------------------------ codes
meta["cellType"] = meta["cellType"].fillna("unannotated")
pat = pd.Categorical(meta["patient_ID"].fillna("NA"))
ct = pd.Categorical(meta["cellType"])
stage = pd.Categorical(meta["tumorStage"].fillna("NA"), categories=["im", "T1.sm1", "T1.sm3"])
pat_code = pat.codes.astype(np.int32)
ct_code = ct.codes.astype(np.int32)
lib = meta["nCount_RNA"].to_numpy(dtype=np.float64)
lib[~np.isfinite(lib) | (lib <= 0)] = np.nan

n_pat, n_ct = len(pat.categories), len(ct.categories)
print(f"{n_pat} patients, {n_ct} cell types")

# normalised value of every stored entry
norm = np.log2(1.0 + 1000.0 * val / lib[col]).astype(np.float32)
good = np.isfinite(norm)
norm = np.where(good, norm, 0.0)

# ------------------------------------------------------------------ pseudobulk
cell_group = pat_code.astype(np.int64) * n_ct + ct_code          # per cell
entry_group = cell_group[col] * n_gene + gene_code               # per stored entry
G = n_pat * n_ct * n_gene

sum_norm = np.bincount(entry_group, weights=norm.astype(np.float64), minlength=G)
n_pos = np.bincount(entry_group, minlength=G)
sum_raw = np.bincount(entry_group, weights=val.astype(np.float64), minlength=G)
cells_per_group = np.bincount(cell_group[np.isfinite(lib)], minlength=n_pat * n_ct)

gi, ci, gg = np.meshgrid(np.arange(n_pat), np.arange(n_ct), np.arange(n_gene), indexing="ij")
flat = pd.DataFrame({
    "patient_ID": np.asarray(pat.categories)[gi.ravel()],
    "cellType": np.asarray(ct.categories)[ci.ravel()],
    "gene": np.asarray(gene_of_row)[gg.ravel()],
    "n_cells": np.repeat(cells_per_group, n_gene),
    "n_positive": n_pos,
    "sum_norm": sum_norm,
    "sum_counts": sum_raw,
})
flat = flat[flat["n_cells"] >= 1].copy()
flat["mean_norm"] = flat["sum_norm"] / flat["n_cells"]
flat["frac_positive"] = flat["n_positive"] / flat["n_cells"]

stage_of_patient = (meta.groupby("patient_ID", observed=True)["tumorStage"]
                    .agg(lambda s: s.dropna().iloc[0] if s.notna().any() else "NA"))
flat["tumorStage"] = flat["patient_ID"].map(stage_of_patient)
flat = flat[["patient_ID", "tumorStage", "cellType", "gene", "n_cells", "n_positive",
             "frac_positive", "mean_norm", "sum_counts"]]
flat.to_csv(f"{CX}/cx_pseudobulk.csv.gz", index=False, float_format="%.6g")
print(f"pseudobulk rows: {len(flat):,}  ->  cx_pseudobulk.csv.gz")

(meta.groupby(["patient_ID", "tumorStage", "cellType"], observed=True)
     .size().rename("n_cells").reset_index()
     .to_csv(f"{CX}/cx_celltype_cells.csv", index=False))
open(f"{CX}/cx_panel_genes.txt", "w").write("\n".join(gene_of_row))

# ------------------------------------------------------------------ cells for the figures
val_df = pd.read_csv(VAL)
on_panel = val_df[val_df["gene"].isin(gene_of_row)].sort_values("lead_FDR")
plot_genes = list(on_panel["gene"].head(N_PLOT_GENES))
print(f"validated genes on the CosMx panel: {len(on_panel)}  plotting: {plot_genes}")

slide_counts = meta["slide_ID_numeric"].value_counts()
slide = int(slide_counts.index[0])
sel_slide = (meta["slide_ID_numeric"].to_numpy() == slide)
gcodes = {g: gene_of_row.index(g) for g in plot_genes}

# dense per-cell values for the plotted genes on this slide
idx = np.flatnonzero(sel_slide)
pos_in_slide = -np.ones(ncell, dtype=np.int64)
pos_in_slide[idx] = np.arange(len(idx))
E = np.zeros((len(idx), len(plot_genes)), dtype=np.float32)
for j, g in enumerate(plot_genes):
    m = (gene_code == gcodes[g]) & sel_slide[col]
    E[pos_in_slide[col[m]], j] = norm[m]

# keep every expressing cell, plus a sample of the rest for tissue context
rng = np.random.default_rng(25)
expressing = (E > 0).any(1)
context = (~expressing) & (rng.random(len(idx)) < 0.22)
keep = expressing | context
print(f"slide {slide}: {len(idx):,} cells, plotting {int(keep.sum()):,} "
      f"({int(expressing.sum()):,} expressing)")

sub = meta.iloc[idx[keep]]
np.savez_compressed(
    f"{CX}/cx_plot_cells.npz",
    x=sub["x_slide_mm"].to_numpy(np.float32), y=sub["y_slide_mm"].to_numpy(np.float32),
    expr=E[keep].astype(np.float16),
    genes=np.array(plot_genes),
    cellType=np.array(sub["cellType"].astype(str)),
    stage=np.array(sub["tumorStage"].astype(str)),
    patient=np.array(sub["patient_ID"].astype(str)),
    slide=np.int32(slide),
    n_cells_slide=np.int64(len(idx)),
    frac_kept=np.float32(keep.mean()),
)
print(f"wrote cx_plot_cells.npz "
      f"({os.path.getsize(f'{CX}/cx_plot_cells.npz')/1e6:.1f} MB)")
print("Done.")
