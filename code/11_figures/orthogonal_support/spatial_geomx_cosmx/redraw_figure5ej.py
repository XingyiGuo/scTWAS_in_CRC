#!/usr/bin/env python3
"""Redraw the six CosMx gene panels (Figure 5e–j).

Pass --workbook SourceData.xlsx, or run it in the folder holding the three
archived source-data CSVs:

    python redraw_figure5ej.py

Writes Figure5ej_CosMx_selected.png / .pdf. Drawing conventions follow the
original code (AI-gene-drug code/48_screened_and_selected.py, lines 320-415).
"""
import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--workbook", type=Path)
parser.add_argument("--out-dir", type=Path, default=Path.cwd())
args = parser.parse_args()

if args.workbook:
    workbook = pd.ExcelFile(args.workbook)
    G = workbook.parse("Fig5ej CosMx Gene Summary")
    PAT = workbook.parse("Fig5ej CosMx Paired Means")
else:
    G = pd.read_csv("Fig30_cosmx_selected_sourcedata.csv")
    PAT = pd.read_csv("Fig30_cosmx_selected_sourcedata_patients.csv")
    CELL = pd.read_csv("Fig30_cosmx_selected_sourcedata_percell.csv.gz")

INK, INK2 = "#0b0b0b", "#52514e"
SURF = "#fcfcfb"
EXPCM = LinearSegmentedColormap.from_list(
    "exp", ["#e9e7e1", "#ffd9a0", "#f59b42", "#e34948", "#8f1d2e"])
NORMC, CRCC, GREEN, GREY = "#2a78d6", "#c3202f", "#1baf7a", "#c9c7c1"
# style block copied verbatim from code/48_screened_and_selected.py, lines 56-66
plt.rcParams.update({
    "figure.facecolor": SURF, "axes.facecolor": SURF, "savefig.facecolor": SURF,
    "font.family": "DejaVu Sans", "font.size": 8.5,
    "axes.edgecolor": "#d8d7d2", "axes.linewidth": 0.8,
    "axes.labelcolor": INK2, "text.color": INK,
    "xtick.color": INK2, "ytick.color": INK2, "xtick.labelsize": 8, "ytick.labelsize": 8,
    "axes.titlesize": 9.5, "axes.titleweight": "bold", "axes.titlelocation": "left",
    "axes.spines.top": False, "axes.spines.right": False,
    "grid.color": "#eceae5", "grid.linewidth": 0.7, "legend.frameon": False,
    "legend.fontsize": 7.6, "pdf.fonttype": 42, "ps.fonttype": 42,
})

ncol, nrow = 3, 2
fig = plt.figure(figsize=(3.45 * ncol, 4.35 * nrow))
gs = fig.add_gridspec(nrow * 3, ncol, height_ratios=[1.0, 0.62, 0.34] * nrow,
                      hspace=0.42, wspace=0.34)
heads = []
for k, r in G.iterrows():
    rr, cc = divmod(k, ncol)
    g = r.gene

    # ---- map: non-expressing cells grey underneath, expressing cells on top in
    # ascending order of expression, colour capped at the per-gene ceiling
    if args.workbook:
        c = workbook.parse(f"Fig5{'efghij'[k]} {g} Cells",
                           usecols=["x_mm", "y_mm_plotted", "expression_log2_norm", "expressing"])
    else:
        c = CELL[CELL.gene == g]
    ax = fig.add_subplot(gs[3 * rr, cc])
    neg = c[~c.expressing]
    ax.scatter(neg.x_mm, neg.y_mm_plotted, s=0.10, c="#eceae5", linewidths=0,
               rasterized=True)
    pos = c[c.expressing].sort_values("expression_log2_norm", kind="stable")
    sc = ax.scatter(pos.x_mm, pos.y_mm_plotted, s=0.30, c=pos.expression_log2_norm,
                    cmap=EXPCM, vmin=0, vmax=r.map_colour_max_log2_norm,
                    linewidths=0, rasterized=True)
    ax.set_aspect("equal"); ax.axis("off")
    cb = fig.colorbar(sc, ax=ax, fraction=0.030, pad=0.005, shrink=0.60,
                      location="left", anchor=(0.0, 0.06))
    cb.ax.tick_params(labelsize=6.0, pad=1); cb.outline.set_visible(False)
    cb.set_label("log$_2$ norm.", fontsize=6.0, labelpad=2)
    heads.append((rr, cc, ax, r))

    # ---- paired dot plot: one line per patient, red if it rises in carcinoma
    p = PAT[PAT.gene == g]
    ax = fig.add_subplot(gs[3 * rr + 1, cc])
    for _, q in p.iterrows():
        ax.plot([0, 1], [q.mean_norm_expression_normal_epithelium,
                         q.mean_norm_expression_carcinoma],
                color="#e6b8b8" if q.rises_in_carcinoma else GREY, lw=1.0, zorder=2)
    ax.scatter(np.zeros(len(p)), p.mean_norm_expression_normal_epithelium, s=28,
               color=NORMC, edgecolor="white", linewidth=0.6, zorder=4)
    ax.scatter(np.ones(len(p)), p.mean_norm_expression_carcinoma, s=28,
               color=CRCC, edgecolor="white", linewidth=0.6, zorder=4)
    for xx, v, col in [(0, p.mean_norm_expression_normal_epithelium, NORMC),
                       (1, p.mean_norm_expression_carcinoma, CRCC)]:
        ax.plot([xx - 0.16, xx + 0.16], [v.mean()] * 2, color=col, lw=2.2, zorder=5)
    ax.set_xticks([0, 1])
    ax.set_xticklabels(["normal\nepithelium", "carcinoma\ncells"], fontsize=7.4)
    ax.set_xlim(-0.35, 1.35)
    if cc == 0:
        ax.set_ylabel("Mean normalised\nexpression per patient", fontsize=7.8)
    pw = r.dotplot_paired_wilcoxon_p
    ax.set_title(f"paired Wilcoxon P = {pw:.3f}", fontsize=8.4, pad=15,
                 color=GREEN if pw < 0.05 else INK2, loc="center")
    ax.text(0.5, 1.012, f"{r.dotplot_n_rise} of {r.dotplot_n_patients} patients rise",
            transform=ax.transAxes, ha="center", va="bottom", fontsize=7.4, color=INK2)
    ax.grid(axis="y"); ax.set_axisbelow(True)

fig.canvas.draw()
for rr, cc, ax, r in heads:
    cell = gs[3 * rr, cc].get_position(fig); m = ax.get_position()
    fig.text(cell.x0 + 0.019, m.y1 + 0.040, r.gene, fontsize=10.2, fontweight="bold",
             color=INK)
    fig.text(cell.x0 + 0.019, m.y1 + 0.024,
             f"CosMx \u0394 {r.map_heading_cosmx_delta:+.3f}, "
             f"P = {r.map_heading_cosmx_p:.3f} (n = {r.dotplot_n_patients})",
             fontsize=7.2, color=INK2)
    fig.text(cell.x0 + 0.019, m.y1 + 0.009, str(r.map_heading_selection_reason),
             fontsize=7, color="#4a3aa7")
    fig.text(cell.x0 - 0.002, m.y1 + 0.055, "efghij"[rr * ncol + cc], fontsize=11,
             fontweight="bold", color=INK, va="top")

if args.workbook:
    workbook.close()
args.out_dir.mkdir(parents=True, exist_ok=True)
fig.savefig(args.out_dir / "Figure5ej_CosMx_selected.png", dpi=220, bbox_inches="tight")
fig.savefig(args.out_dir / "Figure5ej_CosMx_selected.pdf", bbox_inches="tight")
