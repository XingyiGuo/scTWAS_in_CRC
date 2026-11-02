#!/usr/bin/env python3
"""
30 - Figure 18: per-cell CosMx expression maps for twelve further risk genes.

Same construction as Figure 16b-e, for twelve further genes chosen on two criteria that
are both checkable rather than editorial:

  statistical  every gene is significant at FDR < 0.05 among the tested risk genes, in
               the TWAS-predicted direction (§2.14). Genes that fail that test are not
               plotted - a spatial map of a gene with no differential expression invites
               the reader to see a pattern the statistics do not support.
  biological   every gene is in direct protein-protein contact with a colorectal-cancer
               somatic driver (§2.4), which is the concrete sense of "CRC-relevant" this
               report can defend.

Rows group them by class and by where the evidence is strongest. Each panel carries two
statistics from two different datasets, and both are named on the panel: above, the CosMx
paired carcinoma-versus-normal-epithelium difference with its nominal p (Wilcoxon signed
rank across the patients who have both populations); below, the GeoMx stage FDR that
qualified the gene for the figure.
"""
import pandas as pd, numpy as np, json, os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.colors import LinearSegmentedColormap

os.makedirs("figures", exist_ok=True)
INK, INK2, INK3 = "#0b0b0b", "#52514e", "#8a8880"
SURF = "#fcfcfb"
EXPCM = LinearSegmentedColormap.from_list("exp", ["#e9e7e1", "#ffd9a0", "#f59b42",
                                                  "#e34948", "#8f1d2e"])
plt.rcParams.update({
    "figure.facecolor": SURF, "axes.facecolor": SURF, "savefig.facecolor": SURF,
    "font.family": "DejaVu Sans", "font.size": 8.5,
    "axes.edgecolor": "#d8d7d2", "axes.linewidth": 0.8,
    "axes.labelcolor": INK2, "text.color": INK,
    "xtick.color": INK2, "ytick.color": INK2, "xtick.labelsize": 8, "ytick.labelsize": 8,
    "axes.titlesize": 9.5, "axes.titleweight": "bold", "axes.titlelocation": "left",
    "axes.spines.top": False, "axes.spines.right": False,
    "grid.color": "#eceae5", "grid.linewidth": 0.7, "legend.frameon": False,
    "legend.fontsize": 7.8, "pdf.fonttype": 42, "ps.fonttype": 42,
})

z = np.load("spatial/cx_plot_cells3.npz", allow_pickle=True)
X = z["x"].astype(float)
Y = z["y"].astype(float)
Y = Y.max() - Y
E = z["expr_q"].astype(np.float32) / float(z["qscale"])
genes = [str(g) for g in z["genes"]]

B = pd.read_csv("results/cosmx_crc_vs_normal_epi.csv").set_index("gene")
V = pd.read_csv("results/validated_risk_genes_all.csv").set_index("gene")

ROWS = [("Oncogene-like, strongest stage evidence", ["SMARCD1", "CSNK2B", "PTGES3", "GTF2I"]),
        ("Oncogene-like, next strongest", ["UBE2V1", "UBE2B", "PSMD3", "LAMC1"]),
        ("Tumour-suppressor-like, strongest stage evidence", ["SMAD9", "BMP2", "WNT4", "HLA-G"])]

fig = plt.figure(figsize=(12.4, 11.0))
gs = fig.add_gridspec(6, 4, height_ratios=[0.16, 1, 0.16, 1, 0.16, 1],
                      hspace=0.22, wspace=0.20)
letters = "abcdefghijkl"
k = 0
for r, (rowlab, gs_row) in enumerate(ROWS):
    for c, g in enumerate(gs_row):
        ax = fig.add_subplot(gs[2 * r + 1, c])
        j = genes.index(g)
        e = E[:, j]
        neg = e <= 0
        ax.scatter(X[neg], Y[neg], s=0.10, c="#eceae5", linewidths=0, rasterized=True)
        pos = ~neg
        o = np.argsort(e[pos])
        vmax = float(np.quantile(e[pos], 0.98)) if pos.sum() > 50 else 1.0
        sc = ax.scatter(X[pos][o], Y[pos][o], s=0.34, c=e[pos][o], cmap=EXPCM,
                        vmin=0, vmax=max(vmax, 0.2), linewidths=0, rasterized=True)
        ax.set_aspect("equal"); ax.axis("off")

        d = B.loc[g] if g in B.index else None
        v = V.loc[g] if g in V.index else None
        # Two statistics from two different datasets sit on every panel, so each is
        # named: the CosMx paired single-cell test above, the GeoMx stage FDR below.
        cos_line = (f"CosMx Δ {d['delta']:+.3f}, P = {d['p']:.3f} (n = {int(d['n_patients'])})"
                    ) if d is not None else ""
        assert v is not None and bool(v["validated"]), f"{g} is not FDR-significant"
        badge = (f"GeoMx FDR {v['lead_FDR']:.1g} · "
                 f"{v['class'].replace('tumour-suppressor-like', 'TSG-like')}")
        bcol = "#1baf7a"
        ax.text(0.0, 1.105, g, transform=ax.transAxes, fontsize=9.6,
                fontweight="bold", color=INK, va="baseline", ha="left")
        ax.text(0.0, 1.035, cos_line, transform=ax.transAxes, fontsize=7.2,
                color=INK2, va="baseline", ha="left")
        ax.text(0.0, -0.015, badge, transform=ax.transAxes, fontsize=7,
                color=bcol, va="top", ha="left")
        ax.text(-0.075, 1.19, letters[k], transform=ax.transAxes, fontsize=11,
                fontweight="bold", color=INK, va="top", ha="left")
        cb = fig.colorbar(sc, ax=ax, fraction=0.030, pad=0.005, shrink=0.62,
                          location="left", anchor=(0.0, 0.06))
        cb.ax.tick_params(labelsize=6.0, pad=1); cb.outline.set_visible(False)
        cb.set_label("log$_2$ norm.", fontsize=6.0, labelpad=2)
        k += 1
    hax = fig.add_subplot(gs[2 * r, :])
    hax.axis("off")
    hax.text(0.0, 0.85, rowlab, transform=hax.transAxes, va="top", ha="left",
             fontsize=9.4, color=INK, fontweight="bold")

fig.text(0.5, 0.975,
         f"Slide {int(z['slide'])}, {int(z['n_cells_slide']):,} cells; "
         f"{len(X):,} plotted (every cell expressing a listed gene, plus context). "
         "Two statistics per panel, from two datasets.",
         ha="center", fontsize=8.2, color=INK2)
fig.text(0.5, 0.955,
         "Above each map: the CosMx paired carcinoma-minus-normal-epithelium difference "
         "across the n patients carrying both populations, with its nominal p. "
         "Below: the GeoMx stage FDR that qualified the gene for this figure.",
         ha="center", fontsize=8.2, color=INK2)
fig.savefig("figures/fig18_cosmx_more_genes.png", dpi=220, bbox_inches="tight")
fig.savefig("figures/fig18_cosmx_more_genes.pdf", bbox_inches="tight")
plt.close(fig)
print("fig18 done")
