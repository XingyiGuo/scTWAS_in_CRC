#!/usr/bin/env python3
"""
47 - Figure 28: per-cell CosMx expression for a requested gene set.

The eight genes asked for were DDX39B, CSNK2B, CLIC1, CTNNB1, UBE2V1, RPS5, TRIM28 and
LIMA1. Four of them can be drawn and four cannot, and the reason is worth stating rather
than quietly dropping them: when all 271 measured risk genes were queried against the
CosMx object, 121 were found on its panel and 150 were not. DDX39B, CLIC1, RPS5 and LIMA1
are in the second group - the CosMx SMI panel carries 6,175 genes, not the transcriptome,
and these four are not among them. No amount of reprocessing recovers them; a different
platform would be needed.

The four that are on the panel - CSNK2B, CTNNB1, UBE2V1, TRIM28 - are all risk-increasing
by TWAS, all oncogene-like by the spatial validation of section 2.14, and all rise in
carcinoma cells. Each gets two panels rather than the one Figure 18 gave:

  top     the per-cell spatial map, as in Figure 18: every cell of slide 2 positioned as
          it sits in the tissue and coloured by its normalised expression of that gene.
  bottom  the patient-level evidence behind the map, which a spatial map cannot show -
          each of the seven patients carrying both populations as one line from its
          normal epithelium to its carcinoma cells, with the paired Wilcoxon p.

The map is the picture; the dot plot is the test. Showing them together stops the map
from being read as evidence on its own.

Outputs
  results/cosmx_requested_genes.csv
  figures/fig28_cosmx_requested_genes.{png,pdf}
"""
import pandas as pd, numpy as np, json, os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.colors import LinearSegmentedColormap
from scipy.stats import wilcoxon

os.makedirs("figures", exist_ok=True)
INK, INK2, INK3 = "#0b0b0b", "#52514e", "#8a8880"
SURF = "#fcfcfb"
EXPCM = LinearSegmentedColormap.from_list("exp", ["#e9e7e1", "#ffd9a0", "#f59b42",
                                                  "#e34948", "#8f1d2e"])
NORMC, CRCC = "#2a78d6", "#c3202f"
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

REQUESTED = ["DDX39B", "CSNK2B", "CLIC1", "CTNNB1", "UBE2V1", "RPS5", "TRIM28", "LIMA1"]
NORMAL, TUMOUR = "normal/undefined epithelium", "CRC"

# ---------------------------------------------------------------- what is available
SRC = {}
for f in ["spatial/cx_plot_cells.npz", "spatial/cx_plot_cells3.npz",
          "spatial/cx_plot_cells2.npz"]:
    z = np.load(f, allow_pickle=True)
    for g in [str(x) for x in z["genes"]]:
        SRC.setdefault(g, f)
panel = {l.strip() for l in open("spatial/cx_panel_genes.txt") if l.strip()}
GENES = [g for g in REQUESTED if g in SRC]
OFF_PANEL = [g for g in REQUESTED if g not in panel]
NO_CELLS = [g for g in REQUESTED if g in panel and g not in SRC]
print("drawable:", GENES)
print("not on the CosMx panel:", OFF_PANEL)
print("on the panel but no per-cell extract:", NO_CELLS)

B = pd.read_csv("results/cosmx_crc_vs_normal_epi.csv").set_index("gene")
V = pd.read_csv("results/validated_risk_genes_all.csv").set_index("gene")
P = pd.read_csv("spatial/cx_pseudobulk.csv.gz")
P = P[P.cellType.isin([NORMAL, TUMOUR])]

# ---------------------------------------------------------------- patient-level pairs
pairs = {}
rows = []
for g in GENES:
    sub = P[P.gene == g].pivot_table(index="patient_ID", columns="cellType",
                                     values="mean_norm")
    sub = sub.dropna(subset=[NORMAL, TUMOUR])
    pairs[g] = sub
    w = wilcoxon(sub[TUMOUR], sub[NORMAL]) if len(sub) >= 3 else None
    d = B.loc[g] if g in B.index else None
    v = V.loc[g] if g in V.index else None
    rows.append({"gene": g, "n_patients": len(sub),
                 "mean_normal_epi": float(sub[NORMAL].mean()),
                 "mean_CRC": float(sub[TUMOUR].mean()),
                 "delta": float((sub[TUMOUR] - sub[NORMAL]).mean()),
                 "p_wilcoxon": float(w.pvalue) if w else np.nan,
                 "n_up": int((sub[TUMOUR] > sub[NORMAL]).sum()),
                 "twas_direction": (v["twas_direction"] if v is not None else ""),
                 "class": (v["class"] if v is not None else ""),
                 "geomx_lead_FDR": (float(v["lead_FDR"]) if v is not None else np.nan),
                 "on_cosmx_panel": True})
for g in OFF_PANEL:
    rows.append({"gene": g, "on_cosmx_panel": False})
R = pd.DataFrame(rows)
R.to_csv("results/cosmx_requested_genes.csv", index=False)
print(R.round(4).to_string(index=False))

# ---------------------------------------------------------------- figure
n = len(GENES)
fig = plt.figure(figsize=(3.35 * n, 7.6))
gs = fig.add_gridspec(2, n, height_ratios=[1.0, 0.72], hspace=0.40, wspace=0.34)
letters = "abcdefgh"

for c, g in enumerate(GENES):
    # ---- top: the spatial map
    z = np.load(SRC[g], allow_pickle=True)
    X = z["x"].astype(float)
    Y = z["y"].astype(float)
    Y = Y.max() - Y
    E = (z["expr_q"].astype(np.float32) / float(z["qscale"])) if "expr_q" in z \
        else z["expr"].astype(np.float32)
    j = [str(x) for x in z["genes"]].index(g)
    e = E[:, j]

    ax = fig.add_subplot(gs[0, c])
    neg = e <= 0
    ax.scatter(X[neg], Y[neg], s=0.10, c="#eceae5", linewidths=0, rasterized=True)
    pos = ~neg
    o = np.argsort(e[pos])
    vmax = float(np.quantile(e[pos], 0.98)) if pos.sum() > 50 else 1.0
    sc = ax.scatter(X[pos][o], Y[pos][o], s=0.34, c=e[pos][o], cmap=EXPCM,
                    vmin=0, vmax=max(vmax, 0.2), linewidths=0, rasterized=True)
    ax.set_aspect("equal"); ax.axis("off")
    d = B.loc[g]
    v = V.loc[g]
    ax.text(0.0, 1.10, g, transform=ax.transAxes, fontsize=10.4, fontweight="bold",
            color=INK, va="baseline", ha="left")
    ax.text(0.0, 1.035,
            f"CosMx Δ {d['delta']:+.3f}, P = {d['p']:.3f} (n = {int(d['n_patients'])})",
            transform=ax.transAxes, fontsize=7.2, color=INK2, va="baseline", ha="left")
    ax.text(0.0, -0.015,
            f"GeoMx FDR {v['lead_FDR']:.1g} · {v['class'].replace('tumour-suppressor-like', 'TSG-like')}",
            transform=ax.transAxes, fontsize=7, color="#1baf7a", va="top", ha="left")
    ax.text(-0.085, 1.20, letters[c], transform=ax.transAxes, fontsize=11,
            fontweight="bold", color=INK, va="top", ha="left")
    cb = fig.colorbar(sc, ax=ax, fraction=0.030, pad=0.005, shrink=0.62,
                      location="left", anchor=(0.0, 0.06))
    cb.ax.tick_params(labelsize=6.0, pad=1); cb.outline.set_visible(False)
    cb.set_label("log$_2$ norm.", fontsize=6.0, labelpad=2)

    # ---- bottom: the patient-level pairs behind it
    ax = fig.add_subplot(gs[1, c])
    ax.text(-0.26, 1.16, letters[c + n], transform=ax.transAxes, fontsize=11,
            fontweight="bold", color=INK, va="top", ha="left")
    sub = pairs[g]
    rng = np.random.default_rng(47 + c)
    for _, r in sub.iterrows():
        up = r[TUMOUR] > r[NORMAL]
        ax.plot([0, 1], [r[NORMAL], r[TUMOUR]], color="#c9c7c1" if not up else "#e6b8b8",
                lw=1.0, zorder=2)
    ax.scatter(np.zeros(len(sub)) + rng.normal(0, 0.012, len(sub)), sub[NORMAL], s=32,
               color=NORMC, edgecolor="white", linewidth=0.6, zorder=4)
    ax.scatter(np.ones(len(sub)) + rng.normal(0, 0.012, len(sub)), sub[TUMOUR], s=32,
               color=CRCC, edgecolor="white", linewidth=0.6, zorder=4)
    for xx, vals, col in [(0, sub[NORMAL], NORMC), (1, sub[TUMOUR], CRCC)]:
        ax.plot([xx - 0.16, xx + 0.16], [vals.mean()] * 2, color=col, lw=2.2, zorder=5)
    rr = R[R.gene == g].iloc[0]
    ax.set_xticks([0, 1])
    ax.set_xticklabels(["normal\nepithelium", "carcinoma\ncells"], fontsize=7.8)
    ax.set_xlim(-0.35, 1.35)
    if c == 0:
        ax.set_ylabel("Mean normalised expression\nper patient", fontsize=8)
    ax.set_title(f"paired Wilcoxon P = {rr.p_wilcoxon:.3f}", fontsize=8.8, pad=16,
                 color="#1baf7a" if rr.p_wilcoxon < 0.05 else INK2, loc="center")
    ax.text(0.5, 1.012, f"{int(rr.n_up)} of {int(rr.n_patients)} patients rise",
            transform=ax.transAxes, ha="center", va="bottom", fontsize=7.6, color=INK2)
    ax.grid(axis="y"); ax.set_axisbelow(True)

fig.text(0.5, 1.005,
         f"Slide {int(z['slide'])}, {int(z['n_cells_slide']):,} cells on the slide. "
         "Top row: every cell plotted where it sits, coloured by expression. "
         "Bottom row: the same contrast as a paired test, one line per patient.",
         ha="center", fontsize=8.2, color=INK2)
miss = ", ".join(OFF_PANEL)
fig.text(0.5, -0.045,
         f"{miss} are not drawn because they are not on the CosMx panel. The panel carries "
         "6,175 genes rather than the transcriptome: of the 271 risk genes measured in "
         "GeoMx, 121 are present on it and 150, including these four,\nare absent. That is "
         "a property of the assay, not a result — nothing about these genes was tested and "
         "found wanting in the single-cell data, and confirming their spatial behaviour "
         "would need a different platform.",
         ha="center", fontsize=7.6, color=INK2)

fig.savefig("figures/fig28_cosmx_requested_genes.png", dpi=220, bbox_inches="tight")
fig.savefig("figures/fig28_cosmx_requested_genes.pdf", bbox_inches="tight")
plt.close(fig)
print("\nfig28 done")
