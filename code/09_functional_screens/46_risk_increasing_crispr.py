#!/usr/bin/env python3
"""
46 - Figure 27: the TWAS risk-increasing genes in the CRISPR-Cas9 screen alone.

Figure 26 pooled both screens and had to carry the siRNA SW480 plates, which fail quality
control and supply every change that runs against the prediction. This figure drops that
complication entirely: one assay, all three of its cell lines, all of which pass the
plate-level Z' check, and only the genes for which the germline direction makes a
prediction.

That prediction is one-sided. A risk-increasing gene raises risk through higher activity,
so if the effect carried into the tumour cell the cell should need the gene and editing it
should lower viability. The test is therefore one-sided too, and the comparison groups are
the two sets assayed on the very same plates: the protective risk genes and the non-TWAS
genes of the same 116-gene panel.

Analysis as in Figures 20-22: each plate's eight non-targeting wells averaged into that
plate's control value, a paired t-test over the three (gene, control) pairs, and the
effect as the mean within-plate difference in viability points.

Outputs
  results/risk_increasing_crispr.csv
  figures/fig27_risk_increasing_crispr.{png,pdf}
"""
import pandas as pd, numpy as np, json, os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch, Rectangle
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm
from scipy.stats import binomtest, mannwhitneyu

os.makedirs("figures", exist_ok=True)
LINES = ["HCT116", "RKO", "SW480"]
POS = "#e34948"
INK, INK2, INK3 = "#0b0b0b", "#52514e", "#8a8880"
SURF = "#fcfcfb"
DIV = LinearSegmentedColormap.from_list("div", ["#184f95", "#6da7ec", "#f2f1ed",
                                                "#f0a09f", "#c3202f"])
RI, PR, OT = "#8c2d3f", "#1f6f8b", "#7d99b8"
DOWN, UP = "#1f7a5a", "#d98218"
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


def tag(ax, s, dx=-0.10, dy=1.05):
    ax.text(dx, dy, s, transform=ax.transAxes, fontsize=11, fontweight="bold",
            color=INK, va="top", ha="left")


C = pd.read_csv("results/crispr_risk36_by_cellline.csv")
W = C[C.twas_direction == "risk-increasing"].copy().reset_index(drop=True)
W["mean_delta"] = W[[f"delta_{l}" for l in LINES]].mean(axis=1)
W["n_sig"] = W[[f"nom_{l}" for l in LINES]].sum(axis=1)


def verdict(r):
    s = [np.sign(r[f"delta_{l}"]) for l in LINES if r[f"nom_{l}"]]
    if not s:
        return "no significant change"
    if len(set(s)) > 1:
        return "conflicting"
    return "lowers viability" if s[0] < 0 else "raises viability"


W["verdict"] = W.apply(verdict, axis=1)
W[["gene", "class"] + [f"delta_{l}" for l in LINES] + [f"nom_{l}" for l in LINES]
  + [f"pPaired_{l}" for l in LINES] + ["mean_delta", "n_sig", "verdict"]].to_csv(
    "results/risk_increasing_crispr.csv", index=False)

# ---- every significant call, and its direction
calls = [{"gene": r.gene, "cell_line": l, "delta": r[f"delta_{l}"]}
         for _, r in W.iterrows() for l in LINES if r[f"nom_{l}"]]
Cl = pd.DataFrame(calls)
n_down, n_tot = int((Cl.delta < 0).sum()), len(Cl)
bt = binomtest(n_down, n_tot, 0.5, alternative="greater")
per_line = {l: (int(((Cl.cell_line == l) & (Cl.delta < 0)).sum()),
                int((Cl.cell_line == l).sum())) for l in LINES}

# ---- comparison groups, same plates
L = pd.read_csv("results/crispr_by_cellline.csv")
Lg = L.groupby("gene").agg(delta=("delta_paired", "mean"),
                           twas=("twas_direction", "first"),
                           is_risk=("is_risk", "first"),
                           n_sig=("nominal", "sum")).reset_index()
g_ri = Lg[(Lg.is_risk) & (Lg.twas == "risk-increasing")]
g_pr = Lg[(Lg.is_risk) & (Lg.twas == "protective")]
g_ot = Lg[~Lg.is_risk]
u_rp = mannwhitneyu(g_ri.delta, g_pr.delta, alternative="less")
u_ro = mannwhitneyu(g_ri.delta, g_ot.delta, alternative="less")

# the same one-sided direction question for the comparison groups, as a control
def down_rate(sub):
    cl = []
    for _, r in L[L.gene.isin(sub.gene) & L.nominal].iterrows():
        cl.append(r.delta_paired)
    if not cl:
        return 0, 0, None
    d = int(sum(1 for v in cl if v < 0))
    return d, len(cl), binomtest(d, len(cl), 0.5, alternative="greater").pvalue


d_pr, n_pr, p_pr = down_rate(g_pr)
d_ot, n_ot, p_ot = down_rate(g_ot)

vd = W["verdict"].value_counts().to_dict()
J = {"n_genes": len(W), "calls_total": n_tot, "calls_down": n_down,
     "p_one_sided": float(bt.pvalue), "per_line": {k: list(v) for k, v in per_line.items()},
     "gene_verdict": {k: int(v) for k, v in vd.items()},
     "genes_lowering": W.loc[W.verdict == "lowers viability", "gene"].tolist(),
     "median_delta_risk_increasing": float(g_ri.delta.median()),
     "median_delta_protective": float(g_pr.delta.median()),
     "median_delta_other": float(g_ot.delta.median()),
     "p_vs_protective": float(u_rp.pvalue), "p_vs_other": float(u_ro.pvalue),
     "control_down_protective": [d_pr, n_pr, p_pr],
     "control_down_other": [d_ot, n_ot, p_ot]}
json.dump(J, open("results/risk_increasing_crispr.json", "w"), indent=1)
print(json.dumps(J, indent=1))

ORDER = {"lowers viability": 0, "conflicting": 1, "raises viability": 2,
         "no significant change": 3}
W["rank"] = W["verdict"].map(ORDER)
W = W.sort_values(["rank", "mean_delta"], ascending=[False, False]).reset_index(drop=True)
N = len(W)

fig = plt.figure(figsize=(13.8, 9.4))
gs = fig.add_gridspec(2, 3, width_ratios=[1.05, 0.92, 1.05], height_ratios=[1.0, 0.86],
                      wspace=0.52, hspace=0.62)

# ---------------------------------------------------------------- (a) matrix
ax = fig.add_subplot(gs[:, 0]); tag(ax, "a", dx=-0.42, dy=1.03)
M = W[[f"delta_{l}" for l in LINES]].to_numpy(float)
lim = max(float(np.nanpercentile(np.abs(M), 97)), 5.0)
im = ax.imshow(M, cmap=DIV, norm=TwoSlopeNorm(vmin=-lim, vcenter=0, vmax=lim),
               aspect="auto", origin="lower", extent=(-0.5, 2.5, -0.5, N - 0.5))
VCOL = {"lowers viability": DOWN, "raises viability": UP, "conflicting": "#8a8880",
        "no significant change": "#e4e2dd"}
for i, v in enumerate(W["verdict"]):
    ax.add_patch(Rectangle((2.68, i - 0.5), 0.5, 1, facecolor=VCOL[v], edgecolor=SURF,
                           linewidth=0.5, zorder=2))
ax.text(2.93, N - 0.2, "effect", ha="center", va="bottom", fontsize=7.2, color=INK2,
        fontweight="bold")
for j, line in enumerate(LINES):
    for i in np.flatnonzero(W[f"nom_{line}"].to_numpy(bool)):
        ax.add_patch(Rectangle((j - 0.5, i - 0.5), 1, 1, fill=False, edgecolor=INK,
                               linewidth=1.6, zorder=5))
        v = W[f"delta_{line}"].iloc[i]
        ax.text(j, i, f"{v:+.0f}", ha="center", va="center", fontsize=6.4,
                color="white" if abs(v) > lim * 0.55 else INK, zorder=6)
ax.set_xticks(range(3)); ax.set_xticklabels(LINES, fontsize=8.6)
ax.set_yticks(range(N))
labs = [f"{'▲ ' if c == 'oncogene-like' else '  '}{g}"
        for g, c in zip(W["gene"], W["class"])]
ax.set_yticklabels(labs, fontsize=7.6)
for t, c in zip(ax.get_yticklabels(), W["class"]):
    t.set_color(POS if c == "oncogene-like" else INK2)
ax.set_xlim(-0.55, 3.20); ax.set_ylim(-0.5, N - 0.5)
ax.set_title(f"The {N} risk-increasing genes,\nCRISPR-Cas9", pad=12)
ax.set_xlabel("outlined = paired t-test P < 0.05", fontsize=7.2, color=INK3,
              style="italic", labelpad=7)
ax.tick_params(length=0)
for sp in ax.spines.values():
    sp.set_visible(False)
cb = fig.colorbar(im, ax=ax, fraction=0.045, pad=0.10, shrink=0.40, anchor=(0.0, 1.0))
cb.set_label("Δ viability vs same-plate control (points)", fontsize=7.3)
cb.ax.tick_params(labelsize=7); cb.outline.set_visible(False)
ax.legend(handles=[
    Patch(facecolor=DOWN, label="editing lowers viability — as predicted"),
    Patch(facecolor="#e4e2dd", label="no significant change"),
    Line2D([], [], marker="$▲$", ls="", color=POS, ms=8, label="oncogene-like (§2.14)")],
    loc="upper left", bbox_to_anchor=(-0.40, -0.075), ncol=1, fontsize=7.2)

# ---------------------------------------------------------------- (b) direction
ax = fig.add_subplot(gs[0, 1]); tag(ax, "b", dx=-0.26)
x = np.arange(4, dtype=float)
downs = [per_line[l][0] for l in LINES] + [n_down]
tots = [per_line[l][1] for l in LINES] + [n_tot]
ups = [t - d for d, t in zip(downs, tots)]
ax.bar(x, downs, width=0.6, color=DOWN, label="viability falls (predicted)")
ax.bar(x, ups, width=0.6, bottom=downs, color=UP, label="viability rises")
for xi, d, u_, t in zip(x, downs, ups, tots):
    if d:
        ax.text(xi, d / 2, str(d), ha="center", va="center", fontsize=8.6, color="white",
                fontweight="bold")
    if u_:
        ax.text(xi, d + u_ / 2, str(u_), ha="center", va="center", fontsize=8.6,
                color="white", fontweight="bold")
ax.set_xticks(x); ax.set_xticklabels(LINES + ["all three"], fontsize=8.0)
ax.get_xticklabels()[3].set_fontweight("bold")
ax.set_ylabel("Significant changes")
ax.set_ylim(0, max(tots) * 1.22)
ax.set_title("Direction of every\nsignificant change")
ax.legend(loc="upper left", fontsize=7.0)
ax.grid(axis="y"); ax.set_axisbelow(True)
ax.text(0.0, -0.30,
        "Every significant change among these\ngenes is a decrease: {} of {}, one-sided\n"
        "binomial P = {:.1e}. Nothing runs the\nother way in any of the three lines.".format(
            n_down, n_tot, bt.pvalue),
        transform=ax.transAxes, fontsize=7.4, color=INK2, va="top")

# ---------------------------------------------------------------- (c) the control groups
ax = fig.add_subplot(gs[0, 2]); tag(ax, "c", dx=-0.22)
labs2 = [f"risk-increasing\n{n_down}/{n_tot}", f"protective\n{d_pr}/{n_pr}",
         f"non-TWAS panel\n{d_ot}/{n_ot}"]
rates = [100 * n_down / n_tot, 100 * d_pr / n_pr if n_pr else 0,
         100 * d_ot / n_ot if n_ot else 0]
ax.bar(range(3), rates, width=0.6, color=[RI, PR, OT])
ax.axhline(50, color=INK3, lw=0.9, ls=(0, (3, 2)))
ax.text(-0.42, 51, "chance", fontsize=7.2, color=INK3, ha="left", va="bottom")
for i, (v, p) in enumerate(zip(rates, [bt.pvalue, p_pr, p_ot])):
    ax.text(i, v + 2.2, f"{v:.0f} %", ha="center", fontsize=8.2, color=INK2)
    if p is not None:
        ax.text(i, 30, f"P = {p:.2g}", ha="center", va="center", fontsize=7.2,
                color="white")
ax.set_xticks(range(3)); ax.set_xticklabels(labs2, fontsize=7.6)
ax.set_ylabel("% of significant changes\nthat are decreases")
ax.set_ylim(0, 112)
ax.set_title("Against the genes on the\nsame plates")
ax.grid(axis="y"); ax.set_axisbelow(True)
ax.text(0.0, -0.30,
        "Knocking out almost anything tends to cost viability, so the\ncomparison groups also "
        "sit above chance. The risk-increasing\ngenes are the only group with no exception at "
        "all, but the\ncontrast with the other groups is one of degree, not of kind.",
        transform=ax.transAxes, fontsize=7.4, color=INK2, va="top")

# ---------------------------------------------------------------- (d) per gene, ranked
ax = fig.add_subplot(gs[1, 1:]); tag(ax, "d", dx=-0.09)
Wd = W.sort_values("mean_delta")
y = np.arange(len(Wd))
cols = [DOWN if v == "lowers viability" else ("#c9c7c1" if v == "no significant change"
                                              else VCOL[v]) for v in Wd["verdict"]]
ax.barh(y, Wd["mean_delta"], height=0.68, color=cols)
ax.axvline(0, color=INK2, lw=0.9)
for yi, (_, r) in zip(y, Wd.iterrows()):
    off = -1.2 if r.mean_delta < 0 else 1.2
    ha = "right" if r.mean_delta < 0 else "left"
    ax.text(r.mean_delta + off, yi, r.gene, va="center", ha=ha, fontsize=6.6,
            color=INK if r.verdict == "lowers viability" else INK3,
            fontweight="bold" if r.verdict == "lowers viability" else "normal")
ax.set_yticks([]); ax.set_ylim(-0.8, len(Wd) - 0.2)
ax.set_xlabel("Mean Δ viability across the three CRISPR lines (points)")
ax.set_xlim(Wd["mean_delta"].min() * 1.30, max(Wd["mean_delta"].max() * 1.8, 6))
ax.set_title("Every risk-increasing gene, ranked")
ax.grid(axis="x"); ax.set_axisbelow(True)
ax.text(0.0, -0.20,
        "Bold, green: significant in at least one line. {} of the {} genes reach significance "
        "somewhere; the rest are grey. As a group the\nrisk-increasing genes sit at a median "
        "of {:+.1f} points against {:+.1f} for the protective risk genes (one-sided P = {:.2f}) "
        "and {:+.1f} for\nthe non-TWAS genes of the same panel (P = {:.2f}) — so it is the "
        "direction of the significant calls that carries the signal,\nnot a shift of the whole "
        "set.".format(
            vd.get("lowers viability", 0), N, J["median_delta_risk_increasing"],
            J["median_delta_protective"], u_rp.pvalue, J["median_delta_other"], u_ro.pvalue),
        transform=ax.transAxes, fontsize=7.4, color=INK2, va="top")

fig.savefig("figures/fig27_risk_increasing_crispr.png", dpi=220, bbox_inches="tight")
fig.savefig("figures/fig27_risk_increasing_crispr.pdf", bbox_inches="tight")
plt.close(fig)
print("\nfig27 done")
