#!/usr/bin/env python3
"""
44 - Figure 25: the 30 risk genes that appear in Figure 20a AND on the siRNA panel,
with both modalities side by side on one row per gene.

Figures 20 and 21 each answer "what happens in this screen". Neither answers the question
a reader actually has, which is whether the two screens say the same thing about the same
gene. Restricting to the intersection makes that comparison exact: 30 of the 36 genes in
Figure 20a were also assayed by siRNA (the six that were not are ABHD11, DIP2B, PLEKHG6,
PYGL, STK39 and TNS3), and every one of the 30 has three CRISPR plates and three siRNA
plates behind it.

Same analysis as Figures 20-22 throughout: each plate's non-targeting wells averaged into
that plate's control value, a paired t-test over the three (gene, control) pairs, and the
effect as the mean within-plate difference in viability points. One colour scale spans
both assays, so the panels are directly comparable and the larger siRNA effects are
visible as larger rather than being normalised away.

Outputs
  results/shared_genes_both_assays.csv
  figures/fig25_shared_genes.{png,pdf}
"""
import pandas as pd, numpy as np, json, os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch, Rectangle
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm
from scipy.stats import spearmanr, fisher_exact

os.makedirs("figures", exist_ok=True)
LINES = ["HCT116", "RKO", "SW480"]
POS, NEG = "#e34948", "#2a78d6"
INK, INK2, INK3 = "#0b0b0b", "#52514e", "#8a8880"
SURF = "#fcfcfb"
DIV = LinearSegmentedColormap.from_list("div", ["#184f95", "#6da7ec", "#f2f1ed",
                                                "#f0a09f", "#c3202f"])
TW = {"risk-increasing": "#8c2d3f", "protective": "#1f6f8b"}
SUPPORTS, OPPOSES, NS = "supports TWAS", "opposes TWAS", "not significant"
VC = {SUPPORTS: "#1f7a5a", OPPOSES: "#d98218", "conflicting": "#8a8880", NS: "#e4e2dd"}
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


C = pd.read_csv("results/crispr_risk36_by_cellline.csv").set_index("gene")
S = pd.read_csv("results/sirna_risk_by_cellline.csv").set_index("gene")
shared = sorted(set(C.index) & set(S.index))
only_crispr = sorted(set(C.index) - set(S.index))

W = pd.DataFrame(index=shared)
for line in LINES:
    W[f"cr_{line}"] = C.loc[shared, f"delta_{line}"]
    W[f"cr_sig_{line}"] = C.loc[shared, f"nom_{line}"].astype(bool)
    W[f"si_{line}"] = S.loc[shared, f"delta_{line}"]
    W[f"si_sig_{line}"] = S.loc[shared, f"nom_{line}"].astype(bool)
W["twas_direction"] = C.loc[shared, "twas_direction"]
W["class"] = C.loc[shared, "class"]
W["cr_verdict"] = C.loc[shared, "twas_support"]
W["si_verdict"] = S.loc[shared, "twas_support"]
# effects for the scatter: the two QC-passing lines, which both assays share
W["cr_mean"] = W[["cr_HCT116", "cr_RKO"]].mean(axis=1)
W["si_mean"] = W[["si_HCT116", "si_RKO"]].mean(axis=1)
W["cr_any"] = W[[f"cr_sig_{l}" for l in LINES]].any(axis=1)
W["si_any"] = W[[f"si_sig_{l}" for l in ["HCT116", "RKO"]]].any(axis=1)


def both_call(r):
    if r.cr_any and r.si_any:
        return "both modalities"
    if r.cr_any:
        return "CRISPR only"
    if r.si_any:
        return "siRNA only"
    return "neither"


W["called_in"] = W.apply(both_call, axis=1)
W["same_sign"] = np.sign(W["cr_mean"]) == np.sign(W["si_mean"])
W = W.reset_index().rename(columns={"index": "gene"})
W.to_csv("results/shared_genes_both_assays.csv", index=False)

rho, prho = spearmanr(W["cr_mean"], W["si_mean"])
counts = W["called_in"].value_counts().to_dict()
tab = [[int(((W.cr_any) & (W.si_any)).sum()), int(((W.cr_any) & (~W.si_any)).sum())],
       [int(((~W.cr_any) & (W.si_any)).sum()), int(((~W.cr_any) & (~W.si_any)).sum())]]
orr, pf = fisher_exact(tab)
print(f"shared genes: {len(W)} (CRISPR-only genes not shown: {', '.join(only_crispr)})")
print(f"Spearman on the two-line mean effect: rho = {rho:+.3f} (P = {prho:.3g}); "
       f"same sign {int(W.same_sign.sum())}/{len(W)}")
print("called in:", counts)
print(f"2x2 CRISPR x siRNA significance: {tab}, OR {orr:.2f}, Fisher P = {pf:.3g}")
print(W[["gene", "twas_direction", "cr_mean", "si_mean", "cr_verdict", "si_verdict",
         "called_in"]].round(2).to_string(index=False))

# order: genes called by both first, then one, then neither; within that by siRNA effect
RANK = {"both modalities": 0, "CRISPR only": 1, "siRNA only": 2, "neither": 3}
W["rank"] = W["called_in"].map(RANK)
W = W.sort_values(["rank", "si_mean"], ascending=[False, False]).reset_index(drop=True)
N = len(W)

fig = plt.figure(figsize=(14.0, 9.8))
gs = fig.add_gridspec(2, 2, width_ratios=[1.45, 1.0], height_ratios=[1.0, 0.92],
                      wspace=0.40, hspace=0.62)

# ---------------------------------------------------------------- (a) the two matrices
ax = fig.add_subplot(gs[:, 0]); tag(ax, "a", dx=-0.30, dy=1.03)
M = W[[f"cr_{l}" for l in LINES] + [f"si_{l}" for l in LINES]].to_numpy(float)
lim = float(np.nanpercentile(np.abs(M), 97))
norm = TwoSlopeNorm(vmin=-lim, vcenter=0, vmax=lim)
im = ax.imshow(M, cmap=DIV, norm=norm, aspect="auto", origin="lower",
               extent=(-0.5, 5.5, -0.5, N - 0.5))
ax.axvline(2.5, color=SURF, lw=5, zorder=6)
ax.axvline(2.5, color=INK3, lw=1.0, zorder=7)

for i, d in enumerate(W["twas_direction"].fillna("")):
    ax.add_patch(Rectangle((-1.30, i - 0.5), 0.55, 1, facecolor=TW.get(d, "#d8d7d2"),
                           edgecolor=SURF, linewidth=0.5, zorder=2))
for k, (col, x0) in enumerate([("cr_verdict", 5.70), ("si_verdict", 6.30)]):
    for i, v in enumerate(W[col].fillna(NS)):
        ax.add_patch(Rectangle((x0, i - 0.5), 0.5, 1, facecolor=VC.get(v, VC[NS]),
                               edgecolor=SURF, linewidth=0.5, zorder=2))
ax.text(-1.02, N - 0.2, "TWAS", ha="center", va="bottom", fontsize=7.2, color=INK2,
        fontweight="bold")
ax.text(6.25, N - 0.2, "verdict", ha="center", va="bottom", fontsize=7.2, color=INK2,
        fontweight="bold")
for xc, lb in [(5.95, "CRISPR"), (6.55, "siRNA")]:
    ax.text(xc, -0.012, lb, ha="center", va="top", fontsize=6.6, color=INK3,
            rotation=90, transform=ax.get_xaxis_transform())

for j, (pre, line) in enumerate([(p, l) for p in ["cr", "si"] for l in LINES]):
    sig = W[f"{pre}_sig_{line}"].to_numpy(bool)
    for i in np.flatnonzero(sig):
        ax.add_patch(Rectangle((j - 0.5, i - 0.5), 1, 1, fill=False, edgecolor=INK,
                               linewidth=1.5, zorder=5))
        v = W[f"{pre}_{line}"].iloc[i]
        ax.text(j, i, f"{v:+.0f}", ha="center", va="center", fontsize=6.2,
                color="white" if abs(v) > lim * 0.55 else INK, zorder=8)

ax.set_xticks(range(6)); ax.set_xticklabels(LINES + LINES, fontsize=7.8)
ax.text(1, -0.085, "CRISPR-Cas9", ha="center", va="top", fontsize=9, color=INK,
        fontweight="bold", transform=ax.get_xaxis_transform())
ax.text(4, -0.085, "siRNA", ha="center", va="top", fontsize=9, color=INK,
        fontweight="bold", transform=ax.get_xaxis_transform())
ax.set_yticks(range(N))
labs = [f"{ {'oncogene-like': '▲', 'tumour-suppressor-like': '▼'}.get(c, ' ') } {g}"
        for g, c in zip(W["gene"], W["class"])]
ax.set_yticklabels(labs, fontsize=7.4)
for t, cls in zip(ax.get_yticklabels(), W["class"]):
    t.set_color(POS if cls == "oncogene-like"
                else (NEG if cls == "tumour-suppressor-like" else INK2))
ax.set_xlim(-1.40, 6.90); ax.set_ylim(-0.5, N - 0.5)
ax.set_title(f"The {N} risk genes assayed by both modalities", pad=14)
ax.tick_params(length=0)
for sp in ax.spines.values():
    sp.set_visible(False)
cb = fig.colorbar(im, ax=ax, fraction=0.034, pad=0.075, shrink=0.40, anchor=(0.0, 1.0))
cb.set_label("Δ viability vs same-plate control (points)", fontsize=7.4)
cb.ax.tick_params(labelsize=7); cb.outline.set_visible(False)
ax.legend(handles=[
    Patch(facecolor=TW["risk-increasing"], label="TWAS risk-increasing"),
    Patch(facecolor=VC[SUPPORTS], label="supports TWAS"),
    Patch(facecolor=TW["protective"], label="TWAS protective"),
    Patch(facecolor=VC[OPPOSES], label="opposes TWAS"),
    Line2D([], [], marker="s", ls="", mfc="none", mec=INK, ms=9, markeredgewidth=1.6,
           label="paired t-test P < 0.05"),
    Patch(facecolor=VC[NS], label="no significant change")],
    loc="upper left", bbox_to_anchor=(-0.33, -0.10), ncol=2, fontsize=7.2,
    columnspacing=1.0)

# ---------------------------------------------------------------- (b) scatter
ax = fig.add_subplot(gs[0, 1]); tag(ax, "b")
ax.axhline(0, color="#d8d7d2", lw=0.9); ax.axvline(0, color="#d8d7d2", lw=0.9)
COL = {"both modalities": "#1f5fae", "CRISPR only": "#8c2d3f", "siRNA only": "#d98218",
       "neither": "#c9c7c1"}
for k in ["neither", "siRNA only", "CRISPR only", "both modalities"]:
    m = W["called_in"] == k
    if m.any():
        ax.scatter(W.loc[m, "si_mean"], W.loc[m, "cr_mean"], s=44 if k != "neither" else 24,
                   color=COL[k], edgecolor="white", linewidth=0.5,
                   zorder=4 if k != "neither" else 3, label=f"{k} ({int(m.sum())})")
lims = [min(W.si_mean.min(), W.cr_mean.min()) - 6, max(W.si_mean.max(), W.cr_mean.max()) + 6]
ax.plot(lims, lims, color=INK3, lw=0.8, ls=(0, (3, 2)), zorder=2)
for _, r in W.iterrows():
    if r.called_in == "both modalities" or abs(r.si_mean) > 22 or abs(r.cr_mean) > 9:
        ax.annotate(r.gene, (r.si_mean, r.cr_mean), textcoords="offset points",
                    xytext=(4, 3), fontsize=6.6, color=INK2)
ax.set_xlabel("siRNA Δ viability (HCT116 + RKO, points)")
ax.set_ylabel("CRISPR Δ viability\n(HCT116 + RKO, points)")
ax.set_title("The two modalities, gene by gene")
ax.text(0.03, 0.06, f"Spearman ρ = {rho:+.2f}, P = {prho:.2f}\n"
                    f"same sign {int(W.same_sign.sum())}/{len(W)}",
        transform=ax.transAxes, fontsize=7.6, color=INK2, va="bottom")
ax.legend(loc="lower right", fontsize=7.0)
ax.grid(True); ax.set_axisbelow(True)

# ---------------------------------------------------------------- (c) agreement
ax = fig.add_subplot(gs[1, 1]); tag(ax, "c")
order = ["both modalities", "CRISPR only", "siRNA only", "neither"]
vals = [counts.get(k, 0) for k in order]
ax.bar(range(4), vals, width=0.62, color=[COL[k] for k in order])
for i, v in enumerate(vals):
    ax.text(i, v + 0.35, str(v), ha="center", fontsize=8.4, color=INK2)
ax.set_xticks(range(4))
ax.set_xticklabels(["both", "CRISPR\nonly", "siRNA\nonly", "neither"], fontsize=8.2)
ax.set_ylabel(f"Genes (of {N})")
ax.set_ylim(0, max(vals) * 1.25)
ax.set_title("Where each gene is called")
ax.grid(axis="y"); ax.set_axisbelow(True)
both = W[W.called_in == "both modalities"]
agree = int((np.sign(both.cr_mean) == np.sign(both.si_mean)).sum())
ax.text(0.0, -0.36,
        "A gene is counted for CRISPR if any of its three plates reaches P < 0.05, and for siRNA if "
        "either QC-passing line does.\nOnly {} of the {} genes are called in both, and {} of those {} "
        "move the same way — so the overlap is small but not\ncontradictory. Across all {} genes the "
        "two modalities correlate at ρ = {:+.2f} (P = {:.2f}) with {}/{} sharing a sign.\n"
        "Knockdown leaves residual protein and knockout removes the gene outright, so partial "
        "agreement is the expected result,\nnot a failure of either screen.".format(
            counts.get("both modalities", 0), N, agree, len(both), N, rho, prho,
            int(W.same_sign.sum()), N),
        transform=ax.transAxes, fontsize=7.4, color=INK2, va="top")

fig.savefig("figures/fig25_shared_genes.png", dpi=220, bbox_inches="tight")
fig.savefig("figures/fig25_shared_genes.pdf", bbox_inches="tight")
plt.close(fig)
print("\nfig25 done")
