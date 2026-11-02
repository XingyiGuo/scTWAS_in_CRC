#!/usr/bin/env python3
"""
43 - Figures 23 and 24: replicate concordance as a correlation heatmap, one figure per
assay. Run with no argument for CRISPR (Figure 23) or "sirna" for siRNA (Figure 24).

A composite of six quality metrics answers the exclusion question but does not let a
reader see the thing being judged. Correlation does: every plate against every other
plate, over the 116 genes they all measure. Agreement between the three plates of one
cell line is what a replicate is for, so the diagonal blocks are the ones to read; the
off-diagonal blocks are a bonus, showing how much of a gene's behaviour carries from one
cell line to another.

Correlations are computed on the effect scale used everywhere else - each gene's value
minus its own plate's non-targeting mean - so the plate offset is removed before
correlating. (Pearson r is unchanged by that subtraction; it is done for clarity.)

Both coefficients are shown because they answer different questions. Pearson is driven by
the large effects, so it reports whether the strong hits reproduce. Spearman uses ranks,
so it reports whether the whole panel, most of which barely moves, reproduces. Where the
two diverge, the plate agrees on its big effects and not on its small ones.

Outputs
  results/<assay>_plate_correlations.csv
  figures/fig23_crispr_concordance.{png,pdf}  /  figures/fig24_sirna_concordance.{png,pdf}
"""
import pandas as pd, numpy as np, json, os, sys
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.colors import LinearSegmentedColormap
from scipy.stats import pearsonr, spearmanr

os.makedirs("figures", exist_ok=True)
LINES = ["HCT116", "RKO", "SW480"]
REPS = ["rep1", "rep2", "rep3"]
INK, INK2, INK3 = "#0b0b0b", "#52514e", "#8a8880"
SURF = "#fcfcfb"
BAD = "#b3472f"
SEQ = LinearSegmentedColormap.from_list(
    "seq", ["#f7f6f3", "#dce8f5", "#9dc1e6", "#4a86c8", "#1f4f8f"])
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


def tag(ax, s, dx=-0.16, dy=1.06):
    ax.text(dx, dy, s, transform=ax.transAxes, fontsize=11, fontweight="bold",
            color=INK, va="top", ha="left")


ASSAY = (sys.argv[1] if len(sys.argv) > 1 else "crispr").lower()
CFG = {
    "crispr": {"src": "results/func_crispr_long.csv", "fig": "fig23_crispr_concordance",
               "out": "results/crispr_plate_correlations.csv", "name": "CRISPR-Cas9",
               "n_nt": 8, "qc_fail_line": None},
    "sirna": {"src": "results/func_sirna_long.csv", "fig": "fig24_sirna_concordance",
              "out": "results/sirna_plate_correlations.csv", "name": "siRNA",
              "n_nt": 12, "qc_fail_line": "SW480"},
}[ASSAY]
C = pd.read_csv(CFG["src"])
# no replicate is excluded in the current analysis; the shading is left available in case
# an exclusion is reinstated
EXCLUDED = set()

# effect per gene per plate: value minus that plate's non-targeting mean
cols, labels = {}, []
for line in LINES:
    d = C[C.cell_line == line]
    ntm = d[d.kind == "negative"].groupby("replicate")["value"].mean()
    G = d[d.kind == "sample"].pivot_table(index="gene", columns="replicate", values="value")
    for rep in REPS:
        cols[(line, rep)] = G[rep] - ntm[rep]
        labels.append((line, rep))
E = pd.DataFrame(cols)
n = len(labels)

P = np.ones((n, n)); S = np.ones((n, n))
rows = []
for i in range(n):
    for j in range(n):
        if i == j:
            continue
        a, b = E.iloc[:, i], E.iloc[:, j]
        P[i, j] = pearsonr(a, b).statistic
        S[i, j] = spearmanr(a, b).statistic
        if i < j:
            rows.append({"plate_a": " ".join(labels[i]), "plate_b": " ".join(labels[j]),
                         "same_cell_line": labels[i][0] == labels[j][0],
                         "pearson_r": P[i, j], "spearman_rho": S[i, j],
                         "n_genes": len(a)})
Cor = pd.DataFrame(rows)
Cor.to_csv(CFG["out"], index=False)

within = Cor[Cor.same_cell_line]
across = Cor[~Cor.same_cell_line]
NG = C[C.kind == "sample"].gene.nunique()
print(f"=== plate-plate correlation, {NG} genes ({CFG['name']}) ===")
print(Cor.round(3).to_string(index=False))
print(f"\nwithin cell line : median Pearson {within.pearson_r.median():.3f} "
      f"(range {within.pearson_r.min():.3f}-{within.pearson_r.max():.3f})")
print(f"across cell lines: median Pearson {across.pearson_r.median():.3f} "
      f"(range {across.pearson_r.min():.3f}-{across.pearson_r.max():.3f})")
for line in LINES:
    w = within[within.plate_a.str.startswith(line)]
    print(f"  {line}: pairwise r = " + ", ".join(f"{v:.2f}" for v in w.pearson_r))

short = [f"{l[0]}{r[-1]}" for l, r in labels]
excl_idx = [k for k, (l, r) in enumerate(labels) if f"{l} {r}" in EXCLUDED]


def heat(ax, M, title, fmt="{:.2f}"):
    im = ax.imshow(M, cmap=SEQ, vmin=0, vmax=1, origin="upper")
    for i in range(n):
        for j in range(n):
            if i == j:
                ax.add_patch(Rectangle((j - 0.5, i - 0.5), 1, 1, facecolor="#eceae5",
                                       edgecolor=SURF, zorder=3))
                continue
            ax.text(j, i, fmt.format(M[i, j]), ha="center", va="center", fontsize=6.6,
                    color="white" if M[i, j] > 0.62 else INK2, zorder=4)
    for k in range(0, n, 3):                       # cell-line blocks
        ax.add_patch(Rectangle((k - 0.5, k - 0.5), 3, 3, fill=False, edgecolor=INK,
                               linewidth=1.8, zorder=5))
    for k in excl_idx:                             # the excluded plate
        ax.add_patch(Rectangle((k - 0.5, -0.5), 1, n, facecolor=BAD, alpha=0.10,
                               edgecolor="none", zorder=2))
        ax.add_patch(Rectangle((-0.5, k - 0.5), n, 1, facecolor=BAD, alpha=0.10,
                               edgecolor="none", zorder=2))
    ax.set_xticks(range(n)); ax.set_xticklabels(short, fontsize=7.6)
    ax.set_yticks(range(n)); ax.set_yticklabels(short, fontsize=7.6)
    for t, k in zip(ax.get_xticklabels(), range(n)):
        if k in excl_idx:
            t.set_color(BAD); t.set_fontweight("bold")
    for t, k in zip(ax.get_yticklabels(), range(n)):
        if k in excl_idx:
            t.set_color(BAD); t.set_fontweight("bold")
    ax.set_title(title, pad=10)
    ax.tick_params(length=0)
    for sp in ax.spines.values():
        sp.set_visible(False)
    for cx, nm in zip((1, 4, 7), LINES):
        ax.text(cx, -0.075, nm, ha="center", va="top", fontsize=7.6, color=INK2,
                transform=ax.get_xaxis_transform())
        ax.text(-0.085, cx, nm, ha="right", va="center", fontsize=7.6, color=INK2,
                rotation=90, transform=ax.get_yaxis_transform())
    return im


fig = plt.figure(figsize=(13.8, 5.2))
gs = fig.add_gridspec(1, 3, width_ratios=[1.0, 1.0, 0.72], wspace=0.30)

ax = fig.add_subplot(gs[0, 0]); tag(ax, "a")
im = heat(ax, P, "Pearson r — do the large effects reproduce?")
ax.set_box_aspect(1.0)
ax.set_xlabel("H = HCT116, R = RKO, S = SW480; the number is the plate",
              fontsize=7.3, color=INK3, labelpad=24)

ax = fig.add_subplot(gs[0, 1]); tag(ax, "b", dx=-0.12)
heat(ax, S, "Spearman ρ — does the whole panel reproduce?")
ax.set_box_aspect(1.0)

cb = fig.colorbar(im, ax=fig.axes[:2], fraction=0.028, pad=0.045, shrink=0.55)
cb.set_label("correlation between plates", fontsize=7.6)
cb.ax.tick_params(labelsize=7); cb.outline.set_visible(False)

# (c) the summary that drives the decision
ax = fig.add_subplot(gs[0, 2]); tag(ax, "c", dx=-0.22)
mean_within = []
for k, (l, r) in enumerate(labels):
    sib = [m for m, (l2, r2) in enumerate(labels) if l2 == l and m != k]
    mean_within.append(np.mean([P[k, m] for m in sib]))
y = np.arange(n)[::-1]
cols_ = [BAD if k in excl_idx else "#4a86c8" for k in range(n)]
ax.barh(y, mean_within, height=0.66, color=cols_)
ax.axvline(0.5, color=INK3, lw=0.9, ls=(0, (3, 2)))
for yi, v in zip(y, mean_within):
    ax.text(v + 0.012, yi, f"{v:.2f}", va="center", fontsize=7.4, color=INK2)
ax.set_yticks(y); ax.set_yticklabels(short, fontsize=7.6)
for t, k in zip(ax.get_yticklabels(), range(n)):
    if k in excl_idx:
        t.set_color(BAD); t.set_fontweight("bold")
ax.set_xlim(0, max(mean_within) * 1.22)
ax.set_xlabel("mean r with the other plates of its own line", fontsize=7.6)
ax.set_title("Concordance per plate")
ax.set_box_aspect(1.0)
ax.grid(axis="x"); ax.set_axisbelow(True)

hct = within[within.plate_a.str.startswith("HCT116")].pearson_r
rko = within[within.plate_a.str.startswith("RKO")].pearson_r
sw = within[within.plate_a.str.startswith("SW480")].pearson_r
sp_med, pe_med = within.spearman_rho.median(), within.pearson_r.median()
if sp_med < pe_med - 0.15:
    spear = (f"Spearman sits well below Pearson (within-line median {sp_med:.2f} against "
             f"{pe_med:.2f}), so the agreement is carried by the large effects while the "
             "bulk of the panel,\nwhich barely moves, does not reproduce.")
else:
    spear = (f"Spearman is as high as Pearson (within-line median {sp_med:.2f} against "
             f"{pe_med:.2f}), so the agreement is not carried by a handful of large effects "
             "alone \u2014 the ranking of the\nwhole panel reproduces too.")
qcnote = ("" if not CFG["qc_fail_line"] else
          f"\n{CFG['qc_fail_line']} nonetheless fails the plate-level Z' check for this assay "
          "(\u00a72.15), which correlations do not reveal: a line can be reproducible and still "
          "unable to tell a killed well from a living one.")
fig.text(0.012, -0.02,
         "Read the three outlined blocks on the diagonal: those are the replicates of one cell line, and their "
         "agreement is what a replicate is for.\n"
         f"HCT116 r {hct.min():.2f}\u2013{hct.max():.2f}, RKO r {rko.min():.2f}\u2013{rko.max():.2f}, "
         f"SW480 r {sw.min():.2f}\u2013{sw.max():.2f}; within-line median {pe_med:.2f} against "
         f"{across.pearson_r.median():.2f} across lines. " + spear + qcnote,
         fontsize=7.6, color=INK2, va="top")

fig.savefig(f"figures/{CFG['fig']}.png", dpi=220, bbox_inches="tight")
fig.savefig(f"figures/{CFG['fig']}.pdf", bbox_inches="tight")
plt.close(fig)
print(f"\n{CFG['fig']} done")
