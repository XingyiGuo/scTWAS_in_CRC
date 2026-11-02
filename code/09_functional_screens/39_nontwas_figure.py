#!/usr/bin/env python3
"""
39 - Figure 22: the same per-cell-line differential analysis as Figure 20, but for the
80 genes on the CRISPR panel that are NOT CRC risk genes from this study.

These genes are the experiment's own comparison group, and they are what makes the risk
gene result interpretable. Figure 20 shows how often a risk gene changes viability;
without this figure there is no way to know whether that rate is high, low or ordinary.
The analysis is identical to Figure 20 - a paired t-test, each replicate plate's eight
non-targeting wells averaged into that plate's control value, the gene paired against its
own plate, three pairs per cell line, nominal P < 0.05. The effect is the mean of the
three within-plate differences in viability points.

Two things replace the TWAS strips of Figure 20, which cannot exist here. On the left,
a mark for the panel's core-machinery genes (ribosomal, spliceosomal, replication and
polymerase subunits), because those are the assay's built-in positive controls and any
working screen must find them. On the right, the number of cell lines in which the gene
changed - the same reproducibility axis Figure 20 uses, with no direction to check it
against.

Outputs
  figures/fig22_nontwas_by_cellline.{png,pdf}
  results/nontwas_by_cellline.csv
  results/nontwas_comparison.json
"""
import pandas as pd, numpy as np, json, os, re
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch, Rectangle
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm
from scipy.stats import fisher_exact, mannwhitneyu

os.makedirs("figures", exist_ok=True)
INK, INK2, INK3 = "#0b0b0b", "#52514e", "#8a8880"
SURF = "#fcfcfb"
DIV = LinearSegmentedColormap.from_list("div", ["#184f95", "#6da7ec", "#f2f1ed",
                                                "#f0a09f", "#c3202f"])
RISK, OTHER, CORE = "#8c2d3f", "#7d99b8", "#4a3aa7"
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


LINES = ["HCT116", "RKO", "SW480"]
ESS = re.compile(r"^(RPS|RPL|PSM|SF3|PRPF|POLR|EIF|SNRP|MCM|TBCE)")
L = pd.read_csv("results/crispr_by_cellline.csv")
L["core"] = L["gene"].str.match(ESS)

O = L[~L.is_risk].copy()
Rk = L[L.is_risk].copy()


def widen(df):
    w = df.pivot_table(index="gene", columns="cell_line",
                       values=["delta_paired", "p_paired", "nominal"])
    f = pd.DataFrame(index=w.index)
    for line in LINES:
        f[f"delta_{line}"] = w[("delta_paired", line)]
        f[f"pT_{line}"] = w[("p_paired", line)]
        f[f"nom_{line}"] = w[("nominal", line)].astype(bool)
    f["n_lines_sig"] = f[[f"nom_{l}" for l in LINES]].sum(axis=1)
    f["sig_any"] = f["n_lines_sig"] > 0
    f["mean_delta"] = f[[f"delta_{l}" for l in LINES]].mean(axis=1)
    f["core"] = [bool(ESS.match(g)) for g in f.index]
    return f.reset_index()


W = widen(O)
RW = widen(Rk)
W.sort_values(["n_lines_sig", "mean_delta"], ascending=[True, False],
              inplace=True, ignore_index=True)
N = len(W)
W.to_csv("results/nontwas_by_cellline.csv", index=False)

# ---------------------------------------------------------------- comparison statistics
n_o, n_r = len(W), len(RW)
sig_o, sig_r = int(W.sig_any.sum()), int(RW.sig_any.sum())
orr, pf = fisher_exact([[sig_r, n_r - sig_r], [sig_o, n_o - sig_o]])
u = mannwhitneyu(RW["mean_delta"], W["mean_delta"], alternative="two-sided")
per_line = {l: {"risk": int(Rk[Rk.cell_line == l]["nominal"].sum()),
                "other": int(O[O.cell_line == l]["nominal"].sum())} for l in LINES}
core = W[W.core]
J = {
    "n_other": n_o, "n_risk": n_r,
    "sig_any_other": sig_o, "sig_any_risk": sig_r,
    "frac_other": sig_o / n_o, "frac_risk": sig_r / n_r,
    "fisher_odds_ratio": float(orr), "fisher_p": float(pf),
    "median_delta_other": float(W["mean_delta"].median()),
    "median_delta_risk": float(RW["mean_delta"].median()),
    "mannwhitney_p": float(u.pvalue),
    "per_line": per_line,
    "per_line_risk_n": n_r, "per_line_other_n": n_o,
    "n_by_lines_other": {int(k): int(v) for k, v in
                         W["n_lines_sig"].value_counts().sort_index().items()},
    "n_by_lines_risk": {int(k): int(v) for k, v in
                        RW["n_lines_sig"].value_counts().sort_index().items()},
    "core_genes": core["gene"].tolist(),
    "n_core": int(len(core)), "n_core_sig": int(core.sig_any.sum()),
    "core_median_delta": float(core["mean_delta"].median()),
    "noncore_median_delta": float(W[~W.core]["mean_delta"].median()),
    "core_p": float(mannwhitneyu(core["mean_delta"], W[~W.core]["mean_delta"],
                                 alternative="less").pvalue),
    "strongest": W.nsmallest(8, "mean_delta")[["gene", "mean_delta", "n_lines_sig"]]
                  .to_dict("records"),
}
json.dump(J, open("results/nontwas_comparison.json", "w"), indent=1, default=float)
print(json.dumps({k: v for k, v in J.items() if k not in ("strongest", "core_genes")},
                 indent=1, default=float))

def draw_matrix(figure, ax):
    """The Figure 20a-style matrix. Used twice: as panel (a) of Figure 22 and as
    a standalone figure, so the two can never drift apart."""
    Dm = W[[f"delta_{l}" for l in LINES]].to_numpy(float)
    lim = max(float(np.nanpercentile(np.abs(Dm), 97)), 5.0)
    norm = TwoSlopeNorm(vmin=-lim, vcenter=0, vmax=lim)
    im = ax.imshow(Dm, cmap=DIV, norm=norm, aspect="auto", origin="lower",
                   extent=(-0.5, 2.5, -0.5, N - 0.5))
    for i, c in enumerate(W["core"]):
        ax.add_patch(Rectangle((-1.15, i - 0.5), 0.5, 1,
                               facecolor=CORE if c else "#e4e2dd",
                               edgecolor=SURF, linewidth=0.5, zorder=2))
    for i, k in enumerate(W["n_lines_sig"]):
        sh = {0: "#e4e2dd", 1: "#bcd2ee", 2: "#6d9fd6", 3: "#1f5fae"}[int(k)]
        ax.add_patch(Rectangle((2.66, i - 0.5), 0.5, 1, facecolor=sh,
                               edgecolor=SURF, linewidth=0.5, zorder=2))
    ax.text(-0.9, N - 0.15, "core", ha="center", va="bottom", fontsize=7.4,
            color=INK2, fontweight="bold")
    ax.text(2.91, N - 0.15, "lines", ha="center", va="bottom", fontsize=7.4,
            color=INK2, fontweight="bold")
    for j, line in enumerate(LINES):
        for i in np.flatnonzero(W[f"nom_{line}"].to_numpy(bool)):
            ax.add_patch(Rectangle((j - 0.5, i - 0.5), 1, 1, fill=False,
                                   edgecolor=INK, linewidth=1.3, zorder=5))
            v = W[f"delta_{line}"].iloc[i]
            ax.text(j, i, f"{v:+.0f}", ha="center", va="center", fontsize=5.7,
                    color="white" if abs(v) > lim * 0.55 else INK, zorder=6)
    ax.set_xticks(range(3)); ax.set_xticklabels(LINES, fontsize=8.8)
    ax.set_yticks(range(N))
    ax.set_yticklabels(W["gene"], fontsize=6.3)
    for t, c in zip(ax.get_yticklabels(), W["core"]):
        t.set_color(CORE if c else INK2)
    ax.set_xlim(-1.25, 3.16); ax.set_ylim(-0.5, N - 0.5)
    ax.set_title(f"The {N} CRISPR panel genes that are not CRC risk genes", pad=20)
    ax.set_xlabel("Outlined cells are significant at P < 0.05 (Welch t-test, 3 replicates vs "
                  "controls)", fontsize=7.2, color=INK3, style="italic", labelpad=8)
    ax.tick_params(length=0)
    for sp in ax.spines.values():
        sp.set_visible(False)
    cb = figure.colorbar(im, ax=ax, fraction=0.038, pad=0.055, shrink=0.26, anchor=(0.0, 1.0))
    cb.set_label("Δ viability vs same-plate control (points)", fontsize=7.4)
    cb.ax.tick_params(labelsize=7); cb.outline.set_visible(False)
    hl = [Patch(facecolor=CORE, label=f"core machinery ({J['n_core']}) — positive controls"),
          Patch(facecolor="#e4e2dd", label="other panel gene"),
          Line2D([], [], marker="s", ls="", mfc="none", mec=INK, ms=9, markeredgewidth=1.4,
                 label="paired t-test P < 0.05")]
    for k, col, lab in [(1, "#bcd2ee", "significant in 1 line"),
                        (2, "#6d9fd6", "2 lines"), (3, "#1f5fae", "all 3 lines")]:
        if J["n_by_lines_other"].get(k, 0):
            hl.append(Patch(facecolor=col, label=lab))
    ax.legend(handles=hl, loc="upper left", bbox_to_anchor=(-0.34, -0.052), ncol=2,
              fontsize=7.2, columnspacing=1.0)



# ---------------------------------------------------------------- figure
fig = plt.figure(figsize=(14.2, 14.6))
gs = fig.add_gridspec(3, 2, width_ratios=[1.06, 1.0], height_ratios=[1, 1, 1],
                      wspace=0.46, hspace=0.80)

# (a) matrix -------------------------------------------------------------------
ax = fig.add_subplot(gs[:, 0]); tag(ax, "a", dx=-0.30, dy=1.016)
draw_matrix(fig, ax)

# (b) per line -----------------------------------------------------------------
ax = fig.add_subplot(gs[0, 1]); tag(ax, "b")
x = np.arange(3, dtype=float)
fr_o = [100 * per_line[l]["other"] / n_o for l in LINES]
fr_r = [100 * per_line[l]["risk"] / n_r for l in LINES]
ax.bar(x - 0.19, fr_r, width=0.36, color=RISK, label=f"CRC risk genes (n = {n_r})")
ax.bar(x + 0.19, fr_o, width=0.36, color=OTHER, label=f"other panel genes (n = {n_o})")
for xi, v, c in zip(x - 0.19, fr_r, [per_line[l]["risk"] for l in LINES]):
    ax.text(xi, v + 0.8, f"{c}", ha="center", fontsize=7.6, color=RISK)
for xi, v, c in zip(x + 0.19, fr_o, [per_line[l]["other"] for l in LINES]):
    ax.text(xi, v + 0.8, f"{c}", ha="center", fontsize=7.6, color="#4a6b8f")
ax.set_xticks(x); ax.set_xticklabels(LINES, fontsize=8.8)
ax.set_ylabel("% of genes with P < 0.05")
ax.set_ylim(0, max(fr_o + fr_r) * 1.30)
ax.set_title("Significant changes per cell line, as a share of each group")
ax.legend(loc="upper left", fontsize=7.2)
ax.grid(axis="y"); ax.set_axisbelow(True)
ax.text(0.0, -0.30, "Counts are printed above the bars; the bars are percentages so the "
        "two groups are comparable despite\ntheir different sizes.",
        transform=ax.transAxes, fontsize=7.4, color=INK2, va="top")

# (c) reproducibility ----------------------------------------------------------
ax = fig.add_subplot(gs[1, 1]); tag(ax, "c")
ks = [0, 1, 2, 3]
b_o = [100 * J["n_by_lines_other"].get(k, 0) / n_o for k in ks]
b_r = [100 * J["n_by_lines_risk"].get(k, 0) / n_r for k in ks]
kk = np.array(ks, float)
ax.bar(kk - 0.19, b_r, width=0.36, color=RISK, label="CRC risk genes")
ax.bar(kk + 0.19, b_o, width=0.36, color=OTHER, label="other panel genes")
for xi, v, c in zip(kk - 0.19, b_r, [J["n_by_lines_risk"].get(k, 0) for k in ks]):
    if c:
        ax.text(xi, v + 1.2, f"{c}", ha="center", fontsize=7.4, color=RISK)
for xi, v, c in zip(kk + 0.19, b_o, [J["n_by_lines_other"].get(k, 0) for k in ks]):
    if c:
        ax.text(xi, v + 1.2, f"{c}", ha="center", fontsize=7.4, color="#4a6b8f")
ax.set_xticks(ks); ax.set_xticklabels(["no line", "1 line", "2 lines", "all 3"], fontsize=8.4)
ax.set_ylabel("% of group")
ax.set_ylim(0, max(b_o + b_r) * 1.26)
ax.set_title("How reproducible the changes are")
ax.legend(loc="upper right", fontsize=7.2)
ax.grid(axis="y"); ax.set_axisbelow(True)
ax.text(0.0, -0.30,
        "Significant in at least one line: {}/{} of the risk genes ({:.0f} %) against {}/{} "
        "of the others ({:.0f} %) —\nodds ratio {:.2f}, Fisher P = {:.2f}. The risk genes are "
        "not more often affected than the genes assayed\nalongside them, which is the "
        "comparison that makes §2.15's negative result meaningful.".format(
            sig_r, n_r, 100 * J["frac_risk"], sig_o, n_o, 100 * J["frac_other"],
            J["fisher_odds_ratio"], J["fisher_p"]),
        transform=ax.transAxes, fontsize=7.4, color=INK2, va="top")

# (d) effect distributions -----------------------------------------------------
ax = fig.add_subplot(gs[2, 1]); tag(ax, "d")
groups = [("CRC risk\ngenes", RW["mean_delta"].to_numpy(float), RISK),
          ("other panel\ngenes", W[~W.core]["mean_delta"].to_numpy(float), OTHER),
          ("core machinery\n(controls)", W[W.core]["mean_delta"].to_numpy(float), CORE)]
for i, (lab, v, col) in enumerate(groups):
    parts = ax.violinplot([v], positions=[i], widths=0.74, showextrema=False)
    for b in parts["bodies"]:
        b.set_facecolor(col); b.set_alpha(0.28); b.set_edgecolor("none")
    ax.scatter(np.full(len(v), i) + np.random.default_rng(39 + i).normal(0, 0.055, len(v)),
               v, s=11, color=col, alpha=0.75, edgecolor="none", zorder=3)
    ax.plot([i - 0.26, i + 0.26], [np.median(v)] * 2, color=INK, lw=1.7, zorder=4)
    ax.text(i, np.median(v), f"{np.median(v):+.1f}", ha="center", va="center",
            fontsize=7.2, color=INK, zorder=6,
            bbox=dict(facecolor=SURF, edgecolor="none", pad=1.4))
ax.axhline(0, color="#d8d7d2", lw=0.9, zorder=1)
ax.set_xticks(range(3)); ax.set_xticklabels([g[0] for g in groups], fontsize=8.2)
ax.set_ylabel("Mean paired Δ viability across\nthe three lines (points)")
ax.set_title("Effect sizes: the controls separate, the two gene sets do not")
ax.grid(axis="y"); ax.set_axisbelow(True)
ax.text(0.0, -0.36,
        "Risk genes against the other panel genes: median {:+.1f} vs {:+.1f} percentage "
        "points, Mann-Whitney P = {:.2f}.\nThe core-machinery genes, which any working "
        "screen must find, sit at {:+.1f} (P = {:.1e} against the rest of\nthe non-risk "
        "panel) — the assay separates what it should, and the risk genes are simply not "
        "separated from\ntheir comparison group.".format(
            J["median_delta_risk"], J["noncore_median_delta"], J["mannwhitney_p"],
            J["core_median_delta"], J["core_p"]),
        transform=ax.transAxes, fontsize=7.4, color=INK2, va="top")

fig.savefig("figures/fig22_nontwas_by_cellline.png", dpi=220, bbox_inches="tight")
fig.savefig("figures/fig22_nontwas_by_cellline.pdf", bbox_inches="tight")
plt.close(fig)
print("fig22 done")

# ---------------------------------------------------------------- standalone matrix
# The direct counterpart of Figure 20a on its own, for use without the summary panels.
h = max(6.0, 0.135 * N + 2.6)
f2 = plt.figure(figsize=(6.6, h))
ax = f2.add_axes([0.20, 0.10, 0.52, 0.84])
draw_matrix(f2, ax)
f2.savefig("figures/fig22a_nontwas_matrix.png", dpi=220, bbox_inches="tight")
f2.savefig("figures/fig22a_nontwas_matrix.pdf", bbox_inches="tight")
plt.close(f2)
print("fig22a (standalone matrix) done")
