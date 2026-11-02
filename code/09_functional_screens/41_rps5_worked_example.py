#!/usr/bin/env python3
"""
41 - RPS5 as a worked example: why the heatmap cell can be dark blue while the paired
t-test says "not significant".

The two things are different quantities. The COLOUR encodes the effect - the mean of the
three within-plate differences. The OUTLINE encodes the p-value, which is that mean
divided by the scatter of the three differences. A large mean with scattered replicates
is a dark cell without an outline; that is not an inconsistency, it is the figure
reporting both of the things a reader needs.

With three pairs the test has 2 degrees of freedom, and the two-sided 5 % critical value
of t on 2 df is 4.303 - not the familiar 1.96. The effect has to be more than 4.3 times
its own standard error before the test will call it, which is a demanding requirement
that a three-point experiment often fails even for a real, large effect.

Outputs
  results/rps5_input_data.csv
  figures/fig_rps5_worked_example.{png,pdf}
"""
import pandas as pd, numpy as np, os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm
from scipy.stats import ttest_rel, t as tdist

os.makedirs("figures", exist_ok=True)
GENE = "RPS5"
LINES = ["HCT116", "RKO", "SW480"]
REPS = ["rep1", "rep2", "rep3"]
INK, INK2, INK3 = "#0b0b0b", "#52514e", "#8a8880"
SURF = "#fcfcfb"
NTC, GC = "#8a8880", "#1f5fae"
DIV = LinearSegmentedColormap.from_list("div", ["#184f95", "#6da7ec", "#f2f1ed",
                                                "#f0a09f", "#c3202f"])
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


def tag(ax, s, dx=-0.13, dy=1.08):
    ax.text(dx, dy, s, transform=ax.transAxes, fontsize=11, fontweight="bold",
            color=INK, va="top", ha="left")


C = pd.read_csv("results/func_crispr_long.csv")
rows, stats = [], {}
for line in LINES:
    d = C[C.cell_line == line]
    for r in REPS:
        nt = d[(d.kind == "negative") & (d.replicate == r)]["value"].to_numpy(float)
        gv = float(d[(d.kind == "sample") & (d.gene == GENE) & (d.replicate == r)]["value"].iloc[0])
        rows.append({"cell_line": line, "replicate": r, "gene": GENE, "gene_value": gv,
                     "n_NT_wells": len(nt), "NT_mean": nt.mean(), "NT_sd": nt.std(ddof=1),
                     "NT_min": nt.min(), "NT_max": nt.max(),
                     "difference": gv - nt.mean(),
                     **{f"NT_well_{i+1}": v for i, v in enumerate(np.sort(nt))}})
    sub = [x for x in rows if x["cell_line"] == line]
    g = np.array([x["gene_value"] for x in sub])
    c = np.array([x["NT_mean"] for x in sub])
    diff = g - c
    tt = ttest_rel(g, c)
    stats[line] = {"diff": diff, "mean": diff.mean(), "sd": diff.std(ddof=1),
                   "se": diff.std(ddof=1) / np.sqrt(3), "t": float(tt.statistic),
                   "p": float(tt.pvalue),
                   "sd_needed": abs(diff.mean()) * np.sqrt(3) / tdist.ppf(0.975, 2)}

D = pd.DataFrame(rows)
D.to_csv("results/rps5_input_data.csv", index=False)
TCRIT = float(tdist.ppf(0.975, 2))
print(D.round(2).to_string(index=False))
print(f"\ncritical |t| on 2 df = {TCRIT:.3f}")
for line in LINES:
    s = stats[line]
    print(f"{line}: differences {np.round(s['diff'],2)} | mean {s['mean']:+.2f} "
          f"sd {s['sd']:.2f} SE {s['se']:.2f} | t {s['t']:+.2f} P {s['p']:.4f} "
          f"| would need sd <= {s['sd_needed']:.2f} to reach P < 0.05")

# ------------------------------------------------------------------ figure
fig = plt.figure(figsize=(12.6, 7.4))
gs = fig.add_gridspec(2, 3, height_ratios=[1.25, 1.0], hspace=0.52, wspace=0.30)

# (a-c) the paired input data, one panel per cell line
for j, line in enumerate(LINES):
    ax = fig.add_subplot(gs[0, j])
    tag(ax, "abc"[j], dx=-0.18 if j == 0 else -0.12)
    sub = D[D.cell_line == line]
    for i, (_, r) in enumerate(sub.iterrows()):
        wells = r[[f"NT_well_{k+1}" for k in range(8)]].to_numpy(float)
        x = i + np.random.default_rng(41 + i).normal(0, 0.035, len(wells))
        ax.scatter(x - 0.16, wells, s=15, color=NTC, alpha=0.65, edgecolor="none", zorder=3)
        ax.plot([i - 0.30, i - 0.02], [r.NT_mean] * 2, color=INK, lw=2.0, zorder=4)
        ax.scatter([i + 0.20], [r.gene_value], s=62, color=GC, edgecolor="white",
                   linewidth=0.7, zorder=5, marker="s")
        ax.annotate("", xy=(i + 0.20, r.gene_value), xytext=(i + 0.20, r.NT_mean),
                    arrowprops=dict(arrowstyle="->", color=GC, lw=1.1))
        ax.text(i + 0.30, (r.gene_value + r.NT_mean) / 2, f"{r.difference:+.1f}",
                fontsize=7.4, color=GC, va="center")
    ax.set_xticks(range(3)); ax.set_xticklabels(["plate 1", "plate 2", "plate 3"], fontsize=8.2)
    ax.set_ylim(40, 125)
    if j == 0:
        ax.set_ylabel("Viability (100 = untreated)")
    ax.set_title(line)
    ax.grid(axis="y"); ax.set_axisbelow(True)
    s = stats[line]
    ax.text(0.02, 0.055, f"mean Δ {s['mean']:+.1f}, sd {s['sd']:.1f}\n"
                         f"t = {s['t']:+.2f}, P = {s['p']:.3f}",
            transform=ax.transAxes, fontsize=7.6,
            color=GC if s["p"] < 0.05 else "#b3472f", va="bottom")

fig.legend(handles=[
    Line2D([], [], marker="o", ls="", color=NTC, ms=6, label="the 8 non-targeting wells of that plate"),
    Line2D([], [], marker="_", ls="", color=INK, ms=14, markeredgewidth=2.5,
           label="their mean — the plate's control value"),
    Line2D([], [], marker="s", ls="", color=GC, ms=7, label=f"{GENE} on that plate")],
    loc="upper center", bbox_to_anchor=(0.5, 0.505), ncol=3, fontsize=7.8)

# (d) the three differences, with the confidence interval that decides significance
ax = fig.add_subplot(gs[1, :2]); tag(ax, "d", dx=-0.075)
ax.axhline(0, color=INK2, lw=1.0, zorder=2)
for j, line in enumerate(LINES):
    s = stats[line]
    sig = s["p"] < 0.05
    col = GC if sig else "#b3472f"
    ax.scatter(np.full(3, j) + np.array([-0.08, 0, 0.08]), s["diff"], s=34, color=col,
               alpha=0.8, edgecolor="none", zorder=4)
    lo, hi = s["mean"] - TCRIT * s["se"], s["mean"] + TCRIT * s["se"]
    ax.plot([j, j], [lo, hi], color=col, lw=2.6, solid_capstyle="round", zorder=3)
    ax.plot([j - 0.13, j + 0.13], [s["mean"]] * 2, color=INK, lw=1.8, zorder=5)
    ax.text(j + 0.22, s["mean"], f"{s['mean']:+.1f}", fontsize=8, color=INK2, va="center")
    ax.text(j, 17.5, "P = %.3f" % s["p"], ha="center", fontsize=8,
            color=col, fontweight="bold" if sig else "normal")
    ax.text(j, 10.5, "interval crosses zero" if not sig else "interval clear of zero",
            ha="center", fontsize=7.3, color=col)
ax.set_xticks(range(3)); ax.set_xticklabels(LINES, fontsize=8.8)
ax.set_xlim(-0.5, 2.6); ax.set_ylim(-70, 24)
ax.set_ylabel("Δ viability vs the same plate's\ncontrol mean (points)")
ax.set_title("The three paired differences, with the 95 % interval that decides the call")
ax.grid(axis="y"); ax.set_axisbelow(True)
ax.text(0.0, -0.30,
        "The bar is mean ± 4.303 × SE. 4.303, not 1.96, is the two-sided 5 % critical value of t on "
        "2 degrees of freedom — all\nthree replicates buy you. In HCT116 the drop is large (−22.7) but "
        "the three plates disagree (−12.2, −35.8, −20.1),\nso the interval reaches past zero and the "
        "gene is not called. It would have needed a replicate sd of "
        f"{stats['HCT116']['sd_needed']:.1f} or less;\nit had {stats['HCT116']['sd']:.1f}. In SW480 a "
        "smaller drop is called easily because its three plates agree to within 2 points.",
        transform=ax.transAxes, fontsize=7.5, color=INK2, va="top")

# (e) what the heatmap cell shows
ax = fig.add_subplot(gs[1, 2]); tag(ax, "e", dx=-0.30)
vals = np.array([[stats[l]["mean"] for l in LINES]]).T
lim = 33.0
ax.imshow(vals, cmap=DIV, norm=TwoSlopeNorm(vmin=-lim, vcenter=0, vmax=lim),
          aspect="auto", origin="lower", extent=(-0.5, 0.5, -0.5, 2.5))
for i, line in enumerate(LINES):
    s = stats[line]
    if s["p"] < 0.05:
        ax.add_patch(plt.Rectangle((-0.5, i - 0.5), 1, 1, fill=False, edgecolor=INK,
                                   linewidth=1.8, zorder=5))
        ax.text(0, i, f"{s['mean']:+.0f}", ha="center", va="center", fontsize=8,
                color="white", zorder=6)
ax.set_yticks(range(3)); ax.set_yticklabels(LINES, fontsize=8.4)
ax.set_xticks([]); ax.set_xlim(-0.5, 2.25)
ax.set_title(f"{GENE} as Figure 20a draws it")
for sp in ax.spines.values():
    sp.set_visible(False)
ax.tick_params(length=0)
ax.text(0.44, 0.5, "colour = the mean difference\n(how big the effect is)\n\n"
                   "outline = the paired t-test\n(how reproducible it is)\n\n"
                   "HCT116 is dark because the\ndrop is large, and unoutlined\n"
                   "because the three plates\ndisagree about its size",
        transform=ax.transAxes, fontsize=7.5, color=INK2, va="center")

fig.savefig("figures/fig_rps5_worked_example.png", dpi=220, bbox_inches="tight")
fig.savefig("figures/fig_rps5_worked_example.pdf", bbox_inches="tight")
plt.close(fig)
print("\nwrote figures/fig_rps5_worked_example.png and results/rps5_input_data.csv")
