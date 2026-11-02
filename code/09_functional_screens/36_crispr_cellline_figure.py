#!/usr/bin/env python3
"""
36 - Figure 20: the CRISPR-Cas9 screen cell line by cell line, for the 36 S1 risk genes,
reported as a plain differential analysis against the controls.

One test only, at the user's request: a PAIRED t-test. Each replicate plate carries eight
non-targeting wells; their mean is that plate's control value, so a gene has one
(gene, control) pair per plate and the three pairs are tested with a paired t-test at
nominal P < 0.05. The effect shown is the mean of the three within-plate differences, in
viability points. Pairing inside the plate removes the plate-to-plate offset without any
normalisation step. The robust-z / empirical-null / FDR analysis
is no longer drawn; it stays in Tables S46 and S47 for anyone who wants it.

The figure also answers the question that motivates the screen. Each gene carries a
germline TWAS direction, and that direction predicts what should happen when the gene is
removed: a risk-increasing gene should be one the tumour cell needs (viability falls), a
protective gene should be one it is better off without (viability rises). The strip on
the right records, per gene, whether the measured change agrees.
"""
import pandas as pd, numpy as np, json, os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch, Rectangle
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm

os.makedirs("figures", exist_ok=True)
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


LINES = ["HCT116", "RKO", "SW480"]
W = pd.read_csv("results/crispr_risk36_by_cellline.csv")
J = json.load(open("results/crispr_by_cellline.json"))
T = J["ttest"]
RANK = {NS: 0, "conflicting": 1, OPPOSES: 2, SUPPORTS: 3}
W["rank"] = W["twas_support"].map(RANK).fillna(0)
W = W.sort_values(["rank", "n_lines_nominal", "mean_delta_all"],
                  ascending=[True, True, False]).reset_index(drop=True)
N = len(W)

fig = plt.figure(figsize=(13.8, 9.8))
gs = fig.add_gridspec(2, 2, width_ratios=[1.32, 1.0], height_ratios=[1.0, 0.92],
                      wspace=0.44, hspace=0.62)

# ---------------------------------------------------------------- (a) matrix
ax = fig.add_subplot(gs[:, 0]); tag(ax, "a", dx=-0.33, dy=1.035)
Dm = W[[f"delta_{l}" for l in LINES]].to_numpy(float)
lim = max(float(np.nanpercentile(np.abs(Dm), 97)), 5.0)
norm = TwoSlopeNorm(vmin=-lim, vcenter=0, vmax=lim)
im = ax.imshow(Dm, cmap=DIV, norm=norm, aspect="auto", origin="lower",
               extent=(-0.5, 2.5, -0.5, N - 0.5))

for i, d in enumerate(W["twas_direction"].fillna("")):
    ax.add_patch(Rectangle((-1.15, i - 0.5), 0.5, 1, facecolor=TW.get(d, "#d8d7d2"),
                           edgecolor=SURF, linewidth=0.6, zorder=2))
for i, v in enumerate(W["twas_support"].fillna(NS)):
    ax.add_patch(Rectangle((2.66, i - 0.5), 0.5, 1, facecolor=VC.get(v, VC[NS]),
                           edgecolor=SURF, linewidth=0.6, zorder=2))
ax.text(-0.9, N - 0.25, "TWAS", ha="center", va="bottom", fontsize=7.4,
        color=INK2, fontweight="bold")
ax.text(2.91, N - 0.25, "verdict", ha="center", va="bottom", fontsize=7.4,
        color=INK2, fontweight="bold")

for j, line in enumerate(LINES):
    nom = W[f"nom_{line}"].to_numpy(bool)
    for i in np.flatnonzero(nom):
        ax.add_patch(Rectangle((j - 0.5, i - 0.5), 1, 1, fill=False,
                               edgecolor=INK, linewidth=1.6, zorder=5))
        v = W[f"delta_{line}"].iloc[i]
        ax.text(j, i, f"{v:+.0f}", ha="center", va="center", fontsize=6.3,
                color="white" if abs(v) > lim * 0.55 else INK, zorder=6)

NPLATE = {r["cell_line"]: int(r["n_rep"]) for _, r in
          pd.read_csv("results/crispr_by_cellline.csv")
            .groupby("cell_line", as_index=False)["n_rep"].max().iterrows()}
ax.set_xticks(range(3))
ax.set_xticklabels([f"{l}\n({NPLATE.get(l, 3)} plates, df = {NPLATE.get(l, 3) - 1})"
                    for l in LINES], fontsize=8.2)
for t, l in zip(ax.get_xticklabels(), LINES):
    if NPLATE.get(l, 3) < 3:
        t.set_color("#b3472f")
ax.set_yticks(range(N))
labs = [f"{ {'oncogene-like': '▲', 'tumour-suppressor-like': '▼'}.get(c, ' ') } {g}"
        for g, c in zip(W["gene"], W["class"])]
ax.set_yticklabels(labs, fontsize=7.4)
for t, cls in zip(ax.get_yticklabels(), W["class"]):
    t.set_color(POS if cls == "oncogene-like"
                else (NEG if cls == "tumour-suppressor-like" else INK2))
ax.set_xlim(-1.25, 3.16); ax.set_ylim(-0.5, N - 0.5)
ax.set_title("CRISPR-Cas9: paired change in viability against the plate controls", pad=14)
ax.set_xlabel("Outlined cells are significant at P < 0.05 (paired t-test, 3 replicates vs "
              "their own plate controls)", fontsize=7.3, color=INK3, style="italic",
              labelpad=7)
ax.tick_params(length=0)
for sp in ax.spines.values():
    sp.set_visible(False)
cb = fig.colorbar(im, ax=ax, fraction=0.040, pad=0.055, shrink=0.40, anchor=(0.0, 1.0))
cb.set_label("Δ viability vs same-plate control (points)", fontsize=7.4)
cb.ax.tick_params(labelsize=7); cb.outline.set_visible(False)
ax.legend(handles=[
    Patch(facecolor=TW["risk-increasing"],
          label=f"TWAS risk-increasing ({J['twas_direction_counts'].get('risk-increasing', 0)})"),
    Patch(facecolor=VC[SUPPORTS], label=f"supports TWAS ({T['gene_verdict'].get(SUPPORTS, 0)})"),
    Patch(facecolor=TW["protective"],
          label=f"TWAS protective ({J['twas_direction_counts'].get('protective', 0)})"),
    Patch(facecolor=VC[OPPOSES], label=f"opposes TWAS ({T['gene_verdict'].get(OPPOSES, 0)})"),
    Line2D([], [], marker="$▲$", ls="", color=POS, ms=8, label="oncogene-like (§2.14)"),
    Patch(facecolor=VC[NS], label=f"no significant change ({T['gene_verdict'].get(NS, 0)})"),
    Line2D([], [], marker="$▼$", ls="", color=NEG, ms=8, label="tumour-suppressor-like"),
    Line2D([], [], marker="s", ls="", mfc="none", mec=INK, ms=9, markeredgewidth=1.7,
           label="paired t-test P < 0.05")],
    loc="upper left", bbox_to_anchor=(-0.37, -0.048), ncol=2, fontsize=7.2,
    columnspacing=1.0)

# ---------------------------------------------------------------- (b) per line
ax = fig.add_subplot(gs[0, 1]); tag(ax, "b")
x = np.arange(3, dtype=float)
sup = np.array([T["per_line_support"][l] for l in LINES], float)
opp = np.array([T["per_line_oppose"][l] for l in LINES], float)
ax.bar(x, sup, width=0.56, color=VC[SUPPORTS], label="change agrees with TWAS direction")
ax.bar(x, opp, width=0.56, bottom=sup, color=VC[OPPOSES], label="change opposes it")
for xi, a, b in zip(x, sup, opp):
    if a:
        ax.text(xi, a / 2, f"{int(a)}", ha="center", va="center", fontsize=8.4,
                color="white", fontweight="bold")
    if b:
        ax.text(xi, a + b / 2, f"{int(b)}", ha="center", va="center", fontsize=8.4,
                color="white", fontweight="bold")
    ax.text(xi, a + b + 0.35, f"{int(a + b)} of 36", ha="center", fontsize=7.6, color=INK2)
ax.set_xticks(x); ax.set_xticklabels(LINES, fontsize=8.8)
ax.set_ylabel("Risk genes with P < 0.05")
ax.set_ylim(0, (sup + opp).max() * 1.30)
ax.set_title("Significant changes per cell line, by agreement with TWAS")
ax.legend(loc="upper left", fontsize=7.2)
ax.grid(axis="y"); ax.set_axisbelow(True)

# ---------------------------------------------------------------- (c) verdict
ax = fig.add_subplot(gs[1, 1]); tag(ax, "c")
ks = sorted(int(k) for k in T["by_lines_verdict"])
bv = T["by_lines_verdict"]
b_sup = [bv[str(k)].get(SUPPORTS, 0) for k in ks]
b_opp = [bv[str(k)].get(OPPOSES, 0) for k in ks]
b_con = [bv[str(k)].get("conflicting", 0) for k in ks]
b_ns = [bv[str(k)].get(NS, 0) for k in ks]
kk = np.array(ks, float)
bot = np.zeros(len(ks))
for vals, col, lab in [(b_ns, VC[NS], "no significant change"),
                       (b_sup, VC[SUPPORTS], "supports TWAS"),
                       (b_opp, VC[OPPOSES], "opposes TWAS"),
                       (b_con, VC["conflicting"], "conflicting between lines")]:
    v = np.array(vals, float)
    if v.sum():
        ax.bar(kk, v, bottom=bot, width=0.6, color=col, label=lab)
        for xi, a, b in zip(kk, v, bot):
            if a:
                ax.text(xi, b + a / 2, f"{int(a)}", ha="center", va="center", fontsize=8,
                        color=INK2 if col == VC[NS] else "white", fontweight="bold")
        bot = bot + v
ax.set_xticks(ks)
ax.set_xticklabels(["no line", "1 line", "2 lines", "all 3"][:len(ks)], fontsize=8.4)
ax.set_xlabel("Cell lines with a significant change", fontsize=8)
ax.set_ylabel("Risk genes (of 36)")
ax.set_ylim(0, bot.max() * 1.22)
ax.set_title("Gene-level verdict")
ax.legend(loc="upper right", fontsize=7.2)
ax.grid(axis="y"); ax.set_axisbelow(True)
cc, cg = T["concordance_calls"], T["concordance_genes"]
lead = ("Agreement runs with the TWAS direction, though on this many calls it is a trend "
        "rather than a result." if cc["p"] >= 0.05 else
        "Agreement is better than chance.")
ax.text(0.0, -0.34,
        "{} Of the {} significant gene × cell-line changes, {} move the way the\ngermline "
        "TWAS direction predicts ({:.0f} %, binomial P = {:.3f}); at the gene level {} of {} "
        "agree ({:.0f} %, P = {:.3f}).\nA risk-increasing gene is expected to lose viability "
        "when edited and a protective gene to gain it; genes with no\nsignificant change in "
        "any line are not counted either way.".format(
            lead, cc["n"], cc["n_support"], 100 * cc["frac"], cc["p"],
            cg["n_support"], cg["n"], 100 * cg["frac"], cg["p"]),
        transform=ax.transAxes, fontsize=7.4, color=INK2, va="top")

fig.savefig("figures/fig20_crispr_by_cellline.png", dpi=220, bbox_inches="tight")
fig.savefig("figures/fig20_crispr_by_cellline.pdf", bbox_inches="tight")
plt.close(fig)
print("fig20 done")
