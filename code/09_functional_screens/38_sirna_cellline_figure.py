#!/usr/bin/env python3
"""
38 - Figure 21: the siRNA screen cell line by cell line, for all 62 S1 risk genes on the
panel, reported as a plain differential analysis against the controls.

Same design as Figure 20 and the same single test: a PAIRED t-test. Each replicate plate's
twelve non-targeting wells are averaged into that plate's control value, a gene has one
(gene, control) pair per plate, and the three pairs are tested together at nominal
P < 0.05. The effect shown is the mean of the three within-plate differences, in viability
points. In this assay the supplied values are already expressed against each plate's own
non-targeting mean, which sits at exactly 100 on every plate, so the pairing is built into
the data and the paired test is a one-sample test against 100. No robust z, no empirical
null, no FDR; those stay in Tables S48 and S49.

The right-hand strip records whether the measured change agrees with the germline TWAS
direction: risk-increasing genes are expected to lose viability when knocked down,
protective genes to gain it. The verdict is taken from HCT116 and RKO only, because
SW480's plates fail quality control - its column is drawn and shaded so the reader can
see both the calls and the reason they are not counted.
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
QC_FAIL = {"SW480"}
W = pd.read_csv("results/sirna_risk_by_cellline.csv")
J = json.load(open("results/sirna_by_cellline.json"))
T = J["ttest"]
RANK = {NS: 0, "conflicting": 1, OPPOSES: 2, SUPPORTS: 3}
W["rank"] = W["twas_support"].map(RANK).fillna(0)
W = W.sort_values(["rank", "n_qc_nominal", "mean_delta_qc"],
                  ascending=[True, True, False]).reset_index(drop=True)
N = len(W)

fig = plt.figure(figsize=(14.2, 13.4))
gs = fig.add_gridspec(3, 2, width_ratios=[1.10, 1.0], height_ratios=[1, 1, 1.05],
                      wspace=0.46, hspace=0.86)

# ---------------------------------------------------------------- (a) matrix
ax = fig.add_subplot(gs[:, 0]); tag(ax, "a", dx=-0.31, dy=1.018)
Dm = W[[f"delta_{l}" for l in LINES]].to_numpy(float)
lim = max(float(np.nanpercentile(np.abs(Dm), 97)), 5.0)
norm = TwoSlopeNorm(vmin=-lim, vcenter=0, vmax=lim)
im = ax.imshow(Dm, cmap=DIV, norm=norm, aspect="auto", origin="lower",
               extent=(-0.5, 2.5, -0.5, N - 0.5))

for i, d in enumerate(W["twas_direction"].fillna("")):
    ax.add_patch(Rectangle((-1.15, i - 0.5), 0.5, 1, facecolor=TW.get(d, "#d8d7d2"),
                           edgecolor=SURF, linewidth=0.5, zorder=2))
for i, v in enumerate(W["twas_support"].fillna(NS)):
    ax.add_patch(Rectangle((2.66, i - 0.5), 0.5, 1, facecolor=VC.get(v, VC[NS]),
                           edgecolor=SURF, linewidth=0.5, zorder=2))
ax.text(-0.9, N - 0.15, "TWAS", ha="center", va="bottom", fontsize=7.4,
        color=INK2, fontweight="bold")
ax.text(2.91, N - 0.15, "verdict", ha="center", va="bottom", fontsize=7.4,
        color=INK2, fontweight="bold")

for j, line in enumerate(LINES):
    nom = W[f"nom_{line}"].to_numpy(bool)
    for i in np.flatnonzero(nom):
        ax.add_patch(Rectangle((j - 0.5, i - 0.5), 1, 1, fill=False,
                               edgecolor=INK, linewidth=1.3, zorder=5))
        v = W[f"delta_{line}"].iloc[i]
        ax.text(j, i, f"{v:+.0f}", ha="center", va="center", fontsize=5.7,
                color="white" if abs(v) > lim * 0.55 else INK, zorder=6)
    if line in QC_FAIL:
        ax.add_patch(Rectangle((j - 0.5, -0.5), 1, N, facecolor="#6f6d68", alpha=0.16,
                               edgecolor="none", zorder=7))

ax.set_xticks(range(3))
ax.set_xticklabels([f"{l}\n(Z' = {J['z_prime_mean'][l]:+.2f})" for l in LINES], fontsize=8.4)
for t, l in zip(ax.get_xticklabels(), LINES):
    t.set_color(INK3 if l in QC_FAIL else INK)
ax.set_yticks(range(N))
labs = [f"{ {'oncogene-like': '▲', 'tumour-suppressor-like': '▼'}.get(c, ' ') } {g}"
        for g, c in zip(W["gene"], W["class"])]
ax.set_yticklabels(labs, fontsize=6.6)
for t, cls in zip(ax.get_yticklabels(), W["class"]):
    t.set_color(POS if cls == "oncogene-like"
                else (NEG if cls == "tumour-suppressor-like" else INK2))
ax.set_xlim(-1.25, 3.16); ax.set_ylim(-0.5, N - 0.5)
ax.set_title(f"siRNA: paired change in viability against the plate controls, "
             f"all {N} risk genes", pad=22)
ax.set_xlabel("Outlined cells are significant at P < 0.05 (paired t-test, 3 replicates vs "
              "their own plate controls).\nSW480 is shaded: its plates fail QC, so its calls "
              "are shown but do not enter the verdict",
              fontsize=7.2, color=INK3, style="italic", labelpad=8)
ax.tick_params(length=0)
for sp in ax.spines.values():
    sp.set_visible(False)
cb = fig.colorbar(im, ax=ax, fraction=0.038, pad=0.055, shrink=0.28, anchor=(0.0, 1.0))
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
    Line2D([], [], marker="s", ls="", mfc="none", mec=INK, ms=9, markeredgewidth=1.5,
           label="paired t-test P < 0.05")],
    loc="upper left", bbox_to_anchor=(-0.36, -0.058), ncol=2, fontsize=7.2,
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
    ax.text(xi, a + b + 0.6, f"{int(a + b)} of {N}", ha="center", fontsize=7.6, color=INK2)
ax.set_xticks(x); ax.set_xticklabels(LINES, fontsize=8.8)
ax.get_xticklabels()[2].set_color(INK3)
ax.set_ylabel(f"Risk genes with P < 0.05")
ax.set_ylim(0, (sup + opp).max() * 1.34)
ax.set_title("Significant changes per cell line, by agreement with TWAS")
ax.legend(loc="upper left", fontsize=7.2)
ax.grid(axis="y"); ax.set_axisbelow(True)
ax.text(0.30, 0.68, "SW480 fails QC, yet returns\nthe most significant changes\nof the "
        "three — and they are no\nbetter than chance",
        transform=ax.transAxes, fontsize=7.2, color=INK2, va="top")

# ---------------------------------------------------------------- (c) verdict
ax = fig.add_subplot(gs[1, 1]); tag(ax, "c")
ks = sorted(int(k) for k in T["by_lines_verdict"])
bv = T["by_lines_verdict"]
kk = np.array(ks, float)
bot = np.zeros(len(ks))
for key, col, lab in [(NS, VC[NS], "no significant change"),
                      (SUPPORTS, VC[SUPPORTS], "supports TWAS"),
                      (OPPOSES, VC[OPPOSES], "opposes TWAS"),
                      ("conflicting", VC["conflicting"], "conflicting between lines")]:
    v = np.array([bv[str(k)].get(key, 0) for k in ks], float)
    if v.sum():
        ax.bar(kk, v, bottom=bot, width=0.58, color=col, label=lab)
        for xi, a, b in zip(kk, v, bot):
            if a:
                ax.text(xi, b + a / 2, f"{int(a)}", ha="center", va="center", fontsize=8,
                        color=INK2 if col == VC[NS] else "white", fontweight="bold")
        bot = bot + v
ax.set_xticks(ks)
ax.set_xticklabels(["neither", "1 line", "both"][:len(ks)], fontsize=8.4)
ax.set_xlabel("QC-passing cell lines with a significant change", fontsize=8)
ax.set_ylabel(f"Risk genes (of {N})")
ax.set_ylim(0, bot.max() * 1.20)
ax.set_title("Gene-level verdict, HCT116 and RKO")
ax.legend(loc="upper right", fontsize=7.2)
ax.grid(axis="y"); ax.set_axisbelow(True)
cq, cg, cs = (T["concordance_calls_qc"], T["concordance_genes"],
              T["concordance_calls_sw480"])
ax.text(0.0, -0.30,
        "siRNA gives no better than chance agreement: {} of the {} significant "
        "changes in the two QC-passing\nlines go the way the TWAS direction predicts ({:.0f} %, "
        "binomial P = {:.2f}), and {} of {} genes ({:.0f} %, P = {:.2f}).\nIn SW480 it is {} of {}. "
        "The verdict strip in (a) is therefore a per-gene record, not evidence of an overall\n"
        "germline-to-dependency carry-over in this assay.".format(
            cq["n_support"], cq["n"], 100 * cq["frac"], cq["p"],
            cg["n_support"], cg["n"], 100 * cg["frac"], cg["p"],
            cs["n_support"], cs["n"]),
        transform=ax.transAxes, fontsize=7.4, color=INK2, va="top")

# ---------------------------------------------------------------- (d) siRNA vs CRISPR
ax = fig.add_subplot(gs[2, 1]); tag(ax, "d")
try:
    Cw = pd.read_csv("results/crispr_risk36_by_cellline.csv").set_index("gene")
    sh = [g for g in W["gene"] if g in Cw.index]
    a = W.set_index("gene").loc[sh, "mean_delta_qc"].to_numpy(float)
    b = Cw.loc[sh, ["delta_HCT116", "delta_RKO"]].mean(axis=1).to_numpy(float)
    ver = W.set_index("gene").loc[sh, "twas_support"].to_numpy()
    from scipy.stats import spearmanr
    rho, prho = spearmanr(a, b)
    ax.axhline(0, color="#d8d7d2", lw=0.8); ax.axvline(0, color="#d8d7d2", lw=0.8)
    for key, lab, s in [(NS, "no significant change", 24),
                        (SUPPORTS, "supports TWAS", 48), (OPPOSES, "opposes TWAS", 48)]:
        m = ver == key
        if m.any():
            ax.scatter(a[m], b[m], s=s, color=VC[key] if key != NS else "#c9c7c1",
                       edgecolor="white", linewidth=0.5, zorder=3 + (key != NS), label=lab)
    for g, xa, yb, v in zip(sh, a, b, ver):
        if v != NS or abs(xa) > 32 or abs(yb) > 16:
            ax.annotate(g, (xa, yb), textcoords="offset points", xytext=(4, 3),
                        fontsize=6.6, color=INK2)
    ax.set_xlabel("siRNA paired Δ viability (HCT116 + RKO, points)")
    ax.set_ylabel("CRISPR paired Δ viability\n(HCT116 + RKO)")
    ax.set_title(f"The two modalities agree only loosely (n = {len(sh)})")
    same = int(np.sum(np.sign(a) == np.sign(b)))
    ax.text(0.02, 0.06, f"Spearman ρ = {rho:+.2f}, P = {prho:.2f}\nsame sign {same}/{len(sh)}",
            transform=ax.transAxes, fontsize=7.6, color=INK2, va="bottom")
    ax.legend(loc="lower right", fontsize=7.0)
    ax.grid(True); ax.set_axisbelow(True)
except Exception as e:                                          # pragma: no cover
    ax.text(0.5, 0.5, f"cross-assay panel unavailable\n{e}", ha="center", va="center")
    ax.axis("off")

fig.savefig("figures/fig21_sirna_by_cellline.png", dpi=220, bbox_inches="tight")
fig.savefig("figures/fig21_sirna_by_cellline.pdf", bbox_inches="tight")
plt.close(fig)
print("fig21 done")
