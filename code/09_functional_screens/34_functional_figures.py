#!/usr/bin/env python3
"""34 - Figure 19: the siRNA and CRISPR viability screens."""
import pandas as pd, numpy as np, json, os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

os.makedirs("figures", exist_ok=True)
C = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#4a3aa7", "#e34948"]
POS, NEG = "#e34948", "#2a78d6"
INK, INK2, INK3 = "#0b0b0b", "#52514e", "#8a8880"
SURF = "#fcfcfb"
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


def tag(ax, s):
    ax.text(-0.11, 1.10, s, transform=ax.transAxes, fontsize=11, fontweight="bold",
            color=INK, va="top", ha="left")


QC = pd.read_csv("results/func_qc.csv")
SI = pd.read_csv("results/func_sirna_gene_stats.csv")
CR = pd.read_csv("results/func_crispr_gene_stats.csv")
RG = pd.read_csv("results/func_risk_gene_results.csv")
MM = pd.read_csv("results/func_sirna_vs_crispr.csv")
F = json.load(open("results/func_summary.json"))
pri = pd.read_csv("results/integrated_gene_priority.csv").set_index("gene")
EXCL = F["excluded_cell_line"]
LINES = ["HCT116", "RKO", "SW480"]
LCOL = {"HCT116": C[0], "RKO": C[2], "SW480": C[1]}

fig = plt.figure(figsize=(15.4, 9.0))
gs = fig.add_gridspec(2, 3, hspace=0.50, wspace=0.30)

# ---------------------------------------------------------------- (a) QC
ax = fig.add_subplot(gs[0, 0]); tag(ax, "a")
x = np.arange(len(LINES))
for k, assay in enumerate(["siRNA", "CRISPR"]):
    for i, line in enumerate(LINES):
        d = QC[(QC.assay == assay) & (QC.cell_line == line)]
        xs = i + (k - 0.5) * 0.36
        ax.plot([xs] * len(d), d["z_prime"], "o", ms=5.5,
                color=LCOL[line] if assay == "siRNA" else "white",
                markeredgecolor=LCOL[line], markeredgewidth=1.3, zorder=3)
        ax.plot([xs - 0.11, xs + 0.11], [d["z_prime"].mean()] * 2,
                color=LCOL[line], lw=2.2)
ax.axhspan(-1, 0, color="#f6e0e0", zorder=0)
ax.axhline(0, color=POS, lw=1, ls=(0, (3, 3)))
ax.axhline(0.5, color=INK3, lw=1, ls=(0, (3, 3)))
ax.text(2.48, 0.52, "0.5 excellent", fontsize=7, color=INK3, ha="right")
ax.text(0.02, 0.055, "below 0 the controls overlap", transform=ax.transAxes,
        fontsize=7, color=POS)
ax.set_xticks(x); ax.set_xticklabels(LINES)
ax.set_ylabel("Z'-factor (per replicate)")
ax.set_ylim(-0.35, 1.0)
ax.set_title("Assay quality: can the plate tell a hit from nothing?")
ax.legend(handles=[Line2D([], [], marker="o", ls="", color=INK3, label="siRNA"),
                   Line2D([], [], marker="o", ls="", mfc="white", mec=INK3, label="CRISPR")],
          loc="upper left", ncol=2)
ax.annotate(f"{EXCL} excluded", xy=(2.18, -0.10), xytext=(1.30, -0.28),
            fontsize=7.8, color=POS, fontweight="bold",
            arrowprops=dict(arrowstyle="->", color=POS, lw=1))
ax.grid(axis="y"); ax.set_axisbelow(True)

# ---------------------------------------------------------------- (b) ranked effects
# A volcano would be uninformative here: the p-value comes from an empirical null on the
# same statistic, so -log10 p is a deterministic function of |z| and the plot is a V.
# The ranking carries the information instead.
ax = fig.add_subplot(gs[0, 1]); tag(ax, "b")
d = SI.sort_values("mean_z").reset_index(drop=True)
d["rank"] = np.arange(len(d))
oth = d[~d.is_risk & ~d.core_machinery]
core = d[d.core_machinery]
rk = d[d.is_risk & ~d.core_machinery]
ax.scatter(oth["rank"], oth["mean_z"], s=10, color="#d8d7d2", edgecolor="none",
           label=f"other genes ({len(oth)})")
ax.scatter(core["rank"], core["mean_z"], s=30, color=C[3], edgecolor="none",
           label=f"ribosome / proteasome / spliceosome ({len(core)})")
ax.scatter(rk["rank"], rk["mean_z"], s=26, color=C[0], edgecolor="none",
           label=f"CRC risk genes ({len(rk)})")
det = d[d.detected]
zc_lo = det[det.mean_z < 0]["mean_z"].max() if (det.mean_z < 0).any() else np.nan
zc_hi = det[det.mean_z > 0]["mean_z"].min() if (det.mean_z > 0).any() else np.nan
for zc in (zc_lo, zc_hi):
    if np.isfinite(zc):
        ax.axhline(zc, color=INK3, lw=0.8, ls=(0, (3, 3)))
ax.axhline(-2, color=POS, lw=0.8, ls=(0, (1, 3)))
ax.text(len(d) - 2, -2.12, "|z| = 2", fontsize=6.8, color=POS, ha="right", va="top")
ax.text(len(d) - 2, zc_lo - 0.02, "FDR 0.05", fontsize=6.8, color=INK3,
        ha="right", va="top")
lab = pd.concat([d[d.detected & (d.is_risk | d.core_machinery)].nsmallest(6, "mean_z"),
                 d[d.detected & d.is_risk].nlargest(2, "mean_z")])
for i, (_, r) in enumerate(lab.iterrows()):
    ax.annotate(r["gene"], (r["rank"], r["mean_z"]), textcoords="offset points",
                xytext=(14, [-16, -4, 8, 20, -28, 32, 6, -10][i % 8]), fontsize=6.8,
                color=INK2,
                arrowprops=dict(arrowstyle="-", lw=0.4, color=INK3, shrinkA=0, shrinkB=2))
ax.set_xlabel("Genes ranked by effect")
ax.set_ylabel("Effect on viability (robust z, siRNA)")
ax.set_title("siRNA screen, 333 genes in HCT116 + RKO")
ax.legend(loc="lower right", fontsize=7)
ax.grid(axis="y"); ax.set_axisbelow(True)

# ---------------------------------------------------------------- (c) calibration
ax = fig.add_subplot(gs[0, 2]); tag(ax, "c")
groups = [("Ribosome etc.", d[d.core_machinery]["mean_z"], C[3]),
          ("CRC risk genes", d[d.is_risk]["mean_z"], C[0]),
          ("Other genes", d[~d.is_risk & ~d.core_machinery]["mean_z"], INK3)]
for i, (nm, vals, col) in enumerate(groups):
    parts = ax.violinplot([vals], positions=[i], widths=0.7, showextrema=False)
    for pc in parts["bodies"]:
        pc.set_facecolor(col); pc.set_alpha(0.38); pc.set_edgecolor("none")
    ax.plot([i - 0.2, i + 0.2], [np.median(vals)] * 2, color=col, lw=2.4)
    ax.text(i, 2.3, f"n={len(vals)}", ha="center", fontsize=7, color=INK3)
ax.axhline(0, color=INK3, lw=0.8)
cal = F["calibration"]["siRNA"]
ax.set_xticks(range(3)); ax.set_xticklabels([g[0] for g in groups], fontsize=8)
ax.set_ylabel("Effect on viability (robust z)")
ax.set_ylim(-5.2, 2.8)
ax.set_title("The screen finds real dependencies")
ax.text(0.03, 0.03, f"core machinery vs rest: P = {cal['p']:.1g}\n"
                    f"{cal['strong_core']}/{cal['n_core']} strong hits vs "
                    f"{cal['strong_rest']}/{len(d)-cal['n_core']}",
        transform=ax.transAxes, fontsize=7.2, color=INK2, va="bottom")
ax.grid(axis="y"); ax.set_axisbelow(True)

# ---------------------------------------------------------------- (d) the risk genes
ax = fig.add_subplot(gs[1, 0]); tag(ax, "d")
h = RG[RG["detected"]].sort_values("mean_pct_NT")
cols = [POS if c == "oncogene-like" else (NEG if c == "tumour-suppressor-like" else INK3)
        for c in h["twas_class"].fillna("")]
y = np.arange(len(h))
ax.barh(y, h["mean_pct_NT"] - 100, left=100, color=cols, height=0.7)
ax.axvline(100, color=INK3, lw=1)
ax.set_yticks(y); ax.set_yticklabels(h["gene"], fontsize=7.4)
for yi, (_, r) in zip(y, h.iterrows()):
    if r["strong"]:
        side = -2.5 if r["mean_pct_NT"] < 100 else 2.5
        ax.text(r["mean_pct_NT"] + side, yi, "strong", va="center",
                ha="right" if side < 0 else "left", fontsize=6.4, color=INK3)
ax.set_xlabel("Viability after knockdown (% of non-targeting control)")
ax.set_title(f"{len(h)} of {len(RG)} risk genes with a detectable effect")
ax.set_xlim(8, 132)
ax.legend(handles=[Patch(color=POS, label="oncogene-like"),
                   Patch(color=NEG, label="tumour-suppressor-like"),
                   Patch(color=INK3, label="not spatially validated")],
          loc="center right", fontsize=7, bbox_to_anchor=(1.0, 0.33))
ax.grid(axis="x"); ax.set_axisbelow(True)

# ---------------------------------------------------------------- (e) predicted vs measured
ax = fig.add_subplot(gs[1, 1]); tag(ax, "e")
for R, nm, col, key in [(SI, "siRNA", C[0], "siRNA"), (CR, "CRISPR", C[3], "CRISPR")]:
    m = R.join(pri[["ko_impact"]], on="gene").dropna(subset=["ko_impact"])
    st = F[key]["ko_impact_vs_measured"]
    ax.scatter(m["ko_impact"], m["mean_z"], s=20, color=col, alpha=0.75, edgecolor="none",
               label=f"{nm}: ρ = {st['rho']:+.2f}, P = {st['p']:.3f} (n = {st['n']})")
    if len(m) > 5:
        b = np.polyfit(m["ko_impact"], m["mean_z"], 1)
        xs = np.linspace(m["ko_impact"].min(), m["ko_impact"].max(), 20)
        ax.plot(xs, np.polyval(b, xs), color=col, lw=1.4, ls=(0, (4, 2)))
ax.axhline(0, color=INK3, lw=0.8)
ax.set_xlabel("Predicted knockout impact (in silico, §2.2)")
ax.set_ylabel("Measured effect on viability (robust z)")
ax.set_title("Does the in-silico knockout predict the real one?")
ax.legend(loc="lower left", fontsize=7)
ax.grid(); ax.set_axisbelow(True)

# ---------------------------------------------------------------- (f) siRNA vs CRISPR
ax = fig.add_subplot(gs[1, 2]); tag(ax, "f")
con = F["concordance"]
ok = MM[~(MM.detected_si | MM.detected_cr)]
hit = MM[MM.detected_si | MM.detected_cr]
ax.axhline(0, color="#e0dfda", lw=0.8); ax.axvline(0, color="#e0dfda", lw=0.8)
ax.scatter(ok["mean_z_si"], ok["mean_z_cr"], s=13, color="#d8d7d2", edgecolor="none",
           label="neither detected")
ag = hit[hit.agree_direction]
dis = hit[~hit.agree_direction]
ax.scatter(ag["mean_z_si"], ag["mean_z_cr"], s=24, color=C[2], edgecolor="none",
           label=f"detected, same direction ({len(ag)})")
ax.scatter(dis["mean_z_si"], dis["mean_z_cr"], s=24, color=POS, edgecolor="none",
           label=f"detected, opposite ({len(dis)})")
for _, r in dis[dis.is_risk_si].iterrows():
    ax.annotate(r["gene"], (r["mean_z_si"], r["mean_z_cr"]), textcoords="offset points",
                xytext=(6, 3), fontsize=6.8, color=INK2)
ax.set_xlabel("siRNA effect (robust z)")
ax.set_ylabel("CRISPR effect (robust z)")
ax.set_title(f"The two modalities disagree ({100*con['frac_direction_agree']:.0f} % same direction)")
ax.text(0.03, 0.03, f"ρ = {con['rho']:+.2f}, P = {con['p']:.2f}\n{con['n_both']} genes in both",
        transform=ax.transAxes, fontsize=7.2, color=INK2, va="bottom")
ax.legend(loc="lower right", fontsize=7)
ax.grid(); ax.set_axisbelow(True)

fig.savefig("figures/fig19_functional_screens.png", dpi=220, bbox_inches="tight")
fig.savefig("figures/fig19_functional_screens.pdf", bbox_inches="tight")
plt.close(fig)
print("fig19 done")
