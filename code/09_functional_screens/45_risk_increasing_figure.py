#!/usr/bin/env python3
"""
45 - Figure 26: the TWAS risk-increasing genes only, across both screens.

This subset is not an arbitrary slice. The germline direction makes a one-sided
prediction for these genes and only for these: a risk-increasing gene raises risk through
higher activity, so if that effect carried into the tumour cell the cell should NEED the
gene and removing it should LOWER viability. Protective genes predict the opposite, and
they are also the less actionable half - restoring a gene is not something a small
molecule does - so the risk-increasing set is both the cleaner test and the one a target
list would be drawn from.

Because the prediction is one-sided, this figure uses a one-sided test where Figures 20
and 21 use two-sided ones: among the significant changes, are more of them decreases than
a coin flip would give? That is the question the germline hypothesis actually poses.

The matrix is the union of the two panels, 40 genes, with the cells of a gene that an
assay did not measure hatched rather than left blank, so a reader is never misled into
reading absence as a null result.

Outputs
  results/risk_increasing_both_assays.csv
  figures/fig26_risk_increasing.{png,pdf}
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
POS, NEG = "#e34948", "#2a78d6"
INK, INK2, INK3 = "#0b0b0b", "#52514e", "#8a8880"
SURF = "#fcfcfb"
DIV = LinearSegmentedColormap.from_list("div", ["#184f95", "#6da7ec", "#f2f1ed",
                                                "#f0a09f", "#c3202f"])
RI = "#8c2d3f"
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


C = pd.read_csv("results/crispr_risk36_by_cellline.csv").set_index("gene")
S = pd.read_csv("results/sirna_risk_by_cellline.csv").set_index("gene")
ci = C[C.twas_direction == "risk-increasing"]
si = S[S.twas_direction == "risk-increasing"]
genes = sorted(set(ci.index) | set(si.index))

rows = []
for g in genes:
    r = {"gene": g,
         "class": (ci.loc[g, "class"] if g in ci.index else si.loc[g, "class"]),
         "in_crispr": g in ci.index, "in_sirna": g in si.index}
    for line in LINES:
        r[f"cr_{line}"] = ci.loc[g, f"delta_{line}"] if g in ci.index else np.nan
        r[f"cr_sig_{line}"] = bool(ci.loc[g, f"nom_{line}"]) if g in ci.index else False
        r[f"si_{line}"] = si.loc[g, f"delta_{line}"] if g in si.index else np.nan
        r[f"si_sig_{line}"] = bool(si.loc[g, f"nom_{line}"]) if g in si.index else False
    rows.append(r)
W = pd.DataFrame(rows)
W["n_sig"] = W[[f"cr_sig_{l}" for l in LINES] + [f"si_sig_{l}" for l in LINES]].sum(axis=1)
W["mean_delta"] = W[[f"cr_{l}" for l in LINES] + [f"si_{l}" for l in LINES]].mean(axis=1)


def direction(r):
    """Sign of the significant changes for this gene, across everything measured."""
    s = []
    for pre in ["cr", "si"]:
        for line in LINES:
            if r[f"{pre}_sig_{line}"] and np.isfinite(r[f"{pre}_{line}"]):
                s.append(np.sign(r[f"{pre}_{line}"]))
    if not s:
        return "no significant change"
    if len(set(s)) > 1:
        return "conflicting"
    return "reduces viability" if s[0] < 0 else "increases viability"


W["verdict"] = W.apply(direction, axis=1)
W.to_csv("results/risk_increasing_both_assays.csv", index=False)

# ---- the one-sided question, over the significant calls
calls = []
for _, r in W.iterrows():
    for pre, assay in [("cr", "CRISPR"), ("si", "siRNA")]:
        for line in LINES:
            if r[f"{pre}_sig_{line}"] and np.isfinite(r[f"{pre}_{line}"]):
                calls.append({"gene": r.gene, "assay": assay, "cell_line": line,
                              "delta": r[f"{pre}_{line}"],
                              # SW480 fails the plate-level Z' check in the siRNA assay
                              "qc_pass": not (assay == "siRNA" and line == "SW480")})
Cl = pd.DataFrame(calls)
QC = Cl[Cl.qc_pass]
qc_down, qc_tot = int((QC.delta < 0).sum()), len(QC)
qc_bt = binomtest(qc_down, qc_tot, 0.5, alternative="greater")
n_down = int((Cl.delta < 0).sum())
n_tot = len(Cl)
bt = binomtest(n_down, n_tot, 0.5, alternative="greater")
by_assay = {a: (int((d.delta < 0).sum()), len(d)) for a, d in Cl.groupby("assay")}
gd = W["verdict"].value_counts().to_dict()
g_down = gd.get("reduces viability", 0)
g_up = gd.get("increases viability", 0)
gbt = binomtest(g_down, g_down + g_up, 0.5, alternative="greater") if (g_down + g_up) else None

# ---- are risk-increasing genes shifted against the rest? CRISPR panel, all genes
L = pd.read_csv("results/crispr_by_cellline.csv")
Lg = L.groupby("gene").agg(delta=("delta_paired", "mean"),
                           twas=("twas_direction", "first"),
                           is_risk=("is_risk", "first")).reset_index()
grp_ri = Lg[(Lg.is_risk) & (Lg.twas == "risk-increasing")]["delta"]
grp_pr = Lg[(Lg.is_risk) & (Lg.twas == "protective")]["delta"]
grp_ot = Lg[~Lg.is_risk]["delta"]
u_rp = mannwhitneyu(grp_ri, grp_pr, alternative="less")
u_ro = mannwhitneyu(grp_ri, grp_ot, alternative="less")

J = {"n_genes": len(W), "n_crispr": int(W.in_crispr.sum()), "n_sirna": int(W.in_sirna.sum()),
     "n_both": int((W.in_crispr & W.in_sirna).sum()),
     "calls_total": n_tot, "calls_down": n_down, "calls_p_one_sided": float(bt.pvalue),
     "by_assay": {k: list(v) for k, v in by_assay.items()},
     "qc_calls_down": qc_down, "qc_calls_total": qc_tot,
     "qc_p_one_sided": float(qc_bt.pvalue),
     "up_calls_all_from": sorted(set(
         f"{r.assay} {r.cell_line}" for _, r in Cl[Cl.delta > 0].iterrows())),
     "gene_verdict": {k: int(v) for k, v in gd.items()},
     "gene_p_one_sided": float(gbt.pvalue) if gbt else None,
     "median_delta_risk_increasing": float(grp_ri.median()),
     "median_delta_protective": float(grp_pr.median()),
     "median_delta_other": float(grp_ot.median()),
     "p_vs_protective": float(u_rp.pvalue), "p_vs_other": float(u_ro.pvalue),
     "genes_reducing": W.loc[W.verdict == "reduces viability", "gene"].tolist(),
     "genes_increasing": W.loc[W.verdict == "increases viability", "gene"].tolist(),
     "genes_conflicting": W.loc[W.verdict == "conflicting", "gene"].tolist()}
json.dump(J, open("results/risk_increasing_summary.json", "w"), indent=1)
print(json.dumps({k: v for k, v in J.items() if not k.startswith("genes_")}, indent=1))
print("reduces:", J["genes_reducing"])
print("increases:", J["genes_increasing"], "| conflicting:", J["genes_conflicting"])

ORDER = {"reduces viability": 0, "conflicting": 1, "increases viability": 2,
         "no significant change": 3}
W["rank"] = W["verdict"].map(ORDER)
W = W.sort_values(["rank", "mean_delta"], ascending=[False, False]).reset_index(drop=True)
N = len(W)

fig = plt.figure(figsize=(14.0, 11.6))
gs = fig.add_gridspec(3, 2, width_ratios=[1.40, 1.0], height_ratios=[1, 1, 1],
                      wspace=0.42, hspace=0.72)

# ---------------------------------------------------------------- (a) matrix
ax = fig.add_subplot(gs[:, 0]); tag(ax, "a", dx=-0.28, dy=1.02)
M = W[[f"cr_{l}" for l in LINES] + [f"si_{l}" for l in LINES]].to_numpy(float)
lim = float(np.nanpercentile(np.abs(M[np.isfinite(M)]), 97))
im = ax.imshow(np.ma.masked_invalid(M), cmap=DIV,
               norm=TwoSlopeNorm(vmin=-lim, vcenter=0, vmax=lim),
               aspect="auto", origin="lower", extent=(-0.5, 5.5, -0.5, N - 0.5))
for i in range(N):
    for j in range(6):
        if not np.isfinite(M[i, j]):
            ax.add_patch(Rectangle((j - 0.5, i - 0.5), 1, 1, facecolor="#f0efec",
                                   edgecolor="#d8d7d2", linewidth=0.4, hatch="////",
                                   zorder=3))
ax.axvline(2.5, color=SURF, lw=5, zorder=6)
ax.axvline(2.5, color=INK3, lw=1.0, zorder=7)
VCOL = {"reduces viability": DOWN, "increases viability": UP,
        "conflicting": "#8a8880", "no significant change": "#e4e2dd"}
for i, v in enumerate(W["verdict"]):
    ax.add_patch(Rectangle((5.72, i - 0.5), 0.5, 1, facecolor=VCOL[v],
                           edgecolor=SURF, linewidth=0.5, zorder=2))
ax.text(5.97, N - 0.2, "effect", ha="center", va="bottom", fontsize=7.2, color=INK2,
        fontweight="bold")
for j, (pre, line) in enumerate([(p, l) for p in ["cr", "si"] for l in LINES]):
    for i in np.flatnonzero(W[f"{pre}_sig_{line}"].to_numpy(bool)):
        ax.add_patch(Rectangle((j - 0.5, i - 0.5), 1, 1, fill=False, edgecolor=INK,
                               linewidth=1.5, zorder=5))
        v = W[f"{pre}_{line}"].iloc[i]
        ax.text(j, i, f"{v:+.0f}", ha="center", va="center", fontsize=6.0,
                color="white" if abs(v) > lim * 0.55 else INK, zorder=8)
ax.set_xticks(range(6)); ax.set_xticklabels(LINES + LINES, fontsize=7.6)
ax.text(1, -0.075, "CRISPR-Cas9", ha="center", va="top", fontsize=9, color=INK,
        fontweight="bold", transform=ax.get_xaxis_transform())
ax.text(4, -0.075, "siRNA", ha="center", va="top", fontsize=9, color=INK,
        fontweight="bold", transform=ax.get_xaxis_transform())
ax.set_yticks(range(N))
labs = [f"{ {'oncogene-like': '▲', 'tumour-suppressor-like': '▼'}.get(c, ' ') } {g}"
        for g, c in zip(W["gene"], W["class"])]
ax.set_yticklabels(labs, fontsize=7.0)
for t, cls in zip(ax.get_yticklabels(), W["class"]):
    t.set_color(POS if cls == "oncogene-like" else INK2)
ax.set_xlim(-0.55, 6.40); ax.set_ylim(-0.5, N - 0.5)
ax.set_title(f"The {N} TWAS risk-increasing genes assayed in either screen", pad=14)
ax.tick_params(length=0)
for sp in ax.spines.values():
    sp.set_visible(False)
cb = fig.colorbar(im, ax=ax, fraction=0.034, pad=0.075, shrink=0.34, anchor=(0.0, 1.0))
cb.set_label("Δ viability vs same-plate control (points)", fontsize=7.4)
cb.ax.tick_params(labelsize=7); cb.outline.set_visible(False)
ax.legend(handles=[
    Patch(facecolor=DOWN, label="removal lowers viability — as predicted"),
    Patch(facecolor=UP, label="removal raises it — against the prediction"),
    Patch(facecolor="#8a8880", label="conflicting between plates"),
    Patch(facecolor="#e4e2dd", label="no significant change"),
    Line2D([], [], marker="s", ls="", mfc="none", mec=INK, ms=9, markeredgewidth=1.5,
           label="paired t-test P < 0.05"),
    Patch(facecolor="#f0efec", hatch="////", edgecolor="#d8d7d2",
          label="not on that panel")],
    loc="upper left", bbox_to_anchor=(-0.31, -0.085), ncol=2, fontsize=7.2,
    columnspacing=1.0)

# ---------------------------------------------------------------- (b) the one-sided test
ax = fig.add_subplot(gs[0, 1]); tag(ax, "b")
groups_b = [("CRISPR\nall 3 lines", Cl[Cl.assay == "CRISPR"]),
            ("siRNA\nHCT116 + RKO", Cl[(Cl.assay == "siRNA") & (Cl.cell_line != "SW480")]),
            ("siRNA\nSW480 (fails QC)", Cl[(Cl.assay == "siRNA") & (Cl.cell_line == "SW480")]),
            ("all QC-passing\nplates", QC)]
x = np.arange(len(groups_b), dtype=float)
downs = [int((d.delta < 0).sum()) for _, d in groups_b]
ups = [int((d.delta > 0).sum()) for _, d in groups_b]
tots = [a + b for a, b in zip(downs, ups)]
ax.bar(x, downs, width=0.6, color=DOWN, label="viability falls (predicted)")
ax.bar(x, ups, width=0.6, bottom=downs, color=UP, label="viability rises")
for xi, d, u_, t in zip(x, downs, ups, tots):
    if d:
        ax.text(xi, d / 2, str(d), ha="center", va="center", fontsize=8.6, color="white",
                fontweight="bold")
    if u_:
        ax.text(xi, d + u_ / 2, str(u_), ha="center", va="center", fontsize=8.6,
                color="white", fontweight="bold")
    if t:
        ax.text(xi, t + 0.6, f"{100 * d / t:.0f} % down", ha="center", fontsize=7.4,
                color=INK2)
ax.set_xticks(x); ax.set_xticklabels([g[0] for g in groups_b], fontsize=7.8)
ax.get_xticklabels()[2].set_color("#b3472f")
ax.set_ylabel("Significant changes")
ax.set_ylim(0, max(tots) * 1.30)
ax.set_title("Direction of every significant change")
ax.legend(loc="upper left", fontsize=7.2)
ax.grid(axis="y"); ax.set_axisbelow(True)
ax.text(0.0, -0.36,
        "The germline direction predicts a fall, so this is the one place a one-sided test is "
        "legitimate. Every single change that\nruns the other way comes from the siRNA SW480 "
        "plates, the set that fails the Z' check. Across all {} calls {} are\ndecreases "
        "(P = {:.1e}); across the {} calls on QC-passing plates, {} are ({}/{}, "
        "P = {:.1e}).".format(
            n_tot, n_down, bt.pvalue, qc_tot, qc_down, qc_down, qc_tot, qc_bt.pvalue),
        transform=ax.transAxes, fontsize=7.4, color=INK2, va="top")

# ---------------------------------------------------------------- (c) shift vs others
ax = fig.add_subplot(gs[1, 1]); tag(ax, "c")
groups = [("risk-increasing\n(n = %d)" % len(grp_ri), grp_ri.to_numpy(), RI),
          ("protective\n(n = %d)" % len(grp_pr), grp_pr.to_numpy(), "#1f6f8b"),
          ("non-TWAS panel\n(n = %d)" % len(grp_ot), grp_ot.to_numpy(), "#7d99b8")]
for i, (lab, v, col) in enumerate(groups):
    parts = ax.violinplot([v], positions=[i], widths=0.74, showextrema=False)
    for b in parts["bodies"]:
        b.set_facecolor(col); b.set_alpha(0.26); b.set_edgecolor("none")
    ax.scatter(np.full(len(v), i) + np.random.default_rng(45 + i).normal(0, 0.055, len(v)),
               v, s=10, color=col, alpha=0.7, edgecolor="none", zorder=3)
    ax.plot([i - 0.26, i + 0.26], [np.median(v)] * 2, color=INK, lw=1.8, zorder=4)
    ax.text(i, np.median(v), f"{np.median(v):+.1f}", ha="center", va="center", fontsize=7.2,
            color=INK, zorder=6, bbox=dict(facecolor=SURF, edgecolor="none", pad=1.4))
ax.axhline(0, color="#d8d7d2", lw=0.9, zorder=1)
ax.set_xticks(range(3)); ax.set_xticklabels([g[0] for g in groups], fontsize=8.0)
ax.set_ylabel("Mean paired Δ viability\nacross the CRISPR plates (points)")
ax.set_title("Are they shifted, as a group?")
ax.grid(axis="y"); ax.set_axisbelow(True)
ax.text(0.0, -0.40,
        "On the CRISPR panel the risk-increasing genes sit at a median of {:+.1f} points against "
        "{:+.1f} for the protective\ngenes (one-sided P = {:.2f}) and {:+.1f} for the non-TWAS "
        "genes assayed beside them (P = {:.2f}). The direction\nof the individual significant "
        "calls is the signal here; the group as a whole is not displaced.".format(
            J["median_delta_risk_increasing"], J["median_delta_protective"], u_rp.pvalue,
            J["median_delta_other"], u_ro.pvalue),
        transform=ax.transAxes, fontsize=7.4, color=INK2, va="top")

# ---------------------------------------------------------------- (d) gene verdicts
ax = fig.add_subplot(gs[2, 1]); tag(ax, "d")
order = ["reduces viability", "increases viability", "conflicting", "no significant change"]
vals = [gd.get(k, 0) for k in order]
ax.bar(range(4), vals, width=0.62, color=[VCOL[k] for k in order])
for i, v in enumerate(vals):
    if v:
        ax.text(i, v + 0.4, str(v), ha="center", fontsize=8.4, color=INK2)
ax.set_xticks(range(4))
ax.set_xticklabels(["lowers\nviability", "raises it", "conflicting", "no change"],
                   fontsize=8.0)
ax.set_ylabel(f"Genes (of {N})")
ax.set_ylim(0, max(vals) * 1.25)
ax.set_title("Per-gene outcome")
ax.grid(axis="y"); ax.set_axisbelow(True)
ax.text(0.0, -0.38,
        "Genes lowering viability: {}.\nRaising it: {}.{}".format(
            ", ".join(J["genes_reducing"]) or "none",
            ", ".join(J["genes_increasing"]) or "none",
            (" Conflicting between plates: " + ", ".join(J["genes_conflicting"]) + ".")
            if J["genes_conflicting"] else ""),
        transform=ax.transAxes, fontsize=7.4, color=INK2, va="top", wrap=True)

fig.savefig("figures/fig26_risk_increasing.png", dpi=220, bbox_inches="tight")
fig.savefig("figures/fig26_risk_increasing.pdf", bbox_inches="tight")
plt.close(fig)
print("\nfig26 done")
