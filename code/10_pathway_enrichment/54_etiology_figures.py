#!/usr/bin/env python3
"""
54 - Figures 35 and 36: the whole risk-gene set, and what it says about CRC aetiology.

Figure 35 is the pathway analysis of all 347 genes at once: what is enriched, whether it
beats both nulls, whether it survives splitting by direction, and how many risk genes sit
inside each canonical CRC family whether or not that family is enriched.

Figure 36 is the aetiology view: gene sets regrouped by the exposure biology colorectal
cancer is known to run on, how much of the risk-gene set the pathway view can see at all,
and the genes themselves.
"""
import pandas as pd, numpy as np, json, os, textwrap
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

os.makedirs("figures", exist_ok=True)
INK, INK2, INK3 = "#0b0b0b", "#52514e", "#8a8880"
SURF = "#fcfcfb"
GREEN, GREY, AMBER = "#1baf7a", "#c9c7c1", "#e8a33d"
UPC, DOWNC = "#b0455a", "#2f6f9f"
BLOCKC = {"MHC class I and the proteasome": "#6a4c93", "BMP / TGF-β": "#2e8b6d",
          "Glucuronidation": "#c07a2c", "Other": "#8a8880"}
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

E = pd.read_csv("results/all347_enrichment.csv")
NUL = pd.DataFrame(json.load(open("results/all347_nulls.json")))
FS = pd.read_csv("results/canonical_membership_summary.csv")
TH = pd.read_csv("results/etiology_themes.csv")
AS = pd.read_csv("results/etiology_gene_assignment.csv")
J = json.load(open("results/all347_etiology.json"))
S = E[E.fdr < 0.1].sort_values("p")


def block_of(genes):
    g = set(str(genes).split("|"))
    if any(x.startswith("HLA") for x in g) or (g & {"PSMA7", "PSMD3", "PSMD13", "UBE2B"}):
        return "MHC class I and the proteasome"
    if g & {"BMP2", "BMP4", "SMAD6", "SMAD7", "SMAD9", "CTDSP2"}:
        return "BMP / TGF-β"
    if any(x.startswith("UGT") for x in g):
        return "Glucuronidation"
    return "Other"


S = S.assign(block=[block_of(g) for g in S.genes])


def tag(ax, s, dx=-0.10, dy=1.05):
    ax.text(dx, dy, s, transform=ax.transAxes, fontsize=11.5, fontweight="bold",
            color=INK, va="bottom", ha="left")


# ==================================================================== FIGURE 35
fig = plt.figure(figsize=(7.9, 9.0))
gs = fig.add_gridspec(3, 1, height_ratios=[1.12, 0.62, 1.14], hspace=1.05)

# ---------------------------------------------------------------- a enrichments
ax = fig.add_subplot(gs[0, 0]); tag(ax, "a", dx=-0.32, dy=1.06)
y = np.arange(len(S))[::-1]
ax.barh(y, -np.log10(S.fdr), height=0.62, color=[BLOCKC[b] for b in S.block], zorder=3)
for yi, r in zip(y, S.itertuples(index=False)):
    ax.text(-np.log10(r.fdr) + 0.03, yi, r.genes.replace("|", ", "), va="center",
            fontsize=6.0, color=INK2)
ax.set_yticks(y)
ax.set_yticklabels([p if len(p) <= 46 else p[:44] + "…" for p in S.pathway], fontsize=6.8)
ax.axvline(-np.log10(0.1), color="#c9c7c1", lw=0.9, ls="--", zorder=2)
ax.set_xlabel("−log$_{10}$ FDR")
ax.set_xlim(0, 5.6)
ax.set_title("Every gene set enriched across all 347 risk genes", pad=10)
ax.grid(axis="x"); ax.set_axisbelow(True)
ax.legend(handles=[Patch(color=c, label=b) for b, c in BLOCKC.items() if b != "Other"],
          loc="upper left", bbox_to_anchor=(0.0, -0.255), ncol=3, fontsize=7.0)
ax.text(0.0, -0.40,
        "Eleven gene sets, three blocks of genes. The dashed line is the FDR < 0.1 "
        "threshold; the genes carrying each set are printed beside it.",
        transform=ax.transAxes, fontsize=7.0, color=INK2, va="top")

# ---------------------------------------------------------------- b nulls
ax = fig.add_subplot(gs[1, 0]); tag(ax, "b", dx=-0.32, dy=1.12)
yy = np.arange(len(NUL))[::-1]
ax.barh(yy, NUL.observed, height=0.56,
        color=[UPC if p < 0.05 else "#c9c7c1" for p in NUL.p_perm_matched], zorder=3,
        label="observed (FDR < 0.1)")
ax.errorbar(NUL.null_mean, yy + 0.15, xerr=NUL.null_sd, fmt="o", color="#9a9893", ms=3.4,
            lw=0.9, capsize=2.2, zorder=4, label="null: size-matched")
ax.errorbar(NUL.matched_mean, yy - 0.15, xerr=NUL.matched_sd, fmt="D", color="#7a3b8f",
            ms=3.4, lw=0.9, capsize=2.2, zorder=4, label="null: matched on annotation depth")
for yi, r in zip(yy, NUL.itertuples(index=False)):
    ax.text(max(r.observed, r.matched_mean + r.matched_sd, r.null_mean + r.null_sd) + 0.3,
            yi, f"P = {r.p_perm_matched:.3f}", va="center", fontsize=7.0,
            color=GREEN if r.p_perm_matched < 0.05 else INK2)
ax.set_yticks(yy)
ax.set_yticklabels([f"{r.gene_list}  ({int(r.n_mapped)} genes)"
                    for r in NUL.itertuples(index=False)], fontsize=7.4)
ax.set_xlabel("Enriched gene sets")
ax.set_xlim(0, 18)
ax.set_title("The whole set beats both nulls; neither half does on its own", pad=10)
ax.legend(loc="upper left", bbox_to_anchor=(0.0, -0.46), ncol=3, fontsize=7.0,
          columnspacing=1.4, handletextpad=0.4)
ax.grid(axis="x"); ax.set_axisbelow(True)
ax.text(0.0, -0.74,
        "Splitting the genes by the direction of risk halves each list and the enrichment "
        "disappears, which is the honest limit of\nthis analysis: the enriched sets contain "
        "risk-increasing and protective genes in equal measure, so none of them can be\n"
        "said to push risk one way.",
        transform=ax.transAxes, fontsize=7.0, color=INK2, va="top")

# ---------------------------------------------------------------- c canonical families
ax = fig.add_subplot(gs[2, 0]); tag(ax, "c", dx=-0.32, dy=1.06)
FSp = FS.sort_values("p")
yy = np.arange(len(FSp))[::-1]
ax.barh(yy, FSp.n_up, height=0.60, color=UPC, zorder=3, label="risk-increasing")
ax.barh(yy, FSp.n_genes - FSp.n_up, left=FSp.n_up, height=0.60, color=DOWNC, zorder=3,
        label="protective")
exp = FSp.family_size * 163 / 8968
ax.scatter(exp, yy, marker="|", s=90, color=INK, zorder=5, linewidths=1.6,
           label="expected by chance")
for yi, r, e_ in zip(yy, FSp.itertuples(index=False), exp):
    ptxt = f"P = {r.p:.1e}" if r.p < 0.001 else f"P = {r.p:.3f}"
    ax.text(max(r.n_genes, e_) + 0.6, yi, ptxt + ("  •" if r.p < 0.05 else ""),
            va="center", fontsize=6.8, color=GREEN if r.p < 0.05 else INK2)
ax.set_yticks(yy); ax.set_yticklabels(FSp.family, fontsize=7.4)
ax.set_xlabel("Risk genes that are members of the family")
ax.set_xlim(0, 27)
ax.set_title("Risk genes inside each canonical CRC family, enriched or not", pad=10)
ax.legend(loc="upper left", bbox_to_anchor=(0.0, -0.215), ncol=3, fontsize=7.0)
ax.grid(axis="x"); ax.set_axisbelow(True)
ax.text(0.0, -0.375,
        "Membership, not enrichment: this asks how many of the 347 are annotated to the "
        "family at all, against how many would be\nthere by chance. The germline-predisposition "
        "pathways come first — TGF-β/BMP and DNA repair — and the somatic\n"
        "driver pathways of established tumours (RTK-MAPK, PI3K, NOTCH) do not separate "
        "from chance.",
        transform=ax.transAxes, fontsize=7.0, color=INK2, va="top")
fig.savefig("figures/fig35_all347_pathways.png", dpi=260, bbox_inches="tight")
fig.savefig("figures/fig35_all347_pathways.pdf", bbox_inches="tight")
plt.close(fig)

# ==================================================================== FIGURE 36
fig = plt.figure(figsize=(7.9, 9.4))
gs = fig.add_gridspec(3, 1, height_ratios=[1.0, 0.52, 1.28], hspace=0.95)

# ---------------------------------------------------------------- a themes
ax = fig.add_subplot(gs[0, 0]); tag(ax, "a", dx=-0.34, dy=1.06)
T = TH.sort_values("fold", ascending=False)
yy = np.arange(len(T))[::-1]
ax.barh(yy, T.fold, height=0.60,
        color=[GREEN if f < 0.1 else "#b9cad8" for f in T.fdr], zorder=3)
ax.axvline(1, color="#a9a7a2", lw=1.0, zorder=2)
for yi, r in zip(yy, T.itertuples(index=False)):
    ax.text(r.fold + 0.12, yi,
            f"{int(r.n_risk_genes)} of {int(r.theme_size)} genes"
            f"   FDR = {r.fdr:.3f}" if r.n_risk_genes else "no risk genes",
            va="center", fontsize=6.8, color=INK2)
ax.set_yticks(yy); ax.set_yticklabels(T.theme, fontsize=7.2)
ax.set_xlabel("Risk genes observed / expected")
ax.set_xlim(0, 13.5)
ax.set_title("Colorectal aetiology, as themes a pathway database is not organised around",
             pad=10)
ax.grid(axis="x"); ax.set_axisbelow(True)
ax.text(0.0, -0.235,
        "Each theme merges the gene sets whose names match a fixed pattern, so the grouping "
        "is reproducible and can be argued with.\nGreen bars clear FDR < 0.1 across the ten "
        "themes. One-carbon and folate metabolism is the informative zero: 108 genes in the\n"
        "universe and not one risk gene among them, in a cancer where folate is one of the "
        "most studied exposures.",
        transform=ax.transAxes, fontsize=7.0, color=INK2, va="top")

# ---------------------------------------------------------------- b coverage
ax = fig.add_subplot(gs[1, 0]); tag(ax, "b", dx=-0.34, dy=1.16)
order = ["in a canonical CRC family", "aetiology theme only", "annotated, no theme",
         "network-adjacent only", "unannotated"]
cnt = [int(J["status"].get(k, 0)) for k in order]
cols = ["#2e8b6d", "#6a4c93", "#7d99b8", "#c07a2c", "#d8d7d2"]
left = 0
for c, col, lab in zip(cnt, cols, order):
    ax.barh([0], [c], left=left, height=0.5, color=col, zorder=3)
    if c > 20:
        ax.text(left + c / 2, 0, str(c), ha="center", va="center", fontsize=7.6,
                color="white", fontweight="bold")
    left += c
ax.set_yticks([]); ax.set_xlim(0, sum(cnt))
ax.set_xlabel("Risk genes")
ax.set_title("What the pathway view can see of the 347", pad=10)
for s_ in ["left", "right", "top"]:
    ax.spines[s_].set_visible(False)
ax.legend(handles=[Patch(color=c, label=l) for c, l in zip(cols, order)],
          loc="upper left", bbox_to_anchor=(0.0, -0.42), ncol=3, fontsize=7.0)
ax.text(0.0, -0.95,
        "Only 163 of the 347 risk genes appear in any gene set of the collection at all, and "
        "64 have no Entrez identifier to map with.\nEvery enrichment in this report is a "
        "statement about that annotated minority; the 182 unannotated genes are not evidence\n"
        "against the pathways named here, they are simply invisible to the method.",
        transform=ax.transAxes, fontsize=7.0, color=INK2, va="top")

# ---------------------------------------------------------------- c the genes
ax = fig.add_subplot(gs[2, 0]); tag(ax, "c", dx=-0.34, dy=1.045)
SH = TH[TH.n_risk_genes > 0].sort_values("p").head(6)
val = set(AS.loc[AS.validated, "gene"])
scr = set(AS.loc[AS.screen_called, "gene"])
dirn = dict(zip(AS.gene, AS.direction))
PER = 5
yrow = 0.0
yticks, ylabels = [], []
for r in SH.itertuples(index=False):
    genes = sorted(str(r.genes).split("|"))
    rows = [genes[i:i + PER] for i in range(0, len(genes), PER)]
    ymid = yrow + (len(rows) - 1) / 2
    for rr in rows:
        for j, g in enumerate(rr):
            up = dirn.get(g) == "risk-increasing"
            ax.scatter([j], [-yrow], s=30, color=UPC if up else DOWNC,
                       edgecolor="black" if g in val else "white",
                       linewidth=0.9 if g in val else 0.4, zorder=3)
            ax.text(j + 0.09, -yrow, g, ha="left", va="center", fontsize=6.0,
                    color=INK if g in val else INK2,
                    fontweight="bold" if g in scr else "normal")
        yrow += 1
    yticks.append(-ymid)
    ylabels.append(f"{r.theme}\n{int(r.n_risk_genes)} genes, FDR = {r.fdr:.3f}")
    yrow += 0.7
ax.set_yticks(yticks); ax.set_yticklabels(ylabels, fontsize=6.8)
ax.set_xlim(-0.3, PER + 0.05); ax.set_ylim(-yrow + 0.4, 0.8)
ax.set_xticks([])
for s_ in ax.spines.values():
    s_.set_visible(False)
ax.tick_params(axis="y", length=0)
ax.set_title("The risk genes behind each theme", pad=10)
ax.legend(handles=[Line2D([], [], marker="o", ls="", color=UPC, markersize=5,
                          label="risk-increasing"),
                   Line2D([], [], marker="o", ls="", color=DOWNC, markersize=5,
                          label="protective"),
                   Line2D([], [], marker="o", ls="", color="#b9b7b1",
                          markeredgecolor="black", markersize=5,
                          label="spatially validated (\u00a72.14)"),
                   Line2D([], [], marker="$\\bf{G}$", ls="", color=INK, markersize=6,
                          label="bold = viability effect in a screen")],
          loc="upper left", bbox_to_anchor=(0.0, -0.01), ncol=4, fontsize=6.6,
          handletextpad=0.3, columnspacing=1.0)
fig.savefig("figures/fig36_etiology_themes.png", dpi=260, bbox_inches="tight")
fig.savefig("figures/fig36_etiology_themes.pdf", bbox_inches="tight")
plt.close(fig)
print("fig35 + fig36 written")
