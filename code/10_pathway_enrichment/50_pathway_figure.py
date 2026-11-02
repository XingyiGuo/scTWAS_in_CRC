#!/usr/bin/env python3
"""
50 - Figure 31: what the cell-type risk genes are enriched in.

Four panels, each answering one part of the question:

  a  which pathways are enriched, in which cell type, and are they canonical CRC routes
     or something else. Rows are merged where several gene sets are carried by exactly
     the same risk genes - Reactome nests its sets, so eight "different" antigen
     presentation pathways can be three HLA genes counted eight times, and a figure that
     hides that overstates the result.
  b  the control: observed significant pathways against 1,000 random modules of the same
     size run through the identical pipeline.
  c  how much of the canonical CRC map is covered, family by family.
  d  which risk genes actually carry the enrichment, and how concentrated it is.

Output: figures/fig31_pathway_celltype.{png,pdf}
"""
import pandas as pd, numpy as np, json, os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

os.makedirs("figures", exist_ok=True)
INK, INK2, INK3 = "#0b0b0b", "#52514e", "#8a8880"
SURF = "#fcfcfb"
GREEN, GREY, AMBER = "#1baf7a", "#c9c7c1", "#e8a33d"
EST_C, NOV_C = "#2f6f9f", "#b0651f"
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
CTLAB = {"ABS": "Absorptive", "ASC": "Absorptive stem-like", "CT": "Crypt top",
         "EE": "Enteroendocrine", "GOB": "Goblet", "SSC": "Secretory stem-like",
         "STM": "Stem", "TAC": "Transit-amplifying", "TUF": "Tuft"}
FAM_C = {"WNT": "#1f6f8b", "TGF-β / BMP": "#2e8b6d", "RTK - MAPK": "#7d5ba6",
         "PI3K - AKT": "#a8577e", "p53": "#b0455a", "Mismatch repair": "#4a7ab8",
         "NOTCH": "#8a7a3c", "Cell cycle": "#356fa8", "Apoptosis": "#a2683a",
         "Chromatin remodelling": "#6a6a6a", "novel": "#c9976a"}

S = pd.read_csv("results/pathway_celltype_significant.csv")
N = pd.read_csv("results/pathway_celltype_null.csv")
G = pd.read_csv("results/pathway_leading_genes.csv")
J = json.load(open("results/pathway_celltype.json"))
FAMS = ["WNT", "TGF-β / BMP", "RTK - MAPK", "PI3K - AKT", "p53", "Mismatch repair",
        "NOTCH", "Cell cycle", "Apoptosis", "Chromatin remodelling"]

S["ct"] = S.tissue + " " + S.cell_type
order = [f"{t} {c}" for t in ["Normal", "Adenoma"]
         for c in ["STM", "TAC", "ASC", "SSC", "ABS", "CT", "GOB", "EE", "TUF"]]
cols = [c for c in order if c in set(S.ct)]

# ---------------------------------------------------------------- merge nested sets
u = (S.groupby("pathway")
     .agg(core=("genes", lambda s: "|".join(sorted({g for x in s for g in x.split("|")}))),
          family=("family", "first"), best_p=("p", "min"))
     .reset_index())
blocks = (u.groupby(["core", "family"])
          .agg(n_sets=("pathway", "size"), best_p=("best_p", "min"),
               label=("pathway", lambda s: min(s, key=len)),
               members=("pathway", lambda s: sorted(s)))
          .reset_index().sort_values("best_p"))
ORDER_KEY = {c: i for i, c in enumerate(cols)}
pw2block = {}
for bi, r in enumerate(blocks.itertuples(index=False)):
    for m in r.members:
        pw2block[m] = bi
S["block"] = S.pathway.map(pw2block)
primary = (S.sort_values("p").groupby("block").first()[["ct", "p"]])
blocks["primary"] = [primary.loc[b, "ct"] for b in range(len(blocks))]
blocks["pkey"] = blocks.primary.map(ORDER_KEY)
blocks = blocks.sort_values(["pkey", "best_p"]).reset_index(drop=True)
pw2block = {m: bi for bi, r in enumerate(blocks.itertuples(index=False)) for m in r.members}
S["block"] = S.pathway.map(pw2block)

# per block × cell type: best FDR and the largest overlap
cell = {}
for (bi, ct), g in S.groupby(["block", "ct"]):
    cell[(bi, ct)] = (float(g.fdr_celltype.min()), int(g.overlap.max()))

nb = len(blocks)
fig = plt.figure(figsize=(7.9, 0.235 * nb + 2.5))
gs = fig.add_gridspec(1, 1)


def tag(ax, s, dx=-0.10, dy=1.05):
    ax.text(dx, dy, s, transform=ax.transAxes, fontsize=11.5, fontweight="bold",
            color=INK, va="bottom", ha="left")


# ============================================================== a  dot matrix
ax = fig.add_subplot(gs[0, 0])
y = np.arange(nb)[::-1]
x = np.arange(len(cols))
for xi in x:
    ax.axvline(xi, color="#f0eee9", lw=9, zorder=0)
for yi, r in zip(y, blocks.itertuples(index=False)):
    ax.axhline(yi, color="#f4f2ed", lw=0.7, zorder=0)
sizes, colours, xs, ys = [], [], [], []
for bi, yi in zip(range(nb), y):
    for xi, ct in enumerate(cols):
        v = cell.get((bi, ct))
        if v is None:
            continue
        xs.append(xi); ys.append(yi)
        sizes.append(11 + 22 * (v[1] - 2))
        colours.append(-np.log10(max(v[0], 1e-6)))
sc = ax.scatter(xs, ys, s=sizes, c=colours, cmap="YlGnBu", vmin=1.0, vmax=3.2,
                edgecolor="white", linewidth=0.5, zorder=3)
ax.xaxis.set_ticks_position("top")
ax.set_xticks(x)
SURV = {f"{r.tissue} {r.cell_type}": r.p_perm_matched_fdr
        for r in N[N.tested].itertuples(index=False)}
ax.set_xticklabels([f"{c.split()[0]} · {CTLAB.get(c.split()[1], c.split()[1])}"
                    + ("" if SURV.get(c, 1) < 0.05 else " †")
                    for c in cols], fontsize=6.0, rotation=55, ha="left",
                   rotation_mode="anchor")
labs = [(r.label if len(r.label) <= 40 else r.label[:38] + "…") +
        ("" if r.n_sets == 1 else f"  (+{r.n_sets - 1})")
        for r in blocks.itertuples(index=False)]
ax.set_yticks(y); ax.set_yticklabels(labs, fontsize=5.9)
for yi, r in zip(y, blocks.itertuples(index=False)):
    ax.text(len(cols) - 0.30, yi, r.core.replace("|", ", "), fontsize=5.4,
            color=INK2, va="center", ha="left")
    ax.add_patch(plt.Rectangle((-0.78, yi - 0.36), 0.20, 0.72,
                               color=FAM_C.get(r.family, GREY), zorder=4, clip_on=False))
ax.set_xlim(-0.85, len(cols) + 5.6)
ax.set_ylim(-0.9, nb - 0.1)
ax.tick_params(axis="y", length=0); ax.tick_params(axis="x", length=0)
for s_ in ax.spines.values():
    s_.set_visible(False)
tag(ax, "", dx=-0.40, dy=1.20)
ax.text(-0.40, 1.163, "Enriched pathways by cell type, with the risk genes that carry them",
        transform=ax.transAxes, fontsize=8.6, fontweight="bold", color=INK,
        va="bottom", ha="left")
cax = ax.inset_axes([0.70, 1.170, 0.26, 0.006])
cb = fig.colorbar(sc, cax=cax, orientation="horizontal")
cb.set_label("−log$_{10}$ FDR within cell type  (capped at 3.2)", fontsize=5.8, labelpad=2)
cb.ax.tick_params(labelsize=5.4, length=2, pad=1); cb.outline.set_visible(False)
cax.xaxis.set_label_position("top"); cax.xaxis.set_ticks_position("bottom")
h = [Line2D([], [], marker="o", ls="", color="#7fb0c8",
            markersize=np.sqrt(11 + 22 * k) / 1.3, label=f"{k + 2} risk genes")
     for k in (0, 1, 2)]
ax.legend(handles=h, loc="upper left", bbox_to_anchor=(0.0, -0.022), ncol=3,
          fontsize=5.8, handletextpad=0.3, columnspacing=1.2)
ax.text(0.0, -0.050,
        "Left bar: canonical CRC family (blue-green) or none of them (tan). Rows carrying "
        "identical gene sets are merged and the\nbracket counts the duplicates; rows are "
        "grouped by the cell type that carries them. † marks a cell type whose count does\n"
        "not survive the annotation-matched null of Figure 32, so its rows are drawn but "
        "should not be read as evidence.",
        transform=ax.transAxes, fontsize=5.8, color=INK2, va="top")

fig.savefig("figures/fig31_pathway_celltype.png", dpi=260, bbox_inches="tight")
fig.savefig("figures/fig31_pathway_celltype.pdf", bbox_inches="tight")
plt.close(fig)

# ==================================================================== FIGURE 32
fig = plt.figure(figsize=(7.9, 8.2))
gs = fig.add_gridspec(3, 1, height_ratios=[1.20, 0.66, 1.08], hspace=0.78)

# ============================================================== b  null control
ax = fig.add_subplot(gs[0, 0]); tag(ax, "a", dx=-0.30, dy=1.06)
T = N[N.tested].copy()
T["ct"] = T.tissue + " " + T.cell_type
T = T.set_index("ct").reindex([c for c in order if c in set(T.ct)]).reset_index()
yy = np.arange(len(T))[::-1]
bars = ax.barh(yy, T.observed_fdr, height=0.60, zorder=3,
               color=[EST_C if p_ < 0.05 else "#b9cad8"
                      for p_ in T.p_perm_matched_fdr], label="observed (FDR < 0.1)")
ax.errorbar(T.null_fdr_mean, yy + 0.16, xerr=T.null_fdr_sd, fmt="o", color="#9a9893",
            ms=3.4, lw=0.9, capsize=2.2, zorder=4, label="null: size-matched")
ax.errorbar(T.matched_fdr_mean, yy - 0.16, xerr=T.matched_fdr_sd, fmt="D", color="#7a3b8f",
            ms=3.4, lw=0.9, capsize=2.2, zorder=4,
            label="null: matched on annotation depth")
for yi, r in zip(yy, T.itertuples(index=False)):
    if r.observed_fdr > 0:
        xmax = max(r.observed_fdr, r.matched_fdr_mean + r.matched_fdr_sd,
                   r.null_fdr_mean + r.null_fdr_sd)
        ax.text(xmax + 0.6, yi, f"P = {r.p_perm_matched_fdr:.3f}", va="center",
                fontsize=6.8, color=GREEN if r.p_perm_matched_fdr < 0.05 else INK2)
ax.set_yticks(yy)
ax.set_yticklabels([f"{r.tissue} {CTLAB.get(r.cell_type, r.cell_type)}  (n = {r.n_mapped})"
                    for r in T.itertuples(index=False)], fontsize=7.0)
ax.set_xlabel("Significant pathways")
ax.set_xlim(0, max(T.observed_fdr.max(),
                   (T.matched_fdr_mean + T.matched_fdr_sd).max()) * 1.36)
ax.set_title("Against two nulls: the second one is the honest one", pad=12)
ax.legend(loc="upper left", bbox_to_anchor=(0.0, -0.215), ncol=3, fontsize=7.0,
          columnspacing=1.6, handletextpad=0.4)
ax.grid(axis="x"); ax.set_axisbelow(True)
ax.text(0.0, -0.335,
        "P is against the annotation-matched null. Adenoma goblet is the instructive failure: "
        "its 17 hits are fewer\nthan genes of the same annotation depth produce by chance, so "
        "that block of WNT and cell-cycle sets\nis not evidence, and the same test is what "
        "keeps tuft, crypt-top and transit-amplifying standing.",
        transform=ax.transAxes, fontsize=7.0, color=INK2, va="top", ha="left")

# ============================================================== c  canonical coverage
ax = fig.add_subplot(gs[1, 0]); tag(ax, "b", dx=-0.30, dy=1.10)
est = S[S.established]
cov = [(f, int(est[est.family == f].pathway.nunique()),
        sorted({g for x in est[est.family == f].genes for g in x.split("|")}))
       for f in FAMS]
yy = np.arange(len(cov))[::-1]
ax.barh(yy, [c[1] for c in cov], height=0.62,
        color=[FAM_C[c[0]] if c[1] else "#e8e6e1" for c in cov], zorder=3)
for yi, (f, n_, gg) in zip(yy, cov):
    if n_:
        ax.text(n_ + 0.2, yi, ", ".join(gg), va="center", fontsize=6.6, color=INK2)
ax.set_yticks(yy); ax.set_yticklabels([c[0] for c in cov], fontsize=7.4)
ax.set_xlabel("Enriched gene sets (FDR < 0.1 in any cell type)")
ax.set_xlim(0, max(c[1] for c in cov) * 2.9)
ax.set_title("The canonical CRC map, family by family")
ax.grid(axis="x"); ax.set_axisbelow(True)

# ============================================================== d  leading genes
ax = fig.add_subplot(gs[2, 0]); tag(ax, "c", dx=-0.30, dy=1.06)
g = G.head(16).iloc[::-1]
yy = np.arange(len(g))
ax.barh(yy, g.n_established, height=0.66, color=EST_C, zorder=3,
        label="canonical CRC family")
ax.barh(yy, g.n_novel, left=g.n_established, height=0.66, color=NOV_C, zorder=3,
        label="outside those families")
ax.set_yticks(yy)
ax.set_yticklabels([f"{r.gene}{' ▲' if r.validated else ''}" for r in g.itertuples(index=False)],
                   fontsize=7.2)
ax.set_xlabel("Enriched gene sets the gene appears in")
ax.set_title("Which risk genes carry the enrichment")
ax.legend(loc="lower right", fontsize=7.0)
ax.grid(axis="x"); ax.set_axisbelow(True)
ax.text(0.0, -0.30,
        "▲ spatially validated (§2.14). Three HLA genes account for the antigen-presentation\n"
        "block, so the count of enriched gene sets is not a count of independent findings.",
        transform=ax.transAxes, fontsize=7.2, color=INK2, va="top")

fig.savefig("figures/fig32_pathway_controls.png", dpi=260, bbox_inches="tight")
fig.savefig("figures/fig32_pathway_controls.pdf", bbox_inches="tight")
plt.close(fig)
print("fig31 +", "fig32 written;", nb, "merged rows,", len(cols), "cell types")
