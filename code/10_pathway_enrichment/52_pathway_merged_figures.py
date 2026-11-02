#!/usr/bin/env python3
"""
52 - Figures 33 and 34: the pooled-by-cell-type enrichment and what pooling did.

Figure 33 is the same dot matrix as Figure 31, drawn on the nine pooled modules.
Figure 34 is the accounting: what each pooled module is made of, how it scores against
the two nulls, and what pooling changed relative to the split analysis - including the
part the request was really about, whether the canonical CRC families come out better.
"""
import pandas as pd, numpy as np, json, os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

os.makedirs("figures", exist_ok=True)
INK, INK2, INK3 = "#0b0b0b", "#52514e", "#8a8880"
SURF = "#fcfcfb"
GREEN, GREY, AMBER = "#1baf7a", "#c9c7c1", "#e8a33d"
EST_C, NOV_C = "#2f6f9f", "#b0651f"
NORM_C, ADEN_C, SHARE_C = "#4a7ab8", "#c98a3c", "#6a9f7a"
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
FAM_C = {"WNT": "#1f6f8b", "TGF-β / BMP": "#2e8b6d", "RTK - MAPK": "#7d5ba6",
         "PI3K - AKT": "#a8577e", "p53": "#b0455a", "Mismatch repair": "#4a7ab8",
         "NOTCH": "#8a7a3c", "Cell cycle": "#356fa8", "Apoptosis": "#a2683a",
         "Chromatin remodelling": "#6a6a6a", "novel": "#c9976a"}
FAMS = ["WNT", "TGF-β / BMP", "RTK - MAPK", "PI3K - AKT", "p53", "Mismatch repair",
        "NOTCH", "Cell cycle", "Apoptosis", "Chromatin remodelling"]

S = pd.read_csv("results/pathway_merged_significant.csv")
N = pd.read_csv("results/pathway_merged_null.csv")
C = pd.read_csv("results/pathway_merged_vs_split.csv")
F = pd.read_csv("results/pathway_merged_families.csv")
SP = pd.read_csv("results/pathway_celltype_significant.csv")
J = json.load(open("results/pathway_merged.json"))

ORDER = ["stem", "transit-amplifying", "absorptive stem-like", "secretory stem-like",
         "absorptive", "crypt top", "goblet", "enteroendocrine", "tuft"]
cols = [c for c in ORDER if c in set(S.label)]
SURV = {r.label: r.p_perm_matched_fdr for r in N[N.tested].itertuples(index=False)}


def tag(ax, s, dx=-0.10, dy=1.05):
    ax.text(dx, dy, s, transform=ax.transAxes, fontsize=11.5, fontweight="bold",
            color=INK, va="bottom", ha="left")


# ==================================================================== FIGURE 33
u = (S.groupby("pathway")
     .agg(core=("genes", lambda s: "|".join(sorted({g for x in s for g in x.split("|")}))),
          family=("family", "first"), best_p=("p", "min"))
     .reset_index())
blocks = (u.groupby(["core", "family"])
          .agg(n_sets=("pathway", "size"), best_p=("best_p", "min"),
               label=("pathway", lambda s: min(s, key=len)),
               members=("pathway", lambda s: sorted(s)))
          .reset_index())
KEY = {c: i for i, c in enumerate(cols)}
pw2b = {m: bi for bi, r in enumerate(blocks.itertuples(index=False)) for m in r.members}
S["block"] = S.pathway.map(pw2b)
primary = S.sort_values("p").groupby("block").first()[["label", "p"]]
blocks["primary"] = [primary.loc[b, "label"] for b in range(len(blocks))]
blocks["pkey"] = blocks.primary.map(KEY)
blocks = blocks.sort_values(["pkey", "best_p"]).reset_index(drop=True)
pw2b = {m: bi for bi, r in enumerate(blocks.itertuples(index=False)) for m in r.members}
S["block"] = S.pathway.map(pw2b)
NEW = set(S.pathway) - set(SP.pathway)
cell, isnew = {}, {}
for (bi, ct), g in S.groupby(["block", "label"]):
    cell[(bi, ct)] = (float(g.fdr_merged.min()), int(g.overlap.max()),
                      str(g.tissues_contributing.iloc[0]))
for bi, r in enumerate(blocks.itertuples(index=False)):
    isnew[bi] = all(m in NEW for m in r.members)

nb = len(blocks)
fig = plt.figure(figsize=(7.9, 0.235 * nb + 2.5))
ax = fig.add_subplot(111)
y = np.arange(nb)[::-1]
x = np.arange(len(cols))
for xi in x:
    ax.axvline(xi, color="#f0eee9", lw=9, zorder=0)
for yi in y:
    ax.axhline(yi, color="#f4f2ed", lw=0.7, zorder=0)
xs, ys, sizes, colours, edges_ = [], [], [], [], []
for bi, yi in zip(range(nb), y):
    for xi, ct in enumerate(cols):
        v = cell.get((bi, ct))
        if v is None:
            continue
        xs.append(xi); ys.append(yi)
        sizes.append(11 + 22 * (v[1] - 2))
        colours.append(-np.log10(max(v[0], 1e-6)))
        edges_.append("#1b1b1b" if v[2] == "normal+adenoma" else "white")
sc = ax.scatter(xs, ys, s=sizes, c=colours, cmap="YlGnBu", vmin=1.0, vmax=3.2,
                edgecolor=edges_, linewidth=0.7, zorder=3)
ax.xaxis.set_ticks_position("top")
ax.set_xticks(x)
ax.set_xticklabels([c[0].upper() + c[1:] + ("" if SURV.get(c, 1) < 0.05 else " †")
                    for c in cols], fontsize=6.4, rotation=55, ha="left",
                   rotation_mode="anchor")
labs = [(r.label if len(r.label) <= 40 else r.label[:38] + "…") +
        ("" if r.n_sets == 1 else f"  (+{r.n_sets - 1})") + ("  ◆" if isnew[i] else "")
        for i, r in enumerate(blocks.itertuples(index=False))]
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
ax.text(-0.40, 1.163,
        "Enriched pathways with normal and adenoma pooled in one module per cell type",
        transform=ax.transAxes, fontsize=8.6, fontweight="bold", color=INK,
        va="bottom", ha="left")
cax = ax.inset_axes([0.70, 1.170, 0.26, 0.006])
cb = fig.colorbar(sc, cax=cax, orientation="horizontal")
cb.set_label("−log$_{10}$ FDR within the pooled module  (capped at 3.2)", fontsize=5.8,
             labelpad=2)
cb.ax.tick_params(labelsize=5.4, length=2, pad=1); cb.outline.set_visible(False)
cax.xaxis.set_label_position("top"); cax.xaxis.set_ticks_position("bottom")
h = [Line2D([], [], marker="o", ls="", color="#7fb0c8",
            markersize=np.sqrt(11 + 22 * k) / 1.3, label=f"{k + 2} risk genes")
     for k in (0, 1, 2)]
h.append(Line2D([], [], marker="o", ls="", color="#7fb0c8", markeredgecolor="#1b1b1b",
                markersize=3.4, label="genes from both tissues"))
ax.legend(handles=h, loc="upper left", bbox_to_anchor=(0.0, -0.022), ncol=4,
          fontsize=5.8, handletextpad=0.3, columnspacing=1.0)
ax.text(0.0, -0.050,
        "◆ marks a pathway that pooling brought out — it is not enriched in either tissue's "
        "module on its own. A dark-rimmed dot means the\\nenrichment draws risk genes from "
        "both the normal and the adenoma list, which is the only way pooling can help. Left "
        "bar: canonical\\nCRC family or none of them. † marks a module whose count does not "
        "survive the annotation-matched null of Figure 34b.".replace("\\n", "\n"),
        transform=ax.transAxes, fontsize=5.8, color=INK2, va="top")
fig.savefig("figures/fig33_pathway_merged.png", dpi=260, bbox_inches="tight")
fig.savefig("figures/fig33_pathway_merged.pdf", bbox_inches="tight")
plt.close(fig)

# ==================================================================== FIGURE 34
fig = plt.figure(figsize=(7.9, 9.4))
gs = fig.add_gridspec(3, 1, height_ratios=[0.92, 1.05, 0.95], hspace=1.02)

# ---------------------------------------------------------------- a  composition
ax = fig.add_subplot(gs[0, 0]); tag(ax, "a", dx=-0.30, dy=1.06)
T = N.set_index("label").reindex([c for c in ORDER if c in set(N.label)]).reset_index()
yy = np.arange(len(T))[::-1]
only_n = T.n_normal - T.n_shared
only_a = T.n_adenoma - T.n_shared
ax.barh(yy, only_n, height=0.62, color=NORM_C, zorder=3, label="normal only")
ax.barh(yy, T.n_shared, left=only_n, height=0.62, color=SHARE_C, zorder=3,
        label="in both lists")
ax.barh(yy, only_a, left=only_n + T.n_shared, height=0.62, color=ADEN_C, zorder=3,
        label="adenoma only")
for yi, r in zip(yy, T.itertuples(index=False)):
    ax.text(r.n_mapped + 1.0, yi, f"{int(r.n_mapped)}", va="center", fontsize=7.0,
            color=INK2)
ax.set_yticks(yy); ax.set_yticklabels(T.label, fontsize=7.4)
ax.set_xlabel("Risk genes in the pathway universe")
ax.set_xlim(0, T.n_mapped.max() * 1.18)
ax.set_title("What pooling actually adds", pad=10)
ax.legend(loc="upper left", bbox_to_anchor=(0.0, -0.235), ncol=3, fontsize=7.0,
          columnspacing=1.8, handletextpad=0.4)
ax.grid(axis="x"); ax.set_axisbelow(True)
ax.text(0.0, -0.40,
        "The two lists barely overlap, so pooling is close to a union of two independent "
        "gene lists rather than a re-count of one.\nIt also raises the bar: a bigger module "
        "needs a bigger overlap to reach the same hypergeometric p, so pooling only pays "
        "when\nboth tissues put genes in the same gene set.",
        transform=ax.transAxes, fontsize=7.0, color=INK2, va="top")

# ---------------------------------------------------------------- b  nulls
ax = fig.add_subplot(gs[1, 0]); tag(ax, "b", dx=-0.30, dy=1.06)
T2 = T[T.tested].copy()
yy = np.arange(len(T2))[::-1]
ax.barh(yy, T2.observed_fdr, height=0.60, zorder=3,
        color=[EST_C if p_ < 0.05 else "#b9cad8" for p_ in T2.p_perm_matched_fdr],
        label="observed (FDR < 0.1)")
ax.errorbar(T2.null_fdr_mean, yy + 0.16, xerr=T2.null_fdr_sd, fmt="o", color="#9a9893",
            ms=3.4, lw=0.9, capsize=2.2, zorder=4, label="null: size-matched")
ax.errorbar(T2.matched_fdr_mean, yy - 0.16, xerr=T2.matched_fdr_sd, fmt="D",
            color="#7a3b8f", ms=3.4, lw=0.9, capsize=2.2, zorder=4,
            label="null: matched on annotation depth")
for yi, r in zip(yy, T2.itertuples(index=False)):
    xm = max(r.observed_fdr, r.matched_fdr_mean + r.matched_fdr_sd,
             r.null_fdr_mean + r.null_fdr_sd)
    ax.text(xm + 0.6, yi, f"P = {r.p_perm_matched_fdr:.3f}", va="center", fontsize=6.8,
            color=GREEN if r.p_perm_matched_fdr < 0.05 else INK2)
ax.set_yticks(yy)
ax.set_yticklabels([f"{r.label}  (n = {int(r.n_mapped)})" for r in T2.itertuples(index=False)],
                   fontsize=7.2)
ax.set_xlabel("Significant pathways")
ax.set_xlim(0, max(T2.observed_fdr.max(),
                   (T2.matched_fdr_mean + T2.matched_fdr_sd).max()) * 1.32)
ax.set_title("The pooled modules against the same two nulls", pad=10)
ax.legend(loc="upper left", bbox_to_anchor=(0.0, -0.215), ncol=3, fontsize=7.0,
          columnspacing=1.6, handletextpad=0.4)
ax.grid(axis="x"); ax.set_axisbelow(True)
ax.text(0.0, -0.335,
        "Pooling brings stem and enteroendocrine cells above the annotation-matched null "
        "for the first time, and it costs goblet\ncells the result they had in Figure 32: "
        "the adenoma block that carried WNT and cell cycle dissolves when its 16 genes are\n"
        "diluted into 55, and what replaces it is carried by CTNNB1, HRAS, ERBB2 and CDH1 — "
        "genes so heavily annotated that the\nmatched null expects more hits than were found.",
        transform=ax.transAxes, fontsize=7.0, color=INK2, va="top")

# ---------------------------------------------------------------- c  canonical families
ax = fig.add_subplot(gs[2, 0]); tag(ax, "c", dx=-0.30, dy=1.06)
yy = np.arange(len(FAMS))[::-1]
w = 0.36
ax.barh(yy + w / 2 + 0.02, [int(F.loc[F.family == f, "split_sets"].iat[0]) for f in FAMS],
        height=w, color="#b9b7b1", zorder=3, label="split by tissue (Figure 31)")
ax.barh(yy - w / 2 - 0.02, [int(F.loc[F.family == f, "merged_sets"].iat[0]) for f in FAMS],
        height=w, color=EST_C, zorder=3, label="pooled per cell type (Figure 33)")
for yi, f in zip(yy, FAMS):
    m = int(F.loc[F.family == f, "merged_sets"].iat[0])
    g = F.loc[F.family == f, "merged_genes"].iat[0]
    if m:
        ax.text(m + 0.25, yi - w / 2 - 0.02, str(g).replace("|", ", "), va="center",
                fontsize=6.4, color=INK2)
ax.set_yticks(yy); ax.set_yticklabels(FAMS, fontsize=7.4)
ax.set_xlabel("Enriched gene sets in that family")
ax.set_xlim(0, 12)
ax.set_title("Did pooling buy the canonical CRC pathways?", pad=10)
ax.legend(loc="lower right", fontsize=7.0)
ax.grid(axis="x"); ax.set_axisbelow(True)
ax.text(0.0, -0.285,
        "Not by count: the canonical total falls from "
        f"{int(F.split_sets.sum() - F.loc[F.family == 'novel', 'split_sets'].iat[0])} sets to "
        f"{int(F.merged_sets.sum() - F.loc[F.family == 'novel', 'merged_sets'].iat[0])}, and "
        "WNT disappears with the adenoma goblet block. What pooling adds is two\ncanonical "
        "gene sets that neither tissue reaches alone — TGF-beta receptor signalling and the "
        "ATM-dependent DNA damage\nresponse, both in goblet cells and both drawing genes from "
        "both lists — and a metabolic axis in stem cells (Figure 33).",
        transform=ax.transAxes, fontsize=7.0, color=INK2, va="top")

fig.savefig("figures/fig34_pathway_merged_controls.png", dpi=260, bbox_inches="tight")
fig.savefig("figures/fig34_pathway_merged_controls.pdf", bbox_inches="tight")
plt.close(fig)
print("fig33 + fig34 written;", nb, "blocks,", len(cols), "pooled modules")
