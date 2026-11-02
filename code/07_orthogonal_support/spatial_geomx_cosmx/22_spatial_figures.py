#!/usr/bin/env python3
"""22 - Spatial figures (12-14)."""
import pandas as pd, numpy as np, json, os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm

os.makedirs("figures", exist_ok=True)
C = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#4a3aa7", "#e34948"]
POS, NEG = "#e34948", "#2a78d6"
SEQ = ["#cde2fb", "#9ec5f4", "#6da7ec", "#3987e5", "#256abf", "#184f95", "#0d366b"]
INK, INK2, INK3 = "#0b0b0b", "#52514e", "#8a8880"
SURF = "#fcfcfb"
SEQCM = LinearSegmentedColormap.from_list("seq", SEQ)
DIV = LinearSegmentedColormap.from_list("div", ["#184f95", "#6da7ec", "#f0efec",
                                                "#f0a09f", "#e34948"])
plt.rcParams.update({
    "figure.facecolor": SURF, "axes.facecolor": SURF, "savefig.facecolor": SURF,
    "font.family": "DejaVu Sans", "font.size": 8.5,
    "axes.edgecolor": "#d8d7d2", "axes.linewidth": 0.8,
    "axes.labelcolor": INK2, "text.color": INK,
    "xtick.color": INK2, "ytick.color": INK2, "xtick.labelsize": 8, "ytick.labelsize": 8,
    "axes.titlesize": 9.5, "axes.titleweight": "bold", "axes.titlelocation": "left",
    "axes.spines.top": False, "axes.spines.right": False,
    "grid.color": "#eceae5", "grid.linewidth": 0.7, "legend.frameon": False,
    "legend.fontsize": 7.8, "pdf.fonttype": 42, "ps.fonttype": 42,
})


def tag(ax, s):
    ax.text(-0.10, 1.09, s, transform=ax.transAxes, fontsize=11, fontweight="bold",
            color=INK, va="top", ha="left")


def save(fig, name):
    fig.savefig(f"figures/{name}.png", dpi=220, bbox_inches="tight")
    fig.savefig(f"figures/{name}.pdf", bbox_inches="tight")
    plt.close(fig)


REGCOL = {"NOR": "#1baf7a", "ADE": "#eda100", "CORE": "#e34948", "INV": "#4a3aa7"}
REGNAME = {"NOR": "Normal mucosa", "ADE": "Adenoma", "CORE": "Tumour core",
           "INV": "Invasive front"}
PRETTY = {"NOR->ADE": "Normal → adenoma", "ADE->CORE": "Adenoma → carcinoma",
          "NOR->CORE": "Normal → carcinoma", "ADE->INV": "Adenoma → invasive front",
          "NOR->INV": "Normal → invasive front", "CORE->INV": "Core → invasive front"}
ALL_T = ["NOR->ADE", "ADE->CORE", "NOR->CORE", "ADE->INV", "NOR->INV", "CORE->INV"]

# ---------------------------------------------------------------- load
R = pd.read_csv("results/spatial_risk_gene_DE.csv")
Cc = pd.read_csv("results/spatial_direction_concordance.csv")
Rep = pd.read_csv("results/spatial_replication.csv")
de2 = pd.concat([pd.read_csv("spatial/de_exp2_epi.tsv.gz", sep="\t"),
                 pd.read_csv("spatial/de_exp2_str.tsv.gz", sep="\t")], ignore_index=True)
de1 = pd.read_csv("spatial/de_exp1.tsv.gz", sep="\t")
ex = pd.read_csv("spatial/geomx_expr_exp2.tsv.gz", sep="\t", index_col=0)
md = pd.read_csv("spatial/exp2_annot_Metadata.tsv", sep="\t")
sp = pd.read_csv("spatial/exp2_annot_SegmentProperties.tsv", sep="\t")
ann = (md.dropna(subset=["SegmentDisplayName"])
       .merge(sp, on="SegmentDisplayName", how="left").set_index("Sample_ID"))
ann = ann.loc[[i for i in ex.index if i in ann.index]]
EX = ex.loc[ann.index]
RISK = set(R["gene"])

# ================================================================== FIGURE 12
GENES = ["CTNNB1", "CEACAM5", "TRIM28", "UGT2B15"]
SLIDE = "slide5"
sel = ann[(ann["slide name"] == SLIDE) & (ann["segment"] == "PanCK+") &
          ann["tissueRegion"].isin(REGCOL)]
print(f"slide {SLIDE}: {len(sel)} epithelial ROIs, regions "
      f"{dict(sel.tissueRegion.value_counts())}")

fig = plt.figure(figsize=(13.8, 8.2))
gs = fig.add_gridspec(2, 3, hspace=0.42, wspace=0.34)
ymax = sel["ROICoordinateY"].max()

def roi_panel(ax, colour_by=None, gene=None):
    y = (ymax - sel["ROICoordinateY"]) / 1000
    x = sel["ROICoordinateX"] / 1000
    if gene is None:
        for reg, g in sel.groupby("tissueRegion"):
            yy = (ymax - g["ROICoordinateY"]) / 1000
            ax.scatter(g["ROICoordinateX"] / 1000, yy, s=60, color=REGCOL[reg],
                       label=REGNAME[reg], edgecolor=SURF, linewidth=0.9)
        return None
    v = EX.loc[sel.index, gene].astype(float)
    return ax.scatter(x, y, c=v, s=64, cmap=SEQCM, edgecolor=SURF, linewidth=0.9,
                      vmin=np.percentile(v, 2), vmax=np.percentile(v, 98))

ax = fig.add_subplot(gs[0, 0])
roi_panel(ax)
ax.set_aspect("equal")
ax.set_xlabel("Scan X (×10³ px)"); ax.set_ylabel("Scan Y (×10³ px)")
ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.17), ncols=2, fontsize=7.2)
ax.set_title(f"GeoMx ROIs, {SLIDE}")
ax.grid(); ax.set_axisbelow(True); tag(ax, "a")

POSN = [(0, 1), (0, 2), (1, 0), (1, 1)]
for k, gene in enumerate(GENES):
    ax = fig.add_subplot(gs[POSN[k]])
    if gene not in EX.columns:
        ax.axis("off"); continue
    sc = roi_panel(ax, gene=gene)
    ax.set_aspect("equal")
    ax.set_xlabel("Scan X (×10³ px)")
    if POSN[k][1] == 0:
        ax.set_ylabel("Scan Y (×10³ px)")
    cb = fig.colorbar(sc, ax=ax, fraction=0.046, pad=0.02)
    cb.set_label("log₂ Q3-norm.", fontsize=7); cb.outline.set_visible(False)
    row = R[(R.gene == gene) & (R.experiment == "exp2") & (R.compartment == "PanCK+") &
            (R.transition == "NOR->CORE")]
    ax.set_title(gene)
    if len(row):
        rr = row.iloc[0]
        ax.text(0.02, -0.30, f"normal → carcinoma  log₂FC {rr.log2FC:+.2f},  FDR {rr.FDR:.1e}",
                transform=ax.transAxes, fontsize=7.2, color=INK2)
    ax.grid(); ax.set_axisbelow(True); tag(ax, "bcde"[k])

ax = fig.add_subplot(gs[1, 2])
sig = R[(R.experiment == "exp2") & (R.compartment == "PanCK+") &
        (R.transition == "NOR->CORE") & (R.FDR < 0.05)]
top = sig.reindex(sig["log2FC"].abs().sort_values(ascending=False).index).head(22)
epi_all = ann[(ann["segment"] == "PanCK+") & ann["tissueRegion"].isin(REGCOL)]
M = np.array([[EX.loc[epi_all[epi_all["tissueRegion"] == reg].index, g].astype(float).mean()
               for reg in ["NOR", "ADE", "CORE", "INV"]] for g in top["gene"]])
Mz = (M - M.mean(axis=1, keepdims=True)) / (M.std(axis=1, keepdims=True) + 1e-9)
im = ax.imshow(Mz, aspect="auto", cmap=DIV, norm=TwoSlopeNorm(0, -2.2, 2.2))
ax.set_xticks(range(4))
ax.set_xticklabels(["Normal", "Adenoma", "Core", "Invasive"], rotation=32, ha="right",
                   fontsize=7.6)
ax.set_yticks(range(len(top)))
ax.set_yticklabels([f"{g}  {'↑' if d=='risk-increasing' else '↓'}"
                    for g, d in zip(top["gene"], top["twas_direction"])], fontsize=7.0)
cb = fig.colorbar(im, ax=ax, fraction=0.04, pad=0.03)
cb.set_label("Row z of region mean", fontsize=7.2); cb.outline.set_visible(False)
ax.set_title("Region means, all 9 patients")
ax.text(0.0, -0.22, "↑ risk-increasing   ↓ protective (TWAS)", transform=ax.transAxes,
        fontsize=7, color=INK2)
tag(ax, "f")
save(fig, "fig12_geomx_spatial_expression")

# ================================================================== FIGURE 13
# Only the three transitions of the adenoma-carcinoma sequence are shown; the three
# invasive-front contrasts are still computed and are in Tables S27 and S40, but they
# are a different question and clutter the panel. Each volcano highlights only the genes
# that are significant AND move the way the TWAS predicts - the genes the report goes on
# to use. Significant genes that move the other way are drawn as open rings and counted
# in the legend, so they are visible without being promoted.
MAIN_T = ["NOR->ADE", "ADE->CORE", "NOR->CORE"]
SHORT = {"NOR->ADE": "N→A", "ADE->CORE": "A→C", "NOR->CORE": "N→C"}
XLAB = {"NOR->ADE": "normal → adenoma", "ADE->CORE": "adenoma → carcinoma",
        "NOR->CORE": "normal → carcinoma"}

fig = plt.figure(figsize=(15.6, 8.6))
gs = fig.add_gridspec(2, 3, hspace=0.46, wspace=0.30)

# ---- (a) counts, three transitions only
ax = fig.add_subplot(gs[0, 0]); tag(ax, "a")
sub = R[(R.experiment == "exp2") & R.transition.isin(MAIN_T)]
cnt = (sub[sub["sig"]].groupby(["compartment", "transition"]).size()
       .unstack(fill_value=0).reindex(columns=MAIN_T).fillna(0))
x = np.arange(len(MAIN_T))
ax.bar(x - 0.19, cnt.loc["PanCK+"], width=0.38, color=C[0], label="Epithelium")
ax.bar(x + 0.19, cnt.loc["PanCK-"], width=0.38, color=C[2], label="Stroma")
for i, v in enumerate(cnt.loc["PanCK+"]):
    ax.text(i - 0.19, v + 2.5, str(int(v)), ha="center", fontsize=7.6, color=INK2)
for i, v in enumerate(cnt.loc["PanCK-"]):
    ax.text(i + 0.19, v + 2.5, str(int(v)), ha="center", fontsize=7.6, color=INK3)
ax.set_xticks(x); ax.set_xticklabels([SHORT[t] for t in MAIN_T], fontsize=9.5)
ax.set_xlabel("N normal · A adenoma · C tumour core", fontsize=7.6)
ax.set_ylabel("Risk genes, FDR < 0.05 within the tested set")
ax.set_ylim(0, max(cnt.values.max() * 1.18, 10))
ax.legend(loc="upper left")
ax.grid(axis="y"); ax.set_axisbelow(True)
ax.set_title("Risk genes that change at each step")
ax.text(0.98, 0.95, f"of {R[R.experiment == 'exp2'].gene.nunique()} measured",
        transform=ax.transAxes, ha="right", fontsize=7.3, color=INK2)

# ---- (b) direction concordance, three transitions only
ax = fig.add_subplot(gs[0, 1]); tag(ax, "b")
cc = Cc[(Cc.experiment == "exp2") & Cc.transition.isin(MAIN_T)].copy()
cc["lab"] = cc["compartment"].map({"PanCK+": "Epithelium", "PanCK-": "Stroma"}) + "  " + \
    cc["transition"].map(SHORT)
cc["ord"] = cc["transition"].map({t: i for i, t in enumerate(MAIN_T)})
cc = cc.sort_values(["compartment", "ord"])
y = np.arange(len(cc))[::-1]
se = np.sqrt(cc["frac_concordant"] * (1 - cc["frac_concordant"]) / cc["n_genes"])
cols = [C[0] if c == "PanCK+" else C[2] for c in cc["compartment"]]
ax.errorbar(cc["frac_concordant"] * 100, y, xerr=1.96 * se * 100, fmt="none",
            ecolor=INK3, elinewidth=1)
ax.scatter(cc["frac_concordant"] * 100, y, s=38, c=cols, zorder=3,
           edgecolor=SURF, linewidth=0.7)
ax.axvline(50, color=INK3, lw=1, ls=(0, (3, 3)))
ax.set_yticks(y); ax.set_yticklabels(cc["lab"], fontsize=7.8)
ax.set_xlabel("Risk genes moving in the TWAS-predicted direction (%)")
ax.set_xlim(35, 68)
ax.legend(handles=[Patch(facecolor=C[0], label="Epithelium"),
                   Patch(facecolor=C[2], label="Stroma")], loc="lower right")
ax.grid(axis="x"); ax.set_axisbelow(True)
ax.set_title("Direction concordance (95 % CI)")
ax.text(0.02, -0.19, "no transition differs from chance", transform=ax.transAxes,
        fontsize=7.3, color=INK2)

# ---- (c) replication
ax = fig.add_subplot(gs[0, 2]); tag(ax, "c")
a = de2[(de2.compartment == "PanCK+") & (de2.transition == "NOR->CORE")]
b = de1[(de1.compartment == "PanCK+") & (de1.transition == "NOR->CORE")]
m = a.merge(b, on="gene", suffixes=("_2", "_1"))
mr = m[m["gene"].isin(RISK)]
ax.scatter(m["log2FC_2"], m["log2FC_1"], s=4, color="#d8d7d2", alpha=0.55,
           edgecolor="none", rasterized=True)
ax.scatter(mr["log2FC_2"], mr["log2FC_1"], s=22, color=C[1], alpha=0.9, edgecolor="none")
lim = np.nanpercentile(np.abs(np.r_[m["log2FC_2"], m["log2FC_1"]]), 99.5)
ax.plot([-lim, lim], [-lim, lim], color=INK3, lw=1, ls=(0, (3, 3)))
ax.axhline(0, color="#e5e3dd", lw=0.8); ax.axvline(0, color="#e5e3dd", lw=0.8)
ax.set_xlim(-lim, lim); ax.set_ylim(-lim, lim)
row = Rep[(Rep.compartment == "PanCK+") & (Rep.transition == "NOR->CORE")].iloc[0]
ax.set_xlabel("log₂FC, cohort 2 (WTA, 9 patients)")
ax.set_ylabel("log₂FC, cohort 1 (CTA, 10 patients)")
ax.text(0.03, 0.95, f"risk genes ρ = {row.rho_risk:+.2f}\n"
                    f"({int(row.n_risk_shared)} genes, {100*row.sign_agreement_risk:.0f} % same sign)\n"
                    f"all genes ρ = {row.rho_all:+.2f}",
        transform=ax.transAxes, va="top", fontsize=7.4, color=INK2)
ax.grid(); ax.set_axisbelow(True)
ax.set_title("Replication, normal → carcinoma")

# ---- (d-f) one volcano per transition, supporting genes highlighted
for k, t in enumerate(MAIN_T):
    ax = fig.add_subplot(gs[1, k]); tag(ax, "def"[k])
    rr = R[(R.experiment == "exp2") & (R.compartment == "PanCK+")
           & (R.transition == t)].copy()
    rr["nlp"] = -np.log10(rr["FDR_risk"].clip(lower=1e-30))
    ns = rr[~rr["sig"]]
    sup = rr[rr["sig"] & rr["concordant"]]
    opp = rr[rr["sig"] & ~rr["concordant"]]
    ax.scatter(ns["log2FC"], ns["nlp"], s=10, color="#dedcd6", edgecolor="none",
               label=f"not significant ({len(ns)})")
    ax.scatter(opp["log2FC"], opp["nlp"], s=17, facecolor="none", edgecolor="#a9a7a0",
               linewidth=0.8, label=f"significant, opposite direction ({len(opp)})")
    ax.scatter(sup["log2FC"], sup["nlp"], s=24, color=POS, edgecolor="none",
               label=f"supports TWAS direction ({len(sup)})")
    ax.axhline(-np.log10(0.05), color=INK3, lw=0.8, ls=(0, (3, 3)))
    ax.axvline(0, color="#e0dfda", lw=0.8)
    ymax = max(rr["nlp"].max() * 1.22, 2.0)
    ax.set_ylim(-0.25, ymax)
    n_lab = min(6, int(rr["sig"].sum() and (rr["sig"] & rr["concordant"]).sum()))
    gutter = 0.30 if n_lab <= 3 else 1.55
    xlo = rr["log2FC"].min() - 0.35
    xhi = rr["log2FC"].max() + gutter
    ax.set_xlim(xlo, xhi)
    if len(sup):
        lab = sup.reindex(sup["log2FC"].abs().sort_values(ascending=False).index).head(6)
        if len(lab) <= 3:
            # few enough to sit beside their own points
            for i, r_ in enumerate(lab.itertuples(index=False)):
                dx = 9 if r_.log2FC >= 0 else -9
                ax.annotate(r_.gene, (r_.log2FC, r_.nlp), textcoords="offset points",
                            xytext=(dx, 7 if i % 2 == 0 else -11), fontsize=7,
                            color=INK2, ha="left" if dx > 0 else "right", va="center")
        else:
            lab = lab.sort_values("nlp", ascending=False)
            ys = np.linspace(ymax * 0.94, ymax * 0.44, len(lab))
            for (yy, r_) in zip(ys, lab.itertuples(index=False)):
                ax.annotate(r_.gene, (r_.log2FC, r_.nlp), xytext=(xhi - 0.12, yy),
                            textcoords="data", fontsize=6.9, color=INK2, ha="right",
                            va="center", arrowprops=dict(arrowstyle="-", lw=0.5,
                                                         color=INK3, shrinkA=1, shrinkB=3))
    ax.set_xlabel(f"log₂ fold change, {XLAB[t]} (epithelium)")
    if k == 0:
        ax.set_ylabel("−log₁₀ FDR (within the tested risk genes)")
    ax.legend(loc="upper left", fontsize=6.9)
    ax.grid(); ax.set_axisbelow(True)
    ax.set_title(f"{XLAB[t][0].upper()}{XLAB[t][1:]}, {len(rr)} tested risk genes")
save(fig, "fig13_spatial_stage_de")

# ================================================================== FIGURE 14
cells = pd.read_csv("spatial/cosmx_cells.tsv.gz", sep="\t")
CT = pd.read_csv("results/cosmx_celltype_trend.csv")
prop = pd.read_csv("results/cosmx_celltype_proportions.csv")

GROUP = {}
for ct in cells["cellType"].unique():
    s = str(ct).lower()
    if ct == "CRC":
        GROUP[ct] = "Carcinoma"
    elif "oncofetal" in s:
        GROUP[ct] = "Oncofetal"
    elif "epitheli" in s:
        GROUP[ct] = "Normal epithelium"
    elif "caf" in s or "fibro" in s or "trophocyte" in s or "telocyte" in s or "stromal" in s:
        GROUP[ct] = "Fibroblast / CAF"
    elif s.startswith("t.") or "cell" in s or "plasma" in s or "macroph" in s or "mast" in s or "lymph" in s:
        GROUP[ct] = "Immune"
    elif "endothel" in s or "pericyte" in s or "smc" in s:
        GROUP[ct] = "Vascular"
    else:
        GROUP[ct] = "Other"
cells["group"] = cells["cellType"].map(GROUP)
GCOL = {"Carcinoma": C[0], "Oncofetal": C[1], "Normal epithelium": C[2],
        "Fibroblast / CAF": C[3], "Immune": C[4], "Vascular": C[6], "Other": "#c9c7c0"}

sl = cells[cells["slide_ID_numeric"] == 1]
print(f"CosMx slide 1: {len(sl):,} cells, {sl.patient_ID.nunique()} patients")

fig, axes = plt.subplots(1, 3, figsize=(15.0, 5.0),
                         gridspec_kw={"width_ratios": [1.05, 1.05, 1.0]})

ax = axes[0]
for g in ["Other", "Immune", "Vascular", "Normal epithelium", "Fibroblast / CAF",
          "Carcinoma", "Oncofetal"]:
    s = sl[sl["group"] == g]
    ax.scatter(s["x_slide_mm"], s["y_slide_mm"], s=0.14, color=GCOL[g], alpha=0.55,
               edgecolor="none", rasterized=True, label=g)
ax.set_aspect("equal"); ax.set_xlabel("Slide x (mm)"); ax.set_ylabel("Slide y (mm)")
leg = ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.13), ncols=4, fontsize=7,
                markerscale=26, handletextpad=0.3)
ax.set_title(f"CosMx slide 1 — {len(sl):,} single cells")
tag(ax, "a")

ax = axes[1]
ax.scatter(sl["x_slide_mm"], sl["y_slide_mm"], s=0.10, color="#e2e0da", alpha=0.6,
           edgecolor="none", rasterized=True)
for g, col, lab in [("oncofetal", C[1], "Oncofetal cells"),
                    ("trophocyte-like CAF", C[7], "Trophocyte-like CAFs")]:
    s = sl[sl["cellType"] == g]
    ax.scatter(s["x_slide_mm"], s["y_slide_mm"], s=0.9, color=col, alpha=0.85,
               edgecolor="none", rasterized=True, label=f"{lab} (n={len(s):,})")
ax.set_aspect("equal"); ax.set_xlabel("Slide x (mm)")
ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.13), ncols=1, fontsize=7.4,
          markerscale=12)
ax.set_title("The states that emerge with invasion")
tag(ax, "b")

ax = axes[2]
order = ["im", "T1.sm1", "T1.sm3"]
xlab = ["Intramucosal", "T1 sm1", "T1 sm3"]
show = ["trophocyte-like CAF", "MMP CAF", "oncofetal", "pericytes"]
cols = [C[7], C[3], C[1], C[6]]
for ct, col in zip(show, cols):
    if ct not in prop.columns:
        continue
    mm = [prop.loc[prop.tumorStage == s, ct].mean() * 100 for s in order]
    for pid, g in prop.groupby("patient_ID"):
        xs, ys = [], []
        for i, s in enumerate(order):
            v = g.loc[g.tumorStage == s, ct]
            if len(v):
                xs.append(i); ys.append(float(v.iloc[0]) * 100)
        ax.plot(xs, ys, color=col, lw=0.7, alpha=0.28)
    r = CT[CT.cellType == ct].iloc[0]
    ax.plot(range(3), mm, color=col, lw=2.4, marker="o", ms=6,
            markeredgecolor=SURF, markeredgewidth=1,
            label=f"{ct}  (ρ={r.spearman_rho:+.2f}, P={r.p_trend:.3f})")
ax.set_xticks(range(3)); ax.set_xticklabels(xlab)
ax.set_ylabel("Percent of cells in the specimen")
ax.legend(loc="upper left", fontsize=7.2)
ax.grid(axis="y"); ax.set_axisbelow(True)
ax.set_title("Composition across invasive progression")
tag(ax, "c")
save(fig, "fig14_cosmx_spatial")

print("written:", sorted(f for f in os.listdir("figures") if f.startswith(("fig12", "fig13", "fig14"))))
