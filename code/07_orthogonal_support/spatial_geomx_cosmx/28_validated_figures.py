#!/usr/bin/env python3
"""28 - Figures 15-17: the validated gene set and the CosMx single-cell evidence."""
import pandas as pd, numpy as np, json, os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from matplotlib.colors import LinearSegmentedColormap

os.makedirs("figures", exist_ok=True)
C = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#4a3aa7", "#e34948"]
POS, NEG = "#e34948", "#2a78d6"
SEQ = ["#cde2fb", "#9ec5f4", "#6da7ec", "#3987e5", "#256abf", "#184f95", "#0d366b"]
INK, INK2, INK3 = "#0b0b0b", "#52514e", "#8a8880"
SURF = "#fcfcfb"
SEQCM = LinearSegmentedColormap.from_list("seq", SEQ)
EXPCM = LinearSegmentedColormap.from_list("exp", ["#e9e7e1", "#ffd9a0", "#f59b42",
                                                  "#e34948", "#8f1d2e"])
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


ONC, TSG = "#e34948", "#2a78d6"
STAGE_ORDER = ["im", "T1.sm1", "T1.sm3"]
STAGE_LAB = {"im": "Intramucosal", "T1.sm1": "T1 sm1", "T1.sm3": "T1 sm3"}
STAGE_COL = {"im": "#1baf7a", "T1.sm1": "#eda100", "T1.sm3": "#e34948"}

V = pd.read_csv("results/validated_risk_genes_all.csv")
val = pd.read_csv("results/validated_priority.csv")
CTv = pd.read_csv("results/validated_by_celltype.csv")
S = json.load(open("results/validated_summary.json"))
rep = pd.read_csv("results/cosmx_geomx_replication.csv")
B = pd.read_csv("results/cosmx_crc_vs_normal_epi.csv")
pb = pd.read_csv("spatial/cx_pseudobulk.csv.gz")
VF = json.load(open("results/validated_vs_failed.json"))
pri = pd.read_csv("results/integrated_gene_priority.csv")

# ============================================================ FIGURE 15
fig = plt.figure(figsize=(11.6, 8.4))
gs = fig.add_gridspec(2, 2, height_ratios=[1.0, 1.05], width_ratios=[1.0, 1.28],
                      hspace=0.52, wspace=0.28)

# (a) the filter
ax = fig.add_subplot(gs[0, 0]); tag(ax, "a")
steps = [("Risk genes (TWAS)", 347, INK3),
         ("Measured in GeoMx epithelium", S["n_measured"], C[0]),
         ("Significant, TWAS direction", S["n_validated"], "#1baf7a")]
ypos = np.arange(len(steps))[::-1]
for y, (lab, n, col) in zip(ypos, steps):
    ax.barh(y, n, color=col, height=0.55, alpha=0.9)
    ax.text(n + 7, y, f"{n}", va="center", ha="left", fontsize=9.5, fontweight="bold")
ax.set_yticks(ypos); ax.set_yticklabels([s0[0] for s0 in steps], fontsize=8)
ax.set_xlim(0, 420); ax.set_xlabel("Genes"); ax.set_ylim(-0.75, 2.6)
ax.set_title("From 347 risk genes to a directionally coherent set", pad=12)
drop = [("never significant", S["n_never_significant"]),
        ("significant but opposite", S["n_sig_opposite_only"]),
        ("significant both ways", S["n_excluded_by_contradiction"])]
ax.text(0.0, -0.30, "Excluded:  " + "   ".join(f"{k} {v}" for k, v in drop),
        transform=ax.transAxes, fontsize=7.4, color=INK3, va="top")

# (b) waterfall of the validated set
ax = fig.add_subplot(gs[0, 1]); tag(ax, "b")
w = val.sort_values("lead_log2FC").reset_index(drop=True)
cols = [ONC if c == "oncogene-like" else TSG for c in w["class"]]
ax.bar(np.arange(len(w)), w["lead_log2FC"], color=cols, width=0.82)
ax.axhline(0, color=INK3, lw=0.8)
lab_idx = list(range(6)) + list(range(len(w) - 8, len(w)))
for i in lab_idx:
    v = w["lead_log2FC"][i]
    ax.text(i, v + (0.13 if v > 0 else -0.13), w["gene"][i], rotation=90,
            ha="center", va="bottom" if v > 0 else "top", fontsize=6.6, color=INK2)
ax.set_xticks([]); ax.set_ylabel("log$_2$ FC, leading transition")
ax.set_ylim(w["lead_log2FC"].min() - 1.3, w["lead_log2FC"].max() + 1.4)
ax.set_title(f"The {len(w)} validated genes: {S['n_tsg_like']} tumour-suppressor-like, "
             f"{S['n_oncogene_like']} oncogene-like")
ax.legend(handles=[Patch(color=ONC, label="oncogene-like (risk-increasing, up)"),
                   Patch(color=TSG, label="tumour-suppressor-like (protective, down)")],
          loc="upper left", fontsize=7.4)

# (c) per TWAS cell type
ax = fig.add_subplot(gs[1, 0]); tag(ax, "c")
ct = CTv[CTv["n_measured"] >= 10].copy()
ct["lab"] = ct["tissue"].str[0] + " " + ct["cell_type"]
ct = ct.sort_values("frac_validated")
y = np.arange(len(ct))
ax.barh(y - 0.19, ct["oncogene_like"] / ct["n_measured"], height=0.36, color=ONC,
        label="oncogene-like")
ax.barh(y + 0.19, ct["tumour_suppressor_like"] / ct["n_measured"], height=0.36, color=TSG,
        label="tumour-suppressor-like")
ax.set_yticks(y); ax.set_yticklabels(ct["lab"], fontsize=7.6)
for yi, (_, r) in zip(y, ct.iterrows()):
    ax.text(r["frac_validated"] + 0.012, yi, f"n={r['n_measured']:.0f}", va="center",
            ha="left", fontsize=6.6, color=INK3)
ax.set_xlim(0, 0.56)
ax.set_xlabel("Fraction of measured risk genes validated")
ax.set_title("Validation rate by the cell type the gene was called in", pad=10)
ax.legend(loc="lower right", fontsize=7.2, bbox_to_anchor=(1.0, 0.02))
ax.text(0.98, 0.30, "no cell type differs after FDR", transform=ax.transAxes,
        fontsize=7.2, color=INK3, style="italic", ha="right")

# (d) validated vs failed
ax = fig.add_subplot(gs[1, 1]); tag(ax, "d")
m = V.merge(pri[["gene", "ko_impact", "degree_ppi", "tractability_z"]], on="gene", how="left")
panels = [("ko_impact", "Knockout impact", VF["ko_impact"]),
          ("degree_ppi", "PPI degree (log$_{10}$)", VF["degree_ppi"]),
          ("tractability_z", "Tractability z", VF["tractability_z"])]
for k, (col, lab, st) in enumerate(panels):
    a = pd.to_numeric(m.loc[m["validated"], col], errors="coerce").dropna()
    b = pd.to_numeric(m.loc[~m["validated"], col], errors="coerce").dropna()
    if col == "degree_ppi":
        a, b = np.log10(a + 1), np.log10(b + 1)
    for j, (vals, colr) in enumerate([(b, INK3), (a, "#1baf7a")]):
        x = k + (j - 0.5) * 0.42
        parts = ax.violinplot([vals], positions=[x], widths=0.36, showextrema=False)
        for pc in parts["bodies"]:
            pc.set_facecolor(colr); pc.set_alpha(0.35); pc.set_edgecolor("none")
        ax.plot([x - 0.13, x + 0.13], [np.median(vals)] * 2, color=colr, lw=2.2)
    ax.text(k, 0.985, f"AUROC {st['auroc']:.2f}   P = {st['p']:.2g}",
            transform=ax.get_xaxis_transform(), ha="center", va="top",
            fontsize=7.2, color=INK2)
ax.set_xticks(range(3)); ax.set_xticklabels([p[1] for p in panels], fontsize=7.8)
ax.set_ylabel("Value"); ax.set_ylim(-2.6, 3.9)
ax.set_title("More connected and more disruptive, not more druggable", pad=10)
ax.legend(handles=[Patch(color=INK3, alpha=.45, label="not validated (209)"),
                   Patch(color="#1baf7a", alpha=.45, label="validated (62)")],
          loc="lower left", fontsize=7.2)
save(fig, "fig15_validated_set")
print("fig15 done")

# ============================================================ FIGURE 16  (CosMx maps)
z = np.load("spatial/cx_plot_cells.npz", allow_pickle=True)
X, Y = z["x"].astype(float), z["y"].astype(float)
E = z["expr"].astype(np.float32)
pgenes = [str(g) for g in z["genes"]]
ctype = np.array([str(s) for s in z["cellType"]])
stage = np.array([str(s) for s in z["stage"]])
Y = Y.max() - Y

SHOW = [g for g in ["CTNNB1", "CEBPB", "TRIM28", "FAS"] if g in pgenes]
fig = plt.figure(figsize=(11.6, 8.8))
gs = fig.add_gridspec(2, 3, height_ratios=[1.0, 1.0], hspace=0.26, wspace=0.36)

GROUP = {"CRC": "Carcinoma", "normal/undefined epithelium": "Normal epithelium",
         "oncofetal": "Oncofetal"}
GCOL = {"Carcinoma": "#e34948", "Normal epithelium": "#1baf7a", "Oncofetal": "#eda100",
        "Fibroblast / stroma": "#8a8880", "Immune": "#9ec5f4", "Other": "#d8d7d2"}


def grp(c):
    if c in GROUP:
        return GROUP[c]
    if "fibro" in c or "CAF" in c or "trophocyte" in c or "telocyte" in c or c == "stromal.1":
        return "Fibroblast / stroma"
    if c.startswith("T.") or c in ("B cell", "plasma", "macrophage", "mast cells",
                                   "MKI67+ lymphocytes"):
        return "Immune"
    return "Other"


g_of = np.array([grp(c) for c in ctype])
ax = fig.add_subplot(gs[0, 0]); tag(ax, "a")
for k in ["Other", "Immune", "Fibroblast / stroma", "Oncofetal", "Normal epithelium",
          "Carcinoma"]:
    m = g_of == k
    ax.scatter(X[m], Y[m], s=0.10, c=GCOL[k], linewidths=0, rasterized=True, label=k)
ax.set_aspect("equal"); ax.axis("off")
npat = len(set(z["patient"].tolist()))
ax.set_title(f"Slide {int(z['slide'])}: {int(z['n_cells_slide']):,} cells, {npat} patients",
             loc="left")
ax.legend(markerscale=26, loc="lower left", bbox_to_anchor=(-0.02, 0.0),
          ncol=1, fontsize=6.6, handletextpad=0.1, labelspacing=0.28)

for k, g in enumerate(SHOW):
    ax = fig.add_subplot(gs[(k + 1) // 3, (k + 1) % 3])
    tag(ax, "bcde"[k])
    j = pgenes.index(g)
    e = E[:, j]
    neg = e <= 0
    ax.scatter(X[neg], Y[neg], s=0.08, c="#e9e7e1", linewidths=0, rasterized=True)
    pos = ~neg
    o = np.argsort(e[pos])
    sc = ax.scatter(X[pos][o], Y[pos][o], s=0.30, c=e[pos][o], cmap=EXPCM,
                    vmin=0, vmax=np.quantile(e[pos], 0.98), linewidths=0, rasterized=True)
    ax.set_aspect("equal"); ax.axis("off")
    row = B[B["gene"] == g]
    # the map is one slide; the statistic is the whole cohort, so both are named
    sub = (f"CosMx Δ {row['delta'].iloc[0]:+.3f}, P = {row['p'].iloc[0]:.3f} "
           f"(n = {int(row['n_patients'].iloc[0])})") if len(row) else ""
    ax.set_title(g, loc="left", fontsize=9.6, pad=12)
    ax.text(0.0, 1.012, sub, transform=ax.transAxes, fontsize=7.2, color=INK2,
            va="baseline", ha="left")
    cb = fig.colorbar(sc, ax=ax, fraction=0.030, pad=0.005, shrink=0.66,
                      location="left", anchor=(0.0, 0.05))
    cb.ax.tick_params(labelsize=6.2, pad=1); cb.outline.set_visible(False)
    cb.set_label("log$_2$ norm. counts", fontsize=6.2, labelpad=2)

# (f) replication CosMx vs GeoMx
ax = fig.add_subplot(gs[1, 2]); tag(ax, "f")
rv = rep[rep["geomx_validated"]]
rn = rep[~rep["geomx_validated"]]
ax.axhline(0, color="#e0dfda", lw=0.8); ax.axvline(0, color="#e0dfda", lw=0.8)
ax.scatter(rn["geomx_log2FC"], rn["delta"], s=13, c=INK3, alpha=0.45, linewidths=0,
           label=f"not validated ({len(rn)})")
ax.scatter(rv["geomx_log2FC"], rv["delta"], s=26, c="#1baf7a", linewidths=0,
           label=f"validated ({len(rv)})")
for _, r in rv.reindex(rv["delta"].abs().sort_values(ascending=False).index).head(6).iterrows():
    ax.annotate(r["gene"], (r["geomx_log2FC"], r["delta"]), fontsize=6.8,
                xytext=(3, 3), textcoords="offset points", color=INK2)
cc = json.load(open("results/cosmx_percell_summary.json"))["replication"]
ax.set_xlabel("GeoMx log$_2$ FC, normal → carcinoma")
ax.set_ylabel("CosMx Δ, carcinoma vs normal epithelium")
ax.set_title("Independent patients, independent platform")
ax.text(0.03, 0.96, f"validated genes: ρ = {cc['rho_validated']:.2f}\n"
                    f"P = {cc['p_validated']:.1e}, "
                    f"{100*cc['sign_agreement_validated']:.0f} % same sign",
        transform=ax.transAxes, va="top", fontsize=7.4, color=INK2)
ax.legend(loc="lower right", fontsize=7.2)
save(fig, "fig16_cosmx_percell_expression")
print("fig16 done")

# ============================================================ FIGURE 17  (stages)
fig = plt.figure(figsize=(11.6, 7.4))
gs = fig.add_gridspec(2, 4, height_ratios=[1.0, 1.06], hspace=0.52, wspace=0.34)

crc = ctype == "CRC"
for k, g in enumerate(SHOW):
    ax = fig.add_subplot(gs[0, k]); tag(ax, "abcd"[k])
    j = pgenes.index(g)
    hi = float(np.quantile(E[crc, j].astype(float), 0.995)) or 1.0
    for si, s in enumerate(STAGE_ORDER):
        m = crc & (stage == s)
        e = E[m, j].astype(float)
        if not m.sum():
            continue
        parts = ax.violinplot([np.clip(e, 0, hi)], positions=[si], widths=0.74,
                              showextrema=False)
        for pc in parts["bodies"]:
            pc.set_facecolor(STAGE_COL[s]); pc.set_alpha(0.45); pc.set_edgecolor("none")
        ax.plot([si - 0.22, si + 0.22], [e.mean()] * 2, color=STAGE_COL[s], lw=2.4)
        ax.text(si, 1.005, f"{m.sum()//1000}k cells", transform=ax.get_xaxis_transform(),
                ha="center", va="bottom", fontsize=6.2, color=INK3)
    # per-patient means on top
    d = pb[(pb["cellType"] == "CRC") & (pb["gene"] == g)]
    for _, r in d.iterrows():
        if r["tumorStage"] in STAGE_ORDER:
            ax.plot(STAGE_ORDER.index(r["tumorStage"]) + 0.30, r["mean_norm"], "o",
                    ms=3.4, color=INK, alpha=0.8)
    ax.set_xticks(range(3))
    ax.set_xticklabels([STAGE_LAB[s].replace("Intramucosal", "im") for s in STAGE_ORDER],
                       fontsize=7.2)
    ax.set_title(g, fontsize=9.5, pad=14)
    if k == 0:
        ax.set_ylabel("log$_2$ normalised counts per cell")
    ax.set_ylim(-0.04 * hi, hi * 1.16)
fig.text(0.5, 0.985, "Per-cell expression in carcinoma cells across the invasive sequence "
                     "(violins: all cells; dots: patient means)",
         ha="center", fontsize=9, color=INK2)

# (e) paired carcinoma vs normal epithelium, per patient
ax = fig.add_subplot(gs[1, :2]); tag(ax, "e")
show2 = [g for g in ["CTNNB1", "CEBPB", "NFAT5", "ZNRF3", "HRAS", "FAS"]
         if g in set(pb["gene"])]
piv = (pb[pb["cellType"].isin(["CRC", "normal/undefined epithelium"])]
       .pivot_table(index=["gene", "patient_ID"], columns="cellType", values="mean_norm")
       .dropna())
for k, g in enumerate(show2):
    d = piv.loc[g] if g in piv.index.get_level_values(0) else None
    if d is None:
        continue
    for _, r in d.iterrows():
        col = POS if r["CRC"] > r["normal/undefined epithelium"] else NEG
        ax.plot([k - 0.17, k + 0.17],
                [r["normal/undefined epithelium"], r["CRC"]], "-o", ms=2.6, lw=0.8,
                color=col, alpha=0.65)
    row = B[B["gene"] == g]
    if len(row):
        ax.text(k, 0.965, f"P = {row['p'].iloc[0]:.3f}", transform=ax.get_xaxis_transform(),
                ha="center", va="top", fontsize=6.8, color=INK2)
ax.set_xticks(range(len(show2))); ax.set_xticklabels(show2, fontsize=8)
ax.set_ylim(top=ax.get_ylim()[1] * 1.20)
ax.set_ylabel("Mean log$_2$ normalised counts")
ax.set_title("Paired within patient: normal epithelium → carcinoma", pad=12)
ax.legend(handles=[Line2D([], [], color=POS, marker="o", ms=3, label="higher in carcinoma"),
                   Line2D([], [], color=NEG, marker="o", ms=3, label="higher in normal")],
          loc="center right", fontsize=7.2, bbox_to_anchor=(1.0, 0.62))

# (f) cell-type heatmap of the validated genes
ax = fig.add_subplot(gs[1, 2:]); tag(ax, "f")
vg = [g for g in val["gene"] if g in set(pb["gene"])][:26]
order = ["CRC", "normal/undefined epithelium", "oncofetal", "myofibroblast.1",
         "myofibroblast.2", "trophocyte-like CAF", "endothelial.stalk", "pericytes",
         "T.CD4", "macrophage", "plasma", "B cell"]
order = [o for o in order if o in set(pb["cellType"])]
M = (pb[pb["gene"].isin(vg) & pb["cellType"].isin(order)]
     .pivot_table(index="gene", columns="cellType", values="mean_norm", aggfunc="mean")
     .reindex(index=vg, columns=order))
Mz = M.sub(M.mean(1), axis=0).div(M.std(1).replace(0, np.nan), axis=0)
DIVCM = LinearSegmentedColormap.from_list("div", ["#184f95", "#9ec5f4", "#f2f1ed",
                                                  "#f0a09f", "#c3202f"])
im = ax.imshow(Mz.values, cmap=DIVCM, aspect="auto", vmin=-1.8, vmax=1.8)
ax.set_yticks(range(len(vg))); ax.set_yticklabels(vg, fontsize=6.2)
ax.set_xticks(range(len(order)))
ax.set_xticklabels([o.replace("normal/undefined epithelium", "normal epi.")
                    .replace("trophocyte-like CAF", "tropho. CAF")
                    .replace("endothelial.stalk", "endoth.") for o in order],
                   rotation=42, ha="right", fontsize=6.8)
cb = fig.colorbar(im, ax=ax, fraction=0.03, pad=0.02)
cb.set_label("row z-score", fontsize=7); cb.ax.tick_params(labelsize=6.5)
cb.outline.set_visible(False)
ax.set_title("Where the validated genes are expressed", pad=10)
save(fig, "fig17_cosmx_stage_celltype")
print("fig17 done")
