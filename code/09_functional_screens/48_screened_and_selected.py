#!/usr/bin/env python3
"""
48 - Figures 29 and 30.

FIGURE 29 is Figure 16f restricted to the high-confidence gene list supplied by the
experimenters (high-confidence-genes.xlsx, 109 genes, every one of them among the 347
TWAS risk genes). Figure 16f asks whether the two spatial platforms agree about a risk
gene: GeoMx measures bulk stage transitions in one set of patients, CosMx measures single
cells in another, and the scatter of one against the other is the only fully independent
replication in this report. Restricting it to the high-confidence list answers the
sharper question - not "do the platforms agree", which section 2.14 settled over all 121
genes, but "do they agree about the genes this project is actually carrying forward". A
shortlist that replicated worse than the genes it was drawn from would be a shortlist
built on noise.

The scatter carries a third layer, because each gene also has a functional verdict from
the viability screens: whether removing it changed the growth of a colorectal cancer cell
line. That makes one set of axes carry genetics, tissue and cell biology at once.

FIGURE 30 applies the Figure 28 treatment - per-cell spatial map above, per-patient paired
test below - to genes chosen by an explicit rule rather than by hand:

  significant   nominal CosMx P < 0.05 in the paired carcinoma-versus-normal-epithelium
                test across the seven patients carrying both populations
  relevant      CRC-related in the concrete sense used in section 2.4 (the gene is itself
                a colorectal-cancer somatic driver, or is in direct protein-protein
                contact with one), OR druggable at tier 1-2 (an approved drug binds it)
  drawable      per-cell expression is available, which restricts the field to the genes
                extracted from the Seurat object

Genes already drawn in Figure 28 are excluded so the figure adds rather than repeats, and
the remaining candidates are ranked by the size of the single-cell difference.

Outputs
  results/screened_replication.csv
  results/cosmx_selected_genes.csv
  figures/fig29_highconf_replication.{png,pdf}
  figures/fig30_cosmx_selected.{png,pdf}
"""
import pandas as pd, numpy as np, json, os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from matplotlib.colors import LinearSegmentedColormap
from scipy.stats import spearmanr, fisher_exact, wilcoxon, binomtest

os.makedirs("figures", exist_ok=True)
INK, INK2, INK3 = "#0b0b0b", "#52514e", "#8a8880"
SURF = "#fcfcfb"
EXPCM = LinearSegmentedColormap.from_list("exp", ["#e9e7e1", "#ffd9a0", "#f59b42",
                                                  "#e34948", "#8f1d2e"])
NORMC, CRCC = "#2a78d6", "#c3202f"
GREEN, GREY, AMBER = "#1baf7a", "#c9c7c1", "#e8a33d"
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


def tag(ax, s, dx=-0.13, dy=1.06):
    ax.text(dx, dy, s, transform=ax.transAxes, fontsize=11, fontweight="bold",
            color=INK, va="top", ha="left")


ALIAS = {"ACAS2": "ACSS2"}
NORMAL, TUMOUR = "normal/undefined epithelium", "CRC"
LINES = ["HCT116", "RKO", "SW480"]

rep = pd.read_csv("results/cosmx_geomx_replication.csv")
HCF = "/mnt/user-data/uploads/AI-gene-drug/high-confidence-genes.xlsx"
HC = set(pd.read_excel(HCF)["gene_name"].astype(str).str.strip())
HC = {ALIAS.get(g, g) for g in HC}

C = pd.read_csv("results/crispr_by_cellline.csv")
S = pd.read_csv("results/sirna_by_cellline.csv")
called = set(C[C.nominal].gene) | set(S[S.nominal & S.qc_pass].gene)
assayed = set(C.gene) | set(S.gene)

rep["screened"] = rep.gene.isin(HC)          # "screened" now means high-confidence
rep["func_called"] = rep.gene.isin(called)
rep["func_assayed"] = rep.gene.isin(assayed)
rep["replicates"] = np.sign(rep.geomx_log2FC) == np.sign(rep.delta)
W = rep[rep.screened].copy()
W.to_csv("results/highconf_replication.csv", index=False)


def rho_of(d):
    r, p = spearmanr(d.geomx_log2FC, d.delta)
    same = int((np.sign(d.geomx_log2FC) == np.sign(d.delta)).sum())
    return float(r), float(p), same, len(d)


R_all = rho_of(rep)
R_scr = rho_of(W)
R_oth = rho_of(rep[~rep.screened])
tab = [[int((W.replicates & W.func_called).sum()), int((W.replicates & ~W.func_called).sum())],
       [int((~W.replicates & W.func_called).sum()), int((~W.replicates & ~W.func_called).sum())]]
orr, pf = fisher_exact(tab)
J29 = {"n_high_conf": len(HC), "n_hc_cosmx": len(W),
       "n_func_assayed": int(W.func_assayed.sum()),
       "n_func_called": int(W.func_called.sum()),
       "n_validated": int(W.geomx_validated.sum()),
       "rho_all": R_all, "rho_screened": R_scr, "rho_other": R_oth,
       "table": tab, "odds_ratio": float(orr), "p_fisher": float(pf)}
json.dump(J29, open("results/highconf_replication.json", "w"), indent=1)
print(json.dumps(J29, indent=1))


def place_labels(ax, xs, ys, names, fontsize=7.0, color=INK2):
    """Greedy non-overlapping point labels, with a leader line when the label is
    pushed away from its point. Offsets are in points; boxes are compared in
    display coordinates against the points and the labels already placed."""
    fig_ = ax.figure
    fig_.canvas.draw()
    rend = fig_.canvas.get_renderer()
    sc = fig_.dpi / 72.0
    inv = ax.transData.inverted()
    pts = [ax.transData.transform((x, y)) for x, y in zip(xs, ys)]
    boxes = [(px - 11, py - 11, px + 11, py + 11) for px, py in pts]
    cands = [(7, 3), (7, -10), (-7, 3), (-7, -10), (7, 13), (7, -20),
             (-7, 13), (-7, -20), (16, 3), (-16, 3), (16, -10), (-16, -10),
             (7, 22), (7, -29), (-7, 22), (-7, -29), (24, 3), (-24, 3)]
    for (x, y), nm in zip(zip(xs, ys), names):
        p0 = ax.transData.transform((x, y))
        for dx, dy in cands:
            t = ax.text(x, y, nm, fontsize=fontsize, color=color, zorder=6,
                        ha="left" if dx > 0 else "right", va="center",
                        transform=ax.transData + matplotlib.transforms.ScaledTranslation(
                            dx / 72.0, dy / 72.0, fig_.dpi_scale_trans))
            bb = t.get_window_extent(renderer=rend)
            box = (bb.x0 - 1.5, bb.y0 - 1.5, bb.x1 + 1.5, bb.y1 + 1.5)
            if not any(not (box[2] < b[0] or box[0] > b[2] or
                            box[3] < b[1] or box[1] > b[3]) for b in boxes):
                boxes.append(box)
                if abs(dx) > 8 or abs(dy) > 8:      # far enough to need a leader
                    ex = p0[0] + dx * sc - (6 * sc if dx > 0 else -6 * sc)
                    ey = p0[1] + dy * sc
                    x1, y1 = inv.transform((ex, ey))
                    ax.plot([x, x1], [y, y1], color="#c4c2bd", lw=0.5, zorder=5,
                            solid_capstyle="butt")
                break
            t.remove()


# ---------------------------------------------------------------- FIGURE 29
# Restricted deliberately: the 41 high-confidence genes carried on the CosMx panel,
# and nothing but the two spatial platforms. No other gene sets, no viability screens.
W = W.copy()
W["agree"] = np.sign(W.geomx_log2FC) == np.sign(W.delta)
bt = binomtest(int(W.agree.sum()), len(W), 0.5)
Wz = W.assign(gz=(W.geomx_log2FC - W.geomx_log2FC.mean()) / W.geomx_log2FC.std(),
              cz=(W.delta - W.delta.mean()) / W.delta.std())
AGREE_C, DISAGREE_C = "#1f7a8c", "#d98b2b"

fig = plt.figure(figsize=(14.6, 7.4))
gs = fig.add_gridspec(1, 3, width_ratios=[1.22, 0.86, 0.92], wspace=0.40)

ax = fig.add_subplot(gs[0, 0]); tag(ax, "a", dx=-0.15, dy=1.07)
xl = float(np.abs(W.geomx_log2FC).max()) * 1.20
pad = (W.delta.max() - W.delta.min()) * 0.12
ylo, yhi = float(W.delta.min() - pad), float(W.delta.max() + pad)
ax.add_patch(plt.Rectangle((0, 0), xl, yhi, color="#eef5f6", zorder=0))
ax.add_patch(plt.Rectangle((-xl, ylo), xl, -ylo, color="#eef5f6", zorder=0))
ax.axhline(0, color="#cfcdc8", lw=0.9, zorder=1)
ax.axvline(0, color="#cfcdc8", lw=0.9, zorder=1)
b, a0 = np.polyfit(W.geomx_log2FC, W.delta, 1)
xs = np.linspace(-xl, xl, 50)
ax.plot(xs, b * xs + a0, color="#8a8880", lw=1.0, ls="--", zorder=2)
for m, col, lab in [(W.agree, AGREE_C, f"same direction on both ({int(W.agree.sum())})"),
                    (~W.agree, DISAGREE_C, f"opposite directions ({int((~W.agree).sum())})")]:
    d = W[m]
    ax.scatter(d.geomx_log2FC, d.delta, s=46, c=col, edgecolor="white", linewidth=0.6,
               zorder=4, label=lab)
lab_idx = set(W.delta.abs().nlargest(8).index) | set(W.geomx_log2FC.abs().nlargest(4).index) \
          | set(W[~W.agree].delta.abs().nlargest(3).index)
LBL = W.loc[sorted(lab_idx)]
ax.set_xlim(-xl, xl); ax.set_ylim(ylo, yhi)
ax.set_xlabel("GeoMx log$_2$ FC, normal → carcinoma (bulk regions)")
place_labels(ax, LBL.geomx_log2FC.tolist(), LBL.delta.tolist(), LBL.gene.tolist())
ax.set_ylabel("CosMx Δ, carcinoma cells vs normal epithelium (single cells)")
ax.set_title("The 41 high-confidence genes, two platforms")
ax.text(0.03, 0.975,
        f"Spearman ρ = {R_scr[0]:.2f}, P = {R_scr[1]:.1e}\n"
        f"Pearson r = {np.corrcoef(W.geomx_log2FC, W.delta)[0, 1]:.2f}\n"
        f"{R_scr[2]} of {R_scr[3]} genes agree in direction ({100 * R_scr[2] / R_scr[3]:.0f} %)",
        transform=ax.transAxes, va="top", fontsize=7.8, color=INK2)
ax.legend(loc="lower right", fontsize=7.4, borderpad=0.25)
ax.grid(True, zorder=0); ax.set_axisbelow(True)
ax.text(0.0, -0.125,
        "Shaded quadrants are where the platforms agree. The dashed\n"
        "line is the least-squares fit, drawn for orientation only.\n"
        "Different patients, different chemistry, bulk regions against\n"
        "single cells: nothing is shared between the two measurements.",
        transform=ax.transAxes, fontsize=7.4, color=INK2, va="top")

ax = fig.add_subplot(gs[0, 1]); tag(ax, "b", dx=-0.19, dy=1.07)
gr = W.geomx_log2FC.rank(); cr = W.delta.rank()
ax.plot([0, len(W) + 1], [0, len(W) + 1], color="#d8d7d2", lw=1.0, zorder=1)
ax.scatter(gr[W.agree], cr[W.agree], s=40, c=AGREE_C, edgecolor="white", linewidth=0.5,
           zorder=3)
ax.scatter(gr[~W.agree], cr[~W.agree], s=40, c=DISAGREE_C, edgecolor="white", linewidth=0.5,
           zorder=3)
ax.set_xlim(0, len(W) + 1); ax.set_ylim(0, len(W) + 1)
ax.set_xlabel("Rank of GeoMx log$_2$ FC (1 = most negative)")
ax.set_ylabel("Rank of CosMx Δ")
ax.set_title("The same data as ranks")
ax.text(0.03, 0.975, f"ρ = {R_scr[0]:.2f} is the correlation\nof these ranks (n = {len(W)})",
        transform=ax.transAxes, va="top", fontsize=7.8, color=INK2)
ax.grid(True); ax.set_axisbelow(True)
ax.text(0.0, -0.125,
        "Spearman is a rank statistic, so this is what it sees. A monotone\n"
        "relationship would sit on the diagonal; the scatter around it is\n"
        "the width of the agreement, and it is not carried by the two\n"
        f"largest genes. Direction alone: {R_scr[2]}/{len(W)}, binomial P = {bt.pvalue:.1e}.",
        transform=ax.transAxes, fontsize=7.4, color=INK2, va="top")

ax = fig.add_subplot(gs[0, 2]); tag(ax, "c", dx=-0.46, dy=1.07)
Wz = Wz.sort_values("cz")
y = np.arange(len(Wz))
for yi, (_, r) in zip(y, Wz.iterrows()):
    ax.plot([r.gz, r.cz], [yi, yi], color=AGREE_C if r.agree else DISAGREE_C,
            lw=1.0, alpha=0.55, zorder=2)
ax.scatter(Wz.gz, y, s=16, c="#6f7f95", edgecolor="white", linewidth=0.4, zorder=3,
           label="GeoMx (bulk stage)")
ax.scatter(Wz.cz, y, s=16, c="#b03a48", edgecolor="white", linewidth=0.4, zorder=3,
           label="CosMx (single cell)")
ax.axvline(0, color="#cfcdc8", lw=0.9, zorder=1)
ax.set_yticks(y); ax.set_yticklabels(Wz.gene, fontsize=5.8)
ax.set_ylim(-1, len(Wz))
ax.set_xlabel("Standardised within platform, across these 41 genes (z)")
ax.set_title("Gene by gene")
ax.legend(loc="lower right", fontsize=7.0, borderpad=0.25)
ax.grid(axis="x"); ax.set_axisbelow(True)
ax.text(0.0, -0.088,
        "Each gene twice, once per platform, on a common scale.\n"
        "Connectors are coloured by whether the two measurements\n"
        "fall on the same side of zero.",
        transform=ax.transAxes, fontsize=7.4, color=INK2, va="top")

fig.savefig("figures/fig29_highconf_replication.png", dpi=220, bbox_inches="tight")
fig.savefig("figures/fig29_highconf_replication.pdf", bbox_inches="tight")
plt.close(fig)
J29["binomial_p"] = float(bt.pvalue)
J29["n_agree"] = int(W.agree.sum())
J29["pearson_r"] = float(np.corrcoef(W.geomx_log2FC, W.delta)[0, 1])
J29["disagree_genes"] = sorted(W.loc[~W.agree, "gene"].tolist())
json.dump(J29, open("results/highconf_replication.json", "w"), indent=1)
print("fig29 done")

# ---------------------------------------------------------------- FIGURE 30 selection
B = pd.read_csv("results/cosmx_crc_vs_normal_epi.csv")
dr = pd.read_csv("results/gene_druggability.csv").set_index("gene")
px = pd.read_csv("results/risk_gene_driver_proximity.csv").set_index("gene")
V = pd.read_csv("results/validated_risk_genes_all.csv").set_index("gene")
ALREADY = ["CSNK2B", "CTNNB1", "UBE2V1", "TRIM28"]          # drawn in Figure 28

SRC = {}
for f in ["spatial/cx_plot_cells.npz", "spatial/cx_plot_cells3.npz",
          "spatial/cx_plot_cells2.npz"]:
    z = np.load(f, allow_pickle=True)
    for g in [str(x) for x in z["genes"]]:
        SRC.setdefault(g, f)

cand = []
for _, r in B[B.p < 0.05].iterrows():
    g = r.gene
    d = dr.loc[g] if g in dr.index else None
    p = px.loc[g] if g in px.index else None
    is_driver = bool(p["is_driver"]) if p is not None else False
    contact = (float(p["dist_to_nearest_driver"]) == 1.0) if p is not None else False
    tier = str(d["tier"]) if d is not None else ""
    cand.append({"gene": g, "cosmx_p": r.p, "delta": r.delta,
                 "twas_direction": r.twas_direction, "concordant": bool(r.concordant),
                 "is_driver": is_driver, "driver_contact": contact,
                 "crc_related": is_driver or contact,
                 "druggable": tier.startswith(("T1", "T2")), "tier": tier,
                 "validated": bool(V.loc[g, "validated"]) if g in V.index else False,
                 "drawable": g in SRC, "already_in_fig28": g in ALREADY})
Cd = pd.DataFrame(cand)
Cd["eligible"] = ((Cd.crc_related | Cd.druggable) & Cd.drawable & ~Cd.already_in_fig28)
Cd = Cd.sort_values(["eligible", "delta"], key=lambda s: s.abs() if s.name == "delta" else s,
                    ascending=[False, False])
Cd.to_csv("results/cosmx_selected_genes.csv", index=False)
PICK = Cd[Cd.eligible].head(6)["gene"].tolist()
print("\neligible:", Cd[Cd.eligible].gene.tolist())
print("picked  :", PICK)
print(Cd[Cd.eligible].head(8).round(4).to_string(index=False))

P = pd.read_csv("spatial/cx_pseudobulk.csv.gz")
P = P[P.cellType.isin([NORMAL, TUMOUR])]
pairs, rows = {}, []
for g in PICK:
    sub = P[P.gene == g].pivot_table(index="patient_ID", columns="cellType",
                                     values="mean_norm").dropna(subset=[NORMAL, TUMOUR])
    pairs[g] = sub
    w = wilcoxon(sub[TUMOUR], sub[NORMAL])
    c = Cd[Cd.gene == g].iloc[0]
    rows.append({"gene": g, "n_patients": len(sub),
                 "mean_normal_epi": float(sub[NORMAL].mean()),
                 "mean_CRC": float(sub[TUMOUR].mean()),
                 "delta": float((sub[TUMOUR] - sub[NORMAL]).mean()),
                 "p_wilcoxon": float(w.pvalue),
                 "n_up": int((sub[TUMOUR] > sub[NORMAL]).sum()),
                 "twas_direction": c.twas_direction, "concordant": c.concordant,
                 "is_driver": c.is_driver, "driver_contact": c.driver_contact,
                 "tier": c.tier, "validated": c.validated})
Sel = pd.DataFrame(rows)
Sel.to_csv("results/cosmx_selected_genes_stats.csv", index=False)
print("\n", Sel.round(4).to_string(index=False))

n = len(PICK)
ncol = 3
nrow = 2
fig = plt.figure(figsize=(3.45 * ncol, 4.35 * nrow))
# three rows per block: map, dot plot, then a spacer so the next block's gene
# title never lands on the previous block's x tick labels
gs = fig.add_gridspec(nrow * 3, ncol, height_ratios=[1.0, 0.62, 0.34] * nrow,
                      hspace=0.42, wspace=0.34)
letters = "abcdefghijkl"
MAPS, DOTS = [], []

for k, g in enumerate(PICK):
    rr, cc = divmod(k, ncol)
    z = np.load(SRC[g], allow_pickle=True)
    X = z["x"].astype(float); Y = z["y"].astype(float); Y = Y.max() - Y
    E = (z["expr_q"].astype(np.float32) / float(z["qscale"])) if "expr_q" in z \
        else z["expr"].astype(np.float32)
    e = E[:, [str(x) for x in z["genes"]].index(g)]

    ax = fig.add_subplot(gs[3 * rr, cc])
    neg = e <= 0
    ax.scatter(X[neg], Y[neg], s=0.10, c="#eceae5", linewidths=0, rasterized=True)
    pos = ~neg
    o = np.argsort(e[pos])
    vmax = float(np.quantile(e[pos], 0.98)) if pos.sum() > 50 else 1.0
    sc = ax.scatter(X[pos][o], Y[pos][o], s=0.30, c=e[pos][o], cmap=EXPCM,
                    vmin=0, vmax=max(vmax, 0.2), linewidths=0, rasterized=True)
    ax.set_aspect("equal"); ax.axis("off")
    c = Cd[Cd.gene == g].iloc[0]
    why = []
    if c.is_driver:
        why.append("CRC driver")
    elif c.driver_contact:
        why.append("contacts a driver")
    if c.druggable:
        why.append(c.tier.split(":")[0])
    # the map axes shrink to the tissue aspect, so axes-fraction coordinates are
    # useless for the headings; they are placed after the draw, in figure
    # coordinates anchored on the grid cell (see below)
    MAPS.append((k, rr, cc, ax, g, c, " · ".join(why)))
    cb = fig.colorbar(sc, ax=ax, fraction=0.030, pad=0.005, shrink=0.60,
                      location="left", anchor=(0.0, 0.06))
    cb.ax.tick_params(labelsize=6.0, pad=1); cb.outline.set_visible(False)
    cb.set_label("log$_2$ norm.", fontsize=6.0, labelpad=2)

    ax = fig.add_subplot(gs[3 * rr + 1, cc])
    sub = pairs[g]
    st = Sel[Sel.gene == g].iloc[0]
    for _, r in sub.iterrows():
        up = r[TUMOUR] > r[NORMAL]
        ax.plot([0, 1], [r[NORMAL], r[TUMOUR]], color="#e6b8b8" if up else "#c9c7c1",
                lw=1.0, zorder=2)
    ax.scatter(np.zeros(len(sub)), sub[NORMAL], s=28, color=NORMC, edgecolor="white",
               linewidth=0.6, zorder=4)
    ax.scatter(np.ones(len(sub)), sub[TUMOUR], s=28, color=CRCC, edgecolor="white",
               linewidth=0.6, zorder=4)
    for xx, vals, col in [(0, sub[NORMAL], NORMC), (1, sub[TUMOUR], CRCC)]:
        ax.plot([xx - 0.16, xx + 0.16], [vals.mean()] * 2, color=col, lw=2.2, zorder=5)
    ax.set_xticks([0, 1])
    ax.set_xticklabels(["normal\nepithelium", "carcinoma\ncells"], fontsize=7.4)
    ax.set_xlim(-0.35, 1.35)
    if cc == 0:
        ax.set_ylabel("Mean normalised\nexpression per patient", fontsize=7.8)
    ax.set_title(f"paired Wilcoxon P = {st.p_wilcoxon:.3f}", fontsize=8.4, pad=15,
                 color=GREEN if st.p_wilcoxon < 0.05 else INK2, loc="center")
    ax.text(0.5, 1.012, f"{int(st.n_up)} of {int(st.n_patients)} patients rise",
            transform=ax.transAxes, ha="center", va="bottom", fontsize=7.4, color=INK2)
    DOTS.append((k, rr, cc, ax))
    ax.grid(axis="y"); ax.set_axisbelow(True)

fig.canvas.draw()
for k, rr, cc, ax, g, c, why in MAPS:
    cell = gs[3 * rr, cc].get_position(fig)
    m = ax.get_position()
    fig.text(cell.x0 + 0.019, m.y1 + 0.040, g, fontsize=10.2, fontweight="bold",
             color=INK, va="baseline", ha="left")
    fig.text(cell.x0 + 0.019, m.y1 + 0.024,
             f"CosMx Δ {c.delta:+.3f}, P = {c.cosmx_p:.3f} (n = 7)",
             fontsize=7.2, color=INK2, va="baseline", ha="left")
    fig.text(cell.x0 + 0.019, m.y1 + 0.009, why, fontsize=7,
             color="#4a3aa7", va="baseline", ha="left")
    fig.text(cell.x0 - 0.002, m.y1 + 0.055, letters[k], fontsize=11,
             fontweight="bold", color=INK, va="top", ha="left")
for k, rr, cc, ax in DOTS:
    cell = gs[3 * rr + 1, cc].get_position(fig)
    d = ax.get_position()
    fig.text(cell.x0 - 0.004, d.y1 + 0.034, letters[k + n], fontsize=11,
             fontweight="bold", color=INK, va="top", ha="left")

fig.text(0.5, 1.012,
         "Selected by rule: nominal CosMx P < 0.05, and either a CRC somatic driver or in direct\n"
         "contact with one, or druggable at tier 1–2; ranked by the size of the single-cell\n"
         "difference. Genes already in Figure 28 are excluded.",
         ha="center", fontsize=8.0, color=INK2)
fig.savefig("figures/fig30_cosmx_selected.png", dpi=220, bbox_inches="tight")
fig.savefig("figures/fig30_cosmx_selected.pdf", bbox_inches="tight")
plt.close(fig)

nd = int(Sel.concordant.sum()) if "concordant" in Sel else 0
json.dump({"picked": PICK, "n_eligible": int(Cd.eligible.sum()),
           "n_cosmx_sig": int((B.p < 0.05).sum()),
           "stats": Sel.to_dict("records"),
           "n_concordant": int(Sel.concordant.sum()),
           "n_drivers": int(Sel.is_driver.sum()),
           "n_contact": int(Sel.driver_contact.sum()),
           "n_validated": int(Sel.validated.sum())},
          open("results/cosmx_selected_genes.json", "w"), indent=1, default=str)
print("\nfig30 done")
