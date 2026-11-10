"""Integration of the 51-TF risk-gene network with the scTWAS EnrichR pathway results.

Produces the tables in ../table and Figure 4e,f in ../figure.

Inputs
------
W:/VUMC_projects/scTWAS/TF_genes_pathways/NL_ADSER_cell_genes_EnrichR.xlsx
L:/CRC51TF-scTWASgenes/TF_network_analysis/data/S1_TF_gene_correlations_all.csv
L:/CRC51TF-scTWASgenes/TF_network_analysis/LQ_review/raw_data/{NL,ADSER}_lead_variant_per_celltype_gene_with_occupying_TFs.csv

Join convention: EnrichR Dataset NL -> tissue 'Normal', ADSER -> tissue 'Adenoma';
match on cell type. A pathway member gene is scored on two layers - TF co-expression
(S1 sig == r > 0.5 and FDR < 0.05 in the SAME cell type) and TF occupancy of that
gene's GWAS lead variant. Genes with both are 'dual evidence'.
"""
import os
import re
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from scipy.stats import hypergeom
from statsmodels.stats.multitest import multipletests

XL = os.environ.get(
    "TF_ENRICHR_XLSX", r"W:\VUMC_projects\scTWAS\TF_genes_pathways\NL_ADSER_cell_genes_EnrichR.xlsx"
)
NET = os.environ.get("TF_NETWORK_DIR", r"L:\CRC51TF-scTWASgenes\TF_network_analysis")
OUT = os.environ.get("TF_PATHWAY_OUT", r"L:\CRC51TF-scTWASgenes\Figure4_TF_pathways")
DS2TIS = {"NL": "Normal", "ADSER": "Adenoma"}

# ---------------------------------------------------------------- load
enr = pd.read_excel(XL)
s1 = pd.read_csv(os.path.join(NET, "data", "S1_TF_gene_correlations_all.csv"))
rd = os.path.join(NET, "LQ_review", "raw_data")
occ = pd.concat(
    [pd.read_csv(os.path.join(rd, "NL_lead_variant_per_celltype_gene_with_occupying_TFs.csv")).assign(Dataset="NL"),
     pd.read_csv(os.path.join(rd, "ADSER_lead_variant_per_celltype_gene_with_occupying_TFs.csv")).assign(Dataset="ADSER")],
    ignore_index=True)
occ["TFset"] = occ.TFs_on_lead_variant.fillna("").apply(lambda s: {x for x in str(s).split(";") if x})
occ_map = {(r.Dataset, r.Cells, r.gene_name): r.TFset for r in occ.itertuples()}

sig = s1[s1.sig]
uni = s1.groupby(["tissue", "cell_type"]).gene.apply(set).to_dict()          # genes tested per cell type
linked = sig.groupby(["tissue", "cell_type"]).gene.apply(set).to_dict()      # genes with >=1 TF link
tf_targets = sig.groupby(["tissue", "cell_type", "TF"]).gene.apply(set).to_dict()
tfs_by_gene = sig.groupby(["tissue", "cell_type", "gene"]).TF.apply(set).to_dict()

# ------------------------------------------------- per-term and per-gene scoring
rows, gene_rows, tp_rows = [], [], []
for r in enr.rename(columns={"Adjusted.P.value": "EnrichR_adjP"}).itertuples(index=False):
    tis, ct = DS2TIS[r.Dataset], r.CellType
    pg = [g for g in str(r.Genes).split(";") if g]
    U, L = uni.get((tis, ct), set()), linked.get((tis, ct), set())
    pg_in = [g for g in pg if g in U]
    lk = [g for g in pg if g in L]
    occ_g = [g for g in pg if occ_map.get((r.Dataset, ct, g), set())]
    dual = []
    for g in pg:
        tfs_corr = tfs_by_gene.get((tis, ct, g), set())
        tfs_occ = occ_map.get((r.Dataset, ct, g), set())
        if tfs_corr & tfs_occ:
            dual.append(g)
        gene_rows.append(dict(Dataset=r.Dataset, CellType=ct, Source=r.Source, Term=r.Term, gene=g,
                              in_universe=g in U, n_TF_corr=len(tfs_corr), TFs_corr=";".join(sorted(tfs_corr)),
                              n_TF_occ=len(tfs_occ), TFs_occ=";".join(sorted(tfs_occ)),
                              TFs_dual=";".join(sorted(tfs_corr & tfs_occ))))
    M, n, N, k = len(U), len(L), len(pg_in), len(set(lk) & U)
    allTF = sorted({t for g in pg for t in tfs_by_gene.get((tis, ct, g), set())})
    rows.append(dict(Dataset=r.Dataset, CellType=ct, Source=r.Source, Term=r.Term,
                     EnrichR_adjP=r.EnrichR_adjP, n_pathway_genes=len(pg), n_in_universe=N,
                     n_TF_linked=len(lk), frac_TF_linked=round(len(lk) / len(pg), 3),
                     n_leadvar_occupied=len(occ_g), n_dual_evidence=len(dual),
                     bg_frac_TF_linked=round(n / M, 3) if M else np.nan,
                     hyper_p=hypergeom.sf(k - 1, M, n, N) if M and N else np.nan,
                     n_TFs=len(allTF), TFs=";".join(allTF),
                     genes_TF_linked=";".join(sorted(set(lk))), genes_dual=";".join(sorted(set(dual)))))
    for tf in allTF:                                   # TF x pathway over-representation
        T = tf_targets.get((tis, ct, tf), set())
        kk = len(set(pg) & T)
        if kk >= 2:
            tp_rows.append(dict(Dataset=r.Dataset, CellType=ct, Source=r.Source, Term=r.Term, TF=tf,
                                n_shared=kk, n_TF_targets=len(T), n_pathway_genes=N,
                                p=hypergeom.sf(kk - 1, len(U), len(T), N),
                                shared_genes=";".join(sorted(set(pg) & T))))

path_tf = pd.DataFrame(rows)
path_tf["hyper_q"] = multipletests(path_tf.hyper_p.fillna(1), method="fdr_bh")[1]
genes_df = pd.DataFrame(gene_rows)
tfp = pd.DataFrame(tp_rows, columns=[
    "Dataset", "CellType", "Source", "Term", "TF", "n_shared",
    "n_TF_targets", "n_pathway_genes", "p", "shared_genes",
])
tfp["q"] = multipletests(tfp.p, method="fdr_bh")[1] if len(tfp) else []

# ---------------------------------------- collapse terms with identical gene sets
g = genes_df.assign(dual=lambda d: d.TFs_dual != "", corr=lambda d: d.n_TF_corr > 0)
sub = path_tf[path_tf.n_pathway_genes >= 3].copy()
key = enr.set_index(["Dataset", "CellType", "Term"]).Genes.to_dict()
sub["members"] = [frozenset(key[(r.Dataset, r.CellType, r.Term)].split(";")) for r in sub.itertuples()]
coll = []
for (ds, ct, ms), d in sub.groupby(["Dataset", "CellType", "members"]):
    d = d.sort_values("EnrichR_adjP")
    lab = d.Term.iloc[0] + (f" (+{len(d) - 1})" if len(d) > 1 else "")
    gg = g[(g.Dataset == ds) & (g.CellType == ct) & (g.gene.isin(ms))].drop_duplicates("gene")
    coll.append(dict(Dataset=ds, CellType=ct, label=lab, n_terms=len(d), n_genes=len(ms),
                     n_dual=int(gg.dual.sum()), n_corr_only=int((gg["corr"] & ~gg.dual).sum()),
                     n_none=int((~gg["corr"]).sum()), minq=d.EnrichR_adjP.min(),
                     TFs=";".join(sorted({t for s in gg.TFs_corr for t in s.split(";") if t}))))
coll = pd.DataFrame(coll).sort_values(["Dataset", "CellType", "n_dual", "n_corr_only"],
                                      ascending=[True, True, False, False])

os.makedirs(os.path.join(OUT, "table"), exist_ok=True)
path_tf.sort_values(["Dataset", "CellType", "EnrichR_adjP"]).to_csv(
    os.path.join(OUT, "table", "T1_pathway_TF_regulation.csv"), index=False)
genes_df.to_csv(os.path.join(OUT, "table", "T2_pathway_gene_TF_evidence.csv"), index=False)
coll.to_csv(os.path.join(OUT, "table", "T3_pathway_modules_collapsed.csv"), index=False)
tfp.sort_values("p").to_csv(os.path.join(OUT, "table", "T4_TF_pathway_pairs_hypergeom.csv"), index=False)

# --------------------------------------------------------- Figure 4e,f
e = enr.copy()
e["members"] = e.Genes.apply(lambda s: frozenset(s.split(";")))
e["nG"] = e.members.apply(len)
e["n_TFlinked"] = e.apply(lambda r: len(r.members & linked.get((DS2TIS[r.Dataset], r.CellType), set())), axis=1)
e["TFs"] = e.apply(lambda r: ";".join(sorted(set().union(
    *[tfs_by_gene.get((DS2TIS[r.Dataset], r.CellType, x), set()) for x in r.members]))), axis=1)

# representative term per gene set: drop non-colon disease names, then prefer a
# canonical process name among terms within 1.5x of the group's best adjusted P
BAD = re.compile(r"allograft|graft-versus|diabetes|thyroid|myocarditis|asthma|leishman|malaria|toxoplasm|"
                 r"amoebiasis|cardiomyopathy|alzheimer|parkinson|huntington|prion|amyotrophic|neurodegener|"
                 r"cocaine|alcohol|diabetic|hepatitis|tuberculos|salmonella|shigell|pertussis|staphyl|"
                 r"influenza|measles|COVID", re.I)
PRIO = re.compile(r"antigen process|interferon|wnt|beta-catenin|cell cycle|oxidative phosphor|senescence|p53|"
                  r"hippo|homologous recomb|ferroptosis|gastrulation|glucuronid|ascorbate|mitotic|"
                  r"reactive oxygen|epithelial mesenchymal|natural killer", re.I)


def pick(d):
    ok = d[~d.Term.str.contains(BAD)]
    ok = ok if len(ok) else d
    near = ok[ok["Adjusted.P.value"] <= ok["Adjusted.P.value"].min() * 1.5]
    pr = near[near.Term.str.contains(PRIO)]
    return (pr if len(pr) else ok).sort_values("Adjusted.P.value").iloc[0]


rows = []
for (ds, ct), d in e.groupby(["Dataset", "CellType"]):
    groups = sorted(d.groupby("members", sort=False),
                    key=lambda kv: (kv[1]["Adjusted.P.value"].min(), -kv[1].nG.max()))
    for _, grp in groups[:5]:
        rep = pick(grp)
        rows.append(dict(Dataset=ds, CellType=ct, Term=rep.Term, n_similar=len(grp) - 1, nG=rep.nG,
                         q=rep["Adjusted.P.value"], n_TFlinked=rep.n_TFlinked, TFs=rep.TFs, Genes=rep.Genes))
P = pd.DataFrame(rows)
P["nlq"] = -np.log10(P.q)
CT_ORDER = ["STM", "TAC", "ABS", "GOB", "EE", "TUF", "CT", "ASC", "SSC"]
PAL = dict(zip(CT_ORDER, [
    "#4C72B0", "#DD8452", "#55A868", "#C44E52", "#8172B3", "#937860", "#DA8BC3",
    "#637939", "#3182BD",
]))
P["cto"] = P.CellType.map({c: i for i, c in enumerate(CT_ORDER)})
P = P.sort_values(["Dataset", "cto", "q"])

# display names: the Enrichr term verbatim is too long for a y tick label, so each
# plotted term gets a short label. The full term stays in the source-data table.
SHORT = {
    "Ascorbate and aldarate metabolism": "Ascorbate/aldarate metabolism",
    "Pentose and glucuronate interconversions": "Pentose\u2013glucuronate interconversion",
    "UV Response Dn": "UV response (down)",
    "Hippo signaling pathway": "Hippo signalling",
    "Human papillomavirus infection": "HPV infection",
    "Human T-cell leukemia virus 1 infection": "HTLV-1 infection",
    "Kaposi sarcoma-associated herpesvirus infection": "KSHV infection",
    "Epithelial Mesenchymal Transition": "EMT",
    "Androgen Response": "Androgen response",
    "Endosomal Vacuolar Pathway": "Endosomal/vacuolar MHC-I pathway",
    "Immunoregulatory Interactions Between a Lymphoid and a non-Lymphoid Cell":
        "Lymphoid\u2013non-lymphoid immunoregulation",
    "Antigen processing and presentation": "Antigen processing & presentation",
    "Degradation of Beta-Catenin by the Destruction Complex": "\u03b2-catenin destruction complex",
    "Beta-catenin Phosphorylation Cascade": "\u03b2-catenin phosphorylation cascade",
    "MASTL Facilitates Mitotic Progression": "MASTL mitotic progression",
    "Inhibition of Replication Initiation of Damaged DNA by RB1 E2F1": "RB1\u2013E2F1 replication block",
    "Mitotic Prophase": "Mitotic prophase",
    "Mitotic Spindle": "Mitotic spindle",
    "Oxidative phosphorylation": "Oxidative phosphorylation",
    "Complex III Assembly": "Complex III assembly",
    "Estrogen Response Late": "Estrogen response (late)",
    "XBP1(S) Activates Chaperone Genes": "XBP1(S) chaperone activation",
    "Homologous recombination": "Homologous recombination",
    "Cellular senescence": "Cellular senescence",
    "Viral carcinogenesis": "Viral carcinogenesis",
    "Endocytosis": "Endocytosis",
    "Ferroptosis": "Ferroptosis",
    "Gastrulation": "Gastrulation",
}
P["short_label"] = P.Term.map(SHORT).fillna(P.Term)

mpl.rcParams.update({"font.size": 8, "axes.titlesize": 8, "axes.labelsize": 8,
                     "axes.spines.top": False, "axes.spines.right": False, "figure.dpi": 110})
DARK, GREY = "#111111", "#8a8a8a"
fig, axes = plt.subplots(1, 2, figsize=(11.4, 5.35), sharex=True)
xmax = P.nlq.max() * 1.20
for ax, ds, ttl in zip(axes, ["NL", "ADSER"], ["Normal colon", "Polyps (adenoma / serrated)"]):
    d = P[P.Dataset == ds].iloc[::-1].reset_index(drop=True)
    yy = np.arange(len(d))
    tf = d.n_TFlinked > 0
    for i, r in d.iterrows():
        ax.plot([0, r.nlq], [i, i], color="#cccccc", lw=.8, zorder=1)
    # terms with TF evidence: tinted band, heavy dark ring, bold dark label, n-linked call-out
    for i in yy[tf.values]:
        ax.axhspan(i - .45, i + .45, color="#f2f2f2", zorder=0)
    ax.scatter(d.nlq, yy, s=22 * d.nG + 18, c=[PAL[c] for c in d.CellType],
               edgecolors=[DARK if t else "#ffffff" for t in tf],
               linewidths=[1.4 if t else .5 for t in tf], zorder=3)
    for i, r in d.iterrows():
        if r.n_TFlinked:
            ax.text(r.nlq + 0.10 + 0.0016 * (22 * r.nG + 18), i, f"{r.n_TFlinked}/{r.nG}",
                    va="center", ha="left", fontsize=6, color=DARK, zorder=4)
    ax.set_yticks(yy)
    # n_similar (further Enrichr terms with an identical risk-gene set) is kept in the
    # source-data table only - printing it on the axis reads as part of the term name
    ax.set_yticklabels(list(d.short_label), fontsize=6.5)
    for lbl, t in zip(ax.get_yticklabels(), tf):
        lbl.set_color(DARK if t else GREY)
        lbl.set_fontweight("bold" if t else "normal")
    ax.axvline(-np.log10(0.05), color="#999999", ls="--", lw=.7, zorder=2)
    ax.set_xlim(0, xmax)
    ax.set_ylim(-0.8, len(d) - 0.2)
    ax.set_xlabel("$-\\log_{10}$ adjusted $P$ (Enrichr)")
    ax.set_title(ttl, loc="left")
    for ct, idx in d.groupby("CellType", sort=False).groups.items():
        idx = list(idx)
        ax.text(1.005, (min(idx) + max(idx)) / 2, ct, transform=ax.get_yaxis_transform(),
                ha="left", va="center", fontsize=7, color=PAL[ct], fontweight="bold")
        if min(idx) > 0:
            ax.axhline(min(idx) - 0.5, color="#e0e0e0", lw=.6)
axes[0].text(-np.log10(0.05) + 0.06, -0.62, "FDR 5%", fontsize=6, color="#777777")
# one figure-level key: it describes both panels, so it sits below them rather than inside g
h = [plt.scatter([], [], s=22 * n + 18, c="#bbbbbb", edgecolors="#ffffff", lw=.5) for n in (1, 3, 5, 7)]
fig.legend(h, ["1", "3", "5", "7"], title="risk genes in term", frameon=False, fontsize=6.5,
           title_fontsize=6.5, loc="lower center", ncol=4, handletextpad=.2, columnspacing=.8,
           borderpad=.2, bbox_to_anchor=(0.26, 0.0))
h2 = [plt.scatter([], [], s=70, c="#bbbbbb", edgecolors=DARK, lw=1.4),
      plt.scatter([], [], s=70, c="#bbbbbb", edgecolors="#ffffff", lw=.5)]
fig.legend(h2, ["TF-regulated term (genes co-expressed with a key TF / genes in term)",
                "no TF-co-expressed gene"], frameon=False, fontsize=6.5, ncol=2,
           loc="lower center", bbox_to_anchor=(0.65, 0.012), handletextpad=.2,
           borderpad=.2, columnspacing=1.6)
for ax, lt in zip(axes, "fg"):
    ax.text(-0.02, 1.04, lt, transform=ax.transAxes, fontsize=11, fontweight="bold", ha="right", va="bottom")
fig.subplots_adjust(left=0.215, right=0.955, top=0.935, bottom=0.155, wspace=0.70)

os.makedirs(os.path.join(OUT, "figure"), exist_ok=True)
# network mounts reject PIL's "w+b" open mode, so hand savefig an explicit binary handle
with open(os.path.join(OUT, "figure", "Figure4FG_pathway_enrichment.png"), "wb") as fh:
    fig.savefig(fh, format="png", dpi=300)
with open(os.path.join(OUT, "figure", "Figure4FG_pathway_enrichment.pdf"), "wb") as fh:
    fig.savefig(fh, format="pdf")
P.to_csv(os.path.join(OUT, "figure", "Figure4FG_sourcedata.csv"), index=False)
