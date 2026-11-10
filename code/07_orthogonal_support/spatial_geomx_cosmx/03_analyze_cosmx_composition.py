#!/usr/bin/env python3
"""
CosMx: what changes across the invasive progression, cell type by cell type.

CosMx experiment 2 profiled 1,257,627 QC-passing cells from 11 pseudotimed early CRCs
spanning intramucosal carcinoma (im) -> earliest submucosal invasion (T1.sm1) -> deeper
T1 (T1.sm3). It deliberately profiled tumour core and invasive front, with little
adjacent normal tissue, so it complements the GeoMx normal->adenoma->carcinoma axis
rather than repeating it.

What this script can and cannot do is worth stating up front. The per-cell gene
expression lives only in a Seurat .rds, and no R interpreter is available in this
environment, so per-cell expression of the risk genes cannot be read. The cell tables
are flat text, so cell identity, niche assignment and true slide coordinates are fully
available. This script therefore tests composition and spatial organisation, and the
per-gene expression question is answered from GeoMx instead.

Trend across the three ordered stages is tested per patient with Spearman correlation
between a cell type's proportion and the stage rank, plus Kruskal-Wallis.
"""
import pandas as pd, numpy as np, json
from scipy.stats import spearmanr, kruskal, mannwhitneyu

np.random.seed(21)
STAGE_RANK = {"im": 0, "T1.sm1": 1, "T1.sm3": 2}
PRETTY = {"im": "Intramucosal", "T1.sm1": "T1 sm1", "T1.sm3": "T1 sm3"}

comp = pd.read_csv("spatial/cosmx_composition.tsv", sep="\t")
niche = pd.read_csv("spatial/cosmx_composition_niche.tsv", sep="\t")
cells = pd.read_csv("spatial/cosmx_cells.tsv.gz", sep="\t")
print(f"cells: {len(cells):,}  patients: {comp.patient_ID.nunique()}  "
      f"cell types: {comp.cellType.nunique()}")
print(comp.groupby("tumorStage")["n"].sum().to_string())


def bh(p):
    p = np.asarray(p, dtype=float); o = np.argsort(p); q = np.empty_like(p)
    q[o] = np.minimum(1, np.minimum.accumulate((p[o] * len(p) / (np.arange(len(p)) + 1))[::-1])[::-1])
    return q


def trend(tbl, label_col):
    tot = tbl.groupby(["patient_ID", "tumorStage"])["n"].sum().rename("total").reset_index()
    t = tbl.merge(tot, on=["patient_ID", "tumorStage"])
    t["prop"] = t["n"] / t["total"]
    full = (t.pivot_table(index=["patient_ID", "tumorStage"], columns=label_col,
                          values="prop", fill_value=0.0).reset_index())
    full["rank"] = full["tumorStage"].map(STAGE_RANK)
    rows = []
    for ct in [c for c in full.columns if c not in ("patient_ID", "tumorStage", "rank")]:
        y = full[ct].values
        r = spearmanr(full["rank"].values, y)
        groups = [full.loc[full["tumorStage"] == s, ct].values for s in STAGE_RANK]
        groups = [g for g in groups if len(g) > 0]
        kw = kruskal(*groups) if len(groups) >= 2 and min(len(g) for g in groups) >= 2 else None
        rows.append({label_col: ct, "n_patients": int(full["patient_ID"].nunique()),
                     "mean_prop_im": float(full.loc[full.tumorStage == "im", ct].mean()),
                     "mean_prop_T1sm1": float(full.loc[full.tumorStage == "T1.sm1", ct].mean()),
                     "mean_prop_T1sm3": float(full.loc[full.tumorStage == "T1.sm3", ct].mean()),
                     "spearman_rho": float(r.statistic), "p_trend": float(r.pvalue),
                     "p_kruskal": float(kw.pvalue) if kw else np.nan})
    out = pd.DataFrame(rows)
    out["fdr_trend"] = bh(out["p_trend"])
    return out.sort_values("p_trend"), full


CT, ct_wide = trend(comp, "cellType")
CT.to_csv("results/cosmx_celltype_trend.csv", index=False)
print("\n=== cell-type proportion across im -> T1.sm1 -> T1.sm3 ===")
print(CT.head(18)[["cellType", "mean_prop_im", "mean_prop_T1sm1", "mean_prop_T1sm3",
                   "spearman_rho", "p_trend", "fdr_trend"]].to_string(index=False))

NI, ni_wide = trend(niche, "neighborhood")
NI.to_csv("results/cosmx_niche_trend.csv", index=False)
print("\n=== niche proportion trend (top 8) ===")
print(NI.head(8).to_string(index=False))

# the paper's own headline states, as a positive control on this reanalysis
print("\n=== positive control: the states the paper reports as emerging with invasion ===")
for ct in ["oncofetal", "trophocyte-like CAF", "trophocyte", "MMP CAF", "telocyte",
           "CRC", "normal/undefined epithelium"]:
    r = CT[CT["cellType"] == ct]
    if len(r):
        r = r.iloc[0]
        print(f"  {ct:28s} im {100*r.mean_prop_im:6.2f}% -> sm1 {100*r.mean_prop_T1sm1:6.2f}% "
              f"-> sm3 {100*r.mean_prop_T1sm3:6.2f}%   rho={r.spearman_rho:+.2f}  P={r.p_trend:.3f}")

# ---------------------------------------------------------------- risk-gene link
risk = pd.read_csv("results/risk_genes_harmonized.csv")
n_by_ct = risk.groupby("cell_type")["gene"].nunique()
print("\nTWAS cell types and how many risk genes each carries:")
print(n_by_ct.sort_values(ascending=False).to_string())
print("\nCosMx exp2 cell-type vocabulary:")
print(sorted(comp["cellType"].unique()))
print("\nThese vocabularies do not overlap: CosMx exp2 annotates tumour, oncofetal and "
      "stromal states, not the normal epithelial lineages (absorptive, goblet, stem, "
      "tuft, enteroendocrine) the TWAS was called in. A gene-level, cell-type-matched "
      "test is therefore not possible in this dataset, and none is reported.")

# what IS testable: are the epithelial compartments the risk genes live in expanding?
epi = pd.read_csv("spatial/cosmx_composition_epi.tsv", sep="\t")
EPI, _ = trend(epi, "cellTypeEpithelium")
EPI.to_csv("results/cosmx_epithelium_trend.csv", index=False)
print("\n=== epithelial compartment trend ===")
print(EPI.to_string(index=False))

json.dump({"n_cells": int(len(cells)),
           "n_patients": int(comp.patient_ID.nunique()),
           "stages": {k: int(v) for k, v in comp.groupby("tumorStage")["n"].sum().items()},
           "celltype_trend": CT.to_dict("records"),
           "niche_trend": NI.head(12).to_dict("records"),
           "epithelium_trend": EPI.to_dict("records"),
           "twas_celltypes": {k: int(v) for k, v in n_by_ct.items()},
           "cosmx_celltypes": sorted(comp["cellType"].unique()),
           "celltype_vocabularies_overlap": False},
          open("results/cosmx_summary.json", "w"), indent=1, default=float)
ct_wide.to_csv("results/cosmx_celltype_proportions.csv", index=False)
print("\nDone.")
