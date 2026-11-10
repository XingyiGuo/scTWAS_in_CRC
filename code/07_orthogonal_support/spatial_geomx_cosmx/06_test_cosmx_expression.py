#!/usr/bin/env python3
"""
What the CosMx single cells say about the risk genes.

Input is the pseudobulk table from 05_summarize_cosmx_expression.py: for every patient, cell type and
gene, the mean normalised expression over that patient's cells of that type. All tests
are across patients, never across cells - 1.26 M cells come from 11 people, and treating
them as 1.26 M independent observations would manufacture significance.

Three questions:

  A. Does the gene change along the invasive sequence (intramucosal -> T1 sm1 -> T1 sm3)
     within a given cell type? Spearman across patients, per cell type.

  B. Within the same patient, is the gene expressed differently in carcinoma cells than
     in the normal/undefined epithelium beside them? This is the single-cell version of
     the oncogene / tumour-suppressor claim, and it is paired, so patient-to-patient
     variation cancels.

  C. Do the genes validated in GeoMx move the same way here? CosMx is a different
     platform, a different gene panel and 11 different patients, so agreement is
     replication rather than reanalysis.
"""
import pandas as pd, numpy as np, json
from scipy.stats import (spearmanr, wilcoxon, mannwhitneyu, binomtest, kruskal,
                         fisher_exact)

np.random.seed(26)
STAGE_RANK = {"im": 0, "T1.sm1": 1, "T1.sm3": 2}
EPI = ["CRC", "normal/undefined epithelium", "oncofetal"]
MIN_CELLS = 50            # a patient x cell type group must have this many cells to count


def bh(p):
    p = np.asarray(p, float)
    ok = np.isfinite(p)
    q = np.full(p.shape, np.nan)
    pp = p[ok]
    o = np.argsort(pp)
    r = np.empty_like(pp)
    r[o] = np.minimum(1, np.minimum.accumulate(
        (pp[o] * len(pp) / (np.arange(len(pp)) + 1))[::-1])[::-1])
    q[ok] = r
    return q


pb = pd.read_csv("spatial/cx_pseudobulk.csv.gz")
pb = pb[pb["n_cells"] >= MIN_CELLS].copy()
pb["rank"] = pb["tumorStage"].map(STAGE_RANK)
allm = pd.read_csv("results/validated_risk_genes_all.csv")
dirmap = dict(zip(allm["gene"], allm["twas_direction"]))
genes = sorted(pb["gene"].unique())
print(f"{len(genes)} risk genes on the CosMx panel; "
      f"{pb['patient_ID'].nunique()} patients; "
      f"{pb['cellType'].nunique()} cell types with >= {MIN_CELLS} cells")

# ------------------------------------------------------------------ A. stage trend
rows = []
for (ct, g), d in pb.groupby(["cellType", "gene"]):
    if d["patient_ID"].nunique() < 6 or d["rank"].nunique() < 3:
        continue
    r = spearmanr(d["rank"], d["mean_norm"])
    grp = [d.loc[d["rank"] == k, "mean_norm"].values for k in (0, 1, 2)]
    grp = [x for x in grp if len(x)]
    kw = kruskal(*grp) if len(grp) >= 2 and max(len(x) for x in grp) >= 2 else None
    rows.append({"cellType": ct, "gene": g, "n_patients": int(d["patient_ID"].nunique()),
                 "mean_im": float(d.loc[d["rank"] == 0, "mean_norm"].mean()),
                 "mean_sm1": float(d.loc[d["rank"] == 1, "mean_norm"].mean()),
                 "mean_sm3": float(d.loc[d["rank"] == 2, "mean_norm"].mean()),
                 "rho": float(r.statistic), "p": float(r.pvalue),
                 "p_kruskal": float(kw.pvalue) if kw else np.nan})
A = pd.DataFrame(rows)
A["twas_direction"] = A["gene"].map(dirmap)
A["fdr"] = bh(A["p"])
A = A.sort_values("p")
A.to_csv("results/cosmx_gene_stage_trend.csv", index=False)
print(f"\n=== A. stage trend: {len(A)} gene x cell-type tests, "
      f"{int((A['fdr'] < 0.05).sum())} at FDR < 0.05 ===")
print(A.head(12)[["gene", "cellType", "mean_im", "mean_sm1", "mean_sm3", "rho", "p", "fdr"]]
      .to_string(index=False))

# ------------------------------------------------------------------ B. carcinoma vs normal epithelium, paired
piv = (pb[pb["cellType"].isin(["CRC", "normal/undefined epithelium"])]
       .pivot_table(index=["gene", "patient_ID"], columns="cellType", values="mean_norm"))
piv = piv.dropna()
rows = []
for g, d in piv.groupby(level=0):
    if len(d) < 6:
        continue
    a = d["CRC"].values
    b = d["normal/undefined epithelium"].values
    w = wilcoxon(a, b)
    diff = float(np.mean(a - b))
    exp_sign = 1.0 if dirmap.get(g) == "risk-increasing" else -1.0
    rows.append({"gene": g, "n_patients": int(len(d)),
                 "mean_CRC": float(a.mean()), "mean_normal_epi": float(b.mean()),
                 "delta": diff, "p": float(w.pvalue),
                 "twas_direction": dirmap.get(g),
                 "concordant": bool(diff * exp_sign > 0)})
B = pd.DataFrame(rows)
B["fdr"] = bh(B["p"])
B = B.sort_values("p")
nsig = int((B["fdr"] < 0.05).sum())
print(f"\n=== B. carcinoma vs adjacent normal epithelium, paired within patient ===")
print(f"{len(B)} genes tested, {nsig} at FDR < 0.05")
print(B.head(15)[["gene", "twas_direction", "mean_CRC", "mean_normal_epi", "delta",
                  "p", "fdr", "concordant"]].to_string(index=False))

# With 7-11 patients the paired Wilcoxon cannot go below 2^-(n-1), so no gene can reach
# FDR < 0.05 across 121 tests however large the effect. The informative quantity is not
# the p-value of any one gene but the DIRECTION the genes take, tested across genes.
n_pat_med = int(B["n_patients"].median())
p_floor = 2.0 ** (-(n_pat_med - 1))
sigB = B[B["p"] < 0.05]
conc_nom = int(sigB["concordant"].sum())
bt = binomtest(conc_nom, len(sigB), 0.5) if len(sigB) else None
print(f"\npaired-Wilcoxon p floor at n = {n_pat_med} patients: {p_floor:.4f}; "
      f"no single gene can clear FDR < 0.05 across {len(B)} tests.")
print(f"nominally significant (P < 0.05): {len(sigB)}; of those {conc_nom} "
      f"({100*conc_nom/max(len(sigB),1):.0f} %) move the way the TWAS predicts"
      + (f", binomial P = {bt.pvalue:.3g}" if bt else ""))

geomx_supported = allm.loc[allm["nor_core_supported"], "gene"]
B["geomx_validated"] = B["gene"].isin(geomx_supported)
conc_all = int(B["concordant"].sum())
bt_all = binomtest(conc_all, len(B), 0.5)
Bv = B[B["geomx_validated"]]
conc_v = int(Bv["concordant"].sum())
bt_v = binomtest(conc_v, len(Bv), 0.5) if len(Bv) else None
print("\nTWAS-direction agreement of the single-cell carcinoma/normal contrast:")
print(f"  all {len(B)} risk genes on the panel: {conc_all} "
      f"({100*conc_all/len(B):.0f} %), binomial P = {bt_all.pvalue:.3g}")
if bt_v:
    print(f"  the {len(Bv)} GeoMx-validated genes:  {conc_v} "
          f"({100*conc_v/len(Bv):.0f} %), binomial P = {bt_v.pvalue:.3g}")
B.to_csv("results/cosmx_crc_vs_normal_epi.csv", index=False)

# ------------------------------------------------------------------ C. replication of the GeoMx-validated set
gm = allm.set_index("gene")
shared = [g for g in B["gene"] if g in gm.index]
rep = B[B["gene"].isin(shared)].copy()
rep["geomx_log2FC"] = rep["gene"].map(gm["log2FC_NOR_CORE"])
rep = rep.dropna(subset=["geomx_log2FC"])
r_all = spearmanr(rep["geomx_log2FC"], rep["delta"])
sign_all = float(np.mean(np.sign(rep["geomx_log2FC"]) == np.sign(rep["delta"])))
sub = rep[rep["geomx_validated"]]
r_val = spearmanr(sub["geomx_log2FC"], sub["delta"]) if len(sub) > 3 else None
sign_val = float(np.mean(np.sign(sub["geomx_log2FC"]) == np.sign(sub["delta"]))) if len(sub) else np.nan
bt_val = binomtest(int(round(sign_val * len(sub))), len(sub), 0.5) if len(sub) else None
print(f"\n=== C. CosMx vs GeoMx, independent patients and platforms ===")
print(f"all {len(rep)} risk genes on both: Spearman rho = {r_all.statistic:.3f} "
      f"(P = {r_all.pvalue:.3g}), sign agreement {100*sign_all:.0f} %")
if r_val is not None:
    print(f"the {len(sub)} GeoMx-validated genes: rho = {r_val.statistic:.3f} "
          f"(P = {r_val.pvalue:.3g}), sign agreement {100*sign_val:.0f} % "
          f"(binomial P = {bt_val.pvalue:.3g})")
rep.to_csv("results/cosmx_geomx_replication.csv", index=False)

# ------------------------------------------------------------------ cell-type specificity of the signal
spec = []
for ct, d in pb.groupby("cellType"):
    sub2 = A[A["cellType"] == ct]
    if not len(sub2):
        continue
    spec.append({"cellType": ct, "n_tests": len(sub2),
                 "n_fdr05": int((sub2["fdr"] < 0.05).sum()),
                 "median_abs_rho": float(sub2["rho"].abs().median()),
                 "n_cells_total": int(pb.loc[pb["cellType"] == ct, "n_cells"].sum()
                                      / max(pb["gene"].nunique(), 1))})
S = pd.DataFrame(spec).sort_values("n_fdr05", ascending=False)
S.to_csv("results/cosmx_celltype_specificity.csv", index=False)
print("\n=== where the stage signal lives ===")
print(S.head(12).to_string(index=False))

json.dump({
    "n_genes_on_panel": len(genes),
    "n_patients": int(pb["patient_ID"].nunique()),
    "min_cells_per_group": MIN_CELLS,
    "stage_trend": {"n_tests": int(len(A)), "n_fdr05": int((A["fdr"] < 0.05).sum()),
                    "top": A.head(15).to_dict("records")},
    "crc_vs_normal": {"n_genes": int(len(B)), "n_fdr05": nsig,
                      "n_patients_median": n_pat_med, "p_floor": float(p_floor),
                      "n_nominal": int(len(sigB)), "n_concordant_of_nominal": conc_nom,
                      "frac_concordant_of_nominal": float(conc_nom / max(len(sigB), 1)),
                      "p_binomial_nominal": float(bt.pvalue) if bt else None,
                      "n_concordant_all": conc_all,
                      "frac_concordant_all": float(conc_all / len(B)),
                      "p_binomial_all": float(bt_all.pvalue),
                      "n_validated_on_panel": int(len(Bv)),
                      "n_concordant_validated": conc_v,
                      "frac_concordant_validated": float(conc_v / max(len(Bv), 1)),
                      "p_binomial_validated": float(bt_v.pvalue) if bt_v else None,
                      "top": B.head(20).to_dict("records")},
    "replication": {"n_shared": int(len(rep)), "rho_all": float(r_all.statistic),
                    "p_all": float(r_all.pvalue), "sign_agreement_all": sign_all,
                    "n_validated": int(len(sub)),
                    "rho_validated": float(r_val.statistic) if r_val is not None else None,
                    "p_validated": float(r_val.pvalue) if r_val is not None else None,
                    "sign_agreement_validated": sign_val,
                    "p_binomial_validated": float(bt_val.pvalue) if bt_val else None},
    "celltype_specificity": S.to_dict("records"),
}, open("results/cosmx_percell_summary.json", "w"), indent=1, default=float)
print("\nDone.")
