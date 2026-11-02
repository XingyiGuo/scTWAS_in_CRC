#!/usr/bin/env python3
"""
37 - The siRNA screen, cell line by cell line, for every risk gene it covers.

Companion to script 35, which did this for CRISPR. The siRNA panel carries 333 genes,
62 of them S1 risk genes - nearly double the CRISPR overlap - so this is the larger
of the two per-line views and the one that covers most of the risk list.

SW480 is kept here but flagged, not silently dropped. Its siRNA plates fail quality
control (Z' = -0.14, -0.08, +0.21, mean -0.00, against 0.32 for HCT116 and 0.30 for
RKO), which is why section 2.15 pooled only HCT116 and RKO. A failing Z' does not mean
the plate is empty, it means the separation between the killing control and the
non-targeting control is no larger than their combined spread - so a call made in SW480
alone is uninterpretable, and the figure marks the column accordingly.

Three significance criteria are carried for every gene x line, because the previous
round showed that two of them are not nested and readers assumed they were:

  emp_nominal   P < 0.05 from the empirical null, uncorrected
  supported     FDR < 0.05 from the same empirical null, corrected across the panel
                within the line, plus all three replicates agreeing in sign
  nominal       P < 0.05 from a Welch t-test of the gene's three replicate values
                (as % of the plate's non-targeting median) against that line's
                non-targeting wells

The first two are nested by construction - the second is always a subset of the first.
The third is a different test with a different variance estimate and is nested in
neither. Reporting all three makes the relationship visible instead of leaving the
reader to infer one that does not hold.

Outputs
  results/sirna_by_cellline.csv        all 333 genes x 3 lines
  results/sirna_risk_by_cellline.csv   the 62 S1 risk genes, wide
  results/sirna_by_cellline.json       counts for the report and the figure
"""
import pandas as pd, numpy as np, json, sys, os
from scipy.stats import (fisher_exact, ttest_ind, ttest_1samp, ttest_rel,
                         mannwhitneyu)

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _ttest_support import (support_call, gene_call, binomial_concordance,
                            SUPPORTS, OPPOSES, NS)

ALIAS = {"ACAS2": "ACSS2"}
LINES = ["HCT116", "RKO", "SW480"]
N_NULL = 20000
QC_FAIL = {"SW480"}          # mean Z' <= 0 in this assay


def bh(p):
    p = np.asarray(p, float)
    o = np.argsort(p)
    q = np.empty_like(p)
    q[o] = np.minimum(1, np.minimum.accumulate(
        (p[o] * len(p) / (np.arange(len(p)) + 1))[::-1])[::-1])
    return q


S = pd.read_csv("results/func_sirna_long.csv")
S["gene"] = S["gene"].replace(ALIAS)
risk = set(pd.read_csv("results/risk_genes_harmonized.csv")["gene"])
allm = pd.read_csv("results/validated_risk_genes_all.csv").set_index("gene")
valset = set(pd.read_csv("results/validated_priority.csv")["gene"])
QC = pd.read_csv("results/func_qc.csv")
QC = QC[QC.assay == "siRNA"]

rows = []
rng = np.random.default_rng(37)
for line in LINES:
    d = S[S.cell_line == line].copy()
    zs = []
    for rep, dd in d.groupby("replicate"):
        smp = dd.loc[dd.kind == "sample", "value"]
        nt = dd.loc[dd.kind == "negative", "value"]
        med, ntmed = smp.median(), nt.median()
        mad = (smp - med).abs().median() * 1.4826 or smp.std(ddof=1)
        e = dd.copy()
        e["z"] = (e["value"] - med) / mad
        e["z_noise"] = (e["value"] - ntmed) / mad
        e["pct_NT"] = 100 * e["value"] / ntmed
        zs.append(e)
    Z = pd.concat(zs, ignore_index=True)
    nt_by_plate = [g["z_noise"].to_numpy(float) - g["z_noise"].median()
                   for _, g in Z[Z.kind == "negative"].groupby("replicate")]
    null = np.sort(np.mean([rng.choice(v, N_NULL, replace=True) for v in nt_by_plate], axis=0))
    zp = QC.loc[QC.cell_line == line, "z_prime"].mean()
    print(f"{line}: {len(nt_by_plate)} plates, {sum(len(v) for v in nt_by_plate)} NT wells, "
          f"null sd {null.std():.3f}, mean Z' {zp:+.3f}"
          f"{'  << fails QC' if line in QC_FAIL else ''}")

    nt_pct = Z.loc[Z.kind == "negative", "pct_NT"].to_numpy(float)

    # --- the paired design, as in the CRISPR screen: each replicate plate's twelve
    # non-targeting wells are averaged into that plate's control value, a gene has one
    # (gene, control) pair per plate, and the three pairs are tested together.
    nt_mean_by_rep = Z[Z.kind == "negative"].groupby("replicate")["value"].mean()
    reps = sorted(nt_mean_by_rep.index)
    print("  plate control means: "
          + ", ".join(f"{r} {nt_mean_by_rep[r]:.1f}" for r in reps))

    out = []
    for g, dd in Z[Z.kind == "sample"].groupby("gene"):
        z = dd["z"].to_numpy(float)
        pct = dd["pct_NT"].to_numpy(float)
        stat = float(np.mean(z))
        p = (int(np.sum(np.abs(null) >= abs(stat))) + 1) / (len(null) + 1)
        tt = ttest_ind(pct, nt_pct, equal_var=False) if len(pct) > 1 else None
        t1 = ttest_1samp(pct, 100.0) if len(pct) > 1 else None

        sv = dd.set_index("replicate")["value"]
        gv = np.array([sv[r] for r in reps], float)
        cv = np.array([nt_mean_by_rep[r] for r in reps], float)
        diff = gv - cv
        tp = ttest_rel(gv, cv) if len(gv) > 1 else None

        out.append({"cell_line": line, "gene": g, "n_rep": len(z), "mean_z": stat,
                    "pct_NT": float(np.mean(pct)),
                    "nt_mean_pct": float(np.mean(nt_pct)),
                    "delta_pct": float(np.mean(pct) - np.mean(nt_pct)),
                    "value_mean": float(gv.mean()), "control_mean": float(cv.mean()),
                    "delta_paired": float(diff.mean()),
                    "sd_diff": float(diff.std(ddof=1)) if len(diff) > 1 else np.nan,
                    "t_paired": float(tp.statistic) if tp else np.nan,
                    "p_paired": float(tp.pvalue) if tp else np.nan,
                    "diffs": ";".join(f"{d:+.2f}" for d in diff),
                    "pct_of_control": float(100 * gv.mean() / cv.mean()),
                    "sd_pct_NT": float(np.std(pct, ddof=1)) if len(pct) > 1 else np.nan,
                    "p": p,
                    "p_ttest_vs_NT": float(tt.pvalue) if tt else np.nan,
                    "p_ttest_vs_100": float(t1.pvalue) if t1 else np.nan,
                    "reps_same_sign": bool(np.all(np.sign(z) == np.sign(z[0]))),
                    "min_z": float(z.min()), "max_z": float(z.max())})
    O = pd.DataFrame(out)
    O["fdr"] = bh(O["p"].to_numpy())
    O["fdr_ttest"] = bh(O["p_ttest_vs_NT"].fillna(1).to_numpy())
    O["emp_nominal"] = O["p"] < 0.05
    O["supported"] = (O["fdr"] < 0.05) & O["reps_same_sign"]
    # the criterion drawn in Figure 21: the paired t-test, as in Figures 20 and 22
    O["nominal"] = O["p_paired"] < 0.05
    O["nominal_welch"] = O["p_ttest_vs_NT"] < 0.05      # the earlier, unpaired version
    O["fdr_paired"] = bh(O["p_paired"].fillna(1).to_numpy())
    O["strong"] = O["supported"] & (O["mean_z"].abs() >= 2)
    O["qc_pass"] = line not in QC_FAIL
    rows.append(O)

L = pd.concat(rows, ignore_index=True)
L["is_risk"] = L["gene"].isin(risk)
L["twas_direction"] = [allm["twas_direction"].get(g, "") for g in L["gene"]]
L["twas_support"] = [support_call(d, dp, s) for d, dp, s
                     in zip(L["twas_direction"], L["delta_paired"], L["nominal"])]
L.to_csv("results/sirna_by_cellline.csv", index=False)
print("\n=== calls per cell line, whole 333-gene panel ===")
print(L.groupby("cell_line")[["emp_nominal", "supported", "nominal", "strong"]].sum().to_string())

# ---------------------------------------------------------------- the risk genes
R = L[L.is_risk].copy()
genes = sorted(R.gene.unique())
print(f"\nS1 risk genes on the siRNA panel: {len(genes)}")

W = R.pivot_table(index="gene", columns="cell_line",
                  values=["mean_z", "fdr", "p", "supported", "emp_nominal", "pct_NT",
                          "delta_pct", "delta_paired", "p_paired", "sd_diff",
                          "p_ttest_vs_NT", "nominal", "sd_pct_NT"])
flat = pd.DataFrame(index=W.index)
for line in LINES:
    flat[f"z_{line}"] = W[("mean_z", line)]
    flat[f"p_{line}"] = W[("p", line)]
    flat[f"fdr_{line}"] = W[("fdr", line)]
    flat[f"pctNT_{line}"] = W[("pct_NT", line)]
    flat[f"delta_{line}"] = W[("delta_paired", line)]        # paired difference, primary
    flat[f"delta_unpaired_{line}"] = W[("delta_pct", line)]
    flat[f"pPaired_{line}"] = W[("p_paired", line)]
    flat[f"sdDiff_{line}"] = W[("sd_diff", line)]
    flat[f"emp_{line}"] = W[("emp_nominal", line)].astype(bool)
    flat[f"sup_{line}"] = W[("supported", line)].astype(bool)
    flat[f"pT_{line}"] = W[("p_ttest_vs_NT", line)]
    flat[f"nom_{line}"] = W[("nominal", line)].astype(bool)
    flat[f"sdpct_{line}"] = W[("sd_pct_NT", line)]

QCOK = [l for l in LINES if l not in QC_FAIL]
for tag, col in [("supported", "sup"), ("nominal", "nom"), ("emp_nominal", "emp")]:
    flat[f"n_lines_{tag}"] = flat[[f"{col}_{l}" for l in LINES]].sum(axis=1)
    flat[f"{tag}_any"] = flat[f"n_lines_{tag}"] > 0
    # the same counts restricted to the two lines that pass QC - the defensible version
    flat[f"n_qc_{tag}"] = flat[[f"{col}_{l}" for l in QCOK]].sum(axis=1)
    flat[f"{tag}_any_qc"] = flat[f"n_qc_{tag}"] > 0

flat["mean_z_all"] = flat[[f"z_{l}" for l in LINES]].mean(axis=1)
flat["mean_z_qc"] = flat[[f"z_{l}" for l in QCOK]].mean(axis=1)
sign = np.sign(flat[[f"z_{l}" for l in LINES]])
flat["all_lines_same_sign"] = (sign.nunique(axis=1) == 1)
flat["qc_lines_same_sign"] = (np.sign(flat[[f"z_{l}" for l in QCOK]]).nunique(axis=1) == 1)


def direction(r, col):
    s = [np.sign(r[f"z_{l}"]) for l in LINES if r[f"{col}_{l}"]]
    if not s:
        return ""
    if len(set(s)) > 1:
        return "conflicting"
    return "reduces viability" if s[0] < 0 else "increases viability"


flat["supported_direction"] = flat.apply(direction, axis=1, col="sup")
flat["nominal_direction"] = flat.apply(direction, axis=1, col="nom")
flat["twas_direction"] = [allm["twas_direction"].get(g, "") for g in flat.index]
flat["class"] = [allm["class"].get(g, "not measured in GeoMx") for g in flat.index]

# ---- the t-test view: per-line TWAS support, and the gene-level verdict.
# QC-passing lines only for the verdict; SW480 is carried per line but excluded from it,
# for the same reason its FDR calls are not counted.
for line in LINES:
    flat[f"sup_twas_{line}"] = [support_call(d, dp, s) for d, dp, s in
                                zip(flat["twas_direction"], flat[f"delta_{line}"],
                                    flat[f"nom_{line}"])]
flat["twas_support"] = [gene_call([r[f"sup_twas_{l}"] for l in QCOK])
                        for _, r in flat.iterrows()]
flat["twas_support_all3"] = [gene_call([r[f"sup_twas_{l}"] for l in LINES])
                             for _, r in flat.iterrows()]
flat["n_lines_supporting_twas"] = sum(
    (flat[f"sup_twas_{l}"] == SUPPORTS).astype(int) for l in QCOK)
flat["n_lines_opposing_twas"] = sum(
    (flat[f"sup_twas_{l}"] == OPPOSES).astype(int) for l in QCOK)
flat["mean_delta_all"] = flat[[f"delta_{l}" for l in LINES]].mean(axis=1)
flat["mean_delta_qc"] = flat[[f"delta_{l}" for l in QCOK]].mean(axis=1)
flat["spatially_validated"] = [g in valset for g in flat.index]
flat = flat.reset_index().sort_values(["n_lines_supported", "mean_z_all"],
                                      ascending=[False, True])
flat.to_csv("results/sirna_risk_by_cellline.csv", index=False)

show = ["gene", "twas_direction", "class", "z_HCT116", "z_RKO", "z_SW480",
        "n_lines_emp_nominal", "n_lines_supported", "n_lines_nominal",
        "supported_direction"]
print("\n=== risk genes with any support, by cell line ===")
print(flat[flat.supported_any | flat.nominal_any][show].round(3).to_string(index=False))

n = len(flat)
print(f"\n--- criterion by criterion, over {n} risk genes ---")
for tag in ["emp_nominal", "supported", "nominal"]:
    print(f"{tag:12s} any line {int(flat[tag + '_any'].sum()):3d} | "
          f"QC-passing lines only {int(flat[tag + '_any_qc'].sum()):3d} | "
          f"by n lines "
          f"{dict(flat['n_lines_' + tag].value_counts().sort_index())}")

print("\nper line (risk genes):")
print(R.groupby("cell_line")[["emp_nominal", "supported", "nominal"]].sum().to_string())

nsup = int(flat["supported_any_qc"].sum())
# TWAS direction vs measured effect, on the two QC-passing lines
ri = flat[flat.twas_direction == "risk-increasing"]["mean_z_qc"]
pr = flat[flat.twas_direction == "protective"]["mean_z_qc"]
u = mannwhitneyu(ri, pr, alternative="less")
print(f"\nTWAS direction vs measured effect (HCT116+RKO): risk-increasing (n={len(ri)}) "
      f"median {ri.median():+.3f} vs protective (n={len(pr)}) {pr.median():+.3f}; "
      f"one-sided P = {u.pvalue:.3f}")

other = L[(~L.is_risk) & (L.cell_line.isin(QCOK))].pivot_table(
    index="gene", columns="cell_line", values="supported")
other_any = (other.fillna(False).astype(bool).sum(axis=1) > 0)
tab = [[nsup, n - nsup], [int(other_any.sum()), len(other_any) - int(other_any.sum())]]
orr, pf = fisher_exact(tab)
print(f"risk genes supported in >=1 QC-passing line: {nsup}/{n} vs "
      f"{int(other_any.sum())}/{len(other_any)} of the other genes "
      f"(OR {orr:.2f}, P = {pf:.3g})")

# concordance with the CRISPR per-line result on the shared genes
try:
    Cw = pd.read_csv("results/crispr_risk36_by_cellline.csv").set_index("gene")
    shared = [g for g in flat["gene"] if g in Cw.index]
    a = flat.set_index("gene").loc[shared, "mean_z_qc"]
    b = Cw.loc[shared, [f"z_{l}" for l in ["HCT116", "RKO"]]].mean(axis=1)
    from scipy.stats import spearmanr
    rho, prho = spearmanr(a, b)
    same = int((np.sign(a) == np.sign(b)).sum())
    print(f"\nshared genes siRNA vs CRISPR (HCT116+RKO mean z): n={len(shared)}, "
          f"rho={rho:+.3f} (P={prho:.3g}), same sign {same}/{len(shared)}")
    cross = {"n": len(shared), "rho": float(rho), "p": float(prho), "same_sign": same}
except Exception as e:                                        # pragma: no cover
    print("cross-assay comparison skipped:", e)
    cross = None

json.dump({
    "n_panel": int(L.gene.nunique()),
    "n_risk_on_panel": n,
    "qc_fail_lines": sorted(QC_FAIL),
    "z_prime_mean": {l: float(QC.loc[QC.cell_line == l, "z_prime"].mean()) for l in LINES},
    "per_line": {c: {l: int(R[(R.cell_line == l)][c].sum()) for l in LINES}
                 for c in ["emp_nominal", "supported", "nominal"]},
    "per_line_all": {c: {l: int(L[(L.cell_line == l)][c].sum()) for l in LINES}
                     for c in ["emp_nominal", "supported", "nominal"]},
    "n_any": {c: int(flat[c + "_any"].sum()) for c in ["emp_nominal", "supported", "nominal"]},
    "n_any_qc": {c: int(flat[c + "_any_qc"].sum())
                 for c in ["emp_nominal", "supported", "nominal"]},
    "n_by_lines": {c: {int(k): int(v) for k, v in
                       flat["n_lines_" + c].value_counts().sort_index().items()}
                   for c in ["emp_nominal", "supported", "nominal"]},
    "nesting_check": {
        "emp_nominal_contains_fdr": bool(
            all(L.loc[L.supported, "emp_nominal"])),
        "ttest_and_fdr_overlap": {
            "both": int((flat["nominal_any"] & flat["supported_any"]).sum()),
            "ttest_only": int((flat["nominal_any"] & ~flat["supported_any"]).sum()),
            "fdr_only": int((~flat["nominal_any"] & flat["supported_any"]).sum()),
            "either": int((flat["nominal_any"] | flat["supported_any"]).sum())},
    },
    "twas_direction_counts": {k: int(v) for k, v in flat["twas_direction"].value_counts().items()},
    # ---- the t-test view used by Figure 21
    "ttest": {
        "per_line_sig": {l: int(R[R.cell_line == l]["nominal"].sum()) for l in LINES},
        "per_line_sig_all": {l: int(L[L.cell_line == l]["nominal"].sum()) for l in LINES},
        "per_line_support": {l: int((R[R.cell_line == l]["twas_support"] == SUPPORTS).sum())
                             for l in LINES},
        "per_line_oppose": {l: int((R[R.cell_line == l]["twas_support"] == OPPOSES).sum())
                            for l in LINES},
        "gene_verdict": {k: int(v) for k, v in flat["twas_support"].value_counts().items()},
        "gene_verdict_all3": {k: int(v) for k, v in
                              flat["twas_support_all3"].value_counts().items()},
        "n_by_lines_sig": {int(k): int(v) for k, v in
                           flat["n_lines_nominal"].value_counts().sort_index().items()},
        # keyed on the QC-passing lines, to match the verdict, which ignores SW480
        "by_lines_verdict": {
            str(k): {v: int(((flat["n_qc_nominal"] == k) &
                             (flat["twas_support"] == v)).sum())
                     for v in [SUPPORTS, OPPOSES, "conflicting", NS]}
            for k in sorted(flat["n_qc_nominal"].unique())},
        "concordance_calls_qc": binomial_concordance(
            int((R[R.cell_line.isin(QCOK)]["twas_support"] == SUPPORTS).sum()),
            int((R[R.cell_line.isin(QCOK)]["twas_support"] == OPPOSES).sum())),
        "concordance_calls_all3": binomial_concordance(
            int((R["twas_support"] == SUPPORTS).sum()),
            int((R["twas_support"] == OPPOSES).sum())),
        "concordance_calls_sw480": binomial_concordance(
            int((R[R.cell_line == "SW480"]["twas_support"] == SUPPORTS).sum()),
            int((R[R.cell_line == "SW480"]["twas_support"] == OPPOSES).sum())),
        "concordance_genes": binomial_concordance(
            int((flat["twas_support"] == SUPPORTS).sum()),
            int((flat["twas_support"] == OPPOSES).sum())),
        "supporting_genes": flat.loc[flat["twas_support"] == SUPPORTS, "gene"].tolist(),
        "opposing_genes": flat.loc[flat["twas_support"] == OPPOSES, "gene"].tolist(),
        "conflicting_genes": flat.loc[flat["twas_support"] == "conflicting", "gene"].tolist(),
    },
    "twas_direction_vs_effect": {
        "n_risk_increasing": int(len(ri)), "n_protective": int(len(pr)),
        "median_z_risk_increasing": float(ri.median()),
        "median_z_protective": float(pr.median()), "p_one_sided": float(u.pvalue)},
    "direction_counts": {k: int(v) for k, v in
                         flat[flat.supported_any]["supported_direction"].value_counts().items()},
    "nominal_direction_counts": {k: int(v) for k, v in
                                 flat[flat.nominal_any]["nominal_direction"].value_counts().items()},
    "n_all_same_sign": int(flat["all_lines_same_sign"].sum()),
    "n_qc_same_sign": int(flat["qc_lines_same_sign"].sum()),
    "vs_other": {"risk": nsup, "n_risk": n, "other": int(other_any.sum()),
                 "n_other": int(len(other_any)), "odds_ratio": float(orr), "p": float(pf)},
    "cross_assay": cross,
    "table": flat[show + ["fdr_HCT116", "fdr_RKO", "fdr_SW480", "supported_any",
                          "supported_any_qc", "nominal_any", "spatially_validated"]
                  ].to_dict("records"),
}, open("results/sirna_by_cellline.json", "w"), indent=1, default=float)
print("\nDone.")
