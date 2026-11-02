#!/usr/bin/env python3
"""
35 - The CRISPR-Cas9 screen, cell line by cell line, for the risk genes it covers.

Section 2.15 pooled HCT116 and RKO because SW480 fails QC in the *siRNA* screen. In the
CRISPR screen SW480 passes (mean Z' 0.59 against 0.50 for HCT116 and 0.78 for RKO), so
all three lines are analysed here and the pooled two-line result of §2.15 is left as it
stands. Keeping the two analyses separate is deliberate: the pooled test asks whether a
gene matters on a common background, this one asks in which lines it matters at all.

For each cell line separately:
  * each of its three replicates is converted to a robust z against that plate's own
    sample median and MAD, exactly as in script 33;
  * a gene's effect is the mean z over its three wells;
  * the p-value is the two-sided position of that mean in a null of 20,000 draws of one
    non-targeting well from each of the line's three plates, so every line is judged
    against its own noise;
  * FDR is computed across the 116 genes within that line.

A gene is "supported" in a line when it clears FDR < 0.05 there and its three replicates
agree in sign. Overall support is support in at least one line, and the number of lines
is reported alongside, because one line out of three is a much weaker claim than three.

Outputs
  results/crispr_by_cellline.csv        all 116 genes x 3 lines
  results/crispr_risk36_by_cellline.csv the 36 S1 risk genes, wide
  results/crispr_by_cellline.json       counts for the report
"""
import pandas as pd, numpy as np, json, sys, os
from scipy.stats import fisher_exact, ttest_ind, ttest_1samp, ttest_rel

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _ttest_support import (support_call, gene_call, binomial_concordance,
                            SUPPORTS, OPPOSES, NS)

ALIAS = {"ACAS2": "ACSS2"}
LINES = ["HCT116", "RKO", "SW480"]
N_NULL = 20000


def bh(p):
    p = np.asarray(p, float)
    o = np.argsort(p)
    q = np.empty_like(p)
    q[o] = np.minimum(1, np.minimum.accumulate(
        (p[o] * len(p) / (np.arange(len(p)) + 1))[::-1])[::-1])
    return q


C = pd.read_csv("results/func_crispr_long.csv")
C["gene"] = C["gene"].replace(ALIAS)

# ---- replicate-level QC (script 42). A plate is dropped when it is the worst of its
# cell line and beyond the screen-wide quartile on three or more of six metrics, or when
# its Z' is at or below zero. The exclusion list is read from the QC output rather than
# written here, so the rule and its application cannot drift apart.
#
# APPLY_REPLICATE_QC is FALSE by decision: dropping one plate of three leaves the paired
# t-test with a single degree of freedom, where the critical value of t is 12.71 instead
# of 4.30, and the p-values collapse while the effects barely move. All nine plates are
# analysed. Script 42 still runs and Table S51 still reports the metrics, and Figure 23
# shows the concordance directly, so the reader can see what was kept and judge it.
APPLY_REPLICATE_QC = False
EXCLUDE = set()
if APPLY_REPLICATE_QC:
    try:
        _qc = json.load(open("results/crispr_replicate_qc.json"))
        EXCLUDE = {tuple(x.split()) for x in _qc["excluded"]}
    except FileNotFoundError:
        print("no replicate QC found - run 42_crispr_replicate_qc.py first")
else:
    print("replicate QC not applied - all plates retained (see script 42 and Table S51)")
if EXCLUDE:
    before = len(C)
    C = C[~pd.Series(list(zip(C.cell_line, C.replicate)), index=C.index).isin(EXCLUDE)]
    print("QC: dropped " + ", ".join(f"{l} {r}" for l, r in sorted(EXCLUDE))
          + f"  ({before - len(C)} wells removed)")
risk = set(pd.read_csv("results/risk_genes_harmonized.csv")["gene"])
allm = pd.read_csv("results/validated_risk_genes_all.csv").set_index("gene")
valset = set(pd.read_csv("results/validated_priority.csv")["gene"])

rows = []
rng = np.random.default_rng(35)
for line in LINES:
    d = C[C.cell_line == line].copy()
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
    print(f"{line}: {len(nt_by_plate)} plates, null sd {null.std():.3f}")

    # the control distribution for the t-test: every non-targeting well of this line,
    # expressed as % of its own plate's NT median so plate level is removed
    nt_pct = Z.loc[Z.kind == "negative", "pct_NT"].to_numpy(float)

    # --- the paired design. Each replicate plate carries 8 non-targeting wells; their
    # mean is that plate's control value, so a gene has one (gene, control) pair per
    # plate and the three pairs are tested with a paired t-test. Pairing inside the
    # plate removes the plate-to-plate offset without any normalisation step, and the
    # difference is already in viability points because the source values are on a
    # 100 = untreated scale.
    nt_mean_by_rep = (Z[Z.kind == "negative"].groupby("replicate")["value"].mean())
    reps = sorted(nt_mean_by_rep.index)
    print(f"  plate control means: "
          + ", ".join(f"{r} {nt_mean_by_rep[r]:.1f}" for r in reps))

    out = []
    for g, dd in Z[Z.kind == "sample"].groupby("gene"):
        z = dd["z"].to_numpy(float)
        pct = dd["pct_NT"].to_numpy(float)
        stat = float(np.mean(z))
        p = (int(np.sum(np.abs(null) >= abs(stat))) + 1) / (len(null) + 1)
        # requested: a t-test of the three replicates against the controls. Welch's is
        # used because the three replicates and the 24 non-targeting wells have no
        # reason to share a variance. The one-sample form against the control mean of
        # 100 is kept alongside; with two degrees of freedom it is the weaker of the two.
        tt = ttest_ind(pct, nt_pct, equal_var=False) if len(pct) > 1 else None
        t1 = ttest_1samp(pct, 100.0) if len(pct) > 1 else None

        # paired: the gene's value on each plate against that plate's own control mean
        s = dd.set_index("replicate")["value"]
        gv = np.array([s[r] for r in reps], float)
        cv = np.array([nt_mean_by_rep[r] for r in reps], float)
        diff = gv - cv
        tp = ttest_rel(gv, cv) if len(gv) > 1 else None

        out.append({"cell_line": line, "gene": g, "n_rep": len(z), "mean_z": stat,
                    "pct_NT": float(np.mean(pct)),
                    "nt_mean_pct": float(np.mean(nt_pct)),
                    "delta_pct": float(np.mean(pct) - np.mean(nt_pct)),
                    "sd_pct_NT": float(np.std(pct, ddof=1)) if len(pct) > 1 else np.nan,
                    # --- the paired analysis, primary for Figures 20 and 22
                    "value_mean": float(gv.mean()),
                    "control_mean": float(cv.mean()),
                    "delta_paired": float(diff.mean()),
                    "sd_diff": float(diff.std(ddof=1)) if len(diff) > 1 else np.nan,
                    "t_paired": float(tp.statistic) if tp else np.nan,
                    "p_paired": float(tp.pvalue) if tp else np.nan,
                    "diffs": ";".join(f"{d:+.2f}" for d in diff),
                    "pct_of_control": float(100 * gv.mean() / cv.mean()),
                    "p": p,
                    "p_ttest_vs_NT": float(tt.pvalue) if tt else np.nan,
                    "p_ttest_vs_100": float(t1.pvalue) if t1 else np.nan,
                    "reps_same_sign": bool(np.all(np.sign(z) == np.sign(z[0]))),
                    "min_z": float(z.min()), "max_z": float(z.max())})
    O = pd.DataFrame(out)
    O["fdr"] = bh(O["p"].to_numpy())
    O["fdr_ttest"] = bh(O["p_ttest_vs_NT"].fillna(1).to_numpy())
    # the uncorrected form of the same empirical test, carried so that the figure can
    # show what correction actually does: it can only remove genes from this set, never
    # add them. Without it a reader sees the FDR count exceed the t-test count and
    # reasonably concludes that correction made the criterion more lenient.
    O["emp_nominal"] = O["p"] < 0.05
    O["supported"] = (O["fdr"] < 0.05) & O["reps_same_sign"]
    # the criterion drawn in Figures 20 and 22: the paired t-test
    O["nominal"] = (O["p_paired"] < 0.05)
    O["nominal_welch"] = (O["p_ttest_vs_NT"] < 0.05)    # the earlier, unpaired version
    O["fdr_paired"] = bh(O["p_paired"].fillna(1).to_numpy())
    O["strong"] = O["supported"] & (O["mean_z"].abs() >= 2)
    rows.append(O)

L = pd.concat(rows, ignore_index=True)
L["is_risk"] = L["gene"].isin(risk)
L["twas_direction"] = [allm["twas_direction"].get(g, "") for g in L["gene"]]
L["twas_support"] = [support_call(d, dp, s) for d, dp, s
                     in zip(L["twas_direction"], L["delta_paired"], L["nominal"])]
L.to_csv("results/crispr_by_cellline.csv", index=False)
print("\n=== supported genes per cell line (of 116) ===")
print(L.groupby("cell_line")[["emp_nominal", "supported", "nominal", "strong"]].sum().to_string())

# ---------------------------------------------------------------- the 36 risk genes
R = L[L.is_risk].copy()
genes = sorted(R.gene.unique())
print(f"\nS1 risk genes on the CRISPR panel: {len(genes)}")

W = R.pivot_table(index="gene", columns="cell_line",
                  values=["mean_z", "fdr", "p", "supported", "emp_nominal", "pct_NT",
                          "delta_pct", "delta_paired", "p_paired", "sd_diff",
                          "p_ttest_vs_NT", "nominal", "sd_pct_NT"])
flat = pd.DataFrame(index=W.index)
for line in LINES:
    flat[f"z_{line}"] = W[("mean_z", line)]
    flat[f"p_{line}"] = W[("p", line)]
    flat[f"emp_{line}"] = W[("emp_nominal", line)].astype(bool)
    flat[f"fdr_{line}"] = W[("fdr", line)]
    flat[f"pctNT_{line}"] = W[("pct_NT", line)]
    flat[f"delta_{line}"] = W[("delta_paired", line)]      # paired difference, primary
    flat[f"delta_unpaired_{line}"] = W[("delta_pct", line)]
    flat[f"pPaired_{line}"] = W[("p_paired", line)]
    flat[f"sdDiff_{line}"] = W[("sd_diff", line)]
    flat[f"sup_{line}"] = W[("supported", line)].astype(bool)
    flat[f"pT_{line}"] = W[("p_ttest_vs_NT", line)]
    flat[f"nom_{line}"] = W[("nominal", line)].astype(bool)
    flat[f"sdpct_{line}"] = W[("sd_pct_NT", line)]
flat["n_lines_supported"] = flat[[f"sup_{l}" for l in LINES]].sum(axis=1)
flat["supported_any"] = flat["n_lines_supported"] > 0
flat["n_lines_nominal"] = flat[[f"nom_{l}" for l in LINES]].sum(axis=1)
flat["nominal_any"] = flat["n_lines_nominal"] > 0
flat["n_lines_emp_nominal"] = flat[[f"emp_{l}" for l in LINES]].sum(axis=1)
flat["emp_nominal_any"] = flat["n_lines_emp_nominal"] > 0
flat["mean_z_all"] = flat[[f"z_{l}" for l in LINES]].mean(axis=1)
sign = np.sign(flat[[f"z_{l}" for l in LINES]])
flat["all_lines_same_sign"] = (sign.nunique(axis=1) == 1)
# direction agreement restricted to the lines that actually support the gene
def nom_dir(r):
    s = [np.sign(r[f"z_{l}"]) for l in LINES if r[f"nom_{l}"]]
    if not s:
        return ""
    if len(set(s)) > 1:
        return "conflicting"
    return "reduces viability" if s[0] < 0 else "increases viability"


def sup_dir(r):
    s = [np.sign(r[f"z_{l}"]) for l in LINES if r[f"sup_{l}"]]
    if not s:
        return ""
    if len(set(s)) > 1:
        return "conflicting"
    return "reduces viability" if s[0] < 0 else "increases viability"


flat["supported_direction"] = flat.apply(sup_dir, axis=1)
flat["nominal_direction"] = flat.apply(nom_dir, axis=1)
flat["twas_direction"] = [allm["twas_direction"].get(g, "") for g in flat.index]
flat["class"] = [allm["class"].get(g, "not measured in GeoMx") for g in flat.index]

# ---- the t-test view: per-line TWAS support, and the gene-level verdict
for line in LINES:
    flat[f"sup_twas_{line}"] = [support_call(d, dp, s) for d, dp, s in
                                zip(flat["twas_direction"], flat[f"delta_{line}"],
                                    flat[f"nom_{line}"])]
flat["twas_support"] = [gene_call([r[f"sup_twas_{l}"] for l in LINES])
                        for _, r in flat.iterrows()]
flat["n_lines_supporting_twas"] = sum(
    (flat[f"sup_twas_{l}"] == SUPPORTS).astype(int) for l in LINES)
flat["n_lines_opposing_twas"] = sum(
    (flat[f"sup_twas_{l}"] == OPPOSES).astype(int) for l in LINES)
flat["mean_delta_all"] = flat[[f"delta_{l}" for l in LINES]].mean(axis=1)
flat["spatially_validated"] = [g in valset for g in flat.index]
flat = flat.reset_index().sort_values(["n_lines_supported", "mean_z_all"],
                                      ascending=[False, True])
flat.to_csv("results/crispr_risk36_by_cellline.csv", index=False)

print("\n=== the 36 risk genes, by cell line ===")
show = ["gene", "twas_direction", "class", "z_HCT116", "z_RKO", "z_SW480",
        "n_lines_nominal", "n_lines_supported", "nominal_direction"]
print(flat[show].round(3).to_string(index=False))

nsup = int(flat["supported_any"].sum())
byn = flat["n_lines_supported"].value_counts().sort_index()
nnom = int(flat["nominal_any"].sum())
print(f"\nnominal t-test P<0.05 in at least one line: {nnom}/{len(flat)}")
print("number of lines nominally significant:")
print(flat["n_lines_nominal"].value_counts().sort_index().to_string())
print("direction among nominally significant genes:")
print(flat[flat.nominal_any]["nominal_direction"].value_counts().to_string())
print(f"\nsupported (FDR) in at least one line: {nsup}/{len(flat)}")
print("number of lines supporting:"); print(byn.to_string())
print(f"all three lines same sign: {int(flat['all_lines_same_sign'].sum())}/{len(flat)}")
print("direction among supported genes:")
print(flat[flat.supported_any]["supported_direction"].value_counts().to_string())

# now that the TWAS direction is on the figure, test whether it predicts anything here:
# a risk-increasing gene should, if the germline effect carried into the tumour cell, be
# more needed for growth than a protective one.
from scipy.stats import mannwhitneyu
ri = flat[flat.twas_direction == "risk-increasing"]["mean_z_all"]
pr = flat[flat.twas_direction == "protective"]["mean_z_all"]
u = mannwhitneyu(ri, pr, alternative="less")
print(f"\nTWAS direction vs measured effect: risk-increasing (n={len(ri)}) median "
      f"{ri.median():+.3f} vs protective (n={len(pr)}) {pr.median():+.3f}; "
      f"one-sided P = {u.pvalue:.3f}")

# is the support rate different for risk genes than the rest of the panel?
other = L[~L.is_risk].pivot_table(index="gene", columns="cell_line", values="supported")
other_any = (other.fillna(False).astype(bool).sum(axis=1) > 0)
tab = [[nsup, len(flat) - nsup], [int(other_any.sum()), len(other_any) - int(other_any.sum())]]
orr, pf = fisher_exact(tab)
print(f"\nrisk genes supported in >=1 line: {nsup}/{len(flat)} vs "
      f"{int(other_any.sum())}/{len(other_any)} of the other genes "
      f"(OR {orr:.2f}, P = {pf:.3g})")

json.dump({
    "n_risk_on_panel": int(len(flat)),
    "per_line_supported": {l: int(L[(L.cell_line == l) & L.is_risk]["supported"].sum())
                           for l in LINES},
    "per_line_supported_all": {l: int(L[L.cell_line == l]["supported"].sum()) for l in LINES},
    "n_supported_any": nsup,
    "n_nominal_any": int(flat["nominal_any"].sum()),
    "per_line_nominal": {l: int(L[(L.cell_line == l) & L.is_risk]["nominal"].sum())
                         for l in LINES},
    "per_line_nominal_all": {l: int(L[L.cell_line == l]["nominal"].sum()) for l in LINES},
    "n_by_lines_nominal": {int(k): int(v) for k, v in
                           flat["n_lines_nominal"].value_counts().sort_index().items()},
    "n_by_lines_emp_nominal": {int(k): int(v) for k, v in
                               flat["n_lines_emp_nominal"].value_counts().sort_index().items()},
    "n_emp_nominal_any": int(flat["emp_nominal_any"].sum()),
    "per_line_emp_nominal": {l: int(L[(L.cell_line == l) & L.is_risk]["emp_nominal"].sum())
                             for l in LINES},
    "per_line_emp_nominal_all": {l: int(L[L.cell_line == l]["emp_nominal"].sum())
                                 for l in LINES},
    "nesting_check": {"every_fdr_gene_is_also_emp_nominal":
                      bool(all(L.loc[L.supported, "emp_nominal"]))},
    # ---- the t-test view used by Figure 20
    "ttest": {
        "per_line_sig": {l: int(R[R.cell_line == l]["nominal"].sum()) for l in LINES},
        "per_line_sig_all": {l: int(L[L.cell_line == l]["nominal"].sum()) for l in LINES},
        "per_line_support": {l: int((R[R.cell_line == l]["twas_support"] == SUPPORTS).sum())
                             for l in LINES},
        "per_line_oppose": {l: int((R[R.cell_line == l]["twas_support"] == OPPOSES).sum())
                            for l in LINES},
        "gene_verdict": {k: int(v) for k, v in flat["twas_support"].value_counts().items()},
        "n_by_lines_sig": {int(k): int(v) for k, v in
                           flat["n_lines_nominal"].value_counts().sort_index().items()},
        "by_lines_verdict": {
            str(k): {v: int(((flat["n_lines_nominal"] == k) &
                             (flat["twas_support"] == v)).sum())
                     for v in [SUPPORTS, OPPOSES, "conflicting", NS]}
            for k in sorted(flat["n_lines_nominal"].unique())},
        "concordance_calls": binomial_concordance(
            int((R["twas_support"] == SUPPORTS).sum()),
            int((R["twas_support"] == OPPOSES).sum())),
        "concordance_genes": binomial_concordance(
            int((flat["twas_support"] == SUPPORTS).sum()),
            int((flat["twas_support"] == OPPOSES).sum())),
        "supporting_genes": flat.loc[flat["twas_support"] == SUPPORTS,
                                     "gene"].tolist(),
        "opposing_genes": flat.loc[flat["twas_support"] == OPPOSES, "gene"].tolist(),
        "opposing_delta": {r["gene"]: float(np.nanmean([r[f"delta_{l}"] for l in LINES]))
                           for _, r in flat[flat["twas_support"] == OPPOSES].iterrows()},
        "conflicting_genes": flat.loc[flat["twas_support"] == "conflicting",
                                      "gene"].tolist(),
    },
    "twas_direction_vs_effect": {
        "n_risk_increasing": int(len(ri)), "n_protective": int(len(pr)),
        "median_z_risk_increasing": float(ri.median()),
        "median_z_protective": float(pr.median()),
        "p_one_sided": float(u.pvalue)},
    "criteria_overlap": {
        "both": int((flat["nominal_any"] & flat["supported_any"]).sum()),
        "nominal_only": int((flat["nominal_any"] & ~flat["supported_any"]).sum()),
        "fdr_only": int((~flat["nominal_any"] & flat["supported_any"]).sum()),
        "either": int((flat["nominal_any"] | flat["supported_any"]).sum()),
    },
    "nominal_direction_counts": {k: int(v) for k, v in
                                 flat[flat.nominal_any]["nominal_direction"]
                                 .value_counts().items()},
    "twas_direction_counts": {k: int(v) for k, v in
                              flat["twas_direction"].value_counts().items()},
    "n_by_lines": {int(k): int(v) for k, v in byn.items()},
    "n_all_same_sign": int(flat["all_lines_same_sign"].sum()),
    "direction_counts": {k: int(v) for k, v in
                         flat[flat.supported_any]["supported_direction"].value_counts().items()},
    "vs_other": {"risk": nsup, "n_risk": int(len(flat)),
                 "other": int(other_any.sum()), "n_other": int(len(other_any)),
                 "odds_ratio": float(orr), "p": float(pf)},
    "table": flat[show + ["fdr_HCT116", "fdr_RKO", "fdr_SW480", "supported_any",
                          "spatially_validated"]].to_dict("records"),
}, open("results/crispr_by_cellline.json", "w"), indent=1, default=float)
print("\nDone.")
