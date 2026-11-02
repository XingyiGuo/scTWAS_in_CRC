#!/usr/bin/env python
"""Do the CRC risk-gene drugs themselves reverse the TWAS profile?

The library-wide scan (scripts 04-06) is a screen over 16,767 signatures and is
null. This script asks the pre-specified version of the question, which is far
better powered: of the 193 approved or phase II/III compounds that target the 63
druggable scTWAS risk genes, which are present in iLINCS, and do they reverse the
cell-type TWAS profiles?

Only compounds with prior evidence enter, so the multiple-testing burden is the
size of the panel rather than the library, and a modest negative connectivity
score carries more weight here than a larger one selected out of 16,767.

Tests
-----
1. Panel versus library: Mann-Whitney U on NCS, panel signatures against all
   others, pooled over models. Asks whether risk-gene drugs as a class shift
   negative. No permutation needed.
2. Per compound: mean NCS over its signatures and models, the fraction of its
   signatures with NCS < 0, and a one-sample Wilcoxon against 0, BH-corrected
   across the compounds in the panel.
3. Two-way evidence: compounds that both reverse the profile and target one of
   the 63 genes are flagged separately - those are the nominations.

Outputs (under <ROOT>/iLINCS_reversal):
  table/ilincs_drug_panel_matching.csv      every one of the 193 drugs, matched or not
  table/ilincs_drug_panel_per_compound.csv  per-compound connectivity and test
  table/ilincs_drug_panel_per_model.csv     compound x model mean NCS
  table/ilincs_drug_panel_summary.csv       panel-vs-library test
"""
import json
import os
import re
import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu
from statsmodels.stats.multitest import multipletests

ROOT = r"W:\VUMC_projects\scTWAS"
OUT = rf"{ROOT}\iLINCS_reversal"
DRUGS = rf"{ROOT}\Drugs\table\63_druggable_genes_drug_records.csv"
CMAP = rf"{OUT}\table\ilincs_cmap_all_signatures.csv"
CIDS = "handoff_panel_cids.json"
META = r"L:\Drugs_Proteins\iLINCS\meta_sigSearch_LIB_5_143374_20231227.csv"

SALTS = (r"\s+(hydrochloride|hcl|mesylate|maleate|citrate|tartrate|sulfate|sulphate|"
         r"phosphate|acetate|fumarate|succinate|besylate|tosylate|sodium|potassium|"
         r"calcium|dihydrate|monohydrate|hydrate|ethanolate|malate|oxalate|bromide|"
         r"chloride|nitrate|lactate|gluconate|carbonate|benzoate|palmitate|"
         r"dimaleate|dimesylate|disodium)\b")


def norm_name(s):
    s = str(s).strip().lower()
    s = re.sub(r"\(.*?\)", " ", s)
    for _ in range(3):
        s = re.sub(SALTS, " ", s)
    return re.sub(r"[^a-z0-9]", "", s)


def main():
    rec = pd.read_csv(DRUGS)
    adv = rec[rec["Development stage (harmonised)"].isin(
        ["Approved drug", "Phase 3 / registration", "Phase 2"])].copy()
    adv["cpd"] = adv["Drug Brand Name"].astype(str).str.strip()
    panel = (adv.groupby(adv.cpd.map(norm_name))
                .agg(drug_name=("cpd", lambda s: sorted(s.unique())[0]),
                     risk_gene_targets=("Genes", lambda s: ";".join(sorted(set(s)))),
                     n_risk_genes=("Genes", "nunique"),
                     stages=("Development stage (harmonised)",
                             lambda s: ";".join(sorted(set(s)))))
                .reset_index().rename(columns={"cpd": "key"}))
    panel = panel[panel.key.ne("")]
    print("approved / phase II-III compounds in the panel:", len(panel))

    cm = pd.read_csv(CMAP)
    # Identity is the PubChem CID where both sides have one, normalised name otherwise.
    # String matching alone silently dropped bortezomib, ixazomib and the ERBB2 code
    # names, and double-counted aliases (HKI-272 = neratinib, BIBW 2992 = afatinib).
    meta = pd.read_csv(META)[["SignatureId", "pubChemID"]]
    cm = cm.merge(meta, left_on="signature", right_on="SignatureId", how="left").drop(columns="SignatureId")
    cm["key"] = np.where(cm.pubChemID.notna(), "CID" + cm.pubChemID.fillna(0).astype("int64").astype(str),
                         cm.Perturbagen.map(norm_name))
    cids = json.load(open(CIDS))
    panel["cid"] = panel.drug_name.map(cids)
    lib_cids = set(cm.loc[cm.pubChemID.notna(), "pubChemID"].astype("int64"))
    panel["key"] = np.where(panel.cid.notna() & panel.cid.fillna(-1).astype("int64").isin(lib_cids),
                            "CID" + panel.cid.fillna(0).astype("int64").astype(str), panel.key)
    # collapse alias rows that resolved to the same compound
    panel = (panel.groupby("key", as_index=False)
                  .agg(drug_name=("drug_name", lambda s: sorted(s, key=lambda x: (any(c.isdigit() for c in x), len(x)))[0]),
                       aliases=("drug_name", lambda s: "; ".join(sorted(s))),
                       risk_gene_targets=("risk_gene_targets",
                                          lambda s: ";".join(sorted({g for x in s for g in str(x).split(";")}))),
                       n_risk_genes=("n_risk_genes", "max"),
                       stages=("stages", lambda s: ";".join(sorted({y for x in s for y in str(x).split(";")}))),
                       cid=("cid", "first")))
    print("panel compounds after collapsing aliases on CID:", len(panel))
    lib_keys = set(cm.key)
    panel["in_iLINCS"] = panel.key.isin(lib_keys)
    panel["n_signatures"] = panel.key.map(cm.drop_duplicates("signature")
                                            .key.value_counts()).fillna(0).astype(int)
    print("matched to an iLINCS perturbagen:", int(panel.in_iLINCS.sum()))

    hit = cm[cm.key.isin(set(panel.loc[panel.in_iLINCS, "key"]))].copy()
    nz = cm[cm.NCS != 0]
    hit_nz = hit[hit.NCS != 0]
    u, p_u = mannwhitneyu(hit_nz.NCS, nz.loc[~nz.index.isin(hit_nz.index), "NCS"],
                          alternative="less")
    summary = pd.DataFrame([{
        "n_panel_compounds": int(panel.in_iLINCS.sum()),
        "n_panel_signatures": hit.signature.nunique(),
        "n_panel_rows": len(hit), "n_library_rows": len(cm),
        "panel_mean_NCS": hit_nz.NCS.mean(), "library_mean_NCS": nz.NCS.mean(),
        "panel_median_NCS": hit_nz.NCS.median(), "library_median_NCS": nz.NCS.median(),
        "mannwhitney_U": u, "P_panel_shifted_negative": p_u}])

    # Compound-level inference. A signature x model row is NOT an independent
    # observation (the 18 models share query genes; a compound contributes several
    # doses and cell lines), so each compound is reduced to ONE number - the mean
    # NCS over its signatures, each signature first averaged across the 18 models -
    # and referenced against a size-matched null of random signature sets drawn
    # from the whole library.
    sig_mean = cm.groupby("signature").NCS.mean()              # one value per signature
    sig_key = cm.drop_duplicates("signature").set_index("signature").key
    vals = sig_mean.to_numpy()
    rng = np.random.default_rng(7)

    def per_cpd(g):
        sigs = g.signature.unique()
        v = sig_mean.loc[sigs].to_numpy()
        return pd.Series({
            "n_signatures": len(sigs), "n_models": g.model.nunique(),
            "mean_NCS": v.mean(), "median_NCS": float(np.median(v)),
            "min_NCS_signature": v.min(), "frac_signatures_negative": float((v < 0).mean()),
            "best_model": g.loc[g.NCS.idxmin(), "model"],
            "best_signature": g.loc[g.NCS.idxmin(), "signature"],
            "cell_lines": ";".join(sorted(g.CellLine.dropna().unique()))[:60]})

    pc = hit.groupby("key").apply(per_cpd, include_groups=False).reset_index()
    pc["n_signatures"] = pc.n_signatures.astype(int)

    # size-matched null: mean NCS of k random signatures, 5,000 draws per distinct k
    null_p, null_mu, null_sd = [], [], []
    cache = {}
    for k, obs in zip(pc.n_signatures, pc.mean_NCS):
        if k not in cache:
            draws = vals[rng.integers(0, len(vals), size=(5000, k))].mean(axis=1)
            cache[k] = draws
        d = cache[k]
        null_p.append((int((d <= obs).sum()) + 1) / (len(d) + 1))
        null_mu.append(d.mean()); null_sd.append(d.std())
    pc["null_mean"], pc["null_sd"], pc["P_vs_sizematched_null"] = null_mu, null_sd, null_p
    pc["FDR_vs_sizematched_null"] = multipletests(pc.P_vs_sizematched_null, method="fdr_bh")[1]
    pc = pc.merge(panel[["key", "drug_name", "aliases", "risk_gene_targets", "n_risk_genes", "stages", "cid"]],
                  on="key", how="left").sort_values("mean_NCS")

    # the same statistic for every library compound, so the panel can be placed in it
    libc = (pd.DataFrame({"key": sig_key.values, "v": sig_mean.loc[sig_key.index].values})
            .groupby("key").v.agg(["mean", "size"]).reset_index()
            .rename(columns={"mean": "mean_NCS", "size": "n_signatures"}))
    libc["in_panel"] = libc.key.isin(set(pc.key))
    libc.to_csv(rf"{OUT}\table\ilincs_all_compounds_mean_NCS.csv", index=False)
    summary["panel_compound_mean_NCS"] = pc.mean_NCS.mean()
    summary["library_compound_mean_NCS"] = libc.loc[~libc.in_panel, "mean_NCS"].mean()
    summary["P_panel_compounds_vs_library_compounds"] = mannwhitneyu(
        pc.mean_NCS, libc.loc[~libc.in_panel, "mean_NCS"], alternative="less").pvalue

    pm = (hit.groupby(["key", "model"]).NCS.mean().reset_index()
             .merge(panel[["key", "drug_name"]], on="key", how="left")
             .pivot(index="drug_name", columns="model", values="NCS"))

    panel.to_csv(rf"{OUT}\table\ilincs_drug_panel_matching.csv", index=False)
    pc.to_csv(rf"{OUT}\table\ilincs_drug_panel_per_compound.csv", index=False)
    pm.to_csv(rf"{OUT}\table\ilincs_drug_panel_per_model.csv")
    summary.to_csv(rf"{OUT}\table\ilincs_drug_panel_summary.csv", index=False)

    print("\n", summary.T.to_string())
    print("\nmost reversing panel compounds:")
    print(pc.head(15)[["drug_name", "risk_gene_targets", "n_signatures", "mean_NCS",
                       "null_mean", "frac_signatures_negative", "P_vs_sizematched_null",
                       "FDR_vs_sizematched_null"]].to_string(index=False))
    print("\nsignificant after BH:", int((pc.FDR_vs_sizematched_null < 0.05).sum()))
    print("\nunmatched panel drugs (first 25):",
          sorted(panel.loc[~panel.in_iLINCS, "drug_name"])[:25])
    return panel, pc, pm, summary


if __name__ == "__main__":
    panel, pc, pm, summary = main()
