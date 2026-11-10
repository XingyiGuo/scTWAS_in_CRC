#!/usr/bin/env python3
"""
Do the CRC risk genes change expression along normal -> adenoma -> carcinoma in
     spatial transcriptomics, and do they move in the direction the TWAS predicts?

Data: GeoMx from Buissant des Amorie et al., Nature 2026 (Zenodo 10.5281/zenodo.17671259).
      exp2 = Whole Transcriptome Atlas, 18,441 genes, 9 T1 patients, 281 segments.
      exp1 = Cancer Transcriptome Atlas, 1,781 genes, 10 independent T1 patients.
      Segments are separated into epithelium (PanCK+) and stroma (PanCK-), and into
      four regions of the same specimen: normal mucosa (NOR), adenoma (ADE), tumour
      core (CORE) and invasive front (INV).

Differential expression is the authors' own linear mixed model,
log2(norm_exp) ~ tissueRegion + (1 + tissueRegion | Patient_ID), which keeps the
comparison within patient. Contrast signs have been harmonised upstream so a positive
log2FC always means higher in the later stage.

Tests
  A  per-gene stage contrasts for every risk gene
  B  are risk genes differentially expressed more often than the rest of the genome?
  C  DIRECTION CONCORDANCE - do risk-increasing genes go up and protective genes go
     down along progression? Against a direction-label permutation null and against
     background genes.
  D  replication in the independent exp1 cohort
  E  compartment specificity - the TWAS genes were called in epithelial cell types, so
     the effect should be stronger in epithelium than in stroma

Also exports the legacy three-transition validation set for follow-up
analyses and the two-cohort normal-to-tumour-core flag for Table 21.
Use --validation-only to refresh these outputs from existing GeoMx results.
"""
import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import binomtest, fisher_exact, mannwhitneyu, spearmanr, wilcoxon

RISK_GENES = Path("results/risk_genes_harmonized.csv")
GEOMX_RESULTS = Path("results/spatial_risk_gene_DE.csv")
MAIN_TRANSITIONS = ["NOR->ADE", "ADE->CORE", "NOR->CORE"]
INVASION_TRANSITIONS = ["ADE->INV", "NOR->INV", "CORE->INV"]
ALL_TRANSITIONS = MAIN_TRANSITIONS + INVASION_TRANSITIONS
PRETTY = {"NOR->ADE": "Normal → adenoma", "ADE->CORE": "Adenoma → carcinoma",
          "NOR->CORE": "Normal → carcinoma", "ADE->INV": "Adenoma → invasive front",
          "NOR->INV": "Normal → invasive front", "CORE->INV": "Tumour core → invasive front"}


def bh(p):
    p = np.asarray(p, dtype=float); o = np.argsort(p); q = np.empty_like(p)
    q[o] = np.minimum(1, np.minimum.accumulate((p[o] * len(p) / (np.arange(len(p)) + 1))[::-1])[::-1])
    return q


def analyze_geomx(risk):
    """Analyze staged GeoMx contrasts and export the support tables."""
    np.random.seed(20)
    gdir = risk.groupby("gene")["z_best"].mean()
    gmeta = (risk.groupby("gene")
             .agg(min_p=("p_best", "min"), n_ct=("cell_type", "size"),
                  tissues=("tissue", lambda s: "|".join(sorted(set(s)))),
                  cell_types=("cell_type", lambda s: "|".join(sorted(set(s)))))
             .reset_index())
    gmeta["twas_direction"] = np.where(gdir[gmeta["gene"]].values > 0,
                                       "risk-increasing", "protective")
    gmeta["expected_sign"] = np.where(gmeta["twas_direction"] == "risk-increasing", 1.0, -1.0)
    RISK = set(gmeta["gene"])
    print(f"risk genes: {len(RISK)} ({(gmeta.twas_direction=='risk-increasing').sum()} "
          f"risk-increasing, {(gmeta.twas_direction=='protective').sum()} protective)")

    de2 = pd.concat([pd.read_csv("spatial/de_exp2_epi.tsv.gz", sep="\t"),
                     pd.read_csv("spatial/de_exp2_str.tsv.gz", sep="\t")], ignore_index=True)
    de1 = pd.read_csv("spatial/de_exp1.tsv.gz", sep="\t")
    de2["experiment"] = "exp2"; de1["experiment"] = "exp1"
    DE = pd.concat([de2, de1], ignore_index=True)
    DE["is_risk"] = DE["gene"].isin(RISK)
    DE["compartment_label"] = DE["compartment"].map({"PanCK+": "Epithelium", "PanCK-": "Stroma"})
    print(f"DE rows: {len(DE)}; genes exp2 {de2.gene.nunique()}, exp1 {de1.gene.nunique()}")
    print(f"risk genes measured: exp2 {de2[de2.gene.isin(RISK)].gene.nunique()}, "
          f"exp1 {de1[de1.gene.isin(RISK)].gene.nunique()}")

    # ============================================================== A  per-gene table
    R = DE[DE["is_risk"]].merge(gmeta, on="gene", how="left").copy()
    R["signed_log2FC"] = R["log2FC"] * R["expected_sign"]
    R["concordant"] = R["signed_log2FC"] > 0
    # Two multiple-testing frames, and the distinction matters. The FDR shipped with the
    # published tables is Benjamini-Hochberg over all ~18,441 analysed genes. The hypothesis
    # space here is not the transcriptome: it is the set of TWAS risk genes that were
    # actually tested. FDR_risk therefore recomputes BH over the risk genes alone, within
    # each experiment x compartment x transition. That is the criterion used for the
    # per-gene calls and for the validated set; the genome-wide FDR is retained because the
    # risk-versus-background comparison below needs both groups on the same footing.
    R["FDR_risk"] = R.groupby(["experiment", "compartment", "transition"])["p"].transform(bh)
    R["sig_genomewide"] = R["FDR"] < 0.05
    R["sig"] = R["FDR_risk"] < 0.05
    R["verdict"] = np.where(~R["sig"], "not significant",
                            np.where(R["concordant"], "supports TWAS direction",
                                     "opposes TWAS direction"))
    R = R.sort_values(["experiment", "compartment", "transition", "FDR"])
    R.to_csv(GEOMX_RESULTS, index=False)
    print(f"\nper-gene rows written: {len(R)}")

    print("\n=== multiple-testing frame ===")
    for t in MAIN_TRANSITIONS:
        sl = R[(R.compartment == "PanCK+") & (R.transition == t)]
        print(f"  {t:11s} epithelium: genome-wide FDR {int((sl.FDR < 0.05).sum()):4d}   "
              f"risk-gene-restricted FDR {int(sl['sig'].sum()):4d}   of {len(sl)} tests")

    print("\n=== A. significant stage changes among risk genes (FDR<0.05, risk-restricted) ===")
    for exp in ["exp2", "exp1"]:
        for comp in ["PanCK+", "PanCK-"]:
            sub = R[(R.experiment == exp) & (R.compartment == comp)]
            if not len(sub):
                continue
            n_genes = sub.gene.nunique()
            line = []
            for t in ALL_TRANSITIONS:
                s = sub[sub.transition == t]
                line.append(f"{t}: {int(s.sig.sum()):4d}/{len(s):4d}")
            print(f"  {exp} {('Epithelium' if comp=='PanCK+' else 'Stroma'):11s} "
                  f"({n_genes} genes)  " + "   ".join(line))

    # ============================================================== B  vs background
    print("\n=== B. are risk genes differentially expressed more often than other genes? ===")
    rowsB = []
    for exp in ["exp2", "exp1"]:
        for comp in ["PanCK+", "PanCK-"]:
            for t in ALL_TRANSITIONS:
                d = DE[(DE.experiment == exp) & (DE.compartment == comp) & (DE.transition == t)]
                if not len(d):
                    continue
                a = int((d.is_risk & (d.FDR < 0.05)).sum())
                b = int((d.is_risk & (d.FDR >= 0.05)).sum())
                c = int((~d.is_risk & (d.FDR < 0.05)).sum())
                e = int((~d.is_risk & (d.FDR >= 0.05)).sum())
                orv, p = fisher_exact([[a, b], [c, e]], alternative="greater")
                u = mannwhitneyu(d.loc[d.is_risk, "log2FC"].abs(),
                                 d.loc[~d.is_risk, "log2FC"].abs(), alternative="greater")
                rowsB.append({"experiment": exp, "compartment": comp, "transition": t,
                              "n_risk": a + b, "frac_risk_sig": a / (a + b) if a + b else np.nan,
                              "n_other": c + e, "frac_other_sig": c / (c + e) if c + e else np.nan,
                              "odds_ratio": orv, "p_fisher": p,
                              "median_abs_lfc_risk": float(d.loc[d.is_risk, "log2FC"].abs().median()),
                              "median_abs_lfc_other": float(d.loc[~d.is_risk, "log2FC"].abs().median()),
                              "p_effectsize": float(u.pvalue)})
    B = pd.DataFrame(rowsB)
    B["fdr_fisher"] = bh(B["p_fisher"])
    B.to_csv("results/spatial_risk_vs_background.csv", index=False)
    for exp in ["exp2", "exp1"]:
        sub = B[B.experiment == exp]
        print(f"\n  {exp}")
        for r in sub.itertuples(index=False):
            print(f"    {('Epi' if r.compartment=='PanCK+' else 'Str')} {r.transition:11s} "
                  f"risk {100*r.frac_risk_sig:5.1f}% vs other {100*r.frac_other_sig:5.1f}%  "
                  f"OR={r.odds_ratio:5.2f}  P={r.p_fisher:.3g}  FDR={r.fdr_fisher:.3g}")

    # ============================================================== C  concordance
    print("\n=== C. direction concordance: do risk-increasing genes go up and protective down? ===")
    NPERM = 10000
    rowsC = []
    for exp in ["exp2", "exp1"]:
        for comp in ["PanCK+", "PanCK-"]:
            for t in ALL_TRANSITIONS:
                s = R[(R.experiment == exp) & (R.compartment == comp) & (R.transition == t)]
                s = s.dropna(subset=["signed_log2FC"])
                if len(s) < 20:
                    continue
                k = int(s.concordant.sum()); n = len(s)
                bt = binomtest(k, n, 0.5, alternative="greater")
                w = wilcoxon(s["signed_log2FC"], alternative="greater")
                # permutation of the TWAS direction labels among the same genes
                lf = s["log2FC"].values; sg = s["expected_sign"].values
                obs = float(np.mean(lf * sg))
                perm = np.array([np.mean(lf * np.random.permutation(sg)) for _ in range(NPERM)])
                p_perm = (np.sum(perm >= obs) + 1) / (NPERM + 1)
                # background genes given random directions, same n
                d = DE[(DE.experiment == exp) & (DE.compartment == comp) & (DE.transition == t)]
                bgl = d.loc[~d.is_risk, "log2FC"].values
                bgd = []
                for _ in range(2000):
                    idx = np.random.choice(len(bgl), n, replace=False)
                    bgd.append(np.mean(bgl[idx] * np.random.choice([1., -1.], n)))
                bgd = np.array(bgd)
                p_bg = (np.sum(bgd >= obs) + 1) / (len(bgd) + 1)
                # significant-only subset
                ss = s[s.sig]
                k2, n2 = int(ss.concordant.sum()), len(ss)
                bt2 = binomtest(k2, n2, 0.5, alternative="greater") if n2 >= 5 else None
                rowsC.append({"experiment": exp, "compartment": comp, "transition": t,
                              "n_genes": n, "n_concordant": k, "frac_concordant": k / n,
                              "p_binomial": bt.pvalue,
                              "median_signed_lfc": float(s["signed_log2FC"].median()),
                              "mean_signed_lfc": obs, "p_wilcoxon": float(w.pvalue),
                              "p_direction_permutation": p_perm,
                              "p_vs_background": p_bg,
                              "n_sig": n2, "n_sig_concordant": k2,
                              "frac_sig_concordant": (k2 / n2) if n2 else np.nan,
                              "p_binomial_sig_only": (bt2.pvalue if bt2 else np.nan)})
    C = pd.DataFrame(rowsC)
    for col in ["p_binomial", "p_direction_permutation", "p_vs_background"]:
        C["fdr_" + col.replace("p_", "")] = bh(C[col])
    C.to_csv("results/spatial_direction_concordance.csv", index=False)
    for exp in ["exp2", "exp1"]:
        print(f"\n  {exp}")
        for r in C[C.experiment == exp].itertuples(index=False):
            star = "  <<<" if (r.p_direction_permutation < 0.05 and r.frac_concordant > 0.5) else ""
            print(f"    {('Epi' if r.compartment=='PanCK+' else 'Str')} {r.transition:11s} "
                  f"concordant {r.n_concordant:3d}/{r.n_genes:3d} = {100*r.frac_concordant:5.1f}%  "
                  f"binom P={r.p_binomial:.3g}  perm P={r.p_direction_permutation:.4f}  "
                  f"vs-bg P={r.p_vs_background:.4f}{star}")

    # ============================================================== D  replication
    print("\n=== D. replication between the two independent GeoMx cohorts ===")
    rowsD = []
    for comp in ["PanCK+", "PanCK-"]:
        for t in ALL_TRANSITIONS:
            a = DE[(DE.experiment == "exp2") & (DE.compartment == comp) & (DE.transition == t)]
            b = DE[(DE.experiment == "exp1") & (DE.compartment == comp) & (DE.transition == t)]
            m = a.merge(b, on="gene", suffixes=("_2", "_1"))
            if len(m) < 50:
                continue
            rho = spearmanr(m["log2FC_2"], m["log2FC_1"])
            same = float(np.mean(np.sign(m["log2FC_2"]) == np.sign(m["log2FC_1"])))
            mr = m[m["gene"].isin(RISK)]
            rho_r = spearmanr(mr["log2FC_2"], mr["log2FC_1"]) if len(mr) >= 15 else None
            rowsD.append({"compartment": comp, "transition": t, "n_shared_genes": len(m),
                          "rho_all": rho.statistic, "p_all": rho.pvalue,
                          "sign_agreement_all": same,
                          "n_risk_shared": len(mr),
                          "rho_risk": rho_r.statistic if rho_r else np.nan,
                          "p_risk": rho_r.pvalue if rho_r else np.nan,
                          "sign_agreement_risk": float(np.mean(
                              np.sign(mr["log2FC_2"]) == np.sign(mr["log2FC_1"]))) if len(mr) else np.nan})
    D = pd.DataFrame(rowsD)
    D.to_csv("results/spatial_replication.csv", index=False)
    for r in D.itertuples(index=False):
        print(f"    {('Epi' if r.compartment=='PanCK+' else 'Str')} {r.transition:11s} "
              f"all genes rho={r.rho_all:+.3f} (n={r.n_shared_genes}, P={r.p_all:.2g}), "
              f"risk genes rho={r.rho_risk:+.3f} (n={r.n_risk_shared}), "
              f"sign agreement {100*r.sign_agreement_risk:.0f}%")

    # ============================================================== E  compartment
    print("\n=== E. epithelium versus stroma (the TWAS genes were called in epithelial cell types) ===")
    rowsE = []
    for exp in ["exp2", "exp1"]:
        for t in ALL_TRANSITIONS:
            e = R[(R.experiment == exp) & (R.compartment == "PanCK+") & (R.transition == t)]
            s = R[(R.experiment == exp) & (R.compartment == "PanCK-") & (R.transition == t)]
            if len(e) < 20 or len(s) < 20:
                continue
            u = mannwhitneyu(e["signed_log2FC"], s["signed_log2FC"], alternative="greater")
            rowsE.append({"experiment": exp, "transition": t,
                          "median_signed_lfc_epithelium": float(e["signed_log2FC"].median()),
                          "median_signed_lfc_stroma": float(s["signed_log2FC"].median()),
                          "frac_concordant_epithelium": float(e.concordant.mean()),
                          "frac_concordant_stroma": float(s.concordant.mean()),
                          "p": float(u.pvalue)})
    E = pd.DataFrame(rowsE)
    if len(E):
        E["fdr"] = bh(E["p"])
        E.to_csv("results/spatial_compartment_specificity.csv", index=False)
        for r in E.itertuples(index=False):
            print(f"    {r.experiment} {r.transition:11s} concordance epi "
                  f"{100*r.frac_concordant_epithelium:5.1f}% vs stroma "
                  f"{100*r.frac_concordant_stroma:5.1f}%,  P={r.p:.4f}")

    # ============================================================== headline genes
    print("\n=== strongest individual risk genes, epithelium, exp2 ===")
    top = (R[(R.experiment == "exp2") & (R.compartment == "PanCK+") & R.sig]
           .sort_values("FDR"))
    for t in MAIN_TRANSITIONS:
        s = top[top.transition == t].head(12)
        print(f"\n  {PRETTY[t]}  ({int((top.transition==t).sum())} significant risk genes)")
        for r in s.itertuples(index=False):
            print(f"    {r.gene:12s} {r.twas_direction:16s} log2FC={r.log2FC:+6.2f} "
                  f"FDR={r.FDR:.2e}  {r.verdict}  [{r.cell_types}]")

    summary = {
        "n_risk_genes": len(RISK),
        "n_measured_exp2": int(de2[de2.gene.isin(RISK)].gene.nunique()),
        "n_measured_exp1": int(de1[de1.gene.isin(RISK)].gene.nunique()),
        "sig_counts": {f"{r.experiment}|{r.compartment}|{r.transition}":
                       {"n_sig": int(R[(R.experiment == r.experiment) &
                                       (R.compartment == r.compartment) &
                                       (R.transition == r.transition)].sig.sum()),
                        "n_tested": int(len(R[(R.experiment == r.experiment) &
                                              (R.compartment == r.compartment) &
                                              (R.transition == r.transition)]))}
                       for r in C.itertuples(index=False)},
        "concordance": C.to_dict("records"),
        "vs_background": B.to_dict("records"),
        "replication": D.to_dict("records"),
        "compartment": E.to_dict("records") if len(E) else [],
    }
    json.dump(summary, open("results/spatial_de_summary.json", "w"), indent=1, default=float)


def select_validated_genes(de, risk):
    """Export validation results from the per-gene GeoMx analysis table."""
    np.random.seed(23)
    transition_names = {"NOR->ADE": "Normal -> adenoma", "ADE->CORE": "Adenoma -> carcinoma",
              "NOR->CORE": "Normal -> carcinoma", "ADE->INV": "Adenoma -> invasive front",
              "NOR->INV": "Normal -> invasive front", "CORE->INV": "Core -> invasive front"}

    epi = de[(de["compartment"] == "PanCK+") & de["is_risk"]].copy()
    print(f"risk-gene epithelial contrasts: {len(epi):,}  genes: {epi.gene.nunique()}")

    # ---------------------------------------------------------------- the rule
    def best_lfc(m, t):
        """Fold change from the most significant panel that measured this transition."""
        sub = m[m["transition"] == t]
        if not len(sub):
            return float("nan")
        return float(sub.sort_values("FDR")["log2FC"].iloc[0])


    rows = []
    for gene, g in epi.groupby("gene"):
        m = g[g["transition"].isin(MAIN_TRANSITIONS)]
        if not len(m):
            continue
        sig = m[m["sig"]]
        sup = sig[sig["concordant"]]
        opp = sig[~sig["concordant"]]
        core = sig[(sig["transition"] == "NOR->CORE") & sig["experiment"].isin(["exp1", "exp2"])]
        twas_dir = g["twas_direction"].iloc[0]
        passes = len(sup) > 0 and len(opp) == 0

        # the transition that carries the claim: strongest supporting one
        lead = sup.sort_values("FDR").iloc[0] if len(sup) else None
        inv = g[g["transition"].isin(INVASION_TRANSITIONS)]
        inv_sig_sup = int(((inv["sig"]) & (inv["concordant"])).sum())
        inv_sig_opp = int(((inv["sig"]) & (~inv["concordant"])).sum())

        # does it move the same way at every measured step of the main sequence?
        monotone = bool((m["signed_log2FC"] > 0).all())

        ct = risk[risk["gene"] == gene]
        rows.append({
            "gene": gene,
            "twas_direction": twas_dir,
            "class": ("oncogene-like" if twas_dir == "risk-increasing" else "tumour-suppressor-like")
                     if passes else "not validated",
            "validated": passes,
            "nor_core_supported": bool(core["concordant"].any() and core["concordant"].all()),
            "n_sig_supporting": int(len(sup)),
            "n_sig_opposing": int(len(opp)),
            "lead_transition": lead["transition"] if lead is not None else "",
            "lead_log2FC": float(lead["log2FC"]) if lead is not None else np.nan,
            "lead_FDR": float(lead["FDR"]) if lead is not None else np.nan,
            "lead_experiment": lead["experiment"] if lead is not None else "",
            "monotone_across_sequence": monotone,
            "invasion_supporting": inv_sig_sup,
            "invasion_opposing": inv_sig_opp,
            # where both panels measured a gene, quote the more significant one rather than
            # an average, so every number in the report refers to one actual measurement
            "log2FC_NOR_ADE": best_lfc(m, "NOR->ADE"),
            "log2FC_ADE_CORE": best_lfc(m, "ADE->CORE"),
            "log2FC_NOR_CORE": best_lfc(m, "NOR->CORE"),
            "FDR_NOR_CORE": float(m.loc[m.transition == "NOR->CORE", "FDR"].min()),
            "twas_min_p": float(g["min_p"].iloc[0]),
            "n_celltypes": int(g["n_ct"].iloc[0]),
            "tissues": g["tissues"].iloc[0],
            "cell_types": g["cell_types"].iloc[0],
            "experiments": "|".join(sorted(g["experiment"].unique())),
        })

    V = pd.DataFrame(rows)
    val = V[V["validated"]].copy()
    print(f"\nmeasured genes: {len(V)}   validated: {len(val)}   "
          f"({100*len(val)/len(V):.1f} %)")
    print(val["class"].value_counts().to_string())
    print(f"monotone across the whole sequence: {int(val['monotone_across_sequence'].sum())}")

    # ---------------------------------------------------------------- enrichment of the rule
    # How many genes would pass by chance? Randomise the TWAS direction labels among the
    # same measured genes and re-apply the rule. This is the honest null: the significance
    # pattern is fixed, only the direction assignment moves.
    mm = epi[epi["transition"].isin(MAIN_TRANSITIONS)].copy()
    # one row per gene x transition: the most significant experiment, and "significant"
    # meaning significant on any panel that measured it - exactly what the rule above used
    best = mm.sort_values("FDR").groupby(["gene", "transition"], as_index=False).first()
    anysig = mm.groupby(["gene", "transition"])["sig"].any().rename("anysig").reset_index()
    best = best.drop(columns=["sig"]).merge(anysig, on=["gene", "transition"])
    genes = np.array(sorted(best["gene"].unique()))
    S = np.zeros((len(genes), len(MAIN_TRANSITIONS)), dtype=bool)
    L = np.zeros((len(genes), len(MAIN_TRANSITIONS)), dtype=float)
    gpos = {g: k for k, g in enumerate(genes)}
    for j, t in enumerate(MAIN_TRANSITIONS):
        sl = best[best["transition"] == t]
        idx = sl["gene"].map(gpos).to_numpy()
        S[idx, j] = sl["anysig"].to_numpy()
        L[idx, j] = sl["log2FC"].to_numpy(float)
    obs_sign = V.set_index("gene").loc[genes, "twas_direction"].map(
        {"risk-increasing": 1.0, "protective": -1.0}).to_numpy()


    def n_pass(sign):
        conc = (L * sign[:, None]) > 0
        sup = (S & conc).any(1)
        opp = (S & ~conc).any(1)
        return int((sup & ~opp).sum())


    obs = n_pass(obs_sign)
    perm = np.array([n_pass(np.random.permutation(obs_sign)) for _ in range(10000)])
    p_perm = float((perm >= obs).mean())
    print(f"\nvalidated under observed directions: {obs}   "
          f"permuted mean {perm.mean():.1f} (sd {perm.std():.1f})   P = {p_perm:.4f}")

    # ---------------------------------------------------------------- per TWAS cell type
    rr = risk[["gene", "tissue", "cell_type", "direction"]].drop_duplicates()
    j = rr.merge(V[["gene", "validated", "class", "lead_log2FC"]], on="gene", how="inner")
    ct_rows = []
    for (tis, ct), g in j.groupby(["tissue", "cell_type"]):
        n, k = len(g), int(g["validated"].sum())
        onc = int((g["class"] == "oncogene-like").sum())
        tsg = int((g["class"] == "tumour-suppressor-like").sum())
        # is this cell type's validation rate different from the rest?
        rest = j[~j.index.isin(g.index)]
        tab = [[k, n - k], [int(rest["validated"].sum()), len(rest) - int(rest["validated"].sum())]]
        orr, pf = fisher_exact(tab)
        ct_rows.append({"tissue": tis, "cell_type": ct, "n_measured": n, "n_validated": k,
                        "frac_validated": k / n if n else np.nan,
                        "oncogene_like": onc, "tumour_suppressor_like": tsg,
                        "odds_ratio_vs_rest": float(orr), "p_fisher": float(pf)})
    CT = pd.DataFrame(ct_rows).sort_values("frac_validated", ascending=False)
    CT["fdr"] = np.minimum(1, CT["p_fisher"] * len(CT))
    print("\n=== validation rate per TWAS cell type (epithelium) ===")
    print(CT.to_string(index=False))

    V.to_csv("results/validated_risk_genes_all.csv", index=False)
    val.sort_values("lead_FDR").to_csv("results/validated_risk_genes.csv", index=False)
    CT.to_csv("results/validated_by_celltype.csv", index=False)

    top = val.sort_values("lead_FDR").head(25)
    print("\n=== the validated set, strongest first ===")
    print(top[["gene", "class", "lead_transition", "lead_log2FC", "lead_FDR",
               "cell_types", "monotone_across_sequence"]].to_string(index=False))

    json.dump({
        "rule": "FDR<0.05 and TWAS-concordant in >=1 of NOR->ADE / ADE->CORE / NOR->CORE, "
                "no significant opposite change in the other two; epithelium (PanCK+)",
        "n_measured": int(len(V)), "n_validated": int(len(val)),
        "frac_validated": float(len(val) / len(V)),
        "n_oncogene_like": int((val["class"] == "oncogene-like").sum()),
        "n_tsg_like": int((val["class"] == "tumour-suppressor-like").sum()),
        "n_monotone": int(val["monotone_across_sequence"].sum()),
        "n_excluded_by_contradiction": int(((V["n_sig_supporting"] > 0) & (V["n_sig_opposing"] > 0)).sum()),
        "n_never_significant": int(((V["n_sig_supporting"] == 0) & (V["n_sig_opposing"] == 0)).sum()),
        "n_sig_opposite_only": int(((V["n_sig_supporting"] == 0) & (V["n_sig_opposing"] > 0)).sum()),
        "permutation": {"observed": obs, "null_mean": float(perm.mean()),
                        "null_sd": float(perm.std()), "p": p_perm, "n_perm": 10000},
        "lead_transition_counts": {transition_names[k]: int(v) for k, v in
                                   val["lead_transition"].value_counts().items()},
        "by_celltype": CT.to_dict("records"),
        "validated_genes": val.sort_values("lead_FDR")[
            ["gene", "class", "lead_transition", "lead_log2FC", "lead_FDR", "cell_types",
             "tissues", "monotone_across_sequence", "twas_min_p"]].to_dict("records"),
    }, open("results/validated_summary.json", "w"), indent=1, default=float)

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--validation-only", action="store_true",
                        help="Read the existing GeoMx result CSV and refresh validation outputs.")
    args = parser.parse_args()

    risk = pd.read_csv(RISK_GENES)
    if not args.validation_only:
        analyze_geomx(risk)
    # Read the exported CSV in both modes so validation uses identical inputs.
    select_validated_genes(pd.read_csv(GEOMX_RESULTS), risk)


if __name__ == "__main__":
    main()
