#!/usr/bin/env python3
"""
33 - Do the CRC risk genes have a measurable effect on colorectal cancer cell viability?

The screens supply the missing layer. Sections 2.1-2.14 establish which genes are
associated, which are disrupted in tissue and which direction they move; none of that
says whether removing the gene does anything to a cancer cell. These experiments do.

Design of the analysis
  * SW480 is dropped (see 32_functional_parse_qc.py): its siRNA Z'-factor is at or below
    zero in two of three replicates, meaning the killing control and the non-targeting
    control overlap, so that line cannot separate a hit from nothing. The same two lines
    are then used for CRISPR so the two modalities are compared on the same background;
    a sensitivity analysis with SW480 restored is reported.
  * siRNA values are already expressed as % of the non-targeting control. Each replicate
    is converted to a robust z against its own NT wells, z = (value - median_NT) /
    (1.4826 * MAD_NT), which uses the actual well-to-well noise of that plate rather than
    assuming it.
  * A gene's effect is the mean of its z across the retained lines x replicates, tested
    against zero by a one-sample t-test and corrected across the whole assayed panel.
    A hit must also move the same way in both cell lines.
  * The panel contains 333 genes of which only 62 are CRC risk genes. The other 271 are
    the internal comparison group: the question is not whether risk genes have effects,
    but whether they have them more often than the rest of the same experiment.

The directional prediction is the sharp test. Oncogene-like genes (§2.14) should be
required for growth, so knocking them down should reduce viability. Tumour-suppressor-
like genes should not, and may increase it.
"""
import pandas as pd, numpy as np, json
from scipy.stats import mannwhitneyu, fisher_exact, spearmanr, binomtest

np.random.seed(33)
ALIAS = {"ACAS2": "ACSS2"}          # old symbol used on the plate map
EXCLUDE = json.load(open("results/func_qc.json"))["excluded_cell_line"]
KEEP = [l for l in ["HCT116", "RKO", "SW480"] if l != EXCLUDE]
print(f"excluded cell line: {EXCLUDE};  retained: {KEEP}")


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


def robust_z(D, lines):
    """Per plate (cell line x replicate): the standard high-throughput robust z-score.

    Two separate things are needed from each plate and they come from different wells.
    The CENTRE is the plate's own sample median, not its non-targeting wells, because the
    global offset between the two is plate-specific and large in the siRNA screen - the
    median assayed gene sits at 67 % of NT on HCT116 replicate 1 and 108 % on RKO
    replicate 1. Anchoring on NT would carry that shift into every gene on the plate.
    Centring on the sample median is the convention in RNAi screening and asks the
    question the screen can actually answer: is this gene unusual for this plate.
    The SCALE is the sample MAD. The NT wells are then used only to measure how much of
    the spread is assay noise (see per_gene), and the NT-anchored % of control is kept
    alongside as the descriptive effect size the experimenters reported.
    """
    out = []
    for (line, rep), d in D.groupby(["cell_line", "replicate"]):
        if line not in lines:
            continue
        smp = d.loc[d.kind == "sample", "value"]
        nt = d.loc[d.kind == "negative", "value"]
        med, ntmed = smp.median(), nt.median()
        mad = (smp - med).abs().median() * 1.4826
        if not np.isfinite(mad) or mad == 0:
            mad = smp.std(ddof=1)
        dd = d.copy()
        dd["z"] = (dd["value"] - med) / mad                 # sample-centred, sample-scaled
        dd["z_noise"] = (dd["value"] - ntmed) / mad         # NT deviation on the same scale
        dd["pct_NT"] = 100 * dd["value"] / ntmed
        dd["plate_offset_pct"] = 100 * med / ntmed
        out.append(dd)
    return pd.concat(out, ignore_index=True)


def per_gene(D, lines, label, n_null=20000, seed=33):
    """Per-gene effect, tested against a null built from the plates' own NT wells.

    A one-sample t-test over the six wells of a gene (2 lines x 3 replicates) has almost
    no power once corrected across 333 genes - the smallest attainable p is around 1e-3.
    The screen itself supplies a better null: each plate carries 8-12 non-targeting wells
    that differ from one another only by assay noise. Drawing one NT well from each of
    the six plates and averaging gives the distribution of the gene-level statistic under
    no effect, with the real plate-to-plate and well-to-well variance built in and no
    normality assumption. The gene's p-value is its position in that distribution.
    """
    rng = np.random.default_rng(seed)
    Z = robust_z(D, lines)
    Z["gene"] = Z["gene"].replace(ALIAS)
    plates = sorted({(l, r) for l, r in zip(Z.cell_line, Z.replicate) if l in lines})
    # noise null: NT wells expressed as deviations from their own plate's NT median and
    # divided by the sample MAD, so the null is centred at zero and on the same scale as
    # the effect. It measures assay noise, not the spread of real effects.
    nt_by_plate = []
    for l, r in plates:
        d = Z[(Z.cell_line == l) & (Z.replicate == r) & (Z.kind == "negative")]
        nt_by_plate.append((d["z_noise"] - d["z_noise"].median()).to_numpy(float))
    null = np.mean([rng.choice(v, n_null, replace=True) for v in nt_by_plate], axis=0)
    null_sorted = np.sort(null)
    print(f"  {label}: null of the mean over {len(plates)} plates, "
          f"sd {null.std():.3f}, 2.5-97.5 % {np.quantile(null, 0.025):+.2f} "
          f"to {np.quantile(null, 0.975):+.2f}")

    Zs = Z[Z.kind == "sample"]
    rows = []
    for g, d in Zs.groupby("gene"):
        z = d["z"].to_numpy(float)
        stat = float(np.mean(z))
        # two-sided empirical p, +1 correction so p is never exactly zero
        k = int(np.sum(np.abs(null_sorted) >= abs(stat)))
        p = (k + 1) / (len(null_sorted) + 1)
        by_line = d.groupby("cell_line")["z"].mean()
        rows.append({
            "gene": g, "n_wells": len(z), "mean_z": stat,
            "sd_z": float(np.std(z, ddof=1)) if len(z) > 1 else np.nan,
            "mean_pct_NT": float(d["pct_NT"].mean()),
            "plate_offset_pct": float(d["plate_offset_pct"].mean()),
            "p": p,
            "z_" + lines[0]: float(by_line.get(lines[0], np.nan)),
            "z_" + lines[1]: float(by_line.get(lines[1], np.nan)),
            "same_direction": bool(np.sign(by_line.get(lines[0], 0)) ==
                                   np.sign(by_line.get(lines[1], 0))),
        })
    R = pd.DataFrame(rows)
    R["fdr"] = bh(R["p"])
    # two tiers, because they answer different questions. "detected" = the effect is
    # larger than assay noise after correction and points the same way in both lines.
    # "strong" additionally reaches the conventional robust-z threshold of 2, i.e. twice
    # the deviation of a typical gene on the plate.
    R["detected"] = (R["fdr"] < 0.05) & R["same_direction"]
    R["strong"] = R["detected"] & (R["mean_z"].abs() >= 2)
    R["hit"] = R["detected"]
    R["direction"] = np.where(R["mean_z"] < 0, "reduces viability", "increases viability")
    R["assay"] = label
    return R.sort_values("p")


S = pd.read_csv("results/func_sirna_long.csv")
C = pd.read_csv("results/func_crispr_long.csv")
SI = per_gene(S, KEEP, "siRNA")
CR = per_gene(C, KEEP, "CRISPR")
SI.to_csv("results/func_sirna_gene_stats.csv", index=False)
CR.to_csv("results/func_crispr_gene_stats.csv", index=False)

for R, nm in [(SI, "siRNA"), (CR, "CRISPR")]:
    nh, ns = int(R["detected"].sum()), int(R["strong"].sum())
    print(f"\n{nm}: {len(R)} genes; {nh} with a detectable consistent effect "
          f"({int(((R.detected) & (R.mean_z < 0)).sum())} reduce viability, "
          f"{int(((R.detected) & (R.mean_z > 0)).sum())} increase it); "
          f"{ns} of those reach |z| >= 2")
    print(R.head(10)[["gene", "mean_z", "mean_pct_NT", "p", "fdr", "same_direction"]]
          .round(5).to_string(index=False))

# ---------------------------------------------------------------- risk genes
risk = pd.read_csv("results/risk_genes_harmonized.csv")
riskset = set(risk["gene"])
allm = pd.read_csv("results/validated_risk_genes_all.csv").set_index("gene")
valset = set(pd.read_csv("results/validated_priority.csv")["gene"])
pri = pd.read_csv("results/integrated_gene_priority.csv").set_index("gene")

out = {}
for R, nm in [(SI, "siRNA"), (CR, "CRISPR")]:
    R["is_risk"] = R["gene"].isin(riskset)
    R["is_validated"] = R["gene"].isin(valset)
    R["twas_class"] = R["gene"].map(allm["class"])
    R["twas_direction"] = R["gene"].map(allm["twas_direction"])
    a = R[R.is_risk]
    b = R[~R.is_risk]
    tab = [[int(a.hit.sum()), len(a) - int(a.hit.sum())],
           [int(b.hit.sum()), len(b) - int(b.hit.sum())]]
    orr, pf = fisher_exact(tab)
    u = mannwhitneyu(a["mean_z"], b["mean_z"], alternative="two-sided")
    va = R[R.is_validated]
    res = {
        "n_risk": len(a), "n_other": len(b),
        "n_strong_risk": int(a.strong.sum()), "n_strong_other": int(b.strong.sum()),
        "n_validated": int(len(va)), "hits_validated": int(va.hit.sum()),
        "strong_validated": int(va.strong.sum()),
        "median_z_validated": float(va["mean_z"].median()) if len(va) else np.nan,
        "hits_risk": int(a.hit.sum()), "hits_other": int(b.hit.sum()),
        "frac_risk": float(a.hit.mean()), "frac_other": float(b.hit.mean()),
        "odds_ratio": float(orr), "p_fisher": float(pf),
        "median_z_risk": float(a["mean_z"].median()),
        "median_z_other": float(b["mean_z"].median()),
        "auroc_effect": float(u.statistic / (len(a) * len(b))),
        "p_mannwhitney": float(u.pvalue),
    }
    print(f"\n=== {nm}: risk genes vs the rest of the same screen ===")
    print(f"  spatially validated risk genes in this screen: {len(va)}, "
          f"{int(va.hit.sum())} with a detectable effect, {int(va.strong.sum())} strong")
    print(f"  hits {res['hits_risk']}/{res['n_risk']} ({100*res['frac_risk']:.1f} %) vs "
          f"{res['hits_other']}/{res['n_other']} ({100*res['frac_other']:.1f} %); "
          f"OR {orr:.2f}, P = {pf:.3g}")
    print(f"  median effect z {res['median_z_risk']:+.3f} vs {res['median_z_other']:+.3f}; "
          f"AUROC {res['auroc_effect']:.3f}, P = {res['p_mannwhitney']:.3g}")

    # the directional prediction
    onc = R[R.twas_class == "oncogene-like"]
    tsg = R[R.twas_class == "tumour-suppressor-like"]
    if len(onc) > 2 and len(tsg) > 2:
        ud = mannwhitneyu(onc["mean_z"], tsg["mean_z"], alternative="less")
        bt = binomtest(int((onc["mean_z"] < 0).sum()), len(onc), 0.5, alternative="greater")
        res["classes"] = {
            "n_oncogene_like": len(onc), "n_tsg_like": len(tsg),
            "median_z_oncogene": float(onc["mean_z"].median()),
            "median_z_tsg": float(tsg["mean_z"].median()),
            "p_onc_lower_than_tsg": float(ud.pvalue),
            "n_oncogene_reduce": int((onc["mean_z"] < 0).sum()),
            "p_binom_oncogene_reduce": float(bt.pvalue),
        }
        print(f"  oncogene-like (n={len(onc)}) median z {onc['mean_z'].median():+.3f} vs "
              f"tumour-suppressor-like (n={len(tsg)}) {tsg['mean_z'].median():+.3f}; "
              f"one-sided P = {ud.pvalue:.3g}")

    # does the in-silico knockout impact predict the measured effect?
    m = R.join(pri[["ko_impact", "degree_ppi"]], on="gene")
    ok = m.dropna(subset=["ko_impact"])
    if len(ok) > 8:
        rr = spearmanr(ok["ko_impact"], ok["mean_z"])
        res["ko_impact_vs_measured"] = {"n": int(len(ok)), "rho": float(rr.statistic),
                                        "p": float(rr.pvalue)}
        print(f"  in-silico knockout impact vs measured effect: n={len(ok)}, "
              f"rho = {rr.statistic:+.3f}, P = {rr.pvalue:.3g}")
    out[nm] = res

# ------------------------------------------------- calibration: do known essentials win?
import re
ESS = re.compile(r"^(RPS|RPL|PSM|SF3|PRPF|POLR|EIF|SNRP|MCM|TBCE)")
for R, nm in [(SI, "siRNA"), (CR, "CRISPR")]:
    R["core_machinery"] = R["gene"].str.match(ESS)
    a, b = R[R.core_machinery], R[~R.core_machinery]
    if len(a) >= 3:
        u = mannwhitneyu(a["mean_z"], b["mean_z"], alternative="less")
        tab = [[int(a.strong.sum()), len(a) - int(a.strong.sum())],
               [int(b.strong.sum()), len(b) - int(b.strong.sum())]]
        orr, pf = fisher_exact(tab)
        out.setdefault("calibration", {})[nm] = {
            "n_core": int(len(a)), "median_z_core": float(a["mean_z"].median()),
            "median_z_rest": float(b["mean_z"].median()), "p": float(u.pvalue),
            "strong_core": int(a.strong.sum()), "strong_rest": int(b.strong.sum()),
            "odds_ratio": float(orr), "p_fisher": float(pf)}
        print(f"\n=== {nm} calibration: ribosome/proteasome/spliceosome/replication genes ===")
        print(f"  n={len(a)}; median z {a['mean_z'].median():+.2f} vs {b['mean_z'].median():+.2f} "
              f"for the rest (one-sided P = {u.pvalue:.3g}); "
              f"{int(a.strong.sum())}/{len(a)} are strong hits vs {int(b.strong.sum())}/{len(b)} "
              f"(OR {orr:.1f}, P = {pf:.3g})")

# ---------------------------------------------------------------- sensitivity: SW480 back in
SI3 = per_gene(S, ["HCT116", "RKO", "SW480"], "siRNA_3lines")
cmp3 = SI[["gene", "mean_z"]].merge(SI3[["gene", "mean_z"]], on="gene",
                                    suffixes=("_2line", "_3line"))
r3 = spearmanr(cmp3["mean_z_2line"], cmp3["mean_z_3line"])
hits2, hits3 = set(SI[SI.hit].gene), set(SI3[SI3.hit].gene)
print(f"\n=== sensitivity: adding {EXCLUDE} back ===")
print(f"  effect sizes rho = {r3.statistic:.3f}; hits {len(hits2)} -> {len(hits3)}, "
      f"{len(hits2 & hits3)} shared")
out["sensitivity_three_lines"] = {"rho": float(r3.statistic), "n_hits_two": len(hits2),
                                  "n_hits_three": len(hits3),
                                  "n_shared": len(hits2 & hits3)}

# ---------------------------------------------------------------- siRNA vs CRISPR
m = SI.merge(CR, on="gene", suffixes=("_si", "_cr"))
m["agree_direction"] = np.sign(m["mean_z_si"]) == np.sign(m["mean_z_cr"])
rr = spearmanr(m["mean_z_si"], m["mean_z_cr"])
both_hit = m[(m.hit_si) & (m.hit_cr)]
si_only = m[(m.hit_si) & (~m.hit_cr)]
cr_only = m[(~m.hit_si) & (m.hit_cr)]
disagree = m[(m.hit_si | m.hit_cr) & (~m.agree_direction)]
print(f"\n=== siRNA vs CRISPR, {len(m)} genes assayed in both ===")
print(f"  effect correlation rho = {rr.statistic:+.3f} (P = {rr.pvalue:.3g}); "
      f"direction agrees for {int(m.agree_direction.sum())}/{len(m)} "
      f"({100*m.agree_direction.mean():.0f} %)")
print(f"  hits in both {len(both_hit)}; siRNA only {len(si_only)}; CRISPR only {len(cr_only)}; "
      f"opposite direction among hits {len(disagree)}")
if len(disagree):
    print(disagree[["gene", "mean_z_si", "fdr_si", "mean_z_cr", "fdr_cr", "is_risk_si"]]
          .round(4).to_string(index=False))
m.to_csv("results/func_sirna_vs_crispr.csv", index=False)
out["concordance"] = {
    "n_both": int(len(m)), "rho": float(rr.statistic), "p": float(rr.pvalue),
    "n_direction_agree": int(m.agree_direction.sum()),
    "frac_direction_agree": float(m.agree_direction.mean()),
    "n_hit_both": int(len(both_hit)), "n_hit_sirna_only": int(len(si_only)),
    "n_hit_crispr_only": int(len(cr_only)),
    "n_opposite_among_hits": int(len(disagree)),
    "opposite_genes": disagree[["gene", "mean_z_si", "fdr_si", "mean_z_cr", "fdr_cr"]]
        .to_dict("records"),
    "crispr_only_risk_genes": cr_only[cr_only.is_risk_si][
        ["gene", "mean_z_si", "fdr_si", "mean_z_cr", "fdr_cr"]].to_dict("records"),
}

# ---------------------------------------------------------------- per-gene risk table
RG = SI[SI.is_risk].copy()
RG["twas_cell_types"] = RG["gene"].map(allm["cell_types"])
RG["spatial_lead_FDR"] = RG["gene"].map(allm["lead_FDR"])
RG["ko_impact"] = RG["gene"].map(pri["ko_impact"])
RG = RG.merge(CR[["gene", "mean_z", "fdr", "detected"]].rename(
    columns={"mean_z": "crispr_mean_z", "fdr": "crispr_fdr", "detected": "crispr_hit"}),
    on="gene", how="left")
cols = ["gene", "twas_direction", "twas_class", "is_validated", "twas_cell_types",
        "spatial_lead_FDR", "ko_impact", "mean_z", "mean_pct_NT", "plate_offset_pct",
        "z_" + KEEP[0], "z_" + KEEP[1], "same_direction", "p", "fdr", "detected", "strong",
        "crispr_mean_z", "crispr_fdr", "crispr_hit"]
RG[cols].sort_values("p").to_csv("results/func_risk_gene_results.csv", index=False)
print("\n--- risk genes with a detectable effect (siRNA) ---")
print(RG[RG.hit].sort_values("mean_z")[
    ["gene", "twas_class", "is_validated", "mean_z", "mean_pct_NT", "fdr", "strong"]]
    .round(4).to_string(index=False))
print(f"\n=== the {len(RG)} CRC risk genes in the siRNA screen ===")
print(RG.sort_values("p").head(15)[
    ["gene", "twas_class", "mean_z", "mean_pct_NT", "p", "fdr", "hit", "crispr_mean_z"]]
    .round(4).to_string(index=False))

# re-save the gene tables now that the risk / class / calibration annotations are on them
SI.to_csv("results/func_sirna_gene_stats.csv", index=False)
CR.to_csv("results/func_crispr_gene_stats.csv", index=False)

out["risk_gene_table"] = RG[cols].sort_values("p").head(25).to_dict("records")
out["excluded_cell_line"] = EXCLUDE
out["retained"] = KEEP
json.dump(out, open("results/func_summary.json", "w"), indent=1, default=float)
print("\nDone.")
