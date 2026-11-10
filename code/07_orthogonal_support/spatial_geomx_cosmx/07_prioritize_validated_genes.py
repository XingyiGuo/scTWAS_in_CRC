#!/usr/bin/env python3
"""
The drug-screening results restricted to the spatially validated genes.

Everything in sections 2.5-2.9 of the report was computed on all 347 risk genes. This
re-runs the ranking on the 62 genes whose somatic direction in tissue agrees with the
genetic direction, which is the set the user asked to report, and asks two further
questions:

  * are the validated genes any more tractable, or any more disruptive when knocked out,
    than the risk genes that failed validation? (If not, validation buys biological
    interpretability, not druggability - and that should be said.)
  * for each validated gene, does a direction-matched drug exist? An oncogene-like gene
    needs an inhibitor; a tumour-suppressor-like gene needs restoration, which small
    molecules rarely achieve, so the two classes are not equally actionable and the
    table says so per gene.

Outputs
  results/validated_priority.csv        the ranked validated set
  results/validated_vs_failed.json      validated against non-validated risk genes
"""
import pandas as pd, numpy as np, json
from scipy.stats import mannwhitneyu, fisher_exact

np.random.seed(27)

V = pd.read_csv("results/validated_risk_genes_all.csv")
pri = pd.read_csv("results/integrated_gene_priority.csv")
cos = pd.read_csv("results/cosmx_crc_vs_normal_epi.csv")

d = V.merge(pri, on="gene", how="left", suffixes=("", "_pri"))
d = d.merge(cos[["gene", "delta", "p", "concordant"]].rename(
    columns={"delta": "cosmx_delta", "p": "cosmx_p", "concordant": "cosmx_concordant"}),
    on="gene", how="left")

val = d[d["validated"]].copy()
fail = d[~d["validated"]].copy()
print(f"validated {len(val)}   not validated {len(fail)}")

# ---------------------------------------------------------------- therapeutic reading
def action(row):
    if row["class"] == "oncogene-like":
        return "inhibit"
    if row["class"] == "tumour-suppressor-like":
        return "restore"
    return ""


val["required_action"] = val.apply(action, axis=1)
val["direction_matched_drug"] = np.where(
    (val["required_action"] == "inhibit") & (val["matched_approved_drugs"].notna())
    & (val["matched_approved_drugs"].astype(str) != ""),
    val["matched_approved_drugs"].astype(str), "")
tier_num = lambda s: pd.to_numeric(s.astype(str).str.extract(r"^T(\d)")[0], errors="coerce")
val["tier_n"] = tier_num(val["tier"])
fail_tier = tier_num(d.loc[~d["validated"], "tier"])
val["actionable_now"] = (val["required_action"] == "inhibit") & (val["tier_n"] <= 2)

# composite score: genetic strength, knockout impact, and tractability, as in section 2.9,
# but the direction term is replaced by the spatial evidence, which is now measured
def z(s):
    s = pd.to_numeric(s, errors="coerce")
    return (s - s.mean()) / (s.std() if s.std() else 1)


val["score_genetic"] = z(-np.log10(val["twas_min_p"].clip(lower=1e-300)))
val["score_ko"] = z(val["ko_impact"])
val["score_spatial"] = z(-np.log10(val["lead_FDR"].clip(lower=1e-300)))
val["score_drug"] = z(val["tractability_z"].fillna(val["tractability_z"].median()))
val["priority"] = (0.30 * val["score_genetic"] + 0.25 * val["score_ko"]
                   + 0.25 * val["score_spatial"] + 0.20 * val["score_drug"])
val = val.sort_values("priority", ascending=False)

cols = ["gene", "class", "required_action", "cell_types", "tissues", "twas_min_p",
        "lead_transition", "lead_log2FC", "lead_FDR", "monotone_across_sequence",
        "cosmx_delta", "cosmx_concordant", "ko_impact", "tier", "n_drugbank_drugs",
        "matched_approved_drugs", "tractability_z", "priority"]
val[cols].to_csv("results/validated_priority.csv", index=False)
print("\n=== validated genes, ranked ===")
print(val[["gene", "class", "required_action", "lead_log2FC", "lead_FDR", "ko_impact",
           "tier", "matched_approved_drugs", "priority"]].head(20).to_string(index=False))

# ---------------------------------------------------------------- validated vs failed
out = {}
for name, col in [("ko_impact", "ko_impact"), ("tractability_z", "tractability_z"),
                  ("degree_ppi", "degree_ppi")]:
    a = pd.to_numeric(val[col], errors="coerce").dropna()
    b = pd.to_numeric(fail[col], errors="coerce").dropna()
    if len(a) > 3 and len(b) > 3:
        u = mannwhitneyu(a, b, alternative="two-sided")
        auroc = float(u.statistic / (len(a) * len(b)))
        out[name] = {"median_validated": float(a.median()), "median_failed": float(b.median()),
                     "auroc": auroc, "p": float(u.pvalue)}
        print(f"{name:16s} validated {a.median():.3f} vs failed {b.median():.3f}  "
              f"AUROC {auroc:.3f}  P = {u.pvalue:.3g}")

t1 = int((val["tier_n"] <= 2).sum()); t1f = int((fail_tier <= 2).sum())
orr, pf = fisher_exact([[t1, len(val) - t1], [t1f, len(fail) - t1f]])
out["tier1_2"] = {"validated": t1, "n_validated": int(len(val)),
                  "failed": t1f, "n_failed": int(len(fail)),
                  "odds_ratio": float(orr), "p": float(pf)}
print(f"\nTier 1-2 druggable: {t1}/{len(val)} validated vs {t1f}/{len(fail)} not "
      f"(OR {orr:.2f}, P = {pf:.3g})")

onc = val[val["class"] == "oncogene-like"]
tsg = val[val["class"] == "tumour-suppressor-like"]
act = val[val["actionable_now"]]
print(f"\noncogene-like {len(onc)}, tumour-suppressor-like {len(tsg)}")
print(f"direction-matched and already druggable (inhibit, tier <= 2): {len(act)}")
if len(act):
    print(act[["gene", "cell_types", "tier", "matched_approved_drugs"]].to_string(index=False))

out["classes"] = {"oncogene_like": int(len(onc)), "tsg_like": int(len(tsg)),
                  "actionable_now": int(len(act)),
                  "actionable_genes": act[["gene", "tier", "matched_approved_drugs"]]
                  .to_dict("records")}
out["top20"] = val[cols].head(20).to_dict("records")
json.dump(out, open("results/validated_vs_failed.json", "w"), indent=1, default=float)
print("\nDone.")
