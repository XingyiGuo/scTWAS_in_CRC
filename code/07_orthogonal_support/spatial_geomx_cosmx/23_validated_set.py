#!/usr/bin/env python3
"""
23 - The spatially validated risk-gene set.

The TWAS gives every risk gene a direction: risk-increasing (higher expression raises
lifetime risk) or protective (higher expression lowers it). The GeoMx stage contrasts
give every measured gene a somatic direction along normal -> adenoma -> carcinoma.
This script keeps only the genes where those two agree, so that each surviving gene
carries one coherent claim:

    risk-increasing AND up in tumour   -> oncogene-like
    protective      AND down in tumour -> tumour-suppressor-like

Inclusion rule (chosen by the user):
  * FDR < 0.05 in at least one of NOR->ADE, ADE->CORE, NOR->CORE, and
  * that significant change is in the TWAS-predicted direction, and
  * no OTHER of those three transitions is significant in the opposite direction.

The compartment used is the epithelium (PanCK+). Every TWAS cell type in this study is
an epithelial lineage (ABS, ASC, CT, EE, GOB, SSC, STM, TAC, TUF), so epithelium is the
matched compartment and the stroma serves as a specificity control rather than a second
test. Both GeoMx experiments are used; a gene measured on both panels must not be
contradicted on either.

Outputs
  results/validated_risk_genes.csv        one row per validated gene
  results/validated_by_celltype.csv       validation rate per TWAS cell type
  results/validated_summary.json          counts and tests for the report
"""
import pandas as pd, numpy as np, json
from scipy.stats import fisher_exact, binomtest, mannwhitneyu

np.random.seed(23)

MAIN = ["NOR->ADE", "ADE->CORE", "NOR->CORE"]          # the sequence the user asked about
INV = ["ADE->INV", "NOR->INV", "CORE->INV"]            # reported, not part of the rule
PRETTY = {"NOR->ADE": "Normal -> adenoma", "ADE->CORE": "Adenoma -> carcinoma",
          "NOR->CORE": "Normal -> carcinoma", "ADE->INV": "Adenoma -> invasive front",
          "NOR->INV": "Normal -> invasive front", "CORE->INV": "Core -> invasive front"}

de = pd.read_csv("results/spatial_risk_gene_DE.csv")
risk = pd.read_csv("results/risk_genes_harmonized.csv")

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
    m = g[g["transition"].isin(MAIN)]
    if not len(m):
        continue
    sig = m[m["sig"]]
    sup = sig[sig["concordant"]]
    opp = sig[~sig["concordant"]]
    twas_dir = g["twas_direction"].iloc[0]
    passes = len(sup) > 0 and len(opp) == 0

    # the transition that carries the claim: strongest supporting one
    lead = sup.sort_values("FDR").iloc[0] if len(sup) else None
    inv = g[g["transition"].isin(INV)]
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
mm = epi[epi["transition"].isin(MAIN)].copy()
# one row per gene x transition: the most significant experiment, and "significant"
# meaning significant on any panel that measured it - exactly what the rule above used
best = mm.sort_values("FDR").groupby(["gene", "transition"], as_index=False).first()
anysig = mm.groupby(["gene", "transition"])["sig"].any().rename("anysig").reset_index()
best = best.drop(columns=["sig"]).merge(anysig, on=["gene", "transition"])
genes = np.array(sorted(best["gene"].unique()))
S = np.zeros((len(genes), len(MAIN)), dtype=bool)
L = np.zeros((len(genes), len(MAIN)), dtype=float)
gpos = {g: k for k, g in enumerate(genes)}
for j, t in enumerate(MAIN):
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
    "lead_transition_counts": {PRETTY[k]: int(v) for k, v in
                               val["lead_transition"].value_counts().items()},
    "by_celltype": CT.to_dict("records"),
    "validated_genes": val.sort_values("lead_FDR")[
        ["gene", "class", "lead_transition", "lead_log2FC", "lead_FDR", "cell_types",
         "tissues", "monotone_across_sequence", "twas_min_p"]].to_dict("records"),
}, open("results/validated_summary.json", "w"), indent=1, default=float)
print("\nDone.")
