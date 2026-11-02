#!/usr/bin/env python3
"""
01 - Parse and harmonize cell-type-specific CRC risk genes from TWAS results.

Inputs
  data/HTAN_NL_WellcomeSanger_NL_7celltypes_peakSNP.EUR_Multi.CauchyMeta.sig.csv   (normal mucosa, 7 cell types)
  data/VUMC_HTAN_ADSER_9CTs_peakSNP.EUR_Multi.sig.csv                              (adenoma, 9 cell types)

Output
  results/risk_genes_harmonized.csv   one row per (tissue, cell type, gene)
"""
import pandas as pd
import numpy as np
import os

os.makedirs("results", exist_ok=True)

NL = "data/HTAN_NL_WellcomeSanger_NL_7celltypes_peakSNP.EUR_Multi.CauchyMeta.sig.csv"
AD = "data/VUMC_HTAN_ADSER_9CTs_peakSNP.EUR_Multi.sig.csv"

nl = pd.read_csv(NL, encoding="utf-8-sig", low_memory=False)
ad = pd.read_csv(AD, encoding="utf-8-sig", low_memory=False)

print("NL shape", nl.shape)
print("AD shape", ad.shape)


def pick(df, names):
    for n in names:
        if n in df.columns:
            return n
    return None


# ---- Normal mucosa -------------------------------------------------------
# Two discovery datasets (HTAN_NL, WellcomeSanger_NL) meta-analysed with Cauchy.
# Use the Multi-ancestry meta p-value as primary; z-score from the driving dataset.
nl_out = pd.DataFrame({
    "tissue": "Normal",
    "cell_type": nl["CellType"],
    "ensembl": nl["gene"],
    "gene": nl["gene_name"],
    "z_primary": np.where(
        nl["P driven_Multi"].astype(str).str.contains("Wellcome", na=False),
        nl["zscore_WellcomeSanger_NL_Multi"],
        nl["zscore_HTAN_NL_Multi"],
    ),
    "z_eur": np.where(
        nl["P driven_EUR"].astype(str).str.contains("Wellcome", na=False),
        nl["zscore_WellcomeSanger_NL_EUR"],
        nl["zscore_HTAN_NL_EUR"],
    ),
    "p_multi": nl["META_P_Multi"],
    "p_eur": nl["META_P_EUR"],
    "direction_concordance": nl["TWAS_zscore_direction_Multi"],
    "sig_resource": nl["Sig_resource"],
})

# ---- Adenoma -------------------------------------------------------------
ad_out = pd.DataFrame({
    "tissue": "Adenoma",
    "cell_type": ad["CellType"],
    "ensembl": ad["gene"],
    "gene": ad["gene_name"],
    "z_primary": ad["zscore_Multi"],
    "z_eur": ad["zscore_EUR"],
    "p_multi": ad["pvalue_Multi"],
    "p_eur": ad["pvalue_EUR"],
    "direction_concordance": np.nan,
    "sig_resource": ad["Sig_resource"],
})

risk = pd.concat([nl_out, ad_out], ignore_index=True)

# Best available p-value / z per row
risk["p_best"] = risk[["p_multi", "p_eur"]].min(axis=1)
risk["z_best"] = np.where(
    risk["p_multi"].fillna(1) <= risk["p_eur"].fillna(1),
    risk["z_primary"], risk["z_eur"]
)
risk["z_best"] = risk["z_best"].fillna(risk["z_primary"]).fillna(risk["z_eur"])

# Directional label:
#   z > 0  -> genetically higher predicted expression increases CRC risk  => "risk-increasing"
#             therapeutic hypothesis = INHIBIT / knock down
#   z < 0  -> higher predicted expression is protective                   => "protective"
#             therapeutic hypothesis = ACTIVATE / restore
risk["direction"] = np.where(risk["z_best"] > 0, "risk-increasing", "protective")
risk["therapeutic_hypothesis"] = np.where(
    risk["z_best"] > 0, "inhibit/degrade", "activate/restore")
risk["abs_z"] = risk["z_best"].abs()
risk["neglog10p"] = -np.log10(risk["p_best"].clip(lower=1e-320))

risk = risk.dropna(subset=["gene"])
risk["gene"] = risk["gene"].astype(str).str.strip()
risk = risk[risk["gene"] != ""]

risk = risk.sort_values(["tissue", "cell_type", "p_best"]).reset_index(drop=True)
risk.to_csv("results/risk_genes_harmonized.csv", index=False)

# ---- Summaries -----------------------------------------------------------
print("\n=== Rows per tissue ===")
print(risk.groupby("tissue").size())

print("\n=== Cell types ===")
for t, g in risk.groupby("tissue"):
    print(f"\n{t} ({g['cell_type'].nunique()} cell types, "
          f"{g['gene'].nunique()} unique genes)")
    print(g.groupby("cell_type").agg(
        n_genes=("gene", "nunique"),
        n_risk_incr=("direction", lambda s: (s == "risk-increasing").sum()),
        n_protective=("direction", lambda s: (s == "protective").sum()),
    ).to_string())

print("\n=== Direction split overall ===")
print(risk.groupby(["tissue", "direction"]).size())

# Gene-level: how many cell types is each gene significant in?
gene_ct = (risk.groupby(["gene"])
           .agg(n_celltype_tissue=("cell_type", "size"),
                n_tissues=("tissue", "nunique"),
                tissues=("tissue", lambda s: "|".join(sorted(set(s)))),
                cell_types=("cell_type", lambda s: "|".join(sorted(set(s)))),
                max_abs_z=("abs_z", "max"),
                min_p=("p_best", "min"),
                directions=("direction", lambda s: "|".join(sorted(set(s)))))
           .sort_values(["n_celltype_tissue", "max_abs_z"], ascending=False))
gene_ct.to_csv("results/gene_celltype_summary.csv")
print("\n=== Top pleiotropic genes (most cell types) ===")
print(gene_ct.head(25).to_string())

print("\nUnique genes overall:", risk["gene"].nunique())
print("Shared between Normal and Adenoma:",
      (gene_ct["n_tissues"] == 2).sum())
