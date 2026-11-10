#!/usr/bin/env python3
"""
Summary statistics for every tested risk gene.

Two tables, both covering all 271 risk genes measured in the GeoMx data.

  risk_gene_summary_statistics.csv   one row per gene. For each stage transition in the
      epithelium: log2 fold change, the NOMINAL p-value from the authors' linear mixed
      model, the risk-restricted FDR (BH over the tested risk genes) and the published
      genome-wide FDR. Then the leading comparison - which transition and which GeoMx
      panel gives that gene its smallest FDR - so every claim about a gene can be traced
      to the comparison it came from. Stroma is summarised as its best nominal p, as a
      specificity check. CosMx columns give the single-cell carcinoma-versus-normal-
      epithelium difference and its paired nominal p where the gene is on the panel.

  risk_gene_all_comparisons.csv      the long form: one row per gene x panel x
      compartment x transition, 3,744 rows, with log2FC, nominal p, both FDRs and the
      direction verdict. This is the table to query for any comparison not summarised
      above, including the three invasion-front contrasts.

Nothing here is filtered. Genes that fail are in the table with their statistics, which
is the point: the reader can see what the evidence was for every gene, not only the
ones that passed.
"""
import pandas as pd, numpy as np, json

MAIN = ["NOR->ADE", "ADE->CORE", "NOR->CORE"]
INV = ["ADE->INV", "NOR->INV", "CORE->INV"]
SHORT = {"NOR->ADE": "N_A", "ADE->CORE": "A_C", "NOR->CORE": "N_C",
         "ADE->INV": "A_I", "NOR->INV": "N_I", "CORE->INV": "C_I"}
PRETTY = {"NOR->ADE": "normal -> adenoma", "ADE->CORE": "adenoma -> carcinoma",
          "NOR->CORE": "normal -> carcinoma", "ADE->INV": "adenoma -> invasive front",
          "NOR->INV": "normal -> invasive front", "CORE->INV": "core -> invasive front"}
PANEL = {"exp2": "WTA (18,441 genes, 9 patients)", "exp1": "CTA (1,781 genes, 10 patients)"}

de = pd.read_csv("results/spatial_risk_gene_DE.csv")
val = pd.read_csv("results/validated_risk_genes_all.csv")
risk = pd.read_csv("results/risk_genes_harmonized.csv")
try:
    cos = pd.read_csv("results/cosmx_crc_vs_normal_epi.csv")
except FileNotFoundError:
    cos = pd.DataFrame(columns=["gene", "delta", "p", "n_patients", "concordant"])

# ---------------------------------------------------------------- long table
long = de.rename(columns={"p": "p_nominal", "FDR": "FDR_genomewide",
                          "FDR_risk": "FDR_risk_restricted",
                          "compartment_label": "compartment_name"}).copy()
long["transition_name"] = long["transition"].map(PRETTY)
long["panel"] = long["experiment"].map(PANEL)
long = long[["gene", "twas_direction", "cell_types", "tissues", "panel", "experiment",
             "compartment", "compartment_name", "transition", "transition_name",
             "log2FC", "p_nominal", "FDR_risk_restricted", "FDR_genomewide",
             "concordant", "sig", "verdict"]].sort_values(
                 ["gene", "compartment", "transition", "p_nominal"])
long.to_csv("results/risk_gene_all_comparisons.csv", index=False)
print(f"long table: {len(long):,} comparisons, {long.gene.nunique()} genes")

# ---------------------------------------------------------------- wide summary
epi = de[de["compartment"] == "PanCK+"].copy()
str_ = de[de["compartment"] == "PanCK-"].copy()
rows = []
for gene, g in epi.groupby("gene"):
    r = {"gene": gene}
    tw = risk[risk["gene"] == gene]
    r["twas_direction"] = g["twas_direction"].iloc[0]
    r["twas_tissues"] = g["tissues"].iloc[0]
    r["twas_cell_types"] = g["cell_types"].iloc[0]
    r["twas_n_celltype_assoc"] = int(g["n_ct"].iloc[0])
    r["twas_min_p"] = float(g["min_p"].iloc[0])
    r["twas_max_abs_z"] = float(tw["abs_z"].max()) if len(tw) else np.nan
    r["geomx_panels"] = "+".join(sorted(g["experiment"].unique()))

    # every main transition, best (most significant) panel for that transition
    for t in MAIN + INV:
        s = g[g["transition"] == t].sort_values("p_nominal" if "p_nominal" in g else "p")
        s = g[g["transition"] == t].sort_values("p")
        k = SHORT[t]
        if len(s):
            b = s.iloc[0]
            r[f"log2FC_{k}"] = float(b["log2FC"])
            r[f"P_{k}"] = float(b["p"])
            r[f"FDR_{k}"] = float(b["FDR_risk"])
            r[f"FDRgw_{k}"] = float(b["FDR"])
            r[f"panel_{k}"] = b["experiment"]
        else:
            for c in (f"log2FC_{k}", f"P_{k}", f"FDR_{k}", f"FDRgw_{k}"):
                r[c] = np.nan
            r[f"panel_{k}"] = ""

    # the leading comparison over the three transitions of the sequence
    m = g[g["transition"].isin(MAIN)].sort_values("FDR_risk")
    if len(m):
        b = m.iloc[0]
        r["lead_comparison"] = PRETTY[b["transition"]]
        r["lead_compartment"] = "epithelium"
        r["lead_panel"] = PANEL[b["experiment"]]
        r["lead_log2FC"] = float(b["log2FC"])
        r["lead_P_nominal"] = float(b["p"])
        r["lead_FDR"] = float(b["FDR_risk"])
        r["lead_FDR_genomewide"] = float(b["FDR"])
        r["lead_direction_verdict"] = b["verdict"]
    r["n_transitions_sig"] = int(m["sig"].sum()) if len(m) else 0
    r["n_transitions_sig_concordant"] = int((m["sig"] & m["concordant"]).sum()) if len(m) else 0
    r["n_transitions_sig_opposing"] = int((m["sig"] & ~m["concordant"]).sum()) if len(m) else 0

    # stroma, as a specificity check
    sg = str_[(str_["gene"] == gene) & (str_["transition"].isin(MAIN))]
    r["stroma_best_P"] = float(sg["p"].min()) if len(sg) else np.nan
    r["stroma_best_FDR"] = float(sg["FDR_risk"].min()) if len(sg) else np.nan
    r["stroma_n_sig"] = int(sg["sig"].sum()) if len(sg) else 0
    rows.append(r)

W = pd.DataFrame(rows)
W = W.merge(val[["gene", "validated", "class", "monotone_across_sequence"]], on="gene", how="left")
W = W.merge(cos.rename(columns={"delta": "cosmx_delta_CRC_minus_normal_epi",
                                "p": "cosmx_P_paired",
                                "n_patients": "cosmx_n_patients",
                                "concordant": "cosmx_concordant"})
            [["gene", "cosmx_delta_CRC_minus_normal_epi", "cosmx_P_paired",
              "cosmx_n_patients", "cosmx_concordant"]], on="gene", how="left")
W["cosmx_on_panel"] = W["cosmx_delta_CRC_minus_normal_epi"].notna()

order = (["gene", "twas_direction", "twas_tissues", "twas_cell_types",
          "twas_n_celltype_assoc", "twas_min_p", "twas_max_abs_z", "geomx_panels"]
         + [f"{p}_{SHORT[t]}" for t in MAIN + INV
            for p in ("log2FC", "P", "FDR", "FDRgw", "panel")]
         + ["lead_comparison", "lead_compartment", "lead_panel", "lead_log2FC",
            "lead_P_nominal", "lead_FDR", "lead_FDR_genomewide", "lead_direction_verdict",
            "n_transitions_sig", "n_transitions_sig_concordant", "n_transitions_sig_opposing",
            "stroma_best_P", "stroma_best_FDR", "stroma_n_sig",
            "validated", "class", "monotone_across_sequence",
            "cosmx_on_panel", "cosmx_delta_CRC_minus_normal_epi", "cosmx_P_paired",
            "cosmx_n_patients", "cosmx_concordant"])
W = W[[c for c in order if c in W.columns]].sort_values("lead_FDR")
W.to_csv("results/risk_gene_summary_statistics.csv", index=False)
print(f"wide table: {len(W)} genes x {W.shape[1]} columns")

# ---------------------------------------------------------------- what it shows
print("\n=== leading comparison, across all 271 genes ===")
print(W["lead_comparison"].value_counts().to_string())
print("\n=== how many genes reach each threshold in the leading comparison ===")
for lab, col, thr in [("nominal P < 0.05", "lead_P_nominal", 0.05),
                      ("nominal P < 0.01", "lead_P_nominal", 0.01),
                      ("FDR (risk-restricted) < 0.05", "lead_FDR", 0.05),
                      ("FDR (genome-wide) < 0.05", "lead_FDR_genomewide", 0.05)]:
    print(f"  {lab:32s} {int((W[col] < thr).sum()):4d} / {len(W)}")
print("\n=== genes with no nominally significant comparison anywhere in the sequence ===")
none_sig = W[W[[f"P_{SHORT[t]}" for t in MAIN]].min(axis=1) >= 0.05]
print(f"  {len(none_sig)} genes")

json.dump({
    "n_genes": int(len(W)),
    "n_comparisons": int(len(long)),
    "lead_comparison_counts": {k: int(v) for k, v in W["lead_comparison"].value_counts().items()},
    "n_lead_nominal_p05": int((W["lead_P_nominal"] < 0.05).sum()),
    "n_lead_nominal_p01": int((W["lead_P_nominal"] < 0.01).sum()),
    "n_lead_fdr_restricted": int((W["lead_FDR"] < 0.05).sum()),
    "n_lead_fdr_genomewide": int((W["lead_FDR_genomewide"] < 0.05).sum()),
    "n_no_nominal_sig": int(len(none_sig)),
    "median_lead_p": float(W["lead_P_nominal"].median()),
}, open("results/summary_statistics_overview.json", "w"), indent=1)
print("\nDone.")
