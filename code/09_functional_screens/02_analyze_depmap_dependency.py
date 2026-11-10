#!/usr/bin/env python3
"""Recompute Figure 7a–b from the archived colorectal DepMap score exports."""

import argparse
from itertools import combinations
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import kruskal, mannwhitneyu
from statsmodels.stats.multitest import multipletests


def summarize_scores(scores):
    scores = scores.dropna(subset=["depmap_gene_effect"]).drop_duplicates(["symbol", "depmap_id"])
    by_gene = scores.groupby("symbol").depmap_gene_effect
    summary = by_gene.agg(n_lines="size", median="median", mean="mean",
                          q25=lambda s: s.quantile(0.25), q75=lambda s: s.quantile(0.75),
                          min="min", max="max")
    summary["frac_lt_0.5"] = by_gene.apply(lambda s: s.lt(-0.5).mean())
    return scores, summary.loc[summary.n_lines.gt(30)]


def analyze_genesets(root):
    membership = pd.read_csv(root / "DepMap_essentiality_by_geneset_sourcedata.csv",
                             usecols=["group", "group_order", "gene"])
    scores = pd.read_csv(root / "depmap_intestine_all_geneset_genes_percellline.csv")
    _, summary = summarize_scores(scores.rename(columns={"gene_effect": "depmap_gene_effect"}))
    columns = {"median": "median_gene_effect", "frac_lt_0.5": "frac_lines_lt_0_5"}
    summary = summary.drop(columns="mean").rename(columns=columns)
    groups = membership.join(summary, on="gene", validate="many_to_one")
    if groups.n_lines.isna().any():
        raise ValueError("A gene in the archived Figure 7a groups lacks more than 30 scores")
    groups = groups.sort_values(["group_order", "median_gene_effect", "gene"])
    values = {label: group.median_gene_effect.to_numpy()
              for label, group in groups.groupby("group", sort=False)}
    omnibus = kruskal(*values.values())
    pairwise = []
    for a, b in combinations(values, 2):
        x, y = values[a], values[b]
        test = mannwhitneyu(x, y, alternative="two-sided")
        pairwise.append({"group_A": a, "group_B": b, "U": test.statistic,
                         "P_two_sided": test.pvalue,
                         "median_difference": np.median(x) - np.median(y),
                         "hodges_lehmann": np.median(x[:, None] - y)})
    tests = pd.DataFrame(pairwise)
    tests["P_two_sided_Holm"] = multipletests(tests.P_two_sided, method="holm")[1]
    tests["significant_Holm_0.05"] = tests.P_two_sided_Holm.lt(0.05)
    return groups, tests, omnibus


def analyze_risk_genes(root):
    raw = pd.read_csv(root / "scTWAS_347genes_DepMap_largeintestine_percellline.csv")
    scores, summary = summarize_scores(raw)
    dependencies = summary.loc[summary["median"].lt(-0.5)].sort_values("median").copy()
    dependencies["class"] = np.where(dependencies["frac_lt_0.5"].ge(0.90), "Pan-essential", "Other")
    values = scores.loc[scores.symbol.isin(dependencies.index),
                        ["symbol", "depmap_id", "cell_line", "depmap_gene_effect"]]
    values = values.join(dependencies[["class", "median", "n_lines"]], on="symbol")
    return dependencies, values


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--experiments", type=Path, required=True, help="Directory with the archived DepMap exports")
    parser.add_argument("--out-dir", type=Path, default=Path("results"))
    args = parser.parse_args()

    groups, tests, omnibus = analyze_genesets(args.experiments)
    dependencies, values = analyze_risk_genes(args.experiments)
    args.out_dir.mkdir(parents=True, exist_ok=True)
    groups.to_csv(args.out_dir / "DepMap_essentiality_by_geneset_sourcedata.csv", index=False)
    tests.to_csv(args.out_dir / "DepMap_geneset_pairwise_tests.csv", index=False)
    dependencies.to_csv(args.out_dir / "DepMap_36dependency_genes_stats.csv")
    values.to_csv(args.out_dir / "DepMap_36dependency_genes_boxplot_sourcedata.csv", index=False)

    print(groups.groupby("group", sort=False).median_gene_effect.agg(["size", "median"]).to_string())
    print(f"Kruskal–Wallis H = {omnibus.statistic:.3f}, P = {omnibus.pvalue:.3g}")
    print(dependencies["class"].value_counts().to_string())


if __name__ == "__main__":
    main()
