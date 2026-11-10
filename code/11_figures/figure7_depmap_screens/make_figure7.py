#!/usr/bin/env python3
"""Draw Figure 7a–c from the archived DepMap and arrayed CRISPR source data."""

import argparse
from itertools import combinations
from pathlib import Path

import pandas as pd
from scipy.stats import mannwhitneyu
from statsmodels.stats.multitest import multipletests

from build_figure7 import GREY, ONC, ONC_L, ONC_M, build_figure7


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-dir", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    args = parser.parse_args()
    data = args.data_dir
    medians = pd.read_csv(data / "DepMap_essentiality_by_geneset_sourcedata.csv")
    order = medians.sort_values("group_order").group.unique().tolist()
    values = {group: medians.loc[medians.group.eq(group), "median_gene_effect"].to_numpy()
              for group in order}
    palette = dict(zip(order, [GREY, ONC_L, ONC_M, ONC], strict=True))
    pairs = list(combinations(order, 2))
    raw_p = [mannwhitneyu(values[a], values[b], alternative="two-sided").pvalue for a, b in pairs]
    adj_p = multipletests(raw_p, method="holm")[1]
    tests = pd.DataFrame(pairs, columns=["group_A", "group_B"])
    tests["P_two_sided"] = raw_p
    tests["P_two_sided_Holm"] = adj_p
    tests["significant_Holm_0.05"] = adj_p < 0.05
    shown = [(order[0], order[2]), (order[1], order[2]), (order[2], order[3])]

    stats = pd.read_csv(data / "DepMap_36dependency_genes_stats.csv", index_col="symbol")
    deps = stats.loc[stats["class"].eq("Pan-essential")].sort_values("median")
    long = pd.read_csv(data / "DepMap_36dependency_genes_boxplot_sourcedata.csv")
    long = long.loc[long.symbol.isin(deps.index)].rename(
        columns={"symbol": "gene", "depmap_gene_effect": "gene_effect"})
    screen = pd.read_csv(data / "Figure7_panelC_sourcedata.csv")
    top = (screen.drop_duplicates("gene").sort_values("panel_rank").set_index("gene")
           .rename(columns={"scTWAS_rep_Z": "rep_Z"})[["rep_Z"]])
    cells = {(r.gene, r.cell_line): (r.delta_viability_points, r.significant_P_lt_0_05)
             for r in screen.itertuples()}
    lines = ["HCT116", "RKO", "SW480"]
    if len(deps) != 24 or len(top) != 24 or len(cells) != 72:
        raise ValueError("Expected 24 DepMap dependencies and 72 arrayed screen results")

    args.out_dir.mkdir(parents=True, exist_ok=True)
    tests.to_csv(args.out_dir / "Figure7a_pairwise_tests.csv", index=False)
    build_figure7(long, deps.index.tolist(), values, order, palette, tests, shown, top,
                  cells, lines, args.out_dir / "Figure7.png", args.out_dir / "Figure7.pdf")


if __name__ == "__main__":
    main()
