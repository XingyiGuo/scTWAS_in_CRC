#!/usr/bin/env python3
"""Draw the comparison of TF occupancy and co-expression (Extended Data 8)."""

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd


TISSUES = ["Normal", "Adenoma"]
COLORS = ["#3978b5", "#ed6b37"]
MARKERS = ["o", "^"]
KEY_TFS = (
    "SP1 RAD21 POLR2A SMC3 MAX CBX3 NR1H3 CEBPB TFAP4 JUND NIPBL MYC NR1H2 IRF1 "
    "GATA6 ELF1 PPARG YY1 SRF MED12 E2F3 RELA CDX2 PAF1 EGR1 USF1 SOX9 FOXA2"
).split()
DUAL_TFS = (
    "RAD21 POLR2A SP1 MAX SMC3 CBX3 ELF1 PPARG NR1H3 CEBPB NIPBL MYC YY1 IRF1 "
    "ETS2 TFAP4 RELA NR1H2 ETV5 ATF2"
).split()


def source_tables(data_dir):
    correlations = pd.read_csv(data_dir / "corr_long.csv")
    occupancy = pd.read_csv(data_dir / "occ_pairs.csv")
    linked = correlations.loc[correlations.r > 0.5]
    edge = ["tissue", "cell_type", "TF", "gene"]
    dual = occupancy[edge].merge(linked[edge], on=edge)

    # Occupancy sums distinct occupied genes per cell type. The tissue-level
    # co-expression and dual-evidence counts count each gene once per TF.
    pairs = pd.MultiIndex.from_product([TISSUES, sorted(correlations.TF.unique())],
                                       names=["tissue", "TF"])
    occupied = occupancy[edge].drop_duplicates().groupby(["tissue", "TF"]).size()
    coexpressed = linked.groupby(["tissue", "TF"]).gene.nunique()
    scatter = pd.DataFrame({"lead_variants_occupied": occupied,
                            "coexpressed_risk_genes": coexpressed})
    scatter = scatter.reindex(pairs).fillna(0).astype(int).reset_index()
    totals = pd.DataFrame({
        label: frame.groupby("tissue").gene.nunique()
        for label, frame in [
            ("Risk genes in expression models", correlations),
            ("Lead variant occupied by ≥1 CRC TF", occupancy),
            ("Co-expressed with ≥1 CRC TF", linked),
            ("Both occupancy and co-expression", dual),
        ]
    }).reindex(TISSUES).T.rename_axis("evidence").reset_index()
    key = scatter.pivot(index="TF", columns="tissue").reindex(KEY_TFS)
    key.columns = [f"{value}_{tissue}" for value, tissue in key.columns]
    dual_counts = (dual.groupby(["TF", "tissue"]).gene.nunique().unstack(fill_value=0)
                   .reindex(index=DUAL_TFS, columns=TISSUES).fillna(0).astype(int))
    return scatter, totals, key.reset_index(), dual_counts.reset_index()


def plot(tables, out_dir):
    scatter, totals, key, dual = tables
    fig = plt.figure(figsize=(12, 10), layout="constrained")
    grid = fig.add_gridspec(2, 4, width_ratios=[1.15, 1.15, 1, 1])
    ax_a = fig.add_subplot(grid[0, :2])
    ax_b = fig.add_subplot(grid[1, :2])
    ax_c1 = fig.add_subplot(grid[0, 2])
    ax_c2 = fig.add_subplot(grid[0, 3], sharey=ax_c1)
    ax_d = fig.add_subplot(grid[1, 2:])
    labels = "GATA6 MAX RAD21 ELF1 SP1 FOXA2 CBX3 SMC3 POLR2A CEBPB JUND".split()
    label_offsets = {"POLR2A": (-39, 5), "CEBPB": (-21, -10)}

    for i, (tissue, color, marker) in enumerate(zip(TISSUES, COLORS, MARKERS)):
        values = scatter.loc[scatter.tissue == tissue]
        ax_a.scatter(values.lead_variants_occupied, values.coexpressed_risk_genes,
                     color=color, marker=marker, label=tissue, s=24, alpha=.85)
        if tissue == "Normal":
            for row in values.loc[values.TF.isin(labels)].itertuples():
                ax_a.annotate(row.TF, (row.lead_variants_occupied, row.coexpressed_risk_genes),
                              xytext=label_offsets.get(row.TF, (3, 4)),
                              textcoords="offset points", fontsize=7)
        bars = ax_b.barh([j + (i - .5) * .36 for j in range(len(totals))],
                         totals[tissue], .34, color=color, label=tissue)
        ax_b.bar_label(bars, padding=3, fontsize=8)
        ax_c1.scatter(key[f"coexpressed_risk_genes_{tissue}"], key.TF,
                      color=color, marker=marker, s=24)
        ax_c2.scatter(key[f"lead_variants_occupied_{tissue}"], key.TF,
                      color=color, marker=marker, label=tissue, s=24)
        ax_d.scatter(dual[tissue], dual.TF, color=color, marker=marker,
                     label=tissue, s=24)

    ax_a.set(xlabel="Risk genes whose lead variant the TF occupies",
             ylabel="Risk genes co-expressed with TF (r > 0.5)")
    ax_b.set_yticks(range(len(totals)), totals.evidence, fontsize=8)
    ax_b.set(xlabel="Number of risk genes", xlim=(0, 300))
    ax_c1.set(title="Co-expressed risk genes", xlabel="Risk genes")
    ax_c2.set(title="Lead variants occupied", xlabel="Risk genes", xlim=(0, 220))
    ax_c2.tick_params(labelleft=False)
    ax_d.set(xlabel="Risk genes with dual evidence", xlim=(0, 20))
    for axis in (ax_b, ax_c1, ax_d):
        axis.invert_yaxis()
    for axis in fig.axes:
        axis.spines[["top", "right"]].set_visible(False)
        axis.grid(alpha=.2)
        axis.set_axisbelow(True)
    for axis in (ax_a, ax_b, ax_c2, ax_d):
        axis.legend(frameon=False, fontsize=8, loc="lower right")
    for panel, axis, left in zip("abcd", (ax_a, ax_b, ax_c1, ax_d), (-.08, -.08, -.18, -.08)):
        axis.text(left, 1.01, panel, transform=axis.transAxes, weight="bold", fontsize=13)

    out_dir.mkdir(parents=True, exist_ok=True)
    for extension in ("pdf", "png"):
        fig.savefig(out_dir / f"ExtendedDataFigure8.{extension}", dpi=200)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-dir", type=Path, default=Path.cwd())
    parser.add_argument("--out-dir", type=Path, default=Path.cwd())
    args = parser.parse_args()
    plot(source_tables(args.data_dir), args.out_dir)
