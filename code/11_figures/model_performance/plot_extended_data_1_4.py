#!/usr/bin/env python3
"""Redraw Extended Data Figures 1 and 4 from staged coordinates and model scores."""

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


COLORS = dict(zip(["STM", "TAC", "ABS", "CT", "GOB", "EE", "TUF", "ASC", "SSC"],
                  ["#00b6eb", "#a58aff", "#f8766d", "#c49a00", "#00c094",
                   "#53b400", "#fb61d7", "#24304c", "#53542b"], strict=True))
TITLES = {"NL": "Normal colon", "ADSER": "Colon polyps"}


def panel(ax, letter):
    ax.text(-0.1, 1.03, letter, transform=ax.transAxes, weight="bold", va="top")
    ax.spines[["top", "right"]].set_visible(False)


def plot_umap(ax, path, letter):
    points = pd.read_csv(path)
    for cell_type, rows in points.groupby("cell_type", sort=False):
        ax.scatter(rows.x, rows.y, s=0.6, lw=0, color=COLORS.get(cell_type, "#aaaaaa"),
                   alpha=0.65, rasterized=True)
    groups = "display_label" if "display_label" in points else "cell_type"
    centers = points.groupby(groups)[["x", "y"]].median().sort_values("y")
    if groups == "cell_type":
        for label, center in centers.iterrows():
            ax.annotate(label, xy=center, xytext=(4, 2), textcoords="offset points", fontsize=7)
    else:
        # Place fine annotations on either side of the UMAP with ordered leader lines.
        for side in (False, True):
            labels = centers.loc[centers.x.gt(centers.x.median()).eq(side)]
            x = points.x.max() if side else points.x.min()
            y = np.linspace(points.y.quantile(0.04), points.y.quantile(0.96), len(labels))
            for (label, center), label_y in zip(labels.iterrows(), y, strict=True):
                ax.annotate(label, xy=center, xytext=(x, label_y), fontsize=4.5,
                            ha="right" if side else "left", va="center",
                            arrowprops={"arrowstyle": "-", "lw": 0.45})
    ax.set(xlabel="UMAP 1", ylabel="UMAP 2")
    ax.set_aspect("equal", adjustable="datalim")
    panel(ax, letter)


def plot_peaks(ax, peaks, group, letter):
    counts = peaks.loc[peaks.Group.eq(group)].set_index("CellType").Num_regions_after_chain_merge
    counts = counts.reindex([ct for ct in COLORS if ct in counts.index]) / 100
    bars = ax.bar(counts.index, counts, color=[COLORS[ct] for ct in counts.index], alpha=0.7)
    ax.bar_label(bars, labels=[f"{round(value):,}" for value in counts], padding=2, fontsize=7)
    ax.set(xlabel=TITLES[group], ylabel=r"scATAC peaks ($10^2$)", ylim=(0, counts.max() * 1.17))
    panel(ax, letter)


def extended_data_1(data):
    fig = plt.figure(figsize=(13, 8), layout="constrained")
    grid = fig.add_gridspec(2, 6)
    plot_umap(fig.add_subplot(grid[0, :3]), data / "ExtendedDataFigure1a.csv", "a")
    plot_umap(fig.add_subplot(grid[0, 3:]), data / "ExtendedDataFigure1b.csv", "b")
    peaks = pd.read_csv(data / "ExtendedDataFigure1_peaks.csv")
    plot_peaks(fig.add_subplot(grid[1, :2]), peaks, "NL", "c")
    plot_umap(fig.add_subplot(grid[1, 2:4]), data / "ExtendedDataFigure1d.csv", "d")
    plot_peaks(fig.add_subplot(grid[1, 4:]), peaks, "ADSER", "e")
    return fig


def extended_data_4(data):
    models = pd.read_csv(data / "ExtendedDataFigure4_R2.csv")
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.7), layout="constrained")
    for ax, group, letter in zip(axes, TITLES, "ab", strict=True):
        rows = models.loc[models.group.eq(group)]
        types = [ct for ct in COLORS if ct in set(rows.cell_type)]
        scores = [rows.loc[rows.cell_type.eq(ct), "R2_global"].to_numpy() for ct in types]
        # The normal-colon panel clips the displayed tail after computing the pooled
        # medians. The archived polyp panel shows the distributions within the cap.
        shown = scores if group == "NL" else [values[values <= 0.1] for values in scores]
        violins = ax.violinplot(shown, widths=0.88, showmedians=True, showextrema=False)
        for body, ct in zip(violins["bodies"], types, strict=True):
            body.set(facecolor=COLORS[ct], edgecolor="#515151", alpha=0.72, linewidth=0.7)
        violins["cmedians"].set(color="white", linewidth=1.4)
        labels = [f"{ct}\n{len(values):,}" if group == "NL" else ct
                  for ct, values in zip(types, scores, strict=True)]
        ax.set_xticks(np.arange(1, len(types) + 1), labels)
        ax.set(xlabel=TITLES[group], ylabel=r"Cross-validation $R^2$", ylim=(0.01, 0.103))
        panel(ax, letter)
    return fig


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-dir", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    args = parser.parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)
    for number, plot in ((1, extended_data_1), (4, extended_data_4)):
        fig = plot(args.data_dir)
        for suffix in ("pdf", "png"):
            fig.savefig(args.out_dir / f"ExtendedDataFigure{number}.{suffix}", dpi=300)
        plt.close(fig)
