#!/usr/bin/env python3
"""Draw the Figure 5d CosMx cell composition map from the archived coordinates."""

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd


COLORS = {"Carcinoma": "#e34948", "Normal epithelium": "#1baf7a", "Oncofetal": "#eda100",
          "Fibroblast / stroma": "#8a8880", "Immune": "#9ec5f4", "Other": "#d8d7d2"}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    args = parser.parse_args()
    columns = ["x", "y", "group"]
    if args.input.suffix.lower() == ".xlsx":
        cells = pd.read_excel(args.input, sheet_name="Fig5d CosMx Composition", usecols=columns)
    else:
        cells = pd.read_csv(args.input, usecols=columns)
    if not set(cells.group.dropna()) <= COLORS.keys():
        raise ValueError("Unexpected CosMx cell composition groups")
    fig, ax = plt.subplots(figsize=(3.4, 4.2), layout="constrained")
    handles = {}
    for group, color in reversed(COLORS.items()):
        subset = cells.loc[cells.group.eq(group)]
        handles[group] = ax.scatter(subset.x, subset.y, s=0.23, c=color, lw=0, rasterized=True)
    ax.set_aspect("equal")
    ax.set_axis_off()
    ax.set_title("Cell composition", loc="left", fontsize=9, weight="bold")
    ax.legend([handles[group] for group in COLORS], COLORS, loc="lower left",
              bbox_to_anchor=(0, -0.12), frameon=False, markerscale=8, fontsize=6)
    args.out_dir.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.out_dir / "Figure5d_CosMx_composition.pdf")
    fig.savefig(args.out_dir / "Figure5d_CosMx_composition.png", dpi=400)


if __name__ == "__main__":
    main()
