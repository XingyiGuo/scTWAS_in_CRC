#!/usr/bin/env python
"""Plot Extended Data Figure 9 from the deposited Table 24 or a panel CSV."""
import argparse
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import pandas as pd

WORKBOOK = (Path(__file__).resolve().parents[2] / "manuscript_text"
            / "Li_etal_Supplementary_Tables_11062026.xlsx")
META_GREY = "#666666"

# publication case for compound names; code names stay upper case
CASE = {
    "AEE-788": "AEE-788", "AT-7519": "AT-7519", "BMS-599626": "BMS-599626",
    "CI-1033": "CI-1033", "CP-724714": "CP-724714", "PX-12": "PX-12",
    "ALPHA-LINOLENIC ACID": "\u03b1-linolenic acid",
    "MYCOPHENOLATE MOFETIL": "Mycophenolate mofetil",
    "VALPROIC ACID": "Valproic acid", "CHOLIC ACID": "Cholic acid",
}
# PSMA7;PSMD13;PSMD3 are the three proteasome subunits in the risk set;
# the compounds bind the assembled particle, so the label names the complex.
PROTEASOME = "PSMA7;PSMD13;PSMD3"


def display_name(s):
    u = str(s).strip().upper()
    return CASE.get(u, u.capitalize())


def target_label(t):
    if t == PROTEASOME:
        return "proteasome"                      # complex, not italicised
    return r"$\mathit{%s}$" % str(t).replace(";", r";\!")


@mpl.rc_context({"font.size": 7, "axes.labelsize": 7, "ytick.labelsize": 6, "pdf.fonttype": 42})
def build(input_csv, out_dir):
    if input_csv is None:
        d = pd.read_excel(WORKBOOK, sheet_name="Supplementary Table 24", header=1)
    else:
        d = pd.read_csv(input_csv)
    d["n_signatures"] = d.n_signatures.astype(int)
    d["compound"] = d.drug_name.map(display_name)
    d = d.sort_values("mean_NCS", ascending=False).reset_index(drop=True)
    d["label"] = [f"{c}  ({target_label(t)}, n = {n})"
                  for c, t, n in zip(d.compound, d.risk_gene_targets, d.n_signatures)]
    d["significant"] = (d.FDR_vs_sizematched_null < 0.05) & (d.mean_NCS < 0)
    significance = f"FDR < 0.05 across {len(d)} compounds"

    fig, ax = plt.subplots(figsize=(7.9, 5.7), layout="constrained")
    y = range(len(d))
    lo = d.null_mean - 2 * d.null_sd
    hi = d.null_mean + 2 * d.null_sd
    ax.barh(list(y), (hi - lo).to_numpy(), left=lo.to_numpy(), height=0.55,
            color="#d9dce1", edgecolor="none", zorder=1,
            label="size-matched null (mean \u00b1 2 s.d.)")
    ax.axvline(0, color=META_GREY, lw=0.6, zorder=2)
    ns = ~d.significant
    ax.scatter(d.mean_NCS[ns], [i for i, f in zip(y, ns) if f], s=13,
               color="#4c5561", zorder=3, label="not significant")
    ax.scatter(d.mean_NCS[d.significant], [i for i, f in zip(y, d.significant) if f],
               s=17, color="#b8232f", zorder=4,
               label=significance)
    ax.set_yticks(list(y))
    ax.set_yticklabels(d.label)
    ax.set_ylim(-0.8, len(d) - 0.2)
    ax.set_xlabel("Mean normalized connectivity score across a compound's signatures")
    ax.tick_params(axis="y", length=0)
    ax.margins(x=0.04)
    ax.set_title(f"Matching on PubChem CID, {d.significant.sum()} of {len(d)} risk-gene drugs "
                 "reverse the TWAS profile", loc="left", pad=23)
    h, l = ax.get_legend_handles_labels()
    order = [l.index(x) for x in ["size-matched null (mean \u00b1 2 s.d.)", "not significant",
                                  significance]]
    ax.legend([h[i] for i in order], [l[i] for i in order],
              loc="lower left", bbox_to_anchor=(0.0, 1.005), frameon=False,
              fontsize=6, ncol=3, handlelength=1.2, labelspacing=0.35, borderpad=0.0)
    ax.spines[["top", "right"]].set_visible(False)

    out_dir.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_dir / "ExtendedDataFig9_ilincs_drug_panel.pdf")
    fig.savefig(out_dir / "ExtendedDataFig9_ilincs_drug_panel.png", dpi=400)
    (d.drop(columns=["label"])
       .assign(target_shown=lambda x: x.risk_gene_targets.where(
           x.risk_gene_targets != PROTEASOME, "proteasome"))
       .sort_values("mean_NCS")
       .to_csv(out_dir / "ExtendedDataFig9_ilincs_drug_panel_sourcedata.csv", index=False))
    return fig, ax, d


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, help="Per-compound CSV. Defaults to deposited Table 24.")
    parser.add_argument("--out-dir", type=Path, default=Path("figures"))
    args = parser.parse_args()
    build(args.input, args.out_dir)
