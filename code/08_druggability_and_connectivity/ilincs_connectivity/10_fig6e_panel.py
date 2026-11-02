#!/usr/bin/env python
"""Figure 6e - iLINCS connectivity of the pre-specified risk-gene drug panel.

Rebuilds panel a of Fig_ilincs_drug_panel.pdf as a stand-alone, assembly-ready
panel in Nature Genetics style: no in-figure title (it moves to the legend), no
panel letter (added at assembly), compound names in publication case, gene
symbols italic, "s.d." rather than "SD".

Panel b of the old figure (per-model NCS for the three credible hits) is not
reproduced here; its numbers are quoted in the Results text instead.

Input : figure/Fig_ilincs_drug_panel_sourcedata.csv  (panel == 'a', 43 rows)
Output: figure/Fig6e_ilincs_drug_panel.pdf / .png
        figure/Fig6e_ilincs_drug_panel_sourcedata.csv
"""
import os
import matplotlib as mpl
import matplotlib.pyplot as plt
import pandas as pd

OUT = r"W:\VUMC_projects\scTWAS\iLINCS_reversal\figure"
SRC = os.path.join(OUT, "Fig_ilincs_drug_panel_sourcedata.csv")

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
    return r"$\mathit{%s}$" % str(t).replace(";", ";\!")


def build():
    d = pd.read_csv(SRC)
    d = d[d.panel == "a"].drop(columns=["panel", "model"]).copy()
    d["n_signatures"] = d.n_signatures.astype(int)
    d["compound"] = d.drug_name.map(display_name)
    d = d.sort_values("mean_NCS", ascending=False).reset_index(drop=True)
    d["label"] = [f"{c}  ({target_label(t)}, n = {n})"
                  for c, t, n in zip(d.compound, d.risk_gene_targets, d.n_signatures)]
    d["significant"] = d.FDR_vs_sizematched_null < 0.05

    apply_figure_style(sizes=(7, 6, 6))
    fig, ax = plt.subplots(figsize=(4.72, 6.10))
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
               label="FDR < 0.05 across the 43 compounds")
    ax.set_yticks(list(y))
    ax.set_yticklabels(d.label)
    ax.set_ylim(-0.8, len(d) - 0.2)
    ax.set_xlabel("Mean normalized connectivity score across a compound's signatures")
    ax.tick_params(axis="y", length=0)
    ax.margins(x=0.04)
    ax.annotate("\u2190 reverses the predicted risk profile", xy=(0.03, -0.075),
                xycoords="axes fraction", ha="left", va="top",
                fontsize=6, color=META_GREY, annotation_clip=False)
    h, l = ax.get_legend_handles_labels()
    order = [l.index(x) for x in ["size-matched null (mean \u00b1 2 s.d.)", "not significant",
                                  "FDR < 0.05 across the 43 compounds"]]
    ax.legend([h[i] for i in order], [l[i] for i in order],
              loc="lower left", bbox_to_anchor=(0.0, 1.005), frameon=False,
              fontsize=6, handlelength=1.2, labelspacing=0.35, borderpad=0.0)
    set_frame(ax, "open")
    ax.tick_params(axis="y", length=0)
    fig.tight_layout()

    os.makedirs(OUT, exist_ok=True)
    fig.savefig(os.path.join(OUT, "Fig6e_ilincs_drug_panel.pdf"))
    fig.savefig(os.path.join(OUT, "Fig6e_ilincs_drug_panel.png"), dpi=400)
    (d.drop(columns=["label"])
       .assign(target_shown=lambda x: x.risk_gene_targets.where(
           x.risk_gene_targets != PROTEASOME, "proteasome"))
       .sort_values("mean_NCS")
       .to_csv(os.path.join(OUT, "Fig6e_ilincs_drug_panel_sourcedata.csv"), index=False))
    return fig, ax, d


if __name__ == "__main__":
    fig, ax, d = build()
