#!/usr/bin/env python3
"""14_build_supptable_S20_cosmx.py

Supplementary Table S20 - single-cell CosMx support for the scTWAS risk genes,
carcinoma versus normal epithelium, written as a STANDALONE workbook:

    W:/VUMC_projects/scTWAS/ms/SupplementaryTable_S20_CosMx_singlecell.xlsx

It backs Figure 5d (the six per-cell CosMx genes) and merges two tables of the
AI-gene-drug report into one sheet, so the selection can be read against its
denominator:
    S36  results/cosmx_crc_vs_normal_epi.csv  - all 121 genes on the CosMx panel
    S58  results/cosmx_selected_genes.csv     - the 35 nominal genes + selection flags

Every column is either copied from those outputs or recomputed with the same rule
as the AI-gene-drug code, and the recomputed values are asserted equal to the
saved ones before anything is written:
    paired test        code/26_cosmx_gene_tests.py, block B (lines 83-104)
    selection flags    code/48_screened_and_selected.py, lines 259-293
scTWAS statistics come from this manuscript's own 347-gene table.

The main workbook ms/SupplementaryTables.xlsx is NOT modified.
"""
from __future__ import annotations

import os
import numpy as np
import pandas as pd
from openpyxl import Workbook
from openpyxl.styles import Alignment, Border, Font, Side
from openpyxl.utils import get_column_letter

AI = r"L:\AI-gene-drug"
LOCAL = "fig30"            # workspace copies, used when L: is not mounted
MS = r"W:\VUMC_projects\scTWAS\ms"
T347 = (r"W:\VUMC_projects\scTWAS\MultiOmics_DE_support\table"
        r"\scTWAS_347genes_multisource_DE_support.csv")
OUT = os.path.join(MS, "SupplementaryTable_S20_CosMx_singlecell.xlsx")

NORMAL, TUMOUR = "normal/undefined epithelium", "CRC"   # as in scripts 26 and 48
MIN_CELLS = 50                                          # script 26, line 31
ALREADY = ["CSNK2B", "CTNNB1", "UBE2V1", "TRIM28"]      # script 48, line 264
FIG_GENES = ["CEBPB", "GTF2I", "HRAS", "SMARCD1", "WNT4", "ZNRF3"]


def src(rel: str) -> str:
    p = os.path.join(AI, rel)
    return p if os.path.exists(p) else os.path.join(LOCAL, os.path.basename(rel))


def npz(name: str) -> str:
    p = os.path.join(AI, "spatial_in", name)
    return p if os.path.exists(p) else name


# ---- inputs ------------------------------------------------------------------
B = pd.read_csv(src("results/cosmx_crc_vs_normal_epi.csv"))
S58 = pd.read_csv(src("results/cosmx_selected_genes.csv"))
ST6 = pd.read_csv(src("results/cosmx_selected_genes_stats.csv"))
dr = pd.read_csv(src("results/gene_druggability.csv")).set_index("gene")
px = pd.read_csv(src("results/risk_gene_driver_proximity.csv")).set_index("gene")
pb = pd.read_csv(src("spatial_in/cx_pseudobulk.csv.gz"))
ms = pd.read_csv(T347).drop_duplicates("gene_name").set_index("gene_name")
assert len(B) == 121 and len(S58) == 35 and sorted(ST6.gene) == sorted(FIG_GENES)

# ---- patients up, with the pairing of script 26 --------------------------------
pbf = pb[(pb.n_cells >= MIN_CELLS) & pb.cellType.isin([NORMAL, TUMOUR])]
pv = (pbf.pivot_table(index=["gene", "patient_ID"], columns="cellType",
                      values="mean_norm").dropna())
n_up = (pv[TUMOUR] > pv[NORMAL]).groupby(level=0).sum().astype(int)
n_pat = pv.groupby(level=0).size()
assert (n_pat.reindex(B.gene).values == B.n_patients.values).all()
assert (n_up.reindex(ST6.gene).values == ST6.n_up.values).all()

# ---- selection flags, rule of script 48 applied to all 121 ---------------------
drawable = set()
for f in ["cx_plot_cells.npz", "cx_plot_cells3.npz", "cx_plot_cells2.npz"]:
    drawable |= {str(x) for x in np.load(npz(f), allow_pickle=True)["genes"]}

T = B.copy()
T["n_up"] = T.gene.map(n_up)
T["is_driver"] = T.gene.map(lambda g: bool(px.loc[g, "is_driver"]) if g in px.index else False)
T["driver_contact"] = T.gene.map(
    lambda g: float(px.loc[g, "dist_to_nearest_driver"]) == 1.0 if g in px.index else False)
T["crc_related"] = T.is_driver | T.driver_contact
T["tier"] = T.gene.map(lambda g: str(dr.loc[g, "tier"]) if g in dr.index else "")
T["druggable"] = T.tier.str.startswith(("T1", "T2"))
T["drawable"] = T.gene.isin(drawable)
T["already"] = T.gene.isin(ALREADY)
T["nominal"] = T.p < 0.05
T["eligible"] = T.nominal & (T.crc_related | T.druggable) & T.drawable & ~T.already
T["shown"] = T.gene.isin(FIG_GENES)

# the recomputed flags must reproduce the saved selection table row for row
chk, ref = T[T.nominal].set_index("gene").sort_index(), S58.set_index("gene").sort_index()
assert list(chk.index) == list(ref.index)
for mine, theirs in [("is_driver", "is_driver"), ("driver_contact", "driver_contact"),
                     ("crc_related", "crc_related"), ("tier", "tier"),
                     ("druggable", "druggable"), ("drawable", "drawable"),
                     ("already", "already_in_fig28"), ("eligible", "eligible")]:
    a = chk[mine].astype(str)
    b = ref[theirs].astype(str).replace("nan", "")
    assert (a == b).all(), mine
# the six drawn genes are the six eligible genes with the largest |delta|
top6 = T[T.eligible].assign(a=lambda d: d.delta.abs()).nlargest(6, "a").gene
assert sorted(top6) == sorted(FIG_GENES)

# ---- scTWAS statistics from the manuscript's 347-gene table --------------------
assert set(T.gene) <= set(ms.index)
for c in ["rep_model", "rep_Z", "rep_P", "n_models", "direction_status"]:
    T[c] = T.gene.map(ms[c])

# ---- assemble ------------------------------------------------------------------
yn = lambda s: s.map({True: "yes", False: "no"})
out = pd.DataFrame({
    "Gene": T.gene,
    "scTWAS representative model": T.rep_model,
    "scTWAS Z": T.rep_Z,
    "scTWAS P": T.rep_P,
    "scTWAS models (n)": T.n_models.astype(int),
    "Direction across models": T.direction_status,
    "TWAS direction used for concordance": T.twas_direction,
    "Patients (n)": T.n_patients.astype(int),
    "Mean expression, normal epithelium": T.mean_normal_epi,
    "Mean expression, carcinoma": T.mean_CRC,
    "Delta (carcinoma - normal)": T.delta,
    "Patients higher in carcinoma (n)": T.n_up.astype(int),
    "P (paired Wilcoxon, nominal)": T.p,
    "FDR (BH, 121 genes)": T.fdr,
    "Direction concordant with TWAS": yn(T.concordant),
    "GeoMx-validated": yn(T.geomx_validated),
    "Nominal P < 0.05": yn(T.nominal),
    "CRC somatic driver": yn(T.is_driver),
    "Direct PPI contact with a driver": yn(T.driver_contact),
    "Druggability tier": T.tier.replace("", "not assessed"),
    "Approved-drug target (tier 1-2)": yn(T.druggable),
    "Per-cell map available": yn(T.drawable),
    "Drawn in an earlier figure (excluded)": yn(T.already),
    "Eligible for Figure 5d": yn(T.eligible),
    "Shown in Figure 5d": yn(T.shown),
})
out = out.sort_values(["P (paired Wilcoxon, nominal)", "Gene"]).reset_index(drop=True)

TITLE = ("Supplementary Table S20. Single-cell CosMx expression of scTWAS risk genes in "
         "carcinoma cells versus normal epithelium, with the rule used to select the genes "
         "shown in Figure 5d.")
n_nom = int(T.nominal.sum()); n_conc_nom = int((T.nominal & T.concordant).sum())
NOTES = [
    "Notes.",
    f"All {len(T)} scTWAS risk genes on the CosMx panel are listed. For each gene, per-patient "
    f"pseudobulk means were computed for carcinoma cells and for normal epithelium, keeping a "
    f"patient x compartment only if it held >= {MIN_CELLS} cells; the {int(T.n_patients.iloc[0])} "
    f"patients with both compartments were compared with a two-sided paired Wilcoxon signed-rank "
    f"test. Delta is the mean within-patient difference, carcinoma minus normal.",
    f"With {int(T.n_patients.iloc[0])} pairs the smallest attainable P is "
    f"{2.0 ** -(int(T.n_patients.iloc[0]) - 1):.4f}, so no gene can pass FDR < 0.05 across "
    f"{len(T)} tests; the smallest FDR is {T.fdr.min():.4f}. P values are nominal and the table "
    f"records direction of change, not genome-wide significance. {n_nom} genes have nominal "
    f"P < 0.05, of which {n_conc_nom} change in the direction the TWAS predicts.",
    "Direction concordant: delta has the sign predicted by the TWAS direction column "
    "(positive for risk-increasing, negative for protective). That direction comes from the "
    "upstream per-gene summary; 'Direction across models' shows whether the gene's scTWAS "
    "models agree in sign. Five genes are discordant across models (BMP4, PREX1, DAB2 among the "
    "nominal genes); for BMP4 the representative-model Z is positive while the direction used "
    "here is protective, so its concordance call depends on the model chosen.",
    "GeoMx-validated: the gene belongs to the set supported by GeoMx DSP compartment-level "
    "differential expression (an independent patient cohort and platform).",
    "Selection for Figure 5d: nominal P < 0.05, AND (CRC somatic driver OR direct "
    "protein-protein contact with one OR target of an approved drug, druggability tier 1-2), "
    "AND a per-cell expression map available, excluding genes already drawn in an earlier "
    f"figure ({', '.join(ALREADY)}). Seven genes are eligible; the six with the largest "
    "|delta| are shown. FAS (delta -0.0092) is eligible but not shown.",
    "WNT4 is shown in Figure 5d although it is direction-discordant: it is protective in the "
    "scTWAS yet higher in carcinoma.",
    "Druggability tiers: T1/T2 approved-drug targets; T3 experimental ligand available; T4 no "
    "known ligand.",
]


def write(path: str) -> None:
    wb = Workbook()
    ws = wb.active
    ws.title = "S20 CosMx single-cell"
    thin = Side(style="thin", color="000000")
    ws.cell(row=1, column=1, value=TITLE).font = Font(bold=True, size=11)
    for j, c in enumerate(out.columns, 1):
        h = ws.cell(row=2, column=j, value=c)
        h.font = Font(bold=True)
        h.alignment = Alignment(wrap_text=True, vertical="bottom")
        h.border = Border(top=thin, bottom=thin)
    fmt = {"scTWAS Z": "0.00", "scTWAS P": "0.00E+00",
           "Mean expression, normal epithelium": "0.0000",
           "Mean expression, carcinoma": "0.0000",
           "Delta (carcinoma - normal)": "0.0000",
           "P (paired Wilcoxon, nominal)": "0.0000", "FDR (BH, 121 genes)": "0.0000"}
    for i, row in enumerate(out.itertuples(index=False), 3):
        for j, (c, v) in enumerate(zip(out.columns, row), 1):
            v = None if (isinstance(v, float) and np.isnan(v)) else v
            cell = ws.cell(row=i, column=j, value=v)
            if c in fmt:
                cell.number_format = fmt[c]
            if c == "Gene":
                cell.font = Font(italic=True)
    last = 2 + len(out)
    for j in range(1, len(out.columns) + 1):
        ws.cell(row=last, column=j).border = Border(bottom=thin)
    for k, note in enumerate(NOTES):
        c = ws.cell(row=last + 2 + k, column=1, value=note)
        c.font = Font(bold=(k == 0), size=9)
    ws.freeze_panes = "B3"
    ws.row_dimensions[2].height = 45
    for j, c in enumerate(out.columns, 1):
        ws.column_dimensions[get_column_letter(j)].width = 11 if j > 1 else 10
    wb.save(path)


if __name__ == "__main__":
    out.to_csv("SupplementaryTable_S20_CosMx_singlecell.csv", index=False)
    write("SupplementaryTable_S20_CosMx_singlecell.xlsx")
    if os.path.isdir(MS):
        write(OUT)
    print(f"rows {len(out)} | nominal {n_nom} (concordant {n_conc_nom}) | "
          f"eligible {int(T.eligible.sum())} | shown {int(T.shown.sum())} | "
          f"min FDR {T.fdr.min():.4f}")
    print("wrote", OUT if os.path.isdir(MS) else "workspace only")
