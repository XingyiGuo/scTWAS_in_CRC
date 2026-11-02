#!/usr/bin/env python
"""Druggability of the 347 scTWAS colorectal risk genes.

Evidence source: W:/VUMC_projects/scTWAS/Drugs/Proteins_Drugs_20250321.xlsx,
sheet 'ChEMBL_DB_TTD' (73,913 gene-drug records pooled from ChEMBL, DrugBank
and TTD).  Each source uses its own development-stage vocabulary; they are
harmonised onto one ordinal ladder (see TIER_NAME / tier_of) and the highest
tier reached by any drug is assigned to the gene.

Outputs (under <ROOT>/Drugs):
  table/scTWAS_347genes_druggability_ChEMBL_DB_TTD.csv   347 rows, one per gene
  table/scTWAS_druggable_gene_drug_pairs.csv             exploded gene x drug records
  table/scTWAS_druggability_tier_summary.csv             tier counts
  figure/Fig_scTWAS_druggability.png / .pdf              two-panel summary
  figure/Fig_scTWAS_druggability_sourcedata.csv
"""
import re
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

ROOT = r"W:\VUMC_projects\scTWAS"
XLSX = rf"{ROOT}\Drugs\Proteins_Drugs_20250321.xlsx"
GENES = rf"{ROOT}\MultiOmics_DE_support\table\scTWAS_347genes_multisource_DE_support.csv"
DEPMAP = rf"{ROOT}\MultiOmics_DE_support\data\scTWAS_347genes_DepMap_largeintestine.csv"
OUT_T = rf"{ROOT}\Drugs\table"
OUT_F = rf"{ROOT}\Drugs\figure"

TIER_NAME = {5: "Approved drug", 4: "Phase 3 / registration", 3: "Phase 2",
             2: "Phase 1 / clinical", 1: "Preclinical / experimental ligand",
             0: "Discontinued or withdrawn only", -1: "No entry in ChEMBL/DrugBank/TTD"}

TTD_MAP = {**{k: 5 for k in ["approved", "approved (orphan drug)", "registered", "phase 4"]},
           **{k: 4 for k in ["phase 3", "phase 2/3", "preregistration", "nda filed",
                             "bla submitted", "application submitted", "approval submitted"]},
           **{k: 3 for k in ["phase 2", "phase 2a", "phase 2b", "phase 1/2",
                             "phase 1b/2a", "phase 1/2a"]},
           **{k: 2 for k in ["phase 1", "phase 1b", "phase 0", "clinical trial"]},
           **{k: 1 for k in ["preclinical", "investigative", "patented"]}}

# action_type values that denote a pharmacological effect on the target, as opposed
# to a binding/cofactor annotation (zinc, NADH, pyridoxal phosphate and the like).
PHARM = {"INHIBITOR", "ANTAGONIST", "AGONIST", "MODULATOR", "ACTIVATOR", "BLOCKER",
         "POTENTIATOR", "INDUCER", "PARTIAL AGONIST", "INVERSE AGONIST",
         "NEGATIVE ALLOSTERIC MODULATOR", "POSITIVE ALLOSTERIC MODULATOR", "ANTIBODY",
         "OPENER", "DISRUPTING AGENT", "DEGRADER", "CROSS-LINKING AGENT",
         "STABILISER", "SEQUESTERING AGENT"}


def tier_of(src, mp):
    """Map one source's development-stage string onto the shared ordinal ladder."""
    s = "" if pd.isna(mp) else str(mp).strip().strip("'\"").lower()
    if src == "ChEMBL":                      # numeric ChEMBL max_phase
        return {"4": 5, "3": 4, "2": 3, "1": 2, "0": 1}.get(s, 1)
    if src == "TTD":
        if s.startswith(("discontinued", "terminated", "withdrawn")):
            return 0
        return TTD_MAP.get(s, 1)
    toks = {t.strip().strip("'\" ") for t in re.split(r",", s)}   # DrugBank group list
    if "approved" in toks:
        return 5
    if "investigational" in toks:
        return 2
    if toks & {"experimental", "nutraceutical", "vet_approved", "illicit"}:
        return 1
    if "withdrawn" in toks:
        return 0
    return 1


# ---------------------------------------------------------------- load + match
dd = pd.read_excel(XLSX, sheet_name="ChEMBL_DB_TTD")
ms = pd.read_csv(GENES)
dep = pd.read_csv(DEPMAP)
genes347 = ms["gene_name"].dropna().unique().tolist()
assert len(genes347) == 347

# 'Genes' holds protein complexes / target families as ';'- or '|'-joined symbols
dd["gene_list"] = dd["Genes"].astype(str).str.split(r"\s*[;|]\s*")
dl = dd.explode("gene_list")
dl["gene_sym"] = dl["gene_list"].str.strip()
dl = dl[dl.gene_sym.ne("") & dl.gene_sym.ne("nan")].copy()
dl["tier"] = [tier_of(s, m) for s, m in zip(dl["Source"], dl["Max Phase"])]

sub = dl[dl.gene_sym.isin(genes347)].copy()
sub["drug"] = sub["Drug Brand Name"].astype(str).str.strip().str.upper().replace("NAN", np.nan)
sub["act"] = (sub["action_type"].astype(str).str.strip().str.strip("'\"").str.upper()
              .replace("NAN", np.nan))
sub["pharmacological_action"] = sub["act"].isin(PHARM)
sub["tier_label"] = sub["tier"].map(TIER_NAME)


# ----------------------------------------------------------------- per gene
def agg(g):
    nd = lambda m: g.loc[m, "drug"].dropna().nunique()
    top = g[g.tier == g.tier.max()]
    return pd.Series({
        "max_tier": int(g.tier.max()),
        "n_records": len(g),
        "n_drugs": g["drug"].dropna().nunique(),
        "n_approved_drugs": nd(g.tier == 5),
        "n_clinical_drugs": nd(g.tier.isin([2, 3, 4])),
        "n_preclinical_drugs": nd(g.tier == 1),
        "n_discontinued_drugs": nd(g.tier == 0),
        "n_NCI_cancer_drugs": nd(g["NCI Cancer Drugs"] == True),  # noqa: E712
        "n_approved_pharmacological": nd((g.tier == 5) & g.pharmacological_action),
        "sources": "|".join(sorted(g["Source"].unique())),
        "action_types": "|".join(sorted(g["act"].dropna().unique()))[:200],
        "top_tier_drug_examples": "; ".join(sorted(top["drug"].dropna().unique())[:6]),
        "NCI_cancer_drug_names": "; ".join(
            sorted(g.loc[g["NCI Cancer Drugs"] == True, "drug"].dropna().unique())[:8]),  # noqa: E712
    })


per_gene = (sub.groupby("gene_sym").apply(agg, include_groups=False)
            .reset_index().rename(columns={"gene_sym": "gene_name"}))

base = ms[["gene", "gene_name", "rep_model", "rep_Dataset", "rep_CellType",
           "rep_Z", "rep_P", "rep_direction", "n_models"]].merge(
    dep[["gene_name", "depmap_n_lines", "depmap_median", "frac_lines_lt_0.5",
         "depmap_dependency", "depmap_selectivity"]], on="gene_name", how="left")
out = base.merge(per_gene, on="gene_name", how="left")
out["max_tier"] = out["max_tier"].fillna(-1).astype(int)
for c in ["n_records", "n_drugs", "n_approved_drugs", "n_clinical_drugs",
          "n_preclinical_drugs", "n_discontinued_drugs", "n_NCI_cancer_drugs",
          "n_approved_pharmacological"]:
    out[c] = out[c].fillna(0).astype(int)
for c in ["sources", "action_types", "top_tier_drug_examples", "NCI_cancer_drug_names"]:
    out[c] = out[c].fillna("")
out["druggability_tier"] = out["max_tier"].map(TIER_NAME)
out["in_ChEMBL_DB_TTD"] = out["max_tier"].ge(0)
out = out.sort_values(["max_tier", "n_approved_pharmacological", "n_drugs"],
                      ascending=False).reset_index(drop=True)
assert len(out) == 347 and out["gene_name"].is_unique

order = [5, 4, 3, 2, 1, 0, -1]
summary = (pd.DataFrame({"tier": order, "tier_label": [TIER_NAME[t] for t in order]})
           .merge(out.groupby("max_tier").size().rename("n_genes").reset_index()
                  .rename(columns={"max_tier": "tier"}), on="tier", how="left"))
summary["n_genes"] = summary["n_genes"].fillna(0).astype(int)
summary["pct_of_347"] = (100 * summary["n_genes"] / 347).round(1)
assert summary["n_genes"].sum() == 347

out.to_csv(rf"{OUT_T}\scTWAS_347genes_druggability_ChEMBL_DB_TTD.csv", index=False)
sub.drop(columns=["gene_list"]).rename(columns={"gene_sym": "gene_name"}).sort_values(
    ["gene_name", "tier"], ascending=[True, False]).to_csv(
    rf"{OUT_T}\scTWAS_druggable_gene_drug_pairs.csv", index=False)
summary.to_csv(rf"{OUT_T}\scTWAS_druggability_tier_summary.csv", index=False)

# ------------------------------------------------------------------- figure
mpl.rcParams.update({"font.size": 8, "axes.titlesize": 8, "axes.labelsize": 8,
                     "xtick.labelsize": 6, "ytick.labelsize": 6, "legend.fontsize": 7,
                     "axes.spines.top": False, "axes.spines.right": False,
                     "axes.linewidth": 0.6, "xtick.major.width": 0.6,
                     "ytick.major.width": 0.6, "figure.dpi": 150,
                     "font.family": "sans-serif", "pdf.fonttype": 42})

DEP_COL = {"pan-essential (>=90% of lines)": "#7a3b8f",
           "selective dependency": "#1b7f5e",
           "not a dependency": "#9aa4ad"}
DEP_LAB = {"pan-essential (>=90% of lines)": "pan-essential in large intestine",
           "selective dependency": "selective dependency",
           "not a dependency": "not a DepMap dependency"}
UNSCORED = "#d9dde1"

fig = plt.figure(figsize=(7.2, 6.6))
gs = fig.add_gridspec(2, 1, height_ratios=[1.0, 1.75], hspace=0.30,
                      left=0.30, right=0.97, top=0.93, bottom=0.10)

# Panel A -- tier distribution over all 347 genes
axA = fig.add_subplot(gs[0])
sA = summary.iloc[::-1].reset_index(drop=True)          # low tier at bottom
ramp = plt.get_cmap("YlGnBu")
cols = [UNSCORED if t == -1 else ramp(0.25 + 0.13 * t) for t in sA["tier"]]
axA.barh(sA["tier_label"], sA["n_genes"], color=cols, edgecolor="white", linewidth=0.5)
for y, n in enumerate(sA["n_genes"]):
    if n == 0:
        axA.plot([0.8], [y], marker="|", color="#444", ms=5, mew=0.8)   # zero stub
    axA.text(n + 5, y, f"{n}" + (f"  ({100*n/347:.0f}%)" if n >= 20 else ""),
             va="center", ha="left", fontsize=7)
axA.set_xlim(0, 330)
axA.set_xlabel("scTWAS risk genes (of 347)")
axA.set_title("Most colorectal risk genes have no ligand on record; 44 have an approved drug",
              loc="left")

# Panel B -- genes with an approved agent carrying a pharmacological action
axB = fig.add_subplot(gs[1])
pb = out[out["n_approved_pharmacological"] > 0].sort_values(
    ["n_approved_pharmacological", "gene_name"], ascending=[True, False])
colsB = [DEP_COL.get(s, UNSCORED) for s in pb["depmap_selectivity"]]
ypos = np.arange(len(pb))
axB.barh(ypos, pb["n_approved_pharmacological"], color=colsB,
         edgecolor="white", linewidth=0.5)
axB.set_yticks(ypos)
axB.set_yticklabels([f"$\\it{{{g}}}$" + ("  \u25b2" if n else "")
                     for g, n in zip(pb["gene_name"], pb["n_NCI_cancer_drugs"])])
for y, (n, ex) in enumerate(zip(pb["n_approved_pharmacological"],
                                pb["top_tier_drug_examples"])):
    first = ex.split(";")[0].strip().title()
    axB.text(n + 0.25, y, f"{n}   {first}", va="center", ha="left", fontsize=6,
             color="#33393f")
axB.set_xlim(0, 27)
axB.set_ylim(-0.8, len(pb) - 0.2)
axB.set_xlabel("Approved agents annotated with a pharmacological action on the target")
n_onc_b = int((pb["n_NCI_cancer_drugs"] > 0).sum())
axB.set_title(f"25 genes are hit by an approved agent; {n_onc_b} of them by a listed "
              "oncology drug (\u25b2)", loc="left")
handles = [Patch(facecolor=DEP_COL[k], label=DEP_LAB[k]) for k in DEP_COL]
handles.append(Patch(facecolor=UNSCORED, label="not scored in DepMap"))
axB.legend(handles=handles, frameon=False, loc="lower right", handlelength=1.1,
           handleheight=0.9, borderpad=0.2, labelspacing=0.35)

for ax, L in ((axA, "a"), (axB, "b")):
    ax.text(-0.30, 1.04, L, transform=ax.transAxes, fontweight="bold", fontsize=11,
            va="bottom", ha="left")

for path, kw in ((rf"{OUT_F}\Fig_scTWAS_druggability.png", dict(format="png", dpi=300)),
                 (rf"{OUT_F}\Fig_scTWAS_druggability.pdf", dict(format="pdf"))):
    with open(path, "wb") as fh:          # network mount rejects PIL 'w+b'
        fig.savefig(fh, **kw)
fig.savefig("Fig_scTWAS_druggability.png", format="png", dpi=300)   # local copy

src = pd.concat([
    summary.assign(panel="a")[["panel", "tier_label", "n_genes", "pct_of_347"]],
    pb.assign(panel="b")[["panel", "gene_name", "n_approved_pharmacological",
                          "n_NCI_cancer_drugs", "depmap_selectivity",
                          "top_tier_drug_examples"]],
], ignore_index=True)
src.to_csv(rf"{OUT_F}\Fig_scTWAS_druggability_sourcedata.csv", index=False)

print(summary.to_string(index=False))
print("genes with >=1 record:", int(out["in_ChEMBL_DB_TTD"].sum()),
      "| approved:", int((out["max_tier"] == 5).sum()),
      "| approved + pharmacological action:", int((out["n_approved_pharmacological"] > 0).sum()),
      "| listed oncology drug:", int((out["n_NCI_cancer_drugs"] > 0).sum()))
