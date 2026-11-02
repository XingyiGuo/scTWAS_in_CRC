#!/usr/bin/env python
"""Support table for the druggable scTWAS risk genes, in 43_druggable_genes.xlsx format.

Rebuilds the drug-evidence table for the 347 scTWAS colorectal risk genes with the two
defects of 43_druggable_genes.xlsx corrected:

  * the 'Genes' column of sheet ChEMBL_DB_TTD is exploded on ';' and '|' first, so
    evidence held in protein-complex / target-family rows is not lost (PSMA7, PSMD3,
    PSMD13, RPS21, MAN2C1, HRAS);
  * the inclusion rule is applied on one harmonised stage ladder across the three
    sources, instead of ChEMBL phases 1-4 but only TTD 'Approved' and only DrugBank
    'approved' (which dropped TTD Phase 1-3 and DrugBank 'investigational' drugs).

Inclusion rule (row level): the drug is approved or has reached clinical development
  ChEMBL   max_phase 1-4
  DrugBank group list contains 'approved' or 'investigational'
  TTD      Approved / Registered / Phase 0-4 / Preregistration / Clinical trial
and, following the convention of the original file, an approved DrugBank row whose
group list contains 'withdrawn' without 'investigational' is dropped.

Layout of Sheet1 reproduces the original exactly: row 1 blank, header on row 2, data
from row 3, same 15 columns in the same order.  Four provenance columns are appended
after them, and a second sheet carries the gene-level summary.

Output: <ROOT>/Drugs/63_druggable_genes.xlsx  (+ .csv of both sheets)
"""
import re
import numpy as np
import pandas as pd

ROOT = r"W:\VUMC_projects\scTWAS"
XLSX = rf"{ROOT}\Drugs\Proteins_Drugs_20250321.xlsx"
OLD = rf"{ROOT}\Drugs\43_druggable_genes.xlsx"
DRUG_TBL = rf"{ROOT}\Drugs\table\scTWAS_347genes_druggability_ChEMBL_DB_TTD.csv"
OUT = rf"{ROOT}\Drugs\63_druggable_genes.xlsx"
OUT_T = rf"{ROOT}\Drugs\table"

COLS15 = ["Source", "Genes", "MolecularID", "Drug Brand Name", "Synonyms", "Max Phase",
          "action_group", "Drug Indications (DailyMed or ATC)",
          "Drug Indications (ClinicalTrials)", "NCI Cancer Drugs", "DrugResouces",
          "Drug Indications", "Approved/Clinical", "Cancer Indication",
          "Approved/Clinical (Cancer)"]

# action_type -> action_group, read off the 157 rows of 43_druggable_genes.xlsx and
# extended to the action_type values that occur elsewhere in the sheet.
GROUP = {**{k: "Inhibitory" for k in ["INHIBITOR", "ANTAGONIST", "NEGATIVE ALLOSTERIC MODULATOR",
                                      "BLOCKER", "INVERSE AGONIST", "DISRUPTING AGENT",
                                      "SEQUESTERING AGENT", "DEGRADER"]},
         **{k: "Activating" for k in ["AGONIST", "ACTIVATOR", "POTENTIATOR", "PARTIAL AGONIST",
                                      "POSITIVE ALLOSTERIC MODULATOR", "INDUCER", "OPENER",
                                      "STABILISER"]},
         **{k: "Binding" for k in ["BINDER", "BINDING", "BINDING AGENT", "ANTIBODY",
                                   "BINDER', 'ANTIBODY", "SUBSTRATE", "CHELATOR"]},
         **{k: "Modulatory" for k in ["MODULATOR", "LIGAND", "COFACTOR", "REGULATOR",
                                      "ALLOSTERIC MODULATOR"]},
         **{k: "Other" for k in ["OTHER", "UNKNOWN", "OTHER/UNKNOWN", "CROSS-LINKING AGENT",
                                 "RELEASING AGENT", "OXIDATIVE ENZYME", "HYDROLYTIC ENZYME",
                                 "PROTEOLYTIC ENZYME"]}}

CANCER_KW = ["NEOPLASM", "CANCER", "CARCINOMA", "TUMOR", "TUMOUR", "LEUKEMIA", "LYMPHOMA",
             "MELANOMA", "SARCOMA", "GLIOMA", "GLIOBLASTOMA", "MYELOMA", "ADENOMA",
             "BLASTOMA", "MALIGNAN", "METASTA", "MESOTHELIOMA", "POLYPS"]

TIER_NAME = {5: "Approved drug", 4: "Phase 3 / registration", 3: "Phase 2",
             2: "Phase 1 / clinical", 1: "Preclinical / experimental ligand",
             0: "Discontinued or withdrawn only"}
TTD_MAP = {**{k: 5 for k in ["approved", "approved (orphan drug)", "registered", "phase 4"]},
           **{k: 4 for k in ["phase 3", "phase 2/3", "preregistration", "nda filed",
                             "bla submitted", "application submitted", "approval submitted"]},
           **{k: 3 for k in ["phase 2", "phase 2a", "phase 2b", "phase 1/2",
                             "phase 1b/2a", "phase 1/2a"]},
           **{k: 2 for k in ["phase 1", "phase 1b", "phase 0", "clinical trial"]},
           **{k: 1 for k in ["preclinical", "investigative", "patented"]}}


def tier_of(src, mp):
    s = "" if pd.isna(mp) else str(mp).strip().strip("'\"").lower()
    if src == "ChEMBL":
        return {"4": 5, "3": 4, "2": 3, "1": 2, "0": 1}.get(s, 1)
    if src == "TTD":
        if s.startswith(("discontinued", "terminated", "withdrawn")):
            return 0
        return TTD_MAP.get(s, 1)
    toks = {t.strip().strip("'\" ") for t in re.split(r",", s)}
    if "approved" in toks:
        return 5
    if "investigational" in toks:
        return 2
    if toks & {"experimental", "nutraceutical", "vet_approved", "illicit"}:
        return 1
    if "withdrawn" in toks:
        return 0
    return 1


def indications(a, b):
    parts = [x for x in (a, b) if isinstance(x, str) and x.strip()]
    return ";".join(parts).upper() if parts else np.nan


dd = pd.read_excel(XLSX, sheet_name="ChEMBL_DB_TTD")
per_gene = pd.read_csv(DRUG_TBL)
old43 = set(pd.read_excel(OLD, sheet_name="Sheet1", header=1)["Genes"].astype(str).str.strip())

dd["record_Genes"] = dd["Genes"].astype(str)
dd["gene_list"] = dd["record_Genes"].str.split(r"\s*[;|]\s*")
dl = dd.explode("gene_list")
dl["gene_sym"] = dl["gene_list"].str.strip()
dl = dl[dl.gene_sym.ne("") & dl.gene_sym.ne("nan")].copy()
dl["tier"] = [tier_of(s, m) for s, m in zip(dl["Source"], dl["Max Phase"])]

genes61 = set(per_gene.loc[per_gene.max_tier >= 2, "gene_name"])
mp = dl["Max Phase"].astype(str).str.lower()
withdrawn_only = (dl.Source == "DB") & mp.str.contains("withdrawn") & ~mp.str.contains("investigational")
keep = dl.gene_sym.isin(genes61) & dl.tier.ge(2) & ~withdrawn_only
r = dl[keep].copy()

r["Genes"] = r["gene_sym"]
r["action_group"] = (r["action_type"].astype(str).str.strip().str.strip("'\"").str.upper()
                     .map(GROUP))
r["Drug Indications"] = [indications(a, b) for a, b in
                         zip(r["Drug Indications (DailyMed or ATC)"],
                             r["Drug Indications (ClinicalTrials)"])]
r["Approved/Clinical"] = True
r["Cancer Indication"] = [bool(isinstance(s, str) and any(k in s for k in CANCER_KW))
                          for s in r["Drug Indications"]]
r["Approved/Clinical (Cancer)"] = r["Approved/Clinical"] & r["Cancer Indication"]
r["Gene record as stored in sheet"] = r["record_Genes"]
r["Evidence scope"] = np.where(r["record_Genes"].str.strip() == r["gene_sym"],
                               "direct symbol row", "complex or family row")
r["Development stage (harmonised)"] = r["tier"].map(TIER_NAME)
r["In 43_druggable_genes.xlsx"] = r["gene_sym"].isin(old43)

EXTRA = ["Gene record as stored in sheet", "Evidence scope",
         "Development stage (harmonised)", "In 43_druggable_genes.xlsx"]

# --- fourth resource: the workbook's OpenTargets sheet -----------------------
# It carries no development phase, only IsApproved; Open Targets lists drugs with
# clinical precedence, so a non-approved entry is recorded as clinical, phase unstated.
ot = pd.read_excel(XLSX, sheet_name="OpenTargets")
ot["record_Genes"] = ot["Genes"].astype(str)
ot["gene_list"] = ot["record_Genes"].str.split(r"\s*[;|]\s*")
ot = ot.explode("gene_list")
ot["gene_sym"] = ot["gene_list"].str.strip()
o = ot[ot.gene_sym.isin(per_gene.gene_name)].copy()
o["Genes"] = o["gene_sym"]
o["Synonyms"] = np.nan
o["Max Phase"] = np.where(o["IsApproved"], "Approved", "Clinical (phase not stated)")
o["action_group"] = np.nan
o["Drug Indications (DailyMed or ATC)"] = np.nan
o["Drug Indications (ClinicalTrials)"] = np.nan
o["Drug Indications"] = np.nan
o["NCI Cancer Drugs"] = o["Commercial CancerDrug"].astype(bool)
o["Approved/Clinical"] = True
o["Cancer Indication"] = o["NCI Cancer Drugs"]
o["Approved/Clinical (Cancer)"] = o["NCI Cancer Drugs"]
o["Source"] = "OpenTargets"
o["Gene record as stored in sheet"] = o["record_Genes"]
o["Evidence scope"] = np.where(o["record_Genes"].str.strip() == o["gene_sym"],
                               "direct symbol row", "complex or family row")
o["Development stage (harmonised)"] = np.where(o["IsApproved"], "Approved drug",
                                               "Clinical, phase not stated")
o["In 43_druggable_genes.xlsx"] = o["gene_sym"].isin(old43)

rows = (pd.concat([r[COLS15 + EXTRA], o[COLS15 + EXTRA]], ignore_index=True)
        .drop_duplicates(subset=["Source", "Genes", "MolecularID", "Drug Brand Name"])
        .sort_values(["Genes", "Source", "Drug Brand Name"])
        .reset_index(drop=True))
genes63 = set(rows["Genes"])
assert len(genes63) == 63, len(genes63)
assert old43 <= genes63

# gene-level fields recomputed so the OpenTargets-only genes are not left blank
ot_only = sorted(genes63 - genes61)
per_gene = per_gene.set_index("gene_name")
for g in ot_only:
    sl = rows[rows.Genes == g]
    per_gene.loc[g, "druggability_tier"] = ("Approved drug" if (sl["Development stage (harmonised)"]
                                                                == "Approved drug").any()
                                            else "Clinical, phase not stated")
    per_gene.loc[g, "n_drugs"] = sl["Drug Brand Name"].str.upper().nunique()
    per_gene.loc[g, "n_approved_drugs"] = sl.loc[sl["Development stage (harmonised)"]
                                                 == "Approved drug", "Drug Brand Name"].str.upper().nunique()
    per_gene.loc[g, "n_approved_pharmacological"] = 0     # OpenTargets carries no action_type
    per_gene.loc[g, "n_clinical_drugs"] = sl.loc[sl["Development stage (harmonised)"]
                                                 != "Approved drug", "Drug Brand Name"].str.upper().nunique()
    per_gene.loc[g, "n_NCI_cancer_drugs"] = sl.loc[sl["NCI Cancer Drugs"], "Drug Brand Name"].str.upper().nunique()
    per_gene.loc[g, "sources"] = "OpenTargets"
    per_gene.loc[g, "top_tier_drug_examples"] = "; ".join(sorted(
        sl.loc[sl["Development stage (harmonised)"] == "Approved drug", "Drug Brand Name"].str.upper().unique())[:6])
    per_gene.loc[g, "NCI_cancer_drug_names"] = "; ".join(sorted(
        sl.loc[sl["NCI Cancer Drugs"], "Drug Brand Name"].str.upper().unique())[:8])
per_gene = per_gene.reset_index()
genes61 = genes63

summ = per_gene[per_gene.gene_name.isin(genes61)][[
    "gene_name", "druggability_tier", "n_drugs", "n_approved_drugs",
    "n_approved_pharmacological", "n_clinical_drugs", "n_NCI_cancer_drugs", "sources",
    "top_tier_drug_examples", "NCI_cancer_drug_names", "rep_model", "rep_Z", "rep_P",
    "rep_direction", "n_models", "depmap_median", "depmap_selectivity"]].copy()
summ.insert(1, "In 43_druggable_genes.xlsx", summ.gene_name.isin(old43))
n_rec = rows.groupby("Genes").size().rename("n_records_in_this_table")
n_cx = (rows["Evidence scope"].eq("complex or family row").groupby(rows["Genes"]).mean()
        .eq(1).rename("evidence_only_from_complex_rows"))
summ = (summ.merge(n_rec, left_on="gene_name", right_index=True)
            .merge(n_cx, left_on="gene_name", right_index=True)
            .sort_values(["druggability_tier", "n_approved_pharmacological", "gene_name"],
                         ascending=[True, False, True]))
assert len(summ) == 63

with pd.ExcelWriter(OUT, engine="openpyxl") as xw:
    rows.to_excel(xw, sheet_name="Sheet1", index=False, startrow=1)   # blank row 1
    summ.to_excel(xw, sheet_name="Gene_summary", index=False)
rows.to_csv(rf"{OUT_T}\63_druggable_genes_drug_records.csv", index=False)
summ.to_csv(rf"{OUT_T}\63_druggable_genes_summary.csv", index=False)

nm = rows["Drug Brand Name"].astype(str).str.upper().str.replace(r"[^A-Z0-9]", "", regex=True)
appr = rows["Development stage (harmonised)"].eq("Approved drug")
print("rows:", len(rows), "genes:", rows['Genes'].nunique(),
      "| new genes vs the 43:", sorted(set(rows['Genes']) - old43))
print(summ["druggability_tier"].value_counts().to_string())
print("distinct compounds:", nm.nunique(), "| approved:", nm[appr].nunique(),
      "| NCI oncology:", nm[rows['NCI Cancer Drugs'].astype(bool)].nunique())
print("rows from complex/family records:", int(rows["Evidence scope"].eq("complex or family row").sum()),
      "| by source:", rows.Source.value_counts().to_dict())
