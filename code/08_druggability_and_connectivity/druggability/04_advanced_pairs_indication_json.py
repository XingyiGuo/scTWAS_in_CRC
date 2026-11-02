#!/usr/bin/env python
"""Radial-tree JSON for the clinically advanced drug-target pairs (Figure 6d).

Supports the manuscript sentence: "we filtered the initial set to compounds that were
approved or had reached phase II or III trials, leaving 55 genes linked to 193 drugs".

The filter is identical to the one in 03_druggability_paragraph.py: rows whose
harmonised development stage is "Approved drug", "Phase 3 / registration" or "Phase 2".
Compounds are counted on the normalised name (upper-case, alphanumerics only), so the
193 in the prose and the 193 distinct drug labels in the JSON are the same quantity.

Hierarchy is built by diseases_protein_drugs_instruction_json() in
src/protein_indication_json.py:

    Colorectal cancer -> gene -> drug (value = "true" if approved/phase 4) -> <=3 cancer indications

Reads   Drugs/63_druggable_genes.xlsx          (sheet Sheet1, header on row 2)
Writes  Drugs/table/63_druggable_genes_advanced_55genes_193drugs.json
        Drugs/table/63_druggable_genes_advanced_records.csv   (the 268 filtered records)
"""
import importlib.util
import json
import re
from pathlib import Path

import pandas as pd

ROOT = Path(r"W:\VUMC_projects\scTWAS\Drugs")
SRC = ROOT / "src"
TAB = ROOT / "table"
XLSX = ROOT / "63_druggable_genes.xlsx"
OUT_JSON = TAB / "63_druggable_genes_advanced_55genes_193drugs.json"
OUT_CSV = TAB / "63_druggable_genes_advanced_records.csv"

ADVANCED_STAGES = ["Approved drug", "Phase 3 / registration", "Phase 2"]

# Figure-only shortenings for names too long to letter into a ring. The full name
# stays the node name in the JSON and the `drug` column of the source data; the
# figure and the source data's `display_label` column use the short form.
# "BATs" = bispecific antibody-armed activated T cells, the usual short form for
# this construct in the HER2 literature.
DISPLAY_OVERRIDES = {
    "Anti-CD3 And anti-Her2/neu Bispecific Antibody-armed Activated T Cells":
        "Anti-CD3 \u00d7 anti-HER2 BATs",
}

# --- import the hierarchy builder from the sibling module --------------------
spec = importlib.util.spec_from_file_location(
    "protein_indication_json", SRC / "protein_indication_json.py"
)
pij = importlib.util.module_from_spec(spec)
spec.loader.exec_module(pij)

# --- load the 63-gene drug records -------------------------------------------
rec = pd.read_excel(XLSX, sheet_name="Sheet1", header=1)

norm = lambda s: s.astype(str).str.upper().str.replace(r"[^A-Z0-9]", "", regex=True)
rec["cpd"] = norm(rec["Drug Brand Name"])

adv = rec[rec["Development stage (harmonised)"].isin(ADVANCED_STAGES)].copy()

# --- one display label per compound ------------------------------------------
# The same compound is spelled differently across sources (DISULFIRAM / Disulfiram,
# KN-026 / KN026). Collapse to a single label per normalised key, otherwise the tree
# would carry 216 drug nodes for 193 compounds.
KEEP_UPPER = re.compile(r"\d")


def _smart_title(name: str) -> str:
    out = []
    for tok in name.split():
        if KEEP_UPPER.search(tok) or (tok.isupper() and len(tok) <= 4):
            out.append(tok)  # codes and acronyms: KN-026, ATP, AMP
        else:
            out.append(tok.capitalize())
    return " ".join(out)


def _canonical(variants) -> str:
    # prefer the spelling that keeps separators, then the longest, then alphabetical
    best = sorted(variants, key=lambda v: (-len(re.sub(r"[A-Za-z0-9]", "", v)), -len(v), v))[0]
    return _smart_title(best.strip())


label = {k: _canonical(v) for k, v in adv.groupby("cpd")["Drug Brand Name"].agg(set).items()}
adv["Drug Brand Name"] = adv["cpd"].map(label)

n_genes = adv["Genes"].nunique()
n_cpd = adv["cpd"].nunique()
assert n_genes == 55, n_genes
assert n_cpd == 193, n_cpd
assert adv["Drug Brand Name"].nunique() == n_cpd

# --- build the hierarchy ------------------------------------------------------
tree = pij.diseases_protein_drugs_instruction_json(adv)

n_gene_nodes = len(tree["children"])
n_drug_nodes = sum(len(g["children"]) for g in tree["children"])
n_drug_labels = len({d["name"] for g in tree["children"] for d in g["children"]})
n_ind_nodes = sum(len(d["children"]) for g in tree["children"] for d in g["children"])
n_drugs_with_ind = sum(
    1 for g in tree["children"] for d in g["children"] if d["children"]
)
n_approved_nodes = sum(
    1 for g in tree["children"] for d in g["children"] if d["value"] == "true"
)

tree["meta"] = {
    "description": (
        "scTWAS CRC risk genes: drug-target pairs restricted to compounds that are "
        "approved or have reached phase II or III trials"
    ),
    "source_table": XLSX.name + " (sheet Sheet1)",
    "filter": "Development stage (harmonised) in " + str(ADVANCED_STAGES),
    "compound_key": "Drug Brand Name upper-cased, non-alphanumerics removed",
    "n_genes": int(n_genes),
    "n_drugs": int(n_cpd),
    "n_gene_drug_records": int(len(adv)),
    "n_gene_drug_pairs": int(n_drug_nodes),
    "n_drug_nodes_approved_or_phase4": int(n_approved_nodes),
    "n_drug_nodes_with_cancer_indication": int(n_drugs_with_ind),
    "n_indication_nodes": int(n_ind_nodes),
    "node_value": "drug node value 'true' = approved or phase 4 in at least one record",
    "display_overrides": DISPLAY_OVERRIDES,
}
assert set(DISPLAY_OVERRIDES) <= set(adv["Drug Brand Name"]), "override name not in the data"

OUT_JSON.write_text(json.dumps(tree, indent=2, ensure_ascii=False), encoding="utf-8")
adv.drop(columns=["cpd"]).to_csv(OUT_CSV, index=False)

print(
    f"genes {n_gene_nodes} | drug labels {n_drug_labels} | gene-drug nodes {n_drug_nodes} "
    f"| records {len(adv)} | approved nodes {n_approved_nodes} "
    f"| drugs with indication {n_drugs_with_ind} | indication nodes {n_ind_nodes}"
)
print(OUT_JSON)
