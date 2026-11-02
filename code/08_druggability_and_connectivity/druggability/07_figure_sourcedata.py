#!/usr/bin/env python
"""Source data for CRC_scTWAS_CellType_Gene_Drug_Indication.{svg,png,pdf}.

    python src/07_figure_sourcedata.py

Writes a reviewer-readable workbook and the same content as a flat CSV:

    Definitions  what the figure encodes, the exact cancer-keyword rule, caveats
    Drugs        220 rows, one per plotted drug node (the outer ring)
    Genes         55 rows, one per plotted gene node (the inner ring)

Every indication column is computed HERE, verbatim from the two source fields in
63_druggable_genes.xlsx, not read back from the figure's display list -- so the
table is an independent check on the figure rather than a restatement of it. The
black/grey call derived here is asserted against the flag the figure actually
drew; the script fails if they disagree.

No cap: all cancer indications are listed (an earlier version showed only 3).

Reads   table/63_druggable_genes_advanced_55genes_193drugs.json   (what was drawn)
        table/63_druggable_genes_advanced_records.csv             (source fields)
Writes  figure/CRC_scTWAS_CellType_Gene_Drug_Indication_sourcedata.xlsx
        figure/CRC_scTWAS_CellType_Gene_Drug_Indication_sourcedata.csv
"""
import importlib.util
import json
import re
from pathlib import Path

import pandas as pd

ROOT = Path(r"W:\VUMC_projects\scTWAS\Drugs")
JSON_PATH = ROOT / "table" / "63_druggable_genes_advanced_55genes_193drugs.json"
REC_PATH = ROOT / "table" / "63_druggable_genes_advanced_records.csv"
OUT_XLSX = ROOT / "figure" / "CRC_scTWAS_CellType_Gene_Drug_Indication_sourcedata.xlsx"
OUT_CSV = OUT_XLSX.with_suffix(".csv")

GREEN, GREY_NODE = "#70c17f", "#808080"
BLACK, GREY_TEXT = "#000000", "#6e6e6e"

# The keyword lists live in the figure's own builder, so the table and the figure
# cannot drift apart.
_spec = importlib.util.spec_from_file_location("pij", ROOT / "src" / "protein_indication_json.py")
pij = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(pij)
TERMS = list(pij.CRC_TERMS) + list(pij.CANCER_TERMS)

IND_COLS = ["Drug Indications (DailyMed or ATC)", "Drug Indications (ClinicalTrials)"]

tree = json.loads(JSON_PATH.read_text(encoding="utf-8"))
rec = pd.read_csv(REC_PATH)
SHORT = tree.get("meta", {}).get("display_overrides", {})


def split(text):
    return [x.strip() for x in re.split(r"[;,|]", str(text)) if x.strip()]


def matched_terms(indication):
    low = indication.lower()
    return [t for t in TERMS if t in low]


# --- indications per gene-drug, verbatim from both source fields --------------
info = {}
for key, sub in rec.groupby(["Genes", "Drug Brand Name"]):
    seen, items = set(), []
    for col in IND_COLS:
        for cell in sub[col].dropna():
            for x in split(cell):
                if x.lower() not in seen:
                    seen.add(x.lower())
                    items.append(x)
    cancer = sorted((x for x in items if matched_terms(x)), key=str.lower)
    other = sorted((x for x in items if not matched_terms(x)), key=str.lower)
    hits = sorted({t for x in cancer for t in matched_terms(x)})
    info[key] = dict(all=items, cancer=cancer, other=other, terms=hits)

prov = rec.groupby(["Genes", "Drug Brand Name"]).agg(
    n_records=("Source", "size"),
    sources=("Source", lambda s: "|".join(sorted(set(s)))),
    max_phase_raw=("Max Phase", lambda s: "; ".join(sorted({str(x) for x in s.dropna()}))),
    development_stage=("Development stage (harmonised)", lambda s: "; ".join(sorted({str(x) for x in s.dropna()}))),
    molecular_ids=("MolecularID", lambda s: "; ".join(sorted({str(x) for x in s.dropna()}))),
    NCI_cancer_drug=("NCI Cancer Drugs", lambda s: bool(s.astype(str).str.upper().eq("TRUE").any())),
)

# --- one row per plotted drug node -------------------------------------------
drug_rows = []
for gene in tree["children"]:
    for drug in gene["children"]:
        key = (gene["name"], drug["name"])
        d = info[key]
        p = prov.loc[key]
        approved = drug["value"] == "true"
        has_cancer = bool(d["cancer"])
        # the figure's own flag: the drug node has >=1 cancer indication child
        assert has_cancer == bool(drug["children"]), f"flag disagrees with figure for {key}"
        drug_rows.append({
            "gene": gene["name"],
            "drug": drug["name"],
            "label_in_figure": SHORT.get(drug["name"], drug["name"]),
            "cancer_indication": "Yes" if has_cancer else "No",
            "label_colour": "Black" if has_cancer else "Grey",
            "approval_status": "Approved" if approved else "Clinical trials (phase II+)",
            "node_colour": "Green" if approved else "Grey",
            "n_indications_total": len(d["all"]),
            "n_cancer_indications": len(d["cancer"]),
            "cancer_indications": "; ".join(d["cancer"]),
            "matched_cancer_terms": "; ".join(d["terms"]),
            "other_indications": "; ".join(d["other"]),
            "development_stage": p.development_stage,
            "max_phase_raw": p.max_phase_raw,
            "sources": p.sources,
            "n_source_records": int(p.n_records),
            "molecular_ids": p.molecular_ids,
            "NCI_cancer_drug": bool(p.NCI_cancer_drug),
            "node_fill_hex": GREEN if approved else GREY_NODE,
            "label_fill_hex": BLACK if has_cancer else GREY_TEXT,
        })
drugs = pd.DataFrame(drug_rows)

# --- one row per plotted gene node -------------------------------------------
genes = (drugs.groupby("gene")
         .agg(n_drugs=("drug", "size"),
              n_approved=("approval_status", lambda s: int((s == "Approved").sum())),
              n_clinical=("approval_status", lambda s: int((s != "Approved").sum())),
              n_with_cancer_indication=("cancer_indication", lambda s: int((s == "Yes").sum())),
              n_without_cancer_indication=("cancer_indication", lambda s: int((s == "No").sum())))
         .reset_index())

# --- definitions sheet --------------------------------------------------------
n_black = int((drugs.cancer_indication == "Yes").sum())
n_nodata = int((drugs.n_indications_total == 0).sum())
defs = [
    ("Figure", "CRC_scTWAS_CellType_Gene_Drug_Indication.{svg,png,pdf}"),
    ("What is plotted", f"{len(genes)} CRC scTWAS risk genes (inner ring) and their "
                        f"{len(drugs)} gene-drug pairs (outer ring; {drugs.drug.nunique()} distinct compounds)"),
    ("Drug set", "Compounds that are approved or have reached phase II/III, from DrugBank, "
                 "ChEMBL, TTD and Open Targets"),
    ("Node colour", f"Green {GREEN} = approved; grey {GREY_NODE} = in clinical trials (phase II+)"),
    ("Label colour", f"Black {BLACK} = the drug has at least one cancer indication; "
                     f"grey {GREY_TEXT} = no cancer indication"),
    ("Cancer rule", "An indication counts as cancer if its lower-cased text CONTAINS any keyword "
                    "below. A drug is black if ANY ONE of its indications matches."),
    ("Cancer keywords", ", ".join(TERMS)),
    ("Indications searched", "Union of 'Drug Indications (DailyMed or ATC)' and "
                             "'Drug Indications (ClinicalTrials)' from 63_druggable_genes.xlsx; "
                             "split on ; , |, de-duplicated case-insensitively"),
    ("Completeness", "All indications are listed - no cap. Columns cancer_indications and "
                     "other_indications together are the full set."),
    ("Caveat: unknown vs negative", f"{n_nodata} of the grey drugs have NO indication text in either "
                                    "source field (n_indications_total = 0). They are unknown, not "
                                    "confirmed non-cancer, but are drawn the same grey."),
    ("Caveat: mixed indications", "Black marks at least one cancer indication; many black drugs also "
                                  "carry non-cancer indications (see other_indications)."),
    ("Caveat: keyword matching", "Free-text keyword matching, not an ontology lookup. Familial "
                                 "adenomatous polyposis is deliberately NOT counted (predisposition "
                                 "syndrome, not a malignancy)."),
    ("Counts", f"{n_black} drug nodes black / {len(drugs) - n_black} grey; "
               f"{int((drugs.node_colour == 'Green').sum())} approved / "
               f"{int((drugs.node_colour == 'Grey').sum())} phase II+"),
    ("Built by", "src/04_advanced_pairs_indication_json.py -> src/05_render_radial_tree.mjs -> "
                 "src/06_export_figure_png_pdf.py -> src/07_figure_sourcedata.py"),
]
definitions = pd.DataFrame(defs, columns=["Item", "Description"])

with pd.ExcelWriter(OUT_XLSX, engine="openpyxl") as xw:
    definitions.to_excel(xw, sheet_name="Definitions", index=False)
    drugs.to_excel(xw, sheet_name="Drugs", index=False)
    genes.to_excel(xw, sheet_name="Genes", index=False)
    widths = {"Definitions": [26, 110], "Drugs": [10, 34, 30, 10, 8, 14, 8, 8, 8, 60, 26, 60],
              "Genes": [10, 10, 12, 12, 16, 18]}
    for sheet, ws in xw.sheets.items():
        for i, w in enumerate(widths.get(sheet, []), start=1):
            ws.column_dimensions[ws.cell(row=1, column=i).column_letter].width = w
        ws.freeze_panes = "A2"

# flat CSV: the same two blocks, for scripted use
flat = pd.concat([genes.assign(element="gene"), drugs.assign(element="drug")], ignore_index=True)
flat = flat[["element"] + [c for c in flat.columns if c != "element"]]
flat.to_csv(OUT_CSV, index=False)

print(f"{len(drugs)} drug rows ({n_black} black / {len(drugs) - n_black} grey; "
      f"{n_nodata} with no indication text), {len(genes)} gene rows")
print(f"cancer indications listed: {int(drugs.n_cancer_indications.sum())} "
      f"(max {int(drugs.n_cancer_indications.max())} for one drug)")
print(OUT_XLSX)
