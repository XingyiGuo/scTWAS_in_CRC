#!/usr/bin/env python3
"""Build GSE201336_atac_manifest.csv joining GEO GSM accessions to Becker et al. sample metadata."""
import csv, re

# GEO GSM -> title, transcribed from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE201336
GEO = """GSM6058770|A001-C-002, Dounce
GSM6058771|A001-C-007, Dounce
GSM6058772|A001-C-014, Dounce
GSM6058773|A001-C-023, Replicate1, Dounce
GSM6058774|A001-C-023, Replicate2, Dounce
GSM6058775|A001-C-104, Replicate1, Dounce
GSM6058776|A001-C-104, Replicate2, Dounce
GSM6058777|A001-C-108, Dounce
GSM6058778|A001-C-119, Dounce
GSM6058779|A001-C-123, Dounce
GSM6058780|A001-C-124, Replicate1, Dounce
GSM6058781|A001-C-124, Replicate2, Dounce
GSM6058782|A001-C-202, Dounce
GSM6058783|A001-C-203, Dounce
GSM6058784|A001-C-207, Dounce
GSM6058785|A001-C-218, Dounce
GSM6058786|A001-C-223, Dounce
GSM6058787|A002-C-010, Replicate1, Dounce
GSM6058788|A002-C-010, Replicate2, Dounce
GSM6058789|A002-C-010, S2
GSM6058790|A002-C-016, Dounce
GSM6058791|A002-C-021, Dounce
GSM6058792|A002-C-024, Dounce
GSM6058793|A002-C-025, Dounce
GSM6058794|A002-C-025, S2
GSM6058795|A002-C-106, Dounce
GSM6058796|A002-C-114, Dounce
GSM6058797|A002-C-116, Dounce
GSM6058798|A002-C-116, S2
GSM6058799|A002-C-121, Dounce
GSM6058800|A002-C-121, S2
GSM6058801|A002-C-201, Dounce
GSM6058802|A002-C-202, Dounce
GSM6058803|A002-C-203, Dounce
GSM6058804|A002-C-204, Dounce
GSM6058805|A002-C-205, Dounce
GSM6058806|A002-C-212, Dounce
GSM6058807|A008-C-008, Dounce
GSM6058808|A008-C-015, Dounce
GSM6058809|A010-C-018, Dounce
GSM6058810|A010-C-023, Dounce
GSM6058811|A014-C-001, Dounce
GSM6058812|A014-C-002, Dounce
GSM6058813|A014-C-008, Dounce
GSM6058814|A014-C-013, Dounce
GSM6058815|A014-C-040, Dounce
GSM6058816|A014-C-043, Dounce
GSM6058817|A014-C-052, Dounce
GSM6058818|A014-C-054, Dounce
GSM6058819|A014-C-101, Dounce
GSM6058820|A014-C-108, Dounce
GSM6058821|A014-C-111, Dounce
GSM6058822|A014-C-114, Dounce
GSM6058823|A014-C-201, Dounce
GSM6058824|A015-C-001, Dounce
GSM6058825|A015-C-002, Dounce
GSM6058826|A015-C-005, Dounce
GSM6058827|A015-C-006, Dounce
GSM6058828|A015-C-008, Dounce
GSM6058829|A015-C-010, Dounce
GSM6058830|A015-C-102, Dounce
GSM6058831|A015-C-104, Dounce
GSM6058832|A015-C-106, Dounce
GSM6058833|A015-C-109, Dounce
GSM6058834|A015-C-202, Dounce
GSM6058835|A015-C-203, Dounce
GSM6058836|A015-C-204, Dounce
GSM6058837|A015-C-206, Dounce
GSM6058838|A015-C-208, Dounce
GSM6058839|A018-C-013, Dounce
GSM6058840|A018-C-020, Dounce
GSM6058841|A022-C-002, Dounce
GSM6058842|CRC-1-8810, Dounce
GSM6058843|CRC-2-15564, Dounce
GSM6058844|CRC-3-11773, Dounce
GSM6058845|CRC-4-8456, Dounce
GSM6058846|F007, Dounce
GSM6058847|F034, Dounce
GSM6058848|F072B, Dounce
GSM6058849|F091, Dounce
GSM6058850|B001-A-301, Dounce
GSM6058851|B001-A-302, Dounce
GSM6058852|B001-A-401, Dounce
GSM6058853|B001-A-406, Dounce
GSM6058854|B001-A-501, Dounce
GSM6058855|B004-A-004, Replicate1, Dounce
GSM6058856|B004-A-004, Replicate2, Dounce
GSM6058857|B004-A-008, Dounce
GSM6058858|B004-A-204, Dounce"""

# SampleNameOnly -> (Donor, FAP, DiseaseState, Location)
# from winstonbecker/scCRC_continuum/hubmap_htan_metadata_atac_and_rna_final.csv
META = {
 "A001-C-002":("A001","Y","Polyp","Descending"), "A001-C-007":("A001","Y","Adenocarcinoma","Descending"),
 "A001-C-014":("A001","Y","Polyp","Descending"), "A001-C-023":("A001","Y","Unaffected","Descending"),
 "A001-C-104":("A001","Y","Polyp","Transverse"), "A001-C-108":("A001","Y","Polyp","Transverse"),
 "A001-C-119":("A001","Y","Polyp","Transverse"), "A001-C-123":("A001","Y","Unaffected","Transverse"),
 "A001-C-124":("A001","Y","Unaffected","Transverse"), "A001-C-202":("A001","Y","Polyp","Ascending"),
 "A001-C-203":("A001","Y","Polyp","Ascending"), "A001-C-207":("A001","Y","Polyp","Ascending"),
 "A001-C-218":("A001","Y","Polyp","Ascending"), "A001-C-223":("A001","Y","Unaffected","Ascending"),
 "A002-C-010":("A002","Y","Polyp","Descending"), "A002-C-016":("A002","Y","Polyp","Descending"),
 "A002-C-021":("A002","Y","Polyp","Descending"), "A002-C-024":("A002","Y","Unaffected","Descending"),
 "A002-C-025":("A002","Y","Unaffected","Descending"), "A002-C-106":("A002","Y","Polyp","Transverse"),
 "A002-C-114":("A002","Y","Polyp","Transverse"), "A002-C-116":("A002","Y","Polyp","Transverse"),
 "A002-C-121":("A002","Y","Unaffected","Transverse"), "A002-C-201":("A002","Y","Polyp","Ascending"),
 "A002-C-202":("A002","Y","Polyp","Ascending"), "A002-C-203":("A002","Y","Polyp","Ascending"),
 "A002-C-204":("A002","Y","Polyp","Ascending"), "A002-C-205":("A002","Y","Polyp","Ascending"),
 "A002-C-212":("A002","Y","Unaffected","Ascending"),
 "A008-E-008":("A008","Y","Polyp","Ascending"), "A008-E-015":("A008","Y","Polyp","Descending"),
 "A010-E-018":("A010","Y","Polyp","Ascending"), "A010-E-023":("A010","Y","Polyp","Ascending"),
 "A014-C-001":("A014","Y","Polyp","Descending"), "A014-C-002":("A014","Y","Polyp","Rectum"),
 "A014-C-008":("A014","Y","Polyp","Rectum"), "A014-C-013":("A014","Y","Polyp","Descending"),
 "A014-C-040":("A014","Y","Polyp","Descending"), "A014-C-043":("A014","Y","Polyp","Descending"),
 "A014-C-052":("A014","Y","Unaffected","Rectum"), "A014-C-054":("A014","Y","Unaffected","Descending"),
 "A014-C-101":("A014","Y","Polyp","Transverse"), "A014-C-108":("A014","Y","Polyp","Transverse"),
 "A014-C-111":("A014","Y","Unaffected","Transverse"), "A014-C-114":("A014","Y","Unaffected","Transverse"),
 "A014-C-201":("A014","Y","Unaffected","Ascending"),
 "A015-C-001":("A015","Y","Adenocarcinoma","Rectum"), "A015-C-002":("A015","Y","Polyp","Descending"),
 "A015-C-005":("A015","Y","Polyp","Rectum"), "A015-C-006":("A015","Y","Polyp","Rectum"),
 "A015-C-008":("A015","Y","Unaffected","Rectum"), "A015-C-010":("A015","Y","Unaffected","Descending"),
 "A015-C-102":("A015","Y","Polyp","Transverse"), "A015-C-104":("A015","Y","Polyp","Transverse"),
 "A015-C-106":("A015","Y","Polyp","Transverse"), "A015-C-109":("A015","Y","Unaffected","Transverse"),
 "A015-C-202":("A015","Y","Polyp","Ascending"), "A015-C-203":("A015","Y","Polyp","Ascending"),
 "A015-C-204":("A015","Y","Polyp","Ascending"), "A015-C-206":("A015","Y","Unaffected","Ascending"),
 "A015-C-208":("A015","Y","Unaffected","Ascending"),
 "A018-E-013":("A018","Y","Polyp","Ascending"), "A018-E-020":("A018","Y","Polyp","Transverse"),
 "A022-E-002":("A022","N","Polyp","Hepatic flexure"),
 "B001-A-301":("B001","N","Normal","Sigmoid"), "B001-A-302":("B001","N","Normal","Sigmoid"),
 "B001-A-401":("B001","N","Normal","Transverse"), "B001-A-406":("B001","N","Normal","Descending"),
 "B001-A-501":("B001","N","Normal","Ascending"),
 "B004-A-004":("B004","N","Normal","Sigmoid"), "B004-A-008":("B004","N","Normal","Descending"),
 "B004-A-204":("B004","N","Normal","Ascending"),
 "F007":("F","Y","Polyp","Ascending"), "F034":("F","Y","Polyp","Transverse"),
 "F072B":("F","Y","Polyp","Descending"), "F091":("F","Y","Polyp","Rectum"),
 "CRC-1-8810":("CRC1","N","Adenocarcinoma",""), "CRC-2-15564":("CRC2","N","Adenocarcinoma",""),
 "CRC-3-11773":("CRC3","N","Adenocarcinoma",""), "CRC-4-8456":("CRC4","N","Adenocarcinoma",""),
}

# Samples annotated Unaffected in the metadata but that the authors treated as their own
# DifferentialGroup because pathology found dysplasia. Kept out of the "normal" condition.
DYSPLASTIC_UNAFFECTED = {"A015-C-109", "A015-C-206", "A015-C-208"}

# The authors' ATAC_Sample IDs do NOT follow a uniform suffix rule, so a naive
# "Replicate1 -> -D-R1" transform produces names that will not match the cell-type
# annotation files. These are the exceptions, taken verbatim from
# hubmap_htan_metadata_atac_and_rna_final.csv.
#
#   * B004-A-004 Replicate1 is "-D" (not "-D-R1"); Replicate2 is "-D-R2".
#   * A002-C-010 Replicate2 is "-D" (no "-D-R2" row exists).
#   * A002-C-202 carries an "-OCT" suffix.
#   * A015-C-010 exists only as "-D-R2" even though GEO lists no replicate.
#
# CAVEAT: GEO's "Replicate1"/"Replicate2" labels cannot be resolved against the authors'
# IDs from names alone for B004-A-004 and A002-C-010. The assignment below is the most
# consistent reading of the metadata + scATAC_construct_epithelial_reference.R Arrow list.
# Verify with `check_annotation_join.R` before trusting downstream peaks.
ARROW_NAME_OVERRIDES = {
    "GSM6058788": "A002-C-010-D",     # A002-C-010, Replicate2
    "GSM6058802": "A002-C-202-D-OCT", # A002-C-202, Dounce
    "GSM6058829": "A015-C-010-D-R2",  # A015-C-010, Dounce
    "GSM6058855": "B004-A-004-D",     # B004-A-004, Replicate1
}

# GEO writes these donors' samples with a "-C-" infix; the authors' metadata uses "-E-".
def normalize(base):
    for d in ("A008", "A010", "A018", "A022"):
        if base.startswith(d + "-C-"):
            return base.replace("-C-", "-E-", 1)
    return base


def parse(title):
    parts = [p.strip() for p in title.split(",")]
    base = parts[0]
    tags = parts[1:]
    if "S2" in tags:
        suffix, diss = "-S2", "S2"
    else:
        diss = "Dounce"
        if "Replicate1" in tags:
            suffix = "-D-R1"
        elif "Replicate2" in tags:
            suffix = "-D-R2"
        else:
            suffix = "-D"
    # authors' IDs use the "-E-" infix for A008/A010/A018/A022
    return base, normalize(base) + suffix, diss


rows = []
for line in GEO.strip().split("\n"):
    gsm, title = line.split("|")
    base, atac_sample, diss = parse(title)
    atac_sample = ARROW_NAME_OVERRIDES.get(gsm, atac_sample)
    key = normalize(base)
    if key not in META:
        raise SystemExit(f"unmatched sample: {base} ({gsm})")
    donor, fap, disease, loc = META[key]

    if disease == "Normal":
        condition = "Normal"
    elif disease == "Unaffected":
        condition = "Unaffected_dysplastic" if key in DYSPLASTIC_UNAFFECTED else "Unaffected"
    elif disease == "Polyp":
        condition = "Adenoma"
    else:
        condition = "CRC"

    rows.append(dict(
        GSM=gsm, ATAC_Sample=atac_sample, SampleNameOnly=key, Donor=donor, FAP=fap,
        DiseaseState=disease, Condition=condition, Location=loc, Dissociation=diss,
        IsReplicate="Y" if re.search(r"-R[12]$", atac_sample) else "N",
        SupplDir=f"https://ftp.ncbi.nlm.nih.gov/geo/samples/{gsm[:7]}nnn/{gsm}/suppl/",
    ))

FIELDS = ["GSM", "ATAC_Sample", "SampleNameOnly", "Donor", "FAP", "DiseaseState",
          "Condition", "Location", "Dissociation", "IsReplicate", "SupplDir"]

with open("GSE201336_atac_manifest.csv", "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=FIELDS)
    w.writeheader()
    w.writerows(rows)

from collections import Counter
print(f"{len(rows)} samples written")
for k, v in sorted(Counter(r["Condition"] for r in rows).items()):
    print(f"  {k:24s} {v:3d}")
print("\ndonors:", len({r['Donor'] for r in rows}))
assert len(rows) == 89, len(rows)
assert len({r["ATAC_Sample"] for r in rows}) == 89, "duplicate ATAC_Sample"
