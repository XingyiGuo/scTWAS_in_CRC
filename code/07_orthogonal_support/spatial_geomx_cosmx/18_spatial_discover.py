#!/usr/bin/env python3
"""
18 - Discovery pass over the Zenodo spatial-omics download.

Runs on the user's machine (pandas/numpy/matplotlib only - no scipy). It finds the
GeoMx and CosMx archives wherever they were saved, extracts them, and reports the
structure: file names, shapes, column names, and the candidate annotation columns that
carry region (normal / adenoma / tumour core / invasive front), patient, and cell type.
Nothing is analysed here; this exists so the analysis scripts are written against the
real schema rather than a guess.

Usage:  python3 18_spatial_discover.py <search-dir> [<search-dir> ...]
"""
import os, sys, zipfile, json, glob

SEARCH = sys.argv[1:] or [os.path.expanduser("~")]
WORK = os.path.join(os.path.expanduser("~"), "spatial_work")
os.makedirs(WORK, exist_ok=True)

WANT = ["Nanostring_GeoMx", "Nanostring_CosMx_exp2", "Nanostring_CosMx_exp1"]

found = {}
for root in SEARCH:
    for dirpath, dirnames, filenames in os.walk(root):
        depth = dirpath[len(root):].count(os.sep)
        if depth > 3:
            dirnames[:] = []
            continue
        for f in filenames:
            for w in WANT:
                if f.startswith(w) and f.endswith(".zip"):
                    p = os.path.join(dirpath, f)
                    found.setdefault(w, p)
            if f.lower().startswith("readme") and "zenodo" in f.lower():
                found.setdefault("README", os.path.join(dirpath, f))

print("=== archives found ===")
for k, v in found.items():
    print(f"  {k:24s} {os.path.getsize(v)/1e6:10.1f} MB  {v}")
if not found:
    print("  none yet")
    sys.exit(0)

for key, zp in found.items():
    if not zp.endswith(".zip"):
        print(f"\n=== {key} (text) ===")
        print(open(zp, encoding="utf-8", errors="replace").read()[:8000])
        continue
    dest = os.path.join(WORK, key)
    if not os.path.isdir(dest):
        print(f"\nextracting {key} -> {dest}")
        with zipfile.ZipFile(zp) as z:
            z.extractall(dest)
    print(f"\n=== {key} contents ===")
    entries = []
    for dirpath, _, filenames in os.walk(dest):
        for f in filenames:
            p = os.path.join(dirpath, f)
            entries.append((os.path.getsize(p), os.path.relpath(p, dest)))
    entries.sort(reverse=True)
    for s, rel in entries[:60]:
        print(f"  {s/1e6:10.2f} MB  {rel}")
    print(f"  ... {len(entries)} files total")

# ------------------------------------------------------------------ peek inside
import pandas as pd

REGION_HINT = ["region", "segment", "roi", "tissue", "stage", "class", "anno", "type",
               "compartment", "location", "area", "group", "label", "niche", "zone"]
CELL_HINT = ["cell_type", "celltype", "cluster", "annotation", "subtype", "ident",
             "leiden", "louvain", "cellclass"]
COORD_HINT = ["x_", "_x", "centroid", "global", "local", "coord", "px", "position"]


def peek(path, nrows=4):
    ext = os.path.splitext(path)[1].lower()
    try:
        if ext in (".csv", ".txt", ".tsv", ".gz"):
            sep = "\t" if ext in (".txt", ".tsv") else ","
            head = pd.read_csv(path, sep=sep, nrows=nrows, low_memory=False)
            n = sum(1 for _ in open(path, "rb")) - 1 if os.path.getsize(path) < 4e8 else None
            return head, n
        if ext in (".xlsx", ".xls"):
            head = pd.read_excel(path, nrows=nrows)
            return head, None
    except Exception as e:
        print(f"    ! could not read: {type(e).__name__}: {e}")
    return None, None


print("\n\n================= SCHEMA =================")
for key in found:
    dest = os.path.join(WORK, key)
    if not os.path.isdir(dest):
        continue
    tabs = []
    for dirpath, _, filenames in os.walk(dest):
        for f in filenames:
            if os.path.splitext(f)[1].lower() in (".csv", ".txt", ".tsv", ".xlsx", ".gz"):
                tabs.append(os.path.join(dirpath, f))
    tabs.sort(key=os.path.getsize, reverse=True)
    print(f"\n########## {key}: {len(tabs)} tabular files")
    for p in tabs[:14]:
        rel = os.path.relpath(p, dest)
        head, n = peek(p)
        print(f"\n--- {rel}   ({os.path.getsize(p)/1e6:.2f} MB"
              + (f", {n} rows)" if n is not None else ")"))
        if head is None:
            continue
        cols = list(head.columns)
        print(f"    {len(cols)} columns")
        print("    columns:", cols[:40], "..." if len(cols) > 40 else "")
        cand_r = [c for c in cols if any(h in str(c).lower() for h in REGION_HINT)]
        cand_c = [c for c in cols if any(h in str(c).lower() for h in CELL_HINT)]
        cand_xy = [c for c in cols if any(h in str(c).lower() for h in COORD_HINT)]
        if cand_r:
            print("    region-like :", cand_r[:12])
        if cand_c:
            print("    celltype-like:", cand_c[:12])
        if cand_xy:
            print("    coord-like  :", cand_xy[:12])
        with pd.option_context("display.width", 200, "display.max_columns", 14):
            print(head.iloc[:3, :12].to_string())
        # value counts for the most promising annotation columns
        for c in (cand_r + cand_c)[:6]:
            try:
                col = pd.read_csv(p, sep=("\t" if p.endswith((".txt", ".tsv")) else ","),
                                  usecols=[c], low_memory=False)[c]
                vc = col.value_counts().head(14)
                if 1 < len(col.unique()) <= 200:
                    print(f"    values of {c!r}:")
                    for k2, v2 in vc.items():
                        print(f"        {str(k2)[:60]:62s} {v2}")
            except Exception:
                pass

json.dump({k: v for k, v in found.items()}, open(os.path.join(WORK, "found.json"), "w"), indent=1)
print(f"\n\nwork dir: {WORK}")
