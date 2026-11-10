#!/usr/bin/env python3
"""
Extract everything needed from the GeoMx archive. Runs on the user's machine
     (pandas + numpy only; no scipy there).

Outputs include:

  geomx_DE_compact.tsv      every gene x compartment x stage contrast from the paper's
                            linear mixed model, with the contrast direction harmonised
                            so that a positive log2 fold change always means HIGHER IN
                            THE LATER STAGE.
  de_exp1.tsv.gz and de_exp2_{epi,str}.tsv.gz
                            staged contrasts for 02_analyze_geomx_support.py.
  geomx_expr_<exp>.tsv      Q3-normalised log2 expression per segment, rebuilt from the
                            raw DCC probe counts and the PKC panel file, for the CRC
                            risk genes plus a matched random background.

The DE tables the paper ships are used as-is rather than refitted: they come from
log2(norm_exp) ~ tissueRegion + (1 + tissueRegion | Patient_ID), which is the right
model (patient random effect with a random slope for region, because segments are
subsampled within patients) and it is the model the published results rest on.
"""
import os, sys, json, re, gzip
import numpy as np
import pandas as pd

ROOT = os.path.expanduser(os.environ.get("SCTWAS_GEOMX_ROOT", "~/spatial_work/GeoMx/Nanostring_GeoMx"))
OUT = os.path.expanduser(os.environ.get("SCTWAS_GEOMX_OUT", "~/spatial_work/out"))
os.makedirs(OUT, exist_ok=True)
RISK = os.path.expanduser(os.environ.get("SCTWAS_RISK_GENES", "~/mnt/AI-gene-drug/results/risk_genes_harmonized.csv"))


def write_staged_contrasts(de, out):
    """Write the analysis inputs at the precision of the archived split tables."""
    exp1 = de.loc[de["experiment"] == "exp1"]
    exp2 = de.loc[de["experiment"] == "exp2"]
    panels = {
        "de_exp1.tsv.gz": exp1,
        "de_exp2_epi.tsv.gz": exp2.loc[exp2["compartment"] == "PanCK+"],
        "de_exp2_str.tsv.gz": exp2.loc[exp2["compartment"] == "PanCK-"],
    }
    for name, panel in panels.items():
        if panel.empty:
            raise ValueError(f"No GeoMx contrast rows for {name}")
    for name, panel in panels.items():
        staged = panel[["gene", "compartment", "transition", "log2FC", "p", "FDR"]].copy()
        staged["log2FC"] = staged["log2FC"].round(4)
        for column in ("p", "FDR"):
            staged[column] = staged[column].map("{:.5g}".format, na_action="ignore")
        staged.to_csv(os.path.join(out, name), sep="\t", index=False)


# ---------------------------------------------------------------- 1. DE tables
DE_FILES = {
    "exp2": f"{ROOT}/exp2/results/DE_analysis/geomx_exp2_LMM_DE_segmentRegion_contrast_results_T1_only.tsv",
    "exp1": f"{ROOT}/exp1/results/DE_analysis/geomx_exp1_LMM_DE_segmentRegion_contrast_results.tsv",
}
STAGE_ORDER = {"NOR": 0, "ADE": 1, "CORE": 2, "INV": 3}

rows = []
for exp, path in DE_FILES.items():
    d = pd.read_csv(path, sep="\t")
    d.columns = [c.strip() for c in d.columns]
    d = d.rename(columns={"Pr(>|t|)": "p"})
    print(f"{exp}: {d.shape[0]} rows, {d['Gene'].nunique()} genes, "
          f"{d['Contrast'].nunique()} contrasts, subsets {sorted(d['Subset'].unique())}")

    pat = re.compile(r"^(\S+?)_([A-Z]+)\s*-\s*(\S+?)_([A-Z]+)$")
    a_stage, b_stage, ok = [], [], []
    for c in d["Contrast"]:
        m = pat.match(str(c).strip())
        if not m:
            a_stage.append(None); b_stage.append(None); ok.append(False); continue
        a_stage.append(m.group(2)); b_stage.append(m.group(4)); ok.append(True)
    d["stage_a"], d["stage_b"], d["parsed"] = a_stage, b_stage, ok
    d = d[d["parsed"] & d["stage_a"].isin(STAGE_ORDER) & d["stage_b"].isin(STAGE_ORDER)].copy()

    # harmonise: log2FC positive = higher in the LATER stage
    later = [a if STAGE_ORDER[a] > STAGE_ORDER[b] else b
             for a, b in zip(d["stage_a"], d["stage_b"])]
    earlier = [b if STAGE_ORDER[a] > STAGE_ORDER[b] else a
               for a, b in zip(d["stage_a"], d["stage_b"])]
    flip = np.where(np.array([STAGE_ORDER[a] for a in d["stage_a"]]) >
                    np.array([STAGE_ORDER[b] for b in d["stage_b"]]), 1.0, -1.0)
    d["from_stage"] = earlier
    d["to_stage"] = later
    d["log2FC"] = d["Estimate"].astype(float) * flip
    d["transition"] = d["from_stage"] + "->" + d["to_stage"]
    d["experiment"] = exp
    d["compartment"] = d["Subset"]
    rows.append(d[["experiment", "Gene", "compartment", "transition", "from_stage",
                   "to_stage", "log2FC", "p", "FDR"]])

DE = pd.concat(rows, ignore_index=True).rename(columns={"Gene": "gene"})
DE.to_csv(f"{OUT}/geomx_DE_compact.tsv", sep="\t", index=False)
write_staged_contrasts(DE, OUT)
print(f"\nwrote geomx_DE_compact.tsv: {DE.shape[0]} rows")
print(DE.groupby(["experiment", "compartment", "transition"]).size().to_string())

# ---------------------------------------------------------------- 2. annotations
for exp in ["exp1", "exp2"]:
    x = f"{ROOT}/{exp}/flatfiles/annotation/geomx_{exp}_sample_annotations.xlsx"
    if os.path.exists(x):
        try:
            sheets = pd.read_excel(x, sheet_name=None)
            for nm, df in sheets.items():
                safe = re.sub(r"\W+", "_", nm)
                df.to_csv(f"{OUT}/{exp}_annot_{safe}.tsv", sep="\t", index=False)
                print(f"{exp} annotation sheet {nm!r}: {df.shape}")
        except Exception as e:
            print(f"!! could not read {x}: {e}")
    p = f"{ROOT}/{exp}/flatfiles/annotation/geomx_{exp}_patient_metadata.tsv"
    if os.path.exists(p):
        pd.read_csv(p, sep="\t").to_csv(f"{OUT}/{exp}_patient_metadata.tsv", sep="\t", index=False)

cms = f"{ROOT}/exp2/results/CMS_prediction/geomx_exp2_cms.tsv"
if os.path.exists(cms):
    pd.read_csv(cms, sep="\t").to_csv(f"{OUT}/exp2_cms.tsv", sep="\t", index=False)

# ---------------------------------------------------------------- 3. expression
risk = pd.read_csv(RISK)
risk_genes = sorted(set(risk["gene"].astype(str)))
print(f"\nrisk genes: {len(risk_genes)}")


def read_pkc(path):
    """RTS probe id -> gene symbol, from the NanoString panel file."""
    with open(path) as fh:
        pk = json.load(fh)
    m = {}
    targets = []
    for grp in pk.get("Targets", []):
        name = grp.get("DisplayName") or grp.get("AnalyteType") or grp.get("TargetName")
        targets.append(grp)
        for pr in grp.get("Probes", []):
            rts = pr.get("RTS_ID") or pr.get("RTS_Seq_ID")
            if rts:
                m[rts] = name
    return m


def read_dcc(path):
    counts, sid = {}, None
    section = None
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith("<") and line.endswith(">"):
                section = line.strip("<>/")
                continue
            if section == "Scan_Attributes" and line.startswith("ID,"):
                sid = line.split(",", 1)[1].strip()
            elif section == "Code_Summary" and "," in line:
                k, v = line.split(",", 1)
                try:
                    counts[k] = int(float(v))
                except ValueError:
                    pass
    return sid, counts


for exp in ["exp2", "exp1"]:
    pkcs = [f"{ROOT}/{exp}/flatfiles/pkcs/{f}"
            for f in os.listdir(f"{ROOT}/{exp}/flatfiles/pkcs")]
    probe2gene = {}
    for p in pkcs:
        probe2gene.update(read_pkc(p))
    print(f"\n{exp}: {len(probe2gene)} probes mapped to {len(set(probe2gene.values()))} targets")

    dcc_dir = f"{ROOT}/{exp}/flatfiles/dccs"
    files = sorted(f for f in os.listdir(dcc_dir) if f.endswith(".dcc"))
    seg_counts, seg_ids = [], []
    for i, f in enumerate(files):
        sid, cnt = read_dcc(os.path.join(dcc_dir, f))
        if sid is None or len(cnt) < 100:
            continue
        g = {}
        for rts, c in cnt.items():
            gn = probe2gene.get(rts)
            if gn is None:
                continue
            g[gn] = g.get(gn, 0) + c
        seg_ids.append(sid); seg_counts.append(g)
        if (i + 1) % 60 == 0:
            print(f"   {i+1}/{len(files)} dcc parsed")
    M = pd.DataFrame(seg_counts, index=seg_ids).fillna(0.0)
    print(f"{exp}: raw matrix {M.shape} (segments x targets)")

    neg = [c for c in M.columns if isinstance(c, str) and
           ("NegProbe" in c or c.startswith("Neg"))]
    q3 = M.drop(columns=neg, errors="ignore").apply(
        lambda r: np.percentile(r[r > 0], 75) if (r > 0).any() else np.nan, axis=1)
    norm = M.div(q3.replace(0, np.nan), axis=0)
    logn = np.log2(norm + 1)
    print(f"{exp}: Q3 median {np.nanmedian(q3):.1f}; {len(neg)} negative probes excluded from Q3")

    keep = [g for g in logn.columns if g in set(risk_genes)]
    rng = np.random.default_rng(19)
    others = [g for g in logn.columns if g not in set(risk_genes) and g not in neg]
    bg = list(rng.choice(others, size=min(1500, len(others)), replace=False))
    sub = logn[keep + bg].copy()
    sub.index.name = "segment_id"
    sub.to_csv(f"{OUT}/geomx_expr_{exp}.tsv", sep="\t")
    json.dump({"risk_genes_present": keep, "background_genes": bg,
               "n_segments": int(sub.shape[0]), "n_targets_total": int(logn.shape[1]),
               "q3_median": float(np.nanmedian(q3))},
              open(f"{OUT}/geomx_expr_{exp}_meta.json", "w"), indent=1)
    print(f"{exp}: wrote expression for {len(keep)} risk genes + {len(bg)} background, "
          f"{sub.shape[0]} segments")

print("\nfiles in", OUT)
for f in sorted(os.listdir(OUT)):
    print(f"  {os.path.getsize(os.path.join(OUT, f))/1e6:8.2f} MB  {f}")
