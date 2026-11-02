#!/usr/bin/env python3
"""
32 - The siRNA and CRISPR-Cas9 viability screens: parsing and quality control.

Input is the experimenters' workbook, siRNA_CRISPR-116.xlsx, with two sheets:

  siRNA        333 genes plus controls, three cell lines (HCT116, RKO, SW480) x three
               replicates, each value a per-well viability read expressed as % of
               control. Controls: DEATH (4 wells, positive = kills cells) and
               NT (12 wells, non-targeting negative).
  crisp-cas9   116 genes plus controls, the same three lines x three replicates, laid
               out as three side-by-side blocks. Each replicate carries the cell number
               and a second column that is a plain linear rescaling of it
               (z = (value - 101)/18, i.e. plate-level, not NT-based), so it adds no
               information and is recomputed here against the NT wells instead.

QC is what decides whether a cell line can be used at all:

  Z'-factor       1 - 3(sd_pos + sd_neg) / |mean_pos - mean_neg|. The standard
                  high-throughput screen acceptance metric: > 0.5 excellent,
                  0 - 0.5 marginal but workable, < 0 means the positive and negative
                  control distributions overlap and the plate cannot separate a hit
                  from nothing.
  SSMD            (mean_neg - mean_pos) / sqrt(sd_neg^2 + sd_pos^2), less sensitive to
                  small control counts than Z'.
  NT CV           spread of the negative controls; a noisy baseline inflates every
                  downstream z-score.
  replicate r     Spearman correlation between replicates over the assayed genes.
  cross-line r    Spearman correlation of per-gene effects between cell lines - a line
                  that is internally reproducible but disagrees with the other two is
                  measuring something else.

Outputs
  results/func_sirna_long.csv, func_crispr_long.csv   tidy per-well values
  results/func_qc.csv, func_qc.json                   the metrics and the exclusion call
"""
import pandas as pd, numpy as np, json
from scipy.stats import spearmanr

XL = "func/siRNA_CRISPR-116.xlsx"
LINES = ["HCT116", "RKO", "SW480"]

# ---------------------------------------------------------------- parse siRNA
si = pd.read_excel(XL, "siRNA", header=None)
line_of_col = {}
cur = None
for c in range(1, si.shape[1]):
    v = si.iat[0, c]
    if isinstance(v, str) and v.strip():
        cur = v.strip()
    line_of_col[c] = cur
rep_of_col = {c: str(si.iat[2, c]).strip() for c in range(1, si.shape[1])}

rows = []
for r in range(4, si.shape[0]):
    g = si.iat[r, 0]
    if not isinstance(g, str) or not g.strip():
        continue
    for c in range(1, si.shape[1]):
        v = pd.to_numeric(si.iat[r, c], errors="coerce")
        if np.isfinite(v):
            rows.append({"assay": "siRNA", "gene": g.strip(), "cell_line": line_of_col[c],
                         "replicate": rep_of_col[c], "value": float(v), "well_row": r})
S = pd.DataFrame(rows)
S["kind"] = np.where(S["gene"] == "NT", "negative",
                     np.where(S["gene"] == "DEATH", "positive", "sample"))
print(f"siRNA wells {len(S):,}  genes {S.loc[S.kind=='sample','gene'].nunique()}  "
      f"lines {sorted(S.cell_line.unique())}")

# ---------------------------------------------------------------- parse CRISPR
cr = pd.read_excel(XL, "crisp-cas9", header=None)
blocks = {}
for c in range(cr.shape[1]):
    v = cr.iat[0, c]
    if isinstance(v, str) and v.strip() in LINES:
        blocks[v.strip()] = c
rows = []
for line, c0 in blocks.items():
    # header row 1: "Gene Symbol", then "cell number r1", blank, "cell number r2", ...
    val_cols = [c for c in range(c0 + 1, min(c0 + 8, cr.shape[1]))
                if isinstance(cr.iat[1, c], str) and "cell number" in str(cr.iat[1, c])]
    for r in range(2, cr.shape[0]):
        g = cr.iat[r, c0]
        if not isinstance(g, str) or not g.strip():
            continue
        for k, c in enumerate(val_cols, start=1):
            v = pd.to_numeric(cr.iat[r, c], errors="coerce")
            if np.isfinite(v):
                rows.append({"assay": "CRISPR", "gene": g.strip(), "cell_line": line,
                             "replicate": f"rep{k}", "value": float(v), "well_row": r})
C = pd.DataFrame(rows)
C["kind"] = np.where(C["gene"] == "NT", "negative",
                     np.where(C["gene"] == "LETHAL", "positive", "sample"))
print(f"CRISPR wells {len(C):,}  genes {C.loc[C.kind=='sample','gene'].nunique()}  "
      f"lines {sorted(C.cell_line.unique())}")

S.to_csv("results/func_sirna_long.csv", index=False)
C.to_csv("results/func_crispr_long.csv", index=False)

# ---------------------------------------------------------------- QC
def zprime(pos, neg):
    d = abs(neg.mean() - pos.mean())
    if d == 0:
        return np.nan
    return 1 - 3 * (pos.std(ddof=1) + neg.std(ddof=1)) / d


def ssmd(pos, neg):
    v = pos.var(ddof=1) + neg.var(ddof=1)
    return (neg.mean() - pos.mean()) / np.sqrt(v) if v > 0 else np.nan


qc_rows = []
for D, assay in [(S, "siRNA"), (C, "CRISPR")]:
    for line in LINES:
        for rep in sorted(D.replicate.unique()):
            d = D[(D.cell_line == line) & (D.replicate == rep)]
            pos = d.loc[d.kind == "positive", "value"]
            neg = d.loc[d.kind == "negative", "value"]
            smp = d.loc[d.kind == "sample", "value"]
            qc_rows.append({
                "assay": assay, "cell_line": line, "replicate": rep,
                "n_samples": len(smp), "n_NT": len(neg), "n_pos": len(pos),
                "NT_mean": neg.mean(), "NT_sd": neg.std(ddof=1),
                "NT_cv_pct": 100 * neg.std(ddof=1) / neg.mean(),
                "pos_mean": pos.mean(), "pos_sd": pos.std(ddof=1),
                "signal_window": neg.mean() / pos.mean() if pos.mean() else np.nan,
                "z_prime": zprime(pos, neg), "ssmd": ssmd(pos, neg),
                "sample_median": smp.median(), "sample_mad": (smp - smp.median()).abs().median(),
            })
QC = pd.DataFrame(qc_rows)

# replicate-replicate agreement over the assayed genes
rep_rows = []
for D, assay in [(S, "siRNA"), (C, "CRISPR")]:
    for line in LINES:
        w = (D[(D.cell_line == line) & (D.kind == "sample")]
             .pivot_table(index="gene", columns="replicate", values="value"))
        reps = list(w.columns)
        for i in range(len(reps)):
            for j in range(i + 1, len(reps)):
                ok = w[[reps[i], reps[j]]].dropna()
                rep_rows.append({"assay": assay, "cell_line": line,
                                 "pair": f"{reps[i]}-{reps[j]}", "n": len(ok),
                                 "spearman": spearmanr(ok[reps[i]], ok[reps[j]]).statistic})
RC = pd.DataFrame(rep_rows)
print("\n=== replicate agreement (Spearman over assayed genes) ===")
print(RC.pivot_table(index=["assay", "cell_line"], values="spearman",
                     aggfunc=["mean", "min"]).round(3).to_string())

# cross-cell-line agreement of per-gene effects
cross_rows = []
for D, assay in [(S, "siRNA"), (C, "CRISPR")]:
    g = (D[D.kind == "sample"].groupby(["cell_line", "gene"])["value"].mean()
         .unstack(0))
    for i in range(len(LINES)):
        for j in range(i + 1, len(LINES)):
            a, b = LINES[i], LINES[j]
            ok = g[[a, b]].dropna()
            cross_rows.append({"assay": assay, "pair": f"{a}-{b}", "n_genes": len(ok),
                               "spearman": spearmanr(ok[a], ok[b]).statistic,
                               "pearson": ok[a].corr(ok[b])})
CX = pd.DataFrame(cross_rows)
print("\n=== cross-cell-line agreement of per-gene means ===")
print(CX.round(3).to_string(index=False))

print("\n=== QC per cell line x replicate ===")
print(QC[["assay", "cell_line", "replicate", "NT_mean", "NT_cv_pct", "pos_mean",
          "signal_window", "z_prime", "ssmd"]].round(3).to_string(index=False))

# ---------------------------------------------------------------- exclusion call
summ = (QC.groupby(["assay", "cell_line"])
        .agg(z_prime=("z_prime", "mean"), ssmd=("ssmd", "mean"),
             NT_cv=("NT_cv_pct", "mean"), window=("signal_window", "mean"))
        .reset_index())
rep_mean = RC.groupby(["assay", "cell_line"])["spearman"].mean().rename("rep_r").reset_index()
summ = summ.merge(rep_mean, on=["assay", "cell_line"])
cross_mean = []
for assay in ["siRNA", "CRISPR"]:
    for line in LINES:
        rs = CX[(CX.assay == assay) & CX.pair.str.contains(line)]["spearman"]
        cross_mean.append({"assay": assay, "cell_line": line, "cross_r": rs.mean()})
summ = summ.merge(pd.DataFrame(cross_mean), on=["assay", "cell_line"])
print("\n=== per cell line summary ===")
print(summ.round(3).to_string(index=False))

si_s = summ[summ.assay == "siRNA"].set_index("cell_line")
worst_z = si_s["z_prime"].idxmin()
worst_cross = si_s["cross_r"].idxmin()
worst_win = si_s["window"].idxmin()
votes = pd.Series([worst_z, worst_cross, worst_win]).value_counts()
excluded = votes.index[0] if votes.iloc[0] >= 2 else None
print(f"\nworst Z' = {worst_z}; worst cross-line agreement = {worst_cross}; "
      f"smallest control window = {worst_win}")
print(f"EXCLUDED: {excluded}" if excluded else "\nno line is consistently worst")

QC.to_csv("results/func_qc.csv", index=False)
RC.to_csv("results/func_qc_replicates.csv", index=False)
CX.to_csv("results/func_qc_crossline.csv", index=False)
json.dump({"per_line": summ.to_dict("records"),
           "replicate_agreement": RC.to_dict("records"),
           "cross_line": CX.to_dict("records"),
           "excluded_cell_line": excluded,
           "exclusion_votes": {k: int(v) for k, v in votes.items()},
           "n_genes_sirna": int(S.loc[S.kind == "sample", "gene"].nunique()),
           "n_genes_crispr": int(C.loc[C.kind == "sample", "gene"].nunique()),
           }, open("results/func_qc.json", "w"), indent=1, default=float)
print("\nDone.")
