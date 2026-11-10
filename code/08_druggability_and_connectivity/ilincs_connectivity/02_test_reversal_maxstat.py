#!/usr/bin/env python
"""Corrected inference for the iLINCS reversal scan, plus perturbagen names.

Script 01 asked, per signature, "is this r more extreme than r for random 100-gene
queries?".  That question is fine for a pre-specified signature but not for the
signature we report, which is the minimum over 16,767 - a selected extreme.  This
script builds the MAX-STATISTIC null: for each of 1,000 random 100-gene queries,
scan all signatures and keep the single most negative r.  The observed best
reversal is then compared with that distribution, which is the null for "the best
reversing perturbagen for this model".

Also adds: Spearman r as a rank-based sensitivity check (Pearson over 100 genes is
sensitive to the few extreme TWAS Z values), and the perturbagen / cell line /
dose / time metadata from the LIB5 signature metadata CSV.

Outputs (under <ROOT>/iLINCS_reversal):
  table/ilincs_reversal_all_signatures_annotated.csv
  table/ilincs_reversal_maxstat_null.csv
  table/ilincs_reversal_top20_per_model.csv
"""
import os
import numpy as np
import pandas as pd
from scipy.stats import rankdata

import importlib.util
from pathlib import Path

spec = importlib.util.spec_from_file_location(
    "reversal_scan", Path(__file__).with_name("01_scan_ilincs_reversal.py")
)
reversal_scan = importlib.util.module_from_spec(spec)
spec.loader.exec_module(reversal_scan)

TWAS_DIR = reversal_scan.TWAS_DIR
OUT = reversal_scan.OUT
META = r"L:\Drugs_Proteins\iLINCS\meta_sigSearch_LIB_5_143374_20231227.csv"
N_PERM, SEED = 1000, 1


def main():
    il = reversal_scan.load_ilincs(reversal_scan.ILINCS)
    sig_ids = np.array(il.columns)
    S = il.to_numpy(dtype=np.float32)
    allr = pd.read_csv(rf"{OUT}\table\ilincs_reversal_all_signatures.csv")
    qdf = pd.read_csv(rf"{OUT}\table\ilincs_query_genes.csv")
    gene_pos = {g: i for i, g in enumerate(il.index)}

    rng = np.random.default_rng(SEED)
    rows, spearman = [], []
    for model, q in qdf.groupby("model"):
        fn = f"{model}_peakSNP.EUR.TWAS.csv"
        tw = pd.read_csv(os.path.join(TWAS_DIR, fn), usecols=["gene_name", "zscore"]).dropna()
        tw = (tw.assign(absz=tw.zscore.abs()).sort_values("absz", ascending=False)
                .drop_duplicates("gene_name"))
        shared = tw[tw.gene_name.isin(il.index)]
        pool_idx = np.array([gene_pos[g] for g in shared.gene_name])
        pool_z = shared.zscore.to_numpy(dtype=np.float64)

        qi = np.array([gene_pos[g] for g in q.gene_name])
        zq = q.zscore.to_numpy(dtype=np.float64)
        r_obs = allr.loc[allr.model.eq(model), ["signature", "pearson_r"]].set_index("signature")
        r_obs = r_obs.loc[sig_ids, "pearson_r"].to_numpy()

        # Spearman on the same query genes: rank the TWAS Z, and the signature values
        # WITHIN the query gene subset
        rs = reversal_scan.corr_against(
            rankdata(zq), np.apply_along_axis(rankdata, 0, S[qi]).astype(np.float32))
        spearman.append(pd.DataFrame({"model": model, "signature": sig_ids, "spearman_r": rs}))

        null_min, null_maxabs = np.empty(N_PERM), np.empty(N_PERM)
        for b in range(N_PERM):
            pos = rng.choice(len(pool_idx), size=len(qi), replace=False)
            rn = reversal_scan.corr_against(pool_z[pos], S[pool_idx[pos]])
            null_min[b], null_maxabs[b] = rn.min(), np.abs(rn).max()
        obs_min = r_obs.min()
        p_max = (int((null_min <= obs_min).sum()) + 1) / (N_PERM + 1)
        rows.append({"model": model, "observed_best_reversal_r": obs_min,
                     "best_signature": sig_ids[r_obs.argmin()],
                     "null_min_r_mean": null_min.mean(), "null_min_r_sd": null_min.std(),
                     "null_min_r_p05": np.percentile(null_min, 5),
                     "P_maxstat": p_max,
                     "n_null_at_least_as_extreme": int((null_min <= obs_min).sum()),
                     "null_maxabs_mean": null_maxabs.mean()})
        print(f"{model:30s} obs {obs_min:+.3f} | null min r {null_min.mean():+.3f}"
              f" +/- {null_min.std():.3f} | P_maxstat {p_max:.3f}")

    mx = pd.DataFrame(rows)
    sp = pd.concat(spearman, ignore_index=True)
    meta = pd.read_csv(META)
    ann = (allr.merge(sp, on=["model", "signature"], how="left")
                .merge(meta, left_on="signature", right_on="SignatureId", how="left")
                .drop(columns=["SignatureId"]))
    ann.to_csv(rf"{OUT}\table\ilincs_reversal_all_signatures_annotated.csv", index=False)
    mx.to_csv(rf"{OUT}\table\ilincs_reversal_maxstat_null.csv", index=False)
    top = (ann.sort_values(["model", "pearson_r"]).groupby("model").head(20)
              .loc[:, ["model", "signature", "Perturbagen", "GeneTargets", "CellLine",
                       "Concentration", "Time", "pearson_r", "spearman_r", "P_parametric",
                       "FDR_parametric", "P_empirical_reversal"]])
    top.to_csv(rf"{OUT}\table\ilincs_reversal_top20_per_model.csv", index=False)
    print("\n", mx.to_string(index=False))
    print("\nmetadata coverage:", int(ann.Perturbagen.notna().sum()), "of", len(ann))
    return ann, mx, top


if __name__ == "__main__":
    ann, mx, top = main()
