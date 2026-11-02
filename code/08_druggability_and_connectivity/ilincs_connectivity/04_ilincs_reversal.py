#!/usr/bin/env python
"""Do iLINCS perturbagen signatures reverse the cell-type TWAS Z profile?

Design
------
For each cell-type TWAS model we build a 100-gene query signature (the 50 most
positive and 50 most negative TWAS Z scores) and correlate it, gene by gene, with
every iLINCS perturbagen signature (log2 differential expression, large intestine
LIB5 library).  Pearson r is computed over the 100 query genes.

  r < 0  the perturbagen pushes risk-increasing genes down and risk-decreasing
         genes up - a REVERSAL of the genetically predicted risk profile, the
         therapeutically interesting direction
  r > 0  the perturbagen mimics the risk profile

Two significance assessments per signature:
  * parametric Pearson P (n = 100), Benjamini-Hochberg across the 16,767 signatures
  * empirical P against 1,000 random 100-gene query sets drawn from the same
    TWAS-by-landmark gene pool, each carrying its own real TWAS Z values - this is
    the "better than 100 random genes" comparison

IMPORTANT constraint: the iLINCS LIB5 matrix covers the 978 LINCS landmark genes
only, of which ~630 have a TWAS Z in a given model.  The top 50 / bottom 50 are
therefore taken among the shared landmark genes, not among all ~16,000 tested
genes - otherwise only a handful of the top 100 genome-wide would be measurable.

Outputs (under <ROOT>/iLINCS_reversal):
  table/ilincs_reversal_all_signatures.csv   every model x signature r, P, FDR, empirical P
  table/ilincs_reversal_top_signatures.csv   FDR<0.05 reversing signatures
  table/ilincs_reversal_model_summary.csv    per-model summary
  table/ilincs_query_genes.csv               the 100 query genes per model
"""
import os
import numpy as np
import pandas as pd
from scipy.stats import t as tdist
from statsmodels.stats.multitest import multipletests

TWAS_DIR = r"W:\VUMC_projects\scTWAS\TWAS\SPrediXcan_peakSNP_association"
ILINCS = r"L:\Drugs_Proteins\iLINCS\LIB5_signature_expr_logDiffExpr_large_intestine.csv"
OUT = r"W:\VUMC_projects\scTWAS\iLINCS_reversal"
N_SIDE, N_PERM, SEED = 50, 1000, 0

os.makedirs(rf"{OUT}\table", exist_ok=True)
os.makedirs(rf"{OUT}\figure", exist_ok=True)
os.makedirs(rf"{OUT}\src", exist_ok=True)


def load_ilincs(path):
    hdr = open(path).readline().rstrip("\n").split(",")
    de = [c for c in hdr if c.endswith("_DiffExp")]
    m = pd.read_csv(path, usecols=["Name_GeneSymbol"] + de,
                    dtype={c: "float32" for c in de}, engine="c").set_index("Name_GeneSymbol")
    m.columns = [c.replace("_DiffExp", "") for c in m.columns]
    return m


def corr_against(z, S):
    """Pearson r between vector z (n,) and every column of S (n x m)."""
    zc = z - z.mean()
    num = zc @ S                                   # S columns are centred below
    n = S.shape[0]
    cm = S.sum(0)
    sq = (S * S).sum(0)
    den = np.sqrt((zc * zc).sum()) * np.sqrt(np.maximum(sq - cm * cm / n, 1e-12))
    return (num - zc.sum() * cm / n) / den


def main():
    il = load_ilincs(ILINCS)
    sig_ids = np.array(il.columns)
    files = sorted(f for f in os.listdir(TWAS_DIR) if f.endswith(".EUR.TWAS.csv"))
    rng = np.random.default_rng(SEED)

    all_rows, qrows, summary = [], [], []
    for fn in files:
        model = fn.replace("_peakSNP.EUR.TWAS.csv", "")
        tw = pd.read_csv(os.path.join(TWAS_DIR, fn), usecols=["gene_name", "zscore"]).dropna()
        # one Z per symbol: keep the most extreme association
        tw = (tw.assign(absz=tw.zscore.abs()).sort_values("absz", ascending=False)
                .drop_duplicates("gene_name").drop(columns="absz"))
        shared = tw[tw.gene_name.isin(il.index)].sort_values("zscore")
        if len(shared) < 4 * N_SIDE:
            continue
        pos = shared.tail(N_SIDE)
        neg = shared.head(N_SIDE)
        q = pd.concat([neg, pos])
        genes = q.gene_name.to_numpy()
        z = q.zscore.to_numpy(dtype=np.float64)

        S = il.loc[genes].to_numpy(dtype=np.float32)
        r = corr_against(z, S).astype(np.float64)

        n = len(z)
        tstat = r * np.sqrt((n - 2) / np.maximum(1 - r ** 2, 1e-12))
        pval = 2 * tdist.sf(np.abs(tstat), n - 2)
        fdr = multipletests(pval, method="fdr_bh")[1]

        # empirical null: random 100-gene query sets from the same shared pool
        pool_g = shared.gene_name.to_numpy()
        pool_z = shared.zscore.to_numpy(dtype=np.float64)
        Spool = il.loc[pool_g].to_numpy(dtype=np.float32)
        ge = np.zeros(len(r), dtype=np.int32)      # |null r| >= |observed r|
        le = np.zeros(len(r), dtype=np.int32)      # null r <= observed r (reversal tail)
        for _ in range(N_PERM):
            idx = rng.choice(len(pool_g), size=n, replace=False)
            rn = corr_against(pool_z[idx], Spool[idx])
            ge += (np.abs(rn) >= np.abs(r))
            le += (rn <= r)
        emp_two = (ge + 1) / (N_PERM + 1)
        emp_rev = (le + 1) / (N_PERM + 1)

        all_rows.append(pd.DataFrame({
            "model": model, "signature": sig_ids, "pearson_r": r, "P_parametric": pval,
            "FDR_parametric": fdr, "P_empirical_2sided": emp_two, "P_empirical_reversal": emp_rev}))
        qrows.append(pd.DataFrame({"model": model, "gene_name": genes, "zscore": z,
                                   "side": ["negative"] * N_SIDE + ["positive"] * N_SIDE}))
        rev = (r < 0) & (fdr < 0.05)
        summary.append({
            "model": model, "n_genes_tested": len(tw), "n_landmark_shared": len(shared),
            "z_range_query": f"{z.min():.2f} to {z.max():.2f}",
            "n_signatures": len(r), "r_min": r.min(), "r_max": r.max(),
            "n_reversing_FDR05": int(rev.sum()),
            "n_reversing_FDR05_emp05": int((rev & (emp_rev < 0.05)).sum()),
            "n_mimicking_FDR05": int(((r > 0) & (fdr < 0.05)).sum()),
            "best_reversal_signature": sig_ids[r.argmin()], "best_reversal_r": r.min(),
            "best_reversal_emp_P": emp_rev[r.argmin()]})
        print(f"{model:34s} shared {len(shared):4d}  r [{r.min():+.3f},{r.max():+.3f}]  "
              f"reversing FDR<0.05 {int(rev.sum()):5d}")

    allr = pd.concat(all_rows, ignore_index=True)
    qdf = pd.concat(qrows, ignore_index=True)
    summ = pd.DataFrame(summary)
    top = (allr[(allr.pearson_r < 0) & (allr.FDR_parametric < 0.05)]
           .sort_values(["model", "pearson_r"]))
    allr.to_csv(rf"{OUT}\table\ilincs_reversal_all_signatures.csv", index=False)
    top.to_csv(rf"{OUT}\table\ilincs_reversal_top_signatures.csv", index=False)
    summ.to_csv(rf"{OUT}\table\ilincs_reversal_model_summary.csv", index=False)
    qdf.to_csv(rf"{OUT}\table\ilincs_query_genes.csv", index=False)
    print("\n", summ.to_string(index=False))
    return allr, qdf, summ, top


if __name__ == "__main__":
    allr, qdf, summ, top = main()
