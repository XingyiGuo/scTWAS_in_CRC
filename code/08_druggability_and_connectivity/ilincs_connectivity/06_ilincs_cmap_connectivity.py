#!/usr/bin/env python
"""CMap weighted connectivity (WTCS/NCS) between cell-type TWAS profiles and iLINCS.

Replaces the Pearson correlation of script 04/05 with the connectivity score of
Subramanian et al. 2017 (Cell 171:1437), the standard CMap statistic:

  * the query is a tagged gene pair - up = the 50 most positive TWAS Z, down = the
    50 most negative, taken among genes shared with the iLINCS landmark space;
  * each perturbagen signature is ranked by log2FC, descending;
  * ES_up and ES_down are weighted Kolmogorov-Smirnov enrichment scores (GSEA
    running sum, weight = |log2FC|, p = 1) of each tag set in that ranked list;
  * WTCS = (ES_up - ES_down) / 2 when the two have opposite signs, else 0.
    WTCS < 0 means the query's up-genes are suppressed and its down-genes induced
    by the perturbagen - a REVERSAL of the genetically predicted risk profile;
  * NCS normalises WTCS within cell line by the mean of the same-sign scores in
    that cell line, so signatures from different lines are comparable.

Inference: NCS and tau_like rank signatures WITHIN a query, which is descriptive
only - the library cannot serve as its own null for the question "is the best of
16,767 unusual?", because the library is the data (its score distribution is also
bimodal with ~50% exact zeros, so a Gaussian z against it is doubly wrong). The
significance test is therefore the same max-statistic permutation as script 05,
at 250 draws: for each random tagged query the most negative WTCS over the whole
library is recorded, and the observed best reversal is read against that.

Outputs (under <ROOT>/iLINCS_reversal):
  table/ilincs_cmap_all_signatures.csv      model x signature WTCS, NCS, ES_up, ES_down
  table/ilincs_cmap_top20_per_model.csv     20 most reversing signatures per model
  table/ilincs_cmap_model_summary.csv       best reversal and library-null stats per model
  table/ilincs_cmap_vs_pearson.csv          agreement between the two statistics
"""
import os
import numpy as np
import pandas as pd
from scipy.stats import norm, rankdata
from statsmodels.stats.multitest import multipletests

import importlib.util
spec = importlib.util.spec_from_file_location("s04", "04_ilincs_reversal.py")
s04 = importlib.util.module_from_spec(spec)
spec.loader.exec_module(s04)

TWAS_DIR, OUT = s04.TWAS_DIR, s04.OUT
META = r"L:\Drugs_Proteins\iLINCS\meta_sigSearch_LIB_5_143374_20231227.csv"
N_SIDE, N_PERM, SEED = 50, 250, 2


def enrichment(rankpos, W, n_genes):
    """Weighted KS enrichment score of a tag set in each ranked signature.

    rankpos : (n_tag, n_sig) int   0-based position of each tag gene in the
                                   descending-|log2FC| ordering of its signature
    W       : (n_tag, n_sig) float weights, |log2FC| of those genes
    Returns (n_sig,) ES, the running-sum deviation of largest magnitude.
    """
    n_tag = rankpos.shape[0]
    o = np.argsort(rankpos, axis=0, kind="stable")           # hits in rank order
    p = np.take_along_axis(rankpos, o, axis=0).astype(np.float64)
    w = np.take_along_axis(W, o, axis=0).astype(np.float64)
    NR = w.sum(0)
    NR[NR == 0] = 1e-12
    cum = np.cumsum(w, axis=0) / NR                          # hit contribution
    i = np.arange(1, n_tag + 1)[:, None]
    miss_after = (p + 1 - i) / (n_genes - n_tag)             # after taking the hit
    miss_before = (p - (i - 1)) / (n_genes - n_tag)
    top = cum - miss_after                                   # running sum at hits
    bot = np.vstack([np.zeros((1, cum.shape[1])), cum[:-1]]) - miss_before
    hi, lo = top.max(0), bot.min(0)
    return np.where(hi > -lo, hi, lo)


def wtcs(up_idx, dn_idx, X, rankpos, n_genes):
    es_u = enrichment(rankpos[up_idx], np.abs(X[up_idx]), n_genes)
    es_d = enrichment(rankpos[dn_idx], np.abs(X[dn_idx]), n_genes)
    out = np.where(np.sign(es_u) != np.sign(es_d), (es_u - es_d) / 2.0, 0.0)
    return out, es_u, es_d


def main():
    il = s04.load_ilincs(s04.ILINCS)
    sig_ids = np.array(il.columns)
    X = il.to_numpy(dtype=np.float32)
    n_genes = X.shape[0]
    order = np.argsort(-X, axis=0, kind="stable")
    rankpos = np.empty_like(order)
    np.put_along_axis(rankpos, order, np.arange(n_genes)[:, None], axis=0)

    rng = np.random.default_rng(SEED)
    meta = pd.read_csv(META)
    cell = (pd.Series(sig_ids).map(meta.set_index("SignatureId").CellLine)
            .fillna("unknown").to_numpy())
    qdf = pd.read_csv(rf"{OUT}\table\ilincs_query_genes.csv")
    pear = pd.read_csv(rf"{OUT}\table\ilincs_reversal_all_signatures.csv",
                       usecols=["model", "signature", "pearson_r"])
    gpos = {g: i for i, g in enumerate(il.index)}

    rows, nulls = [], []
    for model, q in qdf.groupby("model"):
        up = np.array([gpos[g] for g in q.loc[q.side.eq("positive"), "gene_name"]])
        dn = np.array([gpos[g] for g in q.loc[q.side.eq("negative"), "gene_name"]])
        w, esu, esd = wtcs(up, dn, X, rankpos, n_genes)

        # NCS: divide by the mean of same-sign scores within each cell line
        ncs = np.zeros_like(w)
        for cl in np.unique(cell):
            m = cell == cl
            for s in (1, -1):
                sel = m & (np.sign(w) == s)
                if sel.sum():
                    mu = np.abs(w[sel]).mean()
                    ncs[sel] = w[sel] / (mu if mu > 0 else 1)

        # Inference without permutation, as CLUE does it: read each signature against
        # the distribution of scores the same query produces over the whole library.
        # tau_like = signed percentile of |NCS| among the same-sign scores (0-100).
        tau = np.zeros_like(ncs)
        for s in (1, -1):
            sel = np.sign(ncs) == s
            if sel.sum():
                tau[sel] = s * 100 * rankdata(np.abs(ncs[sel])) / sel.sum()

        # Valid null for "is the BEST of 16,767 unusual?": random tagged queries.
        # The library's own distribution cannot answer this - the library is the data.
        tw = pd.read_csv(os.path.join(TWAS_DIR, f"{model}_peakSNP.EUR.TWAS.csv"),
                         usecols=["gene_name", "zscore"]).dropna()
        tw = (tw.assign(a=tw.zscore.abs()).sort_values("a", ascending=False)
                .drop_duplicates("gene_name"))
        pool = np.array([gpos[g] for g in tw.gene_name if g in gpos])
        null_min = np.empty(N_PERM)
        for b in range(N_PERM):
            pick = rng.choice(pool, size=2 * N_SIDE, replace=False)
            null_min[b] = wtcs(pick[:N_SIDE], pick[N_SIDE:], X, rankpos, n_genes)[0].min()

        obs = w.min()
        p_max = (int((null_min <= obs).sum()) + 1) / (N_PERM + 1)
        nulls.append({"model": model, "observed_best_WTCS": obs,
                      "null_min_WTCS_mean": null_min.mean(),
                      "null_min_WTCS_sd": null_min.std(), "P_maxstat": p_max,
                      "n_perm": N_PERM,
                      "best_signature": sig_ids[w.argmin()],
                      "best_NCS": ncs[w.argmin()], "best_tau": tau[w.argmin()],
                      "library_NCS_p01": np.percentile(ncs[ncs != 0], 1),
                      "n_WTCS_lt_0": int((w < 0).sum()),
                      "frac_zero_WTCS": float((w == 0).mean())})
        rows.append(pd.DataFrame({"model": model, "signature": sig_ids, "ES_up": esu,
                                  "ES_down": esd, "WTCS": w, "NCS": ncs, "tau_like": tau, "CellLine": cell}))
        print(f"{model:30s} best WTCS {obs:+.3f} NCS {ncs[w.argmin()]:+.2f} | null {null_min.mean():+.3f}"
              f" +/- {null_min.std():.3f} | P {p_max:.3f} | zero-WTCS {100*(w==0).mean():.0f}%")

    allc = pd.concat(rows, ignore_index=True).merge(pear, on=["model", "signature"], how="left")
    allc = allc.merge(meta[["SignatureId", "Perturbagen", "GeneTargets", "Concentration", "Time"]],
                      left_on="signature", right_on="SignatureId", how="left").drop(columns="SignatureId")
    nl = pd.DataFrame(nulls).merge(
        meta[["SignatureId", "Perturbagen", "GeneTargets", "CellLine"]],
        left_on="best_signature", right_on="SignatureId", how="left").drop(columns="SignatureId")
    top = (allc.sort_values(["model", "WTCS"]).groupby("model").head(20)
              .loc[:, ["model", "signature", "Perturbagen", "GeneTargets", "CellLine",
                       "Concentration", "Time", "WTCS", "NCS", "tau_like", "ES_up", "ES_down", "pearson_r"]])
    agree = (allc.groupby("model")
             .apply(lambda g: pd.Series({
                 "pearson_vs_WTCS_r": g.WTCS.corr(g.pearson_r),
                 "spearman_of_ranks": g.WTCS.corr(g.pearson_r, method="spearman"),
                 "overlap_top20": len(set(g.nsmallest(20, "WTCS").signature) &
                                      set(g.nsmallest(20, "pearson_r").signature))}),
                    include_groups=False).reset_index())

    allc.to_csv(rf"{OUT}\table\ilincs_cmap_all_signatures.csv", index=False)
    top.to_csv(rf"{OUT}\table\ilincs_cmap_top20_per_model.csv", index=False)
    nl.to_csv(rf"{OUT}\table\ilincs_cmap_model_summary.csv", index=False)
    agree.to_csv(rf"{OUT}\table\ilincs_cmap_vs_pearson.csv", index=False)
    print("\n", nl[["model", "observed_best_WTCS", "best_NCS", "P_maxstat", "null_min_WTCS_mean",
                    "Perturbagen", "CellLine"]].to_string(index=False))
    print("\n", agree.to_string(index=False))
    return allc, top, nl, agree


if __name__ == "__main__":
    allc, top, nl, agree = main()
