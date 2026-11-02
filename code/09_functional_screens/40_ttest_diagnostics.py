#!/usr/bin/env python3
"""
40 - How the t-test p-values in Figures 20-22 are derived, and what they cost.

This script does not change any result; it documents and audits the test, and evaluates
the one alternative worth considering.

THE NORMALISATION, AND WHY NON-TARGETING WELLS SIT ABOVE 100
Every well on a plate is expressed as a percentage of that plate's own non-targeting
median. The non-targeting wells go through the same division, so about half of them land
above 100 by construction - in this screen exactly 4 of 8 on every plate. A value of 104
is not an anomaly to be removed, it IS the control noise, and it is the only direct
estimate of how much a well can move for no biological reason.

THE TEST
For cell line L and gene g:
    x = the gene's 3 replicate values, as % of each plate's NT median
    y = that line's 24 non-targeting wells (8 per plate x 3 plates), same scale
    Welch two-sample t-test of x against y:
        t  = (mean(x) - mean(y)) / sqrt(s_x^2/3 + s_y^2/24)
        df = (s_x^2/3 + s_y^2/24)^2 / [ (s_x^2/3)^2/2 + (s_y^2/24)^2/23 ]
Welch rather than Student because the two samples have neither the same size nor any
reason to share a variance. The effect reported in the figures is mean(x) - mean(y).

THE COST
Because 24 control wells are averaged and only 3 replicates are, s_y^2/24 is tiny and
s_x^2/3 dominates: the gene's three replicates contribute a median 96 % of the squared
standard error, so the Welch df collapses to about 2. Four fifths of gene x line
combinations have df < 2.5. The 24 non-targeting wells set the reference level and
contribute almost nothing to the power - the test is effectively a one-sample t-test on
three numbers, and it can miss very large effects: RPL5 loses 39 percentage points of
viability in HCT116 and returns P = 0.063.

THE ALTERNATIVE
Moderated t (Smyth 2004, the test behind limma), which keeps the same statistic but
replaces each gene's variance with a posterior that borrows from all 116 genes on the
panel:
    s~^2 = (d0 * s0^2 + 2 * s_g^2) / (d0 + 2),  df = d0 + 2
d0 and s0 are estimated from the spread of the observed variances. On this panel d0 is
about 17 in HCT116 and RKO, so each gene gains roughly 17 degrees of freedom, and genes
whose three replicates happen to agree closely are no longer rewarded for it.
"""
import pandas as pd, numpy as np
from scipy.special import digamma, polygamma
from scipy.optimize import brentq
from scipy.stats import t as tdist

LINES = ["HCT116", "RKO", "SW480"]
N_REP = 3


def normalise(line, C):
    """Every well as a percentage of its own plate's non-targeting median."""
    parts = []
    for _, dd in C[C.cell_line == line].groupby("replicate"):
        nt = dd.loc[dd.kind == "negative", "value"]
        e = dd.copy()
        e["pct"] = 100 * e["value"] / nt.median()
        parts.append(e)
    return pd.concat(parts, ignore_index=True)


def fit_prior(s2, d):
    """Smyth's moment estimator for the prior degrees of freedom and variance."""
    z = np.log(s2)
    e = z - digamma(d / 2) + np.log(d / 2)
    ebar = e.mean()
    target = ((e - ebar) ** 2).sum() / (len(e) - 1) - polygamma(1, d / 2)
    if target <= 0:                       # boundary: shrink to one common variance
        d0 = 1e4
    else:
        d0 = 2 * brentq(lambda x: polygamma(1, x / 2) - target, 1e-3, 1e4)
    return d0, float(np.exp(ebar + digamma(d0 / 2) - np.log(d0 / 2)))


C = pd.read_csv("results/func_crispr_long.csv")
C["gene"] = C["gene"].replace({"ACAS2": "ACSS2"})
risk = set(pd.read_csv("results/risk_genes_harmonized.csv")["gene"])

print("=" * 78)
print("1. THE CONTROL WELLS")
print("=" * 78)
for line in LINES:
    Z = normalise(line, C)
    for rep, dd in Z[Z.kind == "negative"].groupby("replicate"):
        p = dd["pct"]
        print(f"{line} {rep}: {len(p)} NT wells, {int((p > 100).sum())} of them above 100, "
              f"range {p.min():.1f}-{p.max():.1f}, sd {p.std(ddof=1):.2f} pp")

rows = []
print("\n" + "=" * 78)
print("2. THE TEST, AND WHERE ITS VARIANCE COMES FROM")
print("=" * 78)
for line in LINES:
    Z = normalise(line, C)
    y = Z.loc[Z.kind == "negative", "pct"].to_numpy(float)
    vy = y.var(ddof=1) / len(y)
    G = Z[Z.kind == "sample"].groupby("gene")["pct"]
    genes = np.array(G.mean().index)
    mean_x = G.mean().to_numpy()
    s_x = G.std(ddof=1).to_numpy()
    vx = s_x ** 2 / N_REP
    se = np.sqrt(vx + vy)
    t = (mean_x - y.mean()) / se
    df = se ** 4 / (vx ** 2 / (N_REP - 1) + vy ** 2 / (len(y) - 1))
    p_welch = 2 * tdist.sf(np.abs(t), df)

    d0, s02 = fit_prior(s_x ** 2, N_REP - 1)
    s2p = (d0 * s02 + (N_REP - 1) * s_x ** 2) / (d0 + N_REP - 1)
    se_m = np.sqrt(s2p / N_REP + vy)
    p_mod = 2 * tdist.sf(np.abs((mean_x - y.mean()) / se_m), d0 + N_REP - 1)

    print(f"\n{line}: NT mean {y.mean():.2f} pp, NT sd {y.std(ddof=1):.2f} pp "
          f"(contributes {vy:.2f} to the squared SE)")
    print(f"  the 3 replicates contribute a median {np.median(vx / (vx + vy)):.1%} of it; "
          f"median Welch df {np.median(df):.2f}; df < 2.5 for {(df < 2.5).mean():.0%} of genes")
    print(f"  moderated prior: d0 = {d0:.1f} extra df, s0 = {np.sqrt(s02):.2f} pp")
    print(f"  significant at P < 0.05, whole panel: Welch {int((p_welch < 0.05).sum())}, "
          f"moderated {int((p_mod < 0.05).sum())} of {len(genes)}")
    rows.append(pd.DataFrame({"cell_line": line, "gene": genes, "delta_pct": mean_x - y.mean(),
                              "sd_reps": s_x, "welch_df": df, "p_welch": p_welch,
                              "p_moderated": p_mod, "is_risk": np.isin(genes, list(risk))}))

R = pd.concat(rows, ignore_index=True)
R.to_csv("results/ttest_diagnostics.csv", index=False)

print("\n" + "=" * 78)
print("3. WHAT THE LOW DEGREES OF FREEDOM COST")
print("=" * 78)
miss = R[(R.delta_pct < -15) & (R.p_welch >= 0.05)].sort_values("delta_pct")
print("large drops the Welch t-test calls non-significant:")
print(miss[["cell_line", "gene", "delta_pct", "sd_reps", "welch_df", "p_welch",
            "p_moderated"]].round(3).to_string(index=False))
small = R[(R.delta_pct.abs() < 7) & (R.p_welch < 0.05)].sort_values("p_welch")
print(f"\nsmall changes it calls significant ({len(small)} of them), the mirror-image error:")
print(small.head(8)[["cell_line", "gene", "delta_pct", "sd_reps", "welch_df", "p_welch",
                     "p_moderated"]].round(3).to_string(index=False))

rr = R[R.is_risk]
print("\n" + "=" * 78)
print("4. EFFECT ON THE 36 RISK GENES IN THE FIGURES")
print("=" * 78)
print(rr.groupby("cell_line").apply(
    lambda d: pd.Series({"Welch P<0.05": int((d.p_welch < 0.05).sum()),
                         "moderated P<0.05": int((d.p_moderated < 0.05).sum())}),
    include_groups=False).to_string())
print(f"\ntotal calls: Welch {int((rr.p_welch < 0.05).sum())}, "
      f"moderated {int((rr.p_moderated < 0.05).sum())} of {len(rr)}")
print("\nwrote results/ttest_diagnostics.csv")
