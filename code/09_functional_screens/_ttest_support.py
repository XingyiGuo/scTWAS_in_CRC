"""
Shared helper for the t-test view of the functional screens (scripts 35 and 37).

Figures 20 and 21 report one test only: a Welch t-test of a gene's three replicate
viability values, each expressed as a percentage of its own plate's non-targeting
median, against that cell line's non-targeting wells. The effect size that goes with
that test is the difference of means in percentage points, so the figures are a plain
differential analysis against the controls - no robust z, no empirical null, no FDR.
Those remain in the tables for anyone who wants them.

On top of that this module answers the question the figures are really for: does the
measured change agree with the direction the TWAS assigned to the gene?

  TWAS risk-increasing  the germline effect raises risk through higher activity, so if
                        that carried into the tumour cell the cell should need the gene
                        and removing it should LOWER viability (delta < 0).
  TWAS protective       the mirror image: removing it should RAISE viability (delta > 0).

A call is made only where the t-test is significant; everywhere else the honest label
is "not significant", not "opposes".
"""
import numpy as np

SUPPORTS, OPPOSES, NS = "supports TWAS", "opposes TWAS", "not significant"
EXPECT = {"risk-increasing": -1.0, "protective": +1.0}


def expected_sign(twas_direction):
    """+1 if knockout should raise viability, -1 if it should lower it, nan if unknown."""
    return EXPECT.get(twas_direction, np.nan)


def support_call(twas_direction, delta_pct, significant):
    """Label one gene x cell line."""
    if not significant:
        return NS
    e = expected_sign(twas_direction)
    if not np.isfinite(e) or not np.isfinite(delta_pct) or delta_pct == 0:
        return NS
    return SUPPORTS if np.sign(delta_pct) == np.sign(e) else OPPOSES


def gene_call(labels):
    """Collapse a gene's per-line labels into one verdict."""
    s = {l for l in labels if l != NS}
    if not s:
        return NS
    if s == {SUPPORTS}:
        return SUPPORTS
    if s == {OPPOSES}:
        return OPPOSES
    return "conflicting"


def binomial_concordance(n_support, n_oppose):
    """Two-sided binomial test of concordant vs discordant significant calls."""
    from scipy.stats import binomtest
    n = n_support + n_oppose
    if n == 0:
        return {"n_support": 0, "n_oppose": 0, "n": 0, "frac": None, "p": None}
    r = binomtest(n_support, n, 0.5, alternative="two-sided")
    return {"n_support": int(n_support), "n_oppose": int(n_oppose), "n": int(n),
            "frac": float(n_support / n), "p": float(r.pvalue)}
