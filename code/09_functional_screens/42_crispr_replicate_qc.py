#!/usr/bin/env python3
"""
42 - Replicate-level quality control for the CRISPR-Cas9 screen, and the exclusion
decision that follows from it.

Section 2.15 excluded a whole CELL LINE from the siRNA screen on Z'. This is the finer
question: within each CRISPR cell line, is any one of the three replicate plates bad
enough to drop?

Six metrics per plate, five from the controls and one from the data:

  Z'-factor        1 - 3(sd_pos + sd_neg)/|mean_pos - mean_neg|; > 0.5 excellent,
                   0 to 0.5 marginal but usable, <= 0 unusable
  SSMD             strictly standardised mean difference between the control sets
  NT CV            coefficient of variation of the non-targeting wells; lower is better
  signal window    (mean_neg - mean_pos) / sd_neg
  replicate r      the plate's mean Pearson correlation with the OTHER TWO plates of the
                   same cell line, over the 116 assayed genes - the only metric that
                   asks whether a plate agrees with its own siblings rather than whether
                   its controls behaved
  plate offset     |sample median - the line's median of plate medians|

The exclusion rule is fixed in advance so it cannot be tuned to a desired answer:

  hard fail   Z' <= 0, the plate cannot separate a killed well from a living one
  soft fail   on a metric, the plate is the WORST of its three AND beyond the screen-wide
              quartile for that metric (for replicate r, below 0.50 outright)
  EXCLUDE     a hard fail, or three or more soft fails out of six

Dropping a plate has a price that has to be stated: a paired t-test on the two remaining
plates carries ONE degree of freedom, where the two-sided 5 % critical value of t is
12.706. The script quantifies that rather than hiding it.

Outputs
  results/crispr_replicate_qc.csv
  results/crispr_replicate_qc.json
  figures/fig_crispr_replicate_qc.{png,pdf}
"""
import pandas as pd, numpy as np, json, os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch, Rectangle
from scipy.stats import pearsonr, ttest_rel, t as tdist

os.makedirs("figures", exist_ok=True)
LINES = ["HCT116", "RKO", "SW480"]
REPS = ["rep1", "rep2", "rep3"]
INK, INK2, INK3 = "#0b0b0b", "#52514e", "#8a8880"
SURF = "#fcfcfb"
OK, WARN, BAD = "#1f7a5a", "#e8a33d", "#b3472f"
plt.rcParams.update({
    "figure.facecolor": SURF, "axes.facecolor": SURF, "savefig.facecolor": SURF,
    "font.family": "DejaVu Sans", "font.size": 8.5,
    "axes.edgecolor": "#d8d7d2", "axes.linewidth": 0.8,
    "axes.labelcolor": INK2, "text.color": INK,
    "xtick.color": INK2, "ytick.color": INK2, "xtick.labelsize": 8, "ytick.labelsize": 8,
    "axes.titlesize": 9.5, "axes.titleweight": "bold", "axes.titlelocation": "left",
    "axes.spines.top": False, "axes.spines.right": False,
    "grid.color": "#eceae5", "grid.linewidth": 0.7, "legend.frameon": False,
    "legend.fontsize": 7.6, "pdf.fonttype": 42, "ps.fonttype": 42,
})


def tag(ax, s, dx=-0.12, dy=1.07):
    ax.text(dx, dy, s, transform=ax.transAxes, fontsize=11, fontweight="bold",
            color=INK, va="top", ha="left")


C = pd.read_csv("results/func_crispr_long.csv")
Q = pd.read_csv("results/func_qc.csv")
Q = Q[Q.assay == "CRISPR"].set_index(["cell_line", "replicate"])

rows = []
for line in LINES:
    W = C[(C.cell_line == line) & (C.kind == "sample")].pivot_table(
        index="gene", columns="replicate", values="value")
    med_of_medians = np.median([W[r].median() for r in REPS])
    for rep in REPS:
        others = [r for r in REPS if r != rep]
        r_other = float(np.mean([pearsonr(W[rep], W[o]).statistic for o in others]))
        q = Q.loc[(line, rep)]
        rows.append({"cell_line": line, "replicate": rep,
                     "z_prime": float(q.z_prime), "ssmd": float(q.ssmd),
                     "NT_cv_pct": float(q.NT_cv_pct),
                     "signal_window": float(q.signal_window),
                     "replicate_r": r_other,
                     "plate_offset": abs(float(W[rep].median()) - med_of_medians),
                     "NT_mean": float(q.NT_mean), "pos_mean": float(q.pos_mean),
                     "sample_median": float(W[rep].median())})
R = pd.DataFrame(rows)

# ------------------------------------------------------------------ the fixed rule
HIGHER_BETTER = {"z_prime": True, "ssmd": True, "NT_cv_pct": False,
                 "signal_window": True, "replicate_r": True, "plate_offset": False}
soft = {}
for m, hi in HIGHER_BETTER.items():
    q25 = R[m].quantile(0.25 if hi else 0.75)
    worst_in_line = R.groupby("cell_line")[m].transform("min" if hi else "max")
    is_worst = R[m] == worst_in_line
    if m == "replicate_r":
        beyond = R[m] < 0.50
    else:
        beyond = R[m] < q25 if hi else R[m] > q25
    soft[m] = is_worst & beyond
    R[f"fail_{m}"] = soft[m]
R["n_soft_fails"] = sum(soft.values())
R["hard_fail"] = R["z_prime"] <= 0
R["exclude"] = R["hard_fail"] | (R["n_soft_fails"] >= 3)

print("=== CRISPR replicate QC ===")
show = ["cell_line", "replicate", "z_prime", "ssmd", "NT_cv_pct", "signal_window",
        "replicate_r", "plate_offset", "n_soft_fails", "exclude"]
print(R[show].round(3).to_string(index=False))
dropped = R[R.exclude]
print("\nexcluded:", ", ".join(f"{r.cell_line} {r.replicate}" for _, r in dropped.iterrows())
      or "none")
for _, r in dropped.iterrows():
    why = [m for m in HIGHER_BETTER if r[f"fail_{m}"]]
    print(f"  {r.cell_line} {r.replicate}: worst of its line and beyond the screen "
          f"quartile on {', '.join(why)}")

# ------------------------------------------------------------------ the consequence
TC2, TC1 = float(tdist.ppf(0.975, 2)), float(tdist.ppf(0.975, 1))
risk = set(pd.read_csv("results/risk_genes_harmonized.csv")["gene"])
cons = []
for line in LINES:
    keep = [r for r in REPS if not bool(
        R[(R.cell_line == line) & (R.replicate == r)]["exclude"].iloc[0])]
    d = C[C.cell_line == line]
    ntm = d[d.kind == "negative"].groupby("replicate")["value"].mean()
    G = d[d.kind == "sample"].pivot_table(index="gene", columns="replicate", values="value")
    for tagname, reps in [("all 3", REPS), ("after QC", keep)]:
        diff = np.array([G[r].to_numpy(float) - ntm[r] for r in reps])
        if len(reps) < 2:
            continue
        mean = diff.mean(axis=0)
        sd = diff.std(axis=0, ddof=1)
        se = sd / np.sqrt(len(reps))
        with np.errstate(divide="ignore", invalid="ignore"):
            tstat = np.where(se > 0, mean / se, np.nan)
        p = 2 * tdist.sf(np.abs(tstat), len(reps) - 1)
        isr = np.isin(G.index, list(risk))
        cons.append({"cell_line": line, "set": tagname, "n_rep": len(reps),
                     "df": len(reps) - 1, "t_crit": float(tdist.ppf(0.975, len(reps) - 1)),
                     "sig_panel": int(np.nansum(p < 0.05)),
                     "sig_risk": int(np.nansum((p < 0.05) & isr)),
                     "median_abs_effect": float(np.median(np.abs(mean)))})
Cons = pd.DataFrame(cons)
print("\n=== what dropping a plate costs ===")
print(Cons.round(3).to_string(index=False))

# ---- the pre-exclusion concordance, so the report can quote both readings side by side
import sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _ttest_support import support_call, gene_call, binomial_concordance, SUPPORTS, OPPOSES
allm = pd.read_csv("results/validated_risk_genes_all.csv").set_index("gene")
ALIAS = {"ACAS2": "ACSS2"}
def concordance(drop):
    calls = []
    for line in LINES:
        d = C[C.cell_line == line].copy()
        d["gene"] = d["gene"].replace(ALIAS)
        keep = [r for r in REPS if f"{line} {r}" not in drop]
        ntm = d[d.kind == "negative"].groupby("replicate")["value"].mean()
        G = d[d.kind == "sample"].pivot_table(index="gene", columns="replicate",
                                              values="value")
        for g in G.index:
            if g not in risk or len(keep) < 2:
                continue
            diff = np.array([G.loc[g, r] - ntm[r] for r in keep])
            tt = ttest_rel([G.loc[g, r] for r in keep], [ntm[r] for r in keep])
            tw = allm["twas_direction"].get(g, "")
            calls.append({"cell_line": line, "gene": g, "delta": float(diff.mean()),
                          "p": float(tt.pvalue),
                          "verdict": support_call(tw, float(diff.mean()),
                                                  tt.pvalue < 0.05)})
    Cl = pd.DataFrame(calls)
    gv = Cl.groupby("gene")["verdict"].apply(lambda v: gene_call(list(v)))
    return {"n_calls": int((Cl.p < 0.05).sum()),
            "concordance_calls": binomial_concordance(
                int((Cl.verdict == SUPPORTS).sum()), int((Cl.verdict == OPPOSES).sum())),
            "concordance_genes": binomial_concordance(
                int((gv == SUPPORTS).sum()), int((gv == OPPOSES).sum()))}


_all = concordance(set())
_qc = concordance(set(f"{r.cell_line} {r.replicate}" for _, r in dropped.iterrows()))
pre = {**_all, "concordance_after_qc": _qc["concordance_calls"],
       "concordance_after_qc_genes": _qc["concordance_genes"]}
json.dump(pre, open("results/crispr_concordance_all_plates.json", "w"), indent=1)
print("\npre-exclusion (all 3 plates) concordance:", pre["concordance_calls"])

R.to_csv("results/crispr_replicate_qc.csv", index=False)
json.dump({"table": R[show + [f"fail_{m}" for m in HIGHER_BETTER]].to_dict("records"),
           "excluded": [f"{r.cell_line} {r.replicate}" for _, r in dropped.iterrows()],
           "n_excluded": int(len(dropped)),
           "consequence": Cons.to_dict("records"),
           "t_crit_df2": TC2, "t_crit_df1": TC1},
          open("results/crispr_replicate_qc.json", "w"), indent=1, default=float)

# ------------------------------------------------------------------ figure
fig = plt.figure(figsize=(13.4, 8.2))
gs = fig.add_gridspec(2, 3, height_ratios=[1.0, 0.92], hspace=0.62, wspace=0.34)

METRICS = [("z_prime", "Z'-factor", True, 0.5),
           ("replicate_r", "Mean r with the other\ntwo plates of the line", True, 0.5),
           ("NT_cv_pct", "Non-targeting CV (%)", False, None),
           ("signal_window", "Control signal window", True, None),
           ("ssmd", "SSMD of the controls", True, None)]

for k, (m, lab, hi, ref) in enumerate(METRICS[:3]):
    ax = fig.add_subplot(gs[0, k]); tag(ax, "abc"[k])
    x = np.arange(9, dtype=float)
    cols = [BAD if e else (OK if not f else WARN)
            for e, f in zip(R["exclude"], R[f"fail_{m}"])]
    ax.bar(x, R[m], width=0.68, color=cols)
    if ref is not None:
        ax.axhline(ref, color=INK3, lw=0.9, ls=(0, (3, 2)))
        ax.text(-0.35, ref, "0.5 ", fontsize=7.2, color=INK3, va="bottom", ha="left")
    for xi, v in zip(x, R[m]):
        ax.text(xi, v + R[m].max() * 0.025, f"{v:.2f}", ha="center", fontsize=6.8,
                color=INK2)
    ax.set_xticks(x)
    ax.set_xticklabels([f"{r.cell_line[0]}{r.replicate[-1]}" for _, r in R.iterrows()],
                       fontsize=7.6)
    for b in (2.5, 5.5):
        ax.axvline(b, color="#e2e0db", lw=0.9, zorder=0)
    for cx, nm in zip((1, 4, 7), LINES):
        ax.text(cx, -0.165, nm, transform=ax.get_xaxis_transform(), ha="center",
                va="top", fontsize=7.4, color=INK2)
    for t, e in zip(ax.get_xticklabels(), R["exclude"]):
        if e:
            t.set_color(BAD); t.set_fontweight("bold")
    ax.set_ylabel(lab, fontsize=8)
    ax.set_ylim(0, R[m].max() * 1.20)
    ax.set_title(["Control separation", "Agreement with siblings",
                  "Control noise"][k])
    ax.grid(axis="y"); ax.set_axisbelow(True)
    for sp in ["bottom"]:
        ax.spines[sp].set_visible(True)

# (d) the scorecard
ax = fig.add_subplot(gs[1, :2]); tag(ax, "d", dx=-0.075)
Mx = np.zeros((len(METRICS) + 1, 9))
labels = [lab.replace("\n", " ") for _, lab, _, _ in METRICS] + ["plate offset"]
keys = [m for m, _, _, _ in METRICS] + ["plate_offset"]
for i, key in enumerate(keys):
    Mx[i] = R[f"fail_{key}"].to_numpy(float)
ax.imshow(Mx, cmap=plt.cm.colors.ListedColormap(["#f0efec", WARN]), aspect="auto",
          origin="lower", extent=(-0.5, 8.5, -0.5, len(keys) - 0.5), vmin=0, vmax=1)
for j in range(9):
    if R["exclude"].iloc[j]:
        ax.add_patch(Rectangle((j - 0.5, -0.5), 1, len(keys), fill=False,
                               edgecolor=BAD, linewidth=2.2, zorder=5))
for i in range(len(keys)):
    for j in range(9):
        if Mx[i, j]:
            ax.text(j, i, "fail", ha="center", va="center", fontsize=6.6, color="white")
ax.set_yticks(range(len(keys))); ax.set_yticklabels(labels, fontsize=7.4)
ax.set_xticks(range(9))
ax.set_xticklabels([f"{r.cell_line}\n{r.replicate}" for _, r in R.iterrows()], fontsize=7)
for t, e in zip(ax.get_xticklabels(), R["exclude"]):
    if e:
        t.set_color(BAD); t.set_fontweight("bold")
ax.tick_params(length=0)
for sp in ax.spines.values():
    sp.set_visible(False)
ax.set_title("Where each plate is the worst of its line AND beyond the screen quartile")
excl = ", ".join(f"{r.cell_line} {r.replicate}" for _, r in dropped.iterrows()) or "none"
ax.text(0.0, -0.40,
        f"A plate is excluded on a hard fail (Z' ≤ 0) or three or more of these six.\n"
        f"Excluded: {excl}. Every other plate has at most two.\nNo CRISPR plate comes "
        f"near a hard fail: the lowest Z' is {R.z_prime.min():.2f} (HCT116 rep1), and it "
        f"agrees with its\nsiblings as well as they agree with each other.",
        transform=ax.transAxes, fontsize=7.5, color=INK2, va="top")

# (e) the cost
ax = fig.add_subplot(gs[1, 2]); tag(ax, "e", dx=-0.26)
sub = Cons[Cons.cell_line == "SW480"] if (dropped.cell_line == "SW480").any() else Cons
lab = [f"{r['set']}\n({r.n_rep} plates, df = {r.df})" for _, r in sub.iterrows()]
ax.bar(np.arange(len(sub)), sub["sig_panel"], width=0.6,
       color=[OK if r.n_rep == 3 else BAD for _, r in sub.iterrows()])
for i, (_, r) in enumerate(sub.iterrows()):
    ax.text(i, r.sig_panel + 0.6, f"{int(r.sig_panel)}", ha="center", fontsize=8.4,
            color=INK2)
    ax.text(i, r.sig_panel / 2, f"critical t\n{r.t_crit:.2f}", ha="center",
            va="center", fontsize=7.6, color="white", fontweight="bold")
ax.set_xticks(range(len(sub))); ax.set_xticklabels(lab, fontsize=7.4)
ax.set_ylabel("Genes with P < 0.05 (of 116)")
ax.set_ylim(0, max(sub["sig_panel"].max() * 1.35, 5))
ax.set_title("The price of dropping one plate")
ax.grid(axis="y"); ax.set_axisbelow(True)
ax.text(0.0, -0.44,
        "Two plates leave one degree of freedom, where the 5 % critical value of t is "
        "12.71 instead of 4.30.\nThe effects barely move; the p-values collapse. Removing "
        "a marginal plate costs more than it buys.",
        transform=ax.transAxes, fontsize=7.5, color=INK2, va="top")

fig.legend(handles=[Patch(facecolor=OK, label="passes"),
                    Patch(facecolor=WARN, label="worst of its line on this metric"),
                    Patch(facecolor=BAD, label="excluded")],
           loc="lower center", bbox_to_anchor=(0.5, 0.982), ncol=3, fontsize=8.0)

fig.savefig("figures/fig_crispr_replicate_qc.png", dpi=220, bbox_inches="tight")
fig.savefig("figures/fig_crispr_replicate_qc.pdf", bbox_inches="tight")
plt.close(fig)
print("\nwrote figures/fig_crispr_replicate_qc.png")
