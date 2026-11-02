#!/usr/bin/env python3
"""
49 - Pathway enrichment of the TWAS risk genes, one cell type at a time.

The question is not "are these genes in pathways" but two sharper ones:

  1. Do the cell-type risk genes land in the pathways colorectal cancer is already known
     to run through - WNT, TGF-beta/BMP, RTK-MAPK, PI3K-AKT, p53, mismatch repair,
     NOTCH, cell cycle, apoptosis, chromatin remodelling?
  2. What else do they land in that is NOT one of those, recurrently enough and with
     enough genes to be worth calling a candidate mechanism?

Design
  universe     genes carried by at least one gene set in the Reactome / WikiPathways /
               PID collection (dhimmel/pathways)
  module       the risk genes of one TWAS cell type, mapped to Entrez, intersected with
               the universe; cell types with fewer than 5 mapped genes are not tested
  test         hypergeometric (over-representation), pathway sizes 5-500, overlap >= 2
  correction   Benjamini-Hochberg within each cell type (the unit a reader interprets),
               and separately across every test in the analysis
  null         1,000 random modules per cell type, drawn from the same universe at the
               same size, each scored through the identical pipeline. This is the control
               that matters: a module of any 90 genes hits some pathways, and the question
               is whether the real module hits more.

The established-vs-novel split is made by matching the pathway NAME against a fixed
family map, so it is reproducible and can be argued with; a pathway that matches none of
the families is called novel, which means "not one of the canonical CRC routes", not
"never before associated with cancer".

Outputs
  results/pathway_celltype_enrichment.csv        every test that reached overlap >= 2
  results/pathway_celltype_significant.csv       FDR < 0.1 within cell type
  results/pathway_celltype_null.csv              permutation summary per cell type
  results/pathway_families.csv                   established-family coverage per cell type
  results/pathway_novel_recurrent.csv            novel pathways, ranked by recurrence
  results/pathway_celltype.json                  headline numbers for the report
"""
import pandas as pd, numpy as np, json, re, os
from scipy.stats import hypergeom

rng = np.random.default_rng(20260919)
os.makedirs("results", exist_ok=True)

risk = pd.read_csv("results/risk_genes_harmonized.csv")
sym2ent = json.load(open("resources/gene_id_map.json"))["sym2ent"]
ent2sym = {v: k for k, v in sym2ent.items()}
risk["entrez"] = risk["gene"].map(sym2ent)
pw = pd.read_pickle("resources/pathways.pkl")

MIN_PW, MAX_PW, MIN_OVERLAP, MIN_MODULE = 5, 500, 2, 5
N_PERM = 1000


def bh(p):
    p = np.asarray(p, float)
    o = np.argsort(p)
    q = np.empty_like(p)
    q[o] = np.minimum(1, np.minimum.accumulate(
        (p[o] * len(p) / (np.arange(len(p)) + 1))[::-1])[::-1])
    return q


# ---------------------------------------------------------------- gene sets
sets, names, sources = [], [], []
for r in pw.itertuples(index=False):
    gs = set()
    for e in r.gene_list:
        try:
            gs.add(int(e))
        except Exception:
            pass
    if MIN_PW <= len(gs) <= MAX_PW:
        sets.append(gs); names.append(r.name); sources.append(r.source)

universe = sorted(set().union(*sets))
U = len(universe)
pos = {g: i for i, g in enumerate(universe)}
M = np.zeros((len(sets), U), dtype=bool)
for i, gs in enumerate(sets):
    M[i, [pos[g] for g in gs]] = True
K = M.sum(1)
DEG = M.sum(0)                     # how many gene sets each universe gene belongs to
edges = np.unique(np.quantile(DEG, np.linspace(0, 1, 11)))
BIN = np.clip(np.searchsorted(edges, DEG, side="right") - 1, 0, len(edges) - 2)
BIN_MEMBERS = [np.where(BIN == b)[0] for b in range(len(edges) - 1)]
print(f"{len(sets)} gene sets ({MIN_PW}-{MAX_PW} genes), universe {U} genes; "
      f"pathway-membership bins: {[len(m) for m in BIN_MEMBERS]}")


def matched_draw(idx):
    """A random module matched gene-for-gene on pathway membership count.

    Well-studied genes sit in more gene sets and are easier to discover, so a uniform
    draw makes any real gene list look enriched. This null removes that.
    """
    out = np.empty(len(idx), int)
    for i, gidx in enumerate(idx):
        pool = BIN_MEMBERS[BIN[gidx]]
        out[i] = pool[rng.integers(len(pool))]
    return np.unique(out)

# ---------------------------------------------------------------- the canonical routes
FAMILIES = {
    "WNT": r"wnt|beta[- ]catenin|destruction complex|\btcf\b|lef1",
    "TGF-β / BMP": r"tgf[- ]?beta|\bbmp\b|smad|activin|nodal",
    "RTK - MAPK": r"\bmapk\b|\berk\b|egfr|erbb|\bras signal|\braf\b|receptor tyrosine|"
                  r"fgfr|\bmet\b signaling|kit signal",
    "PI3K - AKT": r"pi3k|\bakt\b|mtor|pten",
    "p53": r"\bp53\b|tp53",
    "Mismatch repair": r"mismatch repair|\bmsh\b|dna damage|base excision|"
                       r"nucleotide excision|homologous recombination",
    "NOTCH": r"notch",
    "Cell cycle": r"cell cycle|mitotic|\bg1/s\b|\bg2/m\b|dna replication",
    "Apoptosis": r"apoptosis|apoptotic|caspase|death receptor",
    "Chromatin remodelling": r"chromatin|swi/snf|histone|hdac|methyltransferase complex",
}


def family_of(name):
    for fam, rx in FAMILIES.items():
        if re.search(rx, str(name), flags=re.I):
            return fam
    return "novel"


FAM = np.array([family_of(n) for n in names])
print("gene sets per canonical family:")
for fam in FAMILIES:
    print(f"  {fam:22s} {(FAM == fam).sum():4d}")
print(f"  {'(none of these)':22s} {(FAM == 'novel').sum():4d}")


# ---------------------------------------------------------------- one enrichment run
def enrich(idx):
    """idx = positions in `universe` of the module's genes.

    Returns (reportable mask, BH q, raw p, overlaps). Every gene set the module touches
    at all counts as a test and enters the correction - that is what over-representation
    tools correct over - while only sets reaching MIN_OVERLAP are reported. The
    permutation below runs the identical rule, so the null absorbs the same convention.
    """
    k = M[:, idx].sum(1)
    touched = k >= 1                      # the tests a module actually generates
    rep = k >= MIN_OVERLAP                # the ones worth reporting
    p = np.ones(len(sets))
    p[touched] = hypergeom.sf(k[touched] - 1, U, K[touched], len(idx))
    q = np.ones(len(sets))
    q[touched] = bh(p[touched])
    return rep, q, p, k


groups = [((t, c), g) for (t, c), g in risk.groupby(["tissue", "cell_type"])]
rows, null_rows, touched_p = [], [], []
modules = {}
for (tissue, ct), grp in groups:
    ents = {int(e) for e in grp["entrez"].dropna()} & set(universe)
    if len(ents) < MIN_MODULE:
        null_rows.append({"tissue": tissue, "cell_type": ct, "n_genes": len(grp),
                          "n_mapped": len(ents), "tested": False})
        continue
    idx = np.array([pos[g] for g in ents])
    modules[(tissue, ct)] = ents
    tested, q, p, k = enrich(idx)
    touched_p.extend(p[k >= 1].tolist())
    for i in np.where(tested)[0]:
        ov = sorted(ents & sets[i])
        rows.append({"tissue": tissue, "cell_type": ct, "pathway": names[i],
                     "source": sources[i], "family": FAM[i],
                     "established": FAM[i] != "novel",
                     "n_pathway": int(K[i]), "n_module": len(ents), "overlap": int(k[i]),
                     "p": float(p[i]), "fdr_celltype": float(q[i]),
                     "genes": "|".join(ent2sym.get(g, str(g)) for g in ov)})
    obs_q = int(((q < 0.1) & tested).sum())
    obs_p = int((p[tested] < 0.01).sum())
    # ---- permutation null: same module size, same pipeline, same two counts
    cq = np.empty(N_PERM, int); cp = np.empty(N_PERM, int)
    mq = np.empty(N_PERM, int); mp = np.empty(N_PERM, int)
    for b in range(N_PERM):
        ridx = rng.choice(U, size=len(idx), replace=False)
        tt, qq, pp, _ = enrich(ridx)
        cq[b] = int(((qq < 0.1) & tt).sum()); cp[b] = int(((pp < 0.01) & tt).sum())
        midx = matched_draw(idx)
        tt, qq, pp, _ = enrich(midx)
        mq[b] = int(((qq < 0.1) & tt).sum()); mp[b] = int(((pp < 0.01) & tt).sum())
    null_rows.append({"tissue": tissue, "cell_type": ct, "n_genes": len(grp),
                      "n_mapped": len(ents), "tested": True,
                      "n_tests": int(tested.sum()),
                      "observed_fdr": obs_q, "observed_p01": obs_p,
                      "null_fdr_mean": float(cq.mean()), "null_fdr_sd": float(cq.std(ddof=1)),
                      "null_p01_mean": float(cp.mean()), "null_p01_sd": float(cp.std(ddof=1)),
                      "null_p01_p95": float(np.quantile(cp, 0.95)),
                      "p_perm_fdr": float((cq >= obs_q).mean()),
                      "p_perm_p01": float((cp >= obs_p).mean()),
                      "matched_fdr_mean": float(mq.mean()),
                      "matched_fdr_sd": float(mq.std(ddof=1)),
                      "matched_p01_mean": float(mp.mean()),
                      "p_perm_matched_fdr": float((mq >= obs_q).mean()),
                      "p_perm_matched_p01": float((mp >= obs_p).mean()),
                      "fold_p01": float(obs_p / max(cp.mean(), 0.5))})
    print(f"{tissue:8s} {ct:4s} n={len(ents):3d} tests={int(tested.sum()):4d}  "
          f"FDR<0.1: {obs_q:3d} (uniform {cq.mean():4.1f} P={(cq >= obs_q).mean():.3f}; "
          f"matched {mq.mean():4.1f} P={(mq >= obs_q).mean():.3f})")

E = pd.DataFrame(rows)
# global correction over every test any module generated, not only the reported ones
tp = np.sort(np.array(touched_p))
gq = bh(tp)                      # monotone in p, so a lookup by position is exact
E["fdr_global"] = gq[np.searchsorted(tp, E["p"].values)]
E = E.sort_values(["tissue", "cell_type", "p"])
E.to_csv("results/pathway_celltype_enrichment.csv", index=False)
N = pd.DataFrame(null_rows)
N.to_csv("results/pathway_celltype_null.csv", index=False)

S = E[E.fdr_celltype < 0.1].copy()
S.to_csv("results/pathway_celltype_significant.csv", index=False)
print(f"\n{len(E)} tests, {len(S)} significant within cell type, "
      f"{int((E.fdr_global < 0.1).sum())} across the whole analysis")

# ---------------------------------------------------------------- pooled, all risk genes
allent = {int(e) for e in risk["entrez"].dropna()} & set(universe)
t_all, q_all, p_all, k_all = enrich(np.array([pos[g] for g in allent]))
POOL = pd.DataFrame({"pathway": names, "source": sources, "family": FAM,
                     "n_pathway": K, "overlap": k_all, "p": p_all, "fdr": q_all})
POOL = POOL[POOL.overlap >= MIN_OVERLAP].sort_values("p")
POOL["genes"] = [
    "|".join(ent2sym.get(g, str(g)) for g in sorted(allent & sets[names.index(n)]))
    for n in POOL.pathway]
POOL.to_csv("results/pathway_pooled_enrichment.csv", index=False)

# ---------------------------------------------------------------- family coverage
fam_rows = []
for (tissue, ct), ents in modules.items():
    sub = S[(S.tissue == tissue) & (S.cell_type == ct)]
    for fam in list(FAMILIES) + ["novel"]:
        f = sub[sub.family == fam]
        fam_rows.append({"tissue": tissue, "cell_type": ct, "family": fam,
                         "n_pathways": len(f),
                         "best_fdr": float(f.fdr_celltype.min()) if len(f) else np.nan,
                         "genes": "|".join(sorted({g for s_ in f.genes for g in s_.split("|")}))
                         if len(f) else ""})
F = pd.DataFrame(fam_rows)
F.to_csv("results/pathway_families.csv", index=False)

# ---------------------------------------------------------------- recurrent novel sets
NOV = (S[~S.established]
       .groupby(["pathway", "source"])
       .agg(n_celltypes=("cell_type", "nunique"),
            n_modules=("cell_type", "size"),
            best_fdr=("fdr_celltype", "min"),
            max_overlap=("overlap", "max"),
            n_pathway=("n_pathway", "first"),
            cell_types=("cell_type", lambda s: "|".join(sorted(set(s)))),
            tissues=("tissue", lambda s: "|".join(sorted(set(s)))),
            genes=("genes", lambda s: "|".join(sorted({g for x in s for g in x.split("|")}))))
       .reset_index()
       .sort_values(["n_celltypes", "best_fdr"], ascending=[False, True]))
NOV.to_csv("results/pathway_novel_recurrent.csv", index=False)

EST = (S[S.established]
       .groupby(["family"])
       .agg(n_pathways=("pathway", "nunique"), n_celltypes=("cell_type", "nunique"),
            best_fdr=("fdr_celltype", "min"),
            genes=("genes", lambda s: "|".join(sorted({g for x in s for g in x.split("|")}))))
       .reset_index().sort_values("best_fdr"))

# ---------------------------------------------------------------- what supports them
val = pd.read_csv("results/validated_risk_genes_all.csv")
validated = set(val.loc[val["validated"], "gene"]) if "validated" in val else set()
C = pd.read_csv("results/crispr_by_cellline.csv")
Sw = pd.read_csv("results/sirna_by_cellline.csv")
called = set(C[C.nominal].gene) | set(Sw[Sw.nominal & Sw.qc_pass].gene)


def support(genes):
    g = set(str(genes).split("|"))
    return len(g & validated), len(g & called), len(g)


for df in (NOV,):
    df[["n_validated", "n_screen_called", "n_genes_total"]] = [
        support(x) for x in df.genes]

gene_rows = []
for gsym in sorted({g for x in S.genes for g in x.split("|")}):
    sub = S[S.genes.str.split("|").apply(lambda L: gsym in L)]
    gene_rows.append({"gene": gsym, "n_pathways": len(sub),
                      "n_celltypes": sub.cell_type.nunique(),
                      "n_established": int(sub.established.sum()),
                      "n_novel": int((~sub.established).sum()),
                      "families": "|".join(sorted(set(sub.family))),
                      "cell_types": "|".join(sorted(set(sub.tissue + " " + sub.cell_type))),
                      "validated": gsym in validated, "screen_called": gsym in called})
G = pd.DataFrame(gene_rows).sort_values(["n_pathways", "gene"], ascending=[False, True])
G.to_csv("results/pathway_leading_genes.csv", index=False)

core = (S.groupby(["tissue", "cell_type"])
        .agg(n_sig=("pathway", "size"),
             n_core_genes=("genes", lambda s_: len({g for x in s_ for g in x.split("|")})),
             core=("genes", lambda s_: "|".join(sorted({g for x in s_ for g in x.split("|")}))))
        .reset_index())
core.to_csv("results/pathway_gene_core.csv", index=False)
print("\ngene core behind the significant pathways, per cell type:")
print(core.to_string(index=False))

J = {
    "n_sets": len(sets), "universe": U,
    "n_celltypes_tested": int(N.tested.sum()), "n_celltypes": len(groups),
    "n_tests": len(E), "n_sig_celltype": len(S),
    "n_sig_global": int((E.fdr_global < 0.1).sum()),
    "n_sig_established": int(S.established.sum()),
    "n_sig_novel": int((~S.established).sum()),
    "n_families_hit": int(S[S.established].family.nunique()),
    "families": EST.to_dict("records"),
    "null": N.to_dict("records"),
    "n_celltypes_above_null": int((N[N.tested].p_perm_p01 < 0.05).sum()),
    "n_celltypes_above_matched": int((N[N.tested].p_perm_matched_fdr < 0.05).sum()),
    "n_celltypes_above_null_fdr": int((N[N.tested].p_perm_fdr < 0.05).sum()),
    "pooled_top": POOL.head(15).to_dict("records"),
    "pooled_n_sig": int((POOL.fdr < 0.1).sum()),
    "novel_top": NOV.head(15).to_dict("records"),
    "leading_genes": G.head(15).to_dict("records"),
    "gene_core": core.to_dict("records"),
    "n_distinct_pathways": int(S.pathway.nunique()),
    "n_distinct_genes": int(len({g for x in S.genes for g in x.split("|")})),
    "per_celltype": S.groupby(["tissue", "cell_type"]).agg(
        n=("pathway", "size"), n_est=("established", "sum")).reset_index().to_dict("records"),
}
json.dump(J, open("results/pathway_celltype.json", "w"), indent=1, default=str)
print("\nestablished families hit:")
print(EST[["family", "n_pathways", "n_celltypes", "best_fdr"]].to_string(index=False))
print("\ntop recurrent novel pathways:")
print(NOV.head(12)[["pathway", "source", "n_celltypes", "max_overlap", "best_fdr",
                    "n_validated", "n_screen_called"]].to_string(index=False))
