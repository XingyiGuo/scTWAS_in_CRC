#!/usr/bin/env python3
"""
51 - The same enrichment as script 49, with normal and adenoma pooled per cell type.

Script 49 tested 16 modules, one per tissue x cell type, and most were too small to say
anything: seven had fewer than 20 genes in the pathway universe. Pooling the normal and
adenoma module of the same cell type is the obvious way to buy power, and the overlap
between them is small (2 to 11 genes), so it is close to a union of two independent gene
lists rather than a re-count of one.

Everything else is held identical to script 49 so the two analyses can be compared
directly: the same 1,626 Reactome / WikiPathways / PID gene sets, the same universe, the
same hypergeometric test, the same correction over every gene set a module touches, and
the same two nulls - size-matched, and matched gene-for-gene on how many gene sets a gene
belongs to. Only the module definition changes.

The question this is asked for is specific: does the extra power bring the canonical CRC
pathways out? The answer has to survive the annotation-matched null to count, because
WNT, cell cycle and p53 are among the most heavily annotated gene sets in any collection
and a larger module hits them more often by construction.

Outputs
  results/pathway_merged_enrichment.csv     every test at overlap >= 2
  results/pathway_merged_significant.csv    FDR < 0.1 within the merged module
  results/pathway_merged_null.csv           both nulls per merged module
  results/pathway_merged_vs_split.csv       merged against script 49, cell type by cell type
  results/pathway_merged.json               headline numbers for the report
"""
import pandas as pd, numpy as np, json, re, os
from scipy.stats import hypergeom

rng = np.random.default_rng(20260920)
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
DEG = M.sum(0)
edges = np.unique(np.quantile(DEG, np.linspace(0, 1, 11)))
BIN = np.clip(np.searchsorted(edges, DEG, side="right") - 1, 0, len(edges) - 2)
BIN_MEMBERS = [np.where(BIN == b)[0] for b in range(len(edges) - 1)]

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


def matched_draw(idx):
    out = np.empty(len(idx), int)
    for i, gidx in enumerate(idx):
        pool = BIN_MEMBERS[BIN[gidx]]
        out[i] = pool[rng.integers(len(pool))]
    return np.unique(out)


def enrich(idx):
    k = M[:, idx].sum(1)
    touched = k >= 1
    rep = k >= MIN_OVERLAP
    p = np.ones(len(sets))
    p[touched] = hypergeom.sf(k[touched] - 1, U, K[touched], len(idx))
    q = np.ones(len(sets))
    q[touched] = bh(p[touched])
    return rep, q, p, k


# ---------------------------------------------------------------- merged modules
CTLAB = {"ABS": "absorptive", "ASC": "absorptive stem-like", "CT": "crypt top",
         "EE": "enteroendocrine", "GOB": "goblet", "SSC": "secretory stem-like",
         "STM": "stem", "TAC": "transit-amplifying", "TUF": "tuft"}
rows, null_rows, touched_p = [], [], []
modules = {}
for ct, grp in risk.groupby("cell_type"):
    ents = {int(e) for e in grp["entrez"].dropna()} & set(universe)
    by_t = {t: {int(e) for e in g["entrez"].dropna()} & set(universe)
            for t, g in grp.groupby("tissue")}
    n_norm, n_aden = len(by_t.get("Normal", set())), len(by_t.get("Adenoma", set()))
    n_both = len(by_t.get("Normal", set()) & by_t.get("Adenoma", set()))
    base = {"cell_type": ct, "label": CTLAB.get(ct, ct), "n_genes": grp.gene.nunique(),
            "n_mapped": len(ents), "n_normal": n_norm, "n_adenoma": n_aden,
            "n_shared": n_both}
    if len(ents) < MIN_MODULE:
        null_rows.append({**base, "tested": False})
        continue
    idx = np.array([pos[g] for g in ents])
    modules[ct] = ents
    tested, q, p, k = enrich(idx)
    touched_p.extend(p[k >= 1].tolist())
    for i in np.where(tested)[0]:
        ov = sorted(ents & sets[i])
        tis = []
        if by_t.get("Normal", set()) & sets[i]:
            tis.append("normal")
        if by_t.get("Adenoma", set()) & sets[i]:
            tis.append("adenoma")
        rows.append({"cell_type": ct, "label": CTLAB.get(ct, ct), "pathway": names[i],
                     "source": sources[i], "family": FAM[i],
                     "established": FAM[i] != "novel", "n_pathway": int(K[i]),
                     "n_module": len(ents), "overlap": int(k[i]),
                     "p": float(p[i]), "fdr_merged": float(q[i]),
                     "tissues_contributing": "+".join(tis),
                     "genes": "|".join(ent2sym.get(g, str(g)) for g in ov)})
    obs_q = int(((q < 0.1) & tested).sum())
    cq = np.empty(N_PERM, int); mq = np.empty(N_PERM, int)
    for b in range(N_PERM):
        tt, qq, pp, _ = enrich(rng.choice(U, size=len(idx), replace=False))
        cq[b] = int(((qq < 0.1) & tt).sum())
        tt, qq, pp, _ = enrich(matched_draw(idx))
        mq[b] = int(((qq < 0.1) & tt).sum())
    null_rows.append({**base, "tested": True, "n_tests": int(tested.sum()),
                      "observed_fdr": obs_q,
                      "null_fdr_mean": float(cq.mean()), "null_fdr_sd": float(cq.std(ddof=1)),
                      "matched_fdr_mean": float(mq.mean()),
                      "matched_fdr_sd": float(mq.std(ddof=1)),
                      "p_perm_fdr": float((cq >= obs_q).mean()),
                      "p_perm_matched_fdr": float((mq >= obs_q).mean())})
    print(f"{CTLAB.get(ct, ct):22s} n={len(ents):3d} (N {n_norm}, A {n_aden}, both {n_both}) "
          f"tests={int(tested.sum()):4d}  FDR<0.1: {obs_q:3d} "
          f"(size {cq.mean():5.1f} P={(cq >= obs_q).mean():.3f}; "
          f"matched {mq.mean():5.1f} P={(mq >= obs_q).mean():.3f})")

E = pd.DataFrame(rows)
tp = np.sort(np.array(touched_p))
gq = bh(tp)
E["fdr_global"] = gq[np.searchsorted(tp, E["p"].values)]
E = E.sort_values(["cell_type", "p"])
E.to_csv("results/pathway_merged_enrichment.csv", index=False)
N = pd.DataFrame(null_rows)
N.to_csv("results/pathway_merged_null.csv", index=False)
S = E[E.fdr_merged < 0.1].copy()
S.to_csv("results/pathway_merged_significant.csv", index=False)

# ---------------------------------------------------------------- merged vs split
SP = pd.read_csv("results/pathway_celltype_significant.csv")
SN = pd.read_csv("results/pathway_celltype_null.csv")
cmp_rows = []
for ct in sorted(set(risk.cell_type)):
    sp = SP[SP.cell_type == ct]
    sn = SN[(SN.cell_type == ct) & (SN.tested)]
    me = S[S.cell_type == ct]
    mn = N[(N.cell_type == ct) & (N.tested)]
    cmp_rows.append({
        "cell_type": ct, "label": CTLAB.get(ct, ct),
        "split_modules": int(len(sn)),
        "split_mapped": int(sn.n_mapped.sum()),
        "split_sig": int(len(sp)),
        "split_established": int(sp.established.sum()) if len(sp) else 0,
        "split_best_matched_p": float(sn.p_perm_matched_fdr.min()) if len(sn) else np.nan,
        "merged_mapped": int(mn.n_mapped.iat[0]) if len(mn) else 0,
        "merged_sig": int(len(me)),
        "merged_established": int(me.established.sum()) if len(me) else 0,
        "merged_matched_p": float(mn.p_perm_matched_fdr.iat[0]) if len(mn) else np.nan,
        "merged_families": "|".join(sorted(set(me[me.established].family))) if len(me) else "",
        "split_families": "|".join(sorted(set(sp[sp.established].family))) if len(sp) else "",
    })
C = pd.DataFrame(cmp_rows)
C.to_csv("results/pathway_merged_vs_split.csv", index=False)
print("\nmerged against split, cell type by cell type:")
print(C[["label", "split_mapped", "merged_mapped", "split_sig", "merged_sig",
         "split_established", "merged_established", "merged_matched_p"]].to_string(index=False))

fam_rows = []
for fam in list(FAMILIES) + ["novel"]:
    m = S[S.family == fam]
    sp = SP[SP.family == fam]
    fam_rows.append({"family": fam,
                     "merged_sets": int(m.pathway.nunique()),
                     "merged_celltypes": int(m.cell_type.nunique()),
                     "merged_best_fdr": float(m.fdr_merged.min()) if len(m) else np.nan,
                     "split_sets": int(sp.pathway.nunique()),
                     "merged_genes": "|".join(sorted({g for x in m.genes for g in x.split("|")})),
                     })
F = pd.DataFrame(fam_rows)
F.to_csv("results/pathway_merged_families.csv", index=False)
print("\ncanonical families, merged vs split:")
print(F[F.family != "novel"][["family", "split_sets", "merged_sets", "merged_celltypes",
                              "merged_best_fdr"]].to_string(index=False))

surv = N[(N.tested) & (N.p_perm_matched_fdr < 0.05)]
J = {
    "n_modules": int(N.tested.sum()), "n_celltypes": int(risk.cell_type.nunique()),
    "n_tests": len(E), "n_sig": len(S),
    "n_sig_established": int(S.established.sum()),
    "n_sig_novel": int((~S.established).sum()),
    "n_distinct_pathways": int(S.pathway.nunique()),
    "n_core_genes": int(len({g for x in S.genes for g in x.split("|")})),
    "n_families_hit": int(S[S.established].family.nunique()),
    "families_hit": sorted(set(S[S.established].family)),
    "families_missed": [f for f in FAMILIES if f not in set(S[S.established].family)],
    "surviving": surv[["label", "n_mapped", "observed_fdr", "matched_fdr_mean",
                       "p_perm_matched_fdr", "p_perm_fdr"]].to_dict("records"),
    "modules": N[N.tested][["label", "n_mapped", "n_normal", "n_adenoma", "n_shared",
                            "n_tests", "observed_fdr", "null_fdr_mean", "matched_fdr_mean",
                            "p_perm_fdr", "p_perm_matched_fdr"]].to_dict("records"),
    "compare": C.to_dict("records"),
    "families": F.to_dict("records"),
    "top": S.nsmallest(20, "p")[["label", "pathway", "family", "overlap", "n_pathway",
                                 "p", "fdr_merged", "tissues_contributing",
                                 "genes"]].to_dict("records"),
    "established_top": S[S.established].nsmallest(12, "p")[
        ["label", "pathway", "family", "overlap", "p", "fdr_merged",
         "tissues_contributing", "genes"]].to_dict("records"),
    "split_n_sig": int(len(SP)), "split_n_established": int(SP.established.sum()),
    "n_new_pathways": int(len(set(S.pathway) - set(SP.pathway))),
    "new_established": sorted(set(S[S.established].pathway) - set(SP[SP.established].pathway)),
    "both_tissues": int((S.tissues_contributing == "normal+adenoma").sum()),
}
json.dump(J, open("results/pathway_merged.json", "w"), indent=1, default=str)
print(f"\n{len(E)} tests, {len(S)} significant, {int(S.established.sum())} canonical; "
      f"{J['n_new_pathways']} pathways not seen in the split analysis; "
      f"{J['both_tissues']} of the significant sets draw genes from both tissues")
print("\ncanonical hits in the merged analysis:")
print(S[S.established][["label", "pathway", "family", "overlap", "fdr_merged",
                        "tissues_contributing", "genes"]].to_string(index=False))
