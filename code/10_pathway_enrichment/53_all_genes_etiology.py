#!/usr/bin/env python3
"""
53 - The whole risk-gene set, and four ways of asking what it means for CRC aetiology.

Sections 2.18 and 2.19 asked the question cell type by cell type. This asks it of all
347 risk genes at once, and then goes past over-representation, because over-representation
on its own has three blind spots that matter here:

  it ignores direction        a gene set counts the same whether its risk genes raise or
                              lower risk, so it cannot say which way a pathway pushes
  it ignores non-membership   a risk gene one step away from the WNT pathway contributes
                              nothing, even though that is exactly how a germline modifier
                              of a somatic pathway would look
  it is blind to exposure     the canonical CRC families are somatic driver pathways; the
                              epidemiology of colorectal cancer is about diet, alcohol,
                              bile acids, inflammation and carcinogen handling, and no
                              pathway database is organised around those

So this script runs five analyses on the same 347 genes:

  A  over-representation of the whole set against both nulls of section 2.18
  B  the same, split into risk-increasing and protective genes, which puts direction back
  C  a membership audit: which risk genes are literally in each canonical CRC family,
     enriched or not, with their direction and their spatial and functional evidence
  D  gene-level network proximity to each canonical family, degree-matched - the question
     "is this gene one step from WNT" rather than "is this gene annotated to WNT"
  E  aetiology themes: gene sets grouped by the exposure biology colorectal cancer is
     known to run on, tested the same way, with the direction balance of each theme

Outputs
  results/all347_enrichment.csv          A: every test, whole set
  results/all347_by_direction.csv        B: the same, per direction
  results/all347_nulls.json              A + B: both nulls for all three gene lists
  results/canonical_membership.csv       C: risk genes inside each canonical family
  results/gene_pathway_proximity.csv     D: gene x family network distance and z
  results/etiology_themes.csv            E: theme-level enrichment and direction
  results/etiology_gene_assignment.csv   E: every risk gene with its themes
  results/all347_etiology.json           headline numbers for the report
"""
import pandas as pd, numpy as np, json, re, os, pickle
import networkx as nx
from scipy.stats import hypergeom, binomtest, fisher_exact

rng = np.random.default_rng(20260921)
np.random.seed(20260921)
os.makedirs("results", exist_ok=True)

risk = pd.read_csv("results/risk_genes_harmonized.csv")
sym2ent = json.load(open("resources/gene_id_map.json"))["sym2ent"]
ent2sym = {v: k for k, v in sym2ent.items()}
risk["entrez"] = risk["gene"].map(sym2ent)
pw = pd.read_pickle("resources/pathways.pkl")

MIN_PW, MAX_PW, MIN_OVERLAP = 5, 500, 2
N_PERM = 1000

# one row per gene: direction is the direction of its strongest association
gene = (risk.sort_values("p_best")
        .groupby("gene")
        .agg(entrez=("entrez", "first"), direction=("direction", "first"),
             tissues=("tissue", lambda s: "|".join(sorted(set(s)))),
             n_ct=("cell_type", "nunique"), p_best=("p_best", "min"))
        .reset_index())
print(f"{len(gene)} risk genes; {int((gene.direction == 'risk-increasing').sum())} "
      f"risk-increasing, {int((gene.direction == 'protective').sum())} protective")


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

# the aetiology of colorectal cancer, as themes a pathway database is not organised around
THEMES = {
    "Carcinogen & drug metabolism": r"glucuronidation|biological oxidation|cytochrome p450|"
                                    r"xenobiotic|phase i|phase ii conjugation|drug metab|"
                                    r"glutathione conjugation|sulfur amino acid",
    "Alcohol & aldehyde handling": r"ethanol|alcohol|aldehyde",
    "Bile acid, lipid & fatty acid": r"bile acid|bile salt|fatty acid|triglyceride|"
                                     r"lipoprotein|cholesterol|sterol",
    "One-carbon, folate & methylation": r"folate|one[- ]carbon|methionine|s-adenosyl|"
                                        r"methylation of|homocysteine",
    "Immune surveillance & inflammation": r"antigen (presentation|processing)|interferon|"
                                          r"\bmhc\b|immunoregulatory|interleukin|"
                                          r"inflammat|toll[- ]like|nf-?kappa|complement",
    "Oxidative stress & redox": r"reactive oxygen|oxidative stress|detoxification of|"
                                r"peroxiredoxin|thioredoxin|glutathione|superoxide",
    "Barrier, mucin & adhesion": r"mucin|o-linked glycosylation|tight junction|"
                                 r"cell[- ]cell junction|adherens|laminin|integrin|"
                                 r"extracellular matrix|focal adhesion|nectin|cadherin",
    "Energy metabolism & insulin": r"insulin|glycolysis|gluconeogenesis|citric acid|"
                                   r"respiratory electron|oxidative phosphorylation|"
                                   r"electron transport|tca cycle|adipo",
    "Transport & absorption": r"\bslc\b|solute carrier|transport of (vitamins|small|"
                              r"nucleosides|inorganic|bile)|absorption|ion channel transport",
    "Protein turnover & stress": r"proteasom|ubiquitin|unfolded protein|autophagy|"
                                 r"chaperone|lysosom|endosomal",
}


def label_of(name, mapping):
    for lab, rx in mapping.items():
        if re.search(rx, str(name), flags=re.I):
            return lab
    return None


FAM = np.array([label_of(n, FAMILIES) or "novel" for n in names])
THEME = np.array([label_of(n, THEMES) or "" for n in names])
print("\ngene sets per aetiology theme:")
for t in THEMES:
    print(f"  {t:36s} {(THEME == t).sum():4d}")
print(f"  {'(no theme)':36s} {(THEME == '').sum():4d}")


def enrich(idx):
    k = M[:, idx].sum(1)
    touched = k >= 1
    rep = k >= MIN_OVERLAP
    p = np.ones(len(sets))
    p[touched] = hypergeom.sf(k[touched] - 1, U, K[touched], len(idx))
    q = np.ones(len(sets))
    q[touched] = bh(p[touched])
    return rep, q, p, k


def matched_draw(idx):
    out = np.empty(len(idx), int)
    for i, g in enumerate(idx):
        pool = BIN_MEMBERS[BIN[g]]
        out[i] = pool[rng.integers(len(pool))]
    return np.unique(out)


def run(label, genes_entrez):
    ents = sorted(set(genes_entrez) & set(universe))
    idx = np.array([pos[g] for g in ents])
    tested, q, p, k = enrich(idx)
    rows = []
    for i in np.where(tested)[0]:
        ov = sorted(set(ents) & sets[i])
        rows.append({"gene_list": label, "pathway": names[i], "source": sources[i],
                     "family": FAM[i], "theme": THEME[i], "established": FAM[i] != "novel",
                     "n_pathway": int(K[i]), "n_list": len(ents), "overlap": int(k[i]),
                     "p": float(p[i]), "fdr": float(q[i]),
                     "genes": "|".join(ent2sym.get(g, str(g)) for g in ov)})
    obs = int(((q < 0.1) & tested).sum())
    cq = np.empty(N_PERM, int); mq = np.empty(N_PERM, int)
    for b in range(N_PERM):
        tt, qq, pp, _ = enrich(rng.choice(U, size=len(idx), replace=False))
        cq[b] = int(((qq < 0.1) & tt).sum())
        tt, qq, pp, _ = enrich(matched_draw(idx))
        mq[b] = int(((qq < 0.1) & tt).sum())
    null = {"gene_list": label, "n_genes": len(genes_entrez), "n_mapped": len(ents),
            "n_tests": int(tested.sum()), "observed": obs,
            "null_mean": float(cq.mean()), "null_sd": float(cq.std(ddof=1)),
            "matched_mean": float(mq.mean()), "matched_sd": float(mq.std(ddof=1)),
            "p_perm": float((cq >= obs).mean()),
            "p_perm_matched": float((mq >= obs).mean())}
    print(f"{label:18s} n={len(ents):3d}  FDR<0.1: {obs:3d}  "
          f"(size {cq.mean():5.1f} P={(cq >= obs).mean():.3f}; "
          f"matched {mq.mean():5.1f} P={(mq >= obs).mean():.3f})")
    return pd.DataFrame(rows), null


# ---------------------------------------------------------------- A and B
ALL = {int(e) for e in gene.entrez.dropna()}
UP = {int(e) for e in gene.loc[gene.direction == "risk-increasing", "entrez"].dropna()}
DOWN = {int(e) for e in gene.loc[gene.direction == "protective", "entrez"].dropna()}
EA, nA = run("all 347", ALL)
EU, nU = run("risk-increasing", UP)
ED, nD = run("protective", DOWN)
EA.sort_values("p").to_csv("results/all347_enrichment.csv", index=False)
pd.concat([EU, ED]).sort_values(["gene_list", "p"]).to_csv(
    "results/all347_by_direction.csv", index=False)
json.dump([nA, nU, nD], open("results/all347_nulls.json", "w"), indent=1)

# ---------------------------------------------------------------- C membership audit
val = pd.read_csv("results/validated_risk_genes_all.csv")
validated = set(val.loc[val["validated"], "gene"]) if "validated" in val else set()
cls = dict(zip(val.gene, val.get("class", pd.Series(index=val.index, dtype=object))))
C = pd.read_csv("results/crispr_by_cellline.csv")
Sw = pd.read_csv("results/sirna_by_cellline.csv")
called = set(C[C.nominal].gene) | set(Sw[Sw.nominal & Sw.qc_pass].gene)
drg = pd.read_csv("results/gene_druggability.csv").set_index("gene")

fam_union = {f: set().union(*[sets[i] for i in np.where(FAM == f)[0]]) for f in FAMILIES}
rows = []
for f, gs in fam_union.items():
    hits = sorted(ALL & gs)
    k, n = len(hits), len(ALL & set(universe))
    p = hypergeom.sf(k - 1, U, len(gs), n)
    up = sum(1 for e in hits
             if gene.loc[gene.entrez == e, "direction"].iat[0] == "risk-increasing")
    for e in hits:
        g = ent2sym.get(e, str(e))
        r = gene[gene.entrez == e].iloc[0]
        rows.append({"family": f, "gene": g, "direction": r.direction,
                     "tissues": r.tissues, "n_celltypes": int(r.n_ct),
                     "p_twas": float(r.p_best),
                     "validated": g in validated, "class": cls.get(g, ""),
                     "screen_called": g in called,
                     "tier": drg.tier.get(g, ""), "family_size": len(gs),
                     "family_p": float(p), "family_hits": k, "family_up": up})
MEM = pd.DataFrame(rows)
MEM.to_csv("results/canonical_membership.csv", index=False)
FAMSUM = (MEM.groupby("family")
          .agg(n_genes=("gene", "nunique"), n_up=("direction", lambda s: int((s == "risk-increasing").sum())),
               n_validated=("validated", "sum"), n_screened=("screen_called", "sum"),
               family_size=("family_size", "first"), p=("family_p", "first"),
               genes=("gene", lambda s: "|".join(sorted(set(s)))))
          .reset_index().sort_values("p"))
FAMSUM["fdr"] = bh(FAMSUM.p.values)
FAMSUM.to_csv("results/canonical_membership_summary.csv", index=False)
print("\nC - risk genes literally inside each canonical CRC family:")
print(FAMSUM[["family", "n_genes", "n_up", "family_size", "p", "n_validated",
              "n_screened"]].to_string(index=False))

# ---------------------------------------------------------------- D network proximity
G = pickle.load(open("resources/ppi_graph.pkl", "rb"))
nodes = sorted(G.nodes())
nidx = {n: i for i, n in enumerate(nodes)}
deg = np.array([G.degree(n) for n in nodes], float)
dedges = np.unique(np.quantile(deg, np.linspace(0, 1, 21)))
dbin = np.clip(np.searchsorted(dedges, deg, side="right") - 1, 0, len(dedges) - 2)
dbins = [np.where(dbin == b)[0] for b in range(len(dedges) - 1)]

prox_rows = []
NPERM_D = 2000
for f, gs in fam_union.items():
    members = sorted(gs & set(nodes))
    if len(members) < 10:
        continue
    d = np.full(len(nodes), np.inf, np.float32)
    for m in members:
        for n, dd in nx.single_source_shortest_path_length(G, m, cutoff=4).items():
            j = nidx[n]
            if dd < d[j]:
                d[j] = dd
    # degree-matched null distribution of the distance, per degree bin
    null_by_bin = {}
    for b in range(len(dbins)):
        pool = dbins[b]
        v = d[pool]
        v = v[np.isfinite(v)]
        null_by_bin[b] = (float(v.mean()), float(v.std(ddof=1)) if len(v) > 2 else np.nan)
    for e in sorted(ALL):
        if e not in nidx:
            continue
        j = nidx[e]
        if not np.isfinite(d[j]):
            continue
        mu, sd = null_by_bin[dbin[j]]
        z = (d[j] - mu) / sd if sd and not np.isnan(sd) else np.nan
        r = gene[gene.entrez == e].iloc[0]
        prox_rows.append({"gene": ent2sym.get(e, str(e)), "family": f,
                          "in_family": bool(e in gs), "distance": float(d[j]),
                          "z": float(z), "direction": r.direction,
                          "degree": int(G.degree(e))})
PROX = pd.DataFrame(prox_rows)
PROX.to_csv("results/gene_pathway_proximity.csv", index=False)
CLOSE = PROX[(~PROX.in_family) & (PROX.z < -1.645)]
print(f"\nD - {CLOSE.gene.nunique()} risk genes are closer to a canonical family than "
      f"degree-matched genes are, without being annotated to it "
      f"({len(CLOSE)} gene-family pairs)")
print(CLOSE.groupby("family").agg(n=("gene", "nunique"),
                                  genes=("gene", lambda s: "|".join(sorted(s)[:8])))
      .to_string())

# ---------------------------------------------------------------- E aetiology themes
theme_union = {t: set().union(*[sets[i] for i in np.where(THEME == t)[0]])
               for t in THEMES if (THEME == t).any()}
rows, assign = [], {}
for t, gs in theme_union.items():
    hits = sorted(ALL & gs)
    k, n = len(hits), len(ALL & set(universe))
    p = hypergeom.sf(k - 1, U, len(gs), n)
    syms = [ent2sym.get(e, str(e)) for e in hits]
    dirs = [gene.loc[gene.entrez == e, "direction"].iat[0] for e in hits]
    up = sum(d == "risk-increasing" for d in dirs)
    bt = binomtest(up, len(hits), 0.5) if hits else None
    rows.append({"theme": t, "n_sets": int((THEME == t).sum()), "theme_size": len(gs),
                 "n_risk_genes": k, "expected": len(gs) * n / U,
                 "fold": k / max(len(gs) * n / U, 1e-9), "p": float(p),
                 "n_up": up, "n_down": len(hits) - up,
                 "p_direction": float(bt.pvalue) if bt else np.nan,
                 "n_validated": sum(s in validated for s in syms),
                 "n_screened": sum(s in called for s in syms),
                 "genes": "|".join(syms)})
    for s in syms:
        assign.setdefault(s, []).append(t)
TH = pd.DataFrame(rows)
TH["fdr"] = bh(TH.p.values)
TH = TH.sort_values("p")
TH.to_csv("results/etiology_themes.csv", index=False)
print("\nE - aetiology themes over all 347 risk genes:")
print(TH[["theme", "theme_size", "n_risk_genes", "expected", "fold", "p", "fdr",
          "n_up", "n_down", "p_direction"]].round(3).to_string(index=False))

AS = gene[["gene", "direction", "tissues", "n_ct", "p_best"]].copy()
AS["themes"] = AS.gene.map(lambda g: "|".join(assign.get(g, [])))
AS["canonical_families"] = AS.gene.map(
    lambda g: "|".join(sorted(MEM.loc[MEM.gene == g, "family"].unique())))
AS["network_close_to"] = AS.gene.map(
    lambda g: "|".join(sorted(CLOSE.loc[CLOSE.gene == g, "family"].unique())))
AS["validated"] = AS.gene.isin(validated)
AS["screen_called"] = AS.gene.isin(called)
AS["in_pathway_universe"] = AS.gene.map(lambda g: sym2ent.get(g, -1) in set(universe))
AS["status"] = np.where(AS.canonical_families != "", "in a canonical CRC family",
                np.where(AS.themes != "", "aetiology theme only",
                 np.where(AS.network_close_to != "", "network-adjacent only",
                  np.where(AS.in_pathway_universe, "annotated, no theme", "unannotated"))))
AS.to_csv("results/etiology_gene_assignment.csv", index=False)
print("\nwhere the 347 risk genes sit:")
print(AS.status.value_counts().to_string())

J = {
    "n_genes": len(gene), "n_up": int((gene.direction == "risk-increasing").sum()),
    "n_down": int((gene.direction == "protective").sum()),
    "nulls": [nA, nU, nD],
    "all_top": EA.nsmallest(20, "p")[["pathway", "source", "family", "theme", "overlap",
                                      "n_pathway", "p", "fdr", "genes"]].to_dict("records"),
    "all_n_sig": int((EA.fdr < 0.1).sum()),
    "up_top": EU.nsmallest(10, "p")[["pathway", "family", "overlap", "p", "fdr",
                                     "genes"]].to_dict("records"),
    "down_top": ED.nsmallest(10, "p")[["pathway", "family", "overlap", "p", "fdr",
                                       "genes"]].to_dict("records"),
    "up_n_sig": int((EU.fdr < 0.1).sum()), "down_n_sig": int((ED.fdr < 0.1).sum()),
    "families": FAMSUM.to_dict("records"),
    "n_in_any_family": int(MEM.gene.nunique()),
    "themes": TH.to_dict("records"),
    "n_in_any_theme": int(sum(1 for g in assign)),
    "proximity": {"n_close": int(CLOSE.gene.nunique()), "n_pairs": int(len(CLOSE)),
                  "by_family": CLOSE.groupby("family").gene.nunique().to_dict(),
                  "top": PROX[~PROX.in_family].nsmallest(15, "z")[
                      ["gene", "family", "distance", "z", "direction"]].to_dict("records")},
    "status": AS.status.value_counts().to_dict(),
}
json.dump(J, open("results/all347_etiology.json", "w"), indent=1, default=str)
print("\nA - top pathways over all 347:")
print(EA.nsmallest(12, "p")[["pathway", "family", "theme", "overlap", "p", "fdr",
                             "genes"]].to_string(index=False))
