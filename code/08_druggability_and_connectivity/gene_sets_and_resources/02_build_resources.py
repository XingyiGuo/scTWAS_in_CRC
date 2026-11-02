#!/usr/bin/env python3
"""
02 - Build gene ID maps, the multi-layer molecular network, and drug resources.

Resources (all obtained as plain files from public GitHub repositories):
  Hetionet v1.0                (hetio/hetionet, dhimmel/integrate)  - nodes, edges
  STRING v10 high-confidence   (hart-lab/bagel Network_STRING_v10)  - PPI
  DrugBank 5.0 derived tables  (dhimmel/drugbank)                   - drugs, targets, chem similarity
  BindingDB derived affinities (dhimmel/bindingdb)
  Pathway gene sets            (dhimmel/pathways: Reactome/WikiPathways/PID)
  LINCS L1000 consensus sigs   (dhimmel/lincs)                      - knockdown / overexpression / drug
  BAGEL reference gene sets    (hart-lab/bagel CEGv2 / NEGv1)
  GDSC IC50 + CCLE annotation  (kimmo1019/DeepCDR)

Outputs -> resources/
"""
import pandas as pd, numpy as np, gzip, bz2, os, json, pickle
import networkx as nx

R = "resources"
os.makedirs(R, exist_ok=True)
HET = "/tmp/het"; INT = "/tmp/integrate"; DB = "/tmp/drugbank"; BD = "/tmp/bindingdb"
PW = "/tmp/pathways"; LI = "/tmp/lincs"; BG = "/tmp/bagel"; CDR = "/tmp/DeepCDR"

# ---------------------------------------------------------------- gene id map
nodes = pd.read_csv(f"{HET}/hetnet/tsv/hetionet-v1.0-nodes.tsv", sep="\t")
genes = nodes[nodes["kind"] == "Gene"].copy()
genes["entrez"] = genes["id"].str.replace("Gene::", "", regex=False).astype(int)
sym2ent = dict(zip(genes["name"], genes["entrez"]))
ent2sym = dict(zip(genes["entrez"], genes["name"]))

li_gi = pd.read_csv(f"{LI}/data/geneinfo/geneinfo.tsv.gz", sep="\t")
li_gi = li_gi.dropna(subset=["pr_gene_symbol"])
for s, e in zip(li_gi["pr_gene_symbol"], li_gi["pr_gene_id"]):
    e = str(e)
    if "///" in e:
        continue
    try:
        e = int(e)
    except Exception:
        continue
    sym2ent.setdefault(s, e); ent2sym.setdefault(e, s)

cg = pd.read_csv(f"{LI}/data/consensi/genes.tsv", sep="\t")
for s, e in zip(cg["symbol"], cg["entrez_gene_id"]):
    sym2ent.setdefault(s, int(e)); ent2sym.setdefault(int(e), s)

json.dump({"sym2ent": {k: int(v) for k, v in sym2ent.items()}}, open(f"{R}/gene_id_map.json", "w"))
print(f"gene id map: {len(sym2ent)} symbols")

# ---------------------------------------------------------------- network
edges = pd.read_csv(f"{INT}/data/edges.sif.gz", sep="\t")
print("hetnet edges:", edges.shape)


def entrez(x):
    return int(x.split("::")[1]) if x.startswith("Gene::") else None


G = nx.Graph()        # undirected functional/physical layer
DG = nx.DiGraph()     # directed regulatory layer

ppi = edges[edges["metaedge"].isin(["GiG", "GcG"])]
for s, m, t in ppi.itertuples(index=False):
    a, b = entrez(s), entrez(t)
    if a and b and a != b:
        w = 0.9 if m == "GiG" else 0.5
        if G.has_edge(a, b):
            G[a][b]["weight"] = max(G[a][b]["weight"], w)
        else:
            G.add_edge(a, b, weight=w)
print(f"after hetionet GiG/GcG: {G.number_of_nodes()} nodes, {G.number_of_edges()} edges")

# STRING v10, high confidence only (score >= 700)
st = pd.read_csv(f"{BG}/Network_STRING_v10.txt", sep="\t", header=None,
                 names=["a", "b", "score"])
st = st[st["score"] >= 700]
n_add = 0
for a, b, s in st.itertuples(index=False):
    ea, eb = sym2ent.get(a), sym2ent.get(b)
    if ea and eb and ea != eb:
        w = float(s) / 1000.0
        if G.has_edge(ea, eb):
            G[ea][eb]["weight"] = max(G[ea][eb]["weight"], w)
        else:
            G.add_edge(ea, eb, weight=w); n_add += 1
print(f"after STRING>=700 (+{n_add} new): {G.number_of_nodes()} nodes, {G.number_of_edges()} edges")

reg = edges[edges["metaedge"] == "Gr>G"]
for s, m, t in reg.itertuples(index=False):
    a, b = entrez(s), entrez(t)
    if a and b and a != b:
        DG.add_edge(a, b, weight=1.0)
print(f"regulatory DiGraph: {DG.number_of_nodes()} nodes, {DG.number_of_edges()} edges")

nx.write_gpickle = None
pickle.dump(G, open(f"{R}/ppi_graph.pkl", "wb"))
pickle.dump(DG, open(f"{R}/reg_graph.pkl", "wb"))

# ---------------------------------------------------------------- drugs
drugs = pd.read_csv(f"{DB}/data/drugbank.tsv", sep="\t")
prot = pd.read_csv(f"{DB}/data/proteins.tsv", sep="\t")
prot = prot[(prot["organism"] == "Human") & (prot["category"] == "target")].copy()
prot["entrez_gene_id"] = pd.to_numeric(prot["entrez_gene_id"], errors="coerce")
prot = prot.dropna(subset=["entrez_gene_id"])
prot["entrez_gene_id"] = prot["entrez_gene_id"].astype(int)
prot["gene_symbol"] = prot["entrez_gene_id"].map(ent2sym)
prot.to_csv(f"{R}/drugbank_targets.tsv", sep="\t", index=False)
drugs.to_csv(f"{R}/drugbank_drugs.tsv", sep="\t", index=False)
print(f"DrugBank: {drugs.shape[0]} drugs, {prot.shape[0]} human target links, "
      f"{prot['drugbank_id'].nunique()} drugs with targets")

bind = pd.read_csv(f"{BD}/data/bindings-drugbank-gene.tsv", sep="\t")
bind.to_csv(f"{R}/bindingdb.tsv", sep="\t", index=False)
print(f"BindingDB: {bind.shape[0]} drug-gene affinity records")

cbg = edges[edges["metaedge"] == "CbG"].copy()
cbg["entrez"] = cbg["target"].map(entrez)
cbg["drugbank_id"] = cbg["source"].str.replace("Compound::", "", regex=False)
cbg[["drugbank_id", "entrez"]].to_csv(f"{R}/hetionet_CbG.tsv", sep="\t", index=False)
print(f"Hetionet compound-binds-gene: {cbg.shape[0]}")

# disease signature for colorectal cancer from Hetionet DuG / DdG
dis = nodes[nodes["kind"] == "Disease"]
crc = dis[dis["name"].str.contains("colorectal|colon", case=False, na=False)]
print("\nCRC-related Hetionet diseases:\n", crc.to_string())
crc_ids = set(crc["id"])
dsig = edges[edges["metaedge"].isin(["DuG", "DdG"]) & edges["source"].isin(crc_ids)].copy()
dsig["entrez"] = dsig["target"].map(entrez)
dsig["dir"] = np.where(dsig["metaedge"] == "DuG", 1, -1)
dsig[["source", "entrez", "dir"]].to_csv(f"{R}/crc_disease_signature.tsv", sep="\t", index=False)
print(f"CRC disease signature genes: {dsig.shape[0]}")

dag = edges[edges["metaedge"] == "DaG"].copy()
dag = dag[dag["source"].isin(crc_ids)]
dag["entrez"] = dag["target"].map(entrez)
dag[["source", "entrez"]].to_csv(f"{R}/crc_DaG.tsv", sep="\t", index=False)

ctd = edges[edges["metaedge"].isin(["CtD", "CpD"])].copy()
ctd["drugbank_id"] = ctd["source"].str.replace("Compound::", "", regex=False)
ctd.to_csv(f"{R}/compound_treats_disease.tsv", sep="\t", index=False)
print(f"CtD/CpD edges: {ctd.shape[0]}; CRC-treating: "
      f"{ctd[ctd['target'].isin(crc_ids)].shape[0]}")

# ---------------------------------------------------------------- pathways
pw = pd.read_csv(f"{PW}/data/pathways.tsv", sep="\t")
pw = pw.dropna(subset=["coding_genes"])
pw["gene_list"] = pw["coding_genes"].astype(str).str.split("|")
pw = pw[pw["gene_list"].apply(lambda g: 5 <= len(g) <= 500)]
pw.to_pickle(f"{R}/pathways.pkl")
print(f"Pathways (5-500 genes): {pw.shape[0]} from {sorted(pw['source'].unique())}")

# ---------------------------------------------------------------- essentiality
ceg = pd.read_csv(f"{BG}/CEGv2.txt", sep="\t")["GENE"].tolist()
neg = pd.read_csv(f"{BG}/NEGv1.txt", sep="\t")["GENE"].tolist()
json.dump({"core_essential": ceg, "nonessential": neg}, open(f"{R}/essentiality_refsets.json", "w"))
print(f"Reference sets: {len(ceg)} core-essential, {len(neg)} non-essential")

print("\nDone.")
