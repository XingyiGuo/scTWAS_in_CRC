"""Cell-type specificity of the 546 scTWAS associations (Results paragraph 1, Figure 4A).

Counts how many of the 347 distinct risk genes are detected in only one cell-type model and
compares that to a permutation null in which each model draws its observed number of genes at
random, without replacement, from the union of all detected genes.

The null deliberately ignores that a gene can only be detected in a cell type where it is
expressed and powered, so it OVERSTATES the overlap expected by chance; the observed excess of
model-unique genes is therefore a conservative statement of cell-type specificity.

Input : the collected scTWAS result table, one row per gene x cell-type association, with
        columns Dataset, CellType, gene_name (546 rows). The Dataset/CellType/gene_name triples
        are identical in the DE-support copy under the same data folder, either may be used.
Output: ../table/T5_gene_celltype_specificity.csv  (per gene: models detected in)
        ../table/T6_celltype_specificity_permutation.csv  (observed vs null summary)
"""
import os
import numpy as np
import pandas as pd

IN = os.environ.get(
    "TF_TWAS_RESULTS",
    r"W:\VUMC_projects\scTWAS\VUMC_HTAN_scTWAS_genes_DE_support\data\scTWAS_ADSER_NL_collected_results.csv",
)
OUT = os.environ.get(
    "TF_SPECIFICITY_OUT", r"W:\VUMC_projects\scTWAS\TF_genes_pathways\Figure4_TF_pathways\table"
)
N_PERM = 10000
SEED = 0

df = pd.read_csv(IN)[["Dataset", "CellType", "gene_name"]].drop_duplicates()
df["model"] = df.Dataset + "_" + df.CellType
os.makedirs(OUT, exist_ok=True)

per = df.groupby("gene_name").agg(n_models_detected=("model", "nunique"),
                                  n_cell_types=("CellType", "nunique"))
per["models"] = df.groupby("gene_name").model.apply(lambda s: ";".join(sorted(s)))
per.reset_index().rename(columns={"gene_name": "gene"}).sort_values(
    ["n_models_detected", "gene"], ascending=[False, True]).to_csv(
    os.path.join(OUT, "T5_gene_celltype_specificity.csv"), index=False)

# per-model counts, and the polyp split into the seven shared populations vs the two
# lesion-restricted states (ASC, SSC)
SHARED7 = ["ABS", "CT", "EE", "GOB", "STM", "TAC", "TUF"]
by_model = df.groupby(["Dataset", "CellType"]).agg(n_associations=("gene_name", "size"),
                                                   n_genes=("gene_name", "nunique")).reset_index()
by_model.to_csv(os.path.join(OUT, "T7_associations_by_model.csv"), index=False)

pol = df[df.Dataset == "ADSER"]
groups = {"NL_7_populations": df[df.Dataset == "NL"],
          "polyp_7_shared_populations": pol[pol.CellType.isin(SHARED7)],
          "polyp_ASC": pol[pol.CellType == "ASC"],
          "polyp_SSC": pol[pol.CellType == "SSC"]}
split = {f"{k}_{what}": (len(v) if what == "assoc" else v.gene_name.nunique())
         for k, v in groups.items() for what in ("assoc", "genes")}
asc, ssc = set(groups["polyp_ASC"].gene_name), set(groups["polyp_SSC"].gene_name)
other = set(groups["polyp_7_shared_populations"].gene_name) | set(groups["NL_7_populations"].gene_name)
split["ASC_genes_in_no_other_model"] = len(asc - ssc - other)
split["SSC_genes_in_no_other_model"] = len(ssc - asc - other)

pool = df.gene_name.unique()
sizes = df.groupby("model").gene_name.nunique().to_dict()
obs = int((per.n_models_detected == 1).sum())

rng = np.random.default_rng(SEED)
sim = np.empty(N_PERM, dtype=int)
for i in range(N_PERM):
    cnt = {}
    for k in sizes.values():
        for g in rng.choice(len(pool), size=k, replace=False):
            cnt[g] = cnt.get(g, 0) + 1
    sim[i] = sum(1 for v in cnt.values() if v == 1)

p = (sim >= obs).mean()
pd.DataFrame([dict(
    n_pairs=len(df), n_genes=len(pool), n_models=len(sizes),
    genes_in_one_model=obs, pct_one_model=round(100 * obs / len(pool), 1),
    genes_in_one_celltype=int((per.n_cell_types == 1).sum()),
    mean_models_per_gene=round(float(per.n_models_detected.mean()), 2),
    max_models_per_gene=int(per.n_models_detected.max()),
    null_mean_one_model=round(float(sim.mean()), 1), null_sd=round(float(sim.std()), 1),
    null_ci_lo=float(np.percentile(sim, 2.5)), null_ci_hi=float(np.percentile(sim, 97.5)),
    n_permutations=N_PERM, empirical_P=f"<{1 / N_PERM:g}" if p == 0 else p,
    NL_pairs=int((df.Dataset == "NL").sum()),
    NL_genes=df[df.Dataset == "NL"].gene_name.nunique(),
    polyp_pairs=int((df.Dataset == "ADSER").sum()),
    polyp_genes=df[df.Dataset == "ADSER"].gene_name.nunique(),
    genes_shared_NL_polyp=len(set(df[df.Dataset == "NL"].gene_name) &
                              set(df[df.Dataset == "ADSER"].gene_name)),
    **split)]).to_csv(
    os.path.join(OUT, "T6_celltype_specificity_permutation.csv"), index=False)
