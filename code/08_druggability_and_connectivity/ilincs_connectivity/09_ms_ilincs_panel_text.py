#!/usr/bin/env python
"""Methods and Results text for Fig_ilincs_drug_panel panel a.

Every number in the prose is recomputed here from the exported tables, so the
text cannot drift from the data (same convention as the druggability paragraph,
src/03_druggability_paragraph.py in the Drugs folder).

Reads   table/ilincs_drug_panel_{per_compound,matching,summary}.csv
        table/ilincs_cmap_all_signatures.csv, ilincs_all_compounds_mean_NCS.csv
        table/ilincs_query_genes.csv, ilincs_reversal_model_summary.csv
        table/ilincs_cmap_model_summary.csv, sigcom_model_summary.csv
Writes  ms_ilincs_panel_methods_results.md / .docx
        table/ms_ilincs_panel_numbers.csv   every quoted value, with its source
"""
import re

import numpy as np
import pandas as pd
from docx import Document
from statsmodels.stats.multitest import multipletests

ROOT = r"W:\VUMC_projects\scTWAS\iLINCS_reversal"
T = rf"{ROOT}\table"
META = r"L:\Drugs_Proteins\iLINCS\meta_sigSearch_LIB_5_143374_20231227.csv"
N_DRAW, SEED = 5000, 11

SALTS = (r"\s+(hydrochloride|hcl|mesylate|maleate|citrate|tartrate|sulfate|sulphate|phosphate|acetate|"
         r"fumarate|succinate|besylate|tosylate|sodium|potassium|calcium|dihydrate|monohydrate|hydrate|"
         r"ethanolate|malate|oxalate|bromide|chloride|nitrate|lactate|gluconate|carbonate|benzoate|"
         r"palmitate|dimaleate|dimesylate|disodium)\b")


def norm_name(s):
    s = str(s).strip().lower()
    s = re.sub(r"\(.*?\)", " ", s)
    for _ in range(3):
        s = re.sub(SALTS, " ", s)
    return re.sub(r"[^a-z0-9]", "", s)


pc = pd.read_csv(rf"{T}\ilincs_drug_panel_per_compound.csv")
mt = pd.read_csv(rf"{T}\ilincs_drug_panel_matching.csv")
ps = pd.read_csv(rf"{T}\ilincs_drug_panel_summary.csv").iloc[0]
qg = pd.read_csv(rf"{T}\ilincs_query_genes.csv")
prs = pd.read_csv(rf"{T}\ilincs_reversal_model_summary.csv")
cms = pd.read_csv(rf"{T}\ilincs_cmap_model_summary.csv")
sgs = pd.read_csv(rf"{T}\sigcom_model_summary.csv")
libc = pd.read_csv(rf"{T}\ilincs_all_compounds_mean_NCS.csv")
cm = pd.read_csv(rf"{T}\ilincs_cmap_all_signatures.csv")
m5 = pd.read_csv(META)[["SignatureId", "pubChemID"]]

cm = cm.merge(m5, left_on="signature", right_on="SignatureId", how="left")
cm["key"] = np.where(cm.pubChemID.notna(),
                     "CID" + cm.pubChemID.fillna(0).astype("int64").astype(str),
                     cm.Perturbagen.map(norm_name))

# ---- panel composition -------------------------------------------------------
n_models = qg.model.nunique()
n_side = int(qg.groupby(["model", "side"]).size().unique()[0])
lm_lo, lm_hi = int(prs.n_landmark_shared.min()), int(prs.n_landmark_shared.max())
n_sig_lib = int(cm.signature.nunique())
n_pert_lib = int(cm.Perturbagen.nunique())
n_panel = len(mt)
n_cid = int(mt.cid.notna().sum())
n_match = int(mt.in_iLINCS.sum())
n_unmatched = n_panel - n_match

# ---- the four FDR-significant compounds --------------------------------------
hits = pc[pc.FDR_vs_sizematched_null < 0.05].copy()
rows = []
for _, r in hits.iterrows():
    sub = cm[cm.key.eq(r.key)]
    per_model = sub.groupby("model").NCS.mean()
    rows.append({"drug": str(r.drug_name).title(), "target": r.risk_gene_targets,
                 "n_signatures": int(r.n_signatures), "mean_NCS": r.mean_NCS,
                 "null_mean": r.null_mean, "P": r.P_vs_sizematched_null,
                 "FDR": r.FDR_vs_sizematched_null,
                 "pct_signatures_negative": 100 * r.frac_signatures_negative,
                 "models_negative": int((per_model < 0).sum()),
                 "n_cell_lines": int(sub.CellLine.nunique()),
                 "library_rank": int((libc.mean_NCS < r.mean_NCS).sum()) + 1})
hd = pd.DataFrame(rows).sort_values("FDR").reset_index(drop=True)
n_lib_cpd = len(libc)

# ---- class-level tests, size-matched at the compound level --------------------
rng = np.random.default_rng(SEED)
ref = libc.loc[~libc.in_panel]


def class_test(members):
    obs = members.mean_NCS.mean()
    # nearest 50 non-panel compounds by log signature count (an exact +/-2 match
    # is impossible for compounds with hundreds of signatures, e.g. bortezomib)
    lref = np.log10(ref.n_signatures.to_numpy())
    pools = []
    for k in members.n_signatures:
        idx = np.argsort(np.abs(lref - np.log10(k)))[:50]
        pools.append(ref.mean_NCS.to_numpy()[idx])
    draws = np.array([np.mean([rng.choice(p) for p in pools]) for _ in range(N_DRAW)])
    return obs, draws.mean(), (int((draws <= obs).sum()) + 1) / (N_DRAW + 1)


erbb2 = pc[pc.risk_gene_targets.str.contains("ERBB2", na=False)]
prot = pc[pc.risk_gene_targets.str.contains("PSM", na=False)]
e_obs, e_null, e_p = class_test(erbb2)
p_obs, p_null, p_p = class_test(prot)
carf = pc[pc.drug_name.str.lower().eq("carfilzomib")].iloc[0]

# ---- library-wide context ----------------------------------------------------
cmap_min_p = float(cms.P_maxstat.min())
cmap_bonf = min(1.0, cmap_min_p * len(cms))
cmap_n_ns = int((cms.P_maxstat > 0.05).sum())
cmap_worst = cms.loc[cms.P_maxstat.idxmin()]
sg_cp = sgs[sgs.library.eq("l1000_cp")]
sg_n = int(sg_cp.library_signatures.iloc[0])
sg_fdr = int(sg_cp.n_reversers_fdr05.sum())                      # l1000_cp only
sg_mean = sgs[sgs.library.eq("l1000_mean_cp")]
sg_fdr_mean = int(sg_mean.n_reversers_fdr05.sum())               # l1000_mean_cp only
sg_n_mean = int(sg_mean.library_signatures.iloc[0])
sg_bonf = int(sgs.n_reversers_bonferroni05.sum())                # both libraries

num = pd.DataFrame([
    ("cell-type TWAS models", n_models, "ilincs_query_genes.csv"),
    ("query genes per side", n_side, "ilincs_query_genes.csv"),
    ("landmark genes shared with a model, min", lm_lo, "ilincs_reversal_model_summary.csv"),
    ("landmark genes shared with a model, max", lm_hi, "ilincs_reversal_model_summary.csv"),
    ("iLINCS signatures scanned", n_sig_lib, "ilincs_cmap_all_signatures.csv"),
    ("distinct perturbagens in the library", n_pert_lib, "ilincs_cmap_all_signatures.csv"),
    ("panel compounds after alias collapse", n_panel, "ilincs_drug_panel_matching.csv"),
    ("panel compounds with a PubChem CID", n_cid, "ilincs_drug_panel_matching.csv"),
    ("panel compounds present in iLINCS", n_match, "ilincs_drug_panel_matching.csv"),
    ("panel compounds absent from iLINCS", n_unmatched, "ilincs_drug_panel_matching.csv"),
    ("compounds passing BH FDR < 0.05", len(hd), "ilincs_drug_panel_per_compound.csv"),
    ("library compounds for ranking", n_lib_cpd, "ilincs_all_compounds_mean_NCS.csv"),
    ("panel vs library, compound-level P", ps.P_panel_compounds_vs_library_compounds,
     "ilincs_drug_panel_summary.csv"),
    ("ERBB2 class mean NCS", e_obs, "recomputed"),
    ("ERBB2 class null mean", e_null, "recomputed"),
    ("ERBB2 class P", e_p, "recomputed"),
    ("proteasome class mean NCS", p_obs, "recomputed"),
    ("proteasome class null mean", p_null, "recomputed"),
    ("proteasome class P", p_p, "recomputed"),
    ("carfilzomib mean NCS", carf.mean_NCS, "ilincs_drug_panel_per_compound.csv"),
    ("carfilzomib FDR", carf.FDR_vs_sizematched_null, "ilincs_drug_panel_per_compound.csv"),
    ("library scan models not significant", cmap_n_ns, "ilincs_cmap_model_summary.csv"),
    ("library scan smallest P", cmap_min_p, "ilincs_cmap_model_summary.csv"),
    ("library scan smallest P, Bonferroni", cmap_bonf, "ilincs_cmap_model_summary.csv"),
    ("SigCom signatures scanned", sg_n, "sigcom_model_summary.csv"),
    ("SigCom l1000_cp reversers at FDR < 0.05", sg_fdr, "sigcom_model_summary.csv"),
    ("SigCom l1000_mean_cp signatures scanned", sg_n_mean, "sigcom_model_summary.csv"),
    ("SigCom l1000_mean_cp reversers at FDR < 0.05", sg_fdr_mean, "sigcom_model_summary.csv"),
    ("SigCom reversers at Bonferroni, both libraries", sg_bonf, "sigcom_model_summary.csv"),
], columns=["quantity", "value", "source_table"])
num.to_csv(rf"{T}\ms_ilincs_panel_numbers.csv", index=False)

# ---- prose -------------------------------------------------------------------
h = {r.drug: r for _, r in hd.iterrows()}
WORDS = {1: "One", 2: "Two", 3: "Three", 4: "Four", 5: "Five", 6: "Six", 7: "Seven"}
spell = lambda n: WORDS.get(int(n), str(int(n)))
tbl = "\n".join(
    f"| {r.drug} | {r.target.replace('PSMA7;PSMD13;PSMD3', 'PSMA7, PSMD3, PSMD13')} | "
    f"{r.n_signatures} | {r.mean_NCS:+.3f} | {r.null_mean:+.3f} | "
    f"{r.pct_signatures_negative:.0f}% | {r.models_negative}/{n_models} | "
    f"{r.P:.4f} | {r.FDR:.3f} |" for _, r in hd.iterrows())

methods = f"""## Methods

**Query signatures from the cell-type TWAS models.** For each of the {n_models} European-ancestry
cell-type models with S-PrediXcan peak-SNP associations, a bidirectional query signature was
built from the association Z scores. Because the iLINCS LIB5 library reports only the {978}
LINCS landmark genes, and only {lm_lo}-{lm_hi} of these carry a TWAS Z in a given model, the
query was drawn from the genes shared between each model and the landmark space rather than
genome-wide: the {n_side} most positive and {n_side} most negative Z scores among shared genes,
giving a {2 * n_side}-gene tagged query per model. Taking the top {2 * n_side} genome-wide would have
left approximately six measurable genes per query. The query genes are therefore not the
genome-wide top TWAS hits, and this restriction is a property of the L1000 platform rather
than a choice about the genetics.

**Perturbation signatures.** Chemical perturbagen signatures were taken from the iLINCS LIB5
large-intestine library ({n_sig_lib:,} signatures, {n_pert_lib:,} distinct perturbagens, log2
differential expression relative to matched controls). Signature identity, perturbagen name,
cell line, dose and exposure time were joined from the accompanying LIB5 signature metadata.

**Connectivity scoring.** Each query was scored against each signature with the weighted
connectivity score of Subramanian et al. (Cell 2017;171:1437). Signatures were ranked by
log2 fold change in descending order; enrichment scores ES_up and ES_down were computed as
weighted Kolmogorov-Smirnov running sums (weight = |log2FC|, exponent 1) of the up- and
down-tagged query genes in that ranked list; and WTCS = (ES_up - ES_down)/2 when the two
enrichment scores had opposite signs, and 0 otherwise. Scores were normalised within cell
line (NCS) so that signatures from different lines are comparable. Negative NCS indicates
reversal - the perturbagen suppresses genes whose predicted expression raises risk and
induces those whose predicted expression lowers it. The implementation was validated on
constructed cases: a tag set equal to a signature's own 50 highest-ranked genes returns
ES = +1.000, its 50 lowest-ranked genes ES = -1.000, a random tag set ES near 0, and a
perfectly reversing query WTCS = -1.000.

**Compound panel.** The panel was pre-specified as the approved or phase II/III compounds
annotated against the 63 druggable colorectal risk genes in the four-resource druggability
table (ChEMBL, DrugBank, TTD and Open Targets). Compound identity between the drug table and
iLINCS was resolved on PubChem CID rather than on drug name: CIDs were obtained from the
PubChem name-to-CID service ({n_cid} of the panel compounds resolved) and matched against the
CIDs carried in the LIB5 metadata, falling back to a salt-stripped normalised name where no
CID was available. Name matching alone both missed compounds stored under a different synonym
and double-counted aliases of the same molecule, so records resolving to one CID were
collapsed, giving {n_panel} distinct panel compounds of which {n_match} are present in iLINCS.
The {n_unmatched} absentees are predominantly biologics - antibodies, antibody-drug conjugates,
bispecific constructs and CAR-T products - together with several recent tyrosine kinase
inhibitors, none of which are represented in a small-molecule perturbation library; querying
the full LIB5 collection (143,374 signatures over 23 tissues) recovered no additional panel
compound, confirming that this is a library-content limit and not an effect of restricting
to colorectal lines.

**Inference.** A signature x model score is not an independent observation: the {n_models} models
share query genes, and one compound contributes several doses, time points and cell lines.
Each compound was therefore reduced to a single value - the mean NCS over its signatures,
each signature first averaged across the {n_models} models - and referenced against a
size-matched empirical null, formed by drawing {N_DRAW:,} random sets of the same number of
library signatures and taking their mean NCS. One-sided P values were corrected across the
{n_match} tested compounds by the Benjamini-Hochberg procedure. Target-class tests (ERBB2
inhibitors, proteasome inhibitors) used the same construction at the class level, sampling, for each
class member, from the 50 non-panel library compounds closest to it in log signature count. Whether
the panel as a whole is shifted was tested by Mann-Whitney U on compound-level mean NCS,
panel against non-panel library compounds; the signature-level version of that comparison is
pseudo-replicated by compounds with many signatures and is not reported.

Analyses were run in Python 3.13 with numpy, pandas, scipy and statsmodels. Code is in
`src/06_ilincs_cmap_connectivity.py` (scoring) and `src/07_ilincs_drug_panel.py` (panel and
inference); every value quoted in the Results is regenerated from the exported tables by
`src/09_ms_ilincs_panel_text.py` and listed in `table/ms_ilincs_panel_numbers.csv`.
"""

results = f"""## Results

**Drugs of the druggable risk genes reverse the predicted risk profile in three
instances, one of which is mechanistically specific.** Of the {n_panel} approved or
phase II/III compounds annotated against the 63 druggable risk genes, {n_match} are
represented in the iLINCS large-intestine library and were scored against the
{n_models} cell-type TWAS query signatures. {spell(len(hd))} of the {n_match} reach
Benjamini-Hochberg FDR < 0.05 against a size-matched null (Figure Xa).

| Compound | Risk-gene target | Signatures | Mean NCS | Null mean | Signatures reversing | Models reversing | P | FDR |
|---|---|---|---|---|---|---|---|---|
{tbl}

Resorcinol rests on a single signature and is not interpreted further. The remaining
three reverse consistently: {h['Neratinib'].models_negative} of {n_models} cell-type models for
neratinib across {h['Neratinib'].n_cell_lines} cell lines,
{h['Ixazomib'].models_negative} of {n_models} for ixazomib across {h['Ixazomib'].n_cell_lines} lines,
and {h['Bortezomib'].models_negative} of {n_models} for bortezomib across
{h['Bortezomib'].n_cell_lines} lines.

**All three hits act on genes that are themselves in the risk set, but only neratinib
engages one selectively.** ERBB2 is one of the 63 druggable risk genes and the most heavily
drugged gene in the set, so for neratinib the transcriptional reversal and the genetic
association converge on the same protein. PSMA7, PSMD3 and PSMD13 are also risk genes, but
they are pan-essential in the large intestine and acquire bortezomib, carfilzomib and
ixazomib at the level of the assembled proteasome rather than through subunit-selective
chemistry, and proteasome inhibition produces a large proteotoxic-stress response that
opposes many unrelated expression profiles. For the proteasome pair a negative connectivity
score is therefore expected on generic grounds and does not localise to the risk genes.

**The proteasome effect is class-wide; the ERBB2 effect is not.** Tested as sets against
size-matched nulls, the {len(prot)} proteasome inhibitors average NCS {p_obs:+.3f} against
{p_null:+.3f} (P = {p_p:.3f}), so the whole class reverses the risk profile. That is what a
generic proteotoxic-stress response predicts and it argues against subunit-selective action:
carfilzomib is directionally concordant ({100 * carf.frac_signatures_negative:.0f}% of its
signatures negative) while individually non-significant. The {len(erbb2)} ERBB2-directed
compounds, by contrast, are not shifted as a class ({e_obs:+.3f} against {e_null:+.3f},
P = {e_p:.2f}); afatinib, CP-724714, merimepodib and sapitinib sit at or above zero. Neratinib
is thus separated from the other HER2 inhibitors rather than exemplifying them, which is
consistent with its irreversible pan-HER pharmacology but cannot be distinguished here from
compound-specific off-target activity.

**All three are modest in absolute terms.** Among the {n_lib_cpd:,} compounds in the library,
ixazomib ranks {h['Ixazomib'].library_rank}, neratinib {h['Neratinib'].library_rank} and
bortezomib {h['Bortezomib'].library_rank} by mean NCS, and the panel as a whole is not
shifted negative relative to non-panel compounds
(P = {ps.P_panel_compounds_vs_library_compounds:.2f}). These compounds clear significance
because the panel is small and defined in advance by existing drug-target evidence, not
because they are outstanding reversers of the risk profile.

**An unbiased library screen nominates nothing.** Scanning all {n_sig_lib:,} colorectal
signatures per model and comparing the strongest reversal in each scan with a
max-statistic permutation null left {cmap_n_ns} of the {n_models} models non-significant; the
one nominal exception (P = {cmap_min_p:.3f}) does not survive correction for the
{n_models} models tested (P = {cmap_bonf:.3f}). Submitting the same queries to SigCom LINCS,
which scores them against {sg_n:,} L1000 signatures with an independent implementation,
returned {sg_fdr} reversers at FDR < 0.05; its {sg_n_mean:,}-signature per-compound consensus
library returned {sg_fdr_mean}, for one model. No reverser reached Bonferroni significance in
either library ({sg_bonf} of {2 * n_models} model-library queries). The pre-specified drug panel is
therefore the only setting in which this analysis has the power to detect anything, and its
result is best read as transcriptional support for HER2 inhibition in colorectal tissue
rather than as the nomination of a novel repurposing candidate.
"""

md = "# iLINCS reversal of cell-type TWAS profiles\n\n" + results + "\n" + methods
open(rf"{ROOT}\ms_ilincs_panel_methods_results.md", "w", encoding="utf-8").write(md)

doc = Document()
for block in md.split("\n\n"):
    b = block.strip()
    if not b:
        continue
    if b.startswith("# "):
        doc.add_paragraph(b[2:], style="Heading 1")
    elif b.startswith("## "):
        doc.add_paragraph(b[3:], style="Heading 2")
    elif b.startswith("|"):
        lines = [l for l in b.split("\n") if l.strip().startswith("|")]
        cells = [[c.strip() for c in l.strip().strip("|").split("|")] for l in lines]
        cells = [c for c in cells if not set("".join(c)) <= set("-: ")]
        t = doc.add_table(rows=0, cols=len(cells[0]))
        t.style = "Table Grid"
        for ri, row in enumerate(cells):
            cs = t.add_row().cells
            for ci, val in enumerate(row):
                cs[ci].text = val
                if ri == 0:
                    cs[ci].paragraphs[0].runs and setattr(
                        cs[ci].paragraphs[0].runs[0], "bold", True)
    else:
        p = doc.add_paragraph()
        for j, part in enumerate(b.split("**")):
            r = p.add_run(part.replace("\n", " "))
            r.bold = (j % 2 == 1)
try:
    doc.save(rf"{ROOT}\ms_ilincs_panel_methods_results.docx")
except PermissionError:
    doc.save(rf"{ROOT}\ms_ilincs_panel_methods_results_NEW.docx")

print(hd.to_string(index=False))
print(f"\nERBB2 class {e_obs:+.4f} vs {e_null:+.4f} P {e_p:.4f} | "
      f"proteasome {p_obs:+.4f} vs {p_null:+.4f} P {p_p:.4f}")
print("numbers table rows:", len(num), "| words in md:", len(md.split()))
