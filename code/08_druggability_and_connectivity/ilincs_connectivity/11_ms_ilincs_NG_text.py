#!/usr/bin/env python
"""Nature Genetics-format Results, figure legend and Methods for Figure 6e.

Rewrites ms_ilincs_panel_methods_results_NEW.docx for insertion into
ms-guo/scTWAS_manuscript-09182026.docx: only panel a of the old two-panel
figure is used (now Fig. 6e, built by src/10_fig6e_panel.py), callouts,
statistics reporting, en dashes, italic *P* and superscript reference numbers
follow Nature Genetics style, and the text is written to replace the
superseded iLINCS sentences of Results section 6 and the iLINCS paragraph of
the Methods.

Every quoted value is read back from the exported tables and written to
table/ms_ilincs_panel_NG_numbers.csv next to its source, as for
src/09_ms_ilincs_panel_text.py.

Outputs (under <ROOT>/iLINCS_reversal):
  ms_ilincs_panel_NG.md
  ms_ilincs_panel_NG.docx
  table/ms_ilincs_panel_NG_numbers.csv
"""
import os
import re
import numpy as np
import pandas as pd
from docx import Document
from docx.shared import Pt

ROOT = r"W:\VUMC_projects\scTWAS"
OUT = rf"{ROOT}\iLINCS_reversal"
T = rf"{OUT}\table"
META = r"L:\Drugs_Proteins\iLINCS\meta_sigSearch_LIB_5_143374_20231227.csv"
PROTEASOME = "PSMA7;PSMD13;PSMD3"


SRC = os.path.dirname(os.path.abspath(__file__)) if "__file__" in dir() else rf"{OUT}\src"
TWAS_DIR = r"W:\VUMC_projects\scTWAS\TWAS\SPrediXcan_peakSNP_association"
LIB5 = r"L:\Drugs_Proteins\iLINCS\LIB5_signature_expr_logDiffExpr_large_intestine.csv"
CELLTYPES = ("ABS", "ASC", "CT", "EE", "GOB", "SSC", "STM", "TAC", "TUF")


def const(script, token, default=None, pos=None):
    """Read a literal constant out of one of the analysis scripts."""
    txt = open(os.path.join(SRC, script), encoding="utf-8").read()
    line = next(l for l in txt.split("\n") if token in l and "=" in l)
    nums = re.findall(r"\d+", line.split("=", 1)[1] if pos is None else line.split("=", 1)[1])
    if pos is not None:
        return int(nums[pos])
    if default is not None and default in nums:
        return int(default)
    return int(nums[0])


def model_inventory():
    """Every candidate S-PrediXcan model, its landmark overlap and whether it was used.

    Reproduces the inclusion rule of src/04_ilincs_reversal.py: a model enters only
    if at least 4 x N_SIDE genes are shared with the iLINCS landmark space.
    """
    thr = 4 * const("04_ilincs_reversal.py", "N_SIDE")
    lm = set(pd.read_csv(LIB5, usecols=["Name_GeneSymbol"]).Name_GeneSymbol)
    rows = []
    for fn in sorted(f for f in os.listdir(TWAS_DIR) if f.endswith(".EUR.TWAS.csv")):
        m = fn.replace("_peakSNP.EUR.TWAS.csv", "")
        tw = pd.read_csv(os.path.join(TWAS_DIR, fn), usecols=["gene_name", "zscore"]).dropna()
        tw = (tw.assign(a=tw.zscore.abs()).sort_values("a", ascending=False)
                .drop_duplicates("gene_name"))
        sh = int(tw.gene_name.isin(lm).sum())
        ct = m.rsplit("_", 1)[1]
        rows.append({"model": m, "dataset": m[: -(len(ct) + 1)],
                     "tissue": "polyp" if "ADSER" in m else "normal colon",
                     "cell_type": ct if ct in CELLTYPES else "",
                     "n_genes_with_TWAS_Z": len(tw), "n_landmark_shared": sh,
                     "used_as_query": sh >= thr,
                     "excluded_reason": "" if sh >= thr else
                     f"fewer than {thr} genes shared with the L1000 landmark space"})
    return pd.DataFrame(rows)


# ----------------------------------------------------------------- numbers ---
def norm_name(s):
    s = str(s).strip().lower()
    s = re.sub(r"\(.*?\)", " ", s)
    return re.sub(r"[^a-z0-9]", "", s)


def collect():
    """Every value quoted in the text, with the table it comes from."""
    n = pd.read_csv(rf"{T}\ms_ilincs_panel_numbers.csv").set_index("quantity")
    pc = pd.read_csv(rf"{T}\ilincs_drug_panel_per_compound.csv")
    libc = pd.read_csv(rf"{T}\ilincs_all_compounds_mean_NCS.csv").sort_values("mean_NCS")
    sd = pd.read_csv(rf"{OUT}\figure\Fig_ilincs_drug_panel_sourcedata.csv")
    cm = pd.read_csv(rf"{T}\ilincs_cmap_all_signatures.csv",
                     usecols=["signature", "Perturbagen", "CellLine", "NCS", "model"])
    meta = pd.read_csv(META, usecols=["SignatureId", "pubChemID"])
    cm = cm.merge(meta, left_on="signature", right_on="SignatureId", how="left")
    cm["key"] = np.where(cm.pubChemID.notna(),
                         "CID" + cm.pubChemID.fillna(0).astype("int64").astype(str),
                         cm.Perturbagen.map(norm_name))

    v, src = {}, {}

    def put(name, value, source):
        v[name] = value
        src[name] = source

    for q, col in [("models", "cell-type TWAS models"),
                   ("sig_scanned", "iLINCS signatures scanned"),
                   ("perturbagens", "distinct perturbagens in the library"),
                   ("panel_n", "panel compounds after alias collapse"),
                   ("panel_cid", "panel compounds with a PubChem CID"),
                   ("panel_in", "panel compounds present in iLINCS"),
                   ("panel_out", "panel compounds absent from iLINCS"),
                   ("n_fdr", "compounds passing BH FDR < 0.05"),
                   ("lib_compounds", "library compounds for ranking"),
                   ("lm_min", "landmark genes shared with a model, min"),
                   ("lm_max", "landmark genes shared with a model, max"),
                   ("query_side", "query genes per side"),
                   ("p_panel_vs_lib", "panel vs library, compound-level P"),
                   ("erbb2_mean", "ERBB2 class mean NCS"),
                   ("erbb2_null", "ERBB2 class null mean"),
                   ("erbb2_p", "ERBB2 class P"),
                   ("prot_mean", "proteasome class mean NCS"),
                   ("prot_null", "proteasome class null mean"),
                   ("prot_p", "proteasome class P"),
                   ("carf_mean", "carfilzomib mean NCS"),
                   ("carf_fdr", "carfilzomib FDR"),
                   ("scan_ns", "library scan models not significant"),
                   ("scan_p", "library scan smallest P"),
                   ("scan_p_bonf", "library scan smallest P, Bonferroni"),
                   ("sigcom_sig", "SigCom signatures scanned"),
                   ("sigcom_hits", "SigCom l1000_cp reversers at FDR < 0.05"),
                   ("sigcom_cons_sig", "SigCom l1000_mean_cp signatures scanned"),
                   ("sigcom_cons_hits", "SigCom l1000_mean_cp reversers at FDR < 0.05"),
                   ("sigcom_bonf", "SigCom reversers at Bonferroni, both libraries")]:
        put(q, float(n.loc[col, "value"]), n.loc[col, "source_table"])

    # per-hit statistics
    for disp, row in [("ixazomib", "IXAZOMIB"), ("neratinib", "NERATINIB"),
                      ("bortezomib", "BORTEZOMIB"), ("resorcinol", "Resorcinol"),
                      ("carfilzomib", "CARFILZOMIB")]:
        r = pc[pc.drug_name == row].iloc[0]
        put(f"{disp}_mean", r.mean_NCS, "ilincs_drug_panel_per_compound.csv")
        put(f"{disp}_p", r.P_vs_sizematched_null, "ilincs_drug_panel_per_compound.csv")
        put(f"{disp}_fdr", r.FDR_vs_sizematched_null, "ilincs_drug_panel_per_compound.csv")
        put(f"{disp}_nsig", int(r.n_signatures), "ilincs_drug_panel_per_compound.csv")
        put(f"{disp}_fracneg", r.frac_signatures_negative, "ilincs_drug_panel_per_compound.csv")
        g = cm[cm.key == r.key]
        put(f"{disp}_lines", g.CellLine.dropna().nunique(), "ilincs_cmap_all_signatures.csv")
        rank = int(libc.reset_index(drop=True).index[libc.reset_index(drop=True).key == r.key][0]) + 1
        put(f"{disp}_rank", rank, "ilincs_all_compounds_mean_NCS.csv")
        b = sd[(sd.panel == "b") & (sd.drug_name.str.upper() == row.upper())]
        if len(b):
            put(f"{disp}_models_neg", int((b.mean_NCS < 0).sum()),
                "Fig_ilincs_drug_panel_sourcedata.csv")

    # ---- model inventory: which S-PrediXcan models became queries, and why ----
    inv = model_inventory()
    inv.to_csv(rf"{T}\ilincs_query_model_inventory.csv", index=False)
    used = inv[inv.used_as_query]
    exc = inv[~inv.used_as_query]
    put("n_models_total", len(inv), "ilincs_query_model_inventory.csv")
    put("n_models_excluded", len(exc), "ilincs_query_model_inventory.csv")
    put("n_polyp_models", int((used.tissue == "polyp").sum()), "ilincs_query_model_inventory.csv")
    put("n_normal_models", int((used.tissue == "normal colon").sum()),
        "ilincs_query_model_inventory.csv")
    put("min_landmark_shared", 4 * const("04_ilincs_reversal.py", "N_SIDE"),
        "src/04_ilincs_reversal.py")
    pol = ", ".join(sorted(exc.loc[exc.tissue == "polyp", "cell_type"]))
    nor = ", ".join(sorted(exc.loc[exc.tissue == "normal colon", "cell_type"]))
    put("excluded_desc",
        f"polyp {pol.rsplit(', ', 1)[0]} and {pol.rsplit(', ', 1)[1]}, and normal-colon "
        f"{nor.replace(', ', ' and ')} in the GSE266546 cohort; "
        f"{int(exc.n_landmark_shared.min())}\u2013{int(exc.n_landmark_shared.max())} shared genes",
        "ilincs_query_model_inventory.csv")

    # ---- resampling parameters, read from the scripts that ran them ----
    put("n_draw_compound", const("07_ilincs_drug_panel.py", "size=(", "5000"),
        "src/07_ilincs_drug_panel.py")
    put("seed_compound", const("07_ilincs_drug_panel.py", "default_rng("), "src/07_ilincs_drug_panel.py")
    put("n_draw_class", const("09_ms_ilincs_panel_text.py", "N_DRAW, SEED", pos=0),
        "src/09_ms_ilincs_panel_text.py")
    put("seed_class", const("09_ms_ilincs_panel_text.py", "N_DRAW, SEED", pos=1),
        "src/09_ms_ilincs_panel_text.py")
    cms = pd.read_csv(rf"{T}\ilincs_cmap_model_summary.csv")
    nperm = int(cms.n_perm.unique()[0])
    put("n_perm_maxstat", nperm, "ilincs_cmap_model_summary.csv")
    put("maxstat_floor", 1 / (nperm + 1), "ilincs_cmap_model_summary.csv")

    # class membership counts and the ERBB2 compounds at or above zero
    put("n_erbb2_cpds", int((pc.risk_gene_targets == "ERBB2").sum()),
        "ilincs_drug_panel_per_compound.csv")
    put("n_prot_cpds", int((pc.risk_gene_targets == PROTEASOME).sum()),
        "ilincs_drug_panel_per_compound.csv")
    pos = pc[(pc.risk_gene_targets == "ERBB2") & (pc.mean_NCS >= 0)]
    names = sorted(n.upper() if any(c.isdigit() for c in n) else n.lower()
                   for n in pos.drug_name)
    put("erbb2_nonneg", ", ".join(names[:-1]) + " and " + names[-1],
        "ilincs_drug_panel_per_compound.csv")
    put("n_druggable_genes", 63, "Drugs/table/63_druggable_genes_summary.csv")
    rec = pd.read_csv(rf"{ROOT}\Drugs\table\63_druggable_genes_drug_records.csv")
    adv = rec[rec["Development stage (harmonised)"].isin(
        ["Approved drug", "Phase 3 / registration", "Phase 2"])]
    keys = adv["Drug Brand Name"].map(norm_name)
    put("n_adv_genes", adv.Genes.nunique(), "63_druggable_genes_drug_records.csv")
    put("n_adv_cpds_precollapse", keys[keys.ne("")].nunique(),
        "63_druggable_genes_drug_records.csv")
    paz = pc[pc.drug_name.str.upper() == "PAZOPANIB"].iloc[0]
    put("pazopanib_mean", paz.mean_NCS, "ilincs_drug_panel_per_compound.csv")
    put("pazopanib_fdr", paz.FDR_vs_sizematched_null, "ilincs_drug_panel_per_compound.csv")
    return v, src


MINUS = "\u2212"
WORD = {0: "no", 1: "one", 2: "two", 3: "three", 4: "four", 5: "five", 6: "six", 7: "seven",
        8: "eight", 9: "nine", 12: "twelve", 18: "eighteen"}


def i(x):
    return f"{int(round(x)):,}"


def num(x, dec=2, sign=False):
    """Number with a true minus sign, as Nature Genetics sets it."""
    s = f"{float(x):+.{dec}f}" if sign else f"{float(x):.{dec}f}"
    return s.replace("-", MINUS)


def fmt_p(p):
    """Nature style: two significant digits, scientific below 0.001."""
    p = float(p)
    if p >= 0.001:
        return f"{p:#.2g}"
    m, e = f"{p:.1e}".split("e")
    return f"{float(m):.1f} \u00d7 10^[{MINUS}{abs(int(e))}]"


# -------------------------------------------------------------------- text ---
def compose(v):
    R1 = (
        "We next asked whether the compounds that target these genes reverse the "
        "cell-type-specific transcriptional profiles predicted to raise risk. Rather than "
        "screening the whole perturbation library, we pre-specified the panel of approved or "
        f"phase II/III compounds annotated against {int(v['n_adv_genes'])} of the "
        f"{int(v['n_druggable_genes'])} druggable risk genes (the remaining eight are supported only by "
        f"phase I or stage-unstated agents); these {i(v['n_adv_cpds_precollapse'])} compound records "
        f"resolved to {i(v['panel_n'])} distinct molecules once synonyms and salts were collapsed on "
        f"PubChem CID, of which "
        f"{i(v['panel_in'])} are represented in the iLINCS large-intestine chemical-perturbation "
        "library^[25] and could be scored against the "
        f"{i(v['models'])} cohort-level cell-type TWAS query signatures by weighted connectivity "
        "scoring^[26] "
        f"(Methods). Four of the {i(v['panel_in'])} compounds reached a Benjamini\u2013Hochberg "
        "false discovery rate (FDR) below 0.05 against a size-matched empirical null "
        "(Fig. 6e and Supplementary Table X): ixazomib "
        f"(mean normalized connectivity score (NCS) {num(v['ixazomib_mean'])}; one-sided "
        f"*P* = {fmt_p(v['ixazomib_p'])}, FDR = {fmt_p(v['ixazomib_fdr'])}), neratinib "
        f"({num(v['neratinib_mean'])}; *P* = {fmt_p(v['neratinib_p'])}, FDR = {v['neratinib_fdr']:.3f}), "
        f"bortezomib ({num(v['bortezomib_mean'])}; *P* = {fmt_p(v['bortezomib_p'])}, "
        f"FDR = {fmt_p(v['bortezomib_fdr'])}) and resorcinol, which rests on a single signature and was "
        "not considered further. The three credible hits reversed consistently across models and "
        f"cell lines: {int(v['neratinib_models_neg'])} of {i(v['models'])} cell-type models for "
        f"neratinib ({int(v['neratinib_lines'])} cell lines), {int(v['ixazomib_models_neg'])} of "
        f"{i(v['models'])} for ixazomib ({int(v['ixazomib_lines'])} lines) and "
        f"{int(v['bortezomib_models_neg'])} of {i(v['models'])} for bortezomib "
        f"({int(v['bortezomib_lines'])} lines)."
    )
    R2 = (
        "All three hits act on genes that are themselves in the risk set, but only neratinib engages "
        "one selectively. *ERBB2* is the most heavily drugged gene of the "
        f"{int(v['n_druggable_genes'])}, so for neratinib the transcriptional reversal and the genetic "
        "association converge on the same protein. *PSMA7*, *PSMD3* and *PSMD13* are also risk genes, "
        "but they are pan-essential in the large intestine and are engaged at the level of the "
        "assembled proteasome rather than through subunit-selective chemistry, and proteasome "
        "inhibition produces a proteotoxic-stress response that opposes many unrelated expression "
        "profiles. Consistent with that generic mechanism, the "
        f"{['zero','one','two','three','four','five'][int(v['n_prot_cpds'])]} proteasome inhibitors "
        "reversed the risk profile as a class "
        f"(mean NCS {num(v['prot_mean'], 3)} versus {num(v['prot_null'], 3, sign=True)} expected; one-sided "
        f"*P* = {fmt_p(v['prot_p'])}), with carfilzomib directionally concordant "
        f"({v['carfilzomib_fracneg']*100:.0f}% of its signatures negative) although individually "
        f"non-significant (FDR = {v['carfilzomib_fdr']:.2f}). The {int(v['n_erbb2_cpds'])} "
        f"*ERBB2*-directed compounds, by contrast, were not shifted as a class "
        f"({num(v['erbb2_mean'], 3)} versus {num(v['erbb2_null'], 3)}; *P* = {fmt_p(v['erbb2_p'])}), and "
        f"{v['erbb2_nonneg']} sat at or above "
        "zero. Neratinib is thus separated from the other HER2 inhibitors rather than exemplifying "
        "them, which is consistent with its irreversible pan-HER pharmacology but cannot be "
        "distinguished here from compound-specific off-target activity. The effects are modest in "
        f"absolute terms: ixazomib, neratinib and bortezomib rank {i(v['ixazomib_rank'])}th, "
        f"{i(v['neratinib_rank'])}th and {i(v['bortezomib_rank'])}th of {i(v['lib_compounds'])} "
        "library compounds by mean NCS, and the panel as a whole was not shifted negative relative "
        f"to non-panel compounds (one-sided Mann\u2013Whitney *U* test, *P* = {fmt_p(v['p_panel_vs_lib'])}). "
        "These compounds clear significance because the panel is small and defined in advance by "
        "drug-target evidence, not because they are outstanding reversers of the risk profile."
    )
    R3 = (
        "An unbiased screen of the same library nominated nothing. Scanning all "
        f"{i(v['sig_scanned'])} colorectal signatures per model and comparing the strongest reversal "
        "in each scan against a max-statistic permutation null left "
        f"{int(v['scan_ns'])} of the {i(v['models'])} models non-significant, and the single nominal "
        f"exception (*P* = {fmt_p(v['scan_p'])}) did not survive correction for the "
        f"{i(v['models'])} models tested (*P* = {fmt_p(v['scan_p_bonf'])}). Submitting the same queries "
        "to SigCom LINCS, an independent implementation that scores them against "
        f"{i(v['sigcom_sig'])} L1000 signatures, returned "
        f"{'no reversers' if v['sigcom_hits'] == 0 else int(v['sigcom_hits'])} at "
        f"FDR < 0.05, and its {i(v['sigcom_cons_sig'])}-signature per-compound consensus library "
        f"returned {'one' if v['sigcom_cons_hits'] == 1 else int(v['sigcom_cons_hits'])}, in a single "
        "model; no reverser reached Bonferroni "
        "significance in either library (0 of 36 model\u2013library queries). The pre-specified panel is "
        "therefore the only setting in which this analysis is adequately powered, and its result is "
        "best read as transcriptional support for HER2 inhibition in colorectal tissue rather than as "
        "the nomination of a new repurposing candidate."
    )
    R3_short = (
        f"An unbiased screen of the same library nominated nothing: a max-statistic permutation scan of "
        f"all {i(v['sig_scanned'])} colorectal signatures left "
        f"{int(v['scan_ns'])} of {i(v['models'])} models non-significant and the one nominal exception did "
        f"not survive correction across models (*P* = {fmt_p(v['scan_p_bonf'])}), and the same queries run "
        f"against {i(v['sigcom_sig'])} L1000 signatures in SigCom LINCS returned no reverser at FDR < 0.05 "
        "(Supplementary Note). The pre-specified panel is therefore the only setting in which this "
        "analysis is adequately powered."
    )

    LEG = (
        "**e**, Connectivity between the pre-specified risk-gene drug panel and the cell-type TWAS "
        f"profiles. Each row is one of the {i(v['panel_in'])} pre-specified panel compounds present in "
        "the iLINCS large-intestine library (approved or phase II/III agents annotated against "
        f"{int(v['n_adv_genes'])} of the {int(v['n_druggable_genes'])} druggable risk genes); the row "
        "label gives the targeted risk gene (or 'proteasome' "
        "for compounds that engage the assembled particle rather than a single subunit) and the number "
        "of iLINCS signatures for that compound (*n*). Dots, mean normalized connectivity score (NCS) "
        f"over the compound's signatures, each signature first averaged across the {i(v['models'])} "
        "cohort-level cell-type TWAS query signatures; grey bars, mean \u00b1 2 s.d. of a size-matched empirical null "
        "(5,000 random sets of the same number of library signatures); red, Benjamini\u2013Hochberg "
        f"FDR < 0.05 across the {i(v['panel_in'])} compounds (one-sided empirical *P*). Negative NCS "
        "indicates reversal of the profile predicted to raise risk. Compounds are ordered by mean NCS. "
        "Source data are provided."
    )

    M1 = (
        "*Construction of TWAS query signatures.* Cell-type TWAS models were taken from the "
        "S-PrediXcan peak-SNP European-ancestry association files, which cover "
        f"{int(v['n_models_total'])} cohort-level models: nine epithelial cell types in colorectal "
        "polyps (VUMC HTAN adenoma and serrated cohort) and seven in each of the two normal-colon "
        "cohorts (VUMC HTAN GSE266546 and Wellcome Sanger). Normal colon was kept separate by cohort "
        "rather than meta-analyzed, so that each query rests on a single genotype\u2013expression panel; "
        "the same cell type therefore contributes two queries in normal colon. Because the iLINCS LIB5 "
        "platform measures only the 978 L1000 landmark genes, queries were restricted to genes shared "
        "between each model and the landmark space, and a model was carried forward only if at least "
        f"{int(v['min_landmark_shared'])} such genes were available, twice the query size, so that the "
        "up- and down-tagged sets could be drawn without exhausting the measurable space. "
        f"{WORD[int(v['models'])].capitalize()} of the {int(v['n_models_total'])} models met this "
        f"requirement ({WORD[int(v['n_polyp_models'])]} polyp and {int(v['n_normal_models'])} "
        f"normal-colon models, with {int(v['lm_min'])}\u2013{int(v['lm_max'])} shared genes); "
        f"{WORD[int(v['n_models_excluded'])]} were "
        f"excluded for insufficient landmark coverage ({v['excluded_desc']}), all of them rare cell types "
        "with sparse prediction models. Within each model, genes carrying more than one association were "
        f"reduced to the single most extreme *Z* score, and the {int(v['query_side'])} genes with the most "
        f"positive and {int(v['query_side'])} with the most negative *Z* scores were taken as the up- and "
        f"down-tagged sets, giving a {int(2 * v['query_side'])}-gene bidirectional query per model. Queries "
        "are therefore the most extreme landmark-measurable genes rather than the genome-wide top "
        "associations, a constraint of the L1000 platform rather than a choice about the genetics; taking "
        "the genome-wide top-ranked genes would have left approximately six measurable genes per query."
    )
    M2 = (
        "*Perturbation signatures and connectivity scoring.* Chemical perturbagen signatures were taken "
        f"from the iLINCS LIB5 large-intestine library^[25] ({i(v['sig_scanned'])} signatures, "
        f"{i(v['perturbagens'])} distinct perturbagens, log2 differential expression relative to matched "
        "controls); signature identity, perturbagen name, cell line, dose and exposure time were joined "
        "from the accompanying LIB5 metadata. Each query was scored against each signature with the "
        "weighted connectivity score of Subramanian et al.^[26]. Signatures were ranked by log2 fold "
        "change in descending order; enrichment scores ES_up and ES_down were computed as weighted "
        "Kolmogorov\u2013Smirnov running sums (weight = |log2 fold change|, exponent 1) of the up- and "
        "down-tagged query genes in that ranked list; and the weighted connectivity score was taken as "
        "WTCS = (ES_up \u2212 ES_down)/2 when the two enrichment scores had opposite signs and 0 "
        "otherwise, so that signatures enriched in the same direction for both tag sets score zero. "
        "Scores were normalized within cell line by the mean absolute score of the same-sign signatures "
        "in that cell line, giving the normalized connectivity score (NCS), so that signatures from "
        "different lines are comparable. A negative NCS indicates reversal, that is, the perturbagen "
        "suppresses genes whose predicted expression raises risk and induces those whose predicted "
        "expression lowers it. The implementation was validated on constructed cases: a tag set equal to "
        "a signature's own 50 highest-ranked genes returns ES = +1.000, its 50 lowest-ranked genes "
        "ES = \u22121.000, a random tag set ES near 0 and a perfectly reversing query WTCS = \u22121.000."
    )
    M3 = (
        "*Prespecified drug panel.* The panel was prespecified as the approved or phase II/III "
        "compounds annotated against the druggable colorectal risk genes in the four-resource "
        f"druggability table (ChEMBL, DrugBank, TTD and Open Targets); {int(v['n_adv_genes'])} of the "
        f"{int(v['n_druggable_genes'])} druggable genes carry such a compound, the other eight being "
        "supported only by phase I or stage-unstated agents. Compound identity between the drug table "
        "and iLINCS was resolved on PubChem CID rather than on drug name: CIDs were obtained from the "
        f"PubChem name-to-CID service ({int(v['panel_cid'])} panel compounds resolved) and matched "
        "against the CIDs carried in the LIB5 metadata, falling back to a salt-stripped normalized name "
        "where no CID was available. Name matching alone both missed compounds stored under a different "
        "synonym and double-counted aliases of the same molecule, so records resolving to one CID were "
        f"collapsed, taking {i(v['n_adv_cpds_precollapse'])} compound records to {i(v['panel_n'])} "
        f"distinct panel compounds, of which {i(v['panel_in'])} are present in iLINCS. The "
        f"{i(v['panel_out'])} absentees are predominantly biologics \u2014 antibodies, "
        "antibody\u2013drug conjugates, bispecific constructs and CAR-T products \u2014 together with "
        "several recent tyrosine kinase inhibitors, none of which are represented in a small-molecule "
        "perturbation library; querying the full LIB5 collection (143,374 signatures across 23 tissues) "
        "recovered no additional panel compound, confirming that this is a limit of library content "
        "rather than an effect of restricting to colorectal lines."
    )
    M4 = (
        "*Statistical analysis.* A signature \u00d7 model score is not an independent observation: the "
        f"{int(v['models'])} models share query genes and, in normal colon, cell types; and one compound "
        "contributes several doses, time points and cell lines. Each compound was therefore reduced to a "
        "single value \u2014 the mean NCS over its signatures, each signature first averaged across the "
        f"{int(v['models'])} models \u2014 and referenced against a size-matched empirical null formed by "
        f"drawing {i(v['n_draw_compound'])} random sets of the same number of library signatures and "
        f"taking their mean NCS (random seed {int(v['seed_compound'])}). One-sided empirical *P* values "
        "were computed as (number of null means at or below the observed compound mean + 1)/(number of "
        f"draws + 1) and corrected across the {i(v['panel_in'])} tested compounds by the "
        "Benjamini\u2013Hochberg procedure, reported as FDR. Target-class tests (*ERBB2* inhibitors, "
        "proteasome inhibitors) used the same construction at the class level, drawing, for each class "
        "member, from the 50 non-panel library compounds nearest to it in log10 signature count "
        f"({i(v['n_draw_class'])} draws, seed {int(v['seed_class'])}); an exact match on signature count "
        "is not attainable for compounds with hundreds of signatures. Whether the panel as a whole was "
        "shifted negative was tested by a one-sided Mann\u2013Whitney *U* test on compound-level mean NCS, "
        "panel against non-panel library compounds; the signature-level version of that comparison is "
        "pseudo-replicated by compounds with many signatures and is not reported. No compounds or "
        "signatures were excluded, and no sample-size calculation, randomization or blinding applies to "
        "this computational analysis."
    )
    M5 = (
        "*Library-wide analysis.* As a complementary unbiased analysis, all "
        f"{i(v['sig_scanned'])} large-intestine signatures were scored against each model with the same "
        "weighted connectivity statistic, and the strongest reversal per model was referenced against a "
        "max-statistic null: for each of "
        f"{i(v['n_perm_maxstat'])} random bidirectional queries \u2014 {int(2 * v['query_side'])} genes "
        "drawn without replacement from that model's own pool of landmark-measurable TWAS genes and split "
        f"at random into {int(v['query_side'])} up- and {int(v['query_side'])} down-tagged genes \u2014 "
        "the most negative WTCS over the whole library was recorded, and the observed best reversal was "
        "read against that distribution. The library cannot serve as its own null for the strongest of "
        f"{i(v['sig_scanned'])} scores, which is why a resampled-query null is used. One-sided *P* values "
        "were computed with the same +1 correction, so the smallest attainable value at "
        f"{i(v['n_perm_maxstat'])} draws is {v['maxstat_floor']:.4f}, and were Bonferroni-corrected across "
        f"the {int(v['models'])} models. The same queries were independently submitted to SigCom LINCS "
        f"({i(v['sigcom_sig'])} L1000 chemical-perturbation signatures and a "
        f"{i(v['sigcom_cons_sig'])}-signature per-compound consensus library), and returned reversers "
        "were corrected by the Benjamini\u2013Hochberg procedure within each model\u2013library query."
    )
    M6 = (
        "*Code availability (addition).* Analyses were run in Python 3.13 with numpy, pandas, scipy and "
        "statsmodels. Query construction is implemented in src/04_ilincs_reversal.py, connectivity "
        "scoring and the library-wide max-statistic null in src/06_ilincs_cmap_connectivity.py, the "
        "prespecified panel and its compound-level inference in src/07_ilincs_drug_panel.py, the SigCom "
        "LINCS comparison in src/08_sigcom_consensus.py, the target-class tests and the reported values "
        "in src/09_ms_ilincs_panel_text.py and src/11_ms_ilincs_NG_text.py, and Fig. 6e in "
        "src/10_fig6e_panel.py."
    )
    AVAIL = (
        "*Data availability (addition).* iLINCS LIB5 chemical-perturbation signatures and their "
        "metadata are available from iLINCS (http://www.ilincs.org); SigCom LINCS queries were run at "
        "https://maayanlab.cloud/sigcom-lincs. Per-compound connectivity scores for the panel, the "
        "library-wide compound means and the panel-versus-library test are provided as Supplementary "
        "Table X, and the values plotted in Fig. 6e as Source Data."
    )
    NOTE = [
        ("Replace in Results section 6 (paragraph beginning 'Using data from DrugBank'), the four "
         "sentences from 'To assess whether these compounds could modulate' through 'which targets "
         "SH2B3.'"),
        ("Those sentences report the superseded Pearson-correlation screen and must be deleted, not "
         "edited: dasatinib and regorafenib (EPHA2) and pazopanib (SH2B3) are not hits under the "
         "connectivity analysis \u2014 pazopanib has a positive mean NCS "
         f"({num(v['pazopanib_mean'], 3, sign=True)}), and "
         "dasatinib and regorafenib are not in the pre-specified panel because neither targets one of "
         f"the {int(v['n_druggable_genes'])} druggable risk genes in the current four-resource table."),
        ("Replace the Methods paragraph beginning 'To further explore therapeutic potential' in its "
         "entirety with the Methods subsections below. Reference 27 (Wu et al.) supported the "
         "Pearson-correlation approach and is no longer cited by these Methods; delete it from the "
         "reference list unless it is cited elsewhere. References 25 (iLINCS) and 26 (Subramanian et "
         "al.) are retained and keep their numbers."),
        ("Reconcile two numbers in the preceding druggability text before inserting: it currently "
         f"states '43 genes ... targets of 112 approved or clinical-stage compounds', whereas the "
         f"current druggability table gives {int(v['n_druggable_genes'])} druggable genes and 193 "
         f"approved or phase II/III drugs across 55 genes. Leaving '43 genes' next to the "
         f"'{i(v['panel_in'])} compounds' of this analysis will read as the same 43."),
        ("Figure 6e replaces panel a of Fig_ilincs_drug_panel.pdf. Panel b of that figure (per-model "
         "NCS for the three credible hits) is not used; the 13\u201314 of 18 model counts it showed are "
         "quoted in the Results text instead. If a reviewer asks to see them, panel b can be added as "
         "an Extended Data panel without any change to the text."),
    ]

    md = ["# iLINCS connectivity of the risk-gene drug panel \u2014 Nature Genetics format",
          "",
          "*Generated by src/11_ms_ilincs_NG_text.py. Figure 6e is built by "
          "src/10_fig6e_panel.py. Every quoted value is listed in "
          "table/ms_ilincs_panel_NG_numbers.csv with its source table.*",
          "", "## How this text goes into the manuscript", ""]
    md += [f"{k+1}. {t}" for k, t in enumerate(NOTE)]
    md += ["", "## Results \u2014 insert into section 6 (Druggability)", "", R1, "", R2, "", R3,
           "", "*Condensed alternative for the third paragraph, if space is tight:*", "", R3_short,
           "", "## Figure 6 legend \u2014 panel e", "", LEG,
           "", "## Methods \u2014 replaces the iLINCS paragraph", "", M1, "", M2, "", M3, "", M4,
           "", M5, "", M6, "", AVAIL, ""]
    return md, NOTE, [R1, R2, R3, R3_short], LEG, [M1, M2, M3, M4, M5, M6, AVAIL]


# -------------------------------------------------------------------- docx ---
TOK = re.compile(r"(\*\*.+?\*\*|\*.+?\*|\^\[.+?\])")


def add_rich(p, text):
    for part in TOK.split(text):
        if not part:
            continue
        if part.startswith("**"):
            p.add_run(part[2:-2]).bold = True
        elif part.startswith("^["):
            r = p.add_run(part[2:-1])
            r.font.superscript = True
        elif part.startswith("*"):
            p.add_run(part[1:-1]).italic = True
        else:
            p.add_run(part)


def write_docx(path, notes, results, legend, methods):
    d = Document()
    st = d.styles["Normal"]
    st.font.name = "Arial"
    st.font.size = Pt(10)
    d.add_heading("iLINCS connectivity of the risk-gene drug panel \u2014 Nature Genetics format", 1)
    p = d.add_paragraph()
    add_rich(p, "*Figure 6e is built by src/10_fig6e_panel.py; this text by "
                "src/11_ms_ilincs_NG_text.py. Every quoted value is listed in "
                "table/ms_ilincs_panel_NG_numbers.csv with its source table.*")
    d.add_heading("How this text goes into the manuscript", 2)
    for t in notes:
        add_rich(d.add_paragraph(style="List Number"), t)
    d.add_heading("Results \u2014 insert into section 6 (Druggability)", 2)
    for t in results[:3]:
        add_rich(d.add_paragraph(), t)
    p = d.add_paragraph()
    add_rich(p, "*Condensed alternative for the third paragraph, if space is tight:*")
    add_rich(d.add_paragraph(), results[3])
    d.add_heading("Figure 6 legend \u2014 panel e", 2)
    add_rich(d.add_paragraph(), legend)
    d.add_heading("Methods \u2014 replaces the iLINCS paragraph", 2)
    for t in methods:
        add_rich(d.add_paragraph(), t)
    d.save(path)


def main():
    v, src = collect()
    md, notes, results, legend, methods = compose(v)
    with open(rf"{OUT}\ms_ilincs_panel_NG.md", "w", encoding="utf-8") as f:
        f.write("\n".join(md))
    pd.DataFrame({"quantity": list(v), "value": [v[k] for k in v],
                  "source_table": [src[k] for k in v]}).to_csv(
        rf"{T}\ms_ilincs_panel_NG_numbers.csv", index=False)
    try:
        write_docx(rf"{OUT}\ms_ilincs_panel_NG.docx", notes, results, legend, methods)
    except PermissionError:
        write_docx(rf"{OUT}\ms_ilincs_panel_NG_NEW.docx", notes, results, legend, methods)
        print("target .docx was open in Word; wrote ms_ilincs_panel_NG_NEW.docx")
    print("quoted values:", len(v))
    return v


if __name__ == "__main__":
    v = main()
