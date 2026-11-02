#!/usr/bin/env python
"""Revised iLINCS connectivity Methods for the manuscript, with Q19/Q20 resolved.

Rewrites the "iLINCS connectivity analysis" Methods section of
Li_etal_scTWAS_10062026.docx. Two editorial queries are carried in the draft:

  Q19  Results, Methods and the figure legend say 18 cell-type models while the
       study describes 16 GReX models. Reconcile.
  Q20  A second, contradictory description (Pearson correlation of TWAS Z against
       perturbation Z, Wu et al.) was deleted. Confirm weighted connectivity is
       the analysis to retain.

Both are answered in the notes, the body text no longer carries them, and the
reconciliation is written into the first paragraph rather than left implicit.
Reference numerals follow the current manuscript (59 iLINCS, 60 Subramanian).

Outputs: ms_iLINCS_Methods_revised.md / .docx
"""
import importlib.util
import pandas as pd
from docx import Document
from docx.shared import Pt

OUT = r"W:\VUMC_projects\scTWAS\iLINCS_reversal"
SRC = rf"{OUT}\src"
spec = importlib.util.spec_from_file_location("s11", rf"{SRC}\11_ms_ilincs_NG_text.py")
s11 = importlib.util.module_from_spec(spec)
spec.loader.exec_module(s11)

N = pd.read_csv(rf"{OUT}\table\ms_ilincs_panel_NG_numbers.csv").set_index("quantity").value
inv = pd.read_csv(rf"{OUT}\table\ilincs_query_model_inventory.csv")
f = lambda k: float(N[k])
n = lambda k: int(round(f(k)))
c = lambda x: f"{x:,}"

n_polyp_cand = int((inv.tissue == "polyp").sum())
n_norm_cand = int((inv.tissue == "normal colon").sum())
ws = inv[(inv.dataset.str.contains("WellcomeSanger")) & inv.used_as_query]
vu = inv[(inv.dataset.str.contains("GSE266546")) & inv.used_as_query]
pol = inv[(inv.tissue == "polyp") & inv.used_as_query]
lst = lambda d: ", ".join(sorted(d.cell_type)[:-1]) + " and " + sorted(d.cell_type)[-1]

P1 = (
    "*Query signatures from the cell-type association models.* Connectivity queries were built from "
    "the cohort-level European-ancestry S-PrediXcan results rather than from the cohort-combined "
    "association statistics, so that each query rests on a single genotype\u2013expression panel; the "
    "two normal-colon cohorts therefore contribute separate queries for the cell types they share. "
    f"This gives {n('n_models_total')} candidate models \u2014 {s11.WORD[n_polyp_cand]} polyp "
    f"populations and "
    f"{n_norm_cand} normal-colon models, seven in each cohort. Because the iLINCS LIB5 library "
    f"reports only the 978 L1000 landmark genes, of which {n('lm_min')}\u2013{n('lm_max')} carry an "
    "association *Z* score in a given model, each query was drawn from the genes shared between that "
    "model and the landmark space rather than genome-wide, and a model was carried forward only if at "
    f"least {n('min_landmark_shared')} shared genes were available, twice the query size, so that the "
    f"up- and down-tagged sets could be drawn without exhausting the measurable space. "
    f"{s11.WORD[n('models')].capitalize()} of the {n('n_models_total')} models met this requirement: "
    f"all seven Wellcome Sanger IBDverse normal-colon models, {s11.WORD[len(vu)]} VUMC normal-colon "
    f"models ({lst(vu)}) and {s11.WORD[len(pol)]} polyp models ({lst(pol)}). The "
    f"{s11.WORD[n('n_models_excluded')]} excluded models \u2014 crypt-top, enteroendocrine and tuft cells in "
    "polyps and enteroendocrine and tuft cells in the VUMC normal-colon cohort \u2014 share only "
    f"{int(inv.loc[~inv.used_as_query, 'n_landmark_shared'].min())}\u2013"
    f"{int(inv.loc[~inv.used_as_query, 'n_landmark_shared'].max())} genes with the landmark space, "
    "reflecting the smaller number of genes with predictable expression in these rarer populations "
    "rather than any property of the tested compounds. Within each retained model, genes carrying more "
    f"than one association were reduced to the single most extreme *Z* score, and the {n('query_side')} "
    f"most positive and {n('query_side')} most negative *Z* scores among shared genes were taken as the "
    f"up- and down-tagged sets, giving a {2 * n('query_side')}-gene bidirectional query per model. "
    f"Taking the genome-wide top {2 * n('query_side')} would have left approximately six measurable "
    "genes per query; the query genes are therefore not the genome-wide strongest associations, a "
    "restriction imposed by the L1000 platform rather than by a choice about the genetics."
)
P2 = (
    "*Perturbation signatures and connectivity scoring.* Chemical perturbagen signatures were taken "
    f"from the iLINCS^[59] LIB5 large-intestine library ({c(n('sig_scanned'))} signatures, "
    f"{c(n('perturbagens'))} distinct perturbagens, log2 differential expression relative to matched "
    "controls), with signature identity, perturbagen name, cell line, dose and exposure time joined "
    "from the accompanying metadata. Each query was scored against each signature with the weighted "
    "connectivity score of Subramanian et al.^[60]: signatures were ranked by log2 fold change in "
    "descending order; enrichment scores for the up- and down-tagged query genes were computed as "
    "weighted Kolmogorov\u2013Smirnov running sums (weight |log2 fold change|, exponent 1); and the "
    "weighted connectivity score was taken as half the difference of the two enrichment scores when "
    "they had opposite signs and zero otherwise, so that signatures enriched in the same direction for "
    "both tag sets, which carry no directional information, score zero. Scores were normalized within "
    "cell line by the mean absolute score of the same-sign signatures in that line, giving the "
    "normalized connectivity score (NCS), so that signatures from different lines are comparable. A "
    "negative NCS indicates reversal, that is, the perturbagen suppresses genes whose predicted "
    "expression raises risk and induces those whose predicted expression lowers it. The implementation "
    "was validated on constructed cases: a tag set equal to a signature's own 50 highest-ranked genes "
    "returns an enrichment score of +1.000, its 50 lowest-ranked genes \u22121.000, a random tag set "
    "approximately 0, and a perfectly reversing query a weighted connectivity score of \u22121.000."
)
P3 = (
    "*Compound panel.* The compound panel was prespecified as the approved or phase II/III compounds "
    "annotated against the druggable risk genes in the four-resource table; "
    f"{n('n_adv_genes')} of the {n('n_druggable_genes')} druggable genes carry such a compound, the "
    "remaining eight being supported only by phase I or stage-unstated agents. Compound identity "
    "between the drug table and iLINCS was resolved on PubChem compound identifier (CID) rather than "
    "on name: CIDs were obtained from the PubChem name-to-CID service "
    f"({n('panel_cid')} panel compounds resolved) and matched against the CIDs in the LIB5 metadata, "
    "falling back to a salt-stripped normalized name where no CID was available. Name matching alone "
    "both missed compounds stored under alternative synonyms and double-counted aliases of the same "
    "molecule, so records resolving to one CID were collapsed, taking "
    f"{c(n('n_adv_cpds_precollapse'))} compound records to {c(n('panel_n'))} distinct panel compounds, "
    f"of which {n('panel_in')} are present in iLINCS. The {n('panel_out')} absentees are predominantly "
    "biologics, including antibodies, antibody\u2013drug conjugates, bispecific constructs and CAR-T "
    "products, together with several recent tyrosine-kinase inhibitors, none of which are represented "
    "in a small-molecule perturbation library; querying the full LIB5 collection (143,374 signatures "
    "across 23 tissues) recovered no additional panel compound, confirming a library-content limit "
    "rather than an effect of restricting to colorectal lines."
)
P4 = (
    "*Statistical analysis.* A signature-by-model score is not an independent observation, because the "
    "models share query genes and, in normal colon, cell types, and because one compound contributes "
    "several doses, time points and cell lines. Each compound was therefore reduced to a single value, "
    "the mean NCS across its signatures with each signature first averaged across the "
    f"{n('models')} models, and referenced against a size-matched empirical null formed by drawing "
    f"{c(n('n_draw_compound'))} random sets of the same number of library signatures and taking their "
    f"mean NCS (random seed {n('seed_compound')}). One-sided empirical *P* values were computed as "
    "(number of null means at or below the observed compound mean + 1)/(number of draws + 1) and "
    f"corrected across the {n('panel_in')} tested compounds by the Benjamini\u2013Hochberg procedure. "
    "Target-class tests for *ERBB2* inhibitors and proteasome inhibitors used the same construction at "
    "class level, drawing for each class member from the 50 non-panel library compounds nearest to it "
    f"in log10 signature count ({c(n('n_draw_class'))} draws, seed {n('seed_class')}); an exact match "
    "on signature count is not attainable for compounds with hundreds of signatures. Whether the panel "
    "as a whole was shifted negative was tested by a one-sided Mann\u2013Whitney *U* test on "
    "compound-level mean NCS, panel against non-panel library compounds; the signature-level version "
    "of that comparison is pseudo-replicated by compounds with many signatures and is not reported. No "
    "compounds or signatures were excluded, and no sample-size calculation, randomization or blinding "
    "applies to this computational analysis."
)

OPTIONAL = (
    "*Library-wide analysis (optional addition).* As a complementary unbiased analysis, all "
    f"{c(n('sig_scanned'))} large-intestine signatures were scored against each model with the same "
    "weighted connectivity statistic, and the strongest reversal per model was referenced against a "
    f"max-statistic null: for each of {n('n_perm_maxstat')} random bidirectional queries, "
    f"{2 * n('query_side')} genes drawn without replacement from that model's own pool of "
    f"landmark-measurable genes and split at random into {n('query_side')} up- and {n('query_side')} "
    "down-tagged genes, the most negative weighted connectivity score over the whole library was "
    "recorded. One-sided *P* values used the same +1 correction, so the smallest attainable value at "
    f"{n('n_perm_maxstat')} draws is {f('maxstat_floor'):.4f}, and were Bonferroni-corrected across "
    f"the {n('models')} models; {n('scan_ns')} of the {n('models')} models were non-significant and "
    f"the single nominal exception did not survive correction (*P* = {f('scan_p_bonf'):.3f})."
)

ANSWERS = [
    ("Q19, 16 versus 18",
     "The two numbers count different things and both are correct once qualified. Sixteen is the "
     "number of cell-type models for which associations are reported: nine polyp populations and seven "
     "normal-colon populations, the latter after combining the two normal-colon cohorts. Eighteen is "
     "the number of connectivity query signatures, which are cohort-level and additionally subject to "
     f"the landmark-overlap requirement. The arithmetic is explicit: on the polyp side, "
     f"{s11.WORD[n_polyp_cand]} reported models of which {s11.WORD[len(pol)]} have at least "
     f"{n('min_landmark_shared')} landmark-measurable genes; on the normal-colon side, seven reported "
     f"cell types each contributing up to two cohort-level models, {n_norm_cand} candidates of which "
     f"{n('n_normal_models')} pass \u2014 {n('n_models_total')} candidates and {n('models')} queries in "
     "total. The revised first paragraph states this, so the reader never meets a bare 'models' count. "
     "Elsewhere the fix is to qualify the noun: use '18 cohort-level cell-type signatures' in the "
     "Results sentence and in the figure legend, and keep '16 cell-type models' for the association "
     "results. Note that the legend referred to in the query, Fig. 6b, is now Extended Data Fig. 7; "
     "its caption already carries the qualifier."),
    ("Q20, which analysis to retain",
     "Confirmed: the weighted-connectivity analysis is the one to retain, and the deleted Pearson "
     "description should stay deleted together with its citation. Every number reported in the "
     "Results and plotted in Extended Data Fig. 7 \u2014 the NCS values, the size-matched nulls, the "
     "four FDR-significant compounds \u2014 comes from the connectivity scoring, not from the "
     "correlation analysis. The two statistics are not interchangeable at the point where nominations "
     "are made: across the 18 models they correlate at r = 0.67\u20130.78, but only 1 to 6 of the 20 "
     "most-reversing signatures overlap between them, so quoting a hit from one while describing the "
     "method of the other would be a real inconsistency rather than a wording slip. The correlation "
     "scan survives only as an internal sensitivity check and is not part of the manuscript."),
]

NOTES = [
    "The library-wide unbiased screen is no longer described anywhere in the manuscript. That is "
    "self-consistent, since the Results claim only the prespecified panel, but it leaves no written "
    "answer to the predictable reviewer question of whether the library was examined without "
    "preselection. The optional paragraph above can be inserted if that answer is wanted.",
    "Results sentence to adjust for Q19: 'Forty-three compounds were represented in the iLINCS "
    "large-intestine perturbation dataset and were evaluated using weighted connectivity scores' "
    "should name the 18 cohort-level cell-type signatures explicitly if the number appears there.",
    "Two typographic defects in the current Results paragraph on this analysis, unrelated to the "
    "Methods: a doubled space in 'whether  targeting prioritized genes', and the fragment 'Target "
    "genes;' before 'Because resorcinol was represented by only a single perturbational signature'.",
    "Reference numerals are kept as in the current manuscript: 59 for iLINCS and 60 for Subramanian "
    "et al. The Wu et al. citation that supported the deleted Pearson description should be removed "
    "from the reference list unless it is cited elsewhere.",
]


def build():
    md = ["# iLINCS connectivity analysis \u2014 revised Methods", "",
          "*Generated by src/13_ms_ilincs_methods_revised.py; values read from "
          "table/ms_ilincs_panel_NG_numbers.csv and table/ilincs_query_model_inventory.csv.*", "",
          "## Methods (drop-in replacement for the iLINCS connectivity section)", "",
          P1, "", P2, "", P3, "", P4, "",
          "## Optional addition", "", OPTIONAL, "",
          "## Answers to the embedded queries", ""]
    for h, t in ANSWERS:
        md += [f"**{h}.** {t}", ""]
    md += ["## Consequential edits elsewhere", ""]
    md += [f"{k + 1}. {t}" for k, t in enumerate(NOTES)]
    with open(rf"{OUT}\ms_iLINCS_Methods_revised.md", "w", encoding="utf-8") as fh:
        fh.write("\n".join(md) + "\n")

    d = Document()
    d.styles["Normal"].font.name = "Arial"
    d.styles["Normal"].font.size = Pt(10)
    d.add_heading("iLINCS connectivity analysis \u2014 revised Methods", 1)
    d.add_heading("Methods (drop-in replacement for the iLINCS connectivity section)", 2)
    for t in (P1, P2, P3, P4):
        s11.add_rich(d.add_paragraph(), t)
    d.add_heading("Optional addition", 2)
    s11.add_rich(d.add_paragraph(), OPTIONAL)
    d.add_heading("Answers to the embedded queries", 2)
    for h, t in ANSWERS:
        s11.add_rich(d.add_paragraph(), f"**{h}.** {t}")
    d.add_heading("Consequential edits elsewhere", 2)
    for t in NOTES:
        s11.add_rich(d.add_paragraph(style="List Number"), t)
    try:
        d.save(rf"{OUT}\ms_iLINCS_Methods_revised.docx")
    except PermissionError:
        d.save(rf"{OUT}\ms_iLINCS_Methods_revised_NEW.docx")
    return md


if __name__ == "__main__":
    build()
