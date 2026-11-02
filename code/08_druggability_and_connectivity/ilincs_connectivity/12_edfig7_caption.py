#!/usr/bin/env python
"""Extended Data Fig. 7 caption and reader's explanation.

The iLINCS connectivity panel was renumbered: it is no longer Fig. 6e but
Extended Data Fig. 7, cited from the Fig. 6 legend together with Supplementary
Table 22 (Li_etal_scTWAS_10062026.docx / Figures_SupFigures_10062026.docx).
This script emits one final caption in the deck's Extended Data style plus the
explanation behind it, with every quoted value read back from
table/ms_ilincs_panel_NG_numbers.csv.

Outputs: ms_EDFig7_caption.md / .docx
"""
import importlib.util
import os
import pandas as pd
from docx import Document
from docx.shared import Pt

OUT = r"W:\VUMC_projects\scTWAS\iLINCS_reversal"
SRC = rf"{OUT}\src"
spec = importlib.util.spec_from_file_location("s11", rf"{SRC}\11_ms_ilincs_NG_text.py")
s11 = importlib.util.module_from_spec(spec)
spec.loader.exec_module(s11)                      # helpers only; main() is guarded

N = pd.read_csv(rf"{OUT}\table\ms_ilincs_panel_NG_numbers.csv").set_index("quantity").value
f = lambda k: float(N[k])
n = lambda k: int(round(f(k)))

CAPTION = (
    "**Extended Data Fig. 7 | Reversal of the cell-type scTWAS signatures by clinically advanced "
    "compounds targeting risk genes.** Connectivity scoring of the "
    f"{n('panel_in')} prespecified panel compounds represented in the iLINCS LIB5 large-intestine "
    f"perturbation library ({n('sig_scanned'):,} signatures) against the {n('models')} cohort-level "
    "cell-type TWAS signatures; compounds are ordered by mean normalized connectivity score (NCS). "
    "For each compound the point gives the mean NCS across its library signatures, each signature "
    f"first averaged over the {n('models')} cell-type models, and the horizontal bar gives the mean "
    "\u00b1 2 s.d. of a size-matched empirical null formed from "
    f"{n('n_draw_compound'):,} random sets of the same number of library signatures. Red, "
    f"Benjamini\u2013Hochberg FDR < 0.05 across the {n('panel_in')} compounds (one-sided empirical "
    "*P*); grey, not significant. Negative NCS indicates reversal of the risk-associated expression "
    "profile, that is, the compound suppresses genes whose genetically predicted expression raises "
    "risk and induces those whose predicted expression lowers it. Row labels give the targeted risk "
    "gene and, in parentheses, *n*, the number of library signatures contributing to that compound's "
    "mean; 'proteasome' denotes compounds that engage the assembled particle rather than a single "
    "subunit (*PSMA7*, *PSMD3* and *PSMD13*). Compound identity between the drug table and iLINCS was "
    "resolved on PubChem CID rather than on drug name (Methods). Resorcinol, the leftmost compound, "
    "rests on a single signature and was not interpreted further. Per-compound statistics are given "
    "in Supplementary Table 22; source data are provided."
)

EXPLAIN = [
    ("What the figure answers",
     "Figure 6 shows which susceptibility genes are druggable. This panel asks the next question: do "
     "the compounds that hit those genes actually push cell-type expression in the direction opposite "
     "to the one predicted to raise risk. Each row is one compound and the x axis is the direction "
     "and strength of that effect, so a reader scans down the left-hand side for candidates and finds "
     f"that only four of {n('panel_in')} separate from chance ({n('n_fdr')} compounds)."),
    ("How to read a row",
     "The point is the compound's summary score and the grey bar is what that compound would have "
     "scored had its signatures been drawn at random from the library. A point to the left of its own "
     "bar is reversing more than chance allows; a point inside the bar is indistinguishable from "
     "noise however far left it sits. Bar width is set by how many signatures the compound "
     f"contributes, which is why bortezomib (n = {n('bortezomib_nsig')}) has a narrow bar and "
     f"sapitinib or mitotane (n = 1) have bars spanning most of the axis. This is the whole point of "
     "the size-matched null: a single-signature compound has to be far more extreme to count."),
    ("The three credible results",
     f"Ixazomib (mean NCS {s11.num(f('ixazomib_mean'))}, FDR = {f('ixazomib_fdr'):.4f}), neratinib "
     f"({s11.num(f('neratinib_mean'))}, FDR = {f('neratinib_fdr'):.3f}) and bortezomib "
     f"({s11.num(f('bortezomib_mean'))}, FDR = {f('bortezomib_fdr'):.4f}) clear the null, and each "
     f"reverses in {n('ixazomib_models_neg')}\u2013{n('neratinib_models_neg')} of the "
     f"{n('models')} cell-type models rather than resting on one. The fourth significant point, "
     "resorcinol, is the most extreme in the figure but rests on a single HeLa signature, which is "
     "why the caption retires it explicitly \u2014 otherwise a reader takes the leftmost red dot as "
     "the headline result."),
    ("What the figure deliberately does not claim",
     "The scores are modest in absolute terms: among the "
     f"{n('lib_compounds'):,} compounds in the library these three rank {n('ixazomib_rank')}th, "
     f"{n('neratinib_rank')}th and {n('bortezomib_rank'):,}th by mean NCS, and the panel as a whole is "
     f"not shifted negative relative to non-panel compounds (one-sided Mann\u2013Whitney *U* test, "
     f"*P* = {f('p_panel_vs_lib'):.3f}). They clear significance because the panel is small and was "
     "defined in advance by drug-target evidence, not because they are outstanding reversers. The "
     "figure supports a mechanism nomination, not a repurposing claim."),
    ("The one comparison a reviewer will make",
     f"Of the {n('n_erbb2_cpds')} *ERBB2*-directed compounds in the panel, only neratinib is "
     "significant, and the class as a whole is not shifted "
     f"(*P* = {f('erbb2_p'):.2f}); the proteasome inhibitors, by contrast, reverse as a class "
     f"(*P* = {f('prot_p'):.4f}), which is what a generic proteotoxic-stress response predicts. Both "
     "patterns are visible in the figure \u2014 the *ERBB2* rows are scattered across zero while the "
     "three proteasome rows all sit left of it \u2014 so the Results text should state the "
     "interpretation rather than leave the reader to infer it."),
]

NOTES = [
    "The artwork currently embedded as Extended Data Fig. 7 (Figures/FigureS_Drugs_S1.png) is the "
    "earlier export. Replace it with figure/Fig6e_ilincs_drug_panel.pdf, which drops the in-figure "
    "title that now duplicates the caption, italicizes the gene symbols, restores the code-name "
    "capitalization (AEE-788, BMS-599626, CI-1033, CP-724714, AT-7519, PX-12, \u03b1-linolenic acid), "
    "and sets 's.d.' rather than 'SD'. Its file name still says Fig6e and should be renamed at "
    "assembly; the artwork itself carries no panel letter.",
    "The manuscript and the figure deck currently carry different Extended Data Fig. 7 captions \u2014 "
    "a three-sentence stub in Li_etal_scTWAS_10062026.docx and a fuller version in "
    "Figures_SupFigures_10062026.docx. Replace both with the caption above.",
    "Spelling: the figure axis, this caption and the Results all use 'normalized'; the deck caption "
    "currently reads 'normalised'.",
    "The Supplementary Table citation is empty in both current captions. Per the table index it is "
    "Supplementary Table 22 (iLINCS connectivity results), which is also what the Fig. 6 legend "
    "cites.",
    "ms_ilincs_panel_NG.docx still refers to this panel as Fig. 6e throughout; its Results and "
    "Methods text is otherwise current and only the callouts need renumbering to Extended Data "
    "Fig. 7.",
]


def build():
    md = ["# Extended Data Fig. 7 \u2014 caption and explanation", "",
          "*Generated by src/12_edfig7_caption.py; values read from "
          "table/ms_ilincs_panel_NG_numbers.csv.*", "", "## Caption", "", CAPTION, "",
          "## Explanation", ""]
    for h, t in EXPLAIN:
        md += [f"**{h}.** {t}", ""]
    md += ["## Notes for assembly", ""]
    md += [f"{k+1}. {t}" for k, t in enumerate(NOTES)]
    with open(rf"{OUT}\ms_EDFig7_caption.md", "w", encoding="utf-8") as fh:
        fh.write("\n".join(md) + "\n")

    d = Document()
    d.styles["Normal"].font.name = "Arial"
    d.styles["Normal"].font.size = Pt(10)
    d.add_heading("Extended Data Fig. 7 \u2014 caption and explanation", 1)
    d.add_heading("Caption", 2)
    s11.add_rich(d.add_paragraph(), CAPTION)
    d.add_heading("Explanation", 2)
    for h, t in EXPLAIN:
        s11.add_rich(d.add_paragraph(), f"**{h}.** {t}")
    d.add_heading("Notes for assembly", 2)
    for t in NOTES:
        s11.add_rich(d.add_paragraph(style="List Number"), t)
    try:
        d.save(rf"{OUT}\ms_EDFig7_caption.docx")
    except PermissionError:
        d.save(rf"{OUT}\ms_EDFig7_caption_NEW.docx")
    return md


if __name__ == "__main__":
    build()
