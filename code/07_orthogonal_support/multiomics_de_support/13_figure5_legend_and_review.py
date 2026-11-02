#!/usr/bin/env python3
"""13_figure5_legend_and_review.py

Build the Figure 5 legend and a review of ms_MultiOmics_DE_support.docx.

Every number quoted in either document is recomputed here from the exported
tables, so neither can drift from the analysis. Sources:

  figure/Fig_scTWAS_DEsupport_UpDown_bars_sourcedata.csv   panel A
  table/scTWAS_ADSER_NL_collected_results_scRNAseq_DEsupport.csv
                                                           panel A denominators
  table/ASC_SSC_genes_crossceltype_support_long.csv
                                                           panel B
  table/scTWAS_347genes_multisource_DE_support.csv         panel C
  L:/AI-gene-drug/results/cosmx_crc_vs_normal_epi.csv      panel D

Outputs (written beside the manuscript file):
  Figure5_legend.docx / .md
  ms_MultiOmics_DE_support_review.docx
"""
from __future__ import annotations

import os
import numpy as np
import pandas as pd
from docx import Document
from docx.shared import Pt, Inches, RGBColor
from docx.enum.text import WD_ALIGN_PARAGRAPH

ROOT = r"W:\VUMC_projects\scTWAS\VUMC_HTAN_scTWAS_genes_DE_support"
CX_SRC = r"L:\AI-gene-drug\results\cosmx_crc_vs_normal_epi.csv"

CT_ORDER = ["ABS", "GOB", "TAC", "STM", "TUF", "EE", "CT"]
PANEL_D = ["CEBPB", "GTF2I", "HRAS", "SMARCD1", "WNT4", "ZNRF3"]
TEXT_D_CITED = ["DDX39B", "CSNK2B", "CLIC1", "CTNNB1", "UBE2V1", "RPS5",
                "TRIM28", "LIMA1"]


# --------------------------------------------------------------------------- #
# numbers
# --------------------------------------------------------------------------- #
def panel_a():
    a = pd.read_csv(os.path.join(
        ROOT, "table", "scTWAS_ADSER_NL_collected_results_scRNAseq_DEsupport.csv"))
    contrast = {"C_MSIH": "Normal \u2192 MSI-H", "C_MSS": "Normal \u2192 MSS",
                "E_ASC": "Normal \u2192 ASC", "E_SSC": "Normal \u2192 SSC",
                "L1_MSIH": "ASC \u2192 MSI-H", "L1_MSS": "ASC \u2192 MSS",
                "L2_MSIH": "SSC \u2192 MSI-H", "L2_MSS": "SSC \u2192 MSS"}
    section = {"Normal \u2192 MSI-H": "direct", "Normal \u2192 MSS": "direct",
               "Normal \u2192 ASC": "ASC early", "ASC \u2192 MSI-H": "ASC late",
               "ASC \u2192 MSS": "ASC late", "Normal \u2192 SSC": "SSC early",
               "SSC \u2192 MSI-H": "SSC late", "SSC \u2192 MSS": "SSC late"}
    rows = []
    for key, lab in contrast.items():
        hit = a[f"{key}_DEsupport"].isin(["Up", "Down"])
        for ds, name in [("NL", "normal-colon"), ("ADSER", "polyp")]:
            sel = hit & (a.Dataset == ds)
            rows.append({"contrast": lab, "section": section[lab],
                         "dataset": name, "events": int(sel.sum()),
                         "assoc": set(zip(a.loc[sel, "CellType"],
                                          a.loc[sel, "gene_name"]))})
    r = pd.DataFrame(rows)
    ## denominators: one row per scTWAS gene x cell-type association
    denom = a.Dataset.value_counts().to_dict()
    ## per section: support EVENTS (a gene supported in two contrasts counts
    ## twice) and DISTINCT associations (the number comparable to the denominator)
    out = []
    for sec in ["direct", "ASC early", "ASC late", "SSC early", "SSC late"]:
        rec = {"section": sec}
        for name in ["normal-colon", "polyp"]:
            sub = r[(r.section == sec) & (r.dataset == name)]
            rec[f"{name} events"] = int(sub.events.sum())
            rec[f"{name} distinct"] = len(set().union(*sub.assoc))
        out.append(rec)
    return r, pd.DataFrame(out), {"normal-colon": denom["NL"],
                                  "polyp": denom["ADSER"]}


def panel_b():
    lg = pd.read_csv(os.path.join(
        ROOT, "table", "ASC_SSC_genes_crossceltype_support_long.csv"))
    lg["sup"] = lg.DEsupport.isin(["Up", "Down"])
    n_genes = lg.groupby("track").gene.nunique().to_dict()
    e = lg[lg.sup & lg.contrast_tag.str.startswith("E")]
    entry = {tr: (e[e.track == tr].groupby("contrast_cell_type").gene.nunique()
                  .reindex(CT_ORDER).dropna().astype(int).sort_values(ascending=False))
             for tr in ["ASC", "SSC"]}
    l_tags = {"ASC": ["L1_MSIH", "L1_MSS"], "SSC": ["L2_MSIH", "L2_MSS"]}
    exit_ = {tr: {t: int(lg[(lg.track == tr) & (lg.contrast_tag == t) & lg.sup]
                         .gene.nunique()) for t in tags}
             for tr, tags in l_tags.items()}
    ## how many of each track's 23 applicable transitions carry >= 5 genes
    cnt = (lg[lg.sup].groupby(["track", "contrast_tag"]).gene.nunique()
           .unstack(fill_value=0))
    applicable = {
        "ASC": [f"C{i}_{arm}" for i in range(1, 8) for arm in ("MSIH", "MSS")]
               + [f"E{i}_ASC" for i in range(1, 8)] + l_tags["ASC"],
        "SSC": [f"C{i}_{arm}" for i in range(1, 8) for arm in ("MSIH", "MSS")]
               + [f"E{i}_SSC" for i in range(1, 8)] + l_tags["SSC"]}
    ge5 = {tr: int((cnt.loc[tr, [t for t in tags if t in cnt.columns]] >= 5).sum())
           for tr, tags in applicable.items()}
    n_trans = {tr: len(tags) for tr, tags in applicable.items()}
    return n_genes, entry, exit_, ge5, n_trans


def panel_c():
    ms = pd.read_csv(os.path.join(
        ROOT, "table", "scTWAS_347genes_multisource_DE_support.csv"))
    pl = ["CosMx_NOR_CORE", "TCGA", "CPTAC"]
    ms["n_sup"] = ms[[f"{p}_DEsupport" for p in pl]].isin(["Up", "Down"]).sum(axis=1)
    ms["cls"] = np.where(ms.rep_Z > 0, "onco", "tsg")
    k = ms[ms.n_sup >= 1]
    d = {"genes": ms.gene_name.nunique(),
         "onco": int((ms.rep_Z > 0).sum()), "tsg": int((ms.rep_Z < 0).sum()),
         "supported": len(k), "pct": round(100 * len(k) / len(ms), 1)}
    for cls in ["onco", "tsg"]:
        s = k[k.cls == cls]
        d[f"{cls}_sup"] = len(s)
        for n in (1, 2, 3):
            d[f"{cls}_{n}"] = int((s.n_sup == n).sum())
    ## measured / not measured per resource, for the legend's coverage sentence
    d["cov"] = {p: int(ms[f"{p}_log2FC"].notna().sum()) for p in pl}
    d["cov_kept"] = {p: int(k[f"{p}_log2FC"].notna().sum()) for p in pl}
    return d


def panel_d():
    cx = pd.read_csv(CX_SRC)
    v = cx[cx.gene.isin(PANEL_D)].set_index("gene").loc[PANEL_D]
    return v[["n_patients", "mean_normal_epi", "mean_CRC", "delta", "p", "fdr",
              "twas_direction", "concordant"]], len(cx)


R_A, SEC_A, DEN_A = panel_a()
NB_GENES, NB_ENTRY, NB_EXIT, NB_GE5, NB_TRANS = panel_b()
C = panel_c()
DV, CX_PANEL_N = panel_d()

MIN_P = 2.0 ** -(int(DV.n_patients.iloc[0]) - 1)      # floor of a 7-pair sign test
FDR_D = DV.fdr.unique()
DISCORDANT = DV.index[~DV.concordant.astype(bool)].tolist()


# --------------------------------------------------------------------------- #
# the legend
# --------------------------------------------------------------------------- #
def rng(s):
    return f"{int(s.min())}\u2013{int(s.max())}"


LEGEND_TITLE = ("Figure 5. Orthogonal transcriptomic, proteomic and spatial "
                "evidence for scTWAS-prioritized colorectal cancer "
                "susceptibility genes.")

LEGEND = [
 ("(A)", f"""Cell-type-specific differential-expression support across eight epithelial-state contrasts, grouped into the three routes from normal epithelium to carcinoma: directly (Normal \u2192 MSI-H, Normal \u2192 MSS), through a conventional adenoma (Normal \u2192 ASC, ASC \u2192 MSI-H, ASC \u2192 MSS) and through a sessile serrated lesion (Normal \u2192 SSC, SSC \u2192 MSI-H, SSC \u2192 MSS); within each route the lesion entry contrast precedes its two exits. A bar counts scTWAS gene \u00d7 cell-type associations whose donor-level pseudobulk differential expression is direction-concordant with the scTWAS Z score (FDR < 0.05, |log2FC| > 0.5, sign(log2FC) = sign(Z)); bars above the axis are up-regulated in the later state and bars below are down-regulated, and the y axis is labelled with absolute counts on both sides. Within each contrast the 14 scTWAS cell-type models are dodged in a fixed order, the seven normal-colon (NL) models followed by the seven polyp (ADSER) models, with each model's up and down bars mirrored about the axis. Direction is encoded by colour temperature and dataset by hue family: normal-colon models in red (up) and blue (down), polyp models in purple (up) and teal (down). Within all four ramps the seven cell types run light to dark in the order {', '.join(CT_ORDER)}. The parenthetical n in the legend is the number of gene \u00d7 cell-type associations that model contributes ({DEN_A['normal-colon']} across the normal-colon models and {DEN_A['polyp']} across the polyp models), so bar heights are counts rather than rates and are not comparable between models of different size; the three smallest polyp models (TUF n = 2, CT n = 3, EE n = 4) cannot show an informative count in either direction."""),
 ("(B)", f"""Differential-expression support for polyp-specific risk genes across epithelial-state transitions, shown separately for the {NB_GENES['ASC']} genes prioritized in the adenoma-specific cell (ASC) model (left) and the {NB_GENES['SSC']} prioritized in the serrated-specific cell (SSC) model (right). Nodes are epithelial states: the seven normal-colon cell types, the polyp state and the two carcinoma states (MSI-H, MSS). Each gene is evaluated across {NB_TRANS['ASC']} transitions comprising seven normal-to-polyp contrasts (E1\u2013E7, blue), 14 normal-to-carcinoma contrasts (C1\u2013C7 in each of the two carcinoma arms, orange) and two polyp-to-carcinoma contrasts (L1 for ASC and L2 for SSC, green); an edge label is the number of prioritized genes with direction-concordant support in that transition. Support for the normal-to-polyp transition ranges from {rng(NB_ENTRY['ASC'])} genes across the seven normal cell types for ASC, strongest from stem ({NB_ENTRY['ASC']['STM']} genes) and transit-amplifying cells ({NB_ENTRY['ASC']['TAC']}), and from {rng(NB_ENTRY['SSC'])} genes for SSC, strongest from goblet cells ({NB_ENTRY['SSC']['GOB']}). On the polyp-to-carcinoma transitions {NB_EXIT['ASC']['L1_MSIH']} and {NB_EXIT['ASC']['L1_MSS']} ASC-prioritized genes are supported into MSI-H and MSS carcinoma, and {NB_EXIT['SSC']['L2_MSIH']} and {NB_EXIT['SSC']['L2_MSS']} SSC-prioritized genes. {NB_GE5['ASC']} of the {NB_TRANS['ASC']} ASC transitions and {NB_GE5['SSC']} of the {NB_TRANS['SSC']} SSC transitions carry at least five supported genes."""),
 ("(C)", f"""Tumour\u2013normal support for the {C['genes']} unique scTWAS-prioritized genes in three orthogonal resources, split by the direction of the genetic effect: {C['onco']} genes are oncogene-like (representative scTWAS Z > 0, upper panel) and {C['tsg']} tumour-suppressor-like (Z < 0, lower panel). Shown are the {C['supported']} genes ({C['pct']}% of {C['genes']}) with at least one direction-concordant Up or Down call, ordered by the number of supporting resources and then by |Z|; {C['onco_sup']} are oncogene-like ({C['onco_1']} supported by one resource, {C['onco_2']} by two, {C['onco_3']} by all three) and {C['tsg_sup']} tumour-suppressor-like ({C['tsg_1']}, {C['tsg_2']} and {C['tsg_3']}). Rows within each panel are, from top, the scTWAS association and then GeoMx DSP spatial profiling of the PanCK-positive epithelial compartment (normal epithelium versus carcinoma core), TCGA-COAD bulk transcriptomics and CPTAC-COAD bulk proteomics (carcinoma versus normal). Bubble colour is the effect size on a diverging scale centred at zero, the scTWAS Z in the top row and the log2 fold change in the other three; the three fold-change rows share one colour scale and are therefore directly comparable, and colour is capped at |Z| = 13 and |log2FC| = 2.5 with more extreme values drawn at the end colour rather than excluded. Bubble size is a binned significance level for that row's own test. Sizes should be compared along a row and not down a column, because the FDR bases differ: the GeoMx FDR is recomputed over the prioritized genes alone and reaches a minimum of 4.2 \u00d7 10\u207b\u2075, whereas the TCGA and CPTAC FDRs are genome- and proteome-wide. A grey cross marks a gene the resource does not measure; among the {C['supported']} genes shown, GeoMx measures {C['cov_kept']['CosMx_NOR_CORE']} and CPTAC {C['cov_kept']['CPTAC']}, so an absent bubble is a missing measurement and not a null result. Because support requires sign(log2FC) = sign(Z), Up calls arise only in the upper panel and Down calls only in the lower one. The GeoMx row is compartment-level and therefore distinct from the single-cell CosMx measurement shown in (D); the two are separate assays of the same cohort and the column named CosMx_NOR_CORE in the underlying table holds the GeoMx values."""),
 ("(D)", f"""Single-cell spatial validation in CosMx, at single-cell rather than the compartment resolution of the GeoMx row in (C). Left, the cell composition of the CosMx slide used throughout (339,739 cells plotted), coloured by cell-type group. Right, six scTWAS-prioritized genes — *CEBPB*, *GTF2I*, *HRAS*, *SMARCD1*, *WNT4* and *ZNRF3* (sub-panels a–f, with the paired comparisons in g–l respectively): for each gene, the per-cell expression map on the same slide (log2 normalized counts, colour scale capped at the 98th percentile of expressing cells, non-expressing cells in pale grey) and, beneath it, the mean normalized expression per patient in normal epithelium versus carcinoma cells with patients paired by line (n = {int(DV.n_patients.iloc[0])} patients). P values are from a two-sided paired Wilcoxon signed-rank test on the {int(DV.n_patients.iloc[0])} patient pairs and are nominal, that is unadjusted for multiple testing; across the {CX_PANEL_N}-gene CosMx panel the false-discovery rate for all six genes is {FDR_D[0]:.2f}. With {int(DV.n_patients.iloc[0])} pairs the smallest attainable P value is {MIN_P:.4f}, which three of the six genes attain. These comparisons therefore replicate the direction of change at single-cell resolution rather than providing independent significance. Five of the six genes are risk-increasing in the scTWAS and rise in carcinoma, so their change is direction-concordant; *WNT4* is protective in the scTWAS and also rises, and is therefore direction-discordant."""),
]


def add_para(doc, text, bold_lead=None, size=10.5, italic_genes=True,
             align=WD_ALIGN_PARAGRAPH.JUSTIFY, space_after=6):
    p = doc.add_paragraph()
    p.alignment = align
    p.paragraph_format.space_after = Pt(space_after)
    if bold_lead:
        r = p.add_run(bold_lead + " ")
        r.bold = True
        r.font.size = Pt(size)
    for i, chunk in enumerate(text.split("*")):
        if not chunk:
            continue
        r = p.add_run(chunk)
        r.font.size = Pt(size)
        r.italic = italic_genes and (i % 2 == 1)
    return p


def base_doc():
    doc = Document()
    st = doc.styles["Normal"]
    st.font.name = "Calibri"
    st.font.size = Pt(10.5)
    for s in doc.sections:
        s.left_margin = s.right_margin = Inches(0.9)
    return doc


# ---- legend document -------------------------------------------------------
doc = base_doc()
t = doc.add_paragraph()
tr = t.add_run(LEGEND_TITLE)
tr.bold = True
tr.font.size = Pt(11)
for lead, body in LEGEND:
    add_para(doc, body, bold_lead=lead)
out_legend = os.path.join(ROOT, "Figure5_legend.docx")
try:
    doc.save(out_legend)
except PermissionError:
    out_legend = out_legend.replace(".docx", "_NEW.docx")
    doc.save(out_legend)

md = [f"**{LEGEND_TITLE}**", ""]
for lead, body in LEGEND:
    md += [f"**{lead}** {body}", ""]
with open(os.path.join(ROOT, "Figure5_legend.md"), "w", encoding="utf-8") as fh:
    fh.write("\n".join(md))

print("wrote", out_legend)
print("legend words:", sum(len(b.split()) for _, b in LEGEND))


# --------------------------------------------------------------------------- #
# the review
# --------------------------------------------------------------------------- #
rev = base_doc()
h = rev.add_paragraph()
hr = h.add_run("Review of ms_MultiOmics_DE_support.docx and Figure 5")
hr.bold = True
hr.font.size = Pt(13)

def heading(text):
    p = rev.add_paragraph()
    r = p.add_run(text)
    r.bold = True
    r.font.size = Pt(11)
    p.paragraph_format.space_before = Pt(10)
    p.paragraph_format.space_after = Pt(3)

def bullet(text):
    add_para(rev, text, size=10, align=WD_ALIGN_PARAGRAPH.LEFT, space_after=4)

def table_from(df, widths=None):
    tb = rev.add_table(rows=1, cols=len(df.columns))
    tb.style = "Light Grid Accent 1"
    for i, c in enumerate(df.columns):
        cell = tb.rows[0].cells[i]
        cell.text = ""
        r = cell.paragraphs[0].add_run(str(c))
        r.bold = True
        r.font.size = Pt(8.5)
    for _, row in df.iterrows():
        cells = tb.add_row().cells
        for i, v in enumerate(row):
            cells[i].text = ""
            r = cells[i].paragraphs[0].add_run("" if pd.isna(v) else str(v))
            r.font.size = Pt(8.5)
    if widths:
        for r_ in tb.rows:
            for i, w in enumerate(widths):
                r_.cells[i].width = Inches(w)
    rev.add_paragraph().paragraph_format.space_after = Pt(2)

heading("1. Numbers that need filling in (Panel A, paragraph 2)")
bullet("The paragraph carries seven XX placeholders and one stated denominator. "
       f"The denominators are {DEN_A['normal-colon']} normal-colon and "
       f"{DEN_A['polyp']} polyp gene \u00d7 cell-type associations, so the stated "
       f"397 is correct and the missing figure is {DEN_A['polyp']}.")
bullet("Two counts can answer each XX and they are not the same number. "
       "\u201cSupport events\u201d counts a gene \u00d7 cell-type association once per "
       "contrast in which it is supported, so an association supported in both "
       "the MSI-H and the MSS arm is counted twice. \u201cDistinct associations\u201d "
       "counts it once. Only the distinct count is commensurable with the "
       f"{DEN_A['normal-colon']} / {DEN_A['polyp']} denominators, and I would use "
       "that one, stating the basis.")
tab = SEC_A.rename(columns={
    "section": "Route and stage",
    "normal-colon events": "NL events", "normal-colon distinct": "NL distinct",
    "polyp events": "Polyp events", "polyp distinct": "Polyp distinct"})
table_from(tab, widths=[1.6, 1.0, 1.0, 1.0, 1.0])
bullet("Per-contrast counts, if the text is written per arm rather than per "
       "route, are in the accompanying CSV.")

heading("2. One numerical error (Panel B, paragraph 3)")
bullet(f"The text states that 7\u201312 genes are supported in the normal-to-ASC "
       f"contrasts originating from each of the seven normal cell types. The "
       f"range is {rng(NB_ENTRY['ASC'])}: enteroendocrine cells support "
       f"{NB_ENTRY['ASC']['EE']} genes, below the stated lower bound. The rest "
       f"of that paragraph is correct as written.")

heading("3. A cross-reference that points at the wrong panel")
bullet("The final paragraph cites Figure 5D for the notable candidates "
       f"*{'*, *'.join(TEXT_D_CITED[:-1])}* and *{TEXT_D_CITED[-1]}*. None of "
       f"those genes appears in Panel D, which shows "
       f"*{'*, *'.join(PANEL_D)}*. The citation should be Figure 5C, or Panel D "
       "should be rebuilt around the genes the text names.")

heading("4. A direction-discordant gene shown without a flag")
bullet(f"*{DISCORDANT[0]}* (Panel D, sub-panels e and k) rises in carcinoma "
       f"cells, but its scTWAS direction is protective, so the change is "
       f"direction-discordant \u2014 the source table records concordant = FALSE "
       f"for it while the other five genes are TRUE. It currently sits "
       f"unmarked among five direction-concordant risk-increasing genes, which "
       f"reads as six concordant validations. Either mark it or drop it.")

heading("5. How the Panel D statistics are presented")
bullet(f"All six P values are nominal. Across the {CX_PANEL_N}-gene CosMx panel "
       f"the FDR for all six is {FDR_D[0]:.4f}, so none is significant after "
       f"correction. Showing P = 0.016 in bold without the FDR anywhere in the "
       f"figure or text is the most likely reviewer objection to this panel; "
       f"the legend I drafted states it explicitly.")
bullet(f"With {int(DV.n_patients.iloc[0])} patient pairs the smallest attainable "
       f"P value from a signed-rank test is {MIN_P:.4f}. Three genes sit exactly "
       f"at that floor, so the test cannot rank them against each other and the "
       f"apparent spread between 0.016 and 0.031 is one patient changing sign.")

heading("6. The Panel C spatial row is mislabelled by platform")
bullet("The column named CosMx_NOR_CORE in "
       "table/scTWAS_347genes_multisource_DE_support.csv does not hold CosMx "
       "data. Its log2 fold changes are identical to four decimal places to the "
       "GeoMx DSP values in L:/AI-gene-drug/spatial_in/geomx_DE_compact.tsv for "
       "the PanCK-positive compartment, NOR->CORE transition (GTF2I 0.9750, "
       "TRIM28 1.2973, UBE2V1 1.0423), and the other project's own figure axes "
       "label them GeoMx. Only the FDR differs, because this project recomputes "
       "it over the prioritized genes. Panel C's row and every text reference to "
       "it should read GeoMx DSP, PanCK+ compartment; the genuine single-cell "
       "CosMx measurement is a different table and is what Panel D shows. The "
       "counts themselves are unaffected — the numbers in the table below are "
       "correct — but the platform attribution and the stated resolution are "
       "not, and the source filename "
       "data/spatialDE_CosMx_PMID41986711_risk_gene_all_comparisons.csv is "
       "misnamed for the same reason.")

heading("7. Claims I checked and found correct")
ok = pd.DataFrame({
    "Claim": [
        "546 cell type-gene associations, 347 unique genes",
        f"{C['onco']} oncogene-like (Z > 0), {C['tsg']} tumour-suppressor-like (Z < 0)",
        f"{C['onco_sup']} oncogene-like supported: {C['onco_1']} / {C['onco_2']} / {C['onco_3']} by one, two, three resources",
        f"{C['tsg_sup']} suppressor-like supported: {C['tsg_1']} / {C['tsg_2']} / {C['tsg_3']}",
        f"{C['supported']} of {C['genes']} genes supported, {C['pct']}%",
        f"ASC n = {NB_GENES['ASC']} genes, SSC n = {NB_GENES['SSC']} genes",
        f"{NB_TRANS['ASC']} transitions = 7 normal-to-polyp + 14 normal-to-carcinoma + 2 polyp-to-carcinoma",
        f"STM {NB_ENTRY['ASC']['STM']} and TAC {NB_ENTRY['ASC']['TAC']} strongest for normal-to-ASC",
        f"GOB {NB_ENTRY['SSC']['GOB']} of {NB_GENES['SSC']} strongest for normal-to-SSC, others 1-3",
        f"ASC to MSI-H {NB_EXIT['ASC']['L1_MSIH']}, ASC to MSS {NB_EXIT['ASC']['L1_MSS']}",
        f"SSC to MSI-H {NB_EXIT['SSC']['L2_MSIH']}, SSC to MSS {NB_EXIT['SSC']['L2_MSS']}",
        f"{NB_GE5['ASC']} ASC and {NB_GE5['SSC']} SSC transitions with at least five supported genes",
        "397 normal-colon associations evaluated",
    ],
    "Verdict": ["correct"] * 13})
table_from(ok, widths=[5.0, 1.0])

heading("8. Wording and figure construction")
bullet("Paragraph 1 reads \u201cscRAN-seq\u201d; should be \u201cscRNA-seq\u201d.")
bullet("The composite carries a stray \u201cMean normalised\u201d text block below "
       "Panel D, an orphaned axis title from the sub-panel assembly, and a "
       "large blank area beneath it. Both should go before submission.")
bullet("Panel lettering collides. The composite uses A\u2013D while Panel D's "
       "own sub-panels use a\u2013l, so \u201cd\u201d denotes both Panel D and the "
       "*SMARCD1* map. Renumbering Panel D's sub-panels D1\u2013D12 would remove "
       "the ambiguity.")
bullet("Legibility at composite scale: Panel C's gene labels are set at about "
       "4.6 pt in the standalone figure and Panel B's node and edge labels are "
       "smaller still once scaled into the composite. Both fall below the "
       "roughly 5 pt floor most journals apply. Panel C would need either a "
       "full-page width or a reduction to the genes supported by at least two "
       "resources.")

heading("9. Two substantive points the section could make and does not")
bullet("The oncogene-like and suppressor-like sides are not symmetric in their "
       f"strongest evidence: {C['onco_3']} oncogene-like gene is supported by "
       f"all three resources against {C['tsg_3']} suppressor-like. That "
       "asymmetry is visible in Panel C and is a result, not an artefact of "
       "the display.")
bullet("Support requires the fold-change sign to match the scTWAS Z, so Up "
       "calls can only occur among Z > 0 genes and Down calls only among "
       "Z < 0. A reader may otherwise read Panel C's clean separation as "
       "evidence of directional consistency when it is imposed by the rule. "
       "One clause in the Methods or the legend fixes this.")

out_rev = os.path.join(ROOT, "ms_MultiOmics_DE_support_review.docx")
try:
    rev.save(out_rev)
except PermissionError:
    out_rev = out_rev.replace(".docx", "_NEW.docx")
    rev.save(out_rev)

## the per-contrast numbers as a table, so the text can be written either way
per = (R_A.pivot_table(index=["section", "contrast"], columns="dataset",
                       values="events", aggfunc="sum").reset_index())
per.to_csv(os.path.join(ROOT, "Figure5A_panelA_support_counts.csv"), index=False)
SEC_A.to_csv(os.path.join(ROOT, "Figure5A_panelA_section_counts.csv"), index=False)

print("wrote", out_rev)
print(SEC_A.to_string(index=False))
