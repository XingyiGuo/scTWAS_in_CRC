#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Reformat "GWAS meta.docx" (W:\\VUMC_projects\\scTWAS\\ms) into a Nature Genetics
Methods subsection. The GWAS meta-analysis is a minor result whose product is the
summary statistics used for the TWAS, so the Methods text is held to three
paragraphs and the two cohort tables are merged into a single table.

Outputs (written next to this script):
  GWAS_meta_Methods_NG.docx / .md      3-paragraph Methods subsection + references
  GWAS_meta_Table_NG.xlsx              single combined cohort table (+ flat sheet)
  GWAS_meta_Tables_NG.docx             same combined table, Word three-line style
  GWAS_meta_reformat_notes.docx / .md  numbered account of every change + open flags

Inline markup understood by add_rich():
  {i:...} italic   {b:...} bold   {s:...} superscript
"""
import os, re, copy
from docx import Document
from docx.shared import Pt, RGBColor, Inches
from docx.enum.text import WD_ALIGN_PARAGRAPH, WD_LINE_SPACING
from docx.enum.table import WD_TABLE_ALIGNMENT
from docx.oxml.ns import qn
from docx.oxml import OxmlElement

OUT = os.environ.get("GWAS_MS_OUT") or os.path.dirname(os.path.abspath(__file__))

# ----------------------------------------------------------------- doc helpers
def new_doc(font="Times New Roman", size=12, double=True, line_numbers=False):
    doc = Document()
    st = doc.styles["Normal"]
    st.font.name = font
    st.font.size = Pt(size)
    st.element.rPr.rFonts.set(qn("w:eastAsia"), font)
    pf = st.paragraph_format
    pf.space_after = Pt(0)
    if double:
        pf.line_spacing_rule = WD_LINE_SPACING.DOUBLE
    for s in doc.sections:
        s.left_margin = s.right_margin = Inches(1.0)
        if line_numbers:
            ln = OxmlElement("w:lnNumType")
            ln.set(qn("w:countBy"), "1")
            ln.set(qn("w:restart"), "continuous")
            s._sectPr.append(ln)
    return doc

RICH = re.compile(r"\{([ibs]):([^{}]*)\}")

def add_rich(par, text, size=None, bold=False, italic=False, color=None):
    pos = 0
    for m in RICH.finditer(text):
        if m.start() > pos:
            _run(par, text[pos:m.start()], size, bold, italic, color)
        kind, val = m.group(1), m.group(2)
        r = _run(par, val, size, bold, italic, color)
        if kind == "i":
            r.italic = True
        elif kind == "b":
            r.bold = True
        elif kind == "s":
            r.font.superscript = True
        pos = m.end()
    if pos < len(text):
        _run(par, text[pos:], size, bold, italic, color)
    return par

def _run(par, txt, size, bold, italic, color):
    r = par.add_run(txt)
    if size:
        r.font.size = Pt(size)
    if bold:
        r.bold = True
    if italic:
        r.italic = True
    if color:
        r.font.color.rgb = color
    return r

def para(doc, text="", style=None, size=None, bold=False, italic=False,
         space_before=0, space_after=0, double=True, align=None, color=None,
         first_line_indent=None):
    p = doc.add_paragraph(style=style)
    pf = p.paragraph_format
    pf.space_before = Pt(space_before)
    pf.space_after = Pt(space_after)
    pf.line_spacing_rule = (WD_LINE_SPACING.DOUBLE if double
                            else WD_LINE_SPACING.SINGLE)
    if align is not None:
        p.alignment = align
    if first_line_indent is not None:
        pf.first_line_indent = Inches(first_line_indent)
    if text:
        add_rich(p, text, size=size, bold=bold, italic=italic, color=color)
    return p

def heading(doc, text, level=1, double=True):
    sizes = {1: 14, 2: 12, 3: 12}
    p = para(doc, "", space_before=12 if level == 1 else 10, space_after=4,
             double=double)
    r = p.add_run(text)
    r.bold = True
    r.font.size = Pt(sizes[level])
    if level == 3:
        r.italic = True
    return p

FLAG = RGBColor(0xC0, 0x00, 0x00)

def flag_para(doc, text, double=True, size=None):
    """[TO CONFIRM] line rendered bold dark red."""
    p = para(doc, "", space_before=2, space_after=2, double=double)
    add_rich(p, text, size=size, bold=True, color=FLAG)
    return p

# --------------------------------------------------------------- table helpers
def _border(tcPr, edge, sz=8, val="single"):
    tag = "w:" + edge
    el = tcPr.find(qn("w:tcBorders"))
    if el is None:
        el = OxmlElement("w:tcBorders")
        tcPr.append(el)
    b = el.find(qn(tag))
    if b is None:
        b = OxmlElement(tag)
        el.append(b)
    b.set(qn("w:val"), val)
    b.set(qn("w:sz"), str(sz))
    b.set(qn("w:space"), "0")
    b.set(qn("w:color"), "000000")

def three_line(table, header_rows=1):
    """Nature-style rule set: rule above header, below header, below last row."""
    n = len(table.rows)
    for i, row in enumerate(table.rows):
        for cell in row.cells:
            tcPr = cell._tc.get_or_add_tcPr()
            for edge in ("top", "bottom", "left", "right"):
                _border(tcPr, edge, val="none", sz=0)
            if i == 0:
                _border(tcPr, "top", sz=12)
            if i == header_rows - 1:
                _border(tcPr, "bottom", sz=8)
            if i == n - 1:
                _border(tcPr, "bottom", sz=12)

def make_table(doc, header, rows, widths=None, size=9, bold_last_row=False):
    t = doc.add_table(rows=1 + len(rows), cols=len(header))
    t.alignment = WD_TABLE_ALIGNMENT.CENTER
    t.autofit = True
    for j, h in enumerate(header):
        c = t.cell(0, j)
        c.text = ""
        p = c.paragraphs[0]
        p.paragraph_format.space_before = Pt(2)
        p.paragraph_format.space_after = Pt(2)
        p.paragraph_format.line_spacing_rule = WD_LINE_SPACING.SINGLE
        add_rich(p, h, size=size, bold=True)
    for i, row in enumerate(rows, start=1):
        for j, v in enumerate(row):
            c = t.cell(i, j)
            c.text = ""
            p = c.paragraphs[0]
            p.paragraph_format.space_before = Pt(1)
            p.paragraph_format.space_after = Pt(1)
            p.paragraph_format.line_spacing_rule = WD_LINE_SPACING.SINGLE
            if j >= len(header) - 2:
                p.alignment = WD_ALIGN_PARAGRAPH.RIGHT
            add_rich(p, str(v), size=size,
                     bold=(bold_last_row and i == len(rows)))
    three_line(t, header_rows=1)
    if widths:
        for j, w in enumerate(widths):
            for row in t.rows:
                row.cells[j].width = Inches(w)
    return t

# =============================================================== SOURCE CONTENT
# ---- Supplementary Table 1: rows grouped into two explicit blocks so that the
#      dataset composition of each final GWAS reads off the table directly.
T1_TITLE = ("Colorectal cancer GWAS datasets contributing to the final European-ancestry and "
            "multi-ancestry meta-analyses.")
T1_GROUPS = [("", 4),
             ("Final EUR GWAS meta-analysis", 2),
             ("Final multi-ancestry GWAS meta-analysis", 2)]
T1_HEADER = ["Study", "Cohort", "Phenotype definition", "Genetic ancestry",
             "Cases ({i:n})", "Controls ({i:n})", "Cases ({i:n})", "Controls ({i:n})"]
DASH = "\u2014"
NOT_INC = "Not included"

# cohort, definition, ancestry, eur_cases, eur_ctrl, multi_cases, multi_ctrl
BLOCK_A_LABEL = ("Datasets included in BOTH the final EUR GWAS and the final multi-ancestry GWAS "
                 "(5 datasets)")
BLOCK_A = [
    ["Previously published CRC GWAS meta-analysis{s:a}", "Study-specific{s:a}", "EUR + EAS",
     78473, 107142, 100204, 154587],
    ["BioVU", "ICD-based algorithm{s:b}", "EUR", 2581, 144746, 2581, 144746],
    ["MVP", "Phecode 153.2 (colon cancer)", "EUR", 9125, 438096, 9125, 438096],
    ["FinnGen (R13)", "C3_COLON endpoint{s:c}", "EUR", 8259, 372159, 8259, 372159],
    ["UCLA ATLAS", "Phecode 153.2 (colon cancer)", "EUR", 822, 57369, 822, 57369],
]
BLOCK_B_LABEL = ("Datasets included in the final multi-ancestry GWAS ONLY (8 datasets; not in the "
                 "final EUR GWAS)")
BLOCK_B = [
    ["BioVU", "ICD-based algorithm{s:b}", "AFR", None, None, 347, 26048],
    ["MVP", "Phecode 153.2 (colon cancer)", "AFR", None, None, 2069, 118718],
    ["MVP", "Phecode 153.2 (colon cancer)", "AMR", None, None, 870, 58633],
    ["UCLA ATLAS", "Phecode 153.2 (colon cancer)", "AMR", None, None, 122, 10381],
    ["UCLA ATLAS", "Phecode 153.2 (colon cancer)", "EAS", None, None, 107, 7347],
    ["UCLA ATLAS", "Phecode 153.2 (colon cancer)", "AFR", None, None, 59, 3686],
    ["MGI (Freeze 6)", "Phecode 153.2 (colon cancer)", "Mixed{s:d}", None, None, 1967, 57860],
    ["TPMI", "Phecode 153.2 (colon cancer)", "EAS", None, None, 5729, 302768],
]
T1_ROWS = BLOCK_A + BLOCK_B
T1_WIDTHS = [0.4, 1.7, 1.55, 0.8, 0.72, 0.85, 0.72, 0.85]

def t1_totals():
    """(eur_cases, eur_ctrl, multi_cases, multi_ctrl) for block A, block B and overall."""
    s = lambda rows, k: sum(r[k] for r in rows if r[k] is not None)
    A = (s(BLOCK_A, 3), s(BLOCK_A, 4), s(BLOCK_A, 5), s(BLOCK_A, 6))
    B = (s(BLOCK_B, 3), s(BLOCK_B, 4), s(BLOCK_B, 5), s(BLOCK_B, 6))
    T = tuple(x + y for x, y in zip(A, B))
    assert B[0] == 0 and B[1] == 0, B
    assert T[:2] == (99260, 1119512), T          # final EUR GWAS, as in the source document
    assert T[2:] == (132261, 1752398), T         # final multi-ancestry GWAS
    return A, B, T

T1_FOOTNOTES = [
    "Each row is one contributing GWAS (cohort \u00d7 genetic ancestry). The two column blocks give "
    "the sample contributed to each final meta-analysis; '" + NOT_INC + "' marks a dataset that does "
    "not enter the final EUR GWAS. Totals equal the exact sum of the rows above them.",
    "{s:a}Previously published multi-ancestry CRC GWAS meta-analysis{s:7}; case and control "
    "definitions follow the contributing studies of that report. Its combined European and East "
    "Asian release (100,204 cases, 154,587 controls) entered the multi-ancestry meta-analysis and "
    "its European subset (78,473 cases, 107,142 controls) entered the EUR meta-analysis; this is the "
    "only dataset whose contribution differs between the two.",
    "{s:b}Cases were ascertained from de-identified electronic health records using ICD diagnosis "
    "codes with the exclusions described in the Methods.",
    "{s:c}FinnGen endpoint C3_COLON, 'Malignant neoplasm of colon (including HILMO)', defined "
    "primarily from ICD-O-3 topography and morphology codes in the Finnish Cancer Registry and "
    "supplemented by harmonized ICD-8, ICD-9 and ICD-10 diagnosis codes.",
    "{s:d}MGI Freeze 6 summary statistics were generated in a combined sample of multiple genetic "
    "ancestry groups and are therefore labelled Mixed.",
    "AFR, African; AMR, admixed American; CRC, colorectal cancer; EAS, East Asian; EUR, European; "
    "ICD, International Classification of Diseases; MGI, Michigan Genomics Initiative; MVP, Million "
    "Veteran Program; TPMI, Taiwan Precision Medicine Initiative; UCLA ATLAS, University of "
    "California, Los Angeles ATLAS Community Health Initiative.",
]
T1_FLAG = ("[TO CONFIRM] The control count for the previously published European-ancestry CRC GWAS "
           "is given as 107,142 in the source document; the published value is 107,143, which would "
           "make the final EUR GWAS control total 1,119,513 and would also reconcile with the "
           "154,587 controls used in the multi-ancestry column (107,143 EUR + 47,444 EAS).")
T1_FLAT_COLS = ["study_no", "block", "cohort", "phenotype_definition", "genetic_ancestry",
                "in_final_eur_gwas", "in_final_multi_gwas",
                "eur_gwas_cases", "eur_gwas_controls", "multi_gwas_cases", "multi_gwas_controls"]

# ---------------------------------------------------------------- METHODS TEXT
# Three paragraphs maximum: the GWAS meta-analysis is a minor result whose
# product is the summary statistics used as input to the TWAS.
# (heading level, kind, text)  kind: p | run-in | flag | h
METHODS = [
 (2, "h", "Colorectal cancer GWAS meta-analysis"),

 (0, "p", "We first performed a CRC GWAS in BioVU, the Vanderbilt University Medical Center DNA "
   "databank, which links DNA samples to de-identified electronic health records (EHRs){s:1}. CRC "
   "cases were ascertained from ICD diagnosis codes ([TO CONFIRM: insert ICD-9-CM and ICD-10-CM code "
   "list]), excluding individuals with a family history of CRC, inflammatory bowel disease, "
   "colorectal polyps, familial adenomatous polyposis or Lynch syndrome, and individuals with "
   "another cancer diagnosed more than 6 months before their CRC diagnosis; controls had no history "
   "of cancer. Samples were genotyped on a custom Illumina Multi-Ethnic Genotyping Array (MEGA-EX) "
   "and sequenced through the Alliance for Genomic Discovery. We excluded samples with less than 90% "
   "of the genome covered at a depth of at least 10\u00d7, transition-to-transversion ratios more "
   "than four s.d. from the cohort mean, or evidence of duplication or sample mismatch, and retained "
   "only sites passing DRAGEN filters with genotype quality \u2265 20, read depth \u2265 10, allele "
   "balance between 0.20 and 0.80 for heterozygous calls and genotype missingness below 10%. "
   "Genotypes were imputed on the Michigan Imputation Server v1.2.4 using Minimac4 (ref. 2) against "
   "the Haplotype Reference Consortium panel v1.1 (ref. 3), following the BioVU imputation "
   "pipeline{s:4}. Ancestry-specific association testing was performed with regenie v4 (ref. 5) "
   "using sex, age and the first five genetic principal components as covariates: step 1 fitted a "
   "whole-genome ridge regression model on LD-pruned biallelic SNPs (1,000-variant window, "
   "100-variant step, {i:r}{s:2} threshold 0.45) with minor allele frequency (MAF) \u2265 1%, minor "
   "allele count \u2265 100, missingness \u2264 10% and Hardy\u2013Weinberg equilibrium {i:P} "
   "\u2265 1 \u00d7 10{s:\u221215} prepared in PLINK 2 (ref. 6), and step 2 tested each variant "
   "conditional on leave-one-chromosome-out predictions, re-testing variants with {i:P} < 0.01 by "
   "approximate Firth logistic regression. This yielded 2,581 cases and 144,746 controls of European "
   "(EUR) ancestry and 347 cases and 26,048 controls of African (AFR) ancestry."),

 (0, "p", "We combined the BioVU GWAS with published and biobank CRC GWAS summary statistics: our "
   "previous multi-ancestry CRC meta-analysis{s:7}, the Million Veteran Program (MVP){s:8}, FinnGen "
   "(R13){s:9}, the University of California, Los Angeles ATLAS Community Health Initiative (UCLA "
   "ATLAS){s:10}, the Michigan Genomics Initiative (MGI) Freeze 6 (ref. 11) and the Taiwan Precision "
   "Medicine Initiative (TPMI){s:12}. In the EHR-linked biobanks (MVP, UCLA ATLAS, MGI and TPMI), "
   "cases were defined by mapping ICD-9-CM and ICD-10-CM diagnosis codes to phecodes and taking "
   "colon cancer phecode 153.2; in FinnGen, cases were defined by the C3_COLON endpoint (malignant "
   "neoplasm of colon, including HILMO) from ICD-O-3 topography and morphology codes in the Finnish "
   "Cancer Registry, supplemented by harmonized ICD-8, ICD-9 and ICD-10 codes. Cohort-specific "
   "sample sizes, phenotype definitions and genetic ancestry groups are given in Supplementary "
   "Table 1."),

 (0, "p", "Summary statistics reported on GRCh37 were converted to GRCh38 with UCSC LiftOver{s:13} "
   "and unmapped variants were discarded; within each dataset, duplicate variants were resolved by "
   "retaining the record with the smallest {i:P} value, and analyses were restricted to autosomal "
   "variants with MAF > 1%. We then performed inverse-variance-weighted fixed-effects meta-analysis "
   "using METAL (release 2011-03-25){s:14}, applying genomic control correction to each contributing "
   "GWAS and to the meta-analysed summary statistics and assessing between-study heterogeneity with "
   "Cochran's {i:Q} test and the {i:I}{s:2} statistic. The EUR-ancestry meta-analysis comprised "
   "99,260 cases and 1,119,512 controls from BioVU, our previous GWAS, MVP, FinnGen (R13) and UCLA "
   "ATLAS; the multi-ancestry meta-analysis additionally incorporated AFR- and admixed American "
   "(AMR)-ancestry participants from BioVU, MVP and UCLA ATLAS, East Asian (EAS)-ancestry "
   "participants from UCLA ATLAS and TPMI, and the mixed-ancestry MGI sample, totalling 132,261 "
   "cases and 1,752,398 controls (Supplementary Table 1). These two sets of meta-analysis summary "
   "statistics provided the GWAS input for the transcriptome-wide association analyses described "
   "below."),
]

ACK = [
 (2, "h", "Cohort acknowledgements (for the manuscript Acknowledgements section, not Methods)"),
 (0, "p", "The authors thank the Million Veteran Program (MVP) staff, researchers and, most "
   "importantly, the veterans whose participation made this research possible."),
 (0, "p", "We acknowledge the participants and investigators of the FinnGen study."),
 (0, "p", "We gratefully acknowledge the support of the Institute for Precision Health, the "
   "participating patients of the UCLA ATLAS Precision Health Biobank, the UCLA David Geffen School "
   "of Medicine, the UCLA Clinical and Translational Science Institute (grant UL1TR001881) and UCLA "
   "Health."),
 (0, "p", "We acknowledge the Michigan Genomics Initiative participants, AI & Digital Health "
   "Innovation at the University of Michigan, the University of Michigan Medical School Central "
   "Biorepository and the University of Michigan Advanced Genomics Core for providing data and "
   "specimen storage, management, processing and distribution services, and the Center for "
   "Statistical Genetics in the Department of Biostatistics at the School of Public Health for "
   "genotype data curation, imputation and management in support of the research reported in this "
   "publication."),
 (0, "p", "We thank all participants and investigators of the Taiwan Precision Medicine Initiative "
   "(TPMI), and we thank Academia Sinica and the cooperating medical centres across Taiwan for their "
   "support in establishing and maintaining the cohort datasets."),
]

# Nature style: first author + et al. for >5 authors; sentence-case titles;
# abbreviated italic journal; bold volume; en-dash page range; year in parentheses.
REFERENCES = [
 "Ritchie, M. D. {i:et al.} Robust replication of genotype\u2013phenotype associations across "
 "multiple diseases in an electronic medical record. {i:Am. J. Hum. Genet.} {b:86}, 560\u2013572 (2010).",
 "Das, S. {i:et al.} Next-generation genotype imputation service and methods. {i:Nat. Genet.} "
 "{b:48}, 1284\u20131287 (2016).",
 "McCarthy, S. {i:et al.} A reference panel of 64,976 haplotypes for genotype imputation. "
 "{i:Nat. Genet.} {b:48}, 1279\u20131283 (2016).",
 "Piekos, J. A. {i:et al.} Evaluating the relationships between genetic ancestry and the clinical "
 "phenome. {i:Pac. Symp. Biocomput.} {b:29}, 389\u2013403 (2024).",
 "Mbatchou, J. {i:et al.} Computationally efficient whole-genome regression for quantitative and "
 "binary traits. {i:Nat. Genet.} {b:53}, 1097\u20131103 (2021).",
 "Chang, C. C. {i:et al.} Second-generation PLINK: rising to the challenge of larger and richer "
 "datasets. {i:Gigascience} {b:4}, 7 (2015).",
 "Fernandez-Rozadilla, C. {i:et al.} Deciphering colorectal cancer genetics through multi-omic "
 "analysis of 100,204 cases and 154,587 controls of European and east Asian ancestries. "
 "{i:Nat. Genet.} {b:55}, 89\u201399 (2023).",
 "Verma, A. {i:et al.} Diversity and scale: genetic architecture of 2,068 traits in the VA Million "
 "Veteran Program. {i:Science} {b:385}, eadj1182 (2024).",
 "Kurki, M. I. {i:et al.} FinnGen provides genetic insights from a well-phenotyped isolated "
 "population. {i:Nature} {b:613}, 508\u2013518 (2023).",
 "Johnson, R. {i:et al.} Leveraging genomic diversity for discovery in an electronic health record "
 "linked biobank: the UCLA ATLAS Community Health Initiative. {i:Genome Med.} {b:14}, 104 (2022).",
 "Zawistowski, M. {i:et al.} The Michigan Genomics Initiative: a biobank linking genotypes and "
 "electronic clinical records in Michigan Medicine patients. {i:Cell Genom.} {b:3}, 100257 (2023).",
 "Yang, H.-C. {i:et al.} The Taiwan Precision Medicine Initiative provides a cohort for large-scale "
 "studies. {i:Nature} {b:648}, 117\u2013127 (2025).",
 "Hinrichs, A. S. {i:et al.} The UCSC Genome Browser Database: update 2006. {i:Nucleic Acids Res.} "
 "{b:34}, D590\u2013D598 (2006).",
 "Willer, C. J., Li, Y. & Abecasis, G. R. METAL: fast and efficient meta-analysis of genomewide "
 "association scans. {i:Bioinformatics} {b:26}, 2190\u20132191 (2010).",
]

REF_MAP = [
 ("1", "Ritchie et al. (BioVU)", "1", "unchanged"),
 ("2", "Piekos et al. (BioVU genetic ancestry / pipeline)", "4", "renumbered"),
 ("3", "Fernandez-Rozadilla et al. (previous CRC GWAS)", "7", "renumbered"),
 ("4", "Verma et al. (MVP)", "8", "renumbered"),
 ("5", "Kurki et al. (FinnGen)", "9", "renumbered"),
 ("6", "Johnson et al. (UCLA ATLAS)", "10", "renumbered"),
 ("7", "Zawistowski et al. (MGI)", "11", "renumbered"),
 ("8", "Yang et al. (TPMI)", "12", "renumbered"),
 ("9", "Hinrichs et al. (UCSC LiftOver)", "13", "renumbered"),
 ("10", "Willer et al. (METAL)", "14", "renumbered"),
 ("\u2014", "Das et al. (Minimac4 / Michigan Imputation Server)", "2", "added"),
 ("\u2014", "McCarthy et al. (HRC reference panel v1.1)", "3", "added"),
 ("\u2014", "Mbatchou et al. (regenie)", "5", "added"),
 ("\u2014", "Chang et al. (PLINK 2)", "6", "added"),
]

# ======================================================== build Methods (.docx)
def build_methods():
    doc = new_doc(line_numbers=True)
    for level, kind, text in METHODS + ACK:
        if kind == "h":
            heading(doc, text, level=level)
        elif kind == "flag":
            flag_para(doc, text)
        elif kind == "run-in":
            lead, body = text.split("|", 1)
            p = para(doc, "", space_after=0)
            r = p.add_run(lead)
            r.bold = True
            add_rich(p, body)
        else:
            para(doc, text)

    heading(doc, "References", level=2)
    for i, ref in enumerate(REFERENCES, start=1):
        p = para(doc, "", space_after=0, double=True)
        p.paragraph_format.left_indent = Inches(0.3)
        p.paragraph_format.first_line_indent = Inches(-0.3)
        add_rich(p, "%d.\t%s" % (i, ref))

    doc.save(os.path.join(OUT, "GWAS_meta_Methods_NG.docx"))

# ========================================= build combined table (.docx + .xlsx)
def _fmt(v, eur_block=False):
    if v is None:
        return NOT_INC if eur_block else DASH
    return format(v, ",d")

def _t1_render_rows():
    """Flat render plan: ('block'|'row'|'sub'|'tot', cells)."""
    A, B, T = t1_totals()
    out = [("block", BLOCK_A_LABEL)]
    n = 0
    for r in BLOCK_A:
        n += 1
        out.append(("row", [str(n), r[0], r[1], r[2],
                            _fmt(r[3], True), _fmt(r[4], True), _fmt(r[5]), _fmt(r[6])]))
    out.append(("sub", ["", "Subtotal, datasets in both GWAS", "", "",
                        _fmt(A[0]), _fmt(A[1]), _fmt(A[2]), _fmt(A[3])]))
    out.append(("block", BLOCK_B_LABEL))
    for r in BLOCK_B:
        n += 1
        out.append(("row", [str(n), r[0], r[1], r[2],
                            _fmt(r[3], True), _fmt(r[4], True), _fmt(r[5]), _fmt(r[6])]))
    out.append(("sub", ["", "Subtotal, multi-ancestry-only datasets", "", "",
                        DASH, DASH, _fmt(B[2]), _fmt(B[3])]))
    out.append(("tot", ["", "Total in each final GWAS", "", "",
                        _fmt(T[0]), _fmt(T[1]), _fmt(T[2]), _fmt(T[3])]))
    return out

def build_tables():
    plan = _t1_render_rows()
    doc = new_doc(double=False)
    for s in doc.sections:
        s.orientation = 1
        s.page_width, s.page_height = s.page_height, s.page_width
        s.left_margin = s.right_margin = Inches(0.7)
    p = para(doc, "", space_after=6, double=False)
    r = p.add_run("Supplementary Table 1 | ")
    r.bold = True
    add_rich(p, T1_TITLE)

    ncol = len(T1_HEADER)
    t = doc.add_table(rows=2 + len(plan), cols=ncol)
    t.alignment = WD_TABLE_ALIGNMENT.CENTER
    col = 0
    for label, span in T1_GROUPS:
        if span > 1:
            merged = t.cell(0, col).merge(t.cell(0, col + span - 1))
            mp = merged.paragraphs[0]
            mp.alignment = WD_ALIGN_PARAGRAPH.CENTER
            mp.paragraph_format.line_spacing_rule = WD_LINE_SPACING.SINGLE
            add_rich(mp, label, size=9, bold=True)
        col += span
    for j, h in enumerate(T1_HEADER):
        cp = t.cell(1, j).paragraphs[0]
        cp.paragraph_format.line_spacing_rule = WD_LINE_SPACING.SINGLE
        if j >= 4:
            cp.alignment = WD_ALIGN_PARAGRAPH.CENTER
        add_rich(cp, h, size=9, bold=True)

    for i, (kind, cells) in enumerate(plan, start=2):
        if kind == "block":
            merged = t.cell(i, 0).merge(t.cell(i, ncol - 1))
            mp = merged.paragraphs[0]
            mp.paragraph_format.space_before = Pt(3)
            mp.paragraph_format.space_after = Pt(2)
            mp.paragraph_format.line_spacing_rule = WD_LINE_SPACING.SINGLE
            add_rich(mp, cells, size=9, bold=True, italic=True)
            continue
        for j, v in enumerate(cells):
            cp = t.cell(i, j).paragraphs[0]
            cp.paragraph_format.space_before = Pt(1)
            cp.paragraph_format.space_after = Pt(1)
            cp.paragraph_format.line_spacing_rule = WD_LINE_SPACING.SINGLE
            if j >= 4:
                cp.alignment = WD_ALIGN_PARAGRAPH.RIGHT
            add_rich(cp, str(v), size=9, bold=(kind in ("sub", "tot")),
                     italic=(v == NOT_INC))
    three_line(t, header_rows=2)
    for j, w in enumerate(T1_WIDTHS):
        for row in t.rows:
            row.cells[j].width = Inches(w)
    for fn in T1_FOOTNOTES:
        para(doc, fn, size=8, space_before=4, double=False)
    flag_para(doc, T1_FLAG, double=False, size=8)
    doc.save(os.path.join(OUT, "GWAS_meta_Tables_NG.docx"))

def build_xlsx():
    from openpyxl import Workbook
    from openpyxl.styles import Font, Alignment, Border, Side, PatternFill
    from openpyxl.utils import get_column_letter
    A, B, T = t1_totals()
    SUP = {"a": "\u1d43", "b": "\u1d47", "c": "\u1d9c", "d": "\u1d48", "n": "n",
           "7": "\u2077", "2": "\u00b2"}
    plain = lambda s: RICH.sub(
        lambda m: ("".join(SUP.get(ch, ch) for ch in m.group(2)) if m.group(1) == "s"
                   else m.group(2)), str(s))
    bare = lambda s: RICH.sub(lambda m: "" if m.group(1) == "s" else m.group(2), str(s))

    wb = Workbook()
    ws = wb.active
    ws.title = "Supplementary Table 1"
    thin = Side(style="thin", color="000000")
    med = Side(style="medium", color="000000")
    band = PatternFill("solid", fgColor="EDEDED")
    ncol = len(T1_HEADER)

    c = ws.cell(row=1, column=1, value="Supplementary Table 1 | " + plain(T1_TITLE))
    c.font = Font(bold=True, size=11)
    ws.merge_cells(start_row=1, start_column=1, end_row=1, end_column=ncol)

    col = 1
    for label, span in T1_GROUPS:
        if span > 1:
            ws.merge_cells(start_row=3, start_column=col, end_row=3, end_column=col + span - 1)
            g = ws.cell(row=3, column=col, value=label)
            g.font = Font(bold=True)
            g.alignment = Alignment(horizontal="center")
            for k in range(col, col + span):
                ws.cell(row=3, column=k).border = Border(top=med, bottom=thin)
                ws.cell(row=3, column=k).fill = band
        else:
            for k in range(col, col + span):
                ws.cell(row=3, column=k).border = Border(top=med)
        col += span
    for j, h in enumerate(T1_HEADER, start=1):
        h2 = ws.cell(row=4, column=j, value=plain(h))
        h2.font = Font(bold=True)
        h2.alignment = Alignment(horizontal="center" if j >= 5 else "left",
                                 vertical="bottom", wrap_text=True)
        h2.border = Border(bottom=thin)

    rn = 5
    for kind, cells in _t1_render_rows():
        if kind == "block":
            ws.merge_cells(start_row=rn, start_column=1, end_row=rn, end_column=ncol)
            bc = ws.cell(row=rn, column=1, value=plain(cells))
            bc.font = Font(bold=True, italic=True)
            for k in range(1, ncol + 1):
                ws.cell(row=rn, column=k).fill = band
            rn += 1
            continue
        for j, v in enumerate(cells, start=1):
            v2 = plain(v)
            if j >= 5 and v2 not in (NOT_INC, DASH, ""):
                cell = ws.cell(row=rn, column=j, value=int(v2.replace(",", "")))
                cell.number_format = "#,##0"
            else:
                cell = ws.cell(row=rn, column=j, value=v2)
            cell.alignment = Alignment(horizontal="right" if j >= 5 else "left")
            cell.font = Font(bold=(kind in ("sub", "tot")),
                             italic=(v2 == NOT_INC),
                             color="808080" if v2 == NOT_INC else "000000")
            if kind == "tot":
                cell.border = Border(top=thin, bottom=med)
        rn += 1

    fr = rn + 1
    for k, fn in enumerate(T1_FOOTNOTES + [T1_FLAG]):
        fc = ws.cell(row=fr + k, column=1, value=plain(fn))
        fc.font = Font(size=9, bold=(k == len(T1_FOOTNOTES)),
                       color="C00000" if k == len(T1_FOOTNOTES) else "000000")
        fc.alignment = Alignment(vertical="top")
        ws.merge_cells(start_row=fr + k, start_column=1, end_row=fr + k, end_column=ncol)
    for j, w in enumerate([7, 42, 30, 16, 14, 15, 14, 15], start=1):
        ws.column_dimensions[get_column_letter(j)].width = w
    ws.row_dimensions[4].height = 30
    ws.freeze_panes = "B5"

    wsf = wb.create_sheet("Table_S1_flat")
    wsf.append(T1_FLAT_COLS)
    for c in wsf[1]:
        c.font = Font(bold=True)
    n = 0
    for blk, rows in (("both", BLOCK_A), ("multi_only", BLOCK_B)):
        for r in rows:
            n += 1
            wsf.append([n, blk, bare(r[0]), bare(r[1]), r[2],
                        "yes" if r[3] is not None else "no", "yes",
                        r[3], r[4], r[5], r[6]])
    wsf.append([None, "total", "Total in each final GWAS", None, None, None, None,
                T[0], T[1], T[2], T[3]])
    for row in wsf.iter_rows(min_row=2, min_col=8, max_col=11):
        for c in row:
            c.number_format = "#,##0"
    for j, w in enumerate([9, 12, 42, 30, 16, 18, 19, 16, 18, 18, 20], start=1):
        wsf.column_dimensions[get_column_letter(j)].width = w
    wsf.freeze_panes = "A2"
    wb.save(os.path.join(OUT, "GWAS_meta_Table_NG.xlsx"))

# ========================================================== build notes (.docx)
NOTES_TABLES = [
 "Made the dataset composition of each final GWAS explicit: rows are grouped into two labelled, "
 "shaded blocks inside the table \u2014 'Datasets included in BOTH the final EUR GWAS and the "
 "final multi-ancestry GWAS (5 datasets)' and 'Datasets included in the final multi-ancestry GWAS "
 "ONLY (8 datasets; not in the final EUR GWAS)' \u2014 each with its own subtotal, above a final "
 "'Total in each final GWAS' row. The two column blocks are headed 'Final EUR GWAS "
 "meta-analysis' and 'Final multi-ancestry GWAS meta-analysis', and cells in the EUR block read "
 "'Not included' in grey italic rather than a dash, so no cell is ambiguous. The generic "
 "'Phenotype' column was dropped (its content is carried by 'Phenotype definition') to make room.",
 "Merged the two cohort tables into a single Supplementary Table 1 with one row per contributing "
 "GWAS (cohort \u00d7 genetic ancestry) and the two meta-analyses side by side under grouped "
 "headers, 'Multi-ancestry meta-analysis' and 'European-ancestry meta-analysis'. An em dash marks "
 "a dataset that does not contribute to a given meta-analysis, so the European column shows at a "
 "glance which samples were dropped. This keeps the previously published GWAS on one row while "
 "still carrying its combined EUR + EAS counts for the multi-ancestry column and its EUR subset for "
 "the European column.",
 "Delivered the merged table as GWAS_meta_Table_NG.xlsx with two sheets: 'Supplementary Table 1' "
 "(formatted for the supplement \u2014 title row, grouped headers, thousands separators, "
 "right-aligned counts, Nature three-line rules, frozen panes and footnotes below the table) and "
 "'Table_S1_flat' (the same values with unique machine-readable column names for downstream code). "
 "The same merged table is also written to GWAS_meta_Tables_NG.docx in landscape for pasting into "
 "Word.",
 "Deleted the internal working note from the original table banner ('The sample info is shown "
 "below: (for BioVU data may need to double check with Jing). Metal conducted by Shuoshuo on Aug "
 "2026'), along with the banner rows whose text repeated across all seven columns and the blank "
 "spacer row.",
 "Removed the stray token 'darkcyan' from the MGI phenotype label ('Colon cancer neoplasms "
 "darkcyan' \u2192 'Colon cancer').",
 "Replaced the 'Disease'/'Code' columns with 'Phenotype' and 'Phenotype definition' so that "
 "'Phecode 153.2' and 'C3_COLON' sit in a column whose header explains them; the long FinnGen "
 "endpoint description moved from a table cell to a footnote.",
 "Renamed 'Meta of CRC (NG)' to 'Previously published CRC GWAS meta-analysis' with a footnote "
 "citing the source publication and stating that the multi-ancestry column uses the combined "
 "European and East Asian release while the European column uses the European subset.",
 "Standardized cohort naming and data freezes ('FinnGen (R13)', 'MGI (Freeze 6)', 'UCLA ATLAS', "
 "'MVP', 'TPMI', 'BioVU'), defined every ancestry abbreviation in a footnote (including 'Mixed' for "
 "the jointly analysed MGI sample), reformatted the count headers as 'Cases ({i:n})' / 'Controls "
 "({i:n})' and added thousands separators throughout.",
 "Verified both sets of totals against the rows in code (the build asserts them): multi-ancestry "
 "132,261 cases / 1,752,398 controls and European-ancestry 99,260 cases / 1,119,512 controls each "
 "equal the exact sum of their columns. The 'Total' row is bold.",
]

NOTES_TEXT = [
 "Condensed the Methods to three paragraphs, on the basis that the GWAS meta-analysis is a minor "
 "result whose product is the summary statistics used for the TWAS: (1) the BioVU GWAS, from "
 "ascertainment through QC, imputation and regenie association testing, ending with its EUR and AFR "
 "sample sizes; (2) the external published and biobank GWAS with their phenotype definitions, "
 "pointing to Supplementary Table 1 for sample sizes; (3) harmonization, the METAL meta-analyses "
 "and the two resulting sample sizes, closing by stating that these summary statistics are the GWAS "
 "input to the TWAS. The subsection is now headed 'Colorectal cancer GWAS meta-analysis'.",
 "Dropped the per-cohort descriptive paragraphs. The MVP, FinnGen, UCLA ATLAS, MGI and TPMI "
 "paragraphs in the source described each biobank's design, participant counts and phenotyping "
 "framework; only the case definition actually used (phecode 153.2, or C3_COLON for FinnGen) is "
 "retained, one clause each, and the sample sizes now live in Supplementary Table 1 rather than "
 "being repeated in prose.",
 "Dropped the Ethics statement, Reporting summary, Data availability and Code availability "
 "sections. These are required once per manuscript rather than per Methods subsection; if this "
 "GWAS is instead submitted as a standalone paper they must be reinstated (see open items).",
 "Moved the five cohort acknowledgements out of the Methods and into a clearly labelled appendix at "
 "the end of the file, since their wording is mandated by the MVP, FinnGen, UCLA ATLAS, MGI and "
 "TPMI data-use agreements and must not be lost in the condensation.",
 "Corrected a factual error in the source TPMI paragraph, which read '5,729 EUR CRC cases and "
 "302,768 controls from FinnGen (R13) were included in the EUR-specific meta-analysis'. TPMI is an "
 "East Asian-ancestry Taiwanese cohort and, per the cohort table, contributes those counts to the "
 "multi-ancestry meta-analysis only. Cohort name, ancestry label and meta-analysis label are all "
 "corrected, and TPMI is now listed among the multi-ancestry-only datasets in paragraph 3.",
 "Added the software and reference-panel citations Nature Genetics requires (regenie, PLINK 2, "
 "Minimac4 / Michigan Imputation Server and the HRC reference panel) and renumbered the reference "
 "list to run in order of first citation; the old \u2192 new mapping is tabulated below for the "
 "EndNote update.",
 "Converted all references to Nature style: surname with spaced initials, first author plus italic "
 "'et al.' where there are more than five authors, sentence-case titles, abbreviated italic journal "
 "names, bold volume, en-dash page ranges and year in parentheses. 'D590-598' was corrected to "
 "'D590\u2013D598' and '2068 traits' given a thousands separator.",
 "Used Nature's '(ref. n)' citation form wherever a superscript reference number would abut a "
 "numeral \u2014 'regenie v4 (ref. 5)', 'PLINK 2 (ref. 6)', 'Minimac4 (ref. 2)', 'HRC panel v1.1 "
 "(ref. 3)', 'MGI Freeze 6 (ref. 11)'. This is the ambiguity that turned into 'v1.2.410' and "
 "'REGENIE 4' in the source document.",
 "Applied Nature house style throughout: '{i:P} value' with italic {i:P} rather than 'P-value'; "
 "italic {i:I}{s:2}, {i:Q} and {i:r}{s:2}; 's.d.' for standard deviation; en dashes in "
 "'Hardy\u2013Weinberg' and numeric ranges; a true minus sign in 1 \u00d7 10{s:\u221215}; "
 "'\u00d7' for sequencing depth; thousands separators in all counts; and 'fixed-effects' added to "
 "the description of the inverse-variance-weighted meta-analysis.",
 "Replaced 'Variants on chromosomes 23\u201327 were excluded' with restriction to biallelic "
 "autosomal SNPs \u2014 the 23\u201327 numbering is PLINK-internal coding and is not informative "
 "to a reader. The step-1 description, which appeared twice in the source, was consolidated, and "
 "'MEGA-ex' \u2192 'MEGA-EX' with the vendor address dropped, as Nature does not carry vendor city "
 "and state.",
 "Set the Methods file for submission: 12-pt Times New Roman, double-spaced, 1-inch margins and "
 "continuous line numbers.",
]

NOTES_FLAGS = [
 "BioVU ICD codes. The source has the placeholder '[xxxxx]' where the CRC case-ascertainment codes "
 "belong. Supply the ICD-9-CM and ICD-10-CM lists (and the phecode, if used).",
 "BioVU WGS versus array data. The population paragraph describes WGS through AGD and WGS-specific "
 "variant filters (DRAGEN, GQ, DP, allele balance), but the pipeline then imputes MEGA-EX array "
 "genotypes against HRC v1.1. State which data the CRC GWAS used, and if both, how they were "
 "combined.",
 "Previously published European control count. This document gives 78,473 cases and 107,142 "
 "controls; the published figure is 107,143 controls, and 107,143 + 47,444 East Asian controls is "
 "exactly the 154,587 used in the multi-ancestry row. If 107,143 is correct, the European-ancestry "
 "control total becomes 1,119,513. The source numbers have been left in place pending your "
 "confirmation; the tables and text are internally consistent either way.",
 "Michigan Imputation Server version. 'v1.2.410' in the source has been rendered as 'v1.2.4'; "
 "confirm the version and whether a citation was intended by the trailing '10'.",
 "Reference 2 of the old list (Piekos et al.) was cited at the HRC panel sentence. It has been kept "
 "in that sentence as the BioVU imputation-pipeline citation, with the primary Minimac4 and HRC "
 "citations added. Confirm this is what was intended.",
 "Ethics, data availability and code availability. These statements were dropped from this "
 "subsection on the assumption that the parent manuscript carries them once for all analyses. Confirm that the parent Methods covers: BioVU IRB approval and consent model plus the data-use "
 "approvals of MVP, FinnGen, UCLA ATLAS, MGI and TPMI; a repository and accession for the CRC "
 "meta-analysis summary statistics generated here and access routes for each contributing "
 "dataset; and the software list (PLINK 2, regenie v4, Michigan Imputation Server v1.2.4 with "
 "Minimac4, METAL release 2011-03-25, UCSC LiftOver) with any in-house harmonization scripts "
 "deposited and cited.",
 "Genomic control and heterogeneity metrics. Report per-study and meta-analysis \u03bb{s:GC} and/or "
 "LDSC intercepts with quantile\u2013quantile plots, and state the heterogeneity threshold applied "
 "to Cochran's {i:Q} / {i:I}{s:2}. Also confirm METAL was run with STDERR weighting, as the text "
 "now states.",
 "UCLA ATLAS data type. The paragraph notes that 117,417 participants had whole-exome sequencing "
 "data; confirm whether the ATLAS CRC GWAS used array genotypes rather than WES.",
 "TPMI phenotype map. Confirm that phecode 153.2 was applied using TPMI's own phecode map, and "
 "confirm the 565,390 participant count for the freeze used.",
 "METAL release. '2011-03-25' is retained; confirm this is the release actually used.",
 "Table placement. The merged table is formatted as Supplementary Table 1; renumber it if the "
 "supplement already has tables, or relabel it Table 1 if it is to sit in the main text.",
 "Cohort acknowledgements are appended to the Methods file as a separate section \u2014 move them "
 "into the manuscript Acknowledgements. Funding, author-contributions and competing-interests "
 "statements are not present in the source document.",
]

def build_notes():
    doc = new_doc(double=False)
    p = para(doc, "", space_after=8, double=False)
    r = p.add_run("Reformatting notes: 'GWAS meta.docx' \u2192 Nature Genetics Methods format")
    r.bold = True
    r.font.size = Pt(14)
    para(doc, "Source: W:\\VUMC_projects\\scTWAS\\ms\\GWAS meta.docx (unmodified). "
              "Deliverables: GWAS_meta_Methods_NG.docx and GWAS_meta_Tables_NG.docx in the same "
              "folder. Every change made to the source is listed below; items flagged "
              "[TO CONFIRM] need an author decision and are flagged in the deliverables as well.",
         size=10, double=False, space_after=8)
    for title, items in (("Tables", NOTES_TABLES), ("Methods text", NOTES_TEXT)):
        heading(doc, title, level=2, double=False)
        for i, it in enumerate(items, start=1):
            p = para(doc, "", double=False, space_after=4)
            p.paragraph_format.left_indent = Inches(0.35)
            p.paragraph_format.first_line_indent = Inches(-0.35)
            add_rich(p, "%d.\t%s" % (i, it), size=10)
    heading(doc, "Open items requiring author input", level=2, double=False)
    for i, it in enumerate(NOTES_FLAGS, start=1):
        p = para(doc, "", double=False, space_after=4)
        p.paragraph_format.left_indent = Inches(0.35)
        p.paragraph_format.first_line_indent = Inches(-0.35)
        add_rich(p, "%d.\t%s" % (i, it), size=10, bold=True, color=FLAG)
    heading(doc, "Reference renumbering (old \u2192 new)", level=2, double=False)
    make_table(doc, ["Old no.", "Reference", "New no.", "Action"], REF_MAP,
               widths=[0.7, 3.6, 0.7, 1.0], size=9)
    doc.save(os.path.join(OUT, "GWAS_meta_reformat_notes.docx"))

# ============================================================= markdown mirrors
def md(text):
    text = RICH.sub(lambda m: ("*%s*" if m.group(1) == "i" else
                               "**%s**" if m.group(1) == "b" else "^%s^") % m.group(2), text)
    return text

def build_md():
    L = ["# Methods (Nature Genetics format)", "",
         "_Reformatted from `GWAS meta.docx`; see `GWAS_meta_reformat_notes.md` for every change._", ""]
    for level, kind, text in METHODS + ACK:
        if kind == "h":
            L += ["", "#" * (level + 1) + " " + text, ""]
        elif kind == "flag":
            L += ["> **" + md(text) + "**", ""]
        elif kind == "run-in":
            lead, body = text.split("|", 1)
            L += ["**" + lead.strip() + "** " + md(body), ""]
        else:
            L += [md(text), ""]
    L += ["", "## References", ""]
    L += ["%d. %s" % (i, md(r)) for i, r in enumerate(REFERENCES, start=1)]
    L += ["", "## Reference renumbering (old -> new)", "",
          "| Old no. | Reference | New no. | Action |", "|---|---|---|---|"]
    L += ["| %s | %s | %s | %s |" % r for r in REF_MAP]

    hdr = [md(h) for h in T1_HEADER]
    L += ["", "---", "", "## Supplementary Table 1 | " + md(T1_TITLE), "",
          "Columns 5-6 = sample contributed to the **final EUR GWAS meta-analysis**; "
          "columns 7-8 = sample contributed to the **final multi-ancestry GWAS meta-analysis**.", "",
          "| " + " | ".join(hdr) + " |", "|" + "---|" * len(hdr)]
    for kind, cells in _t1_render_rows():
        if kind == "block":
            L += ["| " + " | ".join(["**_" + md(cells) + "_**"] + [""] * (len(hdr) - 1)) + " |"]
        elif kind == "row":
            L += ["| " + " | ".join(md(str(v)) for v in cells) + " |"]
        else:
            L += ["| " + " | ".join("**" + md(str(v)) + "**" if v else "" for v in cells) + " |"]
    L += [""] + [md(f) for f in T1_FOOTNOTES] + ["", "**" + md(T1_FLAG) + "**", ""]
    open(os.path.join(OUT, "GWAS_meta_Methods_NG.md"), "w", encoding="utf-8").write("\n".join(L))

    N = ["# Reformatting notes: `GWAS meta.docx` -> Nature Genetics Methods format", "",
         "Source: `W:\\VUMC_projects\\scTWAS\\ms\\GWAS meta.docx` (unmodified). Deliverables: "
         "`GWAS_meta_Methods_NG.docx`, `GWAS_meta_Tables_NG.docx`.", ""]
    N += ["## Tables", ""] + ["%d. %s" % (i, md(t)) for i, t in enumerate(NOTES_TABLES, 1)]
    N += ["", "## Methods text", ""] + ["%d. %s" % (i, md(t)) for i, t in enumerate(NOTES_TEXT, 1)]
    N += ["", "## Open items requiring author input", ""]
    N += ["%d. **%s**" % (i, md(t)) for i, t in enumerate(NOTES_FLAGS, 1)]
    N += ["", "## Reference renumbering (old -> new)", "",
          "| Old no. | Reference | New no. | Action |", "|---|---|---|---|"]
    N += ["| %s | %s | %s | %s |" % r for r in REF_MAP]
    open(os.path.join(OUT, "GWAS_meta_reformat_notes.md"), "w", encoding="utf-8").write("\n".join(N))

if __name__ == "__main__":
    build_methods()
    build_tables()
    build_xlsx()
    build_notes()
    build_md()
    for f in ("GWAS_meta_Methods_NG.docx", "GWAS_meta_Methods_NG.md",
              "GWAS_meta_Table_NG.xlsx", "GWAS_meta_Tables_NG.docx",
              "GWAS_meta_reformat_notes.docx", "GWAS_meta_reformat_notes.md"):
        p = os.path.join(OUT, f)
        print("%-36s %8d bytes" % (f, os.path.getsize(p)))
