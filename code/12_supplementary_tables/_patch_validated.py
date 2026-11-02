#!/usr/bin/env python3
"""Splice the validated-gene section (2.14) and its methods into the report builder."""
import re

P = "code/12_build_report.js"
s = open(P).read()

# ---------------------------------------------------------------- 1. methods 1.9
methods = '''body.push(H2("1.9  The validated gene set and CosMx per-cell expression"));
body.push(P("A risk gene is called validated when, in the GeoMx epithelium, it is significant at FDR < 0.05 in at least one of normal \\u2192 adenoma, adenoma \\u2192 carcinoma or normal \\u2192 carcinoma, with the fold change in the direction the TWAS z-score predicts, and no significant change the opposite way in the other two transitions; a gene measured on both panels must not be contradicted on either. Validated genes are labelled oncogene-like (risk-increasing and higher in the later stage) or tumour-suppressor-like (protective and lower). The number passing was compared with 10,000 permutations of the TWAS direction labels among the same genes, holding the significance pattern fixed. CosMx per-cell counts were read from the Seurat object by parsing R\\u2019s serialization format directly \\u2014 no R interpreter was available and no RDS-reading package could be installed \\u2014 streaming the compressed object and extracting only the rows of the requested genes from the sparse count matrix (script 24). Counts were normalised to the cell\\u2019s total transcripts from the panel, log2(1 + 1000 \\u00d7 count / nCount_RNA). Cells were reduced to per-patient, per-cell-type means before testing; the carcinoma-versus-normal-epithelium contrast is a Wilcoxon signed-rank test paired within patient, and the invasion trend is a Spearman correlation of patient means with stage rank."));
'''
methods += '''body.push(H2("1.10  Functional screens"));
body.push(P("Two independent high-throughput viability screens were supplied by the experimenters: siRNA knockdown of 333 genes and CRISPR-Cas9 editing of 116, each in HCT116, RKO and SW480 with three replicates and both carrying non-targeting negative and killing positive controls. Quality control per plate used the Z\u2019-factor, 1 \u2212 3(sd_pos + sd_neg)/|mean_pos \u2212 mean_neg|, the strictly standardised mean difference, the non-targeting coefficient of variation and the control signal window, plus replicate and cross-cell-line rank correlation. SW480 was excluded from the siRNA analysis because its mean Z\u2019 is at or below zero, and CRISPR was restricted to the same two lines so the modalities share a background; a sensitivity analysis restores SW480. Effects are robust z-scores against each plate\u2019s own sample median and MAD rather than its non-targeting wells, because the offset between the two is plate-specific and large; the non-targeting wells instead furnish the null. A gene\u2019s statistic is the mean z over its six wells, and its p-value is the two-sided position of that statistic in a null of 20,000 draws of one non-targeting well from each plate, corrected across the assayed panel. A detected gene clears FDR < 0.05 and moves the same way in both lines; a strong gene additionally reaches |z| \u2265 2. Enrichment of risk genes was tested against the other genes of the same screen by Fisher exact test and Mann-Whitney U. The per-cell-line figures report a different and simpler analysis of the same wells, and the same one for both assays. The non-targeting wells of each replicate plate are averaged into that plate’s control value, giving one (gene, control) pair per plate, and the three pairs are compared by a paired t-test at a nominal P < 0.05, uncorrected; the effect reported is the mean of the three within-plate differences in viability points, so pairing removes the plate-to-plate offset in the controls without a normalisation step. In the CRISPR screen those control means run from 102.0 to 105.4 within one cell line; in the siRNA screen the supplied values are already expressed against each plate’s non-targeting mean, which is exactly 100 on every plate, so the paired test there reduces to a one-sample test against the control level. With three pairs the test carries two degrees of freedom, so it is underpowered against a gene whose replicates disagree. Replicate plates were scored on six quality measures (script 42, Table S51) and none was removed: the one plate meeting the pre-specified exclusion rule would have left that cell line with a single degree of freedom, collapsing its p-values while barely changing its effects, so the concordance between plates is reported directly instead, as Pearson and Spearman correlation over the whole panel. A significant change is scored as agreeing with the TWAS when a risk-increasing gene loses viability or a protective gene gains it, and as opposing it otherwise; genes with no significant change are left uncounted, and agreement is tested against 0.5 by a two-sided binomial test, both over gene × cell-line calls and over genes. For the siRNA screen the verdict uses HCT116 and RKO only, SW480 having failed quality control."));
'''
methods += '''body.push(H2("1.11  The high-confidence list and the rule-based single-cell selection"));
body.push(P("A list of 109 high-confidence genes was supplied separately (high-confidence-genes.xlsx, one column of symbols); all 109 fall within the TWAS risk set, 93 were measured on the GeoMx panel and 41 are carried on the 6,175-gene CosMx panel. The cross-platform comparison of \\u00a72.17 restricts the analysis of \\u00a72.14 to those 41: the GeoMx bulk stage log2 fold change for the transition carrying each gene\\u2019s call, against the CosMx mean per-patient difference between carcinoma cells and the same patient\\u2019s normal epithelium, related by Spearman correlation with sign agreement counted separately. The same two statistics were computed for the remaining CosMx-tested risk genes and for the pooled set, and are reported side by side rather than tested against each other, since 41 against 80 genes cannot distinguish the two coefficients. The evidence stack counts, for the same 41 genes, how many agree in direction across platforms, how many are validated under the \\u00a72.14 rule, how many reach nominal CosMx P < 0.05, how many appear on either screening panel, and how many show a nominal viability effect by the paired t-test of \\u00a71.10. Genes drawn in Figure 30 were selected by a rule fixed before inspection: nominal CosMx P < 0.05 in the paired carcinoma-versus-normal-epithelium test; CRC-related in the sense of \\u00a72.4 (a colorectal-cancer somatic driver, or in direct protein-protein contact with one) or druggable at tier 1\\u20132; per-cell expression extractable from the Seurat object; and not already drawn in Figure 28. Seven genes qualified and the six with the largest absolute difference were plotted. Because the paired Wilcoxon over seven patients cannot return a value below 0.016, no across-gene correction is applied to these p-values and none is claimed."));
'''
methods += '''body.push(H2("1.12  Pathway enrichment by cell type"));
body.push(P("Each TWAS cell type\u2019s risk genes were tested for over-representation against the Reactome, WikiPathways and PID gene sets distributed by dhimmel/pathways, restricted to sets of 5 to 500 genes (1,626 sets over a universe of 8,968 genes). Symbols were mapped to Entrez through the same identifier map used throughout; a cell type was tested when at least five of its genes mapped into the universe, which left 12 of the 16 cell types. The statistic is the hypergeometric survival function with the pathway universe as background; every gene set a module touches at all counts as a test and enters the Benjamini-Hochberg correction, which is applied within the cell type, and a separate correction across all tests of the analysis is reported alongside. Only sets sharing at least two risk genes are reported. Two null models were run for each cell type, 1,000 draws each, scored through the identical pipeline: random modules of the same size drawn uniformly from the universe, and random modules matched gene-for-gene on the number of gene sets a gene belongs to, binned by decile of that count. The second null is the one quoted, because annotation depth is a confounder for any gene list assembled from the literature-facing end of a genetic study. A pathway was classified as canonical when its name matches a fixed family map (WNT, TGF-beta/BMP, RTK-MAPK, PI3K-AKT, p53, mismatch and other DNA repair, NOTCH, cell cycle, apoptosis, chromatin remodelling) and novel otherwise, where novel means outside those families rather than previously unstudied. Because Reactome nests its gene sets, enriched pathways sharing an identical set of risk genes are merged into one row for display and counted once when describing the result. The pooled analysis of \u00a72.19 repeats this exactly, changing only the module definition: the normal and adenoma risk genes of a cell type are combined into one module, giving nine modules of which eight are testable. The two source lists share between one and eleven genes per cell type, so a pooled module is close to a union of two independent scans; because the hypergeometric test scales the overlap against the module size, pooling raises the threshold a gene set must clear and helps only where both tissues contribute genes to the same set, which is recorded per enrichment. \u00a72.20 tests the whole risk-gene set at once and adds four analyses to the same genes: the same test split by the direction of risk; a membership audit counting how many risk genes are annotated to each canonical family whether or not any of its gene sets is enriched, tested by hypergeometric against the same universe; gene-level network proximity, the shortest path in the protein interaction graph from each risk gene to the nearest member of a canonical family, z-scored against genes of the same degree in twenty degree bins; and ten aetiology themes, each merging the gene sets whose names match a fixed pattern for a colorectal exposure mechanism (carcinogen and drug metabolism, alcohol and aldehyde handling, bile acid and lipid, one-carbon and folate, immune surveillance, oxidative stress, barrier and adhesion, energy metabolism and insulin, transport, protein turnover), tested the same way with a binomial test of the direction balance of the risk genes that fall in them. The driver annotation of \u00a72.21 uses two reference sets: the colorectal somatic driver list of \u00a72.4 (54 genes from TCGA colorectal, the COSMIC Cancer Gene Census and IntOGen), and a pan-cancer annotation curated by hand from the cancer-gene literature, because the Census requires a login and the analysis environment blocks the file hosts that mirror it. Each pan-cancer entry carries the tumour context, the somatic role and a confidence label of established, probable or contextual. Where a gene acts through its own expression level rather than through a fusion, the direction of the germline association is scored against that role: an oncogene is consistent when risk-increasing and a tumour suppressor when protective."));
'''

anchor = 'body.push(new Paragraph({ children: [new PageBreak()] }));\nbody.push(H1("2  Results"));'
assert anchor in s
s = s.replace(anchor, methods + anchor, 1)

# ---------------------------------------------------------------- 2. sections 2.14, 2.15
sec = open("code/_insert_validated.js").read()
sec += "\n\n" + open("code/_insert_screens.js").read()
sec += "\n\n" + open("code/_insert_pathways.js").read()
sec += "\n\n" + open("code/_insert_drivers.js").read()
anchor2 = 'body.push(H1("3  Limitations"));'
assert anchor2 in s
s = s.replace(anchor2, sec + "\n\n" + anchor2, 1)

# ---------------------------------------------------------------- 3. retire the stale claim in 2.13
old = ('body.push(P("What CosMx cannot do here has to be said plainly, because it was one of the '
       'planned analyses. The per-cell gene expression exists only inside a Seurat R object, and no R '
       'interpreter is available in this environment, so the risk genes cannot be read out cell by cell. '
       'Separately, CosMx experiment 2 annotates tumour, oncofetal and stromal states \\u2014 carcinoma, '
       'oncofetal, trophocyte-like CAF, myofibroblast, T cell and so on \\u2014 and not the normal epithelial '
       'lineages (absorptive, goblet, stem, tuft, enteroendocrine) in which the TWAS genes were called. Even '
       'with the expression matrix open, a gene-level test matched to each gene\\u2019s assigned cell type would '
       'not be possible in this dataset. The cell-type-resolved check therefore remains open, and \\u00a74 names '
       'the data that would close it."));')
idx = s.find('body.push(P("What CosMx cannot do here')
assert idx >= 0, "2.13 CosMx-limitation paragraph not found"
end = s.index('\n', s.index('"));', idx))
new = ('body.push(P("One limit of this dataset is structural rather than technical. CosMx experiment 2 '
       'annotates tumour, oncofetal and stromal states — carcinoma, oncofetal, trophocyte-like CAF, '
       'myofibroblast, T cell and so on — and not the normal epithelial lineages (absorptive, goblet, '
       'stem, tuft, enteroendocrine) in which the TWAS genes were called, so a gene-level test matched to '
       'each gene’s own lineage is not possible here whatever is done with the expression matrix. The '
       'per-cell expression itself is available: §2.14 reads it out of the Seurat object and uses it to '
       'test the validated genes cell by cell."));')
s = s[:idx] + new + s[end:]

# ---------------------------------------------------------------- 4. limitations
old_lim = ('  "Per-cell CosMx expression sits inside a Seurat R object and no R interpreter was available, '
           'so the risk genes could not be read out cell by cell. CosMx experiment 2 also annotates tumour '
           'and stromal states rather than normal epithelial lineages, so even with the matrix open a '
           'lineage-matched gene-level test would not be possible in that dataset.",')
assert old_lim in s
new_lim = ('  "CosMx experiment 2 annotates tumour and stromal states rather than normal epithelial '
           'lineages, so the single-cell tests compare carcinoma cells with normal epithelium and cannot '
           'be matched to the lineage each gene was called in. The cohort also starts at intramucosal '
           'carcinoma, so it cannot speak to the normal → adenoma step.",\n'
           '  "The validated set is a filter, not an enrichment: the number of genes passing is no greater '
           'than permuting the direction labels would give. Its justification is the independent CosMx '
           'replication, not the count itself, and a gene excluded for never reaching significance may '
           'simply be under-powered on the panel rather than genuinely unchanged.",\n'
           '  "Eleven CosMx patients across three invasion stages cannot support per-gene multiple-testing '
           'correction; the paired test floors at P = 0.016. Only the across-gene direction tests are '
           'interpretable, and no single-gene single-cell claim should be read as established.",\n'
           '  "ERBB2 fails spatial validation while remaining the best pharmacological candidate. The two '
           'lines of evidence are independent and should be reported as such rather than reconciled.",')
s = s.replace(old_lim, new_lim, 1)

# ---------------------------------------------------------------- 5. summary bullet
anchor3 = '''  [{ text: "The pharmacology converges on ERBB2. ", bold: true },'''
assert anchor3 in s
add = '''  [{ text: "Sixty-two genes carry a coherent, replicated direction. ", bold: true },
   `Requiring the somatic change in tissue to agree with the genetic direction leaves ${D.validated.n_validated} of the ${D.validated.n_measured} measured risk genes \\u2014 ${D.validated.n_oncogene_like} oncogene-like and ${D.validated.n_tsg_like} tumour-suppressor-like. Read out of the CosMx Seurat object cell by cell, ${n2(100 * D.cosmx_percell.crc_vs_normal.frac_concordant_validated, 0)} % of these move the predicted way between carcinoma cells and the patient\\u2019s own normal epithelium (P = ${fmtP(D.cosmx_percell.crc_vs_normal.p_binomial_validated)}) against ${n2(100 * D.cosmx_percell.crc_vs_normal.frac_concordant_all, 0)} % for the unfiltered list \\u2014 eleven patients on a different platform, in neither GeoMx cohort. CTNNB1 gains genetics, driver status and tissue expression at once; ERBB2 does not validate.`],
'''
s = s.replace(anchor3, add + anchor3, 1)


# ---------------------------------------------------------------- 6. cross-reference 2.8
old28 = ('body.push(P(`Pathway enrichment recovers coherent, cell-type-appropriate biology '
         '(${D.n_pathway_sig} enrichments at FDR < 0.1), which is a useful sanity check on '
         'the cell-type assignments themselves.`));')
assert old28 in s
new28 = ('body.push(P(`Pathway enrichment recovers coherent, cell-type-appropriate biology '
         '(${D.n_pathway_sig} enrichments at FDR < 0.1), which is a useful sanity check on '
         'the cell-type assignments themselves. \u00a72.18 runs this analysis again with the '
         'controls it needs \u2014 a null matched on how heavily each gene is annotated \u2014 '
         'and the table below should be read against that section, which removes the adenoma '
         'Wnt and cell-cycle block and leaves the antigen-presentation and detoxification '
         'signals standing.`));')
s = s.replace(old28, new28, 1)

# ---------------------------------------------------------------- 7. pathway limitations
anchor_lim = ('  "CosMx experiment 2 annotates tumour and stromal states rather than normal '
              'epithelial ')
assert anchor_lim in s
pw_lim = ('  "Twenty-five of the 347 risk genes carry a cancer-driver annotation, but only '
          'three are on the colorectal somatic driver list. The pan-cancer annotations for '
          'the other twenty-two were curated by hand from the cancer-gene literature rather '
          'than pulled from a COSMIC Cancer Gene Census release, which this environment could '
          'not reach, and the ten labelled probable or contextual should be re-checked against '
          'a current release.",\n'
          '  "Pathway enrichment of the cell-type risk genes is carried by 35 genes in total, '
          'and three of them (HLA-A, HLA-F, HLA-G) sit in the MHC, a long-range '
          'linkage-disequilibrium region where a TWAS cannot separate neighbouring genes. '
          'The antigen-presentation result should be read as one locus, not three genes, '
          'until fine-mapping says otherwise. Pooling the normal and adenoma genes of each '
          'cell type (\u00a72.19) was run as a sensitivity analysis and did not change this: '
          'larger modules produced fewer canonical hits, not more.",\n'
          '  "The pathway background is the set of pathway-annotated genes rather than the '
          'set of genes the TWAS could have discovered. The annotation-matched null corrects '
          'for how heavily a gene is studied but not for which genes had usable expression '
          'models, so a residual ascertainment bias remains. Over-representation also '
          'ignores effect size and direction: nothing in \u00a72.18 says which way a pathway '
          'pushes risk.",\n')
s = s.replace(anchor_lim, pw_lim + anchor_lim, 1)

open(P, "w").write(s)
print("patched", P)
