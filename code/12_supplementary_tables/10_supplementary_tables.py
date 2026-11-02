#!/usr/bin/env python3
"""10 - Assemble the supplementary workbook."""
import pandas as pd, numpy as np, os
from openpyxl.styles import Font, PatternFill, Alignment
from openpyxl.utils import get_column_letter

OUT = "deliverables/Supplementary_Tables.xlsx"
os.makedirs("deliverables", exist_ok=True)

risk = pd.read_csv("results/risk_genes_harmonized.csv")
imp = pd.read_csv("results/gene_ko_impact.csv")
rev = pd.read_csv("results/signature_reversal_all.csv.gz")

sheets = [
    ("S1 Risk genes", risk.drop(columns=["ensembl"], errors="ignore"),
     "Harmonised cell-type-specific CRC risk genes from the two TWAS files, with the "
     "directional call used throughout (risk-increasing vs protective)."),
    ("S2 Gene KO impact", imp[imp["is_risk_gene"]].sort_values("ko_impact", ascending=False),
     "Network-propagation knockout impact per risk gene, with degree, influence reach and "
     "the number of genes experimentally dysregulated in the matching LINCS knockdown."),
    ("S3 Module KO impact", pd.read_csv("results/module_knockout_impact.csv"),
     "Module-dismantling virtual knockout: total-variation shift in the cell-type risk "
     "module's diffusion profile when each gene is deleted from the network."),
    ("S4 KO vs CRC state", pd.read_csv("results/ko_direction_vs_crc_state.csv"),
     "Correlation of each gene's experimental knockdown / overexpression signature with an "
     "independent colorectal-cancer disease signature."),
    ("S5 Pathways", pd.read_csv("results/pathway_enrichment.csv"),
     "Hypergeometric pathway enrichment of each cell-type risk module "
     "(Reactome, WikiPathways, PID)."),
    ("S6 Reversal top20", pd.read_csv("results/signature_reversal_top20.csv"),
     "Top signature-reversal candidates per cell type (all three criteria satisfied)."),
    ("S7 Reversal recurrent", pd.read_csv("results/signature_reversal_recurrent.csv"),
     "Compounds nominated as reversers in more than one cell type."),
    ("S8 Positive control", pd.read_csv("results/reversal_positive_control.csv").head(300),
     "Positive-control screen against an independent CRC disease signature "
     "(top 300 of 1,170 compounds)."),
    ("S9 Druggability", pd.read_csv("results/gene_druggability.csv"),
     "Ligandability tier of every risk gene, with the direction of pharmacological action "
     "the genetics calls for."),
    ("S10 Proximity", pd.read_csv("results/network_proximity_significant.csv"),
     "Compounds whose target set is significantly proximal to a cell-type risk module "
     "(z < -1.645 and FDR < 0.05)."),
    ("S11 DTI predictions", pd.read_csv("results/dti_predictions_validated_regime.csv"),
     "NRLMF drug-target predictions restricted to the regime the model was validated in "
     "(risk genes that already have ligand evidence)."),
    ("S12 Tractability", pd.read_csv("results/gene_tractability_score.csv"),
     "Network distance from each risk gene to the multi-ligand druggable proteome, "
     "degree-matched z-score."),
    ("S13 Gene-drug pairs", pd.read_csv("results/gene_drug_candidates.csv"),
     "Direction-annotated gene-drug pairs with the number of independent evidence layers."),
    ("S14 Target priority", pd.read_csv("results/integrated_gene_priority.csv"),
     "Integrated priority score across genetic, knockout, direction-consistency and "
     "druggability evidence."),
    ("S15 GDSC", pd.read_csv("results/gdsc_proximity_bydrug.csv"),
     "GDSC compounds: proximity to the risk modules and selectivity for the 44 colorectal "
     "cell lines (delta ln IC50, CRC minus all other lineages)."),
    ("S16 KO impact spectrum", pd.read_csv("results/ko_impact_spectrum_full.csv"),
     "Network knockout impact on one common ranking for CRC somatic drivers, core-essential "
     "and non-essential reference genes, the risk genes and a random background."),
    ("S17 Driver proximity", pd.read_csv("results/risk_gene_driver_proximity.csv"),
     "Network distance from each risk gene to the nearest CRC somatic driver, with a "
     "degree-matched z-score."),
    ("S18 Driver perturbation", pd.read_csv("results/driver_perturbation_by_ko.csv"),
     "How often each CRC driver is significantly dysregulated by risk-gene knockdowns "
     "versus background knockdowns."),
    ("S19 Module-pathway prox", pd.read_csv("results/module_pathway_proximity.csv"),
     "Degree-matched proximity of each cell-type risk module to the canonical CRC "
     "pathways."),
    ("S20 KO downstream paths", pd.read_csv("results/ko_downstream_pathways.csv"),
     "Pathways enriched among the genes each risk-gene knockdown actually dysregulates."),
    ("S21 KO vs background paths", pd.read_csv("results/ko_downstream_vs_background.csv"),
     "The control: hit rate of every pathway across 81 risk-gene knockdowns versus 713 "
     "background knockdowns."),
    ("S22 Module convergence", pd.read_csv("results/ko_module_convergence.csv"),
     "Within-module versus outside-module similarity of knockdown signatures, per module."),
    ("S23 Lineage proximity", pd.read_csv("results/lineage_network_proximity.csv"),
     "Degree-matched proximity of each risk module to each intestinal lineage-identity "
     "programme."),
    ("S24 Lineage own vs other", pd.read_csv("results/lineage_own_vs_other.csv"),
     "Each module's proximity to its own lineage programme versus the others."),
    ("S25 Stage contrast", pd.read_csv("results/stage_contrast.csv"),
     "Normal-only, adenoma-only and shared risk genes compared on every knockout property."),
    ("S26 Epistasis", pd.read_csv("results/epistasis_controlled.csv"),
     "Double-deletion synergy for risk-gene pairs within each module, z-scored against "
     "degree-matched random pairs deleted from the same module."),
    ("S27 Spatial DE risk genes", pd.read_csv("results/spatial_risk_gene_DE.csv"),
     "Stage-transition differential expression of every CRC risk gene in the GeoMx "
     "spatial data, both cohorts, both compartments, with a per-gene verdict on whether "
     "the change matches the TWAS direction."),
    ("S28 Spatial concordance", pd.read_csv("results/spatial_direction_concordance.csv"),
     "Direction-concordance test per cohort, compartment and transition, against a "
     "TWAS-direction permutation and against background genes."),
    ("S29 Spatial vs background", pd.read_csv("results/spatial_risk_vs_background.csv"),
     "Are risk genes differentially expressed more often than the rest of the genome?"),
    ("S30 Spatial replication", pd.read_csv("results/spatial_replication.csv"),
     "Agreement of stage-transition fold changes between the two independent GeoMx "
     "cohorts (different patients, different probe panels)."),
    ("S31 CosMx cell-type trend", pd.read_csv("results/cosmx_celltype_trend.csv"),
     "Cell-type proportion across intramucosal -> T1sm1 -> T1sm3 in 1,257,627 CosMx "
     "single cells from 11 patients."),
    ("S32 CosMx niche trend", pd.read_csv("results/cosmx_niche_trend.csv"),
     "Spatial niche proportion across the same invasive progression."),
    ("S33 Validation verdict", pd.read_csv("results/validated_risk_genes_all.csv"),
     "Every risk gene measured in GeoMx epithelium, its fold change at each stage "
     "transition, and whether it passes the validation rule: significant and in the "
     "TWAS-predicted direction, with no significant change the opposite way."),
    ("S34 Validated priority", pd.read_csv("results/validated_priority.csv"),
     "The validated genes ranked, each classified oncogene-like or tumour-suppressor-like, "
     "with the action that classification implies and the druggability it carries."),
    ("S35 Validation by cell type", pd.read_csv("results/validated_by_celltype.csv"),
     "Validation rate for the risk genes of each TWAS cell type, split by class."),
    ("S36 CosMx carcinoma vs normal", pd.read_csv("results/cosmx_crc_vs_normal_epi.csv"),
     "Per-cell CosMx expression of each risk gene in carcinoma cells against the normal "
     "epithelium of the same patient, paired across patients."),
    ("S37 CosMx stage by cell type", pd.read_csv("results/cosmx_gene_stage_trend.csv"),
     "Risk-gene expression across intramucosal -> T1sm1 -> T1sm3, tested across patients "
     "within each annotated cell type."),
    ("S38 CosMx-GeoMx replication", pd.read_csv("results/cosmx_geomx_replication.csv"),
     "The same genes measured on two platforms in two patient cohorts: GeoMx "
     "normal->carcinoma fold change against the CosMx single-cell carcinoma/normal "
     "difference."),
    ("S39 Summary statistics", pd.read_csv("results/risk_gene_summary_statistics.csv"),
     "Summary statistics for every one of the 271 tested risk genes: log2 fold change, "
     "nominal p, risk-restricted FDR and genome-wide FDR for each stage transition in "
     "the epithelium, the leading comparison that gives the gene its smallest FDR, the "
     "stroma as a specificity check, and the CosMx single-cell difference where the gene "
     "is on that panel. No gene is filtered out."),
    ("S40 All comparisons", pd.read_csv("results/risk_gene_all_comparisons.csv"),
     "The long form of S39: one row per gene x GeoMx panel x compartment x transition, "
     "with log2 fold change, nominal p, both FDRs and the direction verdict."),
    ("S41 Screen QC", pd.read_csv("results/func_qc.csv"),
     "Quality control of the siRNA and CRISPR-Cas9 viability screens, per cell line and "
     "replicate: non-targeting control mean and CV, killing-control mean, signal window, "
     "Z'-factor and SSMD. SW480 is excluded from the siRNA analysis on these numbers."),
    ("S42 siRNA gene results", pd.read_csv("results/func_sirna_gene_stats.csv"),
     "Per-gene siRNA result over HCT116 and RKO: robust z against the plate's own sample "
     "distribution, viability as % of non-targeting control, empirical p from the plates' "
     "non-targeting wells, FDR, per-line effects and the detected / strong calls."),
    ("S43 CRISPR gene results", pd.read_csv("results/func_crispr_gene_stats.csv"),
     "The same for the CRISPR-Cas9 screen of 116 genes."),
    ("S44 Risk genes in the screens", pd.read_csv("results/func_risk_gene_results.csv"),
     "The 62 CRC risk genes assayed by siRNA, with their TWAS direction and spatial "
     "class, the in-silico knockout impact, the measured siRNA effect and the matching "
     "CRISPR result where the gene was also assayed that way."),
    ("S45 siRNA vs CRISPR", pd.read_csv("results/func_sirna_vs_crispr.csv"),
     "The 88 genes assayed by both modalities, side by side, with the direction-agreement "
     "flag used to decide where CRISPR is reported."),
    ("S46 CRISPR by cell line", pd.read_csv("results/crispr_by_cellline.csv"),
     "The CRISPR-Cas9 screen analysed separately in HCT116, RKO and SW480: robust z, "
     "viability as % of non-targeting control, empirical p and FDR computed within each "
     "line across the 116 panel genes, for every gene."),
    ("S47 Risk genes by cell line", pd.read_csv("results/crispr_risk36_by_cellline.csv"),
     "The 36 S1 risk genes on the CRISPR panel, wide by cell line, with the number of "
     "lines supporting each gene and the direction those lines agree on. SW480 is "
     "included here because it passes QC in the CRISPR screen, unlike in siRNA."),
    ("S48 siRNA by cell line", pd.read_csv("results/sirna_by_cellline.csv"),
     "The siRNA screen analysed separately in HCT116, RKO and SW480: robust z, viability "
     "as % of non-targeting control, the empirical p and FDR within each line across the "
     "333 panel genes, and the Welch t-test of the three replicates against that line's "
     "non-targeting wells, for every gene. qc_pass marks the two lines whose plates pass "
     "quality control; SW480 does not."),
    ("S49 siRNA risk genes by line", pd.read_csv("results/sirna_risk_by_cellline.csv"),
     "The 62 S1 risk genes on the siRNA panel, wide by cell line, carrying all three "
     "criteria (uncorrected empirical p, empirical FDR, nominal t-test) with counts "
     "restricted to the QC-passing lines alongside the counts over all three."),
    ("S50 Non-risk CRISPR genes", pd.read_csv("results/nontwas_by_cellline.csv"),
     "The 80 genes on the CRISPR panel that are not CRC risk genes from this study - the "
     "experiment's own comparison group - wide by cell line, with the change in viability "
     "against the non-targeting controls, the Welch t-test p-value per line, the number of "
     "lines showing a significant change, and a flag for the core-machinery genes that act "
     "as the assay's positive controls."),
    ("S51 CRISPR replicate QC", pd.read_csv("results/crispr_replicate_qc.csv"),
     "Quality control for each of the nine CRISPR replicate plates: Z'-factor, SSMD, "
     "non-targeting CV, control signal window, the plate's mean correlation with the other "
     "two plates of its own cell line across the 116 genes, and its offset from the line's "
     "median plate, with the per-metric failure flags and the exclusion decision. A plate "
     "is excluded on a hard fail (Z' at or below zero) or on three or more soft fails. "
     "No plate is removed from the reported analysis; see Figure 23."),
    ("S52 Plate correlations", pd.concat([
        pd.read_csv("results/crispr_plate_correlations.csv").assign(assay="CRISPR"),
        pd.read_csv("results/sirna_plate_correlations.csv").assign(assay="siRNA")]),
     "Pearson and Spearman correlation between every pair of replicate plates within each "
     "assay, computed on the effect scale (each gene's value minus its own plate's "
     "non-targeting mean) over the whole assayed panel, with a flag for pairs from the "
     "same cell line. Drawn as Figures 23 and 24."),
    ("S53 Genes in both screens", pd.read_csv("results/shared_genes_both_assays.csv"),
     "The 30 risk genes assayed by both modalities, with the paired change in viability and "
     "the significance flag for each of the six plates, the per-assay TWAS verdict, whether "
     "the gene is called in both screens, one or neither, and whether the two effects share "
     "a sign. Drawn as Figure 25."),
    ("S54 Risk-increasing genes", pd.read_csv("results/risk_increasing_both_assays.csv"),
     "The 40 TWAS risk-increasing genes assayed in either screen, with the paired change in "
     "viability and the significance flag for every plate of both assays, a flag for whether "
     "each assay carried the gene, and the per-gene verdict on whether its significant "
     "changes lower viability as the germline direction predicts. Drawn as Figure 26."),
    ("S55 Risk-increasing CRISPR", pd.read_csv("results/risk_increasing_crispr.csv"),
     "The 24 TWAS risk-increasing genes on the CRISPR-Cas9 panel: the paired change in "
     "viability, the paired t-test p-value and the significance flag for each of the three "
     "cell lines, the mean effect and the per-gene outcome. Drawn as Figure 27."),
    ("S56 Requested CosMx genes", pd.read_csv("results/cosmx_requested_genes.csv"),
     "The eight genes requested for single-cell display: for the four carried by the CosMx "
     "panel, the per-patient mean expression in normal epithelium and in carcinoma cells, "
     "the paired difference, the Wilcoxon p, how many of the seven patients rise, the TWAS "
     "direction, the §2.14 class and the GeoMx lead FDR; the other four are listed with "
     "on_cosmx_panel = False because the panel does not carry them. Drawn as Figure 28."),
    ("S57 High-confidence replication", pd.read_csv("results/highconf_replication.csv"),
     "The 41 genes of the supplied high-confidence list that are carried on the CosMx panel: "
     "the GeoMx bulk stage log2 fold change, the CosMx per-patient means in normal epithelium "
     "and carcinoma cells, their paired difference, the Wilcoxon p and FDR, the TWAS "
     "direction, the §2.14 validation call, whether the two platforms agree in sign, and "
     "whether the gene was assayed in a viability screen and showed a nominal effect. "
     "Drawn as Figure 29."),
    ("S58 Rule-selected CosMx genes", pd.read_csv("results/cosmx_selected_genes.csv"),
     "The 35 CosMx-tested risk genes reaching nominal P < 0.05, scored against the rest of "
     "the Figure 30 selection rule: the paired "
     "carcinoma-versus-normal-epithelium difference and its nominal p, the TWAS direction and "
     "whether the single-cell change agrees with it, CRC somatic-driver status and "
     "driver-contact status, druggability tier, whether per-cell expression could be "
     "extracted, whether the gene was already drawn in Figure 28, and the resulting "
     "eligibility flag. The six plotted genes are the eligible ones with the largest "
     "absolute difference. Drawn as Figure 30."),
    ("S59 Pathway tests", pd.read_csv("results/pathway_celltype_enrichment.csv"),
     "Every pathway over-representation test that reached at least two risk genes: cell "
     "type, gene set and source, the canonical CRC family it was assigned to (or novel), "
     "the size of the gene set and of the cell type's module, the overlap, the "
     "hypergeometric p, the Benjamini-Hochberg q within the cell type and across the "
     "whole analysis, and the genes shared. Drawn as Figure 31."),
    ("S60 Pathway nulls", pd.read_csv("results/pathway_celltype_null.csv"),
     "Per cell type: the number of risk genes mapped into the pathway universe, the "
     "number of gene sets the module generates as tests, the observed counts of enriched "
     "sets at FDR < 0.1 and at nominal p < 0.01, and the mean, spread and permutation "
     "p-value of each count under two nulls of 1,000 random modules - size-matched, and "
     "matched gene-for-gene on how many gene sets a gene belongs to. Drawn as Figure 32a."),
    ("S61 Pathway leading genes", pd.read_csv("results/pathway_leading_genes.csv"),
     "Every risk gene appearing in at least one enriched gene set, with how many sets and "
     "cell types it appears in, how many of those sets are canonical CRC families, which "
     "families, and whether the gene is spatially validated or called in a viability "
     "screen. Drawn as Figure 32c."),
    ("S62 Novel pathways", pd.read_csv("results/pathway_novel_recurrent.csv"),
     "Enriched gene sets outside the canonical CRC families, ranked by how many cell types "
     "carry them, with the best FDR, the largest overlap, the cell types and tissues "
     "involved, the genes, and how many of those genes are spatially validated or called "
     "in a viability screen."),
    ("S63 Pathways pooled", pd.read_csv("results/pathway_pooled_enrichment.csv"),
     "The same over-representation test run once on all 347 risk genes together, ignoring "
     "cell type: gene set, source, canonical family, set size, overlap, p and FDR."),
    ("S64 Pooled pathway tests", pd.read_csv("results/pathway_merged_enrichment.csv"),
     "The same over-representation analysis with the normal and adenoma risk genes of a "
     "cell type pooled into one module: every test reaching at least two risk genes, with "
     "the gene set, its canonical family, the overlap, the hypergeometric p, the "
     "Benjamini-Hochberg q within the pooled module and across the analysis, which "
     "tissue's list contributed the overlapping genes, and the genes. Drawn as Figure 33."),
    ("S65 Pooled pathway nulls", pd.read_csv("results/pathway_merged_null.csv"),
     "Per pooled module: the risk genes it contains and how they split between the normal "
     "list, the adenoma list and both, the number of gene sets it tests, the observed "
     "count of enriched sets, and the mean, spread and permutation p-value of that count "
     "under the size-matched and annotation-matched nulls. Drawn as Figure 34a and 34b."),
    ("S66 Pooled vs split", pd.read_csv("results/pathway_merged_vs_split.csv"),
     "Cell type by cell type: risk genes, enriched gene sets and canonical-family hits in "
     "the split analysis of §2.18 against the pooled analysis of §2.19, with the "
     "annotation-matched null p-value of each and the families reached."),
    ("S67 All 347 pathways", pd.read_csv("results/all347_enrichment.csv"),
     "Over-representation of all 347 risk genes together against the Reactome, "
     "WikiPathways and PID collection: gene set, source, canonical family, aetiology "
     "theme, set size, overlap, hypergeometric p, Benjamini-Hochberg q and the genes. "
     "Drawn as Figure 35a."),
    ("S68 By direction", pd.read_csv("results/all347_by_direction.csv"),
     "The same analysis run separately on the risk-increasing and the protective genes. "
     "Neither half yields an enriched gene set at FDR < 0.1. Drawn as Figure 35b."),
    ("S69 Canonical membership", pd.read_csv("results/canonical_membership.csv"),
     "Every risk gene that is a member of a canonical CRC pathway family, with the "
     "family, the direction of risk, the tissues and cell types carrying the "
     "association, the TWAS p-value, whether the gene is spatially validated and its "
     "class, whether it changed viability in a screen, its druggability tier, and the "
     "hypergeometric test of the family's membership count. Drawn as Figure 35c."),
    ("S70 Gene-pathway proximity", pd.read_csv("results/gene_pathway_proximity.csv"),
     "Shortest path in the protein interaction network from each risk gene to the "
     "nearest member of each canonical CRC family, with a z-score against genes of the "
     "same network degree, and a flag for genes that are themselves family members."),
    ("S71 Aetiology themes", pd.read_csv("results/etiology_themes.csv"),
     "The ten colorectal aetiology themes: the gene sets merged into each, the theme "
     "size, the risk genes found in it against the number expected, the fold, the "
     "hypergeometric p and FDR, the split between risk-increasing and protective genes "
     "with a binomial test of that split, how many are spatially validated or called in "
     "a screen, and the genes. Drawn as Figures 36a and 36c."),
    ("S72 Gene assignment", pd.read_csv("results/etiology_gene_assignment.csv"),
     "Every one of the 347 risk genes with the direction of risk, the tissues and number "
     "of cell types carrying it, its aetiology themes, the canonical families it belongs "
     "to, any family it is network-adjacent to without being a member, whether it is "
     "spatially validated or called in a screen, whether it appears in the pathway "
     "universe at all, and the resulting status. Drawn as Figure 36b."),
    ("S73 Driver genes", pd.read_csv("results/driver_overlap.csv"),
     "Every risk gene with a reported cancer-driver annotation: whether it is on the "
     "colorectal somatic driver list used in §2.4, the confidence and somatic role of the "
     "pan-cancer annotation with its tumour context, the TWAS direction, tissues, cell "
     "types and p-value, spatial validation and class, viability-screen effect, "
     "druggability tier, the canonical pathway families the gene belongs to, and whether "
     "the germline direction agrees with the somatic role. Drawn as Figure 37."),
]

with pd.ExcelWriter(OUT, engine="openpyxl") as xw:
    idx = pd.DataFrame([{"Sheet": n, "Rows": len(d), "Contents": c} for n, d, c in sheets])
    idx.to_excel(xw, sheet_name="Index", index=False)
    for name, df, _ in sheets:
        df.to_excel(xw, sheet_name=name[:31], index=False)

    wb = xw.book
    hdr_fill = PatternFill("solid", fgColor="1F3864")
    hdr_font = Font(color="FFFFFF", bold=True, size=10)
    for ws in wb.worksheets:
        ws.freeze_panes = "A2"
        for c in ws[1]:
            c.fill = hdr_fill; c.font = hdr_font
            c.alignment = Alignment(vertical="center", wrap_text=True)
        ws.row_dimensions[1].height = 28
        for j, col in enumerate(ws.iter_cols(min_row=1, max_row=min(ws.max_row, 200)), start=1):
            w = max((len(str(c.value)) for c in col if c.value is not None), default=10)
            ws.column_dimensions[get_column_letter(j)].width = min(max(w + 2, 10), 46)
        ws.auto_filter.ref = ws.dimensions

print("wrote", OUT)
for n, d, _ in sheets:
    print(f"  {n:24s} {len(d):6d} rows")
