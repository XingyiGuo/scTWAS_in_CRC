# Viz - figure scripts

Figure and panel scripts kept together under their original working-directory name. They
are the canonical copies; the numbered sections do not duplicate them.

| File | Figure it draws |
|---|---|
| `Manhattan_NL_EUR.R`, `Manhattan_NL_Multi_EURnotSig_MultiSig_only.R` | Stacked Manhattan panels, normal colorectal mucosa (7 cell types), European and multi-ancestry |
| `Manhattan_ADSER_EUR.R`, `Manhattan_ADSER_Multi_EURnotSig_MultiSig_only.R` | The same for the polyp (adenoma/serrated) arm, 9 cell types |
| `Violin_ADSER_NL_EUR_R2gt0.01_with_peak_barplots.R` | Model R-squared violins with the accessible-chromatin peak-count barplots |
| `CRC51TF_JASPAR2024_motif_logos.R` | JASPAR2024 motif logos for the 51 CRC susceptibility-associated TFs |
| `Plot_Gene_SNP_TF_Cistrome_ggplot2_updated_36x14_TSS_TFcategory_topRSlabels.R` | Gene-SNP-TF Cistrome panel, ordered by TSS and TF category |
| `Final_Viz_4panel.R` | Combined four-panel figure |
| `08_prep_cosmx_highlight_8panel_linux.py`, `08_run_cosmx_8panel_linux.sh`, `09_fig_NOR_CORE_highlight_8panel_linux.R` | CosMx highlight panels, 8-panel version. The 4-panel siblings are in `07_orthogonal_support/multiomics_de_support/`. |

The knockout-validation DE figures are in `09_functional_screens/experimental_validation/`.

Related figure code that stays with its analysis: `11_figures/` (workflow schematic,
annotation UMAPs, Figure 7 layout, Manhattan source-data extraction) and
`06_tf_network/network_figures/` (the network diagrams).
