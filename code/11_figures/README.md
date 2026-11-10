# Figure reproduction

Use the current [figure captions](../manuscript_text/Li_etal_Figures_SupFigures_11062026.docx)
for numbering. The requested supplementary figures are labelled **Extended Data**
in that document. Older filenames use different numbers.

## Source data

[SourceData.xlsx](SourceData.xlsx) contains the source values for Figures 2–7
and Extended Data Figures 1 and 4–9. Its **Index** links to each data sheet
and records the archive source and expected CSV path where applicable.
Use the repository's
[Python](../../env/python_requirements.txt) and [R](../../env/r_packages.txt)
dependencies.

| Figure | Renderer and input |
|---|---|
| 2, 3 | [Manhattan scripts](manhattan) from the archived plotted-gene tables with significance flags and Bonferroni cutoffs. |
| 4a | [Model overlap triangle](orthogonal_support/multiomics_de_support/01_model_overlap_triangle.R) from the 546 gene–cell-type associations. |
| 4b–i | [TF networks and enrichment](network_figures/published_panels) from the correlation and lead-variant occupancy tables. |
| 5a–c | [Grouped bars](orthogonal_support/multiomics_de_support/04_de_support_grouped_bars.R), [ASC/SSC linkage](orthogonal_support/multiomics_de_support/02_asc_ssc_linkage.R), and [support bubbles](orthogonal_support/multiomics_de_support/05_de_support_bubble_heatmap.R). Panel c uses the WTA GeoMx results. |
| 5d | [CosMx composition map](orthogonal_support/spatial_geomx_cosmx/plot_figure5d.py) from the archived coordinates and cell groups. |
| 5e–j | [CosMx gene maps and paired plots](orthogonal_support/spatial_geomx_cosmx/redraw_figure5ej.py). The workbook has six per-cell sheets, a paired-patient table and a gene summary. |
| 6 | [Radial tree](druggability) from the drug–target records. |
| 7 | [DepMap and arrayed CRISPR panels](figure7_depmap_screens/make_figure7.py) from the archived analysis inputs. [Section 09](../09_functional_screens/README.md) documents the upstream analyses. |
| Extended Data 1 and 4 | [UMAPs, peak counts and model R²](model_performance/plot_extended_data_1_4.py). The R² source includes both normal-colon cohorts and the polyp cohort. |
| Extended Data 5 and 6 | [Multi-ancestry Manhattan scripts](manhattan) from the archived association tables with multi-ancestry-only highlights. |
| Extended Data 7 | [TF evidence heatmaps](network_figures/published_panels/mkfigS1.py) using the Figure 4b–i inputs. |
| Extended Data 8 | [TF occupancy and co-expression comparison](network_figures/published_panels/plot_extended_data_8.py) using the same inputs. The workbook identifies which panel values were reconstructed. |
| Extended Data 9 | [iLINCS compound plot](druggability/plot_ilincs_compound_panel.py) from Supplementary Table 24. |

## Render the figures

Run from the repository root. Figures 5d–j and 6 can read the workbook
directly. For the other scripts, export the required worksheets as **CSV UTF-8**
to the paths listed in the workbook Index under your working directory. Keep
the worksheet headers when exporting. The scripts for Figures 4 and 5 render
individual panels. The final page layouts are archived as
`Figure4_TF_Genes_Pathways.ai` and `Figure5_MultiOmics_Support.ai` under
`W:/VUMC_projects/scTWAS/ms-guo/Figures`.

```powershell
$code = (Resolve-Path code/11_figures).Path
$book = Join-Path $code 'SourceData.xlsx'
$stage = (New-Item -ItemType Directory -Force scratch/figure_reproduction).FullName
$out = Join-Path $stage 'rendered'
```

**Figures 2, 3 and Extended Data 5, 6.** Set `SCTWAS_MANHATTAN_DATA` to the
corresponding CSV under `$stage/manhattan`. This route uses the saved plot
coordinates and flags, and skips raw association and EnsDb processing.

| CSV | R script under `$code/manhattan` |
|---|---|
| `Figure2.csv` | `Manhattan_NL_EUR.R` |
| `Figure3.csv` | `Manhattan_ADSER_EUR.R` |
| `ExtendedDataFigure5.csv` | `Manhattan_NL_Multi_EURnotSig_MultiSig_only.R` |
| `ExtendedDataFigure6.csv` | `Manhattan_ADSER_Multi_EURnotSig_MultiSig_only.R` |

```powershell
$env:SCTWAS_MANHATTAN_DATA = "$stage/manhattan/Figure2.csv"
$env:SCTWAS_FIGURE_OUT = $out
Rscript "$code/manhattan/Manhattan_NL_EUR.R"
```

**Figures 4a and 5a–c.** The plots and source-data exports go under
`$stage/multiomics/figure`.

```powershell
$env:SC_TWAS_SUPPORT_ROOT = "$stage/multiomics"
Rscript "$code/orthogonal_support/multiomics_de_support/01_model_overlap_triangle.R"
Rscript "$code/orthogonal_support/multiomics_de_support/04_de_support_grouped_bars.R"
Rscript "$code/orthogonal_support/multiomics_de_support/02_asc_ssc_linkage.R"
Rscript "$code/orthogonal_support/multiomics_de_support/05_de_support_bubble_heatmap.R"
```

**Figure 4b–i and Extended Data 7–8.** Run from the directory containing the TF CSVs.
The script names retain the archive numbering and their output names follow
the current manuscript.

```powershell
Push-Location "$stage/tf"
python "$code/network_figures/published_panels/mkfig1.py"   # Figure 4b,c
python "$code/network_figures/published_panels/mkfig2.py"   # Figure 4d,e
python "$code/network_figures/published_panels/mkfig3.py"   # Figure 4f–i
python "$code/network_figures/published_panels/mkfigS1.py"  # Extended Data 7
python "$code/network_figures/published_panels/plot_extended_data_8.py" --out-dir $out
Pop-Location
```

**Figure 5d–j.**

```powershell
python "$code/orthogonal_support/spatial_geomx_cosmx/plot_figure5d.py" `
  --input $book --out-dir $out
python "$code/orthogonal_support/spatial_geomx_cosmx/redraw_figure5ej.py" `
  --workbook $book --out-dir $out
```

**Figure 6.** Step 01 builds the tree from the workbook. The Node.js renderer
draws the SVG and step 02 exports PDF and PNG.

```powershell
$env:SCTWAS_DRUGS_ROOT = "$stage/drugs"
$env:SCTWAS_FIGURE_WORKBOOK = $book
npm install --prefix "$stage/drugs" d3@7 jsdom@26
python "$code/druggability/01_prepare_radial_tree.py"
node "$code/druggability/05_render_radial_tree.mjs"
python "$code/druggability/02_export_radial_tree.py"
```

**Figure 7 and Extended Data 1, 4, 9.**

```powershell
python "$code/figure7_depmap_screens/make_figure7.py" --data-dir "$stage/figure7" --out-dir $out
python "$code/model_performance/plot_extended_data_1_4.py" --data-dir "$stage/models" --out-dir $out
python "$code/druggability/plot_ilincs_compound_panel.py" `
  --input "$stage/ilincs/ExtendedDataFigure9_Table24.csv" --out-dir $out
```
