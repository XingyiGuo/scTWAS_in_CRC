# CRISPR-Cas9 functional screens

## Arrayed viability screen

Reproduce the arrayed viability results for Figure 7c from the staged
`siRNA_CRISPR-116.xlsx` workbook. The analysis uses its `crisp-cas9` sheet.
The [supplementary workbook](../manuscript_text/Li_etal_Supplementary_Tables_11062026.xlsx)
supplies the 24 risk-increasing genes, their order, and the displayed scTWAS
Z scores from Table 26.

## Run the analysis

Use the [Python dependencies](../../env/python_requirements.txt) and run
[01_analyze_arrayed_crispr.py](01_analyze_arrayed_crispr.py) with
`--workbook` pointing to the staged input. For example, from the repo root
on the project workstation:

```powershell
python code/09_functional_screens/01_analyze_arrayed_crispr.py `
  --workbook 'L:/AI-gene-drug/siRNA_CRISPR-116.xlsx' `
  --out-dir results
```

| Output | Contents |
|---|---|
| `crispr_plate_qc.csv` | Control counts, means, and Z-prime for the nine replicate plates |
| `crispr_by_cellline.csv` | Effect sizes and paired tests for all 116 targets in each of the three cell lines |
| `Figure7_panelC_sourcedata.csv` | The 72 Figure 7c effects and nominal calls, with paired P values and standard deviations |

The archived Figure 7c analysis subtracts each plate's mean of eight
non-targeting wells from the target well on the workbook's staged viability
scale. A two-sided paired t-test across the three within-plate differences
uses nominal P < 0.05. All nine plates pass the positive Z-prime check.
The archived source in the project `Experiments/Archive` folder agrees on
all 72 effects and calls. The 18 significant decreases involve 13 genes
and include five in HCT116, two in RKO, and 11 in SW480.

## Manuscript and figure handoffs

Table 26 contains raw nuclei counts and percentages of dead cells. Dividing
those nuclei counts by the same-plate non-targeting mean and testing against
100%, as described in the dedicated manuscript Methods paragraph, gives 17
significant decreases across the same 13 genes. The HLA-E call in HCT116
accounts for the difference. The manuscript's general statistics paragraph
also specifies BH correction, whereas Figure 7c uses nominal P values.

## DepMap dependency analysis

Run [02_analyze_depmap_dependency.py](02_analyze_depmap_dependency.py) with
`--experiments 'W:/VUMC_projects/scTWAS/Experiments' --out-dir results`.
The `--experiments` directory supplies these inputs.

| Input | Role |
|---|---|
| `depmap_intestine_all_geneset_genes_percellline.csv` | Per-cell-line Chronos scores for the Figure 7a comparisons |
| `DepMap_essentiality_by_geneset_sourcedata.csv` | Original gene-set membership and display order. The script recalculates the statistics from per-cell-line scores. |
| `scTWAS_347genes_DepMap_largeintestine_percellline.csv` | Per-cell-line Chronos scores for the Figure 7b dependency calls |

The script recalculates the four gene-set distributions, Kruskal–Wallis
and Holm-adjusted pairwise tests, and dependency classifications. The archived
run retained 274 scTWAS genes with scores in more than 30 lines. Of these,
36 had a median score below −0.5 and 24 were below −0.5 in at least 90%
of evaluable lines.

[Section 11](../11_figures/figure7_depmap_screens/build_figure7.py) contains
the Figure 7 renderer. Panels a and b use the DepMap analysis above and
Supplementary Table 25 contains the deposited gene-effect scores. The
script starts from exported scores and the original comparison groups.
