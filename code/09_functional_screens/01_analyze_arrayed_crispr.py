#!/usr/bin/env python3
"""Reproduce Figure 7c from the staged CRISPR workbook and Table 26 gene list."""

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import ttest_1samp


LINES = ("HCT116", "RKO", "SW480")
REPLICATES = ("cell number r1", "cell number r2", "cell number r3")
TABLES = (Path(__file__).resolve().parents[1] / "manuscript_text"
          / "Li_etal_Supplementary_Tables_11062026.xlsx")


def analyze_screen(workbook):
    sheet = pd.read_excel(workbook, sheet_name="crisp-cas9", header=None)
    stats, qc = [], []
    for line in LINES:
        starts = sheet.columns[sheet.iloc[0].eq(line)]
        if len(starts) != 1:
            raise ValueError(f"Expected one CRISPR block for {line}")
        start = starts[0]
        headers = sheet.iloc[1, start + 1:start + 8]
        counts = [headers.index[headers.eq(rep)].item() for rep in REPLICATES]
        block = sheet.iloc[2:, [start, *counts]].dropna(subset=[start])
        genes = block.pop(start).str.strip().replace({"ACAS2": "ACSS2"})
        values = block.apply(pd.to_numeric, errors="raise").set_axis(genes, axis=0)
        values.columns = ["rep1", "rep2", "rep3"]
        nt = values.loc[values.index == "NT"]
        lethal = values.loc[values.index == "LETHAL"]
        samples = values.drop(index=["NT", "LETHAL"])
        if (len(nt) != 8 or len(lethal) != 2 or not samples.index.is_unique
                or not np.isfinite(values.to_numpy()).all()):
            raise ValueError(f"Invalid gene wells or incomplete controls in {line}")

        control = nt.mean()
        separation = control - lethal.mean()
        z_prime = 1 - 3 * (nt.std() + lethal.std()) / separation
        if not ((separation > 0) & (z_prime > 0)).all():
            raise ValueError(f"A {line} plate fails the positive Z-prime check")
        qc.append(pd.DataFrame({"cell_line": line, "replicate": values.columns,
                                "n_samples": len(samples), "n_NT": len(nt),
                                "n_lethal": len(lethal), "NT_mean": control.to_numpy(),
                                "lethal_mean": lethal.mean().to_numpy(),
                                "z_prime": z_prime.to_numpy()}))

        # The archived figure used the staged viability scale for both wells and controls.
        # Keep the within-plate differences on that scale for the paired test.
        diff = samples - control
        test = ttest_1samp(diff.to_numpy(), 0, axis=1)
        stats.append(pd.DataFrame({"cell_line": line, "gene": samples.index,
                                   "n_rep": len(REPLICATES), "value_mean": samples.mean(axis=1),
                                   "control_mean": control.mean(), "delta_paired": diff.mean(axis=1),
                                   "sd_diff": diff.std(axis=1), "t_paired": test.statistic,
                                   "p_paired": test.pvalue, "nominal": test.pvalue < 0.05}))
    return pd.concat(stats, ignore_index=True), pd.concat(qc, ignore_index=True)


def figure7c_source(stats, workbook):
    table = pd.read_excel(workbook, sheet_name="Supplementary Table 26", header=1)
    panel = table.loc[table["Well type"].eq("Risk-increasing scTWAS gene"),
                      ["Gene", "scTWAS Z (CRC risk score)"]].rename(
                          columns={"Gene": "gene", "scTWAS Z (CRC risk score)": "scTWAS_rep_Z"})
    if len(panel) != 24 or not panel.gene.is_unique:
        raise ValueError("Expected 24 unique risk-increasing genes in Supplementary Table 26")
    panel["panel_rank"] = range(1, len(panel) + 1)
    source = panel.merge(stats, on="gene", how="left", validate="one_to_many")
    if len(source) != len(panel) * len(LINES) or source.cell_line.isna().any():
        raise ValueError("A Figure 7c gene is missing a cell-line result")
    source = source.rename(columns={"delta_paired": "delta_viability_points",
                                    "nominal": "significant_P_lt_0_05"})
    source["plate_qc_pass"] = True
    source["mean_delta_qc_passing"] = source.groupby("gene").delta_viability_points.transform("mean")
    return source[["panel_rank", "gene", "cell_line", "delta_viability_points",
                   "significant_P_lt_0_05", "plate_qc_pass", "scTWAS_rep_Z",
                   "mean_delta_qc_passing", "sd_diff", "p_paired"]]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--workbook", type=Path, required=True, help="Staged siRNA_CRISPR-116.xlsx")
    parser.add_argument("--panel-workbook", type=Path, default=TABLES)
    parser.add_argument("--out-dir", type=Path, default=Path("results"))
    args = parser.parse_args()

    stats, qc = analyze_screen(args.workbook)
    source = figure7c_source(stats, args.panel_workbook)
    args.out_dir.mkdir(parents=True, exist_ok=True)
    stats.to_csv(args.out_dir / "crispr_by_cellline.csv", index=False)
    qc.to_csv(args.out_dir / "crispr_plate_qc.csv", index=False)
    source.to_csv(args.out_dir / "Figure7_panelC_sourcedata.csv", index=False)

    calls = source[source.significant_P_lt_0_05]
    print(f"{len(calls)} of {len(source)} tests have nominal P < 0.05 across {calls.gene.nunique()} genes")
    print(calls.groupby("cell_line").size().reindex(LINES, fill_value=0).to_string())


if __name__ == "__main__":
    main()
