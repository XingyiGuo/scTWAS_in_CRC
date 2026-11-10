#!/usr/bin/env python3

import argparse
import os

import numpy as np
import pandas as pd

DEFAULT_ASSOCIATION_DIR = (
    "/data/l2_bioinfo1/liq17/scTWAS/"
    "TWAS_HTAN_GSE266546_WellcomeSanger_res/"
    "SPrediXcan_peakSNP_association"
)

CELL_TYPES = ["ABS", "CT", "EE", "GOB", "STM", "TAC", "TUF"]


def bool_to_label(x):
    return "Yes" if bool(x) else "No"


def process_cell_type(base_dir, cohort, ct, summary_row):
    in_file = os.path.join(
        base_dir,
        f"HTAN_NL_WellcomeSanger_NL_{ct}_peakSNP.{cohort}.CauchyMeta.csv"
    )

    if not os.path.exists(in_file):
        raise FileNotFoundError(f"Missing meta file: {in_file}")

    df = pd.read_csv(in_file)

    meta_col = "META_P" if "META_P" in df.columns else f"META_{cohort}_P"
    required_cols = [
        "gene",
        meta_col,
        "pvalue_HTAN_NL",
        "pvalue_WellcomeSanger_NL",
    ]
    missing = [x for x in required_cols if x not in df.columns]
    if missing:
        raise ValueError(
            f"{ct}: missing required columns: {', '.join(missing)}\n"
            f"Available columns: {', '.join(df.columns)}"
        )

    # ------------------------------------------------------------------
    # Denominators
    # ------------------------------------------------------------------
    n_meta = int(summary_row["Genes_in_merged_output"])
    n_htan = int(summary_row["HTAN_NL_kept"])
    n_wellcome = int(summary_row["WellcomeSanger_NL_kept"])

    if len(df) != n_meta:
        print(
            f"WARNING [{ct}]: rows in meta file ({len(df)}) != "
            f"Genes_in_merged_output ({n_meta})"
        )

    # Bonferroni thresholds
    meta_threshold = 0.05 / n_meta if n_meta else np.nan
    htan_threshold = 0.05 / n_htan if n_htan else np.nan
    wellcome_threshold = 0.05 / n_wellcome if n_wellcome else np.nan

    # ------------------------------------------------------------------
    # Per-gene significance annotation
    # ------------------------------------------------------------------
    df["META_Bonferroni_threshold"] = meta_threshold
    df["VUMC_HTAN_GSE266546_Bonferroni_threshold"] = htan_threshold
    df["WellcomeSanger_Bonferroni_threshold"] = wellcome_threshold

    df["meta_sig"] = (
        df[meta_col].notna()
        & (df[meta_col] < meta_threshold)
    )

    df["VUMC_HTAN_GSE266546_sig"] = (
        df["pvalue_HTAN_NL"].notna()
        & (df["pvalue_HTAN_NL"] < htan_threshold)
    )

    df["WellcomeSanger_sig"] = (
        df["pvalue_WellcomeSanger_NL"].notna()
        & (df["pvalue_WellcomeSanger_NL"] < wellcome_threshold)
    )

    # Requested interpretation:
    # for genes that are NOT meta-significant, check whether either
    # individual dataset is significant using its own kept-gene denominator.
    conditions = [
        df["meta_sig"],
        (~df["meta_sig"])
        & df["VUMC_HTAN_GSE266546_sig"]
        & df["WellcomeSanger_sig"],
        (~df["meta_sig"])
        & df["VUMC_HTAN_GSE266546_sig"]
        & (~df["WellcomeSanger_sig"]),
        (~df["meta_sig"])
        & (~df["VUMC_HTAN_GSE266546_sig"])
        & df["WellcomeSanger_sig"],
    ]

    choices = [
        "Meta_sig",
        "Meta_not_sig_Both_individual_sig",
        "Meta_not_sig_VUMC_HTAN_GSE266546_sig",
        "Meta_not_sig_WellcomeSanger_sig",
    ]

    df["significance_category"] = np.select(
        conditions,
        choices,
        default="Not_significant"
    )

    # User-friendly Yes/No columns
    df["meta_sig"] = df["meta_sig"].map(bool_to_label)
    df["VUMC_HTAN_GSE266546_sig"] = (
        df["VUMC_HTAN_GSE266546_sig"].map(bool_to_label)
    )
    df["WellcomeSanger_sig"] = df["WellcomeSanger_sig"].map(bool_to_label)

    # ------------------------------------------------------------------
    # Write annotated per-gene file
    # ------------------------------------------------------------------
    out_file = os.path.join(
        base_dir,
        f"HTAN_NL_WellcomeSanger_NL_{ct}_peakSNP.{cohort}."
        f"CauchyMeta.significance_annotated.csv"
    )
    df.to_csv(out_file, index=False)

    # ------------------------------------------------------------------
    # Summary counts
    # ------------------------------------------------------------------
    category_counts = df["significance_category"].value_counts()

    n_meta_sig = int(category_counts.get("Meta_sig", 0))
    n_both_individual = int(
        category_counts.get("Meta_not_sig_Both_individual_sig", 0)
    )
    n_htan_only = int(
        category_counts.get("Meta_not_sig_VUMC_HTAN_GSE266546_sig", 0)
    )
    n_wellcome_only = int(
        category_counts.get("Meta_not_sig_WellcomeSanger_sig", 0)
    )
    n_not_sig = int(category_counts.get("Not_significant", 0))

    result = {
        "Cell_type": ct,

        "Genes_in_merged_output": n_meta,
        "VUMC_HTAN_GSE266546_genes_kept": n_htan,
        "WellcomeSanger_genes_kept": n_wellcome,

        "Meta_Bonferroni_threshold": meta_threshold,
        "VUMC_HTAN_GSE266546_Bonferroni_threshold": htan_threshold,
        "WellcomeSanger_Bonferroni_threshold": wellcome_threshold,

        "Meta_sig_n": n_meta_sig,

        "Meta_not_sig_Both_individual_sig_n": n_both_individual,
        "Meta_not_sig_VUMC_HTAN_GSE266546_sig_n": n_htan_only,
        "Meta_not_sig_WellcomeSanger_sig_n": n_wellcome_only,

        "Not_significant_n": n_not_sig,

        # Useful aggregate:
        "Meta_not_sig_but_any_individual_sig_n":
            n_both_individual + n_htan_only + n_wellcome_only,
    }

    print("=" * 70)
    print(f"Cell type: {ct}")
    print(f"Genes_in_merged_output                : {n_meta}")
    print(f"Meta threshold                        : {meta_threshold:.6g}")
    print(f"VUMC kept / threshold                 : {n_htan} / {htan_threshold:.6g}")
    print(f"Wellcome kept / threshold             : {n_wellcome} / {wellcome_threshold:.6g}")
    print(f"Meta significant                      : {n_meta_sig}")
    print(f"Meta nonsig, both individual sig      : {n_both_individual}")
    print(f"Meta nonsig, VUMC only sig            : {n_htan_only}")
    print(f"Meta nonsig, WellcomeSanger only sig  : {n_wellcome_only}")
    print(f"Not significant                       : {n_not_sig}")
    print(f"Wrote: {out_file}")

    return result


def main():
    parser = argparse.ArgumentParser(description="Annotate normal-colon ACAT significance.")
    parser.add_argument("cohort", choices=("EUR", "Multi"), nargs="?", default="EUR")
    args = parser.parse_args()
    association_dir = os.environ.get("ASSO_OUTPUT_DIR", DEFAULT_ASSOCIATION_DIR)
    base_dir = os.path.join(association_dir, f"Cauchy_Meta_NL_WellcomeSanger_NL_{args.cohort}")
    summary_file = os.path.join(
        base_dir,
        f"HTAN_NL_WellcomeSanger_NL_{args.cohort}_CauchyMeta_summary_7celltypes_R2global_gt0.01.csv",
    )
    if args.cohort == "EUR" and not os.path.exists(summary_file):
        legacy_dir = os.path.join(association_dir, "Cauchy_Meta_NL_WellcomeSanger_NL")
        legacy_summary = os.path.join(
            legacy_dir, "HTAN_NL_WellcomeSanger_NL_CauchyMeta_summary_7celltypes_R2global_gt0.01.csv"
        )
        if os.path.exists(legacy_summary):
            base_dir, summary_file = legacy_dir, legacy_summary

    if not os.path.exists(summary_file):
        raise FileNotFoundError(f"Missing summary file: {summary_file}")

    summary = pd.read_csv(summary_file)

    required_summary_cols = [
        "Cell_type",
        "Genes_in_merged_output",
        "HTAN_NL_kept",
        "WellcomeSanger_NL_kept",
    ]
    missing = [x for x in required_summary_cols if x not in summary.columns]
    if missing:
        raise ValueError(
            "Summary file missing required columns: "
            + ", ".join(missing)
        )

    summary["Cell_type"] = summary["Cell_type"].astype(str).str.upper()
    summary = summary.set_index("Cell_type")

    results = []

    for ct in CELL_TYPES:
        if ct not in summary.index:
            raise ValueError(f"Cell type {ct} not found in summary table.")

        row = summary.loc[ct]
        if isinstance(row, pd.DataFrame):
            raise ValueError(
                f"More than one summary row found for cell type {ct}."
            )

        results.append(process_cell_type(base_dir, args.cohort, ct, row))

    out_summary = os.path.join(
        base_dir, "HTAN_NL_WellcomeSanger_NL_significance_summary_7celltypes.csv"
    )
    out_summary_df = pd.DataFrame(results)
    out_summary_df.to_csv(out_summary, index=False)

    print("\n" + "=" * 70)
    print("7-cell-type significance summary")
    print("=" * 70)
    print(out_summary_df.to_string(index=False))
    print(f"\nSummary written to: {out_summary}")


if __name__ == "__main__":
    main()
