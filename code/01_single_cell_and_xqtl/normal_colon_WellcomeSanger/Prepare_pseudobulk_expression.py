#!/usr/bin/env python3
"""
prepare_pseudobulk_expression.py

Input
-----
Gene x sample pseudobulk expression/count matrix:
    gene, sample1, sample2, ...

Example:
    /data/l2_bioinfo1/liq17/scTWAS/data/WellcomeSanger/expression/
    Healthy_TI_then_rectum_278samples_outputs_LauLabel/
    pseudobulk_Epithelial_LauLabel_x_individual/
    pseudobulk_expression_ABS.csv.gz

Processing
----------
1. Read GENCODE GTF with annotation.py.
2. Match Ensembl IDs (version-insensitive) to gene annotation / HUGO-style gene_name.
3. Calculate exon-union gene length from the GTF.
4. Counts -> TPM.
5. If multiple Ensembl IDs map to the same gene symbol, sum their TPM values
   AFTER per-Ensembl length normalization, yielding one row per gene symbol.
6. log2(TPM + 1).
7. Quantile normalization across samples using phenotype.py.
8. Row-wise inverse normal transformation across samples using phenotype.py.

No low-expression filtering is applied in this script.

Example
-------
python prepare_pseudobulk_expression.py \
  --gtf /path/to/gencode.vXX.annotation.gtf.gz
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd


DEFAULT_QTLS_DIR = Path("/data/l2_bioinfo1/liq17/LQ_frameworks/QTLS")
if not DEFAULT_QTLS_DIR.is_dir():
    DEFAULT_QTLS_DIR = Path(__file__).resolve().parents[2] / "lib" / "QTLS"

DEFAULT_INPUT = (
    "/data/l2_bioinfo1/liq17/scTWAS/data/WellcomeSanger/expression/"
    "Healthy_TI_then_rectum_278samples_outputs_LauLabel/"
    "pseudobulk_Epithelial_LauLabel_x_individual/"
    "pseudobulk_expression_GOB.csv.gz" #ABS, CT, EE, GOB, STM, TAC, TUF
)


def parse_args():
    p = argparse.ArgumentParser(
        description=(
            "Annotate pseudobulk Ensembl genes with GENCODE gene names, "
            "then TPM -> log2(TPM+1) -> quantile normalization -> "
            "inverse normal transformation."
        )
    )

    p.add_argument(
        "--input",
        default=DEFAULT_INPUT,
        help="Gene x sample pseudobulk CSV/CSV.GZ. Default: GOB file.",
    )
    p.add_argument(
        "--gtf",
        default="/data/l2_bioinfo1/liq17/LQ_references/gencode/gencode.v26.genes.gtf",
        help="GENCODE GTF or GTF.GZ used for gene name and exon-union length.",
    )
    p.add_argument(
        "--qtls-dir",
        default=DEFAULT_QTLS_DIR,
        help="Directory containing annotation.py and phenotype.py.",
    )
    p.add_argument(
        "--output-dir",
        default=None,
        help="Output directory. Default: <input_dir>/processed_expression",
    )
    p.add_argument(
        "--prefix",
        default=None,
        help="Output prefix. Default inferred from input filename, e.g. ABS.",
    )
    return p.parse_args()


def infer_prefix(path: Path) -> str:
    name = path.name
    for suffix in [".csv.gz", ".csv", ".tsv.gz", ".tsv"]:
        if name.endswith(suffix):
            name = name[: -len(suffix)]
            break
    if name.startswith("pseudobulk_expression_"):
        name = name[len("pseudobulk_expression_") :]
    return name


def write_matrix(
    matrix: pd.DataFrame,
    output_file: Path,
    index_name: str = "gene",
):
    out = matrix.copy()
    out.index.name = index_name
    out.reset_index().to_csv(output_file, index=False, compression="gzip")


def main():
    args = parse_args()

    input_file = Path(args.input)
    gtf_file = Path(args.gtf)
    qtls_dir = Path(args.qtls_dir)

    if not input_file.exists():
        raise FileNotFoundError(f"Input expression file not found: {input_file}")
    if not gtf_file.exists():
        raise FileNotFoundError(f"GTF not found: {gtf_file}")
    if not (qtls_dir / "annotation.py").exists():
        raise FileNotFoundError(f"Missing: {qtls_dir / 'annotation.py'}")
    if not (qtls_dir / "phenotype.py").exists():
        raise FileNotFoundError(f"Missing: {qtls_dir / 'phenotype.py'}")

    # Import the user's QTLS helper modules.
    sys.path.insert(0, str(qtls_dir))

    from annotation import build_gene_annotation, annotate_expression_genes
    from phenotype import (
        quantile_normalize,
        inverse_normal_transform,
        check_quantile_normalization,
        check_inverse_normal_transform,
    )

    prefix = args.prefix or infer_prefix(input_file)

    output_dir = (
        Path(args.output_dir)
        if args.output_dir
        else input_file.parent / "processed_expression"
    )
    output_dir.mkdir(parents=True, exist_ok=True)

    print("=" * 72)
    print("PSEUDOBULK EXPRESSION PREPROCESSING")
    print("=" * 72)
    print(f"Input       : {input_file}")
    print(f"GTF         : {gtf_file}")
    print(f"QTLS module : {qtls_dir}")
    print(f"Output      : {output_dir}")
    print(f"Prefix      : {prefix}")
    print()

    # --------------------------------------------------------------------------
    # STEP 1. Read G x N expression/count matrix
    # --------------------------------------------------------------------------
    print("STEP 1: Read pseudobulk matrix")

    expr = pd.read_csv(input_file)

    if "gene" not in expr.columns:
        raise ValueError(
            f"Expected a 'gene' column; found columns: {list(expr.columns[:10])}"
        )

    sample_cols = [c for c in expr.columns if c != "gene"]

    if len(sample_cols) == 0:
        raise ValueError("No sample columns detected.")

    if expr["gene"].duplicated().any():
        dup = expr.loc[expr["gene"].duplicated(keep=False), "gene"]
        dup_file = output_dir / f"{prefix}.input_duplicate_ENSG.txt"
        dup.drop_duplicates().to_csv(dup_file, index=False, header=False)
        raise ValueError(
            f"Input contains duplicated gene IDs. See: {dup_file}"
        )

    # Enforce numeric expression values.
    expr[sample_cols] = expr[sample_cols].apply(pd.to_numeric, errors="raise")

    x = expr[sample_cols].to_numpy(dtype=float)

    if not np.isfinite(x).all():
        raise ValueError("Expression matrix contains NaN/Inf values.")

    if (x < 0).any():
        raise ValueError("TPM cannot be calculated from negative expression values.")

    integer_like = np.isclose(x, np.round(x), atol=1e-6)
    integer_fraction = float(integer_like.mean())

    print(f"  genes   : {expr.shape[0]:,}")
    print(f"  samples : {len(sample_cols):,}")
    print(f"  integer-like values: {integer_fraction:.4%}")

    if integer_fraction < 0.999:
        print(
            "  WARNING: expression values are not almost entirely integer-like. "
            "TPM is normally calculated from raw/pseudobulk counts; verify that "
            "the source matrix represents count-scale expression."
        )

    # --------------------------------------------------------------------------
    # STEP 2. Build GTF annotation + exon-union gene lengths
    # --------------------------------------------------------------------------
    print("\nSTEP 2: Build gene annotation from GTF")

    gene_annot = build_gene_annotation(
        str(gtf_file),
        autosomes_only=True,
    )

    annot_file = output_dir / f"{prefix}.GENCODE_gene_annotation.csv.gz"
    gene_annot.to_csv(annot_file, index=False, compression="gzip")

    print(f"  GTF autosomal genes with length: {gene_annot.shape[0]:,}")

    # --------------------------------------------------------------------------
    # STEP 3. Match ENSG -> gene_name and gene length
    # --------------------------------------------------------------------------
    print("\nSTEP 3: Match expression genes to GENCODE")

    annotated = annotate_expression_genes(
        expr,
        gene_annot,
    )

    # annotate_expression_genes preserves original sample column names.
    annotated_sample_cols = [
        c for c in sample_cols if c in annotated.columns
    ]

    if annotated_sample_cols != sample_cols:
        raise RuntimeError("Sample columns changed unexpectedly during annotation.")

    annotated["geneid_nover"] = (
        annotated["geneid"]
        .astype(str)
        .str.replace(r"\.\d+$", "", regex=True)
    )

    # Explicitly recover matched/unmatched input genes.
    matched_ensg = set(annotated["gene"].astype(str))
    unmatched = expr.loc[
        ~expr["gene"].astype(str).isin(matched_ensg),
        ["gene"],
    ].copy()

    unmatched_file = output_dir / f"{prefix}.unmatched_genes.csv.gz"
    unmatched.to_csv(unmatched_file, index=False, compression="gzip")

    mapping_cols = [
        "gene",
        "geneid",
        "geneid_nover",
        "genename",
        "gene_type",
        "CHR",
        "start",
        "end",
        "direction",
        "Length",
    ]

    mapping_file = output_dir / f"{prefix}.ENSG_to_gene_name_mapping.csv.gz"
    annotated[mapping_cols].to_csv(
        mapping_file,
        index=False,
        compression="gzip",
    )

    print(f"  matched genes   : {annotated.shape[0]:,}")
    print(f"  unmatched genes : {unmatched.shape[0]:,}")

    # --------------------------------------------------------------------------
    # STEP 4. Counts -> TPM using exon-union gene length
    #
    # RPK = counts / length_kb
    # TPM = RPK / sum(RPK per sample) * 1e6
    # --------------------------------------------------------------------------
    print("\nSTEP 4: TPM normalization")

    lengths_bp = annotated["Length"].to_numpy(dtype=float)

    if np.any(lengths_bp <= 0):
        raise ValueError("Zero/negative gene lengths detected after annotation.")

    counts = annotated[sample_cols].astype(float)

    rpk = counts.div(lengths_bp / 1000.0, axis=0)
    rpk_sum = rpk.sum(axis=0)

    if (rpk_sum <= 0).any():
        bad = list(rpk_sum.index[rpk_sum <= 0])
        raise ValueError(
            "Samples with zero total RPK detected: " + ", ".join(bad[:10])
        )

    tpm_ensg = rpk.div(rpk_sum, axis=1) * 1_000_000.0

    # Use gene symbol as requested.
    # If multiple Ensembl genes map to the same symbol, TPM is summed AFTER
    # each Ensembl gene has been length-normalized.
    tpm_ensg.index = annotated["genename"].astype(str).values

    symbol_counts = pd.Series(
        tpm_ensg.index,
        name="gene",
    ).value_counts()

    duplicate_symbols = symbol_counts[symbol_counts > 1]

    duplicate_symbol_file = output_dir / f"{prefix}.duplicate_gene_symbols.csv.gz"
    if len(duplicate_symbols):
        duplicate_map = annotated[
            annotated["genename"].astype(str).isin(duplicate_symbols.index)
        ][mapping_cols].copy()
        duplicate_map.to_csv(
            duplicate_symbol_file,
            index=False,
            compression="gzip",
        )
    else:
        pd.DataFrame(
            columns=mapping_cols
        ).to_csv(
            duplicate_symbol_file,
            index=False,
            compression="gzip",
        )

    tpm = tpm_ensg.groupby(level=0, sort=False).sum()

    # Sanity check: TPM sums should be ~1e6 per sample.
    tpm_sums = tpm.sum(axis=0)
    max_tpm_sum_error = float(np.max(np.abs(tpm_sums - 1_000_000.0)))

    print(f"  Ensembl genes entering TPM : {tpm_ensg.shape[0]:,}")
    print(f"  unique gene symbols        : {tpm.shape[0]:,}")
    print(f"  duplicated symbols merged  : {len(duplicate_symbols):,}")
    print(f"  max |sum(TPM)-1e6|         : {max_tpm_sum_error:.6f}")

    tpm_file = output_dir / f"{prefix}.TPM.csv.gz"
    write_matrix(tpm, tpm_file)

    # --------------------------------------------------------------------------
    # STEP 5. log2(TPM + 1)
    # --------------------------------------------------------------------------
    print("\nSTEP 5: log2(TPM + 1)")

    log2_tpm = np.log2(tpm + 1.0)

    log2_file = output_dir / f"{prefix}.TPM.log2.csv.gz"
    write_matrix(log2_tpm, log2_file)

    # --------------------------------------------------------------------------
    # STEP 6. Quantile normalization
    # phenotype.quantile_normalize expects rows=features, columns=samples.
    # --------------------------------------------------------------------------
    print("\nSTEP 6: Quantile normalization")

    qn = quantile_normalize(log2_tpm)

    check_quantile_normalization(qn)

    qn_file = output_dir / f"{prefix}.TPM.log2.QN.csv.gz"
    write_matrix(qn, qn_file)

    # --------------------------------------------------------------------------
    # STEP 7. Row-wise inverse normal transformation across samples
    # --------------------------------------------------------------------------
    print("\nSTEP 7: Inverse normal transformation")

    inv = inverse_normal_transform(qn)

    inv_check = check_inverse_normal_transform(inv)
    inv_check.index = inv.index
    inv_check.index.name = "gene"

    inv_file = output_dir / f"{prefix}.TPM.log2.QN.INV.csv.gz"
    write_matrix(inv, inv_file)

    inv_check_file = output_dir / f"{prefix}.TPM.log2.QN.INV.check.csv.gz"
    inv_check.reset_index().to_csv(
        inv_check_file,
        index=False,
        compression="gzip",
    )

    # --------------------------------------------------------------------------
    # STEP 8. Processing summary
    # --------------------------------------------------------------------------
    summary_file = output_dir / f"{prefix}.preprocessing_summary.txt"

    with open(summary_file, "w") as fh:
        fh.write(f"input={input_file}\n")
        fh.write(f"gtf={gtf_file}\n")
        fh.write(f"qtls_dir={qtls_dir}\n")
        fh.write(f"n_input_genes={expr.shape[0]}\n")
        fh.write(f"n_samples={len(sample_cols)}\n")
        fh.write(f"integer_like_fraction={integer_fraction:.8f}\n")
        fh.write(f"n_annotated_ensembl_genes={annotated.shape[0]}\n")
        fh.write(f"n_unmatched_genes={unmatched.shape[0]}\n")
        fh.write(f"n_duplicate_gene_symbols={len(duplicate_symbols)}\n")
        fh.write(f"n_final_gene_symbols={tpm.shape[0]}\n")
        fh.write(f"max_abs_TPM_sum_error={max_tpm_sum_error:.8f}\n")
        fh.write("processing=TPM->log2(TPM+1)->QN->rowwise_INV\n")
        fh.write(f"final={inv_file}\n")

    print("\n" + "=" * 72)
    print("DONE")
    print("=" * 72)
    print(f"Mapping : {mapping_file}")
    print(f"TPM     : {tpm_file}")
    print(f"log2    : {log2_file}")
    print(f"QN      : {qn_file}")
    print(f"INV     : {inv_file}")
    print(f"Summary : {summary_file}")


if __name__ == "__main__":
    main()
