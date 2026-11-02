#!/usr/bin/env python3
"""
08_prepare_pseudobulk_expression.py

Prepare sample-level cell-type pseudobulk count matrices for PEER / QTL analysis.

This version reuses the shared QTLS framework:
    annotation.py
        - get_gene_annotation()

    phenotype.py
        - quantile_normalize()
        - inverse_normal_transform()
        - check_quantile_normalization()
        - check_inverse_normal_transform()

Input
-----
pseudobulk_counts_<CT>.csv

Rows    = genes
Columns = samples
First column = gene identifier / gene symbol
Values  = raw summed pseudobulk counts

Workflow
--------
raw pseudobulk counts
    -> match GENCODE v26 annotation
    -> autosomal genes (chr1-22)
    -> calculate exon-union gene length
    -> raw TPM
    -> expression QC:
         TPM >= 0.1 in >=20% samples
         AND
         raw count >=6 in >=20% samples
    -> log2(TPM + 1)
    -> quantile normalization
    -> inverse-normal transform per gene
    -> PEER-ready gene x sample matrix

Final output
------------
pseudobulk_counts_<CT>.tpmlog2.qc.qn.rmlowcounts.inv.csv

NOTE
----
No normalization by cell number is performed.
"""

from __future__ import annotations

import argparse
import importlib
import math
import os
import re
import shlex
import sys
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pandas as pd


# ==============================================================================
# Shared QTLS framework
# ==============================================================================

QTLS_FRAMEWORK_CANDIDATES = [
    "/v5000/data/l2_bioinfo1/liq17/LQ_frameworks/QTLS",
    "/data/l2_bioinfo1/liq17/LQ_frameworks/QTLS",
]

QTLS_FRAMEWORK = None

for framework_path in QTLS_FRAMEWORK_CANDIDATES:
    if os.path.isdir(framework_path):
        QTLS_FRAMEWORK = framework_path
        break

if QTLS_FRAMEWORK is None:
    raise FileNotFoundError(
        "Could not locate the QTLS framework. Checked:\n  "
        + "\n  ".join(QTLS_FRAMEWORK_CANDIDATES)
    )

if QTLS_FRAMEWORK not in sys.path:
    sys.path.insert(0, QTLS_FRAMEWORK)

import annotation
import phenotype

print(f"Using QTLS framework: {QTLS_FRAMEWORK}")
print(f"annotation.py: {annotation.__file__}")
print(f"phenotype.py:  {phenotype.__file__}")


# ==============================================================================
# Helpers specific to pseudobulk TPM
# ==============================================================================

def _parse_bash_scalar(value: str):
    """Parse a simple bash KEY=VALUE scalar used by the tensorQTL config."""
    value = value.strip()
    if not value:
        return ""
    if len(value) >= 2 and value[0] == value[-1] and value[0] in {"'", '"'}:
        return value[1:-1]
    if re.fullmatch(r"[-+]?\d+", value):
        return int(value)
    if re.fullmatch(r"[-+]?(?:\d+\.\d*|\.\d+)", value):
        return float(value)
    return value


def _parse_bash_config(config_path: str):
    """Read the simple bash-style tensorQTL config used for VUMC_HTAN ADSER."""
    path = Path(config_path)
    if not path.is_file():
        raise FileNotFoundError(f"Config file not found: {path}")

    raw = path.read_text()
    values = {}

    # Parse bash arrays such as CELL_TYPES=("ASC" "SSC" ...).
    array_pattern = re.compile(r'(?ms)^\s*([A-Za-z_][A-Za-z0-9_]*)\s*=\s*\((.*?)\)\s*$')
    spans = []
    for match in array_pattern.finditer(raw):
        key = match.group(1)
        values[key] = shlex.split(match.group(2), comments=True, posix=True)
        spans.append(match.span())

    # Blank array blocks before scalar parsing.
    chars = list(raw)
    for start, end in spans:
        for i in range(start, end):
            if chars[i] != '\n':
                chars[i] = ' '
    scalar_text = ''.join(chars)

    for line in scalar_text.splitlines():
        line = line.strip()
        if not line or line.startswith('#') or '=' not in line:
            continue
        key, value = line.split('=', 1)
        key = key.strip()
        if not re.fullmatch(r'[A-Za-z_][A-Za-z0-9_]*', key):
            continue
        try:
            pieces = shlex.split(value, comments=True, posix=True)
            parsed = '' if not pieces else (_parse_bash_scalar(pieces[0]) if len(pieces) == 1 else ' '.join(pieces))
        except ValueError:
            parsed = _parse_bash_scalar(value)
        values[key] = parsed

    cfg = SimpleNamespace(**values)
    if not hasattr(cfg, 'WORK_DIR'):
        raise ValueError(f"{path} must define WORK_DIR.")

    dataset_prefix = getattr(cfg, 'DATASET_PREFIX', 'VUMC_HTAN')
    tissue = getattr(cfg, 'TISSUE', 'ADSER')

    # Derive Step-2 expression paths from the eQTL config.
    if not hasattr(cfg, 'EXPRESSION_PREP_DIR'):
        cfg.EXPRESSION_PREP_DIR = os.path.join(
            str(cfg.WORK_DIR), 'data', str(dataset_prefix), 'expression',
            f'{tissue}_pseudobulk_expression_for_peer'
        )
    if not hasattr(cfg, 'PSEUDOBULK_INPUT_PATTERN'):
        cfg.PSEUDOBULK_INPUT_PATTERN = os.path.join(
            cfg.EXPRESSION_PREP_DIR, '{ct}', 'pseudobulk_counts_{ct}.csv'
        )

    # tensorQTL config name -> preprocessing-script name.
    if not hasattr(cfg, 'GENE_ANNOT_FILE'):
        if hasattr(cfg, 'GENE_ANNOTATION'):
            cfg.GENE_ANNOT_FILE = cfg.GENE_ANNOTATION
        else:
            raise ValueError(f"{path} must define GENE_ANNOTATION or GENE_ANNOT_FILE.")

    # Step-2 defaults used by the VUMC_HTAN preprocessing pipeline.
    defaults = {
        'AUTOSOMES_ONLY': True,
        'TPM_THRESHOLD': 0.1,
        'COUNT_THRESHOLD': 6,
        'SAMPLE_FRAC_THRESHOLD': 0.20,
        'LOG2_PSEUDOCOUNT': 1.0,
        'MIN_SAMPLE_TOTAL_COUNTS': 0,
        'SAMPLE_KEEP_FILE': None,
        'SAMPLE_KEEP_COLUMN': 'Sample',
        'SAVE_INTERMEDIATES': True,
    }
    for key, value in defaults.items():
        if not hasattr(cfg, key):
            setattr(cfg, key, value)

    if not hasattr(cfg, 'CELL_TYPES'):
        raise ValueError(f"{path} must define CELL_TYPES.")

    cfg.CONFIG_SOURCE = str(path)
    return cfg


def load_config(config_spec: str):
    """Load a bash-style config file path or, for compatibility, a Python module."""
    config_path = Path(config_spec).expanduser()
    if config_path.is_file() or config_path.suffix.lower() in {'.txt', '.conf', '.sh'}:
        return _parse_bash_config(str(config_path))
    return importlib.import_module(config_spec)


def read_sample_keep_file(
    path,
    sample_column="Sample",
):
    """Optional explicit sample whitelist."""
    if not path:
        return None

    sep = (
        "\t"
        if str(path).endswith(
            (".tsv", ".txt")
        )
        else ","
    )

    df = pd.read_csv(
        path,
        sep=sep,
    )

    if sample_column not in df.columns:
        raise ValueError(
            f"{path} must contain column "
            f"'{sample_column}'."
        )

    return set(
        df[sample_column].astype(str)
    )


# ==============================================================================
# Cell-type processing
# ==============================================================================

def process_cell_type(
    ct: str,
    cfg,
    gene_annot: pd.DataFrame,
):

    input_file = Path(
        cfg.PSEUDOBULK_INPUT_PATTERN.format(
            ct=ct
        )
    )

    if not input_file.exists():
        print(
            f"[{ct}] SKIP: input not found: "
            f"{input_file}"
        )
        return None

    out_dir = (
        Path(cfg.EXPRESSION_PREP_DIR)
        / ct
    )

    out_dir.mkdir(
        parents=True,
        exist_ok=True,
    )

    print(
        f"\n[{ct}] Reading: "
        f"{input_file}"
    )

    expr = pd.read_csv(
        input_file
    )

    if expr.shape[1] < 2:
        raise ValueError(
            f"{input_file} contains no samples."
        )

    input_gene_col = expr.columns[0]

    sample_cols = list(
        expr.columns[1:]
    )

    expr[sample_cols] = (
        expr[sample_cols]
        .apply(
            pd.to_numeric,
            errors="raise",
        )
    )

    # If gene names are duplicated in the Seurat input,
    # sum them before annotation.
    if expr[
        input_gene_col
    ].duplicated().any():

        print(
            f"[{ct}] Collapsing duplicated "
            "input genes by SUM."
        )

        expr = (
            expr.groupby(
                input_gene_col,
                as_index=False,
            )[sample_cols]
            .sum()
        )

    n_input_genes = expr.shape[0]

    n_input_samples = (
        len(sample_cols)
    )

    # --------------------------------------------------------------------------
    # GENCODE annotation
    # --------------------------------------------------------------------------

    ann_expr = annotation.annotate_expression_genes(
        expr,
        gene_annot,
    )

    metadata_cols = [
        "gene",
        "geneid",
        "genename",
        "gene_type",
        "CHR",
        "start",
        "end",
        "anno",
        "direction",
        "Length",
    ]

    sample_cols = [
        c for c in ann_expr.columns
        if c not in metadata_cols
    ]

    # --------------------------------------------------------------------------
    # Optional sample whitelist
    # --------------------------------------------------------------------------

    keep_set = read_sample_keep_file(
        getattr(
            cfg,
            "SAMPLE_KEEP_FILE",
            None,
        ),
        getattr(
            cfg,
            "SAMPLE_KEEP_COLUMN",
            "Sample",
        ),
    )

    if keep_set is not None:

        sample_cols = [
            sample
            for sample in sample_cols
            if sample in keep_set
        ]

    if len(sample_cols) == 0:
        raise RuntimeError(
            f"[{ct}] No samples remain."
        )

    counts = (
        ann_expr[sample_cols]
        .astype(float)
        .copy()
    )

    # --------------------------------------------------------------------------
    # Optional sample-level library-size filter
    # --------------------------------------------------------------------------

    min_library_counts = float(
        getattr(
            cfg,
            "MIN_SAMPLE_TOTAL_COUNTS",
            0,
        )
    )

    library_sizes = counts.sum(
        axis=0
    )

    keep_samples = (
        library_sizes[
            library_sizes
            >= min_library_counts
        ]
        .index
        .tolist()
    )

    dropped_samples = [
        sample
        for sample in sample_cols
        if sample not in keep_samples
    ]

    sample_cols = keep_samples

    counts = counts[
        sample_cols
    ]

    if counts.shape[1] == 0:
        raise RuntimeError(
            f"[{ct}] All samples removed "
            "by library-size QC."
        )

    # --------------------------------------------------------------------------
    # TPM
    # --------------------------------------------------------------------------

    gene_length_kb = (
        ann_expr["Length"]
        .to_numpy(dtype=float)
        / 1000.0
    )

    rate = (
        counts.to_numpy(dtype=float)
        / gene_length_kb[:, None]
    )

    rate_sum = rate.sum(
        axis=0
    )

    if np.any(rate_sum <= 0):

        bad_samples = np.array(
            sample_cols
        )[rate_sum <= 0].tolist()

        raise RuntimeError(
            f"[{ct}] Zero TPM denominator: "
            f"{bad_samples}"
        )

    tpm_raw = (
        rate
        / rate_sum[None, :]
        * 1e6
    )

    tpm_raw_df = pd.DataFrame(
        tpm_raw,
        index=ann_expr.index,
        columns=sample_cols,
    )

    # --------------------------------------------------------------------------
    # Gene expression QC
    #
    # TPM threshold is applied to RAW TPM, not log2(TPM+1).
    # --------------------------------------------------------------------------

    n_samples = len(
        sample_cols
    )

    min_samples_for_qc = int(
        math.ceil(
            float(
                cfg.SAMPLE_FRAC_THRESHOLD
            )
            * n_samples
        )
    )

    tpm_pass = (
        (
            tpm_raw_df
            >= float(
                cfg.TPM_THRESHOLD
            )
        )
        .sum(axis=1)
        >= min_samples_for_qc
    )

    count_pass = (
        (
            counts
            >= float(
                cfg.COUNT_THRESHOLD
            )
        )
        .sum(axis=1)
        >= min_samples_for_qc
    )

    qc_mask = (
        tpm_pass
        & count_pass
    )

    meta_qc = (
        ann_expr.loc[
            qc_mask,
            metadata_cols,
        ]
        .reset_index(drop=True)
    )

    counts_qc = (
        counts.loc[
            qc_mask
        ]
        .reset_index(drop=True)
    )

    tpm_qc = (
        tpm_raw_df.loc[
            qc_mask
        ]
        .reset_index(drop=True)
    )

    if meta_qc.shape[0] == 0:
        raise RuntimeError(
            f"[{ct}] No genes passed QC."
        )

    # --------------------------------------------------------------------------
    # log2(TPM + 1)
    # --------------------------------------------------------------------------

    log2_tpm_qc = np.log2(
        tpm_qc
        + float(
            cfg.LOG2_PSEUDOCOUNT
        )
    )

    log2_tpm_qc = pd.DataFrame(
        log2_tpm_qc,
        index=meta_qc.index,
        columns=sample_cols,
    )

    # --------------------------------------------------------------------------
    # Quantile normalization
    #
    # Reuse shared phenotype.quantile_normalize().
    # --------------------------------------------------------------------------

    print(
        f"[{ct}] Quantile normalization..."
    )

    qn = phenotype.quantile_normalize(
        log2_tpm_qc
    )

    phenotype.check_quantile_normalization(
        qn
    )

    # --------------------------------------------------------------------------
    # Inverse normal transformation
    #
    # Reuse shared phenotype.inverse_normal_transform().
    # This transforms each gene across samples.
    # --------------------------------------------------------------------------

    print(
        f"[{ct}] Inverse-normal transformation..."
    )

    inv = phenotype.inverse_normal_transform(
        qn
    )

    inv_diagnostic = (
        phenotype.check_inverse_normal_transform(
            inv
        )
    )

    # --------------------------------------------------------------------------
    # Outputs
    # --------------------------------------------------------------------------

    prefix = (
        f"pseudobulk_counts_{ct}"
    )

    final_file = (
        out_dir
        / (
            f"{prefix}."
            "tpmlog2.qc.qn."
            "rmlowcounts.inv.csv"
        )
    )

    final_df = pd.concat(
        [
            meta_qc[["gene"]]
            .reset_index(drop=True),
            inv.reset_index(drop=True),
        ],
        axis=1,
    )

    final_df.to_csv(
        final_file,
        index=False,
    )

    # GENCODE information for downstream BED / TensorQTL conversion.
    annotation_file = (
        out_dir
        / (
            f"{prefix}."
            "gene_annotation.qc.csv"
        )
    )

    meta_qc.to_csv(
        annotation_file,
        index=False,
    )

    # INT diagnostic.
    inv_diag_file = (
        out_dir
        / (
            f"{prefix}."
            "INV.diagnostic.csv"
        )
    )

    inv_diagnostic.insert(
        0,
        "gene",
        meta_qc["gene"].values,
    )

    inv_diagnostic.to_csv(
        inv_diag_file,
        index=False,
    )

    # --------------------------------------------------------------------------
    # Optional intermediates
    # --------------------------------------------------------------------------

    if getattr(
        cfg,
        "SAVE_INTERMEDIATES",
        True,
    ):

        # Raw TPM before gene filtering.
        pd.concat(
            [
                ann_expr[
                    ["gene"]
                ].reset_index(drop=True),
                tpm_raw_df.reset_index(
                    drop=True
                ),
            ],
            axis=1,
        ).to_csv(
            out_dir
            / f"{prefix}.tpm.csv",
            index=False,
        )

        # Raw counts after gene QC.
        pd.concat(
            [
                meta_qc[
                    ["gene"]
                ],
                counts_qc,
            ],
            axis=1,
        ).to_csv(
            out_dir
            / (
                f"{prefix}."
                "counts.qc.csv"
            ),
            index=False,
        )

        # log2 TPM after gene QC.
        pd.concat(
            [
                meta_qc[
                    ["gene"]
                ],
                log2_tpm_qc,
            ],
            axis=1,
        ).to_csv(
            out_dir
            / (
                f"{prefix}."
                "tpmlog2.qc.csv"
            ),
            index=False,
        )

        # Quantile normalized expression.
        pd.concat(
            [
                meta_qc[
                    ["gene"]
                ],
                qn.reset_index(
                    drop=True
                ),
            ],
            axis=1,
        ).to_csv(
            out_dir
            / (
                f"{prefix}."
                "tpmlog2.qc.qn.csv"
            ),
            index=False,
        )

    summary = {
        "cell_type": ct,
        "input_genes": n_input_genes,
        "annotated_autosomal_genes": ann_expr.shape[0],
        "genes_after_qc": meta_qc.shape[0],
        "input_samples": n_input_samples,
        "samples_after_filter": len(sample_cols),
        "samples_removed_low_library": len(
            dropped_samples
        ),
        "min_samples_required_for_gene_qc": (
            min_samples_for_qc
        ),
        "tpm_threshold": (
            cfg.TPM_THRESHOLD
        ),
        "count_threshold": (
            cfg.COUNT_THRESHOLD
        ),
        "sample_frac_threshold": (
            cfg.SAMPLE_FRAC_THRESHOLD
        ),
        "gene_annotation": (
            cfg.GENE_ANNOT_FILE
        ),
        "quantile_normalization_function": (
            "phenotype.quantile_normalize"
        ),
        "inverse_normal_function": (
            "phenotype.inverse_normal_transform"
        ),
        "final_file": str(
            final_file
        ),
    }

    print(
        f"[{ct}] "
        f"{n_input_genes} input genes -> "
        f"{ann_expr.shape[0]} "
        "annotated autosomal genes -> "
        f"{meta_qc.shape[0]} "
        "genes after QC"
    )

    print(
        f"[{ct}] "
        f"{n_input_samples} input samples -> "
        f"{len(sample_cols)} retained"
    )

    print(
        f"[{ct}] Wrote: "
        f"{final_file}"
    )

    return summary


# ==============================================================================
# Main
# ==============================================================================

def main():

    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--config",
        default=(
            "/v5000/data/l2_bioinfo1/liq17/scTWAS/data/"
            "VUMC_HTAN/src/HTAN_ADSER_9cts.eqtl.config_"
            "COLONMAP95_PCAfixed_pathfix_PEERlocal.txt"
        ),
        help=(
            "Bash-style tensorQTL config path, or legacy Python config "
            "module name without .py"
        ),
    )

    parser.add_argument(
        "--cell-types",
        nargs="*",
        default=None,
        help=(
            "Optional subset, e.g. "
            "--cell-types ABS STM TAC"
        ),
    )

    args = parser.parse_args()

    cfg = load_config(
        args.config
    )

    print("\nConfiguration")
    print("-------------")
    print(
        "CONFIG_SOURCE            =",
        getattr(cfg, "CONFIG_SOURCE", args.config),
    )

    print(
        "GENE_ANNOT_FILE          =",
        cfg.GENE_ANNOT_FILE,
    )

    print(
        "PSEUDOBULK_INPUT_PATTERN =",
        cfg.PSEUDOBULK_INPUT_PATTERN,
    )

    print(
        "EXPRESSION_PREP_DIR      =",
        cfg.EXPRESSION_PREP_DIR,
    )

    print(
        "TPM_THRESHOLD            =",
        cfg.TPM_THRESHOLD,
    )

    print(
        "COUNT_THRESHOLD          =",
        cfg.COUNT_THRESHOLD,
    )

    print(
        "SAMPLE_FRAC_THRESHOLD    =",
        cfg.SAMPLE_FRAC_THRESHOLD,
    )

    print(
        "LOG2_PSEUDOCOUNT         =",
        cfg.LOG2_PSEUDOCOUNT,
    )

    print(
        "\nBuilding GENCODE annotation "
        "and exon-union gene lengths..."
    )

    gene_annot = annotation.build_gene_annotation(
        cfg.GENE_ANNOT_FILE,
        autosomes_only=bool(
            getattr(cfg, "AUTOSOMES_ONLY", True)
        ),
    )

    print(
        "Annotated autosomal genes:",
        gene_annot.shape[0],
    )

    cell_types = (
        args.cell_types
        if args.cell_types
        else cfg.CELL_TYPES
    )

    summaries = []

    for ct in cell_types:

        summary = process_cell_type(
            ct,
            cfg,
            gene_annot,
        )

        if summary is not None:
            summaries.append(
                summary
            )

    if summaries:

        summary_df = pd.DataFrame(
            summaries
        )

        summary_file = (
            Path(
                cfg.EXPRESSION_PREP_DIR
            )
            / (
                "pseudobulk_expression_"
                "preprocessing_summary.csv"
            )
        )

        summary_file.parent.mkdir(
            parents=True,
            exist_ok=True,
        )

        summary_df.to_csv(
            summary_file,
            index=False,
        )

        print(
            "\nWrote preprocessing summary:",
            summary_file,
        )

    print("\nDone.")


if __name__ == "__main__":
    main()
