#!/usr/bin/env python3
"""
PEER_residual_to_tensorqtl_bed.py

Convert PEER residual expression to tensorQTL BED format.

Input expression:
    rows    = genes
    columns = expression samples
    first column = gene symbol or Ensembl gene ID

Output BED:
    #chr  start  end  gene_id  genotype_id_1  genotype_id_2 ...

Sample names are mapped:
    expression_id -> genotype_id

Gene positions use the same convention as the supplied annotation.py /
tensorQTL example:
    + strand: start = GTF start
    - strand: start = GTF end
    end = start + 1

By default, versioned GENCODE gene IDs are retained in the BED.
"""

import argparse
import gzip
import re
from pathlib import Path

import numpy as np
import pandas as pd


def open_text(path):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path, "rt")


def parse_gtf(gtf_file):
    rows = []
    with open_text(gtf_file) as fh:
        for line in fh:
            if not line or line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] != "gene":
                continue

            attrs = dict(re.findall(r'(\S+)\s+"([^"]+)"', f[8]))
            gene_id = attrs.get("gene_id")
            gene_name = attrs.get("gene_name")
            if not gene_id or not gene_name:
                continue

            chrom = f[0]
            if not re.fullmatch(r"chr(?:[1-9]|1[0-9]|2[0-2])", chrom):
                continue

            start = int(f[3])
            end = int(f[4])
            strand = f[6]
            tss = start if strand == "+" else end

            rows.append({
                "gene_id": gene_id,
                "gene_id_nover": re.sub(r"\.\d+$", "", gene_id),
                "gene_name": gene_name,
                "#chr": chrom,
                "start": tss,
                "end": tss + 1,
                "strand": strand,
                "gene_length": end - start,
            })

    gtf = pd.DataFrame(rows)
    if gtf.empty:
        raise RuntimeError(f"No autosomal gene records were read from {gtf_file}")

    # Unique Ensembl ID; for ambiguous symbols, prefer the longest gene record.
    gtf = gtf.sort_values("gene_length", ascending=False)
    return gtf


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--expression", required=True)
    ap.add_argument("--mapping", required=True)
    ap.add_argument("--gtf", required=True)
    ap.add_argument("--output", required=True)
    ap.add_argument("--audit-prefix", default=None)
    ap.add_argument("--expression-id-col", default="expression_id")
    ap.add_argument("--genotype-id-col", default="genotype_id")
    ap.add_argument(
        "--gene-id-style",
        choices=["versioned", "unversioned"],
        default="versioned",
    )
    args = ap.parse_args()

    out = Path(args.output)
    out.parent.mkdir(parents=True, exist_ok=True)

    audit_prefix = (
        Path(args.audit_prefix)
        if args.audit_prefix
        else out.with_suffix("").with_suffix("")
    )
    audit_prefix.parent.mkdir(parents=True, exist_ok=True)

    # ------------------------------------------------------------------
    # Sample mapping
    # ------------------------------------------------------------------
    smap = pd.read_csv(args.mapping, dtype=str)

    required = {args.expression_id_col, args.genotype_id_col}
    missing = required.difference(smap.columns)
    if missing:
        raise RuntimeError(
            "Mapping file missing column(s): " + ", ".join(sorted(missing))
        )

    smap = smap[[args.expression_id_col, args.genotype_id_col]].copy()
    smap.columns = ["expression_id", "genotype_id"]
    smap["expression_id"] = smap["expression_id"].astype(str).str.strip()
    smap["genotype_id"] = smap["genotype_id"].astype(str).str.strip()

    if smap["expression_id"].duplicated().any():
        raise RuntimeError("Mapping file contains duplicated expression_id.")
    if smap["genotype_id"].duplicated().any():
        raise RuntimeError("Mapping file contains duplicated genotype_id.")

    expr_to_geno = dict(zip(smap["expression_id"], smap["genotype_id"]))
    geno_set = set(smap["genotype_id"])

    # ------------------------------------------------------------------
    # Expression
    # ------------------------------------------------------------------
    expr = pd.read_csv(args.expression, compression="infer")
    if expr.shape[1] < 2:
        raise RuntimeError("Expression file must have gene + sample columns.")

    gene_col = expr.columns[0]
    input_samples = [str(x) for x in expr.columns[1:]]

    renamed = []
    keep_cols = []
    sample_audit = []

    for sid in input_samples:
        if sid in expr_to_geno:
            gid = expr_to_geno[sid]
            source = "expression_id"
        elif sid in geno_set:
            gid = sid
            source = "genotype_id"
        else:
            gid = None
            source = "unmapped"

        sample_audit.append({
            "input_expression_label": sid,
            "genotype_id": gid,
            "mapping_type": source,
        })

        if gid is not None:
            keep_cols.append(sid)
            renamed.append(gid)

    sample_audit = pd.DataFrame(sample_audit)
    sample_audit.to_csv(
        str(audit_prefix) + ".sample_mapping.tsv",
        sep="\t",
        index=False,
    )

    if not renamed:
        raise RuntimeError("No expression samples mapped to genotype_id.")
    if len(set(renamed)) != len(renamed):
        raise RuntimeError("Sample mapping produced duplicated genotype_id.")

    expr = expr[[gene_col] + keep_cols].copy()
    expr.columns = [gene_col] + renamed

    # Preserve covariate/mapping-file order.
    ordered_samples = [
        x for x in smap["genotype_id"].tolist()
        if x in set(renamed)
    ]
    expr = expr[[gene_col] + ordered_samples]

    # ------------------------------------------------------------------
    # Gene mapping
    # ------------------------------------------------------------------
    gtf = parse_gtf(args.gtf)

    by_id = (
        gtf.drop_duplicates("gene_id_nover", keep="first")
        .set_index("gene_id_nover", drop=False)
    )
    by_name = (
        gtf.drop_duplicates("gene_name", keep="first")
        .set_index("gene_name", drop=False)
    )

    mapped_rows = []
    gene_audit = []

    for i, raw in enumerate(expr[gene_col].astype(str)):
        raw = raw.strip()
        nover = re.sub(r"\.\d+$", "", raw)

        ann = None
        mapping_type = None

        if nover in by_id.index:
            ann = by_id.loc[nover]
            mapping_type = "ENSG"
        elif raw in by_name.index:
            ann = by_name.loc[raw]
            mapping_type = "gene_symbol"

        if ann is None:
            gene_audit.append({
                "input_gene": raw,
                "gene_id": None,
                "gene_id_nover": None,
                "mapping_type": "unmapped",
            })
            continue

        if isinstance(ann, pd.DataFrame):
            ann = ann.iloc[0]

        gene_id_out = (
            ann["gene_id"]
            if args.gene_id_style == "versioned"
            else ann["gene_id_nover"]
        )

        mapped_rows.append({
            "expr_row": i,
            "#chr": ann["#chr"],
            "start": int(ann["start"]),
            "end": int(ann["end"]),
            "gene_id": gene_id_out,
            "gene_id_nover": ann["gene_id_nover"],
            "mapping_type": mapping_type,
        })

        gene_audit.append({
            "input_gene": raw,
            "gene_id": ann["gene_id"],
            "gene_id_nover": ann["gene_id_nover"],
            "mapping_type": mapping_type,
        })

    gene_audit = pd.DataFrame(gene_audit)
    gene_audit.to_csv(
        str(audit_prefix) + ".gene_mapping.tsv.gz",
        sep="\t",
        index=False,
        compression="gzip",
    )

    mapped = pd.DataFrame(mapped_rows)
    if mapped.empty:
        raise RuntimeError("No expression genes mapped to the GTF.")

    # One phenotype per output gene ID.
    mapped = mapped.drop_duplicates("gene_id", keep="first").copy()

    values = expr.iloc[mapped["expr_row"].to_numpy()][ordered_samples].copy()
    values = values.apply(pd.to_numeric, errors="raise")

    if not np.isfinite(values.to_numpy(dtype=float)).all():
        raise RuntimeError("Expression BED contains NA/NaN/Inf.")

    bed = pd.concat(
        [
            mapped[["#chr", "start", "end", "gene_id"]].reset_index(drop=True),
            values.reset_index(drop=True),
        ],
        axis=1,
    )

    chrom_order = {f"chr{i}": i for i in range(1, 23)}
    bed["_chr_order"] = bed["#chr"].map(chrom_order)
    bed = (
        bed.sort_values(["_chr_order", "start", "end", "gene_id"])
        .drop(columns="_chr_order")
        .reset_index(drop=True)
    )

    bed.to_csv(
        out,
        sep="\t",
        index=False,
        compression="gzip" if str(out).endswith(".gz") else None,
    )

    print("=" * 72)
    print("tensorQTL BED created")
    print("=" * 72)
    print("Expression :", args.expression)
    print("Mapping    :", args.mapping)
    print("GTF        :", args.gtf)
    print("Output     :", out)
    print("Samples    :", len(ordered_samples))
    print("Genes      :", bed.shape[0])
    print("Gene ID    :", args.gene_id_style)
    print("Header     :", ", ".join(bed.columns[:6]) + " ...")


if __name__ == "__main__":
    main()
