"""
gene_annotation.py
================================================================
解析 GTF 基因注释文件。

旧代码里有两份几乎一样但不完全一致的 get_gene_annotation()：
  - 一份提取 geneid / transid / genename / gene_type
  - 一份只提取 geneid / genename
这里合并成一个函数，transid/gene_type 默认提取，不需要的话忽略那两列即可，
避免维护两份重复逻辑。
"""

import re
import gzip
from collections import defaultdict
import pandas as pd
import numpy as np
def get_gene_annotation(gene_annot_file_name: str) -> pd.DataFrame:
    """
    读取 GTF 文件，只保留 feature == "gene" 的行，提取常用字段。
    Returns
    -------
    DataFrame with columns:
        geneid, transid, genename, gene_type, chr, start, end, anno, direction
    """
    gene_df = pd.read_csv(
        gene_annot_file_name, comment="#", header=None, sep="\t",
        dtype=str, engine="python",
    )
    gene_df1 = gene_df[gene_df[2] == "gene"]

    def _extract(pattern, attr_series, group=0):
        return attr_series.apply(
            lambda x: re.search(pattern, x).group(group) if re.search(pattern, x) else None
        )
    attrs = gene_df1[8]
    geneid = _extract(r"ENSG\d+\.\d+", attrs)
    transid = _extract(r"ENST\d+\.\d+", attrs)
    genename = _extract(r"gene_name (\S+);", attrs, group=1)
    genetype = _extract(r"gene_type (\S+);", attrs, group=1)
    gene_used = pd.DataFrame({
        "geneid": geneid,
        "transid": transid,
        "genename": genename,
        "gene_type": genetype,
        "chr": gene_df1[0],
        "start": gene_df1[3],
        "end": gene_df1[4],
        "anno": gene_df1[2],
        "direction": gene_df1[6],
    })
    for col in ["transid", "genename", "gene_type"]:
        gene_used[col] = gene_used[col].str.replace('"', "", regex=False)
    gene_used = gene_used.drop_duplicates(subset=["geneid"])
    gene_used = gene_used[~gene_used["chr"].isin(['chrX', 'chrY', 'chrM'])]
    return gene_used


def get_tss_table(gene_annot_file_name: str) -> pd.DataFrame:
    """
    生成以 TSS 为中心点的基因坐标表，用于 gene_tss 类型的表型（如 bulk 表达量）。
    + 链用 start 作为 TSS，- 链用 end 作为 TSS；输出区间为 [tss, tss+1)。

    Returns columns: #chr, start, end, gene_id
    """
    gtf_df = get_gene_annotation(gene_annot_file_name)
    gtf_df["tss"] = np.where(gtf_df["direction"] == "+", gtf_df["start"], gtf_df["end"])
    gtf_short = gtf_df[["chr", "tss", "geneid"]].rename(
        columns={"chr": "#chr", "geneid": "gene_id"}
    )
    gtf_short["tss"] = gtf_short["tss"].astype(int)
    gtf_short["start"] = gtf_short["tss"]
    gtf_short["end"] = gtf_short["tss"] + 1
    return gtf_short[["#chr", "start", "end", "gene_id"]]


# ==============================================================================
# Expression / TPM annotation helpers
# ==============================================================================

def _parse_gtf_attributes(attr_string: str) -> dict:
    """Parse standard GTF key-value attributes."""
    return dict(re.findall(r'(\S+)\s+"([^"]+)"', attr_string))


def _merge_intervals(intervals) -> int:
    """
    Return the union length of 1-based closed genomic intervals.

    Overlapping exons from multiple transcripts are counted only once.
    """
    if not intervals:
        return 0

    intervals = sorted(intervals)
    cur_start, cur_end = intervals[0]
    total = 0

    for start, end in intervals[1:]:
        if start <= cur_end + 1:
            cur_end = max(cur_end, end)
        else:
            total += cur_end - cur_start + 1
            cur_start, cur_end = start, end

    total += cur_end - cur_start + 1
    return int(total)


def get_exon_union_lengths(gtf_file: str) -> pd.DataFrame:
    """
    Calculate exon-union length for each GENCODE gene.

    TPM requires a gene-length denominator. For each Ensembl gene ID, this
    function unions all exon intervals across transcripts so overlapping exons
    are counted only once.

    Parameters
    ----------
    gtf_file : str
        GENCODE GTF file (.gtf or .gtf.gz).

    Returns
    -------
    pandas.DataFrame
        Columns:
            geneid_nover
            Length
    """
    exon_intervals = defaultdict(list)

    opener = gzip.open if str(gtf_file).endswith(".gz") else open

    with opener(gtf_file, "rt") as fh:
        for line in fh:
            if not line or line.startswith("#"):
                continue

            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9 or fields[2] != "exon":
                continue

            attrs = _parse_gtf_attributes(fields[8])
            gene_id = attrs.get("gene_id")
            if gene_id is None:
                continue

            gene_id_nover = gene_id.split(".")[0]
            start = int(fields[3])
            end = int(fields[4])

            exon_intervals[gene_id_nover].append((start, end))

    rows = []
    for gene_id_nover, intervals in exon_intervals.items():
        length = _merge_intervals(intervals)
        if length > 0:
            rows.append({
                "geneid_nover": gene_id_nover,
                "Length": length,
            })

    length_df = pd.DataFrame(rows)

    if length_df.empty:
        raise RuntimeError(
            f"No exon lengths could be calculated from {gtf_file}"
        )

    return length_df


def build_gene_annotation(gtf_file: str, autosomes_only: bool = True) -> pd.DataFrame:
    """
    Build a gene-level annotation table for expression / TPM preprocessing.

    Reuses get_gene_annotation() for GENCODE gene metadata, then adds exon-union
    gene length.

    Parameters
    ----------
    gtf_file : str
        GENCODE GTF file.
    autosomes_only : bool, default True
        If True, retain chr1-chr22 only.

    Returns
    -------
    pandas.DataFrame
        Columns:
            geneid
            geneid_nover
            genename
            gene_type
            CHR
            start
            end
            anno
            direction
            Length
    """
    gene_annot = get_gene_annotation(gtf_file).copy()

    required = {
        "geneid",
        "genename",
        "gene_type",
        "chr",
        "start",
        "end",
        "anno",
        "direction",
    }

    missing = required.difference(gene_annot.columns)
    if missing:
        raise ValueError(
            "get_gene_annotation() is missing columns: "
            + ", ".join(sorted(missing))
        )

    gene_annot["geneid"] = gene_annot["geneid"].astype(str)
    gene_annot["geneid_nover"] = (
        gene_annot["geneid"]
        .str.replace(r"\.\d+$", "", regex=True)
    )
    gene_annot["genename"] = gene_annot["genename"].astype(str)

    if autosomes_only:
        autosomes = {f"chr{i}" for i in range(1, 23)}
        gene_annot = gene_annot[
            gene_annot["chr"].isin(autosomes)
        ].copy()

    length_df = get_exon_union_lengths(gtf_file)

    gene_annot = gene_annot.merge(
        length_df,
        on="geneid_nover",
        how="left",
    )

    gene_annot = gene_annot.dropna(
        subset=["Length"]
    ).copy()

    gene_annot["Length"] = gene_annot["Length"].astype(int)
    gene_annot["start"] = gene_annot["start"].astype(int)
    gene_annot["end"] = gene_annot["end"].astype(int)

    gene_annot = gene_annot.rename(
        columns={"chr": "CHR"}
    )

    gene_annot = gene_annot.sort_values(
        ["CHR", "start"]
    ).reset_index(drop=True)

    return gene_annot


def annotate_expression_genes(
    expr: pd.DataFrame,
    gene_annot: pd.DataFrame,
) -> pd.DataFrame:
    """
    Match expression genes to GENCODE annotation.

    Input
    -----
    expr
        First column = gene identifier / gene symbol.
        Remaining columns = sample expression/count values.

    Matching order
    --------------
    1. Ensembl gene ID without version.
    2. Gene symbol (genename).

    If multiple GENCODE IDs share the same gene symbol, the annotation with the
    largest exon-union length is retained.

    Returns
    -------
    pandas.DataFrame
        Metadata columns followed by sample columns:
            gene, geneid, genename, gene_type, CHR, start, end,
            anno, direction, Length, sample1, sample2, ...
    """
    required_annot_cols = {
        "geneid",
        "geneid_nover",
        "genename",
        "gene_type",
        "CHR",
        "start",
        "end",
        "anno",
        "direction",
        "Length",
    }

    missing_annot_cols = required_annot_cols.difference(gene_annot.columns)

    if missing_annot_cols:
        raise ValueError(
            "gene_annot is missing required columns: "
            + ", ".join(sorted(missing_annot_cols))
            + f"; available columns: {list(gene_annot.columns)}"
        )

    input_gene_col = expr.columns[0]

    x = expr.rename(
        columns={input_gene_col: "gene"}
    ).copy()

    x["gene"] = x["gene"].astype(str)
    x["gene_nover"] = (
        x["gene"]
        .str.replace(r"\.\d+$", "", regex=True)
    )

    sample_cols = [
        c for c in x.columns
        if c not in {"gene", "gene_nover"}
    ]

    by_id = (
        gene_annot
        .drop_duplicates("geneid_nover", keep="first")
        .set_index("geneid_nover", drop=False)
    )

    by_name = (
        gene_annot
        .sort_values("Length", ascending=False)
        .drop_duplicates("genename", keep="first")
        .set_index("genename", drop=False)
    )

    rows = []

    for _, row in x.iterrows():
        gene_nover = row["gene_nover"]
        gene_name = row["gene"]

        ann = None

        if gene_nover in by_id.index:
            ann = by_id.loc[gene_nover]

        elif gene_name in by_name.index:
            ann = by_name.loc[gene_name]

        if ann is None:
            continue

        # Both lookup tables are deduplicated above, so each lookup should
        # return exactly one Series. Keep this guard in case future annotation
        # changes introduce duplicated index entries.
        if isinstance(ann, pd.DataFrame):
            ann = ann.iloc[0]

        rec = {
            "gene": row["gene"],
            "geneid": ann["geneid"],
            "genename": ann["genename"],
            "gene_type": ann["gene_type"],
            "CHR": ann["CHR"],
            "start": int(ann["start"]),
            "end": int(ann["end"]),
            "anno": ann["anno"],
            "direction": ann["direction"],
            "Length": int(ann["Length"]),
        }

        for sample in sample_cols:
            rec[sample] = row[sample]

        rows.append(rec)

    annotated = pd.DataFrame(rows)

    if annotated.empty:
        raise RuntimeError(
            "No expression genes matched GENCODE annotation."
        )

    return annotated

