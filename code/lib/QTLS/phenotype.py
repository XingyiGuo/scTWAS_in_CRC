"""
phenotype_prep.py
================================================================
把各种表型数据（基因表达 / 单细胞表达 / APA / 可变剪接）转换成
tensorQTL 需要的 BED 格式 （#chr, start, end, gene_id, sample1, sample2, ...）。
#chr    start   end     gene_id VM1003  
chr1    29570   29571   ENSG00000227232.5       -0.8654488973804301
"""

import pandas as pd
import numpy as np
from annotation import *
import scipy.stats as stats

def deconv_sc_expr_to_bed_tss(df: pd.DataFrame, gene_annot_file: str,
                              id_col: str = "geneid") -> pd.DataFrame:
    """
    UVA / ACCC 单细胞表达转 BED（按 geneid 而不是 genename 合并）。
    对应旧代码里 uva_sc_df_annot / accc_sc_df_annot 那两段。

    df 需要已经 drop 掉 GeneSymbol / Length 列，只剩 geneid + 样本列。
    """
    tss_table = get_tss_table(gene_annot_file)
    df=df.drop(columns=["GeneSymbol","Length"])
    df.rename(columns={id_col:"gene_id"}, inplace=True)
    merged = pd.merge(tss_table, df, on="gene_id", how="right")
    merged = merged.dropna(subset=["start"])
    merged["start"]=merged["start"].astype("int")
    merged["end"]=merged["end"].astype("int")
    return merged.sort_values(["#chr", "start", "end"])


def expr_to_bed_tss(expr_df: pd.DataFrame, gene_annot_file: str,
                     sample_id_row: str = "sampleID") -> pd.DataFrame:
    """
    把"行=样本，列=基因"的表达矩阵转置后，用 TSS 坐标生成 BED。
    对应旧代码：UVA_Expr / ACCC_Expr 那两段（uva_expr_t_tensorqtl / accc_expr_t_tensorqtl）。

    expr_df: 原始读入的表达矩阵，第一行（或某一行）是样本 ID，
             转置后行变成基因，列变成样本。
    """
    t = expr_df.transpose()
    t.columns = t.iloc[0, :]
    t = t.drop(index=(sample_id_row))
    t["gene_id"] = t.index
    t = t.reset_index(drop=True)
    tss_table = get_tss_table(gene_annot_file)
    merged = pd.merge(tss_table, t, on="gene_id", how="right")
    merged = merged.dropna(subset=["start"])
    merged["start"]=merged["start"].astype("int")
    merged["end"]=merged["end"].astype("int")
    return merged.sort_values(["#chr", "start", "end"])


def apa_to_bed(apa_df: pd.DataFrame, id_col: str = "GeneID") -> pd.DataFrame:
    """
    APA（poly-A 使用率）表型转 BED。
    GeneID 格式形如 "...chr1:12345-12399"，从中解析出 chr/start/end。
    对应旧代码：uva_apa_tensorqtl / accc_apa_tensorqtl。
    """
    content = apa_df.iloc[:, 1:]
    chrom_list, start_list, end_list = [], [], []
    for apa_id in apa_df[id_col]:
        region = apa_id.split("_")[-1]          # "chr1:12345-12399"
        chrom, pos = region.split(":")
        start, end = pos.split("-")
        chrom_list.append(chrom)
        start_list.append(start)
        end_list.append(end)

    prefix = apa_df[[id_col]].copy()
    prefix["#chr"] = chrom_list
    prefix["start"] = start_list
    prefix["end"] = end_list
    prefix = prefix.rename(columns={id_col: "gene_id"})[["#chr", "start", "end", "gene_id"]]

    out = pd.concat([prefix, content], axis=1)
    out["start"] = out["start"].astype(int)
    out["end"] = out["end"].astype(int)
    return out.sort_values(["#chr", "start", "end"])


def altersplice_to_bed(alters_df: pd.DataFrame, sample_id_row: str = "sampleID",
                       add_chr_prefix: bool = False) -> pd.DataFrame:
    """
    altersplice 可变剪接定量（intron cluster ratio）转 BED。
    gene_id 格式形如 "1:12345:12399:clu_1_+"，从中解析 chr/start/end。
    对应旧代码：uva_alters_t_tensorqtl / accc_alters_t_tensorqtl。

    add_chr_prefix: ACCC 那一版额外加了 "chr" 前缀（旧代码不一致，UVA 没加），
                    根据实际基因型文件里 chr 的格式来定。
    """
    t = alters_df.transpose()
    t.columns = t.iloc[0, :]
    t = t.drop(index=(sample_id_row))
    t["gene_id"] = t.index
    t = t.reset_index(drop=True)
    content = t.iloc[:, :-1]
    t["#chr"] = t["gene_id"].str.split(":").str[0]
    t["start"] = t["gene_id"].str.split(":").str[1]
    t["end"] = t["gene_id"].str.split(":").str[2]
    if add_chr_prefix:
        t["#chr"] = "chr" + t["#chr"]
    out = pd.concat([t[["#chr", "start", "end", "gene_id"]], content], axis=1)
    out["start"] = out["start"].astype(int)
    out["end"] = out["end"].astype(int)
    return out.sort_values(["#chr", "start", "end"])
 

"""
Quantile normalization to the average empirical distribution

Note: replicates behavior of R function normalize.quantiles
      from library("preprocessCore")

Reference:
 [1] Bolstad et al., Bioinformatics 19(2), pp. 185-193, 2003

Adapted from https://github.com/andrewdyates/quantile_normalize
""" 
def quantile_normalize(df: pd.DataFrame) -> pd.DataFrame:
    M = df.to_numpy(copy=True)
    Q = M.argsort(axis=0)
    m, n = M.shape
    # compute quantile vector
    quantiles = np.zeros(m)
    for i in range(n):
        quantiles += M[Q[:, i], i]
    quantiles = quantiles / n
    for i in range(n):
        # Get equivalence classes; unique values == 0
        dupes = np.zeros(m, dtype=np.int64)
        for j in range(m - 1):
            if M[Q[j, i], i] == M[Q[j + 1, i], i]:
                dupes[j + 1] = dupes[j] + 1
        # Replace column with quantile ranks
        M[Q[:, i], i] = quantiles
        # Average together equivalence classes
        j = m - 1
        while j >= 0:
            if dupes[j] == 0:
                j -= 1
            else:
                idxs = Q[j - dupes[j]:j + 1, i]
                M[idxs, i] = np.median(M[idxs, i])
                j -= 1 + dupes[j]
        assert j == -1
    return pd.DataFrame(M, index=df.index, columns=df.columns)


def inverse_normal_transform(df: pd.DataFrame) -> pd.DataFrame:
    M = df.to_numpy(copy=True)
    """Transform rows to a standard normal distribution"""
    R = stats.rankdata(M, axis=1)  # ties are averaged
    Q = stats.norm.ppf(R / (M.shape[1] + 1))
    Q = pd.DataFrame(Q, index=df.index, columns=df.columns)
    return Q
    
    
def check_quantile_normalization(df):
    sorted_values = np.sort(df.to_numpy(), axis=0)
    ref = sorted_values[:, 0]
    max_diff = 0
    for i in range(1, sorted_values.shape[1]):
        d = np.max(np.abs(sorted_values[:, i] - ref))
        max_diff = max(max_diff, d)
    print(f"Maximum quantile difference: {max_diff:.3e}")
    
    
def check_inverse_normal_transform(df: pd.DataFrame) -> pd.DataFrame:
    """
    Check whether each feature (row) approximately follows N(0,1).

    Returns
    -------
    DataFrame
        Per-feature mean and standard deviation.
    """
    return pd.DataFrame({
        "mean": df.mean(axis=1),
        "sd": df.std(axis=1)
    })