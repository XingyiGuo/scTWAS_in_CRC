#!/usr/bin/env python3
"""06_test_cptac_proteome_wide.py

Proteome-wide tumour vs adjacent-normal differential abundance for CPTAC-COAD.

Where this sits relative to the existing code
    05_annotate_cptac_protein_support.py answers "what does CPTAC
    say about these 347 scTWAS genes". This script answers the prior question
    "what does CPTAC say about every protein it measured", so the result is a
    reusable reference table rather than a gene-set annotation. The statistics
    are deliberately identical to that script (paired t-test, BH correction) so
    the two agree where they overlap; that agreement is asserted at the end.

Cohort
    CPTAC-2 Prospective Colon Cancer. Vasaikar et al., Cell 2019,
    PMID 31031003, doi 10.1016/j.cell.2019.03.030.
    Accessed through the cptac Python package, umich TMT proteomics
    (gene/protein level, median-normalized log2 abundance).

Output
    data/CPTAC_COAD_proteome_wide_DE_tumor_vs_normal.csv
    data/CPTAC_COAD_proteome_wide_DE_README.md   (methods, counts, validation)

Columns
    gene_symbol, protein_id            protein identity as distributed
    log2FC_tumor_vs_normal             mean paired difference, + = up in tumour
    p_value, padj_BH                   paired t-test, BH over tested proteins
    t_statistic
    n_pairs, pct_pairs_quantified      how much data each protein rests on
    mean_log2_tumor, mean_log2_normal
    tested                             n_pairs >= MIN_PAIRS; padj_BH is NaN otherwise
    gene_level_representative          one row per symbol, for gene-level joins
    dataset .. doi, data_access        provenance, constant down the column

Install notes, all three learned the hard way on this machine
    1. cptac needs pkg_resources, which setuptools >= 81 no longer ships:
           pip install "setuptools<81"
    2. cptac writes its download cache next to its own package files. If those
       are read-only, install cptac somewhere writable or set
       cptac.tools.download_tools.DATA_DIR before first use.
    3. cptac 1.5.14's threaded downloader divides by a file size that Zenodo
       does not always report, raising
           TypeError: unsupported operand type(s) for //: 'NoneType' and 'int'
       patch_cptac_downloader() below replaces it with a plain stream. This is
       the same failure the streamfix in the older script works around.
"""

from __future__ import annotations

from importlib.metadata import version as pkg_version
from pathlib import Path
import warnings

import numpy as np
import pandas as pd
import requests
from scipy import stats
from statsmodels.stats.multitest import multipletests

warnings.filterwarnings("ignore")

PROJ = Path(r"W:/VUMC_projects/scTWAS/VUMC_HTAN_scTWAS_genes_DE_support")
OUT_CSV = PROJ / "data" / "CPTAC_COAD_proteome_wide_DE_tumor_vs_normal.csv"

SOURCE = "umich"      # 'bcm' fails inside the package's own loader, see README
MIN_PAIRS = 20        # a protein needs this many matched pairs to be tested
FDR_CUT = 0.05

PUBLICATION = {
    "publication": "Vasaikar S et al. Cell 2019",
    "pmid": "31031003",
    "doi": "10.1016/j.cell.2019.03.030",
}


def patch_cptac_downloader() -> None:
    """Replace cptac's threaded Zenodo downloader with a single stream.

    The shipped implementation computes ``chunk_size = file_size //
    num_threads`` from Zenodo record metadata that can carry a null size. The
    URL it has is correct, so streaming it directly avoids needing the size at
    all. Relative destinations are resolved against the package DATA_DIR,
    which is what the original does internally.
    """
    import cptac.tools.download_tools as dt

    def _stream(url: str, output_file: str = "", *args, **kwargs) -> str:
        dest = Path(output_file)
        if not dest.is_absolute():
            dest = Path(dt.DATA_DIR) / dest
        dest.parent.mkdir(parents=True, exist_ok=True)
        with requests.get(url, stream=True, timeout=600) as resp:
            resp.raise_for_status()
            with open(dest, "wb") as fh:
                for chunk in resp.iter_content(1 << 20):
                    if chunk:
                        fh.write(chunk)
        return str(dest)

    dt.get_data = _stream


def load_proteomics(source: str = SOURCE) -> pd.DataFrame:
    """Return the samples x proteins log2 abundance matrix."""
    import cptac

    patch_cptac_downloader()
    return cptac.Coad().get_proteomics(source)


def paired_de(prot: pd.DataFrame, min_pairs: int = MIN_PAIRS) -> pd.DataFrame:
    """Paired tumour-vs-normal t-test per protein.

    Normal samples carry a ``.N`` suffix on the patient id, which is how the
    source distinguishes the two arms; stripping it recovers the patient and
    therefore the pairing. Only patients present in BOTH arms are used, so the
    unmatched tumours and normals are dropped rather than quietly turning the
    comparison into a two-sample test.
    """
    names = prot.columns.get_level_values("Name").astype(str)
    dbid = prot.columns.get_level_values("Database_ID").astype(str)

    idx = prot.index.astype(str)
    is_normal = idx.str.endswith(".N")
    patient = np.where(is_normal, idx.str.replace(r"\.N$", "", regex=True), idx)

    tumour = prot.loc[~is_normal].copy()
    tumour.index = patient[~is_normal]
    normal = prot.loc[is_normal].copy()
    normal.index = patient[is_normal]

    pairs = sorted(set(tumour.index) & set(normal.index))
    if not pairs:
        raise RuntimeError("no patient has both a tumour and a normal sample")
    print(f"{prot.shape[0]} samples, {int((~is_normal).sum())} tumour / "
          f"{int(is_normal.sum())} normal, {len(pairs)} matched pairs")

    tm = tumour.loc[pairs].to_numpy(float)
    nm = normal.loc[pairs].to_numpy(float)
    diff = tm - nm

    n_pairs = np.sum(~np.isnan(diff), axis=0)
    with np.errstate(invalid="ignore"):
        t_stat, p_val = stats.ttest_rel(tm, nm, axis=0, nan_policy="omit")

    res = pd.DataFrame({
        "gene_symbol": names,
        "protein_id": dbid,
        "log2FC_tumor_vs_normal": np.nanmean(diff, axis=0),
        "t_statistic": np.asarray(t_stat, float),
        "p_value": np.asarray(p_val, float),
        "n_pairs": n_pairs,
        "pct_pairs_quantified": np.round(100 * n_pairs / len(pairs), 1),
        "mean_log2_tumor": np.nanmean(tm, axis=0),
        "mean_log2_normal": np.nanmean(nm, axis=0),
    })

    ## BH is applied over the tested set only. Including proteins with no p
    ## value would inflate the number of hypotheses and deflate every FDR.
    tested = (res.n_pairs >= min_pairs) & res.p_value.notna()
    res["tested"] = tested
    res.loc[tested, "padj_BH"] = multipletests(
        res.loc[tested, "p_value"], method="fdr_bh")[1]

    ## one row per symbol for gene-level joins: the best-supported entry
    rep = (res[tested]
           .sort_values(["n_pairs", "p_value"], ascending=[False, True])
           .drop_duplicates("gene_symbol"))
    rep_keys = set(zip(rep.gene_symbol, rep.protein_id))
    res["gene_level_representative"] = [
        (g, p) in rep_keys for g, p in zip(res.gene_symbol, res.protein_id)]
    return res


def annotate_provenance(res: pd.DataFrame, source: str = SOURCE) -> pd.DataFrame:
    res = res.copy()
    res["dataset"] = "CPTAC-2 Prospective Colon Cancer (COAD)"
    res["assay"] = ("TMT global proteomics, gene/protein level, "
                    f"MD-normalized ({source} pipeline)")
    res["contrast"] = "primary tumor vs patient-matched adjacent normal colon"
    res["test"] = "paired two-sided t-test on log2 abundance differences"
    res["padj_method"] = "Benjamini-Hochberg over tested proteins"
    for k, v in PUBLICATION.items():
        res[k] = v
    res["data_access"] = (f"cptac Python package v{pkg_version('cptac')}, "
                          f"Coad().get_proteomics('{source}')")
    return res


def check_against_existing(res: pd.DataFrame) -> None:
    """Assert agreement with the 347-gene table from the older script.

    Both analyses use the same cohort and the same test, so the log2FC values
    should be near-identical. This is the guard that catches a silently changed
    upstream file or a broken pairing.
    """
    old_path = PROJ / "data" / "TCGA_CPTAC_support.csv"
    if not old_path.exists():
        print("existing table absent; skipping cross-check")
        return
    old = (pd.read_csv(old_path, skiprows=1)
           .drop_duplicates("gene_name")
           .dropna(subset=["CPTAC_log2FC"]))
    rep = res[res.gene_level_representative]
    m = old.merge(rep, left_on="gene_name", right_on="gene_symbol")
    r = m.CPTAC_log2FC.corr(m.log2FC_tumor_vs_normal)
    sign_ok = int((np.sign(m.CPTAC_log2FC) == np.sign(m.log2FC_tumor_vs_normal)).sum())
    print(f"cross-check on {len(m)} genes: log2FC r = {r:.4f}, "
          f"sign agreement {sign_ok}/{len(m)}")
    assert r > 0.95, f"log2FC disagrees with the existing table (r={r:.3f})"
    assert sign_ok == len(m), "a gene changed direction versus the existing table"


def main() -> None:
    prot = load_proteomics()
    res = annotate_provenance(paired_de(prot))

    order = ["gene_symbol", "protein_id", "log2FC_tumor_vs_normal", "p_value",
             "padj_BH", "t_statistic", "n_pairs", "pct_pairs_quantified",
             "mean_log2_tumor", "mean_log2_normal", "tested",
             "gene_level_representative", "dataset", "assay", "contrast",
             "test", "padj_method", "publication", "pmid", "doi",
             "data_access"]
    res = res[order].sort_values(["padj_BH", "p_value"], na_position="last")

    check_against_existing(res)
    OUT_CSV.parent.mkdir(parents=True, exist_ok=True)
    res.to_csv(OUT_CSV, index=False)

    t = res[res.tested]
    print(f"\n{len(res):,} protein entries, {res.gene_symbol.nunique():,} symbols, "
          f"{len(t):,} tested")
    print(f"FDR < {FDR_CUT}: {int((t.padj_BH < FDR_CUT).sum()):,} "
          f"({int(((t.padj_BH < FDR_CUT) & (t.log2FC_tumor_vs_normal > 0)).sum()):,} up, "
          f"{int(((t.padj_BH < FDR_CUT) & (t.log2FC_tumor_vs_normal < 0)).sum()):,} down)")
    print(f"wrote {OUT_CSV}")


if __name__ == "__main__":
    main()
