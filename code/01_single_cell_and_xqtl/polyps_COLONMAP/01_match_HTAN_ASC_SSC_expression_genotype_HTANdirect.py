#!/usr/bin/env python3
"""
08_match_HTAN_ASC_SSC_expression_genotype_HTANdirect.py

Prepare HTAN ASC/SSC disease-specific pseudobulk count matrices for PEER/TWAS
WITHOUT restricting samples to the legacy N103 expression/genotype set.

The authoritative genotype set is the actual processed HTAN_GSE266546 genotype
PGEN/PSAM. HTAN donor -> genotype sample mapping comes from
VUMC_HTAN/575samples_meta_20240904.csv (STUDYID_HTAN -> VCFHead).

ASC: expected Sample_Classification = AD
SSC: expected Sample_Classification = SER

For donors with multiple disease specimens, retain one specimen using the same
segment priority used by the established NL workflow.

For donors with multiple VCFHead technical/plate genotypes:
  1) restrict candidates to samples actually present in the genotype PSAM;
  2) if >1 remain, prefer the sample with the lowest mean PLINK F_MISS across
     available chromosome post-QC .smiss files;
  3) if missingness information is unavailable/tied, use a deterministic
     lexicographic tie-break and record it in the audit.

Age/Gender/SEGMENT come from the CURRENT disease specimen in HTAN_metadata.csv.
Genotype PC1-PC5 come from the genome-wide PCA generated from the actual
processed genotype set.

Outputs per cell type:
  <expression_root>/<CT>/pseudobulk_counts_<CT>.csv
  <expression_root>/<CT>/<CT>.age_sex_tissue_geno5PCs.csv
  <expression_root>/<CT>/matching_<CT>/*audit files

Example:
  python 08_match_HTAN_ASC_SSC_expression_genotype_HTANdirect.py ASC \
    --count-matrix /path/count_matrix_batch12_ASC.chr1_22.csv
"""

from __future__ import annotations

import argparse
import glob
import re
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path("/v5000/data/l2_bioinfo1/liq17/scTWAS/data/VUMC_HTAN")
EXPRESSION_ROOT = ROOT / "expression" / "ADSER_pseudobulk_expression_for_peer"
HTAN_META_DIR = ROOT / "genotype"
GENO_ROOT = ROOT / "genotype" 

DEFAULT_HTAN_META = HTAN_META_DIR / "575samples_meta_20240904.csv"
DEFAULT_HTAN_SAMPLE_META = HTAN_META_DIR / "HTAN_metadata.csv"
DEFAULT_GENOTYPE_PSAM = (
    GENO_ROOT / "05_pgen_rsid" / "chr1.HTAN_GSE266546.qc.rsid.psam"
)
DEFAULT_PCA_EIGENVEC = (
    GENO_ROOT / "09_jointPCA" / "HTAN_GSE266546.chr1_22.PCA.eigenvec"
)
DEFAULT_SMISS_GLOB = str(GENO_ROOT / "qc" / "chr*.HTAN_GSE266546.postQC_stats.smiss")

CELL_TYPES = ("ASC", "SSC", "ABS", "EE", "GOB", "CT", "TUF", "TAC", "STM")
DEFAULT_CLASSIFICATIONS = ("AD", "SER")

SEGMENT_PRIORITY = [
    "Ascending colon",
    "Terminal ileum",
    "Cecum",
    "Descending colon",
    "Rectum",
    "Sigmoid colon",
    "Hepatic flexure",
    "Splenic flexure",
    "Transverse colon",
]
SEGMENT_PRIORITY_MAP = {
    x.casefold(): i for i, x in enumerate(SEGMENT_PRIORITY, start=1)
}

# Use Ascending colon as the fixed reference category when present. This is
# natural for the established colon workflow and keeps dummy naming stable.
SEGMENT_REFERENCE = "Ascending colon"


def clean_string(x):
    if pd.isna(x):
        return None
    x = str(x).strip()
    return x if x else None


def require_columns(df: pd.DataFrame, required, label: str):
    missing = sorted(set(required) - set(df.columns))
    if missing:
        raise ValueError(f"{label} missing column(s): {', '.join(missing)}")


def canonicalize_htan_expression_id(expression_id: str):
    """Normalize batch-specific HTAN specimen spelling to HTAN_metadata SPECID.htan.

    Examples:
      HTA11_10034_200000_1011 -> HTA11_10034_2000001011
      HTA11_104_2000001011    -> unchanged
    """
    x = str(expression_id).strip()
    x = re.sub(r"_(\d{6})_(\d+)$", r"_\1\2", x)
    return x


def donor_from_expression(expression_id: str):
    m = re.match(r"^(HTA11_[^_]+)_", str(expression_id).strip(), flags=re.I)
    return m.group(1) if m else None


def parse_sample_columns(columns, cell_type: str):
    prefix = f"{cell_type}-"
    rows = []
    for order, col in enumerate(columns):
        c = str(col)
        if not c.upper().startswith(prefix.upper()):
            continue

        m = re.match(
            rf"^{re.escape(cell_type)}-(HTA11_.+)-([A-Za-z0-9_]+)$",
            c,
            flags=re.I,
        )
        if not m:
            raise ValueError(
                f"Could not parse {cell_type} sample column: {c}\n"
                f"Expected {cell_type}-HTA11_...-<classification>."
            )

        raw_expression_id = m.group(1)
        expression_id = canonicalize_htan_expression_id(raw_expression_id)
        rows.append({
            "raw_column": c,
            "raw_expression_id": raw_expression_id,
            "expression_id": expression_id,
            "id_canonicalized": raw_expression_id != expression_id,
            "suffix_class": m.group(2).upper(),
            "_expression_order": order,
            "donor_id": donor_from_expression(expression_id),
        })

    out = pd.DataFrame(rows)
    if out.empty:
        raise RuntimeError(f"No {cell_type}-... sample columns found.")

    if out["expression_id"].duplicated().any():
        bad = out.loc[
            out["expression_id"].duplicated(keep=False),
            ["raw_column", "expression_id"],
        ]
        raise RuntimeError(
            "Duplicated bare expression IDs after stripping prefix/suffix:\n"
            + bad.to_string(index=False)
        )

    if out["donor_id"].isna().any():
        bad = out.loc[out["donor_id"].isna(), "expression_id"].tolist()
        raise RuntimeError("Could not derive HTAN donor IDs:\n" + "\n".join(bad))

    return out


def load_htan_sample_metadata(path: Path):
    meta = pd.read_csv(path)
    require_columns(
        meta,
        ["SPECID.htan", "STUDYID.htan", "Sample_Classification", "Age", "Gender", "SEGMENT"],
        "HTAN sample metadata",
    )
    meta = meta[
        ["SPECID.htan", "STUDYID.htan", "Sample_Classification", "Age", "Gender", "SEGMENT"]
    ].copy()

    for c in ["SPECID.htan", "STUDYID.htan", "Sample_Classification", "Gender", "SEGMENT"]:
        meta[c] = meta[c].map(clean_string)

    # Exact specimen annotations should be unique/consistent.
    check = (
        meta.dropna(subset=["SPECID.htan"])
        .assign(
            _combo=lambda d: (
                d["STUDYID.htan"].fillna("<NA>").astype(str)
                + "||" + d["Sample_Classification"].fillna("<NA>").astype(str)
                + "||" + d["Age"].fillna(-999999).astype(str)
                + "||" + d["Gender"].fillna("<NA>").astype(str)
                + "||" + d["SEGMENT"].fillna("<NA>").astype(str)
            )
        )
        .groupby("SPECID.htan")["_combo"]
        .nunique(dropna=False)
    )
    conflict = check[check > 1]
    if not conflict.empty:
        bad = meta[meta["SPECID.htan"].isin(conflict.index)].drop_duplicates()
        raise RuntimeError(
            "Conflicting HTAN annotations for SPECID.htan:\n" + bad.to_string(index=False)
        )

    return (
        meta.drop_duplicates("SPECID.htan")
        .rename(columns={"SPECID.htan": "expression_id", "STUDYID.htan": "metadata_donor_id"})
    )


def read_psam_ids(path: Path):
    if not path.exists():
        raise FileNotFoundError(f"Genotype PSAM not found: {path}")

    psam = pd.read_csv(path, sep=r"\s+", dtype=str)
    iid_col = None
    for c in ["IID", "#IID", "ID", "#ID"]:
        if c in psam.columns:
            iid_col = c
            break
    if iid_col is None:
        # PLINK2 PSAM normally has #FID IID; provide a robust fallback.
        if psam.shape[1] >= 2:
            iid_col = psam.columns[1]
        else:
            iid_col = psam.columns[0]

    ids = psam[iid_col].astype(str).str.strip()
    if ids.duplicated().any():
        dup = ids[ids.duplicated(keep=False)].drop_duplicates().tolist()
        raise RuntimeError(f"Duplicated genotype IDs in PSAM. Examples: {dup[:10]}")

    return set(ids), iid_col


def load_mean_fmiss(smiss_glob: str):
    files = sorted(glob.glob(smiss_glob))
    if not files:
        print(f"NOTE: no .smiss files matched: {smiss_glob}")
        print("      duplicate VCFHead ties will use deterministic genotype-ID order.")
        return pd.Series(dtype=float), []

    pieces = []
    used = []
    for f in files:
        try:
            x = pd.read_csv(f, sep=r"\s+", dtype={"IID": str, "#IID": str})
        except Exception as exc:
            print(f"WARNING: could not read {f}: {exc}")
            continue

        iid_col = "IID" if "IID" in x.columns else ("#IID" if "#IID" in x.columns else None)
        if iid_col is None or "F_MISS" not in x.columns:
            print(f"WARNING: {f} lacks IID/#IID or F_MISS; skipped")
            continue

        y = x[[iid_col, "F_MISS"]].copy()
        y.columns = ["genotype_id", "F_MISS"]
        y["genotype_id"] = y["genotype_id"].astype(str)
        y["F_MISS"] = pd.to_numeric(y["F_MISS"], errors="coerce")
        pieces.append(y)
        used.append(f)

    if not pieces:
        return pd.Series(dtype=float), []

    all_miss = pd.concat(pieces, ignore_index=True)
    mean_fmiss = all_miss.groupby("genotype_id")["F_MISS"].mean()
    return mean_fmiss, used


def load_genotype_map(path: Path, relevant_donors, genotype_ids, mean_fmiss):
    meta = pd.read_csv(path)
    require_columns(meta, ["STUDYID_HTAN", "VCFHead"], "HTAN genotype metadata")
    meta = meta[["STUDYID_HTAN", "VCFHead"]].copy()
    meta["STUDYID_HTAN"] = meta["STUDYID_HTAN"].map(clean_string)
    meta["VCFHead"] = meta["VCFHead"].map(clean_string)
    meta = meta[
        meta["STUDYID_HTAN"].isin(set(relevant_donors)) & meta["VCFHead"].notna()
    ].drop_duplicates()

    resolved = []
    audit = []

    for donor in sorted(set(relevant_donors)):
        grp = meta[meta["STUDYID_HTAN"] == donor]
        all_vcfs = grp["VCFHead"].dropna().astype(str).drop_duplicates().tolist()
        present = [v for v in all_vcfs if v in genotype_ids]

        chosen = None
        resolution = None
        chosen_fmiss = np.nan

        if len(present) == 0:
            resolution = "no_VCFHead_present_in_actual_genotype_PSAM"
        elif len(present) == 1:
            chosen = present[0]
            chosen_fmiss = float(mean_fmiss.get(chosen, np.nan))
            resolution = "single_VCFHead_in_actual_genotype_PSAM"
        else:
            ranked = []
            for v in present:
                fm = mean_fmiss.get(v, np.nan)
                # Missing F_MISS sorts after observed values; genotype ID provides
                # deterministic final tie-break.
                rank_fm = float(fm) if pd.notna(fm) else float("inf")
                ranked.append((rank_fm, v, fm))
            ranked.sort(key=lambda z: (z[0], z[1]))
            _, chosen, fm = ranked[0]
            chosen_fmiss = float(fm) if pd.notna(fm) else np.nan
            if any(pd.notna(z[2]) for z in ranked):
                resolution = "multiple_VCFHead_resolved_by_lowest_mean_F_MISS"
            else:
                resolution = "multiple_VCFHead_resolved_by_genotype_ID_tiebreak_no_F_MISS"
            print(
                f"NOTE: {donor} has {len(present)} genotype samples in PSAM; "
                f"selected {chosen} ({resolution})."
            )

        audit.append({
            "donor_id": donor,
            "n_vcfhead_metadata": len(all_vcfs),
            "vcfhead_all": ";".join(all_vcfs),
            "n_vcfhead_in_genotype": len(present),
            "vcfhead_in_genotype": ";".join(present),
            "genotype_id": chosen,
            "chosen_mean_F_MISS": chosen_fmiss,
            "genotype_resolution": resolution,
        })

        if chosen is not None:
            resolved.append({
                "donor_id": donor,
                "genotype_id": chosen,
                "chosen_mean_F_MISS": chosen_fmiss,
                "genotype_resolution": resolution,
            })

    return pd.DataFrame(resolved), pd.DataFrame(audit)


def load_pca(path: Path, n_pcs=5):
    if not path.exists():
        raise FileNotFoundError(
            f"Genome-wide PCA eigenvec not found: {path}\n"
            "Run the genotype PCA step first or supply --pca-eigenvec."
        )

    pca = pd.read_csv(path, sep=r"\s+", dtype=str)
    iid_col = None
    for c in ["IID", "#IID"]:
        if c in pca.columns:
            iid_col = c
            break
    if iid_col is None:
        if pca.shape[1] < n_pcs + 1:
            raise ValueError(f"Could not identify IID column in PCA file: {path}")
        # PLINK eigenvec without header can be FID IID PC1...; however current
        # PLINK2 output normally has a header. Fail rather than silently guess if
        # no PC labels are available.
        raise ValueError(
            f"PCA file has no IID/#IID header: {path}. "
            "Please use PLINK2 .eigenvec output with header."
        )

    pc_cols = [f"PC{i}" for i in range(1, n_pcs + 1)]
    require_columns(pca, [iid_col, *pc_cols], "PCA eigenvec")
    out = pca[[iid_col, *pc_cols]].copy()
    out = out.rename(columns={iid_col: "genotype_id"})
    out["genotype_id"] = out["genotype_id"].astype(str)
    out[pc_cols] = out[pc_cols].apply(pd.to_numeric, errors="raise")

    if out["genotype_id"].duplicated().any():
        raise RuntimeError("Duplicated genotype_id values in PCA eigenvec.")
    return out


def encode_sex(gender: pd.Series):
    norm = gender.fillna("").astype(str).str.strip().str.casefold()
    mapping = {
        "female": 0,
        "f": 0,
        "male": 1,
        "m": 1,
    }
    out = norm.map(mapping)
    if out.isna().any():
        bad = sorted(gender[out.isna()].astype(str).unique().tolist())
        raise ValueError(f"Unrecognized Gender values: {bad}")
    return out.astype(int)


def add_segment_dummies(df: pd.DataFrame):
    """Add full segment dummy set with Ascending colon as reference when present."""
    segments = df["SEGMENT"].fillna("").astype(str).str.strip()
    if (segments == "").any():
        bad = df.loc[segments == "", ["expression_id", "donor_id"]]
        raise RuntimeError("Missing SEGMENT in selected disease specimens:\n" + bad.to_string(index=False))

    observed = sorted(segments.unique().tolist())
    reference = SEGMENT_REFERENCE if SEGMENT_REFERENCE in observed else observed[0]
    dummy_cols = []
    for seg in observed:
        if seg == reference:
            continue
        col = f"segment_{seg}"
        df[col] = segments.eq(seg).astype(int)
        dummy_cols.append(col)
    return df, reference, dummy_cols


def set_exclusion_reasons(audit: pd.DataFrame, target_classes):
    target_label = ",".join(target_classes)
    audit["exclusion_reason"] = ""

    wrong_suffix = ~audit["suffix_class"].str.upper().isin(target_classes)
    no_meta = ~audit["specid_in_metadata"]
    donor_mismatch = audit["specid_in_metadata"] & ~audit["metadata_donor_match"]
    wrong_meta_class = audit["specid_in_metadata"] & ~audit["classification_match"]
    no_segment = (
        audit["specid_in_metadata"]
        & audit["classification_match"]
        & ~audit["segment_in_priority"]
    )
    no_gt = (
        audit["specid_in_metadata"]
        & audit["classification_match"]
        & audit["segment_in_priority"]
        & ~audit["has_genotype"]
    )
    donor_not_selected = audit["selection_candidate"] & ~audit["selected_for_donor"]
    no_pca = audit["selected_for_donor"] & ~audit["has_pca"]

    audit.loc[wrong_suffix, "exclusion_reason"] = (
        "sample_suffix_class_not_expected:" + audit.loc[wrong_suffix, "suffix_class"].astype(str)
    )
    audit.loc[no_meta, "exclusion_reason"] = "HTAN_SPECID_not_found_in_metadata"
    audit.loc[donor_mismatch, "exclusion_reason"] = "expression_donor_mismatch_with_HTAN_metadata"
    audit.loc[wrong_meta_class, "exclusion_reason"] = (
        "HTAN_classification_not_" + target_label + ":"
        + audit.loc[wrong_meta_class, "Sample_Classification"].fillna("NA").astype(str)
    )
    audit.loc[no_segment, "exclusion_reason"] = (
        "HTAN_segment_not_in_priority:"
        + audit.loc[no_segment, "SEGMENT"].fillna("NA").astype(str)
    )
    audit.loc[no_gt, "exclusion_reason"] = "HTAN_donor_without_genotype_in_actual_PSAM"
    audit.loc[donor_not_selected, "exclusion_reason"] = (
        "HTAN_not_selected_by_segment_priority;selected="
        + audit.loc[donor_not_selected, "selected_expression_id_for_donor"].fillna("NA").astype(str)
    )
    audit.loc[no_pca, "exclusion_reason"] = "selected_genotype_missing_from_PCA_eigenvec"
    return audit


def main():
    p = argparse.ArgumentParser()
    p.add_argument("cell_type", choices=CELL_TYPES)
    p.add_argument(
        "--count-matrix", default=None,
        help=("Default: <expression-root>/<CT>/"
              "count_matrix_batch12_<CT>.chr1_22.csv"),
    )
    p.add_argument(
        "--classifications", nargs="+", default=list(DEFAULT_CLASSIFICATIONS),
        help="Accepted Sample_Classification/suffix values (default: AD SER)",
    )
    p.add_argument("--htan-meta", default=str(DEFAULT_HTAN_META))
    p.add_argument("--htan-sample-meta", default=str(DEFAULT_HTAN_SAMPLE_META))
    p.add_argument("--genotype-psam", default=str(DEFAULT_GENOTYPE_PSAM), help="Actual genotype PSAM used to define available genotype samples")
    p.add_argument("--pca-eigenvec", default=str(DEFAULT_PCA_EIGENVEC), help="Genome-wide genotype PCA eigenvec")
    p.add_argument("--smiss-glob", default=DEFAULT_SMISS_GLOB, help="Glob for per-chromosome PLINK sample missingness files")
    p.add_argument("--expression-root", default=str(EXPRESSION_ROOT))
    args = p.parse_args()

    ct = args.cell_type.upper()
    target_classes = tuple(dict.fromkeys(x.upper() for x in args.classifications))
    if not target_classes:
        p.error("--classifications must contain at least one value")
    cell_dir = Path(args.expression_root) / ct
    count_file = (Path(args.count_matrix) if args.count_matrix else
                  cell_dir / f"count_matrix_batch12_{ct}.chr1_22.csv")
    match_dir = cell_dir / f"matching_{ct}"
    cell_dir.mkdir(parents=True, exist_ok=True)
    match_dir.mkdir(parents=True, exist_ok=True)

    print("=" * 78)
    print(f"HTAN {ct}: direct HTAN genotype matching (NO N103 restriction)")
    print("=" * 78)
    print(f"Count matrix       : {count_file}")
    print(f"Accepted classes   : {','.join(target_classes)}")
    print(f"Genotype PSAM      : {args.genotype_psam}")
    print(f"PCA eigenvec       : {args.pca_eigenvec}")
    print(f"Cell output dir    : {cell_dir}")

    header = pd.read_csv(count_file, nrows=0)
    samples = parse_sample_columns(header.columns, ct)
    print(f"Detected specimens : {len(samples)}")
    print(f"Unique donors      : {samples['donor_id'].nunique()}")
    print(f"Suffix classes     : {samples['suffix_class'].value_counts().to_dict()}")

    sample_meta = load_htan_sample_metadata(Path(args.htan_sample_meta))
    audit = samples.merge(
        sample_meta,
        on="expression_id",
        how="left",
        validate="one_to_one",
        indicator="_meta_merge",
    )
    audit["specid_in_metadata"] = audit["_meta_merge"].eq("both")
    audit = audit.drop(columns="_meta_merge")
    audit["metadata_donor_match"] = (
        audit["metadata_donor_id"].fillna("").astype(str).eq(audit["donor_id"].fillna("").astype(str))
    )
    audit["classification_match"] = (
        audit["Sample_Classification"].fillna("").astype(str).str.strip().str.upper().isin(target_classes)
    )
    audit["segment_priority"] = (
        audit["SEGMENT"].fillna("").astype(str).str.strip().str.casefold().map(SEGMENT_PRIORITY_MAP)
    )
    audit["segment_in_priority"] = audit["segment_priority"].notna()

    genotype_ids, psam_iid_col = read_psam_ids(Path(args.genotype_psam))
    print(f"Actual genotype IDs : {len(genotype_ids)} (PSAM column {psam_iid_col})")

    mean_fmiss, smiss_files = load_mean_fmiss(args.smiss_glob)
    print(f"Missingness files   : {len(smiss_files)}")

    geno, genotype_resolution = load_genotype_map(
        Path(args.htan_meta),
        audit["donor_id"].dropna().astype(str).unique().tolist(),
        genotype_ids,
        mean_fmiss,
    )
    audit = audit.merge(geno, on="donor_id", how="left", validate="many_to_one")
    audit["has_genotype"] = audit["genotype_id"].notna()

    pca = load_pca(Path(args.pca_eigenvec), n_pcs=5)
    pca_ids = set(pca["genotype_id"])
    audit["has_pca"] = audit["genotype_id"].astype("string").isin(pca_ids)

    audit["selection_candidate"] = (
        audit["suffix_class"].str.upper().isin(target_classes)
        & audit["specid_in_metadata"]
        & audit["metadata_donor_match"]
        & audit["classification_match"]
        & audit["segment_in_priority"]
        & audit["has_genotype"]
    )
    audit["selected_for_donor"] = False

    candidates = audit[audit["selection_candidate"]].copy()
    if not candidates.empty:
        selected_idx = (
            candidates.sort_values(
                ["donor_id", "segment_priority", "_expression_order"],
                kind="stable",
            )
            .drop_duplicates("donor_id", keep="first")
            .index
        )
        audit.loc[selected_idx, "selected_for_donor"] = True

    selected_rows = audit[audit["selected_for_donor"]][
        ["donor_id", "expression_id", "SEGMENT", "segment_priority"]
    ]
    audit["selected_expression_id_for_donor"] = audit["donor_id"].map(
        dict(zip(selected_rows["donor_id"], selected_rows["expression_id"]))
    )
    audit["selected_segment_for_donor"] = audit["donor_id"].map(
        dict(zip(selected_rows["donor_id"], selected_rows["SEGMENT"]))
    )

    audit["eligible_for_peer"] = audit["selected_for_donor"] & audit["has_pca"]
    audit = set_exclusion_reasons(audit, target_classes)
    final = audit[audit["eligible_for_peer"]].copy()

    if final.empty:
        raise RuntimeError(
            f"No {ct} samples remain after disease-class, donor, actual-genotype, and PCA matching."
        )

    # ------------------------------------------------------------------
    # Build disease-specimen-specific PEER covariates.
    # ------------------------------------------------------------------
    final_cov = final[
        [
            "expression_id", "genotype_id", "donor_id", "Sample_Classification",
            "Age", "Gender", "SEGMENT",
        ]
    ].merge(pca, on="genotype_id", how="left", validate="one_to_one")

    final_cov["Age"] = pd.to_numeric(final_cov["Age"], errors="raise")
    final_cov["Sex"] = encode_sex(final_cov["Gender"])
    final_cov, segment_reference, segment_dummy_cols = add_segment_dummies(final_cov)

    numeric_cols = ["Age", "Sex", *segment_dummy_cols, "PC1", "PC2", "PC3", "PC4", "PC5"]
    final_cov[numeric_cols] = final_cov[numeric_cols].apply(pd.to_numeric, errors="raise")
    if final_cov[numeric_cols].isna().any().any():
        bad = final_cov.columns[final_cov.isna().any()].tolist()
        raise RuntimeError(f"Missing values in final {ct} PEER covariates: {bad}")

    cov_out_cols = [
        "expression_id", "genotype_id", "donor_id", "Sample_Classification",
        "Age", "Gender", "Sex", "SEGMENT", *segment_dummy_cols,
        "PC1", "PC2", "PC3", "PC4", "PC5",
    ]
    cov_file = cell_dir / f"{ct}.age_sex_tissue_geno5PCs.csv"
    final_cov[cov_out_cols].to_csv(cov_file, index=False)

    # Save PEER covariate column manifest so the PEER script can audit the exact
    # model matrix while still auto-detecting segment_* columns.
    manifest = pd.DataFrame({
        "covariate": ["Age", "Sex", *segment_dummy_cols, "PC1", "PC2", "PC3", "PC4", "PC5"]
    })
    manifest.to_csv(cell_dir / f"{ct}.PEER_covariate_columns.txt", index=False, header=False)

    # ------------------------------------------------------------------
    # Clean raw count matrix.
    # ------------------------------------------------------------------
    raw = pd.read_csv(count_file)
    gene_col = "geneid" if "geneid" in raw.columns else raw.columns[0]
    keep_raw_cols = final["raw_column"].tolist()

    counts = raw[[gene_col, *keep_raw_cols]].copy()
    counts[keep_raw_cols] = counts[keep_raw_cols].apply(pd.to_numeric, errors="coerce")
    n_missing = int(counts[keep_raw_cols].isna().sum().sum())
    if n_missing:
        print(f"Replacing missing sample-count entries with 0: {n_missing:,}")
        counts[keep_raw_cols] = counts[keep_raw_cols].fillna(0)

    arr = counts[keep_raw_cols].to_numpy(dtype=float)
    if not np.isfinite(arr).all():
        raise RuntimeError("Non-finite counts remain after NaN handling.")
    if (arr < 0).any():
        raise RuntimeError("Negative counts detected.")

    rename_map = dict(zip(final["raw_column"], final["expression_id"]))
    counts = counts.rename(columns={gene_col: "gene", **rename_map})
    if counts["gene"].duplicated().any():
        counts = counts.groupby("gene", as_index=False).sum(numeric_only=True)

    raw_out = cell_dir / f"pseudobulk_counts_{ct}.csv"
    counts.to_csv(raw_out, index=False)

    # ------------------------------------------------------------------
    # Audits.
    # ------------------------------------------------------------------
    audit.drop(columns="_expression_order").to_csv(
        match_dir / f"01_{ct}_expression_genotype_audit.csv", index=False
    )
    genotype_resolution.to_csv(
        match_dir / f"00_{ct}_donor_genotype_resolution_audit.csv", index=False
    )
    final.drop(columns="_expression_order").to_csv(
        match_dir / f"02_{ct}_expression_genotype_matched.csv", index=False
    )
    final[["raw_column", "expression_id", "donor_id", "genotype_id", "SEGMENT"]].to_csv(
        match_dir / f"03_{ct}_sample_rename_map.csv", index=False
    )
    final[["expression_id"]].to_csv(
        match_dir / f"04_{ct}_expression_samples_to_keep.txt", index=False, header=False
    )
    final[["genotype_id"]].drop_duplicates().to_csv(
        match_dir / f"05_{ct}_genotype_samples_to_keep.txt", index=False, header=False
    )

    summary = pd.DataFrame([{
        "cell_type": ct,
        "accepted_classifications": ";".join(target_classes),
        "input_specimens": len(audit),
        "input_unique_donors": audit["donor_id"].nunique(),
        "exact_metadata_matches": int(audit["specid_in_metadata"].sum()),
        "classification_matches": int(audit["classification_match"].sum()),
        "specimens_with_actual_genotype": int(audit["has_genotype"].sum()),
        "unique_donors_with_actual_genotype": audit.loc[audit["has_genotype"], "donor_id"].nunique(),
        "donor_selected": int(audit["selected_for_donor"].sum()),
        "selected_with_PCA": int(audit["eligible_for_peer"].sum()),
        "final_PEER_samples": len(final),
        "segment_reference": segment_reference,
        "segment_dummy_columns": ";".join(segment_dummy_cols),
        "missing_count_entries_filled_zero": n_missing,
        "genotype_psam": str(args.genotype_psam),
        "pca_eigenvec": str(args.pca_eigenvec),
        "raw_pseudobulk_file": str(raw_out),
        "cell_specific_covariate_file": str(cov_file),
    }])
    summary.to_csv(match_dir / f"06_{ct}_matching_summary.csv", index=False)

    print()
    print("=" * 78)
    print("DONE")
    print("=" * 78)
    print(f"Final PEER samples : {len(final)}")
    print(f"Raw pseudobulk     : {raw_out}")
    print(f"Covariates         : {cov_file}")
    print(f"Segment reference  : {segment_reference}")
    print(f"Segment dummies    : {segment_dummy_cols}")
    print(f"Audit directory    : {match_dir}")
    print()
    print("Next:")
    print(f"  python 08_prepare_pseudobulk_expression_ASC_SSC.py --cell-types {ct}")
    print(f"  python Process_peer_HTAN_ASC_SSC_HTANdirect.py {ct}")


if __name__ == "__main__":
    main()
