#!/usr/bin/env python3
"""
08_match_N157_expression_genotype_segment_priority_NL.py

Goal
----
1. Determine how many samples in N157_Expr_head.csv have genotype data.
2. Generate one master ID-matching table between N157 expression sample IDs
   and exact genotype/VCF sample IDs for BOTH:
       - HTAN
       - GSE266546

Inputs
------
A. N157 expression header:
   N157_Expr_head.csv

   Expected:
       gene,<sample1>,<sample2>,...

   HTAN expression IDs look like:
       HTA11_10034_2000002011
       HTA11_104_2000002021

   GSE266546 expression IDs look like:
       5222-YX-10
       5314-YX-24


B. HTAN genotype metadata:
   575samples_meta_20240904.csv

   Required columns:
       STUDYID_HTAN
       VCFHead

   Genotype matching:
       HTA11_10034_2000002011
           -> HTA11_10034
           -> VCFHead


C. HTAN expression phenotype metadata:
   HTAN_metadata.csv

   Required columns:
       SPECID.htan
       Sample_Classification
       SEGMENT

   HTAN inclusion rule:
       N157 expression_id must EXACTLY match SPECID.htan
       AND Sample_Classification == "NL".

       If a donor has multiple eligible NL expression specimens, retain ONE
       specimen using this SEGMENT priority (highest to lowest):
           1. Ascending colon
           2. Terminal ileum
           3. Cecum
           4. Descending colon
           5. Rectum
           6. Sigmoid colon
           7. Hepatic flexure
           8. Splenic flexure
           9. Transverse colon

       If two specimens from the same donor have the same priority segment,
       the one appearing first in the N157 expression header is retained.

   Examples:
       donor with Ascending colon + Rectum -> keep Ascending colon
       donor with Terminal ileum + Cecum   -> keep Terminal ileum
       donor with one NL Rectum specimen   -> keep Rectum
       AD/SER specimen                     -> excluded
       no exact SPECID.htan match          -> excluded


D. GSE266546 genotype/scRNA mapping:


D. GSE266546 genotype/scRNA mapping:
   GDA_genotype_sc_exp_sample_matched.csv

   Required columns:
       IID
       tissue_id

   Example:
       IID = GCA017
       tissue_id = GSM8249577


E. GSE266546 VCF-head sample list:
   GSE266546_geno_N144_samples.txt

   One exact VCF sample ID per line, e.g.:
       16_GCA017DNA
       50_GCA054
       68_gca072dna

   Matching is done through a canonical GCA donor ID, while the ORIGINAL
   VCF sample name is preserved for downstream bcftools/plink use.


F. GSE266546 GEO sample-description mapping

   The expression IDs such as "5314-YX-24" are not directly present in the
   supplied GDA mapping table. The bridge is:

       tissue_id (GSM...)
           -> GEO Sample_description
           -> N157 expression ID

   The script uses a local GEO cache if available. If the cache does not
   exist, it queries NCBI GEO and creates the cache automatically.

   Expected cache columns:
       tissue_id
       geo_description
       expression_id

Outputs
-------
01_N157_expression_genotype_audit.csv
    One row per N157 expression sample, including matched/unmatched status.

02_N157_expression_genotype_matched.csv
    Final analysis set:
      - HTAN: exact SPECID.htan match + Sample_Classification == NL
              + one specimen per donor selected by SEGMENT priority + genotype
      - GSE266546: genotype matched

03_HTAN_expression_genotype_matched.csv
04_GSE266546_expression_genotype_matched.csv

05_HTAN_vcf_samples_to_keep.txt
06_GSE266546_vcf_samples_to_keep.txt

07_N157_expression_samples_with_genotype.txt
08_N157_expression_excluded_from_analysis.csv

09_N157_expression_genotype_summary.csv

10_GSE266546_GEO_sample_description_cache.csv

Important
---------
HTAN expression specimens are not necessarily independent donors.
This version resolves that issue during matching by retaining at most ONE
NL HTAN expression specimen per donor according to the SEGMENT priority rule.
GSE266546 matching behavior is unchanged.
"""

from __future__ import annotations

import argparse
import re
import sys
import time
from pathlib import Path
from urllib.error import HTTPError, URLError
from urllib.request import Request, urlopen

import pandas as pd


# =============================================================================
# Generic helpers
# =============================================================================

def clean_string(x):
    if pd.isna(x):
        return None
    x = str(x).strip()
    return x if x else None


def require_columns(df, required, label):
    missing = sorted(set(required) - set(df.columns))
    if missing:
        raise ValueError(
            f"{label} is missing required column(s): {', '.join(missing)}"
        )


def read_expression_samples(expression_head):
    """
    Read only the CSV header.

    First column is assumed to be gene/gene_id.
    Remaining columns are expression sample IDs.
    """
    cols = list(pd.read_csv(expression_head, nrows=0).columns)

    if len(cols) < 2:
        raise ValueError(
            f"No expression sample columns found in {expression_head}"
        )

    samples = [str(x).strip() for x in cols[1:]]

    if len(samples) != len(set(samples)):
        dup = pd.Series(samples)
        dup = dup[dup.duplicated(keep=False)].unique().tolist()
        raise ValueError(
            "Duplicated expression sample IDs in N157 header:\n"
            + "\n".join(map(str, dup))
        )

    return samples


# =============================================================================
# HTAN helpers
# =============================================================================

# For HTAN donors with multiple NL expression specimens, keep one specimen
# according to this priority (1 = highest priority).
HTAN_SEGMENT_PRIORITY = [
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

HTAN_SEGMENT_PRIORITY_MAP = {
    segment.casefold(): rank
    for rank, segment in enumerate(HTAN_SEGMENT_PRIORITY, start=1)
}

def htan_donor_from_expression(expression_id):
    """
    HTA11_10034_2000002011 -> HTA11_10034
    """
    x = str(expression_id).strip()

    m = re.match(r"^(HTA11_[^_]+)_", x, flags=re.IGNORECASE)

    if not m:
        return None

    return m.group(1)


def match_htan(
    expression_ids,
    htan_genotype_meta_file,
    htan_sample_meta_file,
):
    """
    Match HTAN N157 expression samples and apply the donor-level selection rule.

    Step 1: expression_id EXACTLY matches HTAN_metadata.csv SPECID.htan.
    Step 2: require Sample_Classification == "NL".
    Step 3: among NL samples with a recognized SEGMENT, retain at most ONE
            expression sample per donor according to HTAN_SEGMENT_PRIORITY.
    Step 4: if two samples from the same donor have the same SEGMENT priority,
            retain the one appearing first in the N157 expression header.
    Step 5: derive donor STUDYID_HTAN from expression_id and map to VCFHead
            using 575samples_meta_20240904.csv.
    Step 6: final eligibility requires the donor-selected sample + genotype.

    All N157 HTAN samples remain in the audit table, including samples excluded
    by classification, segment priority, donor selection, or genotype status.
    """

    expr = pd.DataFrame({
        "expression_id": expression_ids
    })

    # Stable tie-breaker when one donor has >1 sample in the same segment.
    expr["_expression_order"] = range(len(expr))

    expr["project"] = "HTAN"
    expr["donor_id"] = expr["expression_id"].map(
        htan_donor_from_expression
    )

    # -------------------------------------------------------------------------
    # HTAN phenotype metadata: exact expression_id -> SPECID.htan
    # -------------------------------------------------------------------------

    sample_meta = pd.read_csv(
        htan_sample_meta_file
    )

    require_columns(
        sample_meta,
        [
            "SPECID.htan",
            "Sample_Classification",
            "SEGMENT",
        ],
        "HTAN sample metadata",
    )

    sample_meta = sample_meta.copy()

    for col in [
        "SPECID.htan",
        "Sample_Classification",
        "SEGMENT",
    ]:
        sample_meta[col] = sample_meta[col].map(clean_string)

    # One SPECID must not have conflicting classification/segment annotations.
    sample_meta["_class_segment"] = (
        sample_meta["Sample_Classification"].fillna("<NA>").astype(str)
        + "||"
        + sample_meta["SEGMENT"].fillna("<NA>").astype(str)
    )

    annotation_check = (
        sample_meta[
            sample_meta["SPECID.htan"].notna()
        ]
        .groupby("SPECID.htan")["_class_segment"]
        .nunique(dropna=False)
    )

    conflicting = annotation_check[
        annotation_check > 1
    ]

    if not conflicting.empty:
        bad = (
            sample_meta[
                sample_meta["SPECID.htan"].isin(
                    conflicting.index
                )
            ][
                [
                    "SPECID.htan",
                    "Sample_Classification",
                    "SEGMENT",
                ]
            ]
            .drop_duplicates()
            .sort_values("SPECID.htan")
        )

        raise RuntimeError(
            "Conflicting HTAN Sample_Classification/SEGMENT values "
            "for the same SPECID.htan:\n"
            + bad.to_string(index=False)
        )

    sample_class = (
        sample_meta[
            [
                "SPECID.htan",
                "Sample_Classification",
                "SEGMENT",
            ]
        ]
        .drop_duplicates("SPECID.htan")
        .rename(
            columns={
                "SPECID.htan": "expression_id",
            }
        )
    )

    # Exact, case-sensitive SPECID.htan match.
    expr = expr.merge(
        sample_class,
        on="expression_id",
        how="left",
        validate="one_to_one",
        indicator="_specid_merge",
    )

    expr["htan_specid_in_metadata"] = (
        expr["_specid_merge"].eq("both")
    )
    expr = expr.drop(columns="_specid_merge")

    expr["htan_is_NL"] = (
        expr["Sample_Classification"]
        .fillna("")
        .astype(str)
        .str.strip()
        .str.upper()
        .eq("NL")
    )

    # Normalize SEGMENT only for matching; preserve original SEGMENT in output.
    segment_norm = (
        expr["SEGMENT"]
        .fillna("")
        .astype(str)
        .str.strip()
        .str.casefold()
    )

    expr["htan_segment_priority"] = (
        segment_norm.map(HTAN_SEGMENT_PRIORITY_MAP)
    )

    expr["htan_segment_in_priority"] = (
        expr["htan_segment_priority"].notna()
    )

    # -------------------------------------------------------------------------
    # HTAN genotype mapping: donor STUDYID_HTAN -> exact VCFHead
    # -------------------------------------------------------------------------

    genotype_meta = pd.read_csv(
        htan_genotype_meta_file
    )

    require_columns(
        genotype_meta,
        ["STUDYID_HTAN", "VCFHead"],
        "HTAN genotype metadata",
    )

    genotype_meta = genotype_meta.copy()

    genotype_meta["STUDYID_HTAN"] = (
        genotype_meta["STUDYID_HTAN"]
        .map(clean_string)
    )

    genotype_meta["VCFHead"] = (
        genotype_meta["VCFHead"]
        .map(clean_string)
    )

    relevant_donors = set(
        expr["donor_id"].dropna()
    )

    genotype_meta = genotype_meta[
        genotype_meta["STUDYID_HTAN"].isin(
            relevant_donors
        )
    ].copy()

    genotype_nonnull = (
        genotype_meta[
            genotype_meta["VCFHead"].notna()
        ][
            ["STUDYID_HTAN", "VCFHead"]
        ]
        .drop_duplicates()
    )

    # One donor must map to only one genotype sample.
    n_gt = (
        genotype_nonnull
        .groupby("STUDYID_HTAN")["VCFHead"]
        .nunique()
    )

    ambiguous = n_gt[
        n_gt > 1
    ]

    if not ambiguous.empty:
        detail = (
            genotype_nonnull[
                genotype_nonnull["STUDYID_HTAN"].isin(
                    ambiguous.index
                )
            ]
            .sort_values(
                [
                    "STUDYID_HTAN",
                    "VCFHead",
                ]
            )
        )

        raise RuntimeError(
            "Ambiguous HTAN genotype mapping: "
            "some N157 donors map to >1 VCFHead.\n"
            + detail.to_string(index=False)
        )

    donor_to_gt = (
        genotype_nonnull
        .drop_duplicates("STUDYID_HTAN")
        .rename(
            columns={
                "STUDYID_HTAN": "donor_id",
                "VCFHead": "genotype_id",
            }
        )
    )

    out = expr.merge(
        donor_to_gt,
        on="donor_id",
        how="left",
        validate="many_to_one",
    )

    out["tissue_id"] = pd.NA

    # Genotype availability is kept separate from phenotype eligibility.
    out["has_genotype"] = (
        out["genotype_id"].notna()
    )

    # -------------------------------------------------------------------------
    # Select ONE HTAN NL expression specimen per donor by SEGMENT priority.
    # -------------------------------------------------------------------------

    out["htan_selection_candidate"] = (
        out["htan_specid_in_metadata"]
        & out["htan_is_NL"]
        & out["htan_segment_in_priority"]
        & out["donor_id"].notna()
    )

    out["htan_selected_for_donor"] = False

    candidates = out[
        out["htan_selection_candidate"]
    ].copy()

    if not candidates.empty:
        selected_idx = (
            candidates
            .sort_values(
                [
                    "donor_id",
                    "htan_segment_priority",
                    "_expression_order",
                ],
                kind="stable",
            )
            .drop_duplicates(
                "donor_id",
                keep="first",
            )
            .index
        )

        out.loc[
            selected_idx,
            "htan_selected_for_donor",
        ] = True

    # Audit columns showing which specimen won for each donor.
    selected_rows = out[
        out["htan_selected_for_donor"]
    ][
        [
            "donor_id",
            "expression_id",
            "SEGMENT",
            "htan_segment_priority",
        ]
    ].copy()

    selected_expr_map = dict(
        zip(
            selected_rows["donor_id"],
            selected_rows["expression_id"],
        )
    )
    selected_segment_map = dict(
        zip(
            selected_rows["donor_id"],
            selected_rows["SEGMENT"],
        )
    )
    selected_priority_map = dict(
        zip(
            selected_rows["donor_id"],
            selected_rows["htan_segment_priority"],
        )
    )

    out["htan_selected_expression_id_for_donor"] = (
        out["donor_id"].map(selected_expr_map)
    )
    out["htan_selected_segment_for_donor"] = (
        out["donor_id"].map(selected_segment_map)
    )
    out["htan_selected_priority_for_donor"] = (
        out["donor_id"].map(selected_priority_map)
    )

    # FINAL HTAN inclusion rule:
    # exact SPECID.htan + NL + recognized segment + donor-selected specimen
    # + genotype.
    out["eligible_for_analysis"] = (
        out["htan_selected_for_donor"]
        & out["has_genotype"]
    )

    out["mapping_source"] = (
        "expression exact SPECID.htan + NL; one sample/donor selected by "
        "SEGMENT priority; donor STUDYID_HTAN -> VCFHead"
    )

    # -------------------------------------------------------------------------
    # Explicit exclusion reason for audit.
    # -------------------------------------------------------------------------

    out["exclusion_reason"] = ""

    no_specid = (
        ~out["htan_specid_in_metadata"]
    )

    non_nl = (
        out["htan_specid_in_metadata"]
        & ~out["htan_is_NL"]
    )

    segment_not_prioritized = (
        out["htan_specid_in_metadata"]
        & out["htan_is_NL"]
        & ~out["htan_segment_in_priority"]
    )

    not_selected = (
        out["htan_selection_candidate"]
        & ~out["htan_selected_for_donor"]
    )

    no_gt = (
        out["htan_selected_for_donor"]
        & ~out["has_genotype"]
    )

    out.loc[
        no_specid,
        "exclusion_reason",
    ] = "HTAN_SPECID_not_found_in_metadata"

    out.loc[
        non_nl,
        "exclusion_reason",
    ] = (
        "HTAN_not_NL:"
        + out.loc[
            non_nl,
            "Sample_Classification",
        ].fillna("NA").astype(str)
    )

    out.loc[
        segment_not_prioritized,
        "exclusion_reason",
    ] = (
        "HTAN_NL_segment_not_in_priority:"
        + out.loc[
            segment_not_prioritized,
            "SEGMENT",
        ].fillna("NA").astype(str)
    )

    out.loc[
        not_selected,
        "exclusion_reason",
    ] = (
        "HTAN_NL_not_selected_by_segment_priority;selected="
        + out.loc[
            not_selected,
            "htan_selected_expression_id_for_donor",
        ].fillna("NA").astype(str)
        + ";selected_segment="
        + out.loc[
            not_selected,
            "htan_selected_segment_for_donor",
        ].fillna("NA").astype(str)
    )

    out.loc[
        no_gt,
        "exclusion_reason",
    ] = "HTAN_selected_NL_sample_without_genotype"

    # The temporary order column was only needed to resolve ties.
    out = out.drop(columns="_expression_order")

    return out


# =============================================================================
# GSE266546 helpers
# =============================================================================

def canonical_gca(x):
    """
    Normalize GSE genotype donor IDs.

    Examples
    --------
    GCA017             -> GCA017
    16_GCA017DNA       -> GCA017
    68_gca072dna       -> GCA072
    283_GCA_125DNA     -> GCA125
    50_GCA054          -> GCA054
    """
    x = clean_string(x)

    if x is None:
        return None

    m = re.search(
        r"GCA_?0*(\d+)",
        x,
        flags=re.IGNORECASE,
    )

    if not m:
        return None

    return f"GCA{int(m.group(1)):03d}"


def read_exact_vcf_samples(sample_file):
    """
    Read exact genotype IDs from a VCF-header sample list.
    """
    with open(sample_file, "r") as fh:
        samples = [
            line.strip()
            for line in fh
            if line.strip()
        ]

    if not samples:
        raise ValueError(
            f"No genotype sample IDs found in {sample_file}"
        )

    if len(samples) != len(set(samples)):
        s = pd.Series(samples)
        dup = s[s.duplicated(keep=False)].unique().tolist()

        raise ValueError(
            "Duplicated exact VCF sample IDs:\n"
            + "\n".join(map(str, dup))
        )

    out = pd.DataFrame({
        "genotype_id": samples
    })

    out["donor_id"] = (
        out["genotype_id"]
        .map(canonical_gca)
    )

    if out["donor_id"].isna().any():
        bad = out.loc[
            out["donor_id"].isna(),
            "genotype_id",
        ].tolist()

        raise ValueError(
            "Could not extract GCA donor ID from VCF sample(s):\n"
            + "\n".join(bad)
        )

    # One donor must correspond to one exact VCF sample.
    duplicated_donor = out[
        out["donor_id"].duplicated(keep=False)
    ]

    if not duplicated_donor.empty:
        raise RuntimeError(
            "Multiple GSE VCF sample IDs map to the same GCA donor:\n"
            + duplicated_donor.to_string(index=False)
        )

    return out


def normalize_geo_expression_id(x):
    """
    Convert GEO sequencing description to the N157 sample ID.

    Examples
    --------
    5314-YX-24_S1_L005 -> 5314-YX-24
    5222-YX-11_S1_L005 -> 5222-YX-11
    """
    x = clean_string(x)

    if x is None:
        return None

    # Remove a standard Illumina sample/lane suffix.
    x = re.sub(
        r"_S\d+_L\d{3}.*$",
        "",
        x,
        flags=re.IGNORECASE,
    )

    # Fallback if a description contains _S1 but no lane token.
    x = re.sub(
        r"_S\d+.*$",
        "",
        x,
        flags=re.IGNORECASE,
    )

    return x


def parse_geo_description(text):
    """
    Parse !Sample_description from one GEO SOFT record.
    """
    values = []

    for line in text.splitlines():
        if line.startswith("!Sample_description"):
            try:
                _, value = line.split("=", 1)
            except ValueError:
                continue

            value = value.strip()

            if value:
                values.append(value)

    if not values:
        return None

    # Prefer the known GSE266546 expression-ID style.
    for value in values:
        if "-YX-" in value:
            return value

    return values[0]


def fetch_geo_description(gsm, retries=3):
    """
    Query one GSM record from NCBI GEO.
    """
    url = (
        "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi"
        f"?acc={gsm}&targ=self&form=text&view=full"
    )

    req = Request(
        url,
        headers={
            "User-Agent":
                "N157-expression-genotype-matcher/1.0"
        },
    )

    last_error = None

    for attempt in range(1, retries + 1):
        try:
            with urlopen(
                req,
                timeout=30,
            ) as response:
                text = response.read().decode(
                    "utf-8",
                    errors="replace",
                )

            description = parse_geo_description(text)

            if description is None:
                raise RuntimeError(
                    f"No GEO Sample_description found for {gsm}"
                )

            return description

        except (
            HTTPError,
            URLError,
            TimeoutError,
            RuntimeError,
        ) as exc:
            last_error = exc

            if attempt < retries:
                time.sleep(attempt)

    raise RuntimeError(
        f"Failed to retrieve GEO metadata for {gsm}: "
        f"{last_error}"
    )


def standardize_geo_cache(cache):
    """
    Accept a few reasonable cache-column names and standardize to:
        tissue_id
        geo_description
        expression_id
    """
    cache = cache.copy()

    require_columns(
        cache,
        ["tissue_id"],
        "GSE GEO cache",
    )

    if "geo_description" not in cache.columns:
        if "Sample_description" in cache.columns:
            cache = cache.rename(
                columns={
                    "Sample_description": "geo_description"
                }
            )
        elif "sample_description" in cache.columns:
            cache = cache.rename(
                columns={
                    "sample_description": "geo_description"
                }
            )

    if "expression_id" not in cache.columns:
        if "expression_sample" in cache.columns:
            cache = cache.rename(
                columns={
                    "expression_sample": "expression_id"
                }
            )
        elif "geo_description" in cache.columns:
            cache["expression_id"] = (
                cache["geo_description"]
                .map(normalize_geo_expression_id)
            )

    require_columns(
        cache,
        [
            "tissue_id",
            "expression_id",
        ],
        "GSE GEO cache",
    )

    if "geo_description" not in cache.columns:
        cache["geo_description"] = pd.NA

    cache["tissue_id"] = (
        cache["tissue_id"]
        .map(clean_string)
    )

    cache["expression_id"] = (
        cache["expression_id"]
        .map(normalize_geo_expression_id)
    )

    return (
        cache[
            [
                "tissue_id",
                "geo_description",
                "expression_id",
            ]
        ]
        .drop_duplicates()
    )


def build_or_update_geo_cache(
    tissue_ids,
    cache_file,
    allow_fetch=True,
):
    """
    Reuse local cache; fetch missing GSM descriptions only when needed.
    """
    cache_file = Path(cache_file)

    if cache_file.exists():
        cache = standardize_geo_cache(
            pd.read_csv(cache_file)
        )
    else:
        cache = pd.DataFrame(
            columns=[
                "tissue_id",
                "geo_description",
                "expression_id",
            ]
        )

    existing = set(
        cache.loc[
            cache["expression_id"].notna(),
            "tissue_id",
        ].astype(str)
    )

    requested = sorted(
        {
            str(x)
            for x in tissue_ids
            if clean_string(x) is not None
        }
    )

    missing = [
        gsm
        for gsm in requested
        if gsm not in existing
    ]

    if missing and not allow_fetch:
        raise RuntimeError(
            f"{len(missing)} GSM IDs are missing from GEO cache, "
            "and --no-geo-fetch was requested.\n"
            f"Cache: {cache_file}\n"
            "Missing examples: "
            + ", ".join(missing[:10])
        )

    new_rows = []

    for i, gsm in enumerate(missing, start=1):
        print(
            f"[GSE GEO] {i}/{len(missing)}: {gsm}",
            flush=True,
        )

        try:
            description = fetch_geo_description(gsm)
        except Exception as exc:
            raise RuntimeError(
                "\nCould not complete the GSE266546 GEO bridge.\n"
                "If your compute node does not have internet access, "
                "run once on a login node with internet or provide a "
                "pre-built --gse-geo-cache containing:\n"
                "  tissue_id,geo_description,expression_id\n\n"
                f"Failure at {gsm}: {exc}"
            ) from exc

        new_rows.append({
            "tissue_id": gsm,
            "geo_description": description,
            "expression_id":
                normalize_geo_expression_id(description),
        })

        # Small delay to avoid hammering NCBI.
        time.sleep(0.35)

    if new_rows:
        cache = pd.concat(
            [
                cache,
                pd.DataFrame(new_rows),
            ],
            ignore_index=True,
        )

    cache = standardize_geo_cache(cache)

    cache = (
        cache
        .drop_duplicates(
            [
                "tissue_id",
                "expression_id",
            ]
        )
        .sort_values(
            [
                "tissue_id",
                "expression_id",
            ]
        )
        .reset_index(drop=True)
    )

    cache_file.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    cache.to_csv(
        cache_file,
        index=False,
    )

    return cache


def match_gse(
    expression_ids,
    gse_mapping_file,
    gse_vcf_samples_file,
    geo_cache_file,
    allow_geo_fetch=True,
):
    """
    Match:
        N157 expression ID
          <- GEO Sample_description
          <- tissue_id/GSM
          <- GDA donor IID
          -> canonical GCA donor
          -> exact VCF-head sample ID
    """
    expr = pd.DataFrame({
        "expression_id": expression_ids
    })

    expr["project"] = "GSE266546"

    gda = pd.read_csv(
        gse_mapping_file
    )

    require_columns(
        gda,
        [
            "IID",
            "tissue_id",
        ],
        "GSE266546 GDA mapping",
    )

    gda = gda.copy()

    gda["tissue_id"] = (
        gda["tissue_id"]
        .map(clean_string)
    )

    gda["donor_id"] = (
        gda["IID"]
        .map(canonical_gca)
    )

    if gda["donor_id"].isna().any():
        bad = (
            gda.loc[
                gda["donor_id"].isna(),
                "IID",
            ]
            .drop_duplicates()
            .tolist()
        )

        raise ValueError(
            "Could not normalize GSE IID donor(s):\n"
            + "\n".join(map(str, bad))
        )

    # Exact genotype IDs from actual VCF-head sample list.
    vcf = read_exact_vcf_samples(
        gse_vcf_samples_file
    )

    # GEO bridge: GSM -> YX expression sample.
    geo = build_or_update_geo_cache(
        tissue_ids=gda["tissue_id"].dropna(),
        cache_file=geo_cache_file,
        allow_fetch=allow_geo_fetch,
    )

    # tissue_id -> donor -> expression ID
    bridge = (
        gda[
            [
                "tissue_id",
                "donor_id",
            ]
        ]
        .drop_duplicates()
        .merge(
            geo[
                [
                    "tissue_id",
                    "geo_description",
                    "expression_id",
                ]
            ],
            on="tissue_id",
            how="left",
            validate="many_to_many",
        )
    )

    # Only expression samples actually present in N157.
    bridge = bridge[
        bridge["expression_id"].isin(
            set(expression_ids)
        )
    ].copy()

    # Verify each N157 expression sample maps to only one biological donor.
    donor_n = (
        bridge[
            bridge["expression_id"].notna()
        ]
        .groupby("expression_id")["donor_id"]
        .nunique()
    )

    ambiguous_donor = donor_n[
        donor_n > 1
    ]

    if not ambiguous_donor.empty:
        bad = bridge[
            bridge["expression_id"].isin(
                ambiguous_donor.index
            )
        ].sort_values(
            [
                "expression_id",
                "donor_id",
            ]
        )

        raise RuntimeError(
            "Ambiguous GSE mapping: an expression sample maps "
            "to >1 GCA donor.\n"
            + bad.to_string(index=False)
        )

    # Collapse repeated identical donor-expression mappings.
    bridge = (
        bridge[
            [
                "expression_id",
                "donor_id",
                "tissue_id",
                "geo_description",
            ]
        ]
        .drop_duplicates()
    )

    # If same expression sample has multiple GSM rows but same donor,
    # retain all GSMs as a semicolon-separated audit field.
    bridge = (
        bridge
        .groupby(
            [
                "expression_id",
                "donor_id",
            ],
            as_index=False,
            dropna=False,
        )
        .agg({
            "tissue_id":
                lambda x: ";".join(
                    sorted(
                        set(
                            str(v)
                            for v in x
                            if pd.notna(v)
                        )
                    )
                ),
            "geo_description":
                lambda x: ";".join(
                    sorted(
                        set(
                            str(v)
                            for v in x
                            if pd.notna(v)
                        )
                    )
                ),
        })
    )

    bridge = bridge.merge(
        vcf,
        on="donor_id",
        how="left",
        validate="many_to_one",
    )

    out = expr.merge(
        bridge,
        on="expression_id",
        how="left",
        validate="one_to_one",
    )

    out["mapping_source"] = (
        "N157 expression -> GEO description -> "
        "tissue_id -> IID/GCA -> exact VCF head"
    )

    out["has_genotype"] = (
        out["genotype_id"].notna()
    )

    # GSE266546 currently has no additional NL/AD filter here.
    out["eligible_for_analysis"] = out["has_genotype"]

    out["exclusion_reason"] = ""
    out.loc[
        ~out["has_genotype"],
        "exclusion_reason",
    ] = "GSE266546_without_genotype"

    return out


# =============================================================================
# Reporting
# =============================================================================

def summarize_project(df, project):
    x = df[
        df["project"] == project
    ]

    genotype_matched = x[
        x["has_genotype"]
    ]

    eligible = x[
        x["eligible_for_analysis"]
    ]

    row = {
        "project": project,
        "expression_samples_total":
            x["expression_id"].nunique(),
        "expression_samples_with_genotype":
            genotype_matched["expression_id"].nunique(),
        "expression_samples_final_analysis":
            eligible["expression_id"].nunique(),
        "expression_samples_excluded":
            x.loc[
                ~x["eligible_for_analysis"],
                "expression_id",
            ].nunique(),
        "unique_final_donors":
            eligible["donor_id"].nunique(),
        "unique_final_genotype_samples":
            eligible["genotype_id"].nunique(),
    }

    if project == "HTAN":
        row["HTAN_NL_expression_samples"] = (
            x.loc[
                x["htan_is_NL"].fillna(False),
                "expression_id",
            ].nunique()
        )
        row["HTAN_NL_priority_segment_expression_samples"] = (
            x.loc[
                (
                    x["htan_is_NL"].fillna(False)
                    & x["htan_segment_in_priority"].fillna(False)
                ),
                "expression_id",
            ].nunique()
        )
        row["HTAN_donor_selected_expression_samples"] = (
            x.loc[
                x["htan_selected_for_donor"].fillna(False),
                "expression_id",
            ].nunique()
        )
    else:
        row["HTAN_NL_expression_samples"] = pd.NA
        row["HTAN_NL_priority_segment_expression_samples"] = pd.NA
        row["HTAN_donor_selected_expression_samples"] = pd.NA

    return row


def write_one_column(values, path):
    pd.Series(
        list(values),
        dtype=str,
    ).to_csv(
        path,
        index=False,
        header=False,
    )


# =============================================================================
# Main
# =============================================================================

def main():
    parser = argparse.ArgumentParser(
        description=(
            "Match N157 expression samples to HTAN and GSE266546 "
            "genotype IDs."
        )
    )

    parser.add_argument(
        "--expression-head",
        default="/data/l2_bioinfo1/liq17/scTWAS/data/VUMC_HTAN_GSE266546/expression/N157_Expr_head.csv",
        help="N157 expression-header CSV",
    )

    parser.add_argument(
        "--htan-meta",
        default="/data/l2_bioinfo1/liq17/scTWAS/data/VUMC_HTAN_GSE266546/genotype/VUMC_HTAN/575samples_meta_20240904.csv",
        help="HTAN genotype metadata CSV",
    )

    parser.add_argument(
        "--htan-sample-meta",
        default="/data/l2_bioinfo1/liq17/scTWAS/data/VUMC_HTAN_GSE266546/genotype/VUMC_HTAN/HTAN_metadata.csv",
        help=(
            "HTAN sample metadata CSV containing SPECID.htan, "
            "Sample_Classification, and SEGMENT."
        ),
    )

    parser.add_argument(
        "--gse-map",
        default="/data/l2_bioinfo1/liq17/scTWAS/data/VUMC_HTAN_GSE266546/genotype/VUMC_GSE266546/GDA_genotype_sc_exp_sample_matched.csv",
        help="GSE266546 genotype/scRNA mapping CSV",
    )

    parser.add_argument(
        "--gse-vcf-samples",
        default="/data/l2_bioinfo1/liq17/scTWAS/data/VUMC_HTAN_GSE266546/genotype/VUMC_GSE266546/GSE266546_geno_N144_samples.txt",
        help="GSE266546 exact VCF sample IDs, one per line",
    )

    parser.add_argument(
        "--gse-geo-cache",
        default=None,
        help=(
            "Optional GSM -> expression-ID cache. "
            "If omitted, OUTDIR/"
            "10_GSE266546_GEO_sample_description_cache.csv "
            "is used."
        ),
    )

    parser.add_argument(
        "--no-geo-fetch",
        action="store_true",
        help=(
            "Do not query NCBI GEO. Requires a complete "
            "--gse-geo-cache."
        ),
    )

    parser.add_argument(
        "--outdir",
        default="/data/l2_bioinfo1/liq17/scTWAS/data/VUMC_HTAN_GSE266546/VUMC_HTAN_GSE266546_metadata",
        help="Output directory",
    )

    args = parser.parse_args()

    outdir = Path(args.outdir)

    outdir.mkdir(
        parents=True,
        exist_ok=True,
    )

    geo_cache = (
        Path(args.gse_geo_cache)
        if args.gse_geo_cache
        else outdir
        / "10_GSE266546_GEO_sample_description_cache.csv"
    )

    # -------------------------------------------------------------------------
    # N157 samples
    # -------------------------------------------------------------------------

    all_expr = read_expression_samples(
        args.expression_head
    )

    htan_expr = [
        x
        for x in all_expr
        if x.upper().startswith("HTA11_")
    ]

    gse_expr = [
        x
        for x in all_expr
        if "-YX-" in x.upper()
    ]

    recognized = set(htan_expr) | set(gse_expr)

    unknown = [
        x
        for x in all_expr
        if x not in recognized
    ]

    print("=" * 72)
    print("N157 INPUT")
    print("=" * 72)
    print(f"Total N157 expression samples : {len(all_expr)}")
    print(f"HTAN expression samples       : {len(htan_expr)}")
    print(f"GSE266546 expression samples  : {len(gse_expr)}")
    print(f"Unclassified expression IDs   : {len(unknown)}")

    if unknown:
        print("\nWARNING: unclassified N157 IDs:")
        for x in unknown:
            print("  ", x)

    if len(recognized) != len(all_expr):
        raise RuntimeError(
            "Not all N157 sample IDs could be classified as "
            "HTAN or GSE266546."
        )

    # -------------------------------------------------------------------------
    # HTAN
    # -------------------------------------------------------------------------

    print()
    print("=" * 72)
    print("MATCH HTAN")
    print("=" * 72)

    htan = match_htan(
        expression_ids=htan_expr,
        htan_genotype_meta_file=args.htan_meta,
        htan_sample_meta_file=args.htan_sample_meta,
    )

    # -------------------------------------------------------------------------
    # GSE266546
    # -------------------------------------------------------------------------

    print()
    print("=" * 72)
    print("MATCH GSE266546")
    print("=" * 72)

    gse = match_gse(
        expression_ids=gse_expr,
        gse_mapping_file=args.gse_map,
        gse_vcf_samples_file=args.gse_vcf_samples,
        geo_cache_file=geo_cache,
        allow_geo_fetch=not args.no_geo_fetch,
    )

    # -------------------------------------------------------------------------
    # Combined audit/master tables
    # -------------------------------------------------------------------------

    wanted_cols = [
        "project",
        "expression_id",
        "donor_id",
        "genotype_id",
        "tissue_id",
        "Sample_Classification",
        "SEGMENT",
        "htan_specid_in_metadata",
        "htan_is_NL",
        "htan_segment_priority",
        "htan_segment_in_priority",
        "htan_selection_candidate",
        "htan_selected_for_donor",
        "htan_selected_expression_id_for_donor",
        "htan_selected_segment_for_donor",
        "htan_selected_priority_for_donor",
        "has_genotype",
        "eligible_for_analysis",
        "exclusion_reason",
        "mapping_source",
    ]

    for col in wanted_cols:
        if col not in htan.columns:
            htan[col] = pd.NA
        if col not in gse.columns:
            gse[col] = pd.NA

    audit = pd.concat(
        [
            htan[wanted_cols],
            gse[wanted_cols],
        ],
        ignore_index=True,
    )

    # Preserve the exact N157 column order.
    order = {
        sample: i
        for i, sample in enumerate(all_expr)
    }

    audit["_order"] = (
        audit["expression_id"]
        .map(order)
    )

    audit = (
        audit
        .sort_values("_order")
        .drop(columns="_order")
        .reset_index(drop=True)
    )

    if audit["expression_id"].nunique() != len(all_expr):
        raise RuntimeError(
            "Internal QC failed: combined audit does not contain "
            "exactly one row per N157 expression sample."
        )

    # Final analysis set:
    #   HTAN      = exact SPECID.htan + NL + one donor-selected segment + genotype
    #   GSE266546 = genotype matched
    matched = audit[
        audit["eligible_for_analysis"]
    ].copy()

    excluded = audit[
        ~audit["eligible_for_analysis"]
    ].copy()

    # -------------------------------------------------------------------------
    # Summaries
    # -------------------------------------------------------------------------

    summary_rows = [
        summarize_project(audit, "HTAN"),
        summarize_project(audit, "GSE266546"),
    ]

    all_summary = {
        "project": "ALL",
        "expression_samples_total":
            audit["expression_id"].nunique(),
        "expression_samples_with_genotype":
            audit.loc[
                audit["has_genotype"],
                "expression_id",
            ].nunique(),
        "expression_samples_final_analysis":
            matched["expression_id"].nunique(),
        "expression_samples_excluded":
            excluded["expression_id"].nunique(),
        "unique_final_donors":
            matched[
                [
                    "project",
                    "donor_id",
                ]
            ]
            .drop_duplicates()
            .shape[0],
        "unique_final_genotype_samples":
            matched[
                [
                    "project",
                    "genotype_id",
                ]
            ]
            .drop_duplicates()
            .shape[0],
        "HTAN_NL_expression_samples":
            audit.loc[
                (
                    (audit["project"] == "HTAN")
                    & audit["htan_is_NL"].fillna(False)
                ),
                "expression_id",
            ].nunique(),
        "HTAN_NL_priority_segment_expression_samples":
            audit.loc[
                (
                    (audit["project"] == "HTAN")
                    & audit["htan_is_NL"].fillna(False)
                    & audit["htan_segment_in_priority"].fillna(False)
                ),
                "expression_id",
            ].nunique(),
        "HTAN_donor_selected_expression_samples":
            audit.loc[
                (
                    (audit["project"] == "HTAN")
                    & audit["htan_selected_for_donor"].fillna(False)
                ),
                "expression_id",
            ].nunique(),
    }

    summary_rows.append(all_summary)

    summary = pd.DataFrame(
        summary_rows
    )

    # -------------------------------------------------------------------------
    # Write files
    # -------------------------------------------------------------------------

    audit_file = (
        outdir
        / "01_N157_expression_genotype_audit.csv"
    )

    matched_file = (
        outdir
        / "02_N157_expression_genotype_matched.csv"
    )

    htan_file = (
        outdir
        / "03_HTAN_expression_genotype_matched.csv"
    )

    gse_file = (
        outdir
        / "04_GSE266546_expression_genotype_matched.csv"
    )

    htan_vcf_file = (
        outdir
        / "05_HTAN_vcf_samples_to_keep.txt"
    )

    gse_vcf_file = (
        outdir
        / "06_GSE266546_vcf_samples_to_keep.txt"
    )

    expr_keep_file = (
        outdir
        / "07_N157_expression_samples_with_genotype.txt"
    )

    excluded_file = (
        outdir
        / "08_N157_expression_excluded_from_analysis.csv"
    )

    summary_file = (
        outdir
        / "09_N157_expression_genotype_summary.csv"
    )

    audit.to_csv(
        audit_file,
        index=False,
    )

    matched.to_csv(
        matched_file,
        index=False,
    )

    matched[
        matched["project"] == "HTAN"
    ].to_csv(
        htan_file,
        index=False,
    )

    matched[
        matched["project"] == "GSE266546"
    ].to_csv(
        gse_file,
        index=False,
    )

    write_one_column(
        matched.loc[
            matched["project"] == "HTAN",
            "genotype_id",
        ].drop_duplicates(),
        htan_vcf_file,
    )

    write_one_column(
        matched.loc[
            matched["project"] == "GSE266546",
            "genotype_id",
        ].drop_duplicates(),
        gse_vcf_file,
    )

    write_one_column(
        matched["expression_id"],
        expr_keep_file,
    )

    excluded.to_csv(
        excluded_file,
        index=False,
    )

    summary.to_csv(
        summary_file,
        index=False,
    )

    # -------------------------------------------------------------------------
    # Screen report
    # -------------------------------------------------------------------------

    print()
    print("=" * 72)
    print("FINAL SUMMARY")
    print("=" * 72)

    print(
        summary.to_string(
            index=False
        )
    )

    print()
    print("Master matched table:")
    print(f"  {matched_file}")

    print()
    print("VCF keep lists:")
    print(f"  HTAN:       {htan_vcf_file}")
    print(f"  GSE266546:  {gse_vcf_file}")

    print()
    print("Important:")
    print(
        "  HTAN final analysis requires exact SPECID.htan match, "
        "Sample_Classification == NL, one sample per donor selected by "
        "SEGMENT priority, and genotype."
    )
    print(
        "  GSE266546 final analysis currently requires genotype match."
    )
    print(
        "  HTAN now contributes at most one selected expression specimen "
        "per donor; GSE266546 matching behavior is unchanged."
    )


if __name__ == "__main__":
    main()
