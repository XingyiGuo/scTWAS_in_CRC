import pandas as pd
import re
import numpy as np
from difflib import SequenceMatcher


CRC_TERMS = [
    "colorectal",
    "colon",
    "rectal",
    "rectum",
    "crc",
]

CANCER_TERMS = [
    "cancer",
    "neoplasm",
    "tumor",
    "tumour",
    "carcinoma",
    "glioma",
    "glioblastoma",
    "leukemia",
    "lymphoma",
    "myeloma",
    "melanoma",
    "sarcoma",
    "blastoma",
    # Malignancies whose names contain none of the words above. Each term
    # was checked against every indication string present in the table:
    # all of its matches are genuine cancers, no false positives.
    "hodgkin",         # Hodgkin / non-Hodgkin disease
    "myelodysplas",    # myelodysplastic syndromes
    "myeloprolifer",   # myeloproliferative neoplasms
    "kaposi",          # Kaposi (written without "sarcoma" in some records)
    "thelioma",        # mesothelioma
    "metasta",         # metastatic / metastasis
]


def diseases_protein_drugs_instruction_json(df: pd.DataFrame) -> dict:
    """
    Build hierarchy for radial tree plots.

    Structure:
        Colorectal cancer
            -> Gene
            -> Drug (value = "true"/"false" approval status)
            -> Top 3 cancer indications

    Indication priority:
        1. Colorectal/colon/rectal indications
        2. Other cancer indications

    Notes:
        - CellType is NOT used in the hierarchy.
        - If the same gene/drug appears in multiple cell types or rows,
          those rows are combined.
    """

    required_cols = [
        "Genes",
        "Drug Brand Name",
        "Max Phase",
        "Drug Indications (DailyMed or ATC)",
        "Drug Indications (ClinicalTrials)",
    ]

    missing = [c for c in required_cols if c not in df.columns]
    if missing:
        raise ValueError(f"Missing required columns: {missing}")

    df = df.copy()

    # ---------------------------------------------------------------
    # Approved / Phase 4 flag
    # ---------------------------------------------------------------
    df["MaxPhase plot"] = np.where(
        df["Max Phase"]
        .astype(str)
        .str.lower()
        .str.contains("approved|4", regex=True),
        "true",
        "false",
    )

    # ---------------------------------------------------------------
    # Indications: UNION of both fields, not DailyMed-first-then-fallback.
    # Coalescing hid trial oncology indications for any drug that already
    # carried a non-cancer approved label, so "has a cancer indication" was
    # really "has a cancer indication in its highest-priority field".
    # ---------------------------------------------------------------
    _dm = df["Drug Indications (DailyMed or ATC)"].fillna("").astype(str).str.strip()
    _ct = df["Drug Indications (ClinicalTrials)"].fillna("").astype(str).str.strip()
    df["Drug Indications"] = (_dm + "; " + _ct).str.strip("; ").replace("", np.nan)

    # ---------------------------------------------------------------
    # Helper functions
    # ---------------------------------------------------------------
    def _deduplicate(items: list) -> list:
        seen = set()
        out = []

        for x in items:
            key = x.lower()

            if key not in seen:
                seen.add(key)
                out.append(x)

        return out

    def _is_cancer(indication: str) -> bool:
        x = indication.lower()

        return (
            any(k in x for k in CRC_TERMS)
            or any(k in x for k in CANCER_TERMS)
        )

    def _score(indication: str) -> tuple:
        """
        CRC-related indications first,
        then all other cancers alphabetically.
        """
        x = indication.lower()

        if any(k in x for k in CRC_TERMS):
            return (0, x)

        return (1, x)

    def _extract_top_indications(
        indication_series: pd.Series,
        max_n: int | None = None,
        similarity_threshold: float = 0.80,
    ) -> list:

        if indication_series.isna().all():
            return []

        merged = ";".join(
            indication_series
            .dropna()
            .astype(str)
        )

        items = [
            x.strip()
            for x in re.split(r"[;,|]", merged)
            if x.strip()
        ]

        # Exact duplicate removal
        items = _deduplicate(items)

        # Keep cancer indications only
        items = [
            x
            for x in items
            if _is_cancer(x)
        ]

        # CRC first
        items = sorted(
            items,
            key=_score,
        )

        # Remove very similar indications
        result = []

        for item in items:

            if max_n is not None and len(result) >= max_n:
                break

            too_similar = any(
                SequenceMatcher(
                    None,
                    item.lower(),
                    kept.lower(),
                ).ratio()
                >= similarity_threshold
                for kept in result
            )

            if not too_similar:
                result.append(item)

        return result

    def _extract_phase(
        per_drug_df: pd.DataFrame,
    ) -> str:

        phases = (
            per_drug_df["MaxPhase plot"]
            .dropna()
            .astype(str)
            .str.lower()
        )

        if phases.empty:
            return "false"

        # If ANY occurrence of this drug is approved / Phase 4,
        # mark the drug node as true.
        return (
            "true"
            if (phases == "true").any()
            else "false"
        )

    # ===============================================================
    # Build hierarchy:
    #
    # CRC
    #   -> Gene
    #       -> Drug
    #           -> Indication
    #
    # CellType layer removed
    # ===============================================================

    gene_list = []

    for gene in sorted(
        df["Genes"]
        .dropna()
        .astype(str)
        .unique()
    ):

        per_gene_df = df[
            df["Genes"].astype(str) == gene
        ]

        drug_list = []

        for drug, per_drug_df in per_gene_df.groupby(
            "Drug Brand Name",
            sort=True,
        ):

            # Skip missing drug names
            if pd.isna(drug):
                continue

            drug = str(drug).strip()

            if drug == "":
                continue

            top_indications = _extract_top_indications(
                per_drug_df["Drug Indications"],
                max_n=None,   # no cap: every cancer indication is kept
                similarity_threshold=0.80,
            )

            drug_list.append({
                "name": drug,

                # Used for drug-node color
                "value": _extract_phase(
                    per_drug_df
                ),

                "children": [
                    {
                        "name": ind.capitalize()
                    }
                    for ind in top_indications
                ],
            })

        # Only include genes that have >=1 drug
        if drug_list:
            gene_list.append({
                "name": gene,
                "children": drug_list,
            })

    return {
        "name": "Colorectal cancer",
        "children": gene_list,
    }