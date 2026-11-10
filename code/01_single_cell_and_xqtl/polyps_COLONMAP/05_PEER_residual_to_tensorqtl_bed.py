#!/usr/bin/env python3
"""Run the shared PEER-to-BED converter from the COLONMAP entry point."""

from pathlib import Path
from runpy import run_path


if __name__ == "__main__":
    converter = (
        Path(__file__).resolve().parents[1]
        / "normal_colon_GSE266546"
        / "xQTL_tensorqtls"
        / "PEER_residual_to_tensorqtl_bed.py"
    )
    run_path(str(converter), run_name="__main__")
