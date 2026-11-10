#!/usr/bin/env bash
set -euo pipefail

# Keep the COLONMAP entry point for existing job commands.
script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
exec bash "${script_dir}/../normal_colon_GSE266546/xQTL_tensorqtls/prepare_tensorqtl_genotype.sh" "$@"
