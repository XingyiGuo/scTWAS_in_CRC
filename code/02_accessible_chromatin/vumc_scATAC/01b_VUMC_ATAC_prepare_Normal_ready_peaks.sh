#!/usr/bin/env bash
set -euo pipefail

###############################################################################
# Populate Normal peak BEDs in ready_full/ and ready_robust/
# without rerunning MACS3 peak calling.
#
# Existing peak calls are read from:
#   MACS3_callpeak_q01/peaks/Normal_<CT>/
#
# Rules:
#   ready_full   = all blacklist-clean MACS3 peak regions
#   ready_robust = >=1M -> full clean peaks
#                  <1M  -> stable2of3 peaks
#
# For the current VUMC data, Normal_TAC is <1M and therefore uses stable2of3
# in ready_robust. Other Normal cell types use full clean peaks.
###############################################################################

ROOT="/data/l2_bioinfo1/scATAC-seq/scATACseq_VUMC/MACS3_callpeak_q01"
PEAK_ROOT="${ROOT}/peaks"
SUMMARY="${ROOT}/peak_calling_summary.tsv"

READY_FULL="${ROOT}/ready_full"
READY_ROBUST="${ROOT}/ready_robust"

mkdir -p "$READY_FULL" "$READY_ROBUST"

CELLTYPES=(ABS CT EE GOB STM TAC TUF)

OUT_SUMMARY="${ROOT}/Normal_ready_peak_summary.tsv"

echo -e \
"Sample\tDepth_class\tFull_clean_peaks\tStable2of3_peaks\tReady_robust_source\tReady_robust_peaks" \
> "$OUT_SUMMARY"

for ct in "${CELLTYPES[@]}"; do

    sample="Normal_${ct}"
    sample_dir="${PEAK_ROOT}/${sample}"

    full="${sample_dir}/${sample}.clean.narrowPeak"
    stable="${sample_dir}/${sample}.stable2of3.bed"

    if [[ ! -s "$full" ]]; then
        echo "ERROR: missing $full" >&2
        exit 1
    fi

    depth_class=$(
        awk -F'\t' -v s="$sample" '
            NR>1 && $1==s {print $3; exit}
        ' "$SUMMARY"
    )

    if [[ -z "$depth_class" ]]; then
        echo "ERROR: cannot find ${sample} in ${SUMMARY}" >&2
        exit 1
    fi

    ###########################################################################
    # ready_full: always use blacklist-clean full MACS3 peaks
    ###########################################################################

    cut -f1-3 "$full" \
        > "${READY_FULL}/${sample}.bed"

    n_full=$(wc -l < "$full")
    n_stable="NA"

    ###########################################################################
    # ready_robust
    ###########################################################################

    if [[ "$depth_class" == "LT1M" ]]; then

        if [[ ! -s "$stable" ]]; then
            echo "ERROR: ${sample} is LT1M but stable2of3 is missing: $stable" >&2
            exit 1
        fi

        cut -f1-3 "$stable" \
            > "${READY_ROBUST}/${sample}.bed"

        n_stable=$(wc -l < "$stable")
        robust_source="stable2of3"
        n_robust="$n_stable"

    else

        cut -f1-3 "$full" \
            > "${READY_ROBUST}/${sample}.bed"

        robust_source="full_clean"
        n_robust="$n_full"

    fi

    printf "%s\t%s\t%s\t%s\t%s\t%s\n" \
        "$sample" \
        "$depth_class" \
        "$n_full" \
        "$n_stable" \
        "$robust_source" \
        "$n_robust" \
        >> "$OUT_SUMMARY"

done

echo
echo "=============================================================="
echo "NORMAL READY PEAK FILES CREATED"
echo "=============================================================="
echo
column -t "$OUT_SUMMARY"
echo
echo "ready_full:"
ls -lh "${READY_FULL}"/Normal_*.bed
echo
echo "ready_robust:"
ls -lh "${READY_ROBUST}"/Normal_*.bed
