#!/usr/bin/env bash
set -euo pipefail

###############################################################################
# 03_VUMC_GSE201336_merge.sh
#
# Merge VUMC robust/chaining peak regions with GSE201336 peak regions.
#
# Design
# ------
# Disease:
#   VUMC Adenoma+Serrated chaining-merged regions
#       +
#   GSE201336 Adenoma peaks
#
#   Cell types:
#       ABS CT EE GOB STM TAC
#
# Normal:
#   VUMC Normal robust peak regions
#       +
#   GSE201336 Normal peaks
#
#   Cell types available in GSE201336:
#       ABS CT GOB STM TAC
#
# NOTE:
#   GSE201336 Normal__EE.narrowPeak.bed was not listed/available, so EE is not
#   included in the Normal cross-study merge.
#
# IMPORTANT:
#   - No summit-centered +/-250 bp conversion.
#   - Original peak/region boundaries are used.
#   - bedtools merge -d 0 is used.
#   - Transitive/chaining merge is intentionally allowed.
#   - GSE201336 narrowPeak files are converted to BED3 before merging.
#
# CPU:
#   Up to 8 merge jobs in parallel.
###############################################################################

export LC_ALL=C

###############################################################################
# Paths
###############################################################################

VUMC_ROOT="/data/l2_bioinfo1/scATAC-seq/scATACseq_VUMC/MACS3_callpeak_q01"

# From 02 pipeline:
VUMC_DISEASE_DIR="${VUMC_ROOT}/Adenoma_Serrated_chaining_merged_peaks"
VUMC_ROBUST_DIR="${VUMC_ROOT}/robust_peak_regions"

GSE_DIR="/data/l2_bioinfo1/scATAC-seq/scATACseq_GSE201336/GSE201336_scATAC_peaks/peaks/per_group"

OUT="${VUMC_ROOT}/VUMC_GSE201336_merged_peaks"

N_PARALLEL=8

mkdir -p \
    "${OUT}/disease" \
    "${OUT}/normal" \
    "${OUT}/qc"

###############################################################################
# Software checks
###############################################################################

for cmd in awk sort cut bedtools xargs wc; do
    if ! command -v "$cmd" >/dev/null 2>&1; then
        echo "ERROR: required command not found: $cmd" >&2
        exit 1
    fi
done

###############################################################################
# Utility: median interval length
###############################################################################

median_bed_length() {
    local bed="$1"

    awk 'NF>=3 {print $3-$2}' "$bed" \
    | sort -n \
    | awk '
        {
            x[NR]=$1
        }
        END {
            if (NR == 0) {
                print "NA"
            } else if (NR % 2 == 1) {
                print x[(NR+1)/2]
            } else {
                printf "%.1f\n", (x[NR/2] + x[NR/2+1]) / 2
            }
        }
    '
}

###############################################################################
# Utility: mean interval length
###############################################################################

mean_bed_length() {
    local bed="$1"

    awk '
        NF>=3 {
            sum += ($3-$2)
            n++
        }
        END {
            if (n == 0)
                print "NA"
            else
                printf "%.2f\n", sum/n
        }
    ' "$bed"
}

###############################################################################
# Merge one VUMC + GSE201336 pair
#
# Arguments:
#   1 mode       disease | normal
#   2 cell type
#   3 VUMC BED
#   4 GSE BED/narrowPeak
###############################################################################

merge_one() {
    set -euo pipefail

    mode="$1"
    ct="$2"
    vumc="$3"
    gse="$4"

    if [[ ! -s "$vumc" ]]; then
        echo "ERROR: missing VUMC input: $vumc" >&2
        exit 1
    fi

    if [[ ! -s "$gse" ]]; then
        echo "ERROR: missing GSE201336 input: $gse" >&2
        exit 1
    fi

    if [[ "$mode" == "disease" ]]; then
        final="${OUT}/disease/VUMC_AdenomaSerrated_GSE201336_Adenoma_${ct}.chaining_merge.bed"
    elif [[ "$mode" == "normal" ]]; then
        final="${OUT}/normal/VUMC_GSE201336_Normal_${ct}.chaining_merge.bed"
    else
        echo "ERROR: unknown mode: $mode" >&2
        exit 1
    fi

    qc="${OUT}/qc/${mode}_${ct}.tsv"

    echo "=============================================================="
    echo "Mode      : $mode"
    echo "Cell type : $ct"
    echo "VUMC      : $vumc"
    echo "GSE201336 : $gse"
    echo "=============================================================="

    ###########################################################################
    # Temporary BED3 versions
    ###########################################################################

    tmp_vumc=$(mktemp)
    tmp_gse=$(mktemp)

    trap 'rm -f "$tmp_vumc" "$tmp_gse"' RETURN

    cut -f1-3 "$vumc" \
    | sort -k1,1 -k2,2n \
    > "$tmp_vumc"

    cut -f1-3 "$gse" \
    | sort -k1,1 -k2,2n \
    > "$tmp_gse"

    ###########################################################################
    # Input counts
    ###########################################################################

    n_vumc=$(wc -l < "$tmp_vumc")
    n_gse=$(wc -l < "$tmp_gse")

    ###########################################################################
    # Cross-study overlap statistics
    #
    # VUMC_overlap_GSE:
    #   number of VUMC regions overlapping >=1 GSE201336 peak
    #
    # GSE_overlap_VUMC:
    #   number of GSE201336 peaks overlapping >=1 VUMC region
    ###########################################################################

    n_vumc_overlap=$(
        bedtools intersect -u \
            -a "$tmp_vumc" \
            -b "$tmp_gse" \
        | wc -l
    )

    n_gse_overlap=$(
        bedtools intersect -u \
            -a "$tmp_gse" \
            -b "$tmp_vumc" \
        | wc -l
    )

    pct_vumc_overlap=$(
        awk -v x="$n_vumc_overlap" -v n="$n_vumc" \
            'BEGIN{if(n>0) printf "%.2f",100*x/n; else print "NA"}'
    )

    pct_gse_overlap=$(
        awk -v x="$n_gse_overlap" -v n="$n_gse" \
            'BEGIN{if(n>0) printf "%.2f",100*x/n; else print "NA"}'
    )

    ###########################################################################
    # Chaining merge:
    # define genomic regions covered by either VUMC or GSE201336 peaks.
    ###########################################################################

    cat "$tmp_vumc" "$tmp_gse" \
    | sort -k1,1 -k2,2n \
    | bedtools merge -d 0 \
    > "$final"

    n_final=$(wc -l < "$final")
    med_final=$(median_bed_length "$final")
    mean_final=$(mean_bed_length "$final")

    ###########################################################################
    # Summary row
    ###########################################################################

    printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
        "$mode" \
        "$ct" \
        "$n_vumc" \
        "$n_gse" \
        "$n_vumc_overlap" \
        "$pct_vumc_overlap" \
        "$n_gse_overlap" \
        "$pct_gse_overlap" \
        "$n_final" \
        "$med_final" \
        "$mean_final" \
        > "$qc"

    echo "DONE: ${mode} ${ct}"
    echo "  VUMC regions                         : $n_vumc"
    echo "  GSE201336 peaks                      : $n_gse"
    echo "  VUMC regions overlapping GSE         : $n_vumc_overlap (${pct_vumc_overlap}%)"
    echo "  GSE peaks overlapping VUMC           : $n_gse_overlap (${pct_gse_overlap}%)"
    echo "  Regions after cross-study chain merge: $n_final"
    echo "  Median merged region length (bp)     : $med_final"
    echo "  Mean merged region length (bp)       : $mean_final"
}

###############################################################################
# Export for xargs workers
###############################################################################

export -f median_bed_length
export -f mean_bed_length
export -f merge_one

export OUT

###############################################################################
# Build job table
#
# Columns:
#   mode
#   cell_type
#   VUMC_input
#   GSE201336_input
###############################################################################

JOBS="${OUT}/merge_jobs.tsv"

cat > "$JOBS" <<EOF
disease	ABS	${VUMC_DISEASE_DIR}/Adenoma_Serrated_ABS.chaining_merge.bed	${GSE_DIR}/Adenoma__ABS.narrowPeak.bed
disease	CT	${VUMC_DISEASE_DIR}/Adenoma_Serrated_CT.chaining_merge.bed	${GSE_DIR}/Adenoma__CT.narrowPeak.bed
disease	EE	${VUMC_DISEASE_DIR}/Adenoma_Serrated_EE.chaining_merge.bed	${GSE_DIR}/Adenoma__EE.narrowPeak.bed
disease	GOB	${VUMC_DISEASE_DIR}/Adenoma_Serrated_GOB.chaining_merge.bed	${GSE_DIR}/Adenoma__GOB.narrowPeak.bed
disease	STM	${VUMC_DISEASE_DIR}/Adenoma_Serrated_STM.chaining_merge.bed	${GSE_DIR}/Adenoma__STM.narrowPeak.bed
disease	TAC	${VUMC_DISEASE_DIR}/Adenoma_Serrated_TAC.chaining_merge.bed	${GSE_DIR}/Adenoma__TAC.narrowPeak.bed
normal	ABS	${VUMC_ROBUST_DIR}/Normal_ABS.robust.bed	${GSE_DIR}/Normal__ABS.narrowPeak.bed
normal	CT	${VUMC_ROBUST_DIR}/Normal_CT.robust.bed	${GSE_DIR}/Normal__CT.narrowPeak.bed
normal	GOB	${VUMC_ROBUST_DIR}/Normal_GOB.robust.bed	${GSE_DIR}/Normal__GOB.narrowPeak.bed
normal	STM	${VUMC_ROBUST_DIR}/Normal_STM.robust.bed	${GSE_DIR}/Normal__STM.narrowPeak.bed
normal	TAC	${VUMC_ROBUST_DIR}/Normal_TAC.robust.bed	${GSE_DIR}/Normal__TAC.narrowPeak.bed
EOF

###############################################################################
# Validate all inputs before launching parallel jobs
###############################################################################

echo
echo "=============================================================="
echo "CHECKING INPUT FILES"
echo "=============================================================="
echo

while IFS=$'\t' read -r mode ct vumc gse; do

    if [[ ! -s "$vumc" ]]; then
        echo "ERROR: missing VUMC file:"
        echo "  $vumc"
        exit 1
    fi

    if [[ ! -s "$gse" ]]; then
        echo "ERROR: missing GSE201336 file:"
        echo "  $gse"
        exit 1
    fi

    echo "OK  ${mode}  ${ct}"

done < "$JOBS"

###############################################################################
# Run merge jobs in parallel
###############################################################################

echo
echo "=============================================================="
echo "VUMC + GSE201336 CROSS-STUDY CHAINING MERGE"
echo "Parallel jobs: ${N_PARALLEL}"
echo "=============================================================="
echo

xargs -r -P "$N_PARALLEL" -n 4 \
    bash -c 'merge_one "$1" "$2" "$3" "$4"' _ \
    < "$JOBS"

###############################################################################
# Combined summary
###############################################################################

SUMMARY="${OUT}/VUMC_GSE201336_merge_summary.tsv"

echo -e \
"Group\tCellType\tNum_VUMC_regions\tNum_GSE201336_peaks\tNum_VUMC_regions_overlapping_GSE201336\tPct_VUMC_regions_overlapping_GSE201336\tNum_GSE201336_peaks_overlapping_VUMC\tPct_GSE201336_peaks_overlapping_VUMC\tNum_regions_after_chain_merge\tMedian_merged_region_length_bp\tMean_merged_region_length_bp" \
> "$SUMMARY"

###############################################################################
# Disease rows
###############################################################################

for ct in ABS CT EE GOB STM TAC; do
    cat "${OUT}/qc/disease_${ct}.tsv" >> "$SUMMARY"
done

###############################################################################
# Normal rows
###############################################################################

for ct in ABS CT GOB STM TAC; do
    cat "${OUT}/qc/normal_${ct}.tsv" >> "$SUMMARY"
done

###############################################################################
# Final report
###############################################################################

echo
echo "=============================================================="
echo "VUMC + GSE201336 MERGE COMPLETE"
echo "=============================================================="
echo

echo "Disease merged peak-region files:"
for ct in ABS CT EE GOB STM TAC; do
    echo "  ${OUT}/disease/VUMC_AdenomaSerrated_GSE201336_Adenoma_${ct}.chaining_merge.bed"
done

echo
echo "Normal merged peak-region files:"
for ct in ABS CT GOB STM TAC; do
    echo "  ${OUT}/normal/VUMC_GSE201336_Normal_${ct}.chaining_merge.bed"
done

echo
echo "Summary:"
column -t "$SUMMARY"

echo
echo "Summary file:"
echo "  $SUMMARY"
