#!/usr/bin/env bash
set -euo pipefail

###############################################################################
# VUMC ATAC: robust peak selection -> summit-centered 500-bp peaks ->
#            Adenoma + Serrated cell-type-specific merged peak sets
#
# CPU strategy:
#   Stage 1: prepare robust summit500 peaks for up to 8 samples in parallel
#   Stage 2: merge 8 cell-type pairs in parallel
#
# Robust peak rule:
#   >=1M usable fragments:
#       use full blacklist-clean MACS3 narrowPeak
#   <1M usable fragments:
#       use peaks supported in >=2/3 80% resampling runs (stable2of3)
#
# Low-depth samples used here:
#   Adenoma_ASC
#   Adenoma_TUF
#   Serrated_SSC
#   Serrated_TUF
#
# Final 8 combined populations:
#   ABS, CT, EE, GOB, STM, TAC, TUF, ASC_SSC
###############################################################################

ROOT="/data/l2_bioinfo1/scATAC-seq/scATACseq_VUMC/MACS3_callpeak_q01"
PEAK_ROOT="${ROOT}/peaks"
ROBUST_DIR="${ROOT}/robust_fixedwidth_peaks"
OUT="${ROOT}/Adenoma_Serrated_combined_peaks"

N_PARALLEL=8
SUMMIT_HALF_WIDTH=250

mkdir -p "$ROBUST_DIR" "$OUT"

###############################################################################
# Check required software
###############################################################################

for cmd in awk sort cut bedtools xargs wc; do
    if ! command -v "$cmd" >/dev/null 2>&1; then
        echo "ERROR: required command not found: $cmd" >&2
        exit 1
    fi
done

###############################################################################
# Determine whether one sample is low depth
###############################################################################

is_low_depth() {
    case "$1" in
        Adenoma_ASC|Adenoma_TUF|Serrated_SSC|Serrated_TUF)
            return 0
            ;;
        *)
            return 1
            ;;
    esac
}

###############################################################################
# Stage 1 function:
# create robust narrowPeak + summit-centered fixed-width BED for one sample
###############################################################################

prepare_summit500() {
    set -euo pipefail

    sample="$1"

    sample_dir="${PEAK_ROOT}/${sample}"
    clean="${sample_dir}/${sample}.clean.narrowPeak"
    stable="${sample_dir}/${sample}.stable2of3.bed"

    robust="${ROBUST_DIR}/${sample}.robust.narrowPeak"
    fixed="${ROBUST_DIR}/${sample}.robust.summit500.bed"
    qc="${ROBUST_DIR}/${sample}.summit500.counts.tsv"

    echo "--------------------------------------------------"
    echo "Preparing robust summit500: ${sample}"
    echo "--------------------------------------------------"

    if [[ ! -s "$clean" ]]; then
        echo "ERROR: missing or empty clean narrowPeak: $clean" >&2
        exit 1
    fi

    n_clean=$(wc -l < "$clean")
    n_stable="NA"
    source_type="full_clean"

    ###########################################################################
    # Low-depth samples: retain only stable2of3 peaks.
    # stable2of3 column 4 stores the original MACS3 peak name.
    # Recover the complete narrowPeak row so MACS3 summit offset (col 10)
    # is preserved.
    ###########################################################################

    if is_low_depth "$sample"; then
        source_type="stable2of3"

        if [[ ! -s "$stable" ]]; then
            echo "ERROR: low-depth sample is missing stable2of3 file: $stable" >&2
            exit 1
        fi

        n_stable=$(wc -l < "$stable")

        awk 'BEGIN{OFS="\t"}
             NR==FNR { keep[$4]=1; next }
             ($4 in keep) { print }' \
            "$stable" "$clean" > "$robust"

        n_robust=$(wc -l < "$robust")

        # Every stable peak should map back to one clean narrowPeak record.
        if [[ "$n_robust" -ne "$n_stable" ]]; then
            echo "ERROR: ${sample}: stable2of3 count (${n_stable}) != recovered robust narrowPeak count (${n_robust})." >&2
            echo "Check whether column 4 peak names match between:" >&2
            echo "  $stable" >&2
            echo "  $clean" >&2
            exit 1
        fi
    else
        cp "$clean" "$robust"
        n_robust="$n_clean"
    fi

    ###########################################################################
    # Validate narrowPeak summit field.
    # narrowPeak column 10 = summit offset from chromStart.
    ###########################################################################

    invalid_summits=$(awk 'NF < 10 || $10 !~ /^[0-9]+$/ {n++} END{print n+0}' "$robust")

    if [[ "$invalid_summits" -gt 0 ]]; then
        echo "ERROR: ${sample}: ${invalid_summits} robust peaks lack a valid MACS3 summit in column 10." >&2
        exit 1
    fi

    ###########################################################################
    # Convert each robust MACS3 peak to an exactly 500-bp interval centered
    # on its summit:
    #
    #   summit = chromStart + summitOffset
    #   BED interval = [summit-250, summit+250)
    #
    # Output columns:
    #   chr, start, end, peak_name, MACS3_score
    ###########################################################################

    awk -v hw="$SUMMIT_HALF_WIDTH" 'BEGIN{OFS="\t"}
        {
            summit = $2 + $10
            start  = summit - hw
            end    = summit + hw

            if (start < 0)
                start = 0

            print $1, start, end, $4, $5
        }' "$robust" \
    | sort -k1,1 -k2,2n \
    > "$fixed"

    n_fixed=$(wc -l < "$fixed")

    if [[ "$n_fixed" -ne "$n_robust" ]]; then
        echo "ERROR: ${sample}: robust peak count (${n_robust}) != summit500 count (${n_fixed})." >&2
        exit 1
    fi

    printf "%s\t%s\t%s\t%s\t%s\t%s\n" \
        "$sample" \
        "$source_type" \
        "$n_clean" \
        "$n_stable" \
        "$n_robust" \
        "$n_fixed" \
        > "$qc"

    echo "DONE summit500: ${sample} | source=${source_type} | peaks=${n_fixed}"
}

###############################################################################
# Stage 2 function:
# combine one Adenoma + Serrated cell-type pair
###############################################################################

merge_one() {
    set -euo pipefail

    label="$1"
    adenoma_ct="$2"
    serrated_ct="$3"

    A="${ROBUST_DIR}/Adenoma_${adenoma_ct}.robust.summit500.bed"
    S="${ROBUST_DIR}/Serrated_${serrated_ct}.robust.summit500.bed"

    FINAL="${OUT}/Adenoma_Serrated_${label}.peaks.bed"
    QC="${OUT}/${label}.counts.tsv"

    echo "=================================================="
    echo "Merging: ${label}"
    echo "Adenoma : ${A}"
    echo "Serrated: ${S}"
    echo "=================================================="

    if [[ ! -s "$A" ]]; then
        echo "ERROR: missing $A" >&2
        exit 1
    fi

    if [[ ! -s "$S" ]]; then
        echo "ERROR: missing $S" >&2
        exit 1
    fi

    ###########################################################################
    # Combine summit500 coordinates and merge overlapping intervals.
    # Final output is BED3 because source-specific scores/names are no longer
    # interpretable after genomic union.
    ###########################################################################

    cat \
        <(cut -f1-3 "$A") \
        <(cut -f1-3 "$S") \
    | sort -k1,1 -k2,2n \
    | bedtools merge \
    > "$FINAL"

    nA=$(wc -l < "$A")
    nS=$(wc -l < "$S")
    nFinal=$(wc -l < "$FINAL")

    printf "%s\t%s\t%s\t%s\n" \
        "$label" "$nA" "$nS" "$nFinal" \
        > "$QC"

    echo "DONE merge: ${label}"
    echo "  Adenoma summit500 : ${nA}"
    echo "  Serrated summit500: ${nS}"
    echo "  Combined peaks     : ${nFinal}"
}

###############################################################################
# Export functions / variables for xargs workers
###############################################################################

export -f is_low_depth
export -f prepare_summit500
export -f merge_one

export ROOT
export PEAK_ROOT
export ROBUST_DIR
export OUT
export SUMMIT_HALF_WIDTH

###############################################################################
# Samples required for the 8 final combined populations
###############################################################################

cat <<'EOF_SAMPLES' > "${ROBUST_DIR}/samples_to_prepare.txt"
Adenoma_ABS
Serrated_ABS
Adenoma_CT
Serrated_CT
Adenoma_EE
Serrated_EE
Adenoma_GOB
Serrated_GOB
Adenoma_STM
Serrated_STM
Adenoma_TAC
Serrated_TAC
Adenoma_TUF
Serrated_TUF
Adenoma_ASC
Serrated_SSC
EOF_SAMPLES

###############################################################################
# Clean only temporary/QC outputs from a previous run.
# Do not delete peak files until each sample is successfully regenerated.
###############################################################################

rm -f "${ROBUST_DIR}"/*.summit500.counts.tsv
rm -f "${OUT}"/*.counts.tsv

###############################################################################
# STAGE 1: prepare robust summit500 peaks, max 8 samples simultaneously
###############################################################################

echo
echo "=============================================================="
echo "STAGE 1: ROBUST PEAK SELECTION + SUMMIT500"
echo "Maximum parallel jobs: ${N_PARALLEL}"
echo "=============================================================="
echo

xargs -r -n 1 -P "$N_PARALLEL" \
    bash -c 'prepare_summit500 "$1"' _ \
    < "${ROBUST_DIR}/samples_to_prepare.txt"

###############################################################################
# Combine summit500 preparation summary
###############################################################################

SUMMIT_SUMMARY="${ROBUST_DIR}/robust_summit500_summary.tsv"

echo -e "Sample\tPeak_source\tClean_peaks\tStable2of3_peaks\tRobust_peaks\tSummit500_peaks" \
    > "$SUMMIT_SUMMARY"

for sample in \
    Adenoma_ABS Serrated_ABS \
    Adenoma_CT Serrated_CT \
    Adenoma_EE Serrated_EE \
    Adenoma_GOB Serrated_GOB \
    Adenoma_STM Serrated_STM \
    Adenoma_TAC Serrated_TAC \
    Adenoma_TUF Serrated_TUF \
    Adenoma_ASC Serrated_SSC; do

    cat "${ROBUST_DIR}/${sample}.summit500.counts.tsv" >> "$SUMMIT_SUMMARY"
done

echo
echo "Summit500 summary:"
column -t "$SUMMIT_SUMMARY"

###############################################################################
# Define the 8 Adenoma/Serrated combinations
###############################################################################

cat <<'EOF_PAIRS' > "${OUT}/peak_pairs.tsv"
ABS	ABS	ABS
CT	CT	CT
EE	EE	EE
GOB	GOB	GOB
STM	STM	STM
TAC	TAC	TAC
TUF	TUF	TUF
ASC_SSC	ASC	SSC
EOF_PAIRS

###############################################################################
# STAGE 2: run all 8 cell-type merges in parallel
###############################################################################

echo
echo "=============================================================="
echo "STAGE 2: ADENOMA + SERRATED PEAK MERGING"
echo "Parallel jobs: ${N_PARALLEL}"
echo "=============================================================="
echo

xargs -r -P "$N_PARALLEL" -n 3 \
    bash -c 'merge_one "$1" "$2" "$3"' _ \
    < "${OUT}/peak_pairs.tsv"

###############################################################################
# Combine final merge summary
###############################################################################

SUMMARY="${OUT}/combined_peak_counts.tsv"

echo -e "CellType\tAdenoma_summit500_peaks\tSerrated_summit500_peaks\tCombined_peaks" \
    > "$SUMMARY"

for label in ABS CT EE GOB STM TAC TUF ASC_SSC; do
    cat "${OUT}/${label}.counts.tsv" >> "$SUMMARY"
done

###############################################################################
# Final report
###############################################################################

echo
echo "=============================================================="
echo "ALL PEAK PREPARATION + MERGING COMPLETE"
echo "=============================================================="
echo

echo "Final combined peak sets:"
for label in ABS CT EE GOB STM TAC TUF ASC_SSC; do
    echo "  ${OUT}/Adenoma_Serrated_${label}.peaks.bed"
done

echo
echo "Combined peak summary:"
column -t "$SUMMARY"

echo
echo "Summit500 summary:"
echo "  $SUMMIT_SUMMARY"
echo
echo "Combined peak summary:"
echo "  $SUMMARY"
