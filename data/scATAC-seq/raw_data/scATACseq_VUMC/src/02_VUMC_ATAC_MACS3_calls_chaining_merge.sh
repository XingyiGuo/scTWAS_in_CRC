#!/usr/bin/env bash
set -euo pipefail

###############################################################################
# VUMC ATAC: robust peak selection -> direct chaining merge of peak regions
#
# Purpose
# -------
# Define genomic regions covered by robust ATAC-seq peaks.
#
# IMPORTANT:
#   - NO summit-centered +/-250 bp conversion is performed.
#   - Original robust peak boundaries are retained.
#   - Overlapping/book-ended regions are merged with bedtools merge.
#   - Transitive/chaining merge is intentionally allowed.
#
# Robust peak rule
# ----------------
# >=1M usable fragments:
#     use full blacklist-clean MACS3 peak regions
#
# <1M usable fragments:
#     use peaks supported in >=2/3 80% resampling runs (stable2of3)
#
# Low-depth samples used here:
#     Adenoma_ASC
#     Adenoma_TUF
#     Serrated_SSC
#     Serrated_TUF
#
# Matched Adenoma + Serrated cell types to chain-merge:
#     ABS CT EE GOB STM TAC TUF
#
# Normal cell types analyzed separately (not chain-merged with disease),
# but their statistics are included in the same final summary rows:
#     ABS CT EE GOB STM TAC TUF
#
# ASC and SSC are NOT merged:
#     Adenoma_ASC -> retain/report original robust BED
#     Serrated_SSC -> retain/report original robust BED
#
# CPU strategy:
#     prepare robust BEDs in parallel
#     chain-merge matched cell types in parallel
###############################################################################

export LC_ALL=C

ROOT="/data/l2_bioinfo1/scATAC-seq/scATACseq_VUMC/MACS3_callpeak_q01"
PEAK_ROOT="${ROOT}/peaks"
ROBUST_DIR="${ROOT}/robust_peak_regions"
OUT="${ROOT}/Adenoma_Serrated_chaining_merged_peaks"

N_PARALLEL=8

mkdir -p "$ROBUST_DIR" "$OUT"

###############################################################################
# Check required software
###############################################################################

for cmd in awk sort cut bedtools xargs wc cp; do
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
        Adenoma_ASC|Adenoma_TUF|Serrated_SSC|Serrated_TUF|Normal_TAC)
            return 0
            ;;
        *)
            return 1
            ;;
    esac
}

###############################################################################
# Median interval length
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
# Mean interval length
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
# Stage 1:
# Prepare one robust BED using ORIGINAL peak boundaries.
#
# >=1M:
#     clean.narrowPeak -> BED3
#
# <1M:
#     stable2of3.bed -> BED3 directly
#
# No summit conversion.
###############################################################################

prepare_robust_bed() {
    set -euo pipefail

    sample="$1"

    sample_dir="${PEAK_ROOT}/${sample}"
    clean="${sample_dir}/${sample}.clean.narrowPeak"
    stable="${sample_dir}/${sample}.stable2of3.bed"

    robust="${ROBUST_DIR}/${sample}.robust.bed"
    qc="${ROBUST_DIR}/${sample}.robust.counts.tsv"

    echo "--------------------------------------------------"
    echo "Preparing robust peak regions: ${sample}"
    echo "--------------------------------------------------"

    if [[ ! -s "$clean" ]]; then
        echo "ERROR: missing or empty clean narrowPeak: $clean" >&2
        exit 1
    fi

    n_clean=$(wc -l < "$clean")
    n_stable="NA"
    source_type="full_clean"

    ###########################################################################
    # Low-depth samples:
    # stable2of3 already contains original robust peak coordinates.
    ###########################################################################

    if is_low_depth "$sample"; then
        source_type="stable2of3"

        if [[ ! -s "$stable" ]]; then
            echo "ERROR: low-depth sample is missing stable2of3 file: $stable" >&2
            exit 1
        fi

        n_stable=$(wc -l < "$stable")

        cut -f1-3 "$stable" \
        | sort -k1,1 -k2,2n \
        > "$robust"

    ###########################################################################
    # >=1M:
    # use original blacklist-clean MACS3 peak boundaries.
    ###########################################################################

    else

        cut -f1-3 "$clean" \
        | sort -k1,1 -k2,2n \
        > "$robust"

    fi

    n_robust=$(wc -l < "$robust")
    medlen=$(median_bed_length "$robust")
    meanlen=$(mean_bed_length "$robust")

    printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
        "$sample" \
        "$source_type" \
        "$n_clean" \
        "$n_stable" \
        "$n_robust" \
        "$medlen" \
        "$meanlen" \
        > "$qc"

    echo "DONE robust BED: ${sample}"
    echo "  Source               : ${source_type}"
    echo "  Robust peak regions  : ${n_robust}"
    echo "  Median length (bp)   : ${medlen}"
    echo "  Mean length (bp)     : ${meanlen}"
}

###############################################################################
# Stage 2:
# Chain-merge one matched Adenoma + Serrated cell type.
#
# bedtools merge intentionally performs transitive/chaining merge.
###############################################################################

chain_merge_one() {
    set -euo pipefail

    ct="$1"

    A="${ROBUST_DIR}/Adenoma_${ct}.robust.bed"
    S="${ROBUST_DIR}/Serrated_${ct}.robust.bed"

    FINAL="${OUT}/Adenoma_Serrated_${ct}.chaining_merge.bed"
    QC="${OUT}/${ct}.chaining_merge.counts.tsv"

    echo "=================================================="
    echo "Chaining merge: ${ct}"
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

    cat "$A" "$S" \
    | sort -k1,1 -k2,2n \
    | bedtools merge -d 0 \
    > "$FINAL"

    nA=$(wc -l < "$A")
    nS=$(wc -l < "$S")

    ###########################################################################
    # Overlap statistics:
    # number of peaks from each condition overlapping >=1 peak in the other.
    ###########################################################################

    nA_overlap=$(
        bedtools intersect -u -a "$A" -b "$S" | wc -l
    )

    nS_overlap=$(
        bedtools intersect -u -a "$S" -b "$A" | wc -l
    )

    pctA_overlap=$(
        awk -v x="$nA_overlap" -v n="$nA" \
            'BEGIN{if(n>0) printf "%.2f",100*x/n; else print "NA"}'
    )

    pctS_overlap=$(
        awk -v x="$nS_overlap" -v n="$nS" \
            'BEGIN{if(n>0) printf "%.2f",100*x/n; else print "NA"}'
    )

    nFinal=$(wc -l < "$FINAL")
    medlen=$(median_bed_length "$FINAL")
    meanlen=$(mean_bed_length "$FINAL")

    printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
        "$ct" \
        "$nA" \
        "$nS" \
        "$nA_overlap" \
        "$pctA_overlap" \
        "$nS_overlap" \
        "$pctS_overlap" \
        "$nFinal" \
        "$medlen" \
        "$meanlen" \
        > "$QC"

    echo "DONE chaining merge: ${ct}"
    echo "  Adenoma peaks                      : ${nA}"
    echo "  Serrated peaks                     : ${nS}"
    echo "  Adenoma overlapping Serrated       : ${nA_overlap} (${pctA_overlap}%)"
    echo "  Serrated overlapping Adenoma       : ${nS_overlap} (${pctS_overlap}%)"
    echo "  Regions after chain-merge          : ${nFinal}"
    echo "  Median merged length (bp)          : ${medlen}"
    echo "  Mean merged length (bp)            : ${meanlen}"
}

###############################################################################
# Export functions / variables for xargs
###############################################################################

export -f is_low_depth
export -f median_bed_length
export -f mean_bed_length
export -f prepare_robust_bed
export -f chain_merge_one

export ROOT
export PEAK_ROOT
export ROBUST_DIR
export OUT

###############################################################################
# Samples required
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
Normal_ABS
Normal_CT
Normal_EE
Normal_GOB
Normal_STM
Normal_TAC
Normal_TUF
EOF_SAMPLES

###############################################################################
# Clean only previous QC outputs from this workflow
###############################################################################

rm -f "${ROBUST_DIR}"/*.robust.counts.tsv
rm -f "${OUT}"/*.chaining_merge.counts.tsv

###############################################################################
# STAGE 1:
# Prepare robust original peak regions, max 8 samples simultaneously
###############################################################################

echo
echo "=============================================================="
echo "STAGE 1: PREPARE ROBUST ORIGINAL PEAK REGIONS"
echo "NO summit +/-250 bp conversion"
echo "Maximum parallel jobs: ${N_PARALLEL}"
echo "=============================================================="
echo

xargs -r -n 1 -P "$N_PARALLEL" \
    bash -c 'prepare_robust_bed "$1"' _ \
    < "${ROBUST_DIR}/samples_to_prepare.txt"

###############################################################################
# Robust peak preparation summary
###############################################################################

ROBUST_SUMMARY="${ROBUST_DIR}/robust_peak_region_summary.tsv"

echo -e \
"Sample\tPeak_source\tClean_peaks\tStable2of3_peaks\tRobust_peaks\tMedian_original_peak_length_bp\tMean_original_peak_length_bp" \
> "$ROBUST_SUMMARY"

for sample in \
    Adenoma_ABS Serrated_ABS \
    Adenoma_CT Serrated_CT \
    Adenoma_EE Serrated_EE \
    Adenoma_GOB Serrated_GOB \
    Adenoma_STM Serrated_STM \
    Adenoma_TAC Serrated_TAC \
    Adenoma_TUF Serrated_TUF \
    Adenoma_ASC Serrated_SSC \
    Normal_ABS Normal_CT Normal_EE Normal_GOB Normal_STM Normal_TAC Normal_TUF; do

    cat "${ROBUST_DIR}/${sample}.robust.counts.tsv" >> "$ROBUST_SUMMARY"
done

###############################################################################
# STAGE 2:
# Chain-merge seven matched Adenoma/Serrated cell types
###############################################################################

MATCHED_CELLTYPES=(ABS CT EE GOB STM TAC TUF)

echo
echo "=============================================================="
echo "STAGE 2: ADENOMA + SERRATED CHAINING MERGE"
echo "Parallel jobs: ${N_PARALLEL}"
echo "=============================================================="
echo

printf "%s\n" "${MATCHED_CELLTYPES[@]}" \
| xargs -r -n 1 -P "$N_PARALLEL" \
    bash -c 'chain_merge_one "$1"' _

###############################################################################
# Final summary
#
# Normal is included in the SAME row as the corresponding cell type.
#
# For ABS/CT/EE/GOB/STM/TAC/TUF:
#   - report Normal robust peak count/source/length statistics
#   - report Adenoma and Serrated peak counts/overlap statistics
#   - report Adenoma+Serrated chaining-merged region statistics
#
# ASC and SSC remain separate and have no matched Normal entry here.
###############################################################################

SUMMARY="${OUT}/chaining_merge_summary.tsv"

echo -e \
"CellType\tNum_peaks_Normal\tNormal_peak_source\tMedian_Normal_peak_length_bp\tMean_Normal_peak_length_bp\tNum_peaks_Adenoma\tNum_peaks_Serrated\tNum_Adenoma_peaks_overlapping_Serrated\tPct_Adenoma_peaks_overlapping_Serrated\tNum_Serrated_peaks_overlapping_Adenoma\tPct_Serrated_peaks_overlapping_Adenoma\tNum_regions_after_chain_merge\tMedian_chain_merged_region_length_bp\tMean_chain_merged_region_length_bp" \
> "$SUMMARY"

###############################################################################
# Seven matched epithelial cell types
###############################################################################

for ct in "${MATCHED_CELLTYPES[@]}"; do

    ###########################################################################
    # Normal statistics
    ###########################################################################

    normal_sample="Normal_${ct}"
    normal_bed="${ROBUST_DIR}/${normal_sample}.robust.bed"
    normal_qc="${ROBUST_DIR}/${normal_sample}.robust.counts.tsv"

    if [[ ! -s "$normal_bed" ]]; then
        echo "ERROR: missing Normal robust BED: $normal_bed" >&2
        exit 1
    fi

    if [[ ! -s "$normal_qc" ]]; then
        echo "ERROR: missing Normal robust QC: $normal_qc" >&2
        exit 1
    fi

    nNormal=$(wc -l < "$normal_bed")
    normal_source=$(awk -F'\t' 'NR==1 {print $2}' "$normal_qc")
    medNormal=$(median_bed_length "$normal_bed")
    meanNormal=$(mean_bed_length "$normal_bed")

    ###########################################################################
    # Adenoma + Serrated chaining-merge statistics
    #
    # Existing per-cell-type file columns:
    #   1 CellType
    #   2 Num_peaks_Adenoma
    #   3 Num_peaks_Serrated
    #   4 Num_Adenoma_peaks_overlapping_Serrated
    #   5 Pct_Adenoma_peaks_overlapping_Serrated
    #   6 Num_Serrated_peaks_overlapping_Adenoma
    #   7 Pct_Serrated_peaks_overlapping_Adenoma
    #   8 Num_regions_after_chain_merge
    #   9 Median_chain_merged_region_length_bp
    #  10 Mean_chain_merged_region_length_bp
    ###########################################################################

    disease_qc="${OUT}/${ct}.chaining_merge.counts.tsv"

    if [[ ! -s "$disease_qc" ]]; then
        echo "ERROR: missing chaining-merge QC: $disease_qc" >&2
        exit 1
    fi

    read -r \
        _ct \
        nA \
        nS \
        nA_overlap \
        pctA_overlap \
        nS_overlap \
        pctS_overlap \
        nFinal \
        medMerged \
        meanMerged \
        < "$disease_qc"

    printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
        "$ct" \
        "$nNormal" \
        "$normal_source" \
        "$medNormal" \
        "$meanNormal" \
        "$nA" \
        "$nS" \
        "$nA_overlap" \
        "$pctA_overlap" \
        "$nS_overlap" \
        "$pctS_overlap" \
        "$nFinal" \
        "$medMerged" \
        "$meanMerged" \
        >> "$SUMMARY"

done

###############################################################################
# ASC: Adenoma only, NOT merged with SSC
###############################################################################

ASC="${ROBUST_DIR}/Adenoma_ASC.robust.bed"

if [[ ! -s "$ASC" ]]; then
    echo "ERROR: missing $ASC" >&2
    exit 1
fi

nASC=$(wc -l < "$ASC")
medASC=$(median_bed_length "$ASC")
meanASC=$(mean_bed_length "$ASC")

printf "ASC\tNA\tNA\tNA\tNA\t%s\tNA\tNA\tNA\tNA\tNA\t%s\t%s\t%s\n" \
    "$nASC" \
    "$nASC" \
    "$medASC" \
    "$meanASC" \
    >> "$SUMMARY"

###############################################################################
# SSC: Serrated only, NOT merged with ASC
###############################################################################

SSC="${ROBUST_DIR}/Serrated_SSC.robust.bed"

if [[ ! -s "$SSC" ]]; then
    echo "ERROR: missing $SSC" >&2
    exit 1
fi

nSSC=$(wc -l < "$SSC")
medSSC=$(median_bed_length "$SSC")
meanSSC=$(mean_bed_length "$SSC")

printf "SSC\tNA\tNA\tNA\tNA\tNA\t%s\tNA\tNA\tNA\tNA\t%s\t%s\t%s\n" \
    "$nSSC" \
    "$nSSC" \
    "$medSSC" \
    "$meanSSC" \
    >> "$SUMMARY"

###############################################################################
# Final report
###############################################################################

echo
echo "=============================================================="
echo "CHAINING MERGE COMPLETE"
echo "=============================================================="
echo

echo "Merged coverage files:"
for ct in "${MATCHED_CELLTYPES[@]}"; do
    echo "  ${OUT}/Adenoma_Serrated_${ct}.chaining_merge.bed"
done

echo
echo "ASC and SSC were NOT merged:"
echo "  ASC: ${ROBUST_DIR}/Adenoma_ASC.robust.bed"
echo "  SSC: ${ROBUST_DIR}/Serrated_SSC.robust.bed"

echo
echo "Final combined Normal + Adenoma/Serrated summary:"
column -t "$SUMMARY"

echo
echo "Robust peak summary:"
echo "  $ROBUST_SUMMARY"

echo
echo "Combined chaining-merge summary:"
echo "  $SUMMARY"
