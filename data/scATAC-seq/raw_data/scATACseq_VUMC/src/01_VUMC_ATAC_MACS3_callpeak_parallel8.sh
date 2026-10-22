#!/usr/bin/env bash

set -euo pipefail

###############################################################################
# VUMC scATAC pseudo-bulk peak calling with MACS3
#
# Strategy
# --------
# 1. Process every cell type x condition BAM independently.
#    Count unique cell barcodes in each original split BAM.
# 2. Run up to 8 BAMs in parallel.
# 3. Filter:
#       properly paired
#       MAPQ >= 30
#       unmapped/mate-unmapped removed
#       secondary removed
#       supplementary removed
#       QC-fail removed
#       PCR duplicates removed
#       canonical chromosomes chr1-22 + chrX only
# 4. MACS3 callpeak:
#       BAMPE
#       q = 0.01
#       --keep-dup all
#       --nolambda
#       --call-summits
#       --cutoff-analysis
# 5. Remove hg38 blacklist.
# 6. For <1M usable fragments:
#       80% random subsampling x 3
#       seeds = 101, 202, 303
#       call peaks independently
#       identify full-data peaks supported in >=2/3 runs
#
# IMPORTANT:
#   Serrated_TAC.sorted.bam is intentionally skipped because
#   Serrated_TAC.bam already represents the same sample.
###############################################################################


###############################################################################
# CONFIGURATION
###############################################################################

BAM_DIR="/data/l2_bioinfo1/scATAC-seq/scATACseq_VUMC/split_bam"

PROJECT_DIR="/data/l2_bioinfo1/scATAC-seq/scATACseq_VUMC"

OUT_ROOT="${PROJECT_DIR}/MACS3_callpeak_q01"

# Put hg38-blacklist.v2.bed here, or change this path.
BLACKLIST="${PROJECT_DIR}/reference/hg38-blacklist.v2.bed"

# Number of BAM files processed simultaneously
N_PARALLEL=8

# One CPU per BAM job.
# This gives approximately 8 CPUs total.
SAMTOOLS_THREADS=1

MAPQ=30

# Cell barcode tags.
# 10x/Cell Ranger BAMs normally use CB:Z:<corrected_cell_barcode>.
# If CB is absent, fall back to CR:Z:<raw_cell_barcode>.
CELL_BARCODE_TAG_PRIMARY="CB"
CELL_BARCODE_TAG_FALLBACK="CR"

# Remove:
#   0x4   unmapped
#   0x8   mate unmapped
#   0x100 secondary
#   0x200 QC fail
#   0x400 duplicate
#   0x800 supplementary
#
# 4 + 8 + 256 + 512 + 1024 + 2048 = 3852?
#
# NOTE:
# We also require -f 2 (properly paired).
#
# To stay compatible with the usual ATAC filtering:
EXCLUDE_FLAGS=3844

LOW_DEPTH_THRESHOLD=1000000

# Low-depth stability analysis
RUN_STABILITY=1
STABILITY_PERCENT=80
STABILITY_MIN_SUPPORT=2

STABILITY_SEEDS=(101 202 303)


###############################################################################
# CHECK SOFTWARE
###############################################################################

for cmd in samtools macs3 bedtools awk sort xargs; do

    if ! command -v "$cmd" >/dev/null 2>&1; then
        echo "ERROR: $cmd not found in PATH"
        exit 1
    fi

done


###############################################################################
# CHECK INPUT
###############################################################################

if [[ ! -d "$BAM_DIR" ]]; then
    echo "ERROR: BAM directory does not exist:"
    echo "$BAM_DIR"
    exit 1
fi

if [[ ! -s "$BLACKLIST" ]]; then
    echo
    echo "ERROR: hg38 blacklist file not found:"
    echo "$BLACKLIST"
    echo
    echo "Please place hg38-blacklist.v2.bed there"
    echo "or modify BLACKLIST= in this script."
    exit 1
fi


###############################################################################
# OUTPUT DIRECTORIES
###############################################################################

mkdir -p \
    "${OUT_ROOT}/clean_bam" \
    "${OUT_ROOT}/peaks" \
    "${OUT_ROOT}/qc" \
    "${OUT_ROOT}/logs" \
    "${OUT_ROOT}/ready_full" \
    "${OUT_ROOT}/ready_robust"


###############################################################################
# RECORD SOFTWARE VERSIONS
###############################################################################

{
    echo "Date: $(date)"
    echo
    echo "samtools:"
    samtools --version | head -2
    echo
    echo "MACS3:"
    macs3 --version
    echo
    echo "bedtools:"
    bedtools --version

} > "${OUT_ROOT}/software_versions.txt"


###############################################################################
# COUNT UNIQUE CELL BARCODES
#
# Cell count is calculated from the ORIGINAL split BAM, before MAPQ/duplicate
# filtering, because the goal is to report the number of cells contributing
# to the pseudo-bulk BAM.
#
# Primary tag:
#     CB:Z:<corrected cell barcode>
#
# Fallback:
#     CR:Z:<raw cell barcode>
#
# If neither tag is present, return:
#     Num_cells = NA
#     Barcode_tag = NA
###############################################################################

count_unique_barcodes_for_tag() {

    local bam="$1"
    local tag="$2"

    samtools view \
        -@ "${SAMTOOLS_THREADS}" \
        "$bam" \
    | awk -v tag="$tag" '
        BEGIN {
            prefix = tag ":Z:"
        }
        {
            for (i=12; i<=NF; i++) {
                if (index($i,prefix) == 1) {
                    bc = substr($i,length(prefix)+1)
                    seen[bc] = 1
                    found = 1
                    break
                }
            }
        }
        END {
            if (!found) {
                print 0
            } else {
                n=0
                for (bc in seen)
                    n++
                print n
            }
        }
    '
}


count_cells_from_bam() {

    local bam="$1"

    local n_cb
    local n_cr

    n_cb=$(
        count_unique_barcodes_for_tag \
            "$bam" \
            "${CELL_BARCODE_TAG_PRIMARY}"
    )

    if [[ "$n_cb" -gt 0 ]]; then

        printf "%s\t%s\n" \
            "$n_cb" \
            "${CELL_BARCODE_TAG_PRIMARY}"

        return 0
    fi

    n_cr=$(
        count_unique_barcodes_for_tag \
            "$bam" \
            "${CELL_BARCODE_TAG_FALLBACK}"
    )

    if [[ "$n_cr" -gt 0 ]]; then

        printf "%s\t%s\n" \
            "$n_cr" \
            "${CELL_BARCODE_TAG_FALLBACK}"

        return 0
    fi

    printf "NA\tNA\n"
}


###############################################################################
# PROCESS ONE BAM
###############################################################################

process_one() {

    set -euo pipefail

    bam="$1"

    sample=$(basename "$bam" .bam)

    echo
    echo "=============================================================="
    echo "START: ${sample}"
    echo "BAM:   ${bam}"
    echo "=============================================================="

    sample_dir="${OUT_ROOT}/peaks/${sample}"
    log="${OUT_ROOT}/logs/${sample}.pipeline.log"

    mkdir -p "$sample_dir"

    ###########################################################################
    # Check BAM
    ###########################################################################

    samtools quickcheck -v "$bam"

    ###########################################################################
    # Count cells represented in this original cell type x condition BAM
    ###########################################################################

    echo "Counting unique cell barcodes..." >> "$log"

    IFS=$'\t' read -r num_cells barcode_tag < <(
        count_cells_from_bam "$bam"
    )

    echo "Num cells           : ${num_cells}" >> "$log"
    echo "Cell barcode tag    : ${barcode_tag}" >> "$log"

    # Index source BAM if necessary
    if [[ ! -f "${bam}.bai" && ! -f "${bam%.bam}.bai" ]]; then

        echo "Indexing original BAM..." >> "$log"

        samtools index \
            -@ "${SAMTOOLS_THREADS}" \
            "$bam"

    fi


    ###########################################################################
    # Determine chromosome naming style
    ###########################################################################

    if samtools idxstats "$bam" | cut -f1 | grep -qx "chr1"; then

        chrom_style="chr"

        chroms=()

        for c in {1..22} X; do
            chroms+=("chr${c}")
        done

    elif samtools idxstats "$bam" | cut -f1 | grep -qx "1"; then

        chrom_style="nochr"

        chroms=()

        for c in {1..22} X; do
            chroms+=("${c}")
        done

    else

        echo "ERROR: cannot determine chromosome naming style for $bam"
        exit 1

    fi


    ###########################################################################
    # Check blacklist chromosome style
    ###########################################################################

    blacklist_first_chr=$(
        awk '
            NF >= 3 && $1 !~ /^#/ {
                print $1
                exit
            }
        ' "$BLACKLIST"
    )

    if [[ "$chrom_style" == "chr" && "$blacklist_first_chr" != chr* ]]; then

        echo "ERROR: BAM uses chr1 but blacklist appears not to use chr prefix."
        exit 1

    fi

    if [[ "$chrom_style" == "nochr" && "$blacklist_first_chr" == chr* ]]; then

        echo "ERROR: BAM does not use chr prefix but blacklist does."
        exit 1

    fi


    ###########################################################################
    # Clean BAM
    #
    # -f 2:
    #     properly paired
    #
    # -F 3844:
    #     remove unmapped, mate-unmapped,
    #     secondary, QC-fail, duplicate,
    #     supplementary alignments
    #
    # -q 30:
    #     MAPQ >= 30
    #
    # canonical:
    #     chr1-22 + chrX
    ###########################################################################

    clean_bam="${OUT_ROOT}/clean_bam/${sample}.clean.bam"

    echo "Filtering BAM..." >> "$log"

    samtools view \
        -@ "${SAMTOOLS_THREADS}" \
        -b \
        -f 2 \
        -F "${EXCLUDE_FLAGS}" \
        -q "${MAPQ}" \
        "$bam" \
        "${chroms[@]}" \
    |
    samtools sort \
        -@ "${SAMTOOLS_THREADS}" \
        -o "$clean_bam" \
        -

    samtools index \
        -@ "${SAMTOOLS_THREADS}" \
        "$clean_bam"


    ###########################################################################
    # Count usable fragments
    ###########################################################################

    reads=$(samtools view -c "$clean_bam")

    fragments=$(( reads / 2 ))

    if (( fragments >= LOW_DEPTH_THRESHOLD )); then
        depth_class="GE1M"
    else
        depth_class="LT1M"
    fi

    echo "Usable paired reads : $reads" >> "$log"
    echo "Usable fragments    : $fragments" >> "$log"
    echo "Depth class         : $depth_class" >> "$log"


    ###########################################################################
    # Primary MACS3 callpeak
    ###########################################################################

    echo "Running MACS3 callpeak..." >> "$log"

    macs3 callpeak \
        -t "$clean_bam" \
        -f BAMPE \
        -g hs \
        -n "$sample" \
        --outdir "$sample_dir" \
        -q 0.01 \
        --keep-dup all \
        --nolambda \
        --call-summits \
        --cutoff-analysis \
        > "${sample_dir}/${sample}.macs3.log" \
        2>&1


    raw_peak="${sample_dir}/${sample}_peaks.narrowPeak"

    if [[ ! -f "$raw_peak" ]]; then

        echo "ERROR: MACS3 failed to generate:"
        echo "$raw_peak"
        exit 1

    fi


    ###########################################################################
    # Blacklist removal
    ###########################################################################

    clean_peak="${sample_dir}/${sample}.clean.narrowPeak"

    bedtools intersect \
        -v \
        -a "$raw_peak" \
        -b "$BLACKLIST" \
        > "$clean_peak"


    ###########################################################################
    # BED3 version
    #
    # Useful later for merging VUMC + GSE201336.
    ###########################################################################

    awk '
        BEGIN {
            OFS="\t"
        }
        {
            print $1,$2,$3
        }
    ' "$clean_peak" \
        > "${sample_dir}/${sample}.clean.bed"


    ###########################################################################
    # Peak statistics
    ###########################################################################

    n_raw=$(
        awk 'END {print NR+0}' "$raw_peak"
    )

    n_clean=$(
        awk 'END {print NR+0}' "$clean_peak"
    )

    total_peak_bp=$(
        awk '
            {
                sum += ($3-$2)
            }
            END {
                print sum+0
            }
        ' "$clean_peak"
    )

    mean_peak_width=$(
        awk '
            {
                sum += ($3-$2)
                n++
            }
            END {
                if(n > 0)
                    printf "%.2f",sum/n
                else
                    printf "NA"
            }
        ' "$clean_peak"
    )


    ###########################################################################
    # LOW-DEPTH STABILITY ANALYSIS
    #
    # Only if usable fragments < 1M.
    #
    # Three independent 80% subsamples:
    #
    # seed 101
    # seed 202
    # seed 303
    #
    # A full-data peak is considered stable when it overlaps a peak in
    # >=2/3 subsampled datasets.
    #
    # Require >=25% of the full peak to overlap the resampled peak.
    ###########################################################################

    stable_n="NA"

    if [[ "$depth_class" == "LT1M" && "$RUN_STABILITY" -eq 1 ]]; then

        echo "Running low-depth stability analysis..." >> "$log"

        stability_dir="${sample_dir}/stability80"

        mkdir -p "$stability_dir"

        full_bed4="${stability_dir}/${sample}.full.clean.bed4"

        awk '
            BEGIN {
                OFS="\t"
            }
            {
                print $1,$2,$3,$4
            }
        ' "$clean_peak" \
            > "$full_bed4"


        #######################################################################
        # Seed 101
        #######################################################################

        for seed in 101 202 303; do

            rep_dir="${stability_dir}/seed${seed}"

            mkdir -p "$rep_dir"

            sampled_bedpe="${rep_dir}/${sample}.80pct.seed${seed}.bedpe"

            macs3 randsample \
                -i "$clean_bam" \
                -f BAMPE \
                -p "${STABILITY_PERCENT}" \
                --seed "${seed}" \
                -o "$sampled_bedpe" \
                > "${rep_dir}/randsample.log" \
                2>&1


            macs3 callpeak \
                -t "$sampled_bedpe" \
                -f BEDPE \
                -g hs \
                -n "${sample}.80pct.seed${seed}" \
                --outdir "$rep_dir" \
                -q 0.01 \
                --keep-dup all \
                --nolambda \
                --call-summits \
                > "${rep_dir}/macs3.log" \
                2>&1


            rep_raw="${rep_dir}/${sample}.80pct.seed${seed}_peaks.narrowPeak"

            rep_clean="${rep_dir}/${sample}.80pct.seed${seed}.clean.narrowPeak"

            bedtools intersect \
                -v \
                -a "$rep_raw" \
                -b "$BLACKLIST" \
                > "$rep_clean"


            ###################################################################
            # Binary support for each full-data peak
            ###################################################################

            bedtools intersect \
                -c \
                -f 0.25 \
                -a "$full_bed4" \
                -b "$rep_clean" \
            |
            awk '
                {
                    if($NF > 0)
                        print 1
                    else
                        print 0
                }
            ' \
                > "${stability_dir}/support.seed${seed}.txt"

        done


        #######################################################################
        # Stable >=2/3
        #
        # Output BED5:
        #
        # chr
        # start
        # end
        # original_peak_name
        # number_of_supporting_subsamples
        #######################################################################

        stable_peak="${sample_dir}/${sample}.stable2of3.bed"

        paste \
            "$full_bed4" \
            "${stability_dir}/support.seed101.txt" \
            "${stability_dir}/support.seed202.txt" \
            "${stability_dir}/support.seed303.txt" \
        |
        awk -v min_support="${STABILITY_MIN_SUPPORT}" '
            BEGIN {
                OFS="\t"
            }
            {
                support=$5+$6+$7

                if(support >= min_support)
                    print $1,$2,$3,$4,support
            }
        ' \
            > "$stable_peak"

        stable_n=$(
            awk 'END {print NR+0}' "$stable_peak"
        )

        echo "Stable >=2/3 peaks : ${stable_n}" >> "$log"

    fi


    ###########################################################################
    # Per-sample QC output
    ###########################################################################

    printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
        "$sample" \
        "$num_cells" \
        "$barcode_tag" \
        "$fragments" \
        "$depth_class" \
        "$n_raw" \
        "$n_clean" \
        "$stable_n" \
        "$mean_peak_width" \
        "$total_peak_bp" \
        > "${OUT_ROOT}/qc/${sample}.tsv"


    echo "DONE: ${sample}"
    echo "Cells: ${num_cells}"
    echo "Barcode tag: ${barcode_tag}"
    echo "Fragments: ${fragments}"
    echo "Clean peaks: ${n_clean}"

    if [[ "$stable_n" != "NA" ]]; then
        echo "Stable >=2/3 peaks: ${stable_n}"
    fi
}


###############################################################################
# EXPORT FUNCTION / VARIABLES FOR xargs WORKERS
###############################################################################

export -f count_unique_barcodes_for_tag
export -f count_cells_from_bam
export -f process_one

export BAM_DIR
export PROJECT_DIR
export OUT_ROOT
export BLACKLIST

export N_PARALLEL
export SAMTOOLS_THREADS

export MAPQ
export CELL_BARCODE_TAG_PRIMARY
export CELL_BARCODE_TAG_FALLBACK
export EXCLUDE_FLAGS
export LOW_DEPTH_THRESHOLD

export RUN_STABILITY
export STABILITY_PERCENT
export STABILITY_MIN_SUPPORT


###############################################################################
# FIND BAM FILES
#
# Skip:
#
#     *.sorted.bam
#
# because Serrated_TAC.sorted.bam duplicates Serrated_TAC.bam.
###############################################################################

echo
echo "=============================================================="
echo "VUMC ATAC MACS3 peak calling"
echo "=============================================================="
echo
echo "BAM directory : $BAM_DIR"
echo "Output        : $OUT_ROOT"
echo "Parallel jobs : $N_PARALLEL"
echo
echo "BAM files:"
echo

find "$BAM_DIR" \
    -maxdepth 1 \
    -type f \
    -name "*.bam" \
    ! -name "*.sorted.bam" \
    -print \
    | sort

echo
echo "Starting peak calling..."
echo


###############################################################################
# RUN 8 BAM FILES IN PARALLEL
###############################################################################

find "$BAM_DIR" \
    -maxdepth 1 \
    -type f \
    -name "*.bam" \
    ! -name "*.sorted.bam" \
    -print0 \
|
sort -z \
|
xargs \
    -0 \
    -r \
    -n 1 \
    -P "${N_PARALLEL}" \
    bash -c 'process_one "$1"' _


###############################################################################
# COMBINE QC RESULTS
###############################################################################

SUMMARY="${OUT_ROOT}/peak_calling_summary.tsv"

echo -e \
"Sample\tNum_cells\tBarcode_tag\tUsable_fragments\tDepth_class\tRaw_peaks\tBlacklist_clean_peaks\tStable2of3_peaks\tMean_peak_width\tTotal_peak_bp" \
    > "$SUMMARY"

cat "${OUT_ROOT}"/qc/*.tsv \
    | sort -k1,1 \
    >> "$SUMMARY"


###############################################################################
# PREPARE PEAK FILES FOR FINAL DISEASE PEAK MERGING
#
# Only:
#
# Adenoma + Serrated
#
# Cell types:
#
# ABS CT EE GOB STM TAC TUF
#
# Two versions:
#
# ready_full/
#     always use full-data blacklist-clean peaks
#
# ready_robust/
#     >=1M  -> full-data blacklist-clean peaks
#     <1M   -> stable2of3 peaks when available
#
# This means TUF will get special low-depth handling.
###############################################################################

TARGET_CELLTYPES=(ABS CT EE GOB STM TAC TUF)
TARGET_CONDITIONS=(Adenoma Serrated)

for condition in "${TARGET_CONDITIONS[@]}"; do

    for ct in "${TARGET_CELLTYPES[@]}"; do

        sample="${condition}_${ct}"

        sample_dir="${OUT_ROOT}/peaks/${sample}"

        full="${sample_dir}/${sample}.clean.narrowPeak"

        stable="${sample_dir}/${sample}.stable2of3.bed"


        #######################################################################
        # Full-data version
        #######################################################################

        if [[ -s "$full" ]]; then

            awk '
                BEGIN {
                    OFS="\t"
                }
                {
                    print $1,$2,$3
                }
            ' "$full" \
                > "${OUT_ROOT}/ready_full/${sample}.bed"

        else

            echo "WARNING: missing full peak file for ${sample}"

        fi


        #######################################################################
        # Robust version
        #######################################################################

        depth_class=$(
            awk -v s="$sample" '
                $1==s {
                    print $5
                }
            ' "$SUMMARY"
        )

        if [[ "$depth_class" == "LT1M" && -s "$stable" ]]; then

            awk '
                BEGIN {
                    OFS="\t"
                }
                {
                    print $1,$2,$3
                }
            ' "$stable" \
                > "${OUT_ROOT}/ready_robust/${sample}.bed"

            echo "${sample}: robust set = stable2of3"

        elif [[ -s "$full" ]]; then

            awk '
                BEGIN {
                    OFS="\t"
                }
                {
                    print $1,$2,$3
                }
            ' "$full" \
                > "${OUT_ROOT}/ready_robust/${sample}.bed"

            echo "${sample}: robust set = full"

        fi

    done

done


###############################################################################
# FINAL REPORT
###############################################################################

echo
echo "=============================================================="
echo "ALL PEAK CALLING FINISHED"
echo "=============================================================="
echo

column -t "$SUMMARY"

echo
echo "Full peak sets:"
echo "${OUT_ROOT}/ready_full"
echo
echo "Robust peak sets:"
echo "${OUT_ROOT}/ready_robust"
echo