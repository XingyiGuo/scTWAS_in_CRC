#!/usr/bin/env bash
set -euo pipefail

###############################################################################
# Add VUMC cell numbers to existing peak_calling_summary.tsv
#
# No MACS3 peak calling is rerun.
#
# Cell number:
#   number of unique cell barcodes in each ORIGINAL split BAM.
#
# Barcode tag:
#   primary  = CB:Z:  (corrected cell barcode)
#   fallback = CR:Z:  (raw barcode)
#
# Existing summary:
#   /data/l2_bioinfo1/scATAC-seq/scATACseq_VUMC/MACS3_callpeak_q01/
#       peak_calling_summary.tsv
#
# Output:
#   1) cell_count_summary.tsv
#   2) peak_calling_summary.tsv updated with Num_cells
#   3) one-time backup:
#      peak_calling_summary.before_cellcounts.tsv
###############################################################################

export LC_ALL=C

BAM_DIR="/data/l2_bioinfo1/scATAC-seq/scATACseq_VUMC/split_bam"
OUT_ROOT="/data/l2_bioinfo1/scATAC-seq/scATACseq_VUMC/MACS3_callpeak_q01"

SUMMARY="${OUT_ROOT}/peak_calling_summary.tsv"
BACKUP="${OUT_ROOT}/peak_calling_summary.before_cellcounts.tsv"
CELL_SUMMARY="${OUT_ROOT}/cell_count_summary.tsv"

N_PARALLEL=8
SAMTOOLS_THREADS=1

###############################################################################
# Checks
###############################################################################

for cmd in samtools awk sort xargs; do
    if ! command -v "$cmd" >/dev/null 2>&1; then
        echo "ERROR: required command not found: $cmd" >&2
        exit 1
    fi
done

if [[ ! -d "$BAM_DIR" ]]; then
    echo "ERROR: BAM directory not found: $BAM_DIR" >&2
    exit 1
fi

if [[ ! -s "$SUMMARY" ]]; then
    echo "ERROR: existing peak summary not found: $SUMMARY" >&2
    exit 1
fi

###############################################################################
# Count unique cell barcodes in one BAM.
#
# One full BAM scan is required because the number of unique cell barcodes
# cannot be inferred from the pseudo-bulk read count.
#
# We count CB directly; if no CB occurs anywhere in the BAM, repeat using CR.
###############################################################################

count_tag() {
    local bam="$1"
    local tag="$2"

    samtools view -@ "$SAMTOOLS_THREADS" "$bam" \
    | awk -v tag="$tag" '
        BEGIN {
            prefix = tag ":Z:"
        }
        {
            for (i=12; i<=NF; i++) {
                if (index($i,prefix) == 1) {
                    bc = substr($i, length(prefix)+1)
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
                n = 0
                for (bc in seen)
                    n++
                print n
            }
        }
    '
}

count_one() {
    set -euo pipefail

    bam="$1"
    sample=$(basename "$bam" .bam)

    # Skip the known duplicate file if ever passed explicitly.
    if [[ "$sample" == "Serrated_TAC.sorted" ]]; then
        return 0
    fi

    echo "Counting cells: ${sample}" >&2

    n_cb=$(count_tag "$bam" "CB")

    if [[ "$n_cb" -gt 0 ]]; then
        num_cells="$n_cb"
        barcode_tag="CB"
    else
        n_cr=$(count_tag "$bam" "CR")

        if [[ "$n_cr" -gt 0 ]]; then
            num_cells="$n_cr"
            barcode_tag="CR"
        else
            num_cells="NA"
            barcode_tag="NA"
        fi
    fi

    printf "%s\t%s\t%s\n" \
        "$sample" "$num_cells" "$barcode_tag"
}

export -f count_tag
export -f count_one
export SAMTOOLS_THREADS

###############################################################################
# Count all samples in parallel
###############################################################################

TMP_DIR=$(mktemp -d "${OUT_ROOT}/cellcount_tmp.XXXXXX")
trap 'rm -rf "$TMP_DIR"' EXIT

find "$BAM_DIR" \
    -maxdepth 1 \
    -type f \
    -name "*.bam" \
    ! -name "*.sorted.bam" \
    -print0 \
| sort -z \
| xargs -0 -r -n 1 -P "$N_PARALLEL" \
    bash -c 'count_one "$1"' _ \
> "${TMP_DIR}/counts.tsv"

{
    echo -e "Sample\tNum_cells\tBarcode_tag"
    sort -k1,1 "${TMP_DIR}/counts.tsv"
} > "$CELL_SUMMARY"

###############################################################################
# Validate that every row in peak_calling_summary.tsv has a cell-count row.
###############################################################################

missing=$(
    awk -F'\t' '
        NR==FNR {
            if (FNR > 1)
                have[$1]=1
            next
        }
        FNR > 1 && !($1 in have) {
            print $1
        }
    ' "$CELL_SUMMARY" "$SUMMARY"
)

if [[ -n "$missing" ]]; then
    echo "ERROR: these peak-summary samples have no cell-count result:" >&2
    echo "$missing" >&2
    exit 1
fi

###############################################################################
# Back up original summary once.
###############################################################################

if [[ ! -e "$BACKUP" ]]; then
    cp -p "$SUMMARY" "$BACKUP"
    echo "Backup created:"
    echo "  $BACKUP"
else
    echo "Backup already exists; leaving it unchanged:"
    echo "  $BACKUP"
fi

###############################################################################
# Add/replace Num_cells column.
#
# Handles either:
#   old summary without Num_cells
# or
#   a summary that already has Num_cells from a previous run.
###############################################################################

UPDATED="${TMP_DIR}/peak_calling_summary.updated.tsv"

awk -F'\t' '
    BEGIN {
        OFS="\t"
    }

    NR==FNR {
        if (FNR > 1) {
            cells[$1]=$2
        }
        next
    }

    FNR==1 {
        # If Num_cells is already column 2, preserve the current structure.
        if ($2 == "Num_cells") {
            print
            already=1
        } else {
            printf "%s\tNum_cells", $1
            for (i=2; i<=NF; i++)
                printf "\t%s", $i
            printf "\n"
            already=0
        }
        next
    }

    {
        if (already) {
            $2=cells[$1]
            print
        } else {
            printf "%s\t%s", $1, cells[$1]
            for (i=2; i<=NF; i++)
                printf "\t%s", $i
            printf "\n"
        }
    }
' "$CELL_SUMMARY" "$SUMMARY" > "$UPDATED"

mv "$UPDATED" "$SUMMARY"

###############################################################################
# Final report
###############################################################################

echo
echo "=============================================================="
echo "CELL COUNTS ADDED TO PEAK CALLING SUMMARY"
echo "=============================================================="
echo

echo "Cell-count table:"
column -t "$CELL_SUMMARY" 2>/dev/null || cat "$CELL_SUMMARY"

echo
echo "Updated peak summary:"
column -t "$SUMMARY" 2>/dev/null || cat "$SUMMARY"

echo
echo "Files:"
echo "  Cell counts : $CELL_SUMMARY"
echo "  Updated     : $SUMMARY"
echo "  Backup      : $BACKUP"
