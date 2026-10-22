#!/usr/bin/env bash
set -euo pipefail

###############################################################################
# GSE201336 / Wellcome Sanger ATAC peak-calling summary
#
# Input BED directory:
#   /data/l2_bioinfo1/scATAC-seq/scATACseq_GSE201336/
#   GSE201336_scATAC_peaks/peaks/per_group/
#
# Output:
#   peak_calling_summary.tsv
#
# Columns:
#   Sample
#   Num_cells
#   Usable_fragments
#   Depth_class
#   Raw_peaks
#   clean_peaks_rm_Blacklist
#   Stable_ge2of3
#   Stability
#   Mean_peak_width
#   Total_peak_bp
#
# Notes:
#   - Raw_peaks = number of intervals in the existing ArchR reproducible BED.
#   - clean_peaks_rm_Blacklist = peaks remaining after hg38 blacklist removal.
#   - Mean_peak_width and Total_peak_bp are calculated from blacklist-clean peaks.
#   - Usable_fragments and Depth_class are NA because they cannot be recovered
#     reliably from these final per-group BEDs.
#   - Stable_ge2of3 is NA because GSE201336 used ArchR reproducible pseudo-bulk
#     peak calling, not the VUMC 3 x 80% resampling workflow.
#   - Stability is recorded as ArchR_reproducible.
###############################################################################

export LC_ALL=C

ROOT="/data/l2_bioinfo1/scATAC-seq/scATACseq_GSE201336/GSE201336_scATAC_peaks"

BED_DIR="${ROOT}/peaks/per_group"

ARCHR_ROOT="${ROOT}/archr"
ARCHR_PROJECT="${ARCHR_ROOT}/proj_epithelial_lau"

# Written by 03_call_peaks.R after setwd(opt$archr)
GROUP_SIZE_CSV="${ARCHR_ROOT}/peak_calling_group_sizes.csv"

SUMMARY="${BED_DIR}/peak_calling_summary.tsv"
CELL_COUNTS="${BED_DIR}/cell_counts_condition_celltype.tsv"

CLEAN_DIR="${BED_DIR}/blacklist_clean"
mkdir -p "$CLEAN_DIR"

###############################################################################
# hg38 blacklist
###############################################################################

BLACKLIST_LOCAL="${ROOT}/reference/hg38-blacklist.v2.bed"
BLACKLIST_VUMC="/data/l2_bioinfo1/scATAC-seq/scATACseq_VUMC/reference/hg38-blacklist.v2.bed"

if [[ -s "$BLACKLIST_LOCAL" ]]; then
    BLACKLIST="$BLACKLIST_LOCAL"
elif [[ -s "$BLACKLIST_VUMC" ]]; then
    BLACKLIST="$BLACKLIST_VUMC"
else
    echo "ERROR: hg38 blacklist not found." >&2
    echo "Checked:" >&2
    echo "  $BLACKLIST_LOCAL" >&2
    echo "  $BLACKLIST_VUMC" >&2
    exit 1
fi

###############################################################################
# Software checks
###############################################################################

for cmd in awk sort bedtools Rscript wc column; do
    if ! command -v "$cmd" >/dev/null 2>&1; then
        echo "ERROR: required command not found: $cmd" >&2
        exit 1
    fi
done

if [[ ! -d "$BED_DIR" ]]; then
    echo "ERROR: BED directory not found: $BED_DIR" >&2
    exit 1
fi

###############################################################################
# Target samples
###############################################################################

SAMPLES=(
    Adenoma__ABS
    Adenoma__CT
    Adenoma__EE
    Adenoma__GOB
    Adenoma__STM
    Adenoma__TAC
    Normal__ABS
    Normal__CT
    Normal__EE
    Normal__GOB
    Normal__STM
    Normal__TAC
)

###############################################################################
# Obtain Condition x CellType cell counts.
#
# Preferred:
#   peak_calling_group_sizes.csv
#
# Fallback:
#   load proj_epithelial_lau and count Condition x CellType directly.
###############################################################################

echo "Obtaining cell counts..."

if [[ -s "$GROUP_SIZE_CSV" ]]; then

    echo "Using:"
    echo "  $GROUP_SIZE_CSV"

    Rscript - "$GROUP_SIZE_CSV" "$CELL_COUNTS" <<'RSCRIPT'
args <- commandArgs(trailingOnly = TRUE)

infile  <- args[1]
outfile <- args[2]

x <- read.csv(
    infile,
    stringsAsFactors = FALSE,
    check.names = FALSE
)

if (!all(c("Group", "nCells") %in% colnames(x))) {
    stop(
        "Expected columns Group and nCells in ",
        infile,
        "; observed: ",
        paste(colnames(x), collapse = ", ")
    )
}

out <- data.frame(
    Sample = as.character(x$Group),
    Num_cells = as.integer(x$nCells),
    stringsAsFactors = FALSE
)

write.table(
    out,
    file = outfile,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    col.names = TRUE
)
RSCRIPT

else

    echo "peak_calling_group_sizes.csv not found."
    echo "Falling back to ArchR project metadata:"
    echo "  $ARCHR_PROJECT"

    if [[ ! -d "$ARCHR_PROJECT" ]]; then
        echo "ERROR: ArchR project not found: $ARCHR_PROJECT" >&2
        exit 1
    fi

    Rscript - "$ARCHR_ROOT" "$ARCHR_PROJECT" "$CELL_COUNTS" <<'RSCRIPT'
args <- commandArgs(trailingOnly = TRUE)

archr_root  <- args[1]
project_dir <- args[2]
outfile     <- args[3]

suppressPackageStartupMessages(
    library(ArchR)
)

setwd(archr_root)

proj <- loadArchRProject(project_dir)

condition <- as.character(proj$Condition)
celltype  <- as.character(proj$CellTypeLau)

ok <- !is.na(condition) &
      !is.na(celltype) &
      nzchar(condition) &
      nzchar(celltype)

sample <- paste(
    condition[ok],
    celltype[ok],
    sep = "__"
)

tab <- sort(table(sample))

out <- data.frame(
    Sample = names(tab),
    Num_cells = as.integer(tab),
    stringsAsFactors = FALSE
)

write.table(
    out,
    file = outfile,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    col.names = TRUE
)
RSCRIPT

fi

###############################################################################
# Create summary header
###############################################################################

echo -e \
"Sample\tNum_cells\tUsable_fragments\tDepth_class\tRaw_peaks\tclean_peaks_rm_Blacklist\tStable_ge2of3\tStability\tMean_peak_width\tTotal_peak_bp" \
> "$SUMMARY"

###############################################################################
# Process each per-group BED
###############################################################################

for sample in "${SAMPLES[@]}"; do

    bed="${BED_DIR}/${sample}.narrowPeak.bed"
    clean="${CLEAN_DIR}/${sample}.blacklist_clean.bed"

    echo
    echo "=============================================================="
    echo "Sample: $sample"
    echo "BED   : $bed"
    echo "=============================================================="

    ###########################################################################
    # Cell count
    ###########################################################################

    num_cells=$(
        awk -F'\t' -v s="$sample" '
            NR > 1 && $1 == s {
                print $2
                exit
            }
        ' "$CELL_COUNTS"
    )

    [[ -n "$num_cells" ]] || num_cells="NA"

    ###########################################################################
    # Missing BED handling
    ###########################################################################

    if [[ ! -s "$bed" ]]; then

        echo "WARNING: missing or empty BED: $bed" >&2

        printf "%s\t%s\tNA\tNA\tNA\tNA\tNA\tMissing_BED\tNA\tNA\n" \
            "$sample" \
            "$num_cells" \
            >> "$SUMMARY"

        continue
    fi

    ###########################################################################
    # Raw peaks
    ###########################################################################

    raw_peaks=$(
        awk 'NF >= 3 {n++} END {print n+0}' "$bed"
    )

    ###########################################################################
    # Blacklist removal
    ###########################################################################

    bedtools intersect \
        -v \
        -a "$bed" \
        -b "$BLACKLIST" \
        > "$clean"

    clean_peaks=$(
        awk 'NF >= 3 {n++} END {print n+0}' "$clean"
    )

    ###########################################################################
    # Peak width statistics from blacklist-clean peaks
    ###########################################################################

    mean_peak_width=$(
        awk '
            NF >= 3 {
                w = $3-$2
                sum += w
                n++
            }
            END {
                if (n > 0)
                    printf "%.2f", sum/n
                else
                    printf "NA"
            }
        ' "$clean"
    )

    total_peak_bp=$(
        awk '
            NF >= 3 {
                sum += ($3-$2)
            }
            END {
                print sum+0
            }
        ' "$clean"
    )

    ###########################################################################
    # Summary row
    ###########################################################################

    printf "%s\t%s\tNA\tNA\t%s\t%s\tNA\tArchR_reproducible\t%s\t%s\n" \
        "$sample" \
        "$num_cells" \
        "$raw_peaks" \
        "$clean_peaks" \
        "$mean_peak_width" \
        "$total_peak_bp" \
        >> "$SUMMARY"

done

###############################################################################
# Final report
###############################################################################

echo
echo "=============================================================="
echo "GSE201336 / WELLCOME SANGER PEAK SUMMARY COMPLETE"
echo "=============================================================="
echo

echo "Blacklist used:"
echo "  $BLACKLIST"

echo
echo "Cell-count table:"
echo "  $CELL_COUNTS"

echo
echo "Blacklist-clean BEDs:"
echo "  $CLEAN_DIR"

echo
echo "Peak-calling summary:"
echo "  $SUMMARY"

echo
column -t "$SUMMARY"
