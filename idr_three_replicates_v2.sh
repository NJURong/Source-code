#!/bin/bash
# =============================================================================
# 3-Replicate IDR Analysis Pipeline (CUT&Tag / ChIP-seq)
# 
# Strategy: Run all 3 pairwise IDR comparisons, then take consensus peaks
#           present in ≥2/3 comparisons via bedtools intersect (overlap-based)
#
# Fixes vs original:
#   1. bedtools intersect replaces uniq -c (handles coordinate near-duplicates)
#   2. --rank changed to p.value (better for CUT&Tag data)
#   3. --soft-idr-threshold used to retain all peaks with IDR scores
#   4. Input file validation before running
#   5. Temp file names scoped by group to avoid collisions
#   6. Full narrowPeak format preserved in consensus output
#   7. Per-pair peak counts logged to summary table
# =============================================================================

set -euo pipefail

# =============================================================================
# Configuration — edit these paths as needed
# =============================================================================
BASE_DIR="/mnt/i/LQ/WR_Chip_seq"
OUT_DIR="${BASE_DIR}/results"
IDR_DIR="${OUT_DIR}/07_idr_three_replicates_v2"
PEAK_DIR="${OUT_DIR}/06_peaks"

GROUPS="A B C D"
THRESHOLDS="0.01 0.05"
unset GROUPS
GROUPS="A B C D"
unset THRESHOLDS
THRESHOLDS="0.01 0.05"

# IDR rank metric: p.value recommended for CUT&Tag; signal.value for ChIP-seq
RANK="p.value"

# =============================================================================
# Environment setup
# =============================================================================
source /home/lq/miniconda3/etc/profile.d/conda.sh
conda activate cuttag

# Verify required tools
for tool in idr bedtools; do
    if ! command -v $tool &>/dev/null; then
        echo "ERROR: '$tool' not found in PATH. Please check your conda environment."
        exit 1
    fi
done

# Create output directories
for threshold in $THRESHOLDS; do
    for group in $GROUPS; do
        mkdir -p "${IDR_DIR}/${threshold}/${group}"
    done
done

# Summary table
SUMMARY="${IDR_DIR}/idr_summary.tsv"
echo -e "Group\tThreshold\tRep1_Peaks\tRep2_Peaks\tRep3_Peaks\tIDR_1v2\tIDR_1v3\tIDR_2v3\tConsensus_2of3\tConsensus_3of3" \
    > "${SUMMARY}"

echo "========================================"
echo "  3-Replicate IDR Analysis v2"
echo "  Rank metric : ${RANK}"
echo "  Started     : $(date)"
echo "========================================"

# =============================================================================
# Main function
# =============================================================================
run_three_rep_idr() {
    local group=$1
    local threshold=$2
    local outdir="${IDR_DIR}/${threshold}/${group}"

    local rep1="${PEAK_DIR}/${group}1_peaks.narrowPeak"
    local rep2="${PEAK_DIR}/${group}2_peaks.narrowPeak"
    local rep3="${PEAK_DIR}/${group}3_peaks.narrowPeak"

    # ------------------------------------------------------------------
    # 1. Input validation
    # ------------------------------------------------------------------
    echo ""
    echo "  [${group}] Threshold=${threshold} — Validating inputs..."
    for f in "$rep1" "$rep2" "$rep3"; do
        if [[ ! -f "$f" ]]; then
            echo "  ERROR: File not found: $f"
            exit 1
        fi
        if [[ ! -s "$f" ]]; then
            echo "  ERROR: File is empty: $f"
            exit 1
        fi
    done

    local n1=$(wc -l < "$rep1")
    local n2=$(wc -l < "$rep2")
    local n3=$(wc -l < "$rep3")
    echo "  [${group}] Input peaks — rep1: ${n1}  rep2: ${n2}  rep3: ${n3}"

    # ------------------------------------------------------------------
    # 2. Pairwise IDR
    #    --soft-idr-threshold: output ALL peaks with IDR scores annotated;
    #    then filter by column 12 (globalIDR) manually for flexibility.
    # ------------------------------------------------------------------
    run_idr_pair() {
        local s1=$1 s2=$2 label=$3
        echo "  [${group}] IDR: ${label}"
        idr --samples "$s1" "$s2" \
            --input-file-type narrowPeak \
            --rank "${RANK}" \
            --output-file "${outdir}/${group}_${label}_all.narrowPeak" \
            --soft-idr-threshold "${threshold}" \
            --plot \
            --log-output-file "${outdir}/${group}_${label}.log" 2>&1 \
            | grep -v "^$" | sed "s/^/    /"

        # Filter by global IDR (column 12 = -log10(globalIDR), convert threshold)
        # globalIDR stored as -log10 value; threshold 0.05 → 12.996; 0.01 → 20
        local idr_col12_cutoff
        idr_col12_cutoff=$(python3 -c "import math; print(-math.log10(${threshold}))")

        awk -v cut="$idr_col12_cutoff" 'BEGIN{OFS="\t"} $12 >= cut {print}' \
            "${outdir}/${group}_${label}_all.narrowPeak" \
            > "${outdir}/${group}_${label}_idr.narrowPeak"

        wc -l < "${outdir}/${group}_${label}_idr.narrowPeak"
    }

    local n_1v2 n_1v3 n_2v3
    n_1v2=$(run_idr_pair "$rep1" "$rep2" "1vs2")
    n_1v3=$(run_idr_pair "$rep1" "$rep3" "1vs3")
    n_2v3=$(run_idr_pair "$rep2" "$rep3" "2vs3")

    echo "  [${group}] IDR peaks — 1v2: ${n_1v2}  1v3: ${n_1v3}  2v3: ${n_2v3}"

    # ------------------------------------------------------------------
    # 3. Extract BED3 for each pair (temp files scoped by group+threshold)
    # ------------------------------------------------------------------
    local tmp_prefix="${outdir}/${group}_${threshold}_tmp"

    cut -f1-3 "${outdir}/${group}_1vs2_idr.narrowPeak" \
        | sort -k1,1 -k2,2n > "${tmp_prefix}_1v2.bed"
    cut -f1-3 "${outdir}/${group}_1vs3_idr.narrowPeak" \
        | sort -k1,1 -k2,2n > "${tmp_prefix}_1v3.bed"
    cut -f1-3 "${outdir}/${group}_2vs3_idr.narrowPeak" \
        | sort -k1,1 -k2,2n > "${tmp_prefix}_2v3.bed"

    # ------------------------------------------------------------------
    # 4. Consensus via bedtools intersect (overlap-based, not string match)
    # ------------------------------------------------------------------

    # 2-of-3 consensus: union of pairwise overlaps, then merge
    bedtools intersect -a "${tmp_prefix}_1v2.bed" -b "${tmp_prefix}_1v3.bed" -u \
        > "${tmp_prefix}_ol_12_13.bed"
    bedtools intersect -a "${tmp_prefix}_1v2.bed" -b "${tmp_prefix}_2v3.bed" -u \
        > "${tmp_prefix}_ol_12_23.bed"
    bedtools intersect -a "${tmp_prefix}_1v3.bed" -b "${tmp_prefix}_2v3.bed" -u \
        > "${tmp_prefix}_ol_13_23.bed"

    cat "${tmp_prefix}_ol_12_13.bed" \
        "${tmp_prefix}_ol_12_23.bed" \
        "${tmp_prefix}_ol_13_23.bed" \
        | sort -k1,1 -k2,2n \
        | bedtools merge -i stdin \
        > "${outdir}/${group}_consensus_2of3.bed"

    # 3-of-3 consensus: overlap across all three pairs
    bedtools intersect -a "${tmp_prefix}_ol_12_13.bed" \
                       -b "${tmp_prefix}_2v3.bed" -u \
        | sort -k1,1 -k2,2n \
        | bedtools merge -i stdin \
        > "${outdir}/${group}_consensus_3of3.bed"

    # ------------------------------------------------------------------
    # 5. Full narrowPeak consensus (retain signal info from 1vs2 as anchor)
    # ------------------------------------------------------------------
    bedtools intersect \
        -a "${outdir}/${group}_1vs2_idr.narrowPeak" \
        -b "${outdir}/${group}_consensus_2of3.bed" \
        -u \
        | sort -k1,1 -k2,2n \
        > "${outdir}/${group}_consensus_2of3.narrowPeak"

    bedtools intersect \
        -a "${outdir}/${group}_1vs2_idr.narrowPeak" \
        -b "${outdir}/${group}_consensus_3of3.bed" \
        -u \
        | sort -k1,1 -k2,2n \
        > "${outdir}/${group}_consensus_3of3.narrowPeak"

    # ------------------------------------------------------------------
    # 6. Clean up temp files
    # ------------------------------------------------------------------
    rm -f "${tmp_prefix}"*.bed

    # ------------------------------------------------------------------
    # 7. Report
    # ------------------------------------------------------------------
    local count_2of3 count_3of3
    count_2of3=$(wc -l < "${outdir}/${group}_consensus_2of3.bed")
    count_3of3=$(wc -l < "${outdir}/${group}_consensus_3of3.bed")

    echo "  [${group}] Consensus peaks — 2/3 pairs: ${count_2of3}  |  3/3 pairs: ${count_3of3}"

    # Append to summary table
    echo -e "${group}\t${threshold}\t${n1}\t${n2}\t${n3}\t${n_1v2}\t${n_1v3}\t${n_2v3}\t${count_2of3}\t${count_3of3}" \
        >> "${SUMMARY}"
}

# =============================================================================
# Run pipeline
# =============================================================================
for threshold in $THRESHOLDS; do
    echo ""
    echo "========================================"
    echo "  IDR Threshold: ${threshold}"
    echo "========================================"
    for group in $GROUPS; do
        run_three_rep_idr "$group" "$threshold"
    done
done

# =============================================================================
# Final summary
# =============================================================================
echo ""
echo "========================================"
echo "  IDR Analysis Complete"
echo "  Finished : $(date)"
echo "  Summary  : ${SUMMARY}"
echo "========================================"
echo ""
echo "--- Summary Table ---"
column -t "${SUMMARY}"
