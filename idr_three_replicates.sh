#!/bin/bash
# Proper 3-replicate IDR analysis
# IDR only supports pairwise comparison, so we run all 3 pairs and take intersection

set -euo pipefail

BASE_DIR="/mnt/i/LQ/WR_Chip_seq"
OUT_DIR="${BASE_DIR}/results"
IDR_DIR="${OUT_DIR}/07_idr_three_replicates"

mkdir -p ${IDR_DIR}/{0.01,0.05}

source /home/lq/miniconda3/etc/profile.d/conda.sh
conda activate cuttag

echo "========================================"
echo "3-Replicate IDR Analysis"
echo "Started: $(date)"
echo "========================================"

run_three_rep_idr() {
    local group=$1
    local threshold=$2
    local outdir="${IDR_DIR}/${threshold}/${group}"
    mkdir -p $outdir
    
    local rep1="${OUT_DIR}/06_peaks/${group}1_peaks.narrowPeak"
    local rep2="${OUT_DIR}/06_peaks/${group}2_peaks.narrowPeak"
    local rep3="${OUT_DIR}/06_peaks/${group}3_peaks.narrowPeak"
    
    echo "Processing ${group} with IDR threshold ${threshold}..."
    
    # Run all 3 pairwise IDR comparisons
    echo "  IDR: ${group}1 vs ${group}2"
    idr --samples $rep1 $rep2 \
        --input-file-type narrowPeak \
        --rank signal.value \
        --output-file ${outdir}/${group}_1vs2_idr.narrowPeak \
        --idr-threshold ${threshold} \
        --plot \
        --log-output-file ${outdir}/${group}_1vs2_idr.log
    
    echo "  IDR: ${group}1 vs ${group}3"
    idr --samples $rep1 $rep3 \
        --input-file-type narrowPeak \
        --rank signal.value \
        --output-file ${outdir}/${group}_1vs3_idr.narrowPeak \
        --idr-threshold ${threshold} \
        --plot \
        --log-output-file ${outdir}/${group}_1vs3_idr.log
    
    echo "  IDR: ${group}2 vs ${group}3"
    idr --samples $rep2 $rep3 \
        --input-file-type narrowPeak \
        --rank signal.value \
        --output-file ${outdir}/${group}_2vs3_idr.narrowPeak \
        --idr-threshold ${threshold} \
        --plot \
        --log-output-file ${outdir}/${group}_2vs3_idr.log
    
    # Create consensus: peaks present in at least 2 of 3 pairwise comparisons
    # Extract peak coordinates (chr, start, end)
    cut -f1-3 ${outdir}/${group}_1vs2_idr.narrowPeak | sort -u > ${outdir}/tmp1.bed
    cut -f1-3 ${outdir}/${group}_1vs3_idr.narrowPeak | sort -u > ${outdir}/tmp2.bed
    cut -f1-3 ${outdir}/${group}_2vs3_idr.narrowPeak | sort -u > ${outdir}/tmp3.bed
    
    # Peaks appearing in at least 2 comparisons (majority vote)
    cat ${outdir}/tmp1.bed ${outdir}/tmp2.bed ${outdir}/tmp3.bed | \
        sort | uniq -c | awk '$1 >= 2 {print $2"\t"$3"\t"$4}' | \
        sort -k1,1 -k2,2n > ${outdir}/${group}_consensus_2of3.bed
    
    # Strict consensus: peaks in all 3 comparisons
    cat ${outdir}/tmp1.bed ${outdir}/tmp2.bed ${outdir}/tmp3.bed | \
        sort | uniq -c | awk '$1 >= 3 {print $2"\t"$3"\t"$4}' | \
        sort -k1,1 -k2,2n > ${outdir}/${group}_consensus_3of3.bed
    
    rm -f ${outdir}/tmp*.bed
    
    local count_2of3=$(wc -l < ${outdir}/${group}_consensus_2of3.bed)
    local count_3of3=$(wc -l < ${outdir}/${group}_consensus_3of3.bed)
    
    echo "  ${group} consensus peaks (2/3 pairs): ${count_2of3}"
    echo "  ${group} consensus peaks (3/3 pairs): ${count_3of3}"
}

# Run for both thresholds
for threshold in 0.01 0.05; do
    echo ""
    echo "=== IDR Threshold: ${threshold} ==="
    for group in A B C D; do
        run_three_rep_idr $group $threshold
    done
done

echo ""
echo "========================================"
echo "3-Replicate IDR Analysis Complete"
echo "Finished: $(date)"
echo "========================================"
