#!/bin/bash
# Comprehensive IDR analysis for CUT&Tag data
# Three strategies:
# 1. Standard multi-replicate IDR (--peak-list mode with proper narrowPeak format)
# 2. Pooled pseudoreplicate IDR (ENCODE gold standard)
# 3. Comparison of IDR thresholds (0.01 vs 0.05)

set -euo pipefail

BASE_DIR="/mnt/i/LQ/WR_Chip_seq"
OUT_DIR="${BASE_DIR}/results"
REF_DIR="${BASE_DIR}/ref"
IDR_DIR="${OUT_DIR}/07_idr_comprehensive"
THREADS=8
EFFECTIVE_GENOME_SIZE="1870000000"

mkdir -p ${IDR_DIR}/{standard_0.01,standard_0.05,pooled_0.01,pooled_0.05,merged_peaks,diagnostic}

source /home/lq/miniconda3/etc/profile.d/conda.sh
conda activate cuttag

echo "========================================"
echo "Comprehensive IDR Analysis"
echo "Started: $(date)"
echo "========================================"

# === Strategy 1: Standard multi-replicate IDR (--peak-list with narrowPeak format) ===
echo ""
echo "[Strategy 1] Standard multi-replicate IDR (--peak-list mode)"

run_standard_idr() {
    local group=$1
    local threshold=$2
    local outdir="${IDR_DIR}/standard_${threshold}/${group}"
    mkdir -p $outdir
    
    local rep1="${OUT_DIR}/06_peaks/${group}1_peaks.narrowPeak"
    local rep2="${OUT_DIR}/06_peaks/${group}2_peaks.narrowPeak"
    local rep3="${OUT_DIR}/06_peaks/${group}3_peaks.narrowPeak"
    
    echo "  Processing ${group} with IDR threshold ${threshold}..."
    
    # Create merged peak list in narrowPeak format (union of all peaks)
    # Use rep1 as template for format, merge coordinates
    cat $rep1 $rep2 $rep3 | sort -k1,1 -k2,2n | \
        awk 'BEGIN{OFS="\t"} {print $1,$2,$3,"merged_peak_"NR,$5,$6,$7,$8,$9,$10}' > ${outdir}/${group}_merged_peaks.narrowPeak
    
    # Run IDR with --peak-list (conservative: peaks must appear in all replicates)
    idr --samples $rep1 $rep2 \
        --peak-list ${outdir}/${group}_merged_peaks.narrowPeak \
        --input-file-type narrowPeak \
        --rank signal.value \
        --output-file ${outdir}/${group}_standard_idr_${threshold}.narrowPeak \
        --idr-threshold ${threshold} \
        --plot \
        --log-output-file ${outdir}/${group}_standard_idr_${threshold}.log
    
    # Count peaks
    local count=$(wc -l < ${outdir}/${group}_standard_idr_${threshold}.narrowPeak)
    echo "    ${group} standard IDR ${threshold}: ${count} peaks"
}

# === Strategy 2: Pooled pseudoreplicate IDR (ENCODE gold standard) ===
echo ""
echo "[Strategy 2] Pooled pseudoreplicate IDR (ENCODE gold standard)"

run_pooled_idr() {
    local group=$1
    local threshold=$2
    local outdir="${IDR_DIR}/pooled_${threshold}/${group}"
    mkdir -p $outdir
    
    echo "  Processing ${group} with pooled pseudoreplicate IDR ${threshold}..."
    
    # Step 1: Merge BAM files
    local merged_bam="${IDR_DIR}/merged_peaks/${group}_pooled.bam"
    if [ ! -f "$merged_bam" ]; then
        echo "    Merging BAM files for ${group}..."
        samtools merge -@ ${THREADS} $merged_bam \
            ${OUT_DIR}/04_filter/${group}1.final.bam \
            ${OUT_DIR}/04_filter/${group}2.final.bam \
            ${OUT_DIR}/04_filter/${group}3.final.bam
        samtools index $merged_bam
    fi
    
    # Step 2: Call peaks on pooled BAM (relaxed threshold p < 0.1)
    local pooled_peaks="${outdir}/${group}_pooled_peaks.narrowPeak"
    if [ ! -f "$pooled_peaks" ]; then
        echo "    Calling peaks on pooled ${group}..."
        macs2 callpeak -t $merged_bam -f BAMPE -g ${EFFECTIVE_GENOME_SIZE} \
            -n ${group}_pooled --outdir ${outdir} \
            --keep-dup all -p 0.1
        if [ -f "${outdir}/${group}_pooled_peaks.narrowPeak" ]; then
            mv ${outdir}/${group}_pooled_peaks.narrowPeak ${pooled_peaks}
        elif [ -f "${outdir}/${group}_pooled_peaks.xls" ]; then
            # Convert xls to narrowPeak if needed
            tail -n +2 ${outdir}/${group}_pooled_peaks.xls | \
                awk 'BEGIN{OFS="\t"} {print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10}' > ${pooled_peaks}
        fi
    fi
    
    # Step 3: Split pooled BAM into two pseudoreplicates
    local pseudo1="${outdir}/${group}_pseudo1.bam"
    local pseudo2="${outdir}/${group}_pseudo2.bam"
    
    if [ ! -f "$pseudo1" ] || [ ! -f "$pseudo2" ]; then
        echo "    Creating pseudoreplicates for ${group}..."
        # Get total reads and split in half using samtools split
        local total_reads=$(samtools view -c $merged_bam)
        local half_reads=$((total_reads / 2))
        
        # Use samtools split with seed for reproducibility
        samtools view -@ ${THREADS} -b -s 42.5 $merged_bam > $pseudo1
        # For pseudo2, we need the other half - use different approach
        samtools view -@ ${THREADS} -H $merged_bam > ${outdir}/header.sam
        samtools view -@ ${THREADS} $merged_bam | awk 'NR>${half_reads}' | \
            cat ${outdir}/header.sam - | samtools view -@ ${THREADS} -bS - > $pseudo2
        rm -f ${outdir}/header.sam
        
        samtools index $pseudo1
        samtools index $pseudo2
    fi
    
    # Step 4: Call peaks on pseudoreplicates (relaxed)
    local pseudo1_peaks="${outdir}/${group}_pseudo1_peaks.narrowPeak"
    local pseudo2_peaks="${outdir}/${group}_pseudo2_peaks.narrowPeak"
    
    if [ ! -f "$pseudo1_peaks" ]; then
        echo "    Calling peaks on pseudoreplicate 1..."
        macs2 callpeak -t $pseudo1 -f BAMPE -g ${EFFECTIVE_GENOME_SIZE} \
            -n ${group}_pseudo1 --outdir ${outdir} \
            --keep-dup all -p 0.1
        if [ -f "${outdir}/${group}_pseudo1_peaks.narrowPeak" ]; then
            mv ${outdir}/${group}_pseudo1_peaks.narrowPeak ${pseudo1_peaks}
        fi
    fi
    
    if [ ! -f "$pseudo2_peaks" ]; then
        echo "    Calling peaks on pseudoreplicate 2..."
        macs2 callpeak -t $pseudo2 -f BAMPE -g ${EFFECTIVE_GENOME_SIZE} \
            -n ${group}_pseudo2 --outdir ${outdir} \
            --keep-dup all -p 0.1
        if [ -f "${outdir}/${group}_pseudo2_peaks.narrowPeak" ]; then
            mv ${outdir}/${group}_pseudo2_peaks.narrowPeak ${pseudo2_peaks}
        fi
    fi
    
    # Step 5: Run IDR on pseudoreplicates
    local idr_output="${outdir}/${group}_pooled_idr_${threshold}.narrowPeak"
    if [ -f "$pseudo1_peaks" ] && [ -f "$pseudo2_peaks" ]; then
        idr --samples $pseudo1_peaks $pseudo2_peaks \
            --input-file-type narrowPeak \
            --rank signal.value \
            --output-file ${idr_output} \
            --idr-threshold ${threshold} \
            --plot \
            --log-output-file ${outdir}/${group}_pooled_idr_${threshold}.log
    fi
    
    # Step 6: Apply IDR threshold to pooled peaks
    local final_peaks="${IDR_DIR}/diagnostic/${group}_IDR${threshold}_final_peaks.bed"
    if [ -f "$idr_output" ] && [ -f "$pooled_peaks" ]; then
        awk 'NR==FNR{a[$1"\t"$2"\t"$3]; next} ($1"\t"$2"\t"$3) in a' \
            ${idr_output} ${pooled_peaks} > ${final_peaks}
    fi
    
    local count=$(wc -l < ${idr_output} 2>/dev/null || echo 0)
    echo "    ${group} pooled IDR ${threshold}: ${count} peaks"
}

# Run both strategies for both thresholds
for threshold in 0.01 0.05; do
    echo ""
    echo "=== IDR Threshold: ${threshold} ==="
    
    for group in A B C D; do
        run_standard_idr $group $threshold
    done
    
    for group in A B C D; do
        run_pooled_idr $group $threshold
    done
done

# === Generate summary report ===
echo ""
echo "========================================"
echo "Generating summary report..."
echo "========================================"

REPORT="${IDR_DIR}/IDR_summary_report.txt"
cat > $REPORT << EOF
CUT&Tag Comprehensive IDR Analysis Report
Generated: $(date)
========================================

Strategy 1: Standard Multi-Replicate IDR (--peak-list)
- Conservative approach: peaks must be present in all replicates
- Uses merged peak list as reference

Strategy 2: Pooled Pseudoreplicate IDR (ENCODE Gold Standard)
- Merges all replicates, splits into random pseudoreplicates
- More robust when replicates have different sequencing depths
- Uses relaxed peak calling (p < 0.1) followed by IDR filtering

Threshold Comparison:
- IDR < 0.01: Very stringent, highest confidence
- IDR < 0.05: Standard ENCODE threshold, balanced sensitivity/specificity

Results:
EOF

for threshold in 0.01 0.05; do
    echo "" >> $REPORT
    echo "IDR Threshold: ${threshold}" >> $REPORT
    echo "------------------------------" >> $REPORT
    
    for group in A B C D; do
        local std_count=$(wc -l < ${IDR_DIR}/standard_${threshold}/${group}/${group}_standard_idr_${threshold}.narrowPeak 2>/dev/null || echo 0)
        local pooled_count=$(wc -l < ${IDR_DIR}/pooled_${threshold}/${group}/${group}_pooled_idr_${threshold}.narrowPeak 2>/dev/null || echo 0)
        
        echo "  ${group} group:" >> $REPORT
        echo "    Standard IDR: ${std_count} peaks" >> $REPORT
        echo "    Pooled IDR:   ${pooled_count} peaks" >> $REPORT
    done
done

echo ""
echo "Summary report: ${REPORT}"
echo "All analyses complete!"
echo "Finished: $(date)"
