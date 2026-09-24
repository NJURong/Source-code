#!/bin/bash
# CUT&Tag Analysis Pipeline for WR_Chip_seq
# Conditions: A=Control, B=TNF-α, C=TNF-α+lactate, D=TNF-α+lactate+oxamate
# 3 biological replicates per condition (A1-A3, B1-B3, C1-C3, D1-D3)
# Target: H3K18la

set -euo pipefail

# === CONFIGURATION ===
BASE_DIR="/mnt/i/LQ/WR_Chip_seq"
DATA_DIR="${BASE_DIR}/CUTTAG"
REF_DIR="${BASE_DIR}/ref"
OUT_DIR="${BASE_DIR}/results"
REF_GENOME="${REF_DIR}/mm10/mm10"  # bowtie2 index prefix
CHROM_SIZE="${REF_DIR}/mm10.chrom.sizes"
BLACKLIST="${REF_DIR}/mm10-blacklist.v2.bed"
TSS_BED="${REF_DIR}/mm10_TSS.bed"
TARGET_GENES="${REF_DIR}/target_genes_mm10.bed"

# Effective genome size for mm10 (macs2)
EFFECTIVE_GENOME_SIZE="1870000000"

# Threads
THREADS=8

# Create output directories
mkdir -p ${OUT_DIR}/{01_fastqc,02_trim,03_align,04_filter,05_bamqc,06_peaks,07_idr,08_correlation,09_bigwig,10_visualization}

# === LOG ===
LOG="${OUT_DIR}/pipeline.log"
exec > >(tee -a "$LOG") 2>&1

echo "========================================"
echo "CUT&Tag Pipeline Started: $(date)"
echo "========================================"

# === STEP 0: Check dependencies ===
echo "[STEP 0] Checking dependencies..."
for cmd in bowtie2 samtools macs2 idr bedtools bamCoverage computeMatrix plotProfile plotHeatmap; do
    if ! command -v $cmd &> /dev/null; then
        echo "ERROR: $cmd not found. Please activate cuttag environment."
        exit 1
    fi
done
echo "All dependencies found."

# === Check reference genome ===
if [ ! -f "${REF_GENOME}.1.bt2" ]; then
    echo "ERROR: Bowtie2 index not found at ${REF_GENOME}"
    echo "Please build index first: bowtie2-build mm10.fa mm10"
    exit 1
fi

# === STEP 1: FastQC (optional, for QC report) ===
# echo "[STEP 1] Running FastQC..."
# for sample_dir in ${DATA_DIR}/*/*; do
#     if [ -d "$sample_dir" ]; then
#         fastq_files=(${sample_dir}/*.fastq.gz)
#         if [ -f "${fastq_files[0]}" ]; then
#             fastqc -t ${THREADS} -o ${OUT_DIR}/01_fastqc ${fastq_files[@]}
#         fi
#     fi
# done

# === STEP 2: Trim Galore (adapter trimming) ===
echo "[STEP 2] Trimming adapters with Trim Galore..."
for group in A B C D; do
    for rep in 1 2 3; do
        SAMPLE="${group}${rep}"
        SAMPLE_DIR="${DATA_DIR}/${group}/${SAMPLE}"
        if [ ! -d "$SAMPLE_DIR" ]; then
            echo "WARNING: $SAMPLE_DIR not found, skipping..."
            continue
        fi
        
        R1=$(ls ${SAMPLE_DIR}/*_R1_*.fastq.gz 2>/dev/null | head -1)
        R2=$(ls ${SAMPLE_DIR}/*_R2_*.fastq.gz 2>/dev/null | head -1)
        
        if [ -z "$R1" ] || [ -z "$R2" ]; then
            echo "WARNING: Fastq files not found for ${SAMPLE}, skipping..."
            continue
        fi
        
        echo "  Processing ${SAMPLE}..."
        trim_galore --paired --cores ${THREADS} \
            --output_dir ${OUT_DIR}/02_trim \
            --basename ${SAMPLE} \
            $R1 $R2
    done
done

# === STEP 3: Alignment with Bowtie2 ===
echo "[STEP 3] Aligning to mm10 with Bowtie2..."
for group in A B C D; do
    for rep in 1 2 3; do
        SAMPLE="${group}${rep}"
        TRIM_R1="${OUT_DIR}/02_trim/${SAMPLE}_val_1.fq.gz"
        TRIM_R2="${OUT_DIR}/02_trim/${SAMPLE}_val_2.fq.gz"
        
        if [ ! -f "$TRIM_R1" ] || [ ! -f "$TRIM_R2" ]; then
            echo "WARNING: Trimmed files not found for ${SAMPLE}, skipping alignment..."
            continue
        fi
        
        echo "  Aligning ${SAMPLE}..."
        bowtie2 -p ${THREADS} -x ${REF_GENOME} \
            -1 $TRIM_R1 -2 $TRIM_R2 \
            --local --very-sensitive-local --no-unal --no-mixed \
            --no-discordant --phred33 -I 10 -X 700 \
            2> ${OUT_DIR}/03_align/${SAMPLE}.bowtie2.log | \
            samtools view -@ ${THREADS} -Sb - > ${OUT_DIR}/03_align/${SAMPLE}.bam
        
        samtools sort -@ ${THREADS} -o ${OUT_DIR}/03_align/${SAMPLE}.sorted.bam ${OUT_DIR}/03_align/${SAMPLE}.bam
        rm ${OUT_DIR}/03_align/${SAMPLE}.bam
        samtools index ${OUT_DIR}/03_align/${SAMPLE}.sorted.bam
    done
done

# === STEP 4: Filter BAM (MAPQ >= 30, remove duplicates, blacklist) ===
echo "[STEP 4] Filtering BAM files..."
for group in A B C D; do
    for rep in 1 2 3; do
        SAMPLE="${group}${rep}"
        BAM="${OUT_DIR}/03_align/${SAMPLE}.sorted.bam"
        
        if [ ! -f "$BAM" ]; then
            echo "WARNING: BAM not found for ${SAMPLE}, skipping..."
            continue
        fi
        
        echo "  Filtering ${SAMPLE}..."
        
        # Remove low quality, sort by coordinate
        samtools view -@ ${THREADS} -b -q 30 $BAM > ${OUT_DIR}/04_filter/${SAMPLE}.q30.bam
        
        # Mark duplicates with Picard or samtools (using coordinate-based approach)
        # For CUT&Tag, we typically remove exact duplicates
        samtools collate -@ ${THREADS} -o - ${OUT_DIR}/04_filter/${SAMPLE}.q30.bam | \
            samtools fixmate -@ ${THREADS} -m - - | \
            samtools sort -@ ${THREADS} - | \
            samtools markdup -@ ${THREADS} -r - ${OUT_DIR}/04_filter/${SAMPLE}.rmdup.bam
        
        # Remove blacklist regions if available
        if [ -f "$BLACKLIST" ]; then
            bedtools intersect -v -a ${OUT_DIR}/04_filter/${SAMPLE}.rmdup.bam -b $BLACKLIST > ${OUT_DIR}/04_filter/${SAMPLE}.filtered.bam
        else
            mv ${OUT_DIR}/04_filter/${SAMPLE}.rmdup.bam ${OUT_DIR}/04_filter/${SAMPLE}.filtered.bam
        fi
        
        samtools sort -@ ${THREADS} -o ${OUT_DIR}/04_filter/${SAMPLE}.final.bam ${OUT_DIR}/04_filter/${SAMPLE}.filtered.bam
        samtools index ${OUT_DIR}/04_filter/${SAMPLE}.final.bam
        
        # Clean up intermediates
        rm -f ${OUT_DIR}/04_filter/${SAMPLE}.q30.bam ${OUT_DIR}/04_filter/${SAMPLE}.rmdup.bam ${OUT_DIR}/04_filter/${SAMPLE}.filtered.bam
        
        # Stats
        echo "  ${SAMPLE} final reads: $(samtools view -c ${OUT_DIR}/04_filter/${SAMPLE}.final.bam)"
    done
done

# === STEP 5: Peak Calling with MACS2 ===
echo "[STEP 5] Calling peaks with MACS2..."
for group in A B C D; do
    for rep in 1 2 3; do
        SAMPLE="${group}${rep}"
        BAM="${OUT_DIR}/04_filter/${SAMPLE}.final.bam"
        
        if [ ! -f "$BAM" ]; then
            echo "WARNING: Final BAM not found for ${SAMPLE}, skipping peak calling..."
            continue
        fi
        
        echo "  Calling peaks for ${SAMPLE}..."
        macs2 callpeak -t $BAM \
            -f BAMPE \
            -g ${EFFECTIVE_GENOME_SIZE} \
            -n ${SAMPLE} \
            --outdir ${OUT_DIR}/06_peaks \
            --keep-dup all \
            --cutoff-analysis \
            2> ${OUT_DIR}/06_peaks/${SAMPLE}.macs2.log
    done
done

# === STEP 6: IDR Analysis (Irreproducible Discovery Rate) ===
echo "[STEP 6] Running IDR analysis for replicates..."

# Function to run IDR for a group
run_idr() {
    local group=$1
    local rep1="${OUT_DIR}/06_peaks/${group}1_peaks.narrowPeak"
    local rep2="${OUT_DIR}/06_peaks/${group}2_peaks.narrowPeak"
    local rep3="${OUT_DIR}/06_peaks/${group}3_peaks.narrowPeak"
    local output_prefix="${OUT_DIR}/07_idr/${group}"
    
    mkdir -p ${output_prefix}
    
    # IDR pairwise comparisons
    if [ -f "$rep1" ] && [ -f "$rep2" ]; then
        echo "  IDR: ${group}1 vs ${group}2"
        idr --samples $rep1 $rep2 \
            --input-file-type narrowPeak \
            --output-file ${output_prefix}/${group}_1vs2_idr.narrowPeak \
            --plot \
            --log-output-file ${output_prefix}/${group}_1vs2_idr.log
    fi
    
    if [ -f "$rep1" ] && [ -f "$rep3" ]; then
        echo "  IDR: ${group}1 vs ${group}3"
        idr --samples $rep1 $rep3 \
            --input-file-type narrowPeak \
            --output-file ${output_prefix}/${group}_1vs3_idr.narrowPeak \
            --plot \
            --log-output-file ${output_prefix}/${group}_1vs3_idr.log
    fi
    
    if [ -f "$rep2" ] && [ -f "$rep3" ]; then
        echo "  IDR: ${group}2 vs ${group}3"
        idr --samples $rep2 $rep3 \
            --input-file-type narrowPeak \
            --output-file ${output_prefix}/${group}_2vs3_idr.narrowPeak \
            --plot \
            --log-output-file ${output_prefix}/${group}_2vs3_idr.log
    fi
    
    # Combine IDR peaks (intersection of all pairwise IDR peaks)
    if [ -f "${output_prefix}/${group}_1vs2_idr.narrowPeak" ] && \
       [ -f "${output_prefix}/${group}_1vs3_idr.narrowPeak" ] && \
       [ -f "${output_prefix}/${group}_2vs3_idr.narrowPeak" ]; then
        
        # Extract peak regions (first 3 columns) and find common peaks
        cut -f1-3 ${output_prefix}/${group}_1vs2_idr.narrowPeak | sort -u > ${output_prefix}/tmp1.bed
        cut -f1-3 ${output_prefix}/${group}_1vs3_idr.narrowPeak | sort -u > ${output_prefix}/tmp2.bed
        cut -f1-3 ${output_prefix}/${group}_2vs3_idr.narrowPeak | sort -u > ${output_prefix}/tmp3.bed
        
        # Find peaks present in at least 2 of 3 comparisons
        cat ${output_prefix}/tmp1.bed ${output_prefix}/tmp2.bed ${output_prefix}/tmp3.bed | \
            sort | uniq -c | awk '$1 >= 2 {print $2"\t"$3"\t"$4}' > ${output_prefix}/${group}_consensus_idr_peaks.bed
        
        rm -f ${output_prefix}/tmp*.bed
        
        echo "  ${group} consensus IDR peaks: $(wc -l < ${output_prefix}/${group}_consensus_idr_peaks.bed)"
    fi
}

for group in A B C D; do
    run_idr $group
done

# === STEP 7: Generate BigWig files for visualization ===
echo "[STEP 7] Generating BigWig files..."
for group in A B C D; do
    for rep in 1 2 3; do
        SAMPLE="${group}${rep}"
        BAM="${OUT_DIR}/04_filter/${SAMPLE}.final.bam"
        
        if [ ! -f "$BAM" ]; then
            continue
        fi
        
        echo "  Creating BigWig for ${SAMPLE}..."
        bamCoverage -b $BAM \
            -o ${OUT_DIR}/09_bigwig/${SAMPLE}.bw \
            --binSize 10 \
            --normalizeUsing RPGC \
            --effectiveGenomeSize ${EFFECTIVE_GENOME_SIZE} \
            --ignoreDuplicates \
            --numberOfProcessors ${THREADS}
    done
done

# === STEP 8: Replicate Correlation Analysis ===
echo "[STEP 8] Computing replicate correlations..."

# Merge replicates per group for comparison
for group in A B C D; do
    BWS=()
    for rep in 1 2 3; do
        BW="${OUT_DIR}/09_bigwig/${group}${rep}.bw"
        [ -f "$BW" ] && BWS+=($BW)
    done
    
    if [ ${#BWS[@]} -ge 2 ]; then
        echo "  Correlation for group ${group}..."
        
        # MultiBigwigSummary
        multiBigwigSummary bins -b ${BWS[@]} \
            -o ${OUT_DIR}/08_correlation/${group}_replicate.npz \
            --binSize 10000 \
            --numberOfProcessors ${THREADS}
        
        # Plot correlation heatmap
        plotCorrelation -in ${OUT_DIR}/08_correlation/${group}_replicate.npz \
            --corMethod spearman \
            --whatToPlot heatmap \
            -o ${OUT_DIR}/08_correlation/${group}_replicate_correlation.pdf \
            --outFileCorMatrix ${OUT_DIR}/08_correlation/${group}_replicate_correlation.txt \
            --plotNumbers
    fi
done

# === STEP 9: TSS Enrichment and Profile Plots ===
echo "[STEP 9] Generating TSS enrichment profiles..."

if [ ! -f "$TSS_BED" ]; then
    echo "WARNING: TSS bed file not found at ${TSS_BED}. Skipping TSS analysis."
else
    for group in A B C D; do
        BWS=()
        for rep in 1 2 3; do
            BW="${OUT_DIR}/09_bigwig/${group}${rep}.bw"
            [ -f "$BW" ] && BWS+=($BW)
        done
        
        if [ ${#BWS[@]} -ge 1 ]; then
            computeMatrix reference-point \
                -S ${BWS[@]} \
                -R $TSS_BED \
                --referencePoint TSS \
                -a 3000 -b 3000 \
                -o ${OUT_DIR}/10_visualization/${group}_TSS_matrix.gz \
                --numberOfProcessors ${THREADS}
            
            plotProfile -m ${OUT_DIR}/10_visualization/${group}_TSS_matrix.gz \
                -o ${OUT_DIR}/10_visualization/${group}_TSS_profile.pdf \
                --perGroup \
                --plotTitle "${group} H3K18la TSS Enrichment"
        fi
    done
fi

# === STEP 10: Gene-specific enrichment (Hk2, Acss2, Acsl4) ===
echo "[STEP 10] Checking Hk2, Acss2, Acsl4 promoter enrichment..."

if [ ! -f "$TARGET_GENES" ]; then
    echo "WARNING: Target genes bed file not found at ${TARGET_GENES}. Skipping target gene analysis."
else
    for group in A B C D; do
        BWS=()
        for rep in 1 2 3; do
            BW="${OUT_DIR}/09_bigwig/${group}${rep}.bw"
            [ -f "$BW" ] && BWS+=($BW)
        done
        
        if [ ${#BWS[@]} -ge 1 ]; then
            computeMatrix scale-regions \
                -S ${BWS[@]} \
                -R $TARGET_GENES \
                -a 3000 -b 3000 \
                -m 5000 \
                -o ${OUT_DIR}/10_visualization/${group}_target_genes_matrix.gz \
                --numberOfProcessors ${THREADS}
            
            plotHeatmap -m ${OUT_DIR}/10_visualization/${group}_target_genes_matrix.gz \
                -o ${OUT_DIR}/10_visualization/${group}_target_genes_heatmap.pdf \
                --colorMap RdBu \
                --whatToShow "heatmap and colorbar"
        fi
    done
fi

echo "========================================"
echo "CUT&Tag Pipeline Completed: $(date)"
echo "========================================"
echo "Results in: ${OUT_DIR}"
