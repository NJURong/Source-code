#!/bin/bash
# Run alignment and downstream analysis for CUT&Tag

set -euo pipefail

BASE_DIR="/mnt/i/LQ/WR_Chip_seq"
REF_GENOME="${BASE_DIR}/ref/mm10/mm10"
OUT_DIR="${BASE_DIR}/results"
THREADS=8

echo "[STEP 3] Aligning to mm10 with Bowtie2..."
for group in A B C D; do
    for rep in 1 2 3; do
        SAMPLE="${group}${rep}"
        TRIM_R1="${OUT_DIR}/02_trim/${SAMPLE}_val_1.fq.gz"
        TRIM_R2="${OUT_DIR}/02_trim/${SAMPLE}_val_2.fq.gz"
        
        if [ ! -f "$TRIM_R1" ] || [ ! -f "$TRIM_R2" ]; then
            echo "WARNING: Trimmed files not found for ${SAMPLE}, skipping..."
            continue
        fi
        
        if [ -f "${OUT_DIR}/03_align/${SAMPLE}.sorted.bam" ]; then
            echo "  ${SAMPLE} already aligned, skipping..."
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
        echo "  ${SAMPLE} alignment done"
    done
done
echo "Alignment complete!"
