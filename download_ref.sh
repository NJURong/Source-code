#!/bin/bash
# Download mm10 reference genome and auxiliary files for CUT&Tag analysis

set -euo pipefail

REF_DIR="/mnt/i/LQ/WR_Chip_seq/ref"
mkdir -p $REF_DIR
cd $REF_DIR

echo "=== Downloading mm10 reference genome ==="

# Option 1: Download pre-built bowtie2 index from AWS
echo "Trying AWS..."
wget -c https://genome-idx.s3.amazonaws.com/bt/UCSC_mm10.tar.gz && \
    tar -xzf UCSC_mm10.tar.gz && \
    mv mm10 mm10_bt2_index && \
    echo "AWS download successful" || echo "AWS failed"

# Option 2: If AWS fails, download from JHU FTP
if [ ! -d "mm10_bt2_index" ]; then
    echo "Trying JHU FTP..."
    wget -c ftp://ftp.ccb.jhu.edu/pub/data/bowtie2_indexes/mm10.zip && \
        unzip mm10.zip && \
        mv mm10 mm10_bt2_index && \
        echo "JHU download successful" || echo "JHU failed"
fi

# Option 3: Download from Illumina iGenomes (most reliable)
if [ ! -d "mm10_bt2_index" ]; then
    echo "Trying Illumina iGenomes..."
    wget -c ftp://igenome:G3nom3s4u@ussd-ftp.illumina.com/Mus_musculus/UCSC/mm10/Mus_musculus_UCSC_mm10.tar.gz && \
        tar -xzf Mus_musculus_UCSC_mm10.tar.gz && \
        mkdir -p mm10_bt2_index && \
        cp Mus_musculus/UCSC/mm10/Sequence/Bowtie2Index/genome* mm10_bt2_index/ && \
        echo "iGenomes download successful" || echo "iGenomes failed"
fi

# Download chromosome sizes
echo "=== Downloading chromosome sizes ==="
wget -q http://hgdownload.soe.ucsc.edu/goldenPath/mm10/bigZips/mm10.chrom.sizes -O mm10.chrom.sizes || \
    echo "WARNING: Could not download chrom.sizes"

# Download blacklist (ENCODE)
echo "=== Downloading blacklist regions ==="
wget -q https://raw.githubusercontent.com/Boyle-Lab/Blacklist/master/lists/mm10-blacklist.v2.bed.gz -O mm10-blacklist.v2.bed.gz && \
    gunzip mm10-blacklist.v2.bed.gz || \
    echo "WARNING: Could not download blacklist"

# Generate TSS bed from refGene (optional, can also use GTF)
echo "=== Note ==="
echo "TSS bed file (mm10_TSS.bed) needs to be generated from UCSC refGene or GTF."
echo "You can download from UCSC Table Browser or use a GTF file with:"
echo "  awk 'BEGIN{OFS=\"\\t\"} \$3==\"transcript\" {print \$1,\$4,\$5,\$10,\$6,\$7}' genes.gtf > mm10_TSS.bed"

echo "=== Reference setup complete ==="
ls -la $REF_DIR
