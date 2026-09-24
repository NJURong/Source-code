#!/bin/bash
# Prepare annotation files for CUT&Tag analysis

set -euo pipefail

REF_DIR="/mnt/i/LQ/WR_Chip_seq/ref"
mkdir -p $REF_DIR
cd $REF_DIR

echo "=== Preparing annotation files ==="

# 1. Chromosome sizes (mm10)
if [ ! -f "mm10.chrom.sizes" ]; then
    echo "Downloading mm10 chromosome sizes..."
    wget -q http://hgdownload.soe.ucsc.edu/goldenPath/mm10/bigZips/mm10.chrom.sizes -O mm10.chrom.sizes || {
        echo "WARNING: Failed to download chrom.sizes. Will generate from FASTA if available."
    }
fi

# 2. Blacklist regions (ENCODE mm10)
if [ ! -f "mm10-blacklist.v2.bed" ]; then
    echo "Downloading mm10 blacklist..."
    wget -q https://raw.githubusercontent.com/Boyle-Lab/Blacklist/master/lists/mm10-blacklist.v2.bed.gz -O mm10-blacklist.v2.bed.gz 2>/dev/null && \
        gunzip mm10-blacklist.v2.bed.gz || {
        echo "WARNING: Failed to download blacklist. Analysis will proceed without blacklist filtering."
    }
fi

# 3. TSS bed file from UCSC refGene
if [ ! -f "mm10_TSS.bed" ]; then
    echo "Generating mm10 TSS bed file..."
    
    # Try to download refGene table from UCSC
    wget -q "http://hgdownload.soe.ucsc.edu/goldenPath/mm10/database/refGene.txt.gz" -O refGene.txt.gz 2>/dev/null && \
        gunzip -f refGene.txt.gz && \
        awk 'BEGIN{OFS="\t"} {
            if ($4 == "+") {
                start = $5
            } else {
                start = $6 - 1
            }
            print $3, start, start+1, $2, $4, $13
        }' refGene.txt > mm10_TSS.bed && \
        echo "TSS bed generated from refGene" || {
        echo "WARNING: Could not generate TSS bed from refGene."
        echo "Please manually provide mm10_TSS.bed or a GTF file."
    }
fi

# 4. Target gene promoter regions (Hk2, Acss2, Acsl4)
# These coordinates are approximate for mm10 - should be verified with current annotation
cat > target_genes_mm10.bed << 'EOF'
# Hk2 (hexokinase 2) - chr6:130,880,000-130,890,000 (approximate)
chr6	130880000	130890000	Hk2
# Acss2 (acetyl-CoA synthetase 2) - chr19:17,445,000-17,455,000 (approximate)
chr19	17445000	17455000	Acss2
# Acsl4 (acyl-CoA synthetase long-chain family member 4) - chr1:75,300,000-75,310,000 (approximate)
chr1	75300000	75310000	Acsl4
EOF

echo "=== Annotation files prepared ==="
ls -la $REF_DIR/*.bed $REF_DIR/*.sizes 2>/dev/null || true
