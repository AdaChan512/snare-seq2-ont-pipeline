#!/bin/bash
# Step 5: minimap2 比对 (师兄参数)
# 用法: bash 05_align.sh <candidates.fastq> <ref_genome> <output.bam> [threads]
set -e
INPUT=$1
REF=$2
OUTPUT=$3
THREADS=${4:-29}
MMI="${REF}.mmi"

echo "===== Step 5: minimap2 比对 ====="
echo "  输入: $INPUT"
echo "  参考: $REF"

if [ ! -f "$MMI" ]; then
    echo "  建索引: minimap2 -d $MMI $REF"
    minimap2 -d "$MMI" "$REF"
fi

minimap2 -x map-ont -k 11 -w 5 -n 2 -m 20 -s 20 \
    --min-dp-len 20 -a -I8G -t $THREADS --eqx --cs=long -Y \
    --secondary=yes -N50 -p0.5 \
    "$MMI" "$INPUT" \
    | samtools sort -@ $THREADS -o "$OUTPUT"

samtools index -@ $THREADS "$OUTPUT"
echo "===== 比对完成: $OUTPUT ====="
samtools flagstat "$OUTPUT"
