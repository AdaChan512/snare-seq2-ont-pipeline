#!/bin/bash
# Step 1: NanoPlot 质控 (不过滤, 师兄做法)
# 用法: bash 01_qc.sh <input.fastq.gz> <output_dir>
set -e
INPUT=$1
OUTDIR=${2:-"./01_qc"}
mkdir -p "$OUTDIR"

echo "===== Step 1: NanoPlot 质控 ====="
echo "  输入: $INPUT"
echo "  不做任何预过滤"

seqkit stats -a "$INPUT" > "${OUTDIR}/seqkit_stats.txt"
echo "  [1/3] seqkit stats -> ${OUTDIR}/seqkit_stats.txt"

NanoPlot --fastq "$INPUT" --N50 --loglength -t 29 \
    -o "${OUTDIR}/nanoplot_full" --named --plots hex dot kde
echo "  [2/3] NanoPlot -> ${OUTDIR}/nanoplot_full/NanoPlot-report.html"

seqkit fx2tab -l -n -q "$INPUT" | awk '{print $2"\t"$3}' > "${OUTDIR}/read_length_quality.tsv"
echo "  [3/3] 读长质量 -> ${OUTDIR}/read_length_quality.tsv"

echo "===== 质控完成 ====="
head -2 "${OUTDIR}/seqkit_stats.txt"
