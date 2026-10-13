#!/bin/bash
# 单样品完整流程
# 用法: bash run_all.sh <input.fastq.gz> <sample_name> <ref_genome.fa>
set -e
INPUT=$1
SAMPLE=$2
REF=$3
THREADS=${4:-29}
OUT="output/${SAMPLE}"

mkdir -p "${OUT}/"{01_qc,02_scan,03_j58,04_candidates,05_align,06_qc,07_classify,08_cluster,09_fragments}

echo "===== [1/9] QC ====="
bash 01_qc.sh "$INPUT" "${OUT}/01_qc"

echo "===== [2/9] 扫描 ====="
python3 02_scan.py "$INPUT" "${OUT}/02_scan/scan.jsonl"

echo "===== [3/9] J58 筛选 ====="
python3 03_j58_filter.py "${OUT}/02_scan/scan.jsonl" "${OUT}/03_j58/strict.jsonl"

echo "===== [4/9] 切分候选 ====="
python3 04_cut_candidates.py "${OUT}/03_j58/strict.jsonl" "$INPUT" "${OUT}/04_candidates/candidates.fastq"

echo "===== [5/9] 比对 ====="
bash 05_align.sh "${OUT}/04_candidates/candidates.fastq" "$REF" "${OUT}/05_align/aligned.bam" $THREADS

echo "===== [6/9] QC 过滤 ====="
python3 06_qc_filter.py "${OUT}/05_align/aligned.bam" "${OUT}/06_qc/qc_passed.bam"

echo "===== [7/9] 末端分类 ====="
python3 07_classify_ends.py "${OUT}/06_qc/qc_passed.bam" "${OUT}/03_j58/strict.jsonl" "${OUT}/07_classify/end_class.tsv"

echo "===== [8/9] 单元聚类 ====="
python3 08_cluster_units.py "${OUT}/03_j58/strict.jsonl" "${OUT}/08_cluster/clustered.jsonl"

echo "===== [9/9] 片段分析 ====="
python3 09_fragment_analysis.py "${OUT}/07_classify/end_class.tsv" "${OUT}/09_fragments"

echo "===== 全部完成: ${OUT} ====="
