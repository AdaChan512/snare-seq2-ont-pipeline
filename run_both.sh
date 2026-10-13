#!/bin/bash
# 双样品 (++ vs +-) 批量运行
# 用法: bash run_both.sh <++.fastq.gz> <+-.fastq.gz> <ref_genome.fa>
set -e
PP=$1
PM=$2
REF=$3
THREADS=${4:-29}

echo "===== 样品 ++ (加 TSO) ====="
bash run_all.sh "$PP" "pp" "$REF" $THREADS

echo ""
echo "===== 样品 +- (不加 TSO) ====="
bash run_all.sh "$PM" "pm" "$REF" $THREADS

echo ""
echo "===== 双样品完成 ====="
echo "++ 结果: output/pp/"
echo "+- 结果: output/pm/"
