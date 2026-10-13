#!/usr/bin/env python3
"""
Step 4: J58 相邻候选片段切分

师兄方法:
  - 仅 full motif 命中可作边界 (C24 和所有 partial 不切分)
  - 边界 ED 上限: L39<=4, J58<=8, E33/T30<=3, PCR23/BC0062<=2, R2/R3<=3
  - 重叠/嵌套的技术命中合并为一个区间
  - 从 J58 末端向外取到最近非重叠技术区间; 无边界取到 read 末端
  - 反向 J58 向左侧截取后反向互补, 统一方向
  - 候选长度 >= 50bp 才输出

用法: python3 04_cut_candidates.py <strict_j58.jsonl> <input.fastq.gz> <output.fastq>
"""
import gzip
import json
import sys
import time

from motifs import BOUNDARY_ED, reverse_complement

# 可作边界的 motif 及其 ED 上限
BOUNDARY_MOTIFS = set(BOUNDARY_ED.keys())


def get_boundary_hits(record):
    """提取可作边界的 full motif 命中 (仅 exact 模式, ED 在上限内)"""
    boundaries = []
    for h in record.get("exact", []):
        base = h["motif"].replace("_rc", "")
        if base in BOUNDARY_MOTIFS and h["ed"] <= BOUNDARY_ED[base]:
            boundaries.append((h["start"], h["end"], base))
    return boundaries


def merge_intervals(intervals):
    """合并重叠/嵌套区间"""
    if not intervals:
        return []
    sorted_ivs = sorted(intervals, key=lambda x: x[0])
    merged = [sorted_ivs[0]]
    for iv in sorted_ivs[1:]:
        if iv[0] <= merged[-1][1]:
            merged[-1] = (merged[-1][0], max(merged[-1][1], iv[1]))
        else:
            merged.append(iv)
    return merged


def cut_candidate(record, seq):
    """切出 J58 相邻的生物片段"""
    j58 = record.get("j58_hit")
    if not j58:
        return None

    j58_start = j58["start"]
    j58_end = j58["end"]
    is_forward = j58["direction"] == "forward"

    boundaries = merge_intervals(get_boundary_hits(record))

    if is_forward:
        # 正向: 片段在 J58 之后, 找最近的下游边界
        downstream = [b for b in boundaries if b[0] >= j58_end]
        if downstream:
            frag_end = min(b[0] for b in downstream)
        else:
            frag_end = len(seq)
        candidate = seq[j58_end:frag_end]
    else:
        # 反向: 片段在 J58 之前, 找最近的上游边界
        upstream = [b for b in boundaries if b[1] <= j58_start]
        if upstream:
            frag_start = max(b[1] for b in upstream)
        else:
            frag_start = 0
        candidate = reverse_complement(seq[frag_start:j58_start])

    if len(candidate) >= 50:
        return candidate
    return None


def process(strict_jsonl, fastq_path, output_fastq):
    open_fn = gzip.open if fastq_path.endswith(".gz") else open
    total = 0
    written = 0
    t0 = time.time()

    with open(strict_jsonl) as fjson, \
         open_fn(fastq_path, "rt") as fq, \
         open(output_fastq, "w") as fout:

        for line in fjson:
            record = json.loads(line)
            total += 1

            # 同步读 FASTQ
            header = fq.readline().strip()
            if not header:
                break
            seq = fq.readline().strip()
            fq.readline()
            fq.readline()

            candidate = cut_candidate(record, seq)
            if candidate:
                fout.write(f"@{record['read_id']}\n{candidate}\n+\n{'I'*len(candidate)}\n")
                written += 1

            if total % 500000 == 0:
                rate = total / (time.time() - t0)
                print(f"  已处理 {total:,} | 输出 {written:,} ({written/total*100:.1f}%) | {rate:.0f} reads/s")

    elapsed = (time.time() - t0) / 60
    print(f"\n切分完成! ({elapsed:.1f} min)")
    print(f"  严格 J58 reads: {total:,}")
    print(f"  输出候选 (>=50bp): {written:,} ({written/total*100:.1f}%)")


if __name__ == "__main__":
    if len(sys.argv) < 4:
        print("用法: python3 04_cut_candidates.py <strict_j58.jsonl> <input.fastq.gz> <output.fastq>")
        sys.exit(1)
    process(sys.argv[1], sys.argv[2], sys.argv[3])
