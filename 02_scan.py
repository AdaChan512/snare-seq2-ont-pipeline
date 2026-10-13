#!/usr/bin/env python3
"""
Step 2: 技术序列两级扫描 (Exact + Partial)

对应师兄 PPTX 方法:
  - Exact:  ED=0, 完整匹配
  - Partial: 端部截断匹配, 动态阈值 m=max(15,ceil(0.8L)), ED<=floor(0.1m)
             blocker 特例 m=12~14, ED<=1

输出 JSONL, 每行: {read_id, read_len, exact:[...], partial:[...]}

用法: python3 02_scan.py <input.fastq.gz> <output.jsonl> [max_reads]
"""
import gzip
import json
import sys
import time
from collections import defaultdict

try:
    import edlib
except ImportError:
    print("错误: pip install edlib"); sys.exit(1)

from motifs import MOTIFS, partial_params, reverse_complement

ALL_MOTIFS = {**MOTIFS, **{f"{k}_rc": reverse_complement(v) for k, v in MOTIFS.items()}}


def scan_exact(seq):
    """Exact = ED=0, str.find"""
    hits = []
    for name, motif in ALL_MOTIFS.items():
        start = 0
        while True:
            pos = seq.find(motif, start)
            if pos == -1:
                break
            hits.append({"motif": name, "start": pos, "end": pos + len(motif),
                         "ed": 0, "mode": "exact"})
            start = pos + 1
    return hits


def scan_partial(seq):
    """Partial = 端部截断匹配, 动态阈值"""
    hits = []
    seq_len = len(seq)
    for name, motif in ALL_MOTIFS.items():
        base = name.replace("_rc", "")
        m, ed_max = partial_params(base, len(motif))

        # 5'端: motif 前缀出现在 read 开头
        head_len = min(m + ed_max, seq_len)
        head = seq[:head_len]
        r = edlib.align(motif, head, task="locations", mode="HW", k=ed_max)
        if 0 <= r["editDistance"] <= ed_max:
            for loc in r["locations"]:
                if loc[0] <= ed_max:
                    hits.append({"motif": name, "start": loc[0], "end": loc[1] + 1,
                                 "ed": r["editDistance"], "mode": "partial", "end": "5'"})

        # 3'端: motif 后缀出现在 read 结尾
        tail_start = max(0, seq_len - m - ed_max)
        tail = seq[tail_start:]
        r = edlib.align(motif, tail, task="locations", mode="HW", k=ed_max)
        if 0 <= r["editDistance"] <= ed_max:
            for loc in r["locations"]:
                abs_end = tail_start + loc[1] + 1
                if abs_end >= seq_len - ed_max:
                    hits.append({"motif": name, "start": tail_start + loc[0],
                                 "end": abs_end, "ed": r["editDistance"],
                                 "mode": "partial", "end": "3'"})
    return hits


def process(fastq_path, output_path, max_reads=None):
    open_fn = gzip.open if fastq_path.endswith(".gz") else open
    total = 0
    exact_counts = defaultdict(int)
    partial_counts = defaultdict(int)
    t0 = time.time()

    with open_fn(fastq_path, "rt") as fin, open(output_path, "w") as fout:
        while True:
            header = fin.readline().strip()
            if not header:
                break
            seq = fin.readline().strip()
            fin.readline()
            fin.readline()
            read_id = header[1:].split()[0]
            total += 1

            exact_hits = scan_exact(seq)
            partial_hits = scan_partial(seq)

            for h in exact_hits:
                exact_counts[h["motif"].replace("_rc", "")] += 1
            for h in partial_hits:
                partial_counts[h["motif"].replace("_rc", "")] += 1

            fout.write(json.dumps({
                "read_id": read_id, "read_len": len(seq),
                "exact": exact_hits, "partial": partial_hits,
            }) + "\n")

            if total % 100000 == 0:
                rate = total / (time.time() - t0)
                print(f"  已处理 {total:,} reads ({rate:.0f} reads/s)")

            if max_reads and total >= max_reads:
                break

    print(f"\n扫描完成: {total:,} reads, 用时 {(time.time()-t0)/60:.1f} min")
    print(f"\n===== Exact 阳性 (占全部reads) =====")
    for name, c in sorted(exact_counts.items(), key=lambda x: -x[1]):
        print(f"  {name:15s}: {c:>10,} ({c/total*100:6.2f}%)")
    print(f"\n===== Partial 阳性 (占全部reads) =====")
    for name, c in sorted(partial_counts.items(), key=lambda x: -x[1]):
        print(f"  {name:15s}: {c:>10,} ({c/total*100:6.2f}%)")


if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("用法: python3 02_scan.py <input.fastq.gz> <output.jsonl> [max_reads]")
        sys.exit(1)
    max_reads = int(sys.argv[3]) if len(sys.argv) > 3 else None
    process(sys.argv[1], sys.argv[2], max_reads)
