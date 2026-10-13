#!/usr/bin/env python3
"""
Step 8: 单元簇聚类

师兄方法:
  重叠 J58 假设先聚类; 簇间组合一致计一次, 冲突或无法归属另列。
  同一 read 可含多个单元 (J58 命中), 需先聚类再汇总末端组合。

用法: python3 08_cluster_units.py <strict_j58.jsonl> <output.jsonl>
"""
import json
import sys
import time
from collections import Counter


def cluster_j58_hits(record, gap=20):
    """
    对一条 read 的所有 J58 命中按位置聚类
    相邻 J58 间距 <= gap 视为同一簇
    返回: [(cluster_id, [j58_hits]), ...]
    """
    hits = [h for h in record.get("exact", []) + record.get("partial", [])
            if h["motif"].replace("_rc", "") == "J58"]
    if not hits:
        return []

    hits.sort(key=lambda h: h["start"])
    clusters = []
    current = [hits[0]]

    for h in hits[1:]:
        if h["start"] - current[-1]["end"] <= gap:
            current.append(h)
        else:
            clusters.append(current)
            current = [h]
    clusters.append(current)
    return clusters


def process(input_jsonl, output_jsonl):
    total = 0
    multi_unit = 0
    cluster_dist = Counter()
    t0 = time.time()

    with open(input_jsonl) as fin, open(output_jsonl, "w") as fout:
        for line in fin:
            record = json.loads(line)
            total += 1
            clusters = cluster_j58_hits(record)
            record["unit_clusters"] = [
                {"cluster_id": i, "n_hits": len(c),
                 "start": c[0]["start"], "end": c[-1]["end"]}
                for i, c in enumerate(clusters)
            ]
            n_clusters = len(clusters)
            cluster_dist[n_clusters] += 1
            if n_clusters > 1:
                multi_unit += 1

            fout.write(json.dumps(record) + "\n")

            if total % 500000 == 0:
                print(f"  已处理 {total:,} | 多单元 {multi_unit:,} ({multi_unit/total*100:.1f}%) | {(time.time()-t0)/60:.1f} min")

    print(f"\n聚类完成! ({(time.time()-t0)/60:.1f} min)")
    print(f"  总 reads: {total:,}")
    print(f"  多单元 reads: {multi_unit:,} ({multi_unit/total*100:.2f}%)")
    print(f"\n  单元簇数分布:")
    for n, c in sorted(cluster_dist.items()):
        print(f"    {n} 簇: {c:,} ({c/total*100:.2f}%)")


if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("用法: python3 08_cluster_units.py <strict_j58.jsonl> <output.jsonl>")
        sys.exit(1)
    process(sys.argv[1], sys.argv[2])
