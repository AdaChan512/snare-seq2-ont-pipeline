#!/usr/bin/env python3
"""
Step 3: J58 双层阈值筛选 + 严格结构确认

师兄方法:
  发现层 (召回):   J58 ED<=8  或  C24 ED<=5
  严格层 (确认):   J58 ED<=6  且  C24<=3 / L15<=2 / ME19<=2
  → 严格层对应师兄 31.27%

输出: j58_strict.jsonl (仅严格层通过的 reads, 带 j58_hit 标注)

用法: python3 03_j58_filter.py <scan.jsonl> <output.jsonl>
"""
import json
import sys
import time
from collections import Counter

from motifs import (J58_DISCOVERY_ED, C24_DISCOVERY_ED,
                    J58_STRICT_ED, C24_STRICT_ED, L15_STRICT_ED, ME19_STRICT_ED)


def get_min_ed(hits, motif_name):
    """获取某 motif 的最小 ED"""
    eds = [h["ed"] for h in hits if h["motif"].replace("_rc", "") == motif_name]
    return min(eds) if eds else 999


def filter_j58_dual(record):
    """
    返回 (is_discovery, is_strict, j58_hit)
    """
    all_hits = record.get("exact", []) + record.get("partial", [])

    j58_hits = [h for h in all_hits if h["motif"].replace("_rc", "") == "J58"]
    if not j58_hits:
        return False, False, None

    best_j58 = min(j58_hits, key=lambda h: h["ed"])
    j58_ed = best_j58["ed"]

    c24_ed = get_min_ed(all_hits, "C24")
    l15_ed = get_min_ed(all_hits, "L15")
    me19_ed = get_min_ed(all_hits, "ME19")

    # 发现层
    is_discovery = (j58_ed <= J58_DISCOVERY_ED) or (c24_ed <= C24_DISCOVERY_ED)

    # 严格层: 必须同时满足
    is_strict = (j58_ed <= J58_STRICT_ED and
                 c24_ed <= C24_STRICT_ED and
                 l15_ed <= L15_STRICT_ED and
                 me19_ed <= ME19_STRICT_ED)

    return is_discovery, is_strict, best_j58


def process(input_jsonl, output_jsonl):
    total = 0
    discovery = 0
    strict = 0
    strict_ed_dist = Counter()
    t0 = time.time()

    with open(input_jsonl) as fin, open(output_jsonl, "w") as fout:
        for line in fin:
            record = json.loads(line)
            total += 1

            is_disc, is_strict, j58_hit = filter_j58_dual(record)

            if is_disc:
                discovery += 1
            if is_strict:
                strict += 1
                strict_ed_dist[j58_hit["ed"]] += 1
                record["j58_hit"] = {
                    "motif": j58_hit["motif"],
                    "start": j58_hit["start"],
                    "end": j58_hit["end"],
                    "ed": j58_hit["ed"],
                    "direction": "forward" if not j58_hit["motif"].endswith("_rc") else "reverse",
                }
                fout.write(json.dumps(record) + "\n")

            if total % 500000 == 0:
                rate = total / (time.time() - t0)
                print(f"  已处理 {total:,} | 发现 {discovery:,} ({discovery/total*100:.1f}%) "
                      f"| 严格 {strict:,} ({strict/total*100:.2f}%) | {rate:.0f} reads/s")

    elapsed = (time.time() - t0) / 60
    print(f"\nJ58 筛选完成! ({elapsed:.1f} min)")
    print(f"  总 reads:     {total:,}")
    print(f"  发现层阳性:   {discovery:,} ({discovery/total*100:.2f}%)")
    print(f"  严格层阳性:   {strict:,} ({strict/total*100:.2f}%)  ← 师兄参考 31.27%")
    print(f"\n  严格层 J58 ED 分布:")
    for ed, c in sorted(strict_ed_dist.items()):
        print(f"    ED={ed}: {c:,} ({c/strict*100:.1f}%)")


if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("用法: python3 03_j58_filter.py <scan.jsonl> <output.jsonl>")
        sys.exit(1)
    process(sys.argv[1], sys.argv[2])
