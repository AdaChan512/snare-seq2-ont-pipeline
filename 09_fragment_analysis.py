#!/usr/bin/env python3
"""
Step 9: 片段大小分布分析

按末端结构分组, 统计片段长度分布, 标注核小体周期 (200/400/600/800/1000bp)

用法: python3 09_fragment_analysis.py <end_class.tsv> <output_dir>
"""
import os
import sys
import csv
from collections import defaultdict

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

NUCLEOSOME = [200, 400, 600, 800, 1000]


def setup_font():
    import platform
    if platform.system() == "Windows":
        plt.rcParams['font.sans-serif'] = ['Microsoft YaHei', 'SimHei']
    else:
        plt.rcParams['font.sans-serif'] = ['Noto Sans CJK SC']
    plt.rcParams['axes.unicode_minus'] = False


def main(tsv_path, output_dir):
    setup_font()
    os.makedirs(output_dir, exist_ok=True)

    sizes_by_cat = defaultdict(list)
    with open(tsv_path) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            # 片段长度来自 alignment length, 这里从 tsv 的 fragment_len 列取
            # 如果没有, 用 read_len 估算
            flen = int(row.get("fragment_len", row.get("aln_len", 0)))
            if flen > 0:
                sizes_by_cat[row["category"]].append(flen)

    # 统计表
    print(f"{'类别':30s} {'n':>10s} {'中位数':>8s} {'均值':>8s} {'<200':>8s} {'200-400':>8s}")
    print("-" * 80)
    stats_rows = []
    for cat, sizes in sorted(sizes_by_cat.items(), key=lambda x: -len(x[1])):
        arr = np.array(sizes)
        row = [cat, len(arr), f"{np.median(arr):.0f}", f"{np.mean(arr):.0f}",
               f"{np.sum(arr<200)/len(arr)*100:.1f}%",
               f"{np.sum((arr>=200)&(arr<400))/len(arr)*100:.1f}%"]
        stats_rows.append(row)
        print(f"  {cat:28s} {row[1]:>10,} {row[2]:>8s} {row[3]:>8s} {row[4]:>8s} {row[5]:>8s}")

    with open(os.path.join(output_dir, "fragment_stats.csv"), "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["category", "n", "median", "mean", "lt200_pct", "200_400_pct"])
        w.writerows(stats_rows)

    # 直方图 (top 6 类别)
    top_cats = sorted(sizes_by_cat, key=lambda c: -len(sizes_by_cat[c]))[:6]
    fig, axes = plt.subplots(2, 3, figsize=(18, 10))
    for idx, cat in enumerate(top_cats):
        ax = axes.flatten()[idx]
        sizes = [s for s in sizes_by_cat[cat] if s <= 2000]
        ax.hist(sizes, bins=100, range=(0, 2000), alpha=0.7, edgecolor="white")
        ax.set_title(f"{cat} (n={len(sizes_by_cat[cat]):,})", fontsize=11)
        ax.set_xlabel("Fragment length (bp)")
        for x in NUCLEOSOME:
            ax.axvline(x, color="red", ls="--", alpha=0.3)
    plt.suptitle("Fragment Size Distribution by End Category", fontsize=14)
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "fragment_histograms.png"), dpi=150)
    plt.close()
    print(f"\n图表: {output_dir}/fragment_histograms.png")


if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("用法: python3 09_fragment_analysis.py <end_class.tsv> <output_dir>")
        sys.exit(1)
    main(sys.argv[1], sys.argv[2])
