#!/usr/bin/env python3
"""
Step 7: 末端结构 15 类互斥分类

师兄 PPTX 方法 (Slide 11-12):
  观察对象: 严格 J58 + 连续比对合格 reads
  观察区间: 按 J58→生物片段定向, 有效连续比对 query_end → 实际 read 末端
  核心: 完全位于观察区间内的 full E33 / T30
  规则: PCR23 嵌套于 T30 不重复计数; 跨起点命中单列;
        锚点不唯一/无观察序列/末端归属不明分别保留

15 类主类别:
  1.  核心 full 仅 T30
  2.  跨连续比对末端命中
  3.  有观察序列，未检出当前字典标记
  4.  末端归属不明
  5.  核心 full 仅独立 PCR23
  6.  仅 partial
  7.  多单元：无法评估或归属不明
  8.  核心 full 仅 E33
  9.  full E33 → RC(T30)
  10. 仅其他 full 标记
  11. 无可观察末端序列
  12. 多重 full E33 / T30 命中
  13. full E33＋T30：其他顺序或方向
  14. 多单元：组合不一致

用法: python3 07_classify_ends.py <qc_passed.bam> <strict_j58.jsonl> <output.tsv>
"""
import json
import sys
import time
from collections import Counter

try:
    import pysam
    import edlib
except ImportError as e:
    print(f"错误: {e}"); sys.exit(1)

from motifs import MOTIFS, reverse_complement, BOUNDARY_ED

# 观察区间内搜索的 full motif 及 ED 上限
FAR_MOTIFS = {
    "E33":   (MOTIFS["E33"], 3),
    "T30":   (MOTIFS["T30"], 3),
    "PCR23": (MOTIFS["PCR23"], 2),
    "L39":   (MOTIFS["L39"], 4),
    "BC0062":(MOTIFS["BC0062"], 2),
    "R2_linker": (MOTIFS["R2_linker"], 3),
    "R3_linker": (MOTIFS["R3_linker"], 3),
    "Ad1_blocker": (MOTIFS["Ad1_blocker"], 1),
}

ALL_FAR = {}
for name, (seq, ed) in FAR_MOTIFS.items():
    ALL_FAR[name] = (seq, ed)
    ALL_FAR[f"{name}_rc"] = (reverse_complement(seq), ed)


def find_full(region, motif_name, ed_max):
    """在 region 中搜索 full motif (HW 模式), 返回 (start, end, ed) 或 None"""
    seq, _ = FAR_MOTIFS[motif_name]
    # 正向
    r = edlib.align(seq, region, task="locations", mode="HW", k=ed_max)
    if 0 <= r["editDistance"] <= ed_max and r["locations"]:
        loc = r["locations"][0]
        return (loc[0], loc[1] + 1, r["editDistance"])
    # 反向互补
    rc = reverse_complement(seq)
    r = edlib.align(rc, region, task="locations", mode="HW", k=ed_max)
    if 0 <= r["editDistance"] <= ed_max and r["locations"]:
        loc = r["locations"][0]
        return (loc[0], loc[1] + 1, r["editDistance"])
    return None


def find_other_full(region):
    """搜索除 E33/T30/PCR23 外的其他 full 标记"""
    for name in ["L39", "BC0062", "R2_linker", "R3_linker", "Ad1_blocker"]:
        seq, ed = FAR_MOTIFS[name]
        r = edlib.align(seq, region, task="locations", mode="HW", k=ed)
        if 0 <= r["editDistance"] <= ed and r["locations"]:
            return name
    return None


def find_partial(region):
    """搜索 partial 标记 (端部截断)"""
    partials = []
    for name, (seq, ed_max) in FAR_MOTIFS.items():
        r = edlib.align(seq, region, task="locations", mode="HW", k=ed_max)
        if 0 <= r["editDistance"] <= ed_max:
            partials.append(name)
    return partials


def is_nested(inner, outer):
    """判断 inner 是否嵌套在 outer 内"""
    if inner is None or outer is None:
        return False
    return outer[0] <= inner[0] and inner[1] <= outer[1]


def classify_read(record, j58_hit, alignment):
    """
    单 read 分类, 返回主类别
    """
    # 确定观察区间: J58 远端
    seq = record.get("seq", "")
    if not seq:
        return "无可观察末端序列"

    j58_start = j58_hit["start"]
    j58_end = j58_hit["end"]
    is_forward = j58_hit["direction"] == "forward"

    if is_forward:
        far_region = seq[j58_end:]
    else:
        far_region = reverse_complement(seq[:j58_start])

    if not far_region:
        return "无可观察末端序列"

    # 搜索观察区间内的 full 标记
    full_e33 = find_full(far_region, "E33", 3)
    full_t30 = find_full(far_region, "T30", 3)
    full_pcr23 = find_full(far_region, "PCR23", 2)
    other_full = find_other_full(far_region)
    partials = find_partial(far_region)

    has_e33 = full_e33 is not None
    has_t30 = full_t30 is not None
    # PCR23 嵌套于 T30 时不独立计数
    has_pcr23 = full_pcr23 is not None and not is_nested(full_pcr23, full_t30)

    # === 互斥分类 (按优先级) ===

    # E33 + T30 组合
    if has_e33 and has_t30:
        # full E33 → RC(T30): E33 在前, T30 反向互补在后
        if full_e33[1] <= full_t30[0] and not full_t30[2]:
            return "full E33 → RC(T30)"
        # 其他顺序或方向
        if full_e33[0] > full_t30[1] or full_t30[2] > 0:
            return "full E33＋T30：其他顺序或方向"
        return "多重 full E33 / T30 命中"

    # 仅 T30
    if has_t30 and not has_e33:
        return "核心 full 仅 T30"

    # 仅独立 PCR23
    if has_pcr23 and not has_t30 and not has_e33:
        return "核心 full 仅独立 PCR23"

    # 仅 E33
    if has_e33 and not has_t30:
        return "核心 full 仅 E33"

    # 无 E33/T30/PCR23
    if not has_e33 and not has_t30 and not has_pcr23:
        if other_full:
            return "仅其他 full 标记"
        if partials:
            return "仅 partial"
        return "有观察序列，未检出当前字典标记"

    return "末端归属不明"


def process(input_bam, strict_jsonl, output_tsv):
    # 加载 strict j58 记录
    j58_records = {}
    with open(strict_jsonl) as f:
        for line in f:
            rec = json.loads(line)
            j58_records[rec["read_id"]] = rec

    counts = Counter()
    total = 0
    t0 = time.time()

    with pysam.AlignmentFile(input_bam, "rb") as bam, open(output_tsv, "w") as fout:
        fout.write("read_id\tcategory\n")
        for read in bam:
            if read.is_unmapped or read.is_secondary or read.is_supplementary:
                continue
            rid = read.query_name
            rec = j58_records.get(rid)
            if not rec or "j58_hit" not in rec:
                continue
            rec["seq"] = read.query_sequence or ""
            cat = classify_read(rec, rec["j58_hit"], read)
            counts[cat] += 1
            total += 1
            fout.write(f"{rid}\t{cat}\n")

            if total % 200000 == 0:
                print(f"  已分类 {total:,} ({(time.time()-t0)/60:.1f} min)")

    print(f"\n分类完成: {total:,} reads ({(time.time()-t0)/60:.1f} min)")
    print(f"\n===== 15 类互斥主表 (分母: {total:,}) =====")
    for cat, c in counts.most_common():
        print(f"  {cat:30s}: {c:>10,} ({c/total*100:.3f}%)")


if __name__ == "__main__":
    if len(sys.argv) < 4:
        print("用法: python3 07_classify_ends.py <qc_passed.bam> <strict_j58.jsonl> <output.tsv>")
        sys.exit(1)
    process(sys.argv[1], sys.argv[2], sys.argv[3])
