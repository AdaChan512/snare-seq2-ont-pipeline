#!/usr/bin/env python3
"""
SNARE-seq2 ONT 分析 — 技术序列定义与阈值
所有步骤共享此模块，避免序列不一致。

来源: 师兄 PPTX (20260914_ONT_AC_plus.pptx) + 原文 SNARE-seq2 Protocol
"""
from math import floor, ceil

# ============================================================
# 技术序列 (5'→3')
# ============================================================
MOTIFS = {
    # ---- J58 复合结构 (C24 + L15 + ME19) ----
    "J58":       "ACGTACTGCAGACTATGTCTACAGCGGTAGCCGTCGCTCAGATGTGTATAAGAGACAG",  # 58nt
    "C24":       "ACGTACTGCAGACTATGTCTACAG",          # 24nt, AC linker 恒定区
    "L15":       "CGGTAGCCGTCGCTC",                    # 15nt, Redesigned adapter 2
    "ME19":      "AGATGTGTATAAGAGACAG",                 # 19nt, Tn5 mosaic end

    # ---- 远端技术序列 (3'端可能出现) ----
    "E33":       "CTGTCTCTTATACACATCTGACGCTGCCGACGA",   # 33nt, RC(Nextera_ad1)
    "T30":       "AAGCAGTGGTATCAACGCAGAGTGAATGGG",      # 30nt, TSO (师兄版本)
    "PCR23":     "AAGCAGTGGTATCAACGCAGAGT",             # 23nt, TSO PCR 引物区

    # ---- 其他技术序列 ----
    "L39":       "CCGAGCCCACGAGACCTGTAGACATAGTCTGCAGTACGT",  # 39nt, AC-R1-linker (BC_0100)
    "BC0062":    "CAGACGTGTGCTCTTCCGATCT",              # 22nt, Round3 5'恒定区
    "R2_linker": "CGAATGCTCTGGCCTCTCAAGCACGTGGAT",      # 30nt
    "R3_linker": "AGTCGTACGCCGATGCGAAACATCGGCCAC",      # 30nt
    "Ad1_blocker":"GACGCTGCCGACGAT",                     # 15nt
}

# ============================================================
# Partial (端部截断) 动态阈值
# 规则: m = max(15, ceil(0.8 * L)); ED_max = floor(0.1 * m)
#       blocker 特例: m = 12~14, ED_max = 1
# ============================================================
def partial_params(name, motif_len):
    """返回 (m, ed_max) 用于 partial 匹配"""
    if name == "Ad1_blocker":
        return 14, 1
    m = max(15, ceil(0.8 * motif_len))
    ed_max = floor(0.1 * m)
    return m, ed_max

# ============================================================
# J58 双层阈值
# ============================================================
J58_DISCOVERY_ED = 8       # 发现层: 召回用
C24_DISCOVERY_ED = 5       # 发现层辅助
J58_STRICT_ED = 6          # 严格层: J58 总 ED
C24_STRICT_ED = 3          # 严格层: C24
L15_STRICT_ED = 2          # 严格层: L15
ME19_STRICT_ED = 2         # 严格层: ME19

# ============================================================
# 切分边界 ED 上限 (仅 full motif 命中可作边界)
# ============================================================
BOUNDARY_ED = {
    "L39": 4, "J58": 8, "E33": 3, "T30": 3,
    "PCR23": 2, "BC0062": 2, "R2_linker": 3, "R3_linker": 3,
}

# ============================================================
# QC 过滤阈值
# ============================================================
MIN_ALN_LEN = 100
MIN_MAPQ = 30
MIN_IDENTITY = 0.90
MIN_COVERAGE = 0.80

# ============================================================
# 辅助函数
# ============================================================
def reverse_complement(seq):
    comp = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'N': 'N'}
    return ''.join(comp.get(b, 'N') for b in reversed(seq))

def build_all_motifs():
    """返回正向 + 反向互补的完整 motif 字典"""
    rc = {f"{k}_rc": reverse_complement(v) for k, v in MOTIFS.items()}
    return {**MOTIFS, **rc}
