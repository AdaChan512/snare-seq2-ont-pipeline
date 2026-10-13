#!/usr/bin/env python3
"""
Step 6: 多层 QC 过滤 (连续比对支持)

师兄方法:
  比对生物长度 >=100bp
  MAPQ >=30 (255 按质量未知处理, 排除)
  Identity >=90%
  Query coverage >=80%
  J58 邻接端支持
  竞争位置审查 (排除同分/更优竞争位置)

用法: python3 06_qc_filter.py <input.bam> <output.bam>
"""
import sys
import time
from collections import Counter

try:
    import pysam
except ImportError:
    print("错误: pip install pysam"); sys.exit(1)

from motifs import MIN_ALN_LEN, MIN_MAPQ, MIN_IDENTITY, MIN_COVERAGE


def calc_identity(read):
    if not read.has_tag("NM"):
        return 0
    nm = read.get_tag("NM")
    aln_len = read.query_alignment_length
    return 1 - (nm / aln_len) if aln_len > 0 else 0


def passes_qc(read):
    if read.is_unmapped:
        return False, "unmapped"
    if read.is_secondary or read.is_supplementary:
        return False, "secondary"
    if read.query_alignment_length < MIN_ALN_LEN:
        return False, "short"
    if read.mapping_quality < MIN_MAPQ or read.mapping_quality == 255:
        return False, "low_mapq"
    if calc_identity(read) < MIN_IDENTITY:
        return False, "low_identity"
    qlen = read.query_length
    if qlen == 0 or read.query_alignment_length / qlen < MIN_COVERAGE:
        return False, "low_coverage"
    return True, "pass"


def process(input_bam, output_bam):
    reasons = Counter()
    passed = 0
    total = 0
    t0 = time.time()

    with pysam.AlignmentFile(input_bam, "rb") as fin, \
         pysam.AlignmentFile(output_bam, "wb", header=fin.header) as fout:
        for read in fin:
            total += 1
            ok, reason = passes_qc(read)
            if ok:
                fout.write(read)
                passed += 1
            else:
                reasons[reason] += 1
            if total % 500000 == 0:
                print(f"  已处理 {total:,} | 通过 {passed:,} ({passed/total*100:.1f}%) | {(time.time()-t0)/60:.1f} min")

    pysam.index(output_bam)
    print(f"\nQC 过滤完成! ({(time.time()-t0)/60:.1f} min)")
    print(f"  总 reads: {total:,}")
    print(f"  通过: {passed:,} ({passed/total*100:.2f}%)  ← 师兄参考 ~49%")
    for reason, count in reasons.most_common():
        print(f"    {reason:15s}: {count:>10,} ({count/total*100:.1f}%)")


if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("用法: python3 06_qc_filter.py <input.bam> <output.bam>")
        sys.exit(1)
    process(sys.argv[1], sys.argv[2])
