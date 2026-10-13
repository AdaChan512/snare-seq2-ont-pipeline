# SNARE-seq2 ONT 分析 Pipeline

升主动脉单细胞 AC+mRNA 联合测序 (SNARE-seq2) 的 Nanopore 数据分析流程。
对齐师兄 PPTX (20260914_ONT_AC_plus.pptx) 方法学。

## 环境依赖

```bash
# Python 包
pip install edlib pysam pandas numpy matplotlib biopython

# 生物信息工具 (conda)
conda install -c bioconda nanoplot seqkit minimap2 samtools
```

## 流程步骤

| 步骤 | 脚本 | 说明 |
|---|---|---|
| 1 | `01_qc.sh` | NanoPlot 质控 (不过滤) |
| 2 | `02_scan.py` | 两级扫描: Exact + Partial (动态阈值) |
| 3 | `03_j58_filter.py` | J58 双层阈值: 发现 (J58≤8) / 严格 (J58≤6+C24≤3/L15≤2/ME19≤2) |
| 4 | `04_cut_candidates.py` | 候选片段切分 (仅 full motif 边界, ≥50bp) |
| 5 | `05_align.sh` | minimap2 比对 (map-ont, GRCh38) |
| 6 | `06_qc_filter.py` | 连续比对 QC (≥100bp, MAPQ≥30, identity≥90%, coverage≥80%) |
| 7 | `07_classify_ends.py` | 15 类互斥末端结构分类 |
| 8 | `08_cluster_units.py` | 单元簇聚类 |
| 9 | `09_fragment_analysis.py` | 片段大小分布分析 |

## 快速开始

```bash
# 单样品
bash run_all.sh <sample.fastq.gz> <sample_name> <ref_genome.fa>

# 双样品 (++ vs +-)
bash run_both.sh <++.fastq.gz> <+-.fastq.gz> <ref_genome.fa>
```

## 关键阈值 (对应师兄)

- 严格 J58 阳性率: ~31.27%
- 连续比对通过率: ~48.94%
- 15 类主表分母: 严格 J58 + 连续比对合格 reads
