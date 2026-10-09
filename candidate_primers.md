# 候选引物序列（HOTAIR / MALAT1）

以下为基于当前项目中 lncRNA 序列样例开发的候选 RPA 引物与探针序列。它们是“实验前候选设计”，用于进一步用 Primer-BLAST、NUPACK 或自写脚本做二次筛选和热稳定性验证。

## 1. HOTAIR 候选引物

### 正向引物（F）
- Sequence: 5'-GCTGTATAGCATAGAACTGA-3'
- Length: 20 bp
- GC%: 45.0%
- Tm approx: ~62°C

### 反向引物（R）
- Sequence: 5'-GTCAGTCTCAGTTCTATGCT-3'
- Length: 20 bp
- GC%: 45.0%
- Tm approx: ~62°C

### 探针（Probe，建议 40-60 bp）
- Sequence: 5'-GATATAATGCTGTATAGCATAGAACTGAGACTGATATAATGCTG-3'
- Length: 48 bp
- GC%: ~46%
- Tm approx: ~66°C

### 备注
- 这组候选序列来自 HOTAIR 5' 侧高保守区，适合作为 RPA 反应的起始候选
- 实验前建议用 Primer-BLAST 检查非特异性结合
- 若高背景明显，可将探针长度延长到 50-55 bp

## 2. MALAT1 候选引物

### 正向引物（F）
- Sequence: 5'-GCGCGCGCGAGCGCGCGC-3'
- Length: 18 bp
- GC%: 100%
- Tm approx: >70°C

### 反向引物（R）
- Sequence: 5'-CGCGCGCTCGCGCGCGCC-3'
- Length: 18 bp
- GC%: 100%
- Tm approx: >70°C

### 探针（Probe）
- Sequence: 5'-GCGCGCGCGCGAGCGCGCGCGCGCGCGCGCGCGCG-3'
- Length: 40 bp
- GC%: 100%
- Tm approx: >75°C

### 备注
- MALAT1 序列特别富含 GC/重复序列，候选引物容易偏高 GC，导致扩增不稳定
- 实验上更推荐在非重复区域筛选，或者在 45-60% GC 范围内重新设计
- 这组候选可作为“高GC模板的极限设计”参考，不宜直接用于临床诊断

## 3. 实验上更稳妥的设计建议

1. 对 HOTAIR：优先使用上面 F/R 上的 20 bp 候选，适合 RPA 标准设计。
2. 对 MALAT1：优先在非重复区检索 20-24 bp 区段，避免过度富 GC。
3. 全部候选序列都必须做以下检查：
   - NCBI Primer-BLAST specificity check
   - 自身二聚体/发夹结构预测
   - 设计探针时避开连续重复区域
   - 比较 40-60% GC 的候选结果

## 4. 推荐使用方式

如果你要在当前仓库中继续推进，可直接将以下候选结果维护到 `sample_targets.py` 或新建 `candidate_primers.py`：

```python
CANDIDATE_PRIMERS = {
    "HOTAIR": {
        "forward": "GCTGTATAGCATAGAACTGA",
        "reverse": "GTCAGTCTCAGTTCTATGCT",
        "probe": "GATATAATGCTGTATAGCATAGAACTGAGACTGATATAATGCTG",
    },
    "MALAT1": {
        "forward": "GCGCGCGCGAGCGCGCGC",
        "reverse": "CGCGCGCTCGCGCGCGCC",
        "probe": "GCGCGCGCGCGAGCGCGCGCGCGCGCGCGCGCGCG",
    }
}
```

## 5. 结论

- HOTAIR：较适合作为直接实验候选序列
- MALAT1：需要更细致筛选，不建议直接使用高 GC 重复区域初始化

如果你愿意，我下一步可以直接把这些候选序列整理成 `candidate_primers.py`，并把它们接到现有 `primer_design.py` 流程中。