# mmc2 Phase 1 数据准备结果

运行脚本：

```bash
python work/mmc2_prepare_dataset.py \
  --source /mnt/workspace/ray/0724/0815zsypre/work/anclaci \
  --out work/results/mmc2_landscape
```

## 结果

| 项目 | 结果 |
|---|---:|
| 统一变体总数 | 2280 |
| LGF phylogenetic library | 1158 |
| EcLacI DMS library | 1122 |
| DMS WT 行 | 1 |
| DMS 非单点突变行（排除 WT） | 0 |
| 无效/空序列行 | 0 |
| kNN 参数 | 34 |
| DMS local peaks | 38 |
| LGF local peaks | 29 |
|

Fitness 范围：

```text
DMS: -9.9658 ~ 0.9322
LGF: -9.9658 ~ 6.9262
```

## 输出文件

```text
work/results/mmc2_landscape/mmc2_unified_variants.csv
work/results/mmc2_landscape/mmc2_dms_single_mutants.csv
work/results/mmc2_landscape/mmc2_lgf_phylogenetic_variants.csv
work/results/mmc2_landscape/mmc2_graph_metrics.csv
work/results/mmc2_landscape/mmc2_validation.json
```

## 数据处理修正

初次运行发现 LGF 序列长度并不统一：

```text
49 aa: 1
58 aa: 22
60 aa: 2
61 aa: 153
62 aa: 45
63 aa: 935
```

因此没有继续强行构造 LGF 全局 consensus，也没有把 consensus 当作实验 WT。当前脚本将 LGF 的 `reference_name` 标记为：

```text
LGF_no_common_reference
```

LGF 图分析按序列长度分层构建 Hamming kNN 图；不同长度的序列不直接计算 Hamming 距离。这个处理避免了把 indel/长度差异误当成普通点突变。

DMS 使用明确的 EcLacI WT reference，并成功解析为 1121 条单点突变 + 1 条 WT。

## 解释边界

- local peak 是当前 34-nearest-neighbor 图上的计算结果，不等于论文原始定义的所有进化峰；
- LGF 与 DMS 的 enrichment 标签保持分开，不能直接比较绝对值；
- DMS 与 LGF 可能存在实验体系和参考尺度差异；
- 下一步需要复现论文使用的图构建、ruggedness/Dirichlet energy 定义，并和当前 kNN 基线区分报告；
- 结构特征与预测模型尚未开始，当前结果仅完成数据准备和景观图基线。
