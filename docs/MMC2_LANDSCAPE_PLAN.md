# mmc2 结构—活性景观计划与当前卡点

> 项目：`zsy-atf-biosensor`  
> 数据来源：`/mnt/workspace/ray/0724/0815zsypre/mmc2.pdf` 及其本地配套数据  
> 状态：Phase 1 数据准备已完成；已生成统一表、分数据集表、kNN 图指标和验证报告。下一步进入 rugged landscape 复现。

## 一、目标

先使用 Meger et al. 2024（Cell Systems）补充材料中的 LacI/GalR 数据，建立可复现的“序列/突变—结构特征—DNA-binding/repression fitness”基线景观；验证方法后，再迁移到 zsy 204 aa aTF 胆酸传感器。

这里的活性标签应准确称为：

- LacO DNA-binding / repression fitness；
- NGS/FACS 得到的 enrichment proxy；
- 不是通用蛋白活性，也不能直接等同于胆酸传感器输出。

## 二、数据范围

### 1. LGF phylogenetic library

约 1158 个 extant/ancestral LacI/GalR family DBD 变体：

```text
LGF_log2_enrichment.csv
LGF_fold_enrichment.csv
LGF_raw_read_counts.csv
LGF_ref_lib.csv
```

主要标签：

```text
avg_log2_enrichment
sd_log2_enrichment
```

LGF 没有一个可当作真实 WT 的共同亲本。当前统一表不构造全局 consensus；LGF 只保留原始序列，并按等长度序列分层做坐标和图分析，不能解释成实验 WT 或突变效应。

### 2. EcLacI DMS library

约 1122 行，包含 WT 和 EcLacI DBD 单点突变：

```text
DMS_log2_normalized_enrichment.csv
DMS_fold_enrichment.csv
DMS_raw_read_counts.csv
DMS_ref_lib.csv
```

主要标签：

```text
avg_log2_normalized_enrichment
sd_log2_normalized_enrichment
```

DMS 以 `DMS_ref_lib.csv` 中的 WT 为参考，逐条解析 `K2A` 等单点突变。

## 三、已完成

- 已核对 `mmc2.pdf`：共 54 页，包含 LGF phylogenetic library、EcLacI DMS、实验方法及补充表图。
- 已确认本地配套数据存在于：

```text
/mnt/workspace/ray/0724/0815zsypre/work/anclaci/
```

- 已确认现有项目已有 AncLacI 初步数据和 baseline 模型，但本次新流程不覆盖旧结果。
- 已在公开项目中写入脚本：

```text
work/mmc2_prepare_dataset.py
work/mmc2_ruggedness.py
```

脚本功能：

1. 自动定位 LGF/DMS enrichment 表；
2. 清理序列；
3. 解析 DMS 单点突变；
4. 为 LGF 建立 consensus coordinate-only reference；
5. 生成统一变体表；
6. 建立序列 Hamming kNN 图；
7. 计算每个变体的局部峰、较高邻居数、节点度数；
8. 计算每个数据集的 graph Dirichlet energy；
9. 写入验证 JSON。

## 四、执行卡点（已解决）

数据准备脚本已成功运行。运行过程中发现 LGF 含有不同长度序列，因此脚本已改为按长度分层建图，并明确不构造 LGF 全局 WT consensus。详细结果见 `work/results/mmc2_landscape/RESULTS.md`。

脚本已实际运行并生成结果表；没有提交 Git，当前只保留工作区变更。

待执行的命令：

```bash
python /volume/schen04/ray/0724/_pub_zsy_atf/work/mmc2_prepare_dataset.py \
  --source /mnt/workspace/ray/0724/0815zsypre/work/anclaci \
  --out /volume/schen04/ray/0724/_pub_zsy_atf/work/results/mmc2_landscape
```

## 五、预期输出

```text
work/results/mmc2_landscape/
├── mmc2_unified_variants.csv
├── mmc2_dms_single_mutants.csv
├── mmc2_lgf_phylogenetic_variants.csv
├── mmc2_graph_metrics.csv
└── mmc2_validation.json
```

重点检查：

- 总行数是否约为 2280；
- DMS 非 WT 行是否均为单点突变；
- DMS WT 是否恰好一行；
- 两类数据的序列长度是否一致；
- 两类 fitness 是否保持独立尺度；
- kNN 图是否成功生成；
- local peak 数量和 Dirichlet energy 是否为有限值。

## 六、后续执行顺序

### Phase 1：数据验证

1. 运行 `mmc2_prepare_dataset.py`；
2. 检查 `mmc2_validation.json`；
3. 检查 DMS 单点解析；
4. 检查重复、缺失和异常长度；
5. 对照论文补充材料中的 1121/1158 规模描述。

### Phase 2：复现 rugged landscape

1. 复现 Hamming/kNN graph；
2. 计算局部峰和 valley 节点；
3. 计算 normalized Dirichlet energy；
4. 分别报告 LGF 和 DMS，不混合原始标签；
5. 生成景观图和可复现结果 Markdown。

### Phase 3：结构特征

优先使用已有 LacI 结构（论文相关结构包括 1EFA）提取：

- 与 DNA 距离；
- recognition helix / hinge helix 区域；
- dimer interface；
- 局部接触数；
- solvent accessibility；
- 位置保守性。

之后再对关键变体使用 ESMFold 或其他结构预测工具。第一版不为全部 LGF 变体盲目生成结构。

### Phase 4：结构—fitness 预测

分别建立：

- DMS 单点模型；
- LGF phylogenetic 模型；
- 联合模型。

必须使用：

- random split；
- position holdout；
- sequence-cluster holdout；
- LGF → DMS transfer；
- DMS → LGF transfer。

模型输出应包括预测均值和不确定性，而不是只有单一排序分数。

### Phase 5：生成式候选搜索

在结构—fitness 预测器通过严格留出验证后，再加入条件 latent diffusion 或 score-based search。生成候选需要经过：

1. 序列合法性；
2. 突变阶数；
3. 结构风险；
4. 训练分布距离；
5. 活性预测；
6. 不确定性；
7. 可实验合成性过滤。

## 七、科学解释边界

- mmc2 标签是 DNA repression fitness proxy；
- LGF consensus 不是实验 WT；
- 结构特征提高预测性能必须通过 held-out 验证证明；
- diffusion 生成的候选只能作为实验优先级排序，不能直接宣称为真实高活性突变；
- 从 LacI/GalR 到 zsy aTF 的迁移需要重新验证，不能默认跨蛋白泛化。
