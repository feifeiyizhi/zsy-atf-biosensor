# AncLacI 参照模型结果（预测 unseen 突变体活性 + ±小分子别构）

日期 2026-09-11。参照 Meger et al. 2024 Cell Systems（mmc2.pdf）的 LacI/GalR DBD 景观数据（与我们的 TetR 同超家族、同 DBD+LBD 架构、同小分子别构调 DNA 结合机制）。
数据在 `work/anclaci/anclaci-ref/`，产物在 `work/results/anclaci_model/`。

## 数据（Zenodo 10652076）
| 集 | n | 目标 | 含义 |
|---|---|---|---|
| LGF 景观 | 1,158 DBD | avg_log2_enrichment (DNA结合fitness) | 稀疏 rugged：仅 22 功能(>0)，892 非功能 |
| DMS 单点 | 1,122 EcLacI 全单点 | log2_normalized_enrichment | 稠密局部景观 |
| Hill 别构 | 54 克隆×8 IPTG 剂量 | GFP（抑制/DNA结合读数） | **±小分子维度**（别构剂量响应） |

## 基线模型（任务#9）：预测 unseen 突变体 DNA 结合活性
特征：ESM2 t30_150M mean-pool 嵌入（640维）。模型：Ridge / GBM / MLP。
- **DMS 内部 5 折 CV**（稠密单点）：Ridge ρ=0.67，**GBM ρ=0.76** —— 稠密景观预测好。
- **LGF 内部 5 折 CV**（稀疏系统发育）：GBM ρ=0.31；但**功能变体 top-k 命中率 0.91–1.00** —— ρ 一般但能精准排出 22 个针（实用信号）。
- **跨集 LGF→DMS**：ρ 0.20–0.26（弱）—— **跨序列区不互迁**，与我们项目"负筛景观↔胆酸克隆不相交"一致。
- 结论：**密度决定质量**；挑针能力强；跨区迁移弱（景观区域特异性坐实）。

## 别构模型（任务#10）：预测 ±小分子 DNA 结合变化（v2 核心）
特征同上。Ridge 多输出 LOO-CV on 54 克隆。
- **别构 shift（IPTG饱-无）Spearman = 0.72** ⭐ —— "小分子有/无 DNA 结合变化"**可预测**，这正是用户要的核心目标。
- 剂量响应曲线形状 ρ=0.54；Hill 下限(bot 基础抑制) 0.70、上限(top 去抑制) 0.43；**EC50=0.03 预测不了**（精确效力是难点，n=54 太小，需更多数据/结构）。
- 图 `fig_allosteric_pred.png`。

## 结论与下一步
1. v2 路线在参照系统（LacI/IPTG）上**端到端跑通**：序列→unseen 突变体 DNA 结合活性（稠密区 0.76）+ 序列→±小分子别构效应（shift 0.72）。
2. **精确 EC50 是难点** → 需要：(a) 更多剂量响应数据；(b) 加结构特征（ESMFold，GPU 已授权待接集群）。
3. 迁移到我们 aTF（胆酸）：用变体面板±胆酸 FACS 数据 fine-tune；ESMFold 自身结构 + ANM 别构通道作结构特征/先验。
4. 与我们负筛景观的关系：本模型补上了我们缺的"±小分子维度"训练范式；我们 338lib 是 DNA 结合（无配体）那一半，正好对应 LGF/DMS；胆酸±对应 Hill。

## 产物
- `lgf_landscape.csv` `dms_landscape.csv` `hill_allosteric.csv`
- `esm2_emb_{lgf,dms,hill}.npy`、`hill_pred.npy` `hill_true.npy`
- `fig_allosteric_pred.png`
- 脚本：`anclaci_eda.py` `baseline_model.py` `allosteric_model.py`
