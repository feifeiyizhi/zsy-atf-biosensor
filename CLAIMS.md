# 主张台账(结论关) — zsy aTF

> 2026-09-14 · 每条主张标状态 + 证据 + 数据版本/特征版本/切分/seed。
> 旧脚本未存盘的历史数字标 **不可复算**,不再靠猜测补解释。
> 建模用数据 = `mutation_activity_full.csv`(已过 `sequence_gate.py` 一致性门)。

## 状态定义
- **已验证(VERIFIED)**:有可复现代码/门证明。
- **复现(REPRODUCED)**:重建脚本跑出同向同量级数字。
- **未建立(NOT ESTABLISHED)**:证据不足或重建不成立。
- **已反驳(REFUTED)**:重建证伪。
- **不可复算(NOT RECOMPUTABLE)**:原脚本未存盘,历史数字无法重算。

## 数据版本登记
| 版本 | 文件 | 参考系 | 用途 | 门 |
|---|---|---|---|---|
| full | `mutation_activity_full.csv` | WT 相对(0810 含 A11) | **建模用** | sequence_gate PASS(28/28 csv==seq) |
| dataset | `mutation_activity_dataset.csv` | 0810 A11 相对(A3="none") | 仅 PPT 对照,**不建模** | 参考系不同,勿混 |
| WT libbackbone | `wt_reference_libbackbone.fasta` | 旧密码子(匹配师兄 338lib reads) | 蛋白骨架基准 | 单独记录 |
| WT alt-codon | `wt_reference_altcodon_userpasted.fasta` | 蛋白等价但 reads 不匹配库 | 仅记录,不用 | 蛋白==libbackbone(等价) |

## 主张台账
| # | 主张 | 状态 | 证据 | 数据/特征/切分/seed |
|---|---|---|---|---|
| 1 | 解析器 R116Q→116/I73N→73 正确 | **已验证** | `sequence_gate.py test_parser` + m[1:-3] 回归测试 | — |
| 2 | 全克隆 csv 突变 == seq aa 逐位一致 | **已验证** | `sequence_gate.py` 28/28 PASS(0810 csv 已含 A11,直接相等) | full.csv / seq 列 / vs WT libbackbone |
| 3 | 建模用 full.csv 参考系统一(无 A11 偏差) | **已验证** | 门显示 0810 "csv==seq" 非"补 A11 后一致"→ full.csv 已统一 | full.csv / 0810 含 A11 |
| 4 | "csv mutations 列两参考系污染建模" | **已反驳(我方错判)** | 修正:该问题在 `dataset.csv`/PPT 层,不在建模用 full.csv;之前对照表"问题4 最严重"过度判断,序列关纠正 | — |
| 5 | 28 克隆 LOO 有信号(+结构+胆酸 ρ≈0.40) | **复现** | `zero_shot_features.py` LOO +结构+胆酸=0.404(旧0.402) | full.csv / ESM2 3B npy / ANM 耦合 / Morgan1024 / Ridge / LOO / seed=见脚本 |
| 6 | LOO ESM2-only=−0.999(过拟合) | **已反驳** | 重建 Ridge(无 PCA)=−0.112;−0.999 是 PCA+Ridge 15 维过拟合假象 | — |
| 7 | 7k-LCA↔CDCA **双向**翻转 −0.37→+0.37 | **已反驳** | 重建仅单向:CDCA→7k-LCA −0.371→+0.296(翻);7k-LCA→CDCA −0.366→−0.070(不翻) | full.csv / 三档特征 / 4 条件切分 |
| 8 | 跨胆酸泛化(结构+配体是泛化关键) | **未建立** | 单向翻转+0.30 不足;**bootstrap CI 全跨0/极宽**:LOO+结构+胆酸 CI[−0.074,+0.722]✗不显著、跨实验CDCA→7kLCA CI[−1,+1]✗;胆酸特征条件内全常数→正ρ必来自条件分离非胆酸特征 | bootstrap 200次/单alpha=100/seed=0 |
| 8b | LOO +结构+胆酸+onehot ρ=0.414 | **复现(point)** | point=0.414 但 bootstrap CI[+0.018,+0.768]下限刚过0 | — |
| 9 | 23k/PCA 退步=origin 偏移+特征偏差双重 | **已反驳(部分)** | 三版位点算法一致→非特征偏差;退步是 origin 偏移(A11/Q136L 拉低 fc)+ 小样本噪声(n=6→6) | — |
| 9b | 胆酸特征能在组内改排名 | **已反驳(机制)** | 胆酸反常检查:4条件内胆酸特征全常数→纯加法线性模型不能靠它改组内排名;组内排名必来自ESM2/结构/位点,跨实验正ρ必来自条件分离 | bile_anomaly_check |
| 10 | 旧表 ρ 数字(双向+0.37 等) | **不可复算** | 原特征脚本未存盘;defensible 重建得单向+0.30 | — |
| 11 | Q164R 6 胆酸共享响应开关 | **已验证** | 6 靶点频率 16.6–99.8% + 口袋壁 3D=0Å + 自身 ESMFold pLDDT94.3 | per_target_six/RESULTS.md |
| 12 | 17 位点全 LBD、pocket↔HTH 33Å 长程别构 | **已验证** | 自身 ESMFold 32.7Å + 3CDL 33.3Å + ANM cross-corr −0.394 | STRUCTURE.md §10 |
| 13 | 7 位组 vs C23 组突变位点零重叠 | **已验证** | 6 胆酸 SMILES 确认后分组 | bile_smiles_confirmed.json |
| 14 | v2 路线在 AncLacI 参照跑通 | **复现(参照)** | 序列→活性 ρ0.76 + 别构 shift 0.72 | anclaci_model/RESULTS.md |
| 15 | A11 origin(0810)系统性拉低响应 | **已验证** | 12 克隆 fc 全挤 3.83–5.04;CDCA 唯一带 Q136L 的 A1 fc=3.85 最低档 | full.csv |

## 可复现性登记
| 脚本 | 存盘? | 产出 |
|---|---|---|
| `sequence_gate.py` | ✓ | 序列关门(解析单测+一致性阻断) |
| `zero_shot_features.py` | ✓ | LOO+跨实验 ρ(三版+bootstrap+反常检查,task10 补) |
| 旧零样本特征脚本 | ✗ | `clone_esm2_3B_emb.npy`+`zero_shot_cross_protein.md`(数字不可复算) |

## 启动门槛(扩展线,引你定的)
跨胆酸预测的启动门槛**未达**:
1. 取得跨骨架、跨配体可比较的**实测标签**(同骨架亲本+子代 × 无配体/目标/非目标胆酸 基础表达+剂量响应+重复)——需师兄湿实验。
2. 冻结基线 + 留出胆酸。
3. 比较"同一预测模型 有/无课题一方法"在新胆酸的**前瞻 top-k 命中率**。

**当前单向 ρ≈0.30 只用于选下一批候选,不充当门槛已过证据。**

## 修正记录
- 2026-09-14:序列关纠正我方错判#4(full.csv 已统一,非"两参考系污染建模");双向翻转#7 降为已反驳(重建仅单向)。
