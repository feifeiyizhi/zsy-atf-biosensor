# zsy aTF 项目文档索引

> 2026-09-14 建,对齐 skill 约定(pi-briefing/paper-graphwalks/ec-target-discovery):
> **根** = 状态真相源(PROJECT_STATE)+ 冻结主张(CLAIMS)+ 总文件(PROJECT_BACKGROUND_PROTOCOL)
> **dashboard/** = 看板(PI 报告 + 流程图)
> **docs/** = 零碎文档归集(planning 背景/规划 + analysis 分析)
> **work/** = 代码 + 数据 + 数据挨着的 RESULTS(脚本依赖,勿动)

## 根(控制面,每次推进先读)
| 文件 | 作用 | 状态 |
|---|---|---|
| `PROJECT_STATE.md` | **状态唯一真相源**(定位/Claims/Decisions/Blockers/Registry/下一步) | 活,每轮更新 |
| `CLAIMS.md` | 冻结主张台账(已验证/复现/未建立/已反驳/不可复算) | 活,只更新 EVIDENCE |
| `PROJECT_BACKGROUND_PROTOCOL.md` | 总文件:背景+操作规程+§9 最新进度+§10 踩坑台账 | 活 |

## dashboard/(看板)
| 文件 | 作用 |
|---|---|
| `dashboard/PI_FULL_REPORT.md` | PI 全景汇报(11 节) |
| `dashboard/experiment_flow.mmd` | Mermaid 流程图 |

## docs/planning/(背景与规划,只读)
| 文件 | 作用 |
|---|---|
| `docs/planning/ZSY_ATF_OVERVIEW_20260912.md` | 对外交流全景版(自足,供张老师/师兄) |
| `docs/planning/LANDSCAPE_PLAN.md` | 数据拆解+三层景观方案落地 |
| `docs/planning/METHODOLOGY.md` | GraphWalks 式过程监督迁到蛋白景观的 step-by-step protocol |
| `docs/planning/BIO_DIRECTION.md` | v1 方向(从 11 篇文献反推) |
| `docs/planning/BIO_DIRECTION_v2.md` | v2(DMS 可验证图,三层递进) |
| `docs/planning/BIO_DIRECTION_v3_minwetlab.md` | v3(湿实验需求压到最小) |

## docs/analysis/(分析文档)
| 文件 | 作用 |
|---|---|
| `docs/analysis/zero_shot_cross_protein.md` | 零样本跨蛋白(28 克隆,旧 claim 已证伪见下) |
| `docs/analysis/zero_shot_rebuild_comparison.md` | 重建对比+bootstrap CI(修正版,以本表为准) |
| `docs/analysis/scaffold_condition_phenotype_table.md` | 骨架-条件-表型对照表(provenance) |
| `docs/analysis/EVOLUTION_GRAPHWALKS.md` | 真实进化图 GraphWalks(862 谷节点破 greedy) |
| `docs/analysis/GRAPHWALKS_FEASIBILITY.md` | GraphWalks 可行性 |

## work/results/(代码 + 数据 + 数据挨着的 RESULTS)
| 类 | 文件 | 说明 |
|---|---|---|
| 代码 | `work/results/sequence_gate.py` | **序列关**(解析单测+一致性阻断,建模前必跑) |
| 代码 | `work/results/zero_shot_features.py` | 零样本特征重建+bootstrap+反常检查(可复现) |
| 数据 | `work/results/mutation_activity_full.csv` | 28 克隆突变+fc+seq(建模用,已过序列关) |
| 数据 | `work/results/clone_esm2_3B_emb.npy` | ESM2 3B 嵌入(2560 维) |
| RESULTS | `work/results/RESULTS.md` | 338lib 负筛景观(109134 变体) |
| RESULTS | `work/results/per_target_six/RESULTS.md` | 6 靶点胆酸全景 |
| RESULTS | `work/results/recommender/RESULTS_annotated.md` | 推荐器+时序+结构→突变分组 |
| RESULTS | `work/results/structure/STRUCTURE.md` | 结构线(含 §10 自身 ESMFold) |
| RESULTS | `work/results/anclaci_model/RESULTS.md` | AncLacI 参照 |
| RESULTS | `work/results/per_target/RESULTS.md` | 2 靶点交叉 |

## 读序建议
1. `PROJECT_STATE.md`(项目在哪)— 2. `CLAIMS.md`(什么立住了)— 3. `PROJECT_BACKGROUND_PROTOCOL.md §9/§10`(最新进度+踩坑)— 4. `dashboard/PI_FULL_REPORT.md`(PI 故事)— 5. 按需读 docs/ + work/results/。

## ⚠ 以修正版为准的几条
- 零样本"双向翻转 −0.37→+0.37"**已证伪** → 以 `docs/analysis/zero_shot_rebuild_comparison.md` 为准(单向 CDCA→7kLCA +0.30,反向不翻,bootstrap 不显著)。
- 旧零样本特征脚本未存盘,其 ρ 数字不可复算;重建脚本 `work/results/zero_shot_features.py` 可复现。
- 建模前必跑 `work/results/sequence_gate.py`(解析单测+csv-vs-seq 一致性,对不上 exit 1)。

## 约定(可复用到其他项目:structgw/paper_gw/EC)
根 PROJECT_STATE.md + CLAIMS.md + 总文件;dashboard/ 看板;docs/{planning,analysis} 归集;work/ 代码+数据+数据挨着 RESULTS 不动。
