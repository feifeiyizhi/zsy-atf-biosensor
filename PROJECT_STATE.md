# PROJECT_STATE — zsy aTF 胆酸生物传感器(轻量控制面)

> 2026-09-13 建(项目此前无控制面);2026-09-14 加进化图 GraphWalks。每项标注状态 + 证据路径 + 最后核实日期。
> 详细见同目录 `PI_FULL_REPORT.md`。各 RESULTS.md 为原始证据。

## 一、项目定位
- 张老师命题、师兄出数据。TetR 家族变构转录因子(aTF,204aa)生物传感器定向进化建模。
- v2 策略(2026-09-11):合成生物学开关元件;顺序两步(定 DNA 结合面板→逐小分子筛别构→AI);建模必须基于结构+别构。
- 与博士主线关系:**待 PI 决策是否并入**(见 `phd-projects-overview`,#4)。

## 二、Claims 状态(见 PI_FULL_REPORT §6)
| Claim | 判定 | 证据 |
|---|---|---|
| Q164R 6 胆酸共享开关 | SUPPORTED | per_target_six/RESULTS.md |
| pocket↔HTH 33Å 长程别构 | SUPPORTED | STRUCTURE.md §10 |
| 7位组 vs C23 组位点零重叠 | SUPPORTED | RESULTS_annotated.md §7 |
| 跨实验翻转=结构+配体泛化关键 | **PRELIMINARY** | zero_shot_cross_protein.md(待 CI) |
| 通用 PLM 预测不了特异性 | SUPPORTED | STRUCTURE.md §6 |
| 23k/PCA 跨实验泛化 | NOT ESTABLISHED | origin 偏移 |
| 靶点特异次级预测 | NOT ESTABLISHED | 时序 top-50 召回 1/14 |
| v2 参照跑通 | SUPPORTED(参照) | anclaci_model/RESULTS.md |
| 真实进化图破 greedy(GraphWalks 底物) | **PRELIMINARY(温和)** | EVOLUTION_GRAPHWALKS.md(862 谷节点,置换 z=2.25/p=0.012,超随机+7%)|
| LacI 结构特征升级预测器 | **REFUTED(收束)** | 防泄漏 ΔR²≈−0.02 CI 跨零(structgw Phase D/F)|

## 三、Decisions
- 2026-09-04:从 GB1 假设数据重定位到师兄真实 17 位点饱和库。
- 2026-09-11:v2 策略修订(顺序两步 + 结构/别构建模)。
- 2026-09-13:自身 ESMFold 结构完成(4090,pLDDT 94.3),4090 可关。
- 2026-09-13:per-target 突变数修正(Q136L origin 归类错误,去骨架后梯度翻转)。
- 2026-09-14:边定义 Hamming-1→突变阶相邻(subset+1)重设计,进化图 80.8% 连通、破 greedy 谷节点超随机(温和)。结构特征升级预测器支线(LacI ESMFold)判 REFUTED 收束,与进化图过程监督分开讲。

## 四、Blockers
| 卡点 | 状态 | 解封条件 |
|---|---|---|
| 28 样本无统计置信度 | 🔴 未做 | bootstrap/permutation(不需 GPU) |
| lib0 起始随机池未测序 | 🔴 湿实验 | 师兄补测 → 真 round-0 基线,重估破 greedy 强度 |
| 破 greedy 信号温和(+7%) | 🟡 已量化 | 精选强 valley 子集(降幅+升幅大)提信噪 + 补 lib0 |
| 7 位组训练数据缺失 | 🟡 诊断完 | 补至少 1 个 7 位修饰胆酸进训练集 |
| 胆酸 SMILES 2/6 缺立体 | 🟡 4/6 已查 | 23K-CDCA + PCA 立体,待 PPT |
| 正/负胆酸 FACS 数据 | 🔴 湿实验 | 师兄侧 |
| 无 git/无控制面 | 🔴 本次部分补 | git init 待 PI 决策 |

## 五、Experiment Registry(简)
| 实验 | 数据 | 主结果 | 日期 | 路径 |
|---|---|---|---|---|
| 338lib 负筛景观 | 109134 变体/6840 可靠 | 17 位点 9-10 阶,PCA 16.1% | 09-08 | results/RESULTS.md |
| 6 靶点胆酸全景 | 6 靶点 reads | Q164R 全共享开关 | 09-12 | per_target_six/ |
| 自身 ESMFold 结构 | WT+3 变体 | pLDDT94.3/ANM-0.394 | 09-13 | structure/STRUCTURE.md §10 |
| 时序突变推荐器 | 2→4 靶点 | AUC0.777 | 09-13 | recommender/ |
| 零样本跨蛋白 | 28 克隆 3 实验 | 7kLCA↔CDCA 翻转 | 09-13 | zero_shot_cross_protein.md |
| AncLacI 参照 | 1158+1122+54 | 别构 shift 0.72 | 09-11 | anclaci_model/ |
| 真实进化图 GraphWalks | 6840 节点/9068 进化边 | 80.8%连通,862破greedy谷节点(z=2.25) | 09-14 | EVOLUTION_GRAPHWALKS.md |

## 六、下一步优先级
1. 翻转 bootstrap CI(不需 GPU/PI)
2. **GraphWalks 任务化**:精选强 valley 子集,greedy vs lookahead 可验证对比(=paper falsifiable 命题,不需 GPU/PI)
3. 7 位组训练数据补 + 重跑时序
4. origin 偏移正式 subsection
5. 胆酸 SMILES 对接
6. **补 lib0**(师兄湿实验)→ 真 round-0 基线,重估破 greedy 强度
7. 建控制面 + git init(待 PI)
