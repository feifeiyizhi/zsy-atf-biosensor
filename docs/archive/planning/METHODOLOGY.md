# 方法论 Protocol：每一步具体做什么（重点写清生物）

日期 2026-08-16。承接 v2/v3。本文件 = 可直接照做的 step-by-step，分 in-silico（零湿实验）与湿实验（一块板）两大块。

---

## 术语先对齐（把 ML 词和生物词钉死）

| 概念 | ML/图 | 生物 |
|---|---|---|
| 节点 node | 图顶点 | 一个**基因型**：在选定 L 个位点上的一组特定氨基酸组合（如 GB1 的 V39/D40/G41/V54 四位点某组合） |
| 边 edge | 有向边 | 一次**单点突变**（在某位点换一个氨基酸，Hamming-1 邻接） |
| 节点值 f(v) | 标量 | 该变体的**实测活性/适应度**（DMS 读数：结合、催化、荧光、存活等） |
| 状态 state | visited 集合 + frontier | **当前基因型 + 迄今获得的突变集合**（因为 epistasis，下一步好坏取决于它） |
| walk / 路径 | BFS/遍历 | 一条**进化轨迹** WT→…→目标多突变体 |
| accessible path | 可达路径 | **每一步单点突变、且活性单调不降**的轨迹（Weinreich 定义：强选择弱突变下可走） |
| checker | 程序重算 | **已测景观直接查表校验**（in-silico）／**湿实验测量**（prospective） |

---

# 第一部分：IN-SILICO（零湿实验，是全部科学主张的主体）

## S1. 建"可验证适应度图"（生物数据 → 图）

**数据源（组合完备 / 近完备的真实 DMS，全部公开、已测全）**：
- **GB1 四位点**（Wu et al. 2016, eLife；Olson 2014）：蛋白 G 的 B1 域，位点 V39/D40/G41/V54，20⁴=160,000 变体测了 ~149k，适应度 = mRNA display 测 IgG 结合的 log-enrichment。**首选**，最经典、最 rugged。
- **TEM-1 β-lactamase 五突变**（Weinreich et al. 2006, Science）：5 个突变、2⁵=32 基因型全测，适应度 = 头孢噻肟 MIC。**accessible-path 的原始铁证集**，小而干净，适合做机制图。
- **avGFP 局部景观**（Sarkisyan et al. 2016, Nature）：~5 万基因型，荧光。
- **ProteinGym 多突变 split**（Notin et al. 2023）：200+ DMS 标准 benchmark，做 baseline 可比性用。

**建图脚本要做的事**（纯确定性、可复现）：
1. 读 DMS 表 → 每行一个基因型 + 实测活性，作为节点 v 和 f(v)。
2. 建边：对任意两节点，若 Hamming 距离=1（差一个单点突变）则连有向边（WT→更远方向）。
3. 归一化 f：以 WT 为 0 基线，或分位标准化（消除 assay 尺度）。
4. 标注：全局最优节点、各局部最优（邻居都不比它高的节点）、ruggedness 指标（局部最优个数、fraction of accessible paths）。
5. **切分防泄漏**：训练/测试用不相交的目标节点或不相交子景观；测试含"训练未见突变阶数/未见位点组合"做外推。

**产物**：与你 GraphWalks `synth_engine` 逐字节同格式的图对象（节点、边、gold 可达路径），只是语义是真实蛋白。

## S2. 定义任务 + 程序 checker（可 100% 精确验证）

**任务 T（序贯规划，不是回归）**：给定图，输出一条从 WT 到目标（全局最优或指定高活性多突变体）的**可达路径** π，要求：
- 每步是单点突变（合法边）；
- 每步 f 单调不降（accessible）；
- 终点 = 目标。

**checker（查表即验，免费）**：
- 逐跳校验：每一跳的 f 是否 ≥ 前一跳（对着 S1 的实测值查表）；
- 终点是否命中目标；
- 输出三档标注（对齐 GraphWalks）：`full-correct`（每跳都可达且到终点）/ `answer-correct-process-wrong`（终点对但中间某跳降了=蒙的）/ `fail`。

**评测指标**：
- **path-accuracy**：路径整体可达率；
- **optimum reachability**：从 WT 能否到全局最优 vs 卡局部最优（rugged 景观的核心指标）；
- **process-acc vs answer-acc**：定位失败是"中途状态失效"还是"终点误判"。

## S3. 合成"可验证状态追踪轨迹"（生物版 synth_bfs_trace）

对每张图，用图算法（BFS/DP over genotype graph）枚举 ground-truth 可达路径，然后**程序合成逐步 CoT**，每步写：
```
Current genotype: {WT+已获得突变集合}, measured fitness = f_k
Candidate next single mutations: {m1, m2, ...}
  - +m1 → fitness f(m1|background); Δ=+/−  (epistasis: 依赖当前背景)
  - +m2 → ...
Accessible move (Δ≥0, 不降): pick m*  → new genotype, fitness f_{k+1}
Updated acquired set: {..., m*}
...
Final path: WT → ... → target (每步活性不降)
```
- 每一行的 fitness 数字都来自 S1 实测景观 → **逐行程序可验证**（对齐你 GraphWalks "每层 frontier 可验证"）。
- 保留推理（不能只留骨架，复现 GraphWalks 的"丢 reasoning 负迁移"）。

**产物**：program-verified SFT 轨迹（100% 通过校验），用于训学生模型学会 epistasis-aware 序贯规划。

## S4. 训练学生模型（沿用你现成 pipeline）

- 输入 = 长上下文（图的全部边 + 已知变体读数）；输出 = 可达路径 + 逐步状态追踪 CoT。
- 监督粒度阶梯消融（论文灵魂，同 GraphWalks）：
  - B1 答案级（只监督终点多突变体）
  - B2 自由 CoT（教师原始推理，不结构化不验证）
  - B3 结构化状态追踪，但含 process-wrong 样本
  - **Ours** = 可验证 full-correct 状态追踪 + 由浅入深课程（先 2 突变路径→再 5+ 突变）
- 复用 GraphWalks 训练 recipe（长文 SFT），只换数据。

## S5. Baseline 头对头（全 in-silico）

在**同一任务（找可达路径 / 到全局最优）**上对比：
- **PLM zero-shot**：ESM log-likelihood ratio 排突变（无路径、无状态）。
- **EVOLVEpro**：PLM embedding + 随机森林回归 + 主动学习，贪心挑变体。
- **MULTI-evolve**：PLM 集成 + FCNN epistasis 回归，外推多突变体。
- **贪心 hill-climbing**：在预测 fitness 上爬山（必陷局部最优）。

**可证伪预测**：低 ruggedness 景观上大家接近；**高 ruggedness / 高阶 epistasis 景观上，序贯规划显著胜出**（回归/贪心卡局部最优）。用 GB1 vs β-lactamase 的 ruggedness 分层验证。

## S6. 把 Layer 3（主动探索）也搬进 in-silico —— reveal-on-query 模拟 oracle（零湿实验）

- 取组合完备景观（GB1），**假装是隐图**：模型初始只见 WT 邻域少数点。
- 模型每"查询"一个变体，就从**已测真实数据**返回其活性（模拟一次湿实验，但不花钱）。
- 模型据此更新状态轨迹、重规划下一批查询。
- **指标 = 达到目标活性所需查询数**（=等效湿实验次数），对比 EVOLVEpro 的主动学习查询数。
- **预测：路径推理更省查询** → 同时验证了算法 + 量化"我们更省湿实验"。

**至此，档 0 完成：一篇完整 ML-for-bio 论文的全部主张，零湿实验、真实生物数据、精确可验证。**

---

# 第二部分：湿实验（唯一一块板，前瞻性确认）

> 目的只有一个：证明"模型规划的路径终点，在真实分子上确实是更好的新蛋白，且赢过 baseline 的 pick"。单轮、无迭代。

## W0. 选蛋白（决定成本≈0 的关键）

**铁律：用实验室已有表达体系 + 已跑通 assay + 630 平台已验证的那个蛋白。** 不为论文新建体系。
- 候选（以实验室现成为准）：**APEX**（工程过氧化物酶，Amplex Red/荧光比色 assay，最简单）、某 Cas 类（切割/编辑读数）、某工程酶（偶联法测 kcat）。
- 从 MULTI-evolve 看 APEX 是他们验证过的靶（>100 倍活性提升），assay 成熟——如实验室有，优先。

## W1. In-silico 出候选（先在电脑上）

1. 用该蛋白的**现有数据**（已有的单/双突变读数，或 PLM zero-shot 先验）建局部图。
2. 我们的模型规划可达路径 → 输出**路径终点的 top-k 多突变体**（k≈8–16）。
3. 同预算下 EVOLVEpro / MULTI-evolve 各出 top-k。
4. 合并去重 → 一份 ≤24（或 ≤48）的变体清单，含 **WT 阳性 + 已知失活突变阴性**对照。

## W2. 造变体（分子克隆，一次性）

- 多突变体用**基因合成**（gBlocks/全长合成，多突变一次到位，省去多轮定点突变）或 **MULTI-assembly**（MULTI-evolve 的多位点组装，一次装 9 突变、~70% 效率）。
- 克隆进实验室现成表达载体（E. coli / 无细胞表达皆可）。

## W3. 表达纯化（96 孔平行 / 630 平台上样）

- 小量平行表达（96 孔或对应板位），粗裂解物或快速纯化（His-tag）。
- **这一步正好是 630 机器人方案里的"孔板取放→上料→启动"**：一块板一次上样、一次启动检测。用"基础版（单台精灵G2+固定物料站）"脚本即可，无需冷链/双机。

## W4. 测活性（一次检测，一块板）

- 用实验室现成 assay 一次读全板：
  - APEX → Amplex Red / luminol，酶标仪读荧光/发光；
  - 酶 → 底物偶联法测初速率 → kcat/KM；
  - Cas → 切割/编辑效率读数。
- 每变体 3 复孔即可（板内足够）。

## W5. 读结论（无迭代）

- 主结论：**我们路径终点的活性 > WT，且 ≥ baseline 的最佳 pick**。
- 样本效率结论：**我们只用 N 次测量拿到该增益；EVOLVEpro 达同等增益的模拟 oracle 曲线需 M≫N 次**（W 的真实点 + S6 的模拟曲线交叉印证）。

**总湿实验量：1 蛋白 × 1 轮 × ≤48 变体 × 1 块板 × 1 次检测。** 落在 630 基础版脚本内。

---

## 论文图对应（湿实验只占 1 张，可放附录）

| 图 | 内容 | 湿实验？ |
|---|---|---|
| Fig 2 | GraphWalks 合成图受控证明（你已有） | 否 |
| Fig 3 | 真实 DMS 图上 路径推理 > 回归 baseline（S5） | 否 |
| Fig 4 | 监督粒度阶梯 + 机制预测（process-acc 先崩、ruggedness 分界） | 否 |
| Fig 5 | 模拟 oracle 主动探索查询效率 > EVOLVEpro（S6） | 否 |
| Fig 6 | 单轮湿实验确认：造出真实更好蛋白 + 赢 baseline（W） | **是，仅此一张** |

---

## 立即可做（今天就能起，全零湿实验）
1. **S1**：下 GB1（Wu 2016）DMS 表 → 建可验证适应度图 + ruggedness 标注。
2. **S2+S3**：写 checker + 迁移 synth_engine 成"可达路径 trace 合成器"。
3. **S5**：跑 PLM zero-shot / EVOLVEpro / 贪心 baseline 打对照。
4. **S6**：实现 reveal-on-query 模拟 oracle，出查询效率曲线。
5. 全部成立 → 拿 Fig3-5 去找张老师谈那"一块板"的 W。

## 待办 / 风险
- [ ] 跟张老师确认 W0 的蛋白（现成 assay + 630 已跑通 = marginal 成本≈0）。
- [ ] 校正张老师方向（本机搜索引擎不可靠；给我他主页/Scholar URL 我用 curl 直取核实）。
- 风险：真实景观 ruggedness 不足 → 序贯优势不显著。**S1+S5 就是最早 go/no-go，在任何湿实验前判定。**
