# 深度方案 v2：把 GraphWalks 可验证过程监督做成"蛋白进化的序贯规划推理"

日期 2026-08-16。作者：陈晨澍。导师：清华药学院 张数一（Shuyi Zhang，方向由 zsyPDF 反推，待校正）。
> v1 的问题：停在"进化像 GraphWalk"的类比层。v2 把它压成一个**精确的计算问题**、一个**可精确验证的 benchmark**、以及**可证伪的机制预测**。

---

## 1. 精确的问题重述（不是类比，是同一个计算问题）

**当前 SOTA 在解什么**（从原文核实）：
- **EVOLVEpro (Science'25)**：PLM embedding → 顶层随机森林回归，iterative active learning，每轮 ~10 实验点，预测**单个变体**的活性排序。
- **MULTI-evolve (Science'26)**：PLM 集成挑单点增益突变 + FCNN 学 epistasis，用**双突变**外推到 8–12 突变体，**单轮**。原文明确："proteins exhibit complex epistasis, where mutational effects combine nonlinearly to create rugged fitness landscapes"。
- **Evolutionary-scale enzymology (Science'25)**：HT-MEK 测 193 个 ADK 直系同源的 kcat/KM，明确解析"**sequence-catalysis landscape 的 topology 与 navigability**：rugged、3 个全局邻域、**path-dependent、stepwise 突变仍 navigable**"。Fig.1A 直接画了 sequence space 上通向不同 optima 的进化路径。

**它们的共同盲区**：全部是**回归/预测单节点的适应度**，然后贪心/排序挑变体。**没有一个在显式地规划"一条路径"、也没有一个在追踪"状态"。**

**我们要解的问题（精确版）**：
> 给定适应度图 G=(V,E)，V=序列变体，E=单点突变（Hamming-1 邻接），每个节点有隐活性 f(v)。求从野生型 WT 到目标高活性多突变体 v* 的一条**可达路径（accessible path）** π=WT→…→v*，使得**每一步 f 单调不降**（或满足给定的路径级约束）。

这为什么是**图上的状态依赖 reachability 问题**、且恰是 GraphWalks 的内核：
- **sign epistasis** 使某突变的正负效应**取决于当前背景已积累的突变集合**——即"边是否可走"依赖 visited 状态。这与 BFS 里 visited 集合决定下一层 frontier 是**同一数学结构**（状态依赖的图遍历）。
- 一个突变的效应不是节点固有属性，而是**(突变, 当前状态) 的函数** → 回归节点值的模型（EVOLVEpro/MULTI-evolve）结构上无法表达路径可达性，会卡在局部最优。
- 这正是你 GraphWalks 论文的论点在生物上的落地："没有显式状态追踪的多跳推理 → 串行深度预算耗尽 → 解不了。"

**理论铁钉**：
- **Weinreich et al. 2006 (Science)**：β-内酰胺酶 5 突变、5!=120 条突变顺序，只有 **18 条**是"每步活性单调上升"的可达路径 → epistasis 把 ~85% 的路径堵死。**可达路径搜索是真问题，不是回归能替代的。**
- **NK / rugged landscape 理论**：ruggedness ⇒ 贪心/单步预测陷局部最优；找可达路径需要序贯前瞻（look-ahead）+ 状态记忆。
- **GraphWalks 论文三定理**（CoT 解串行问题、图连通性需 log-depth、CoT=serial-circuit）：给"为什么必须显式序贯状态追踪推理、而非并行回归"提供计算复杂度层面的理由。

---

## 2. 真正的桥梁：组合完备 DMS = 一张已知且可精确验证的适应度图

这是 v1 完全缺失、也是把工作从"类比"变"硬核"的**关键中间层**。

**洞察**：组合完备（combinatorially-complete）的深度突变扫描数据集 = **每个节点活性都测过的完整图**。于是"在这张测好的图上找可达路径/全局最优"就变成**程序可 100% 精确验证**的任务——和 GraphWalks 的程序 checker 逐字对应，只是语义换成真实蛋白。

**现成的完整/近完整景观（可直接建图，无需任何湿实验）**：
| 数据集 | 规模 | 为何是完整图 |
|---|---|---|
| **GB1**（Wu et al. 2016） | 4 位点，20⁴≈160k 变体，测了 ~149k | 组合近完整，经典 rugged landscape |
| **β-lactamase**（Weinreich 2006） | 5 突变，2⁵=32 节点全测 | accessible-path 的原始证据集 |
| **Entacmaea/avGFP**（Poelwijk, Sarkisyan） | 局部组合景观 | 高阶 epistasis 完整测量 |
| **ProteinGym** DMS + multi-mutant split（Notin et al.） | 200+ DMS，含多突变 | 标准化 benchmark，可比性强 |

**把 GraphWalks 那套原封搬过来**：
- 图已知 → 逐步路径每一跳"活性是否不降 / 是否走向 v*"可对着测好的景观**逐跳程序校验** → 生成 100% 正确的**状态追踪 gold 轨迹**（= 你 synth_engine 里 `synth_bfs_trace` 的生物版）。
- 训练：可验证过程监督教模型做 epistasis-aware 的序贯规划（维护"已获得突变=状态"，逐跳推理下一个可达突变，避开 sign-epistasis 陷阱）。
- **评测（可证伪，全 in-silico）**：
  - **path-accuracy**：模型给的路径是否真可达（对景观校验）。
  - **optimum reachability**：从 WT 出发能否到全局最优 vs 卡局部最优。
  - **对比 baseline**：EVOLVEpro / MULTI-evolve / 贪心 / 纯 PLM zero-shot——预测**这些回归方法在高 ruggedness / 高阶 epistasis 区显著劣于序贯规划**。
  - **process-acc 先于 answer-acc 崩**（沿用 GraphWalks 的探针）：把"陷局部最优"定位成"状态追踪在某突变阶数 d* 失效"。

**这一层的价值**：审稿人最认——真实生物数据、精确可验证、可复现、无需湿实验就能打赢现有回归范式。它是 Paper-2 的地基。

---

## 3. 三层递进（完整故事线）

1. **GraphWalks（合成图，免费验证）** = 纯机制证明：可验证过程监督能教会状态追踪推理。（Paper 1，已在推进）
2. **组合完备 DMS（已知真实图，可精确验证）** = 迁到真实生物语义，正面打赢 EVOLVEpro/MULTI-evolve 的回归范式，**全 in-silico**。（Paper 2 地基，新增的关键层）
3. **630 平台前瞻验证（隐图，湿实验做 oracle）** = 图未知：模型在**部分观测**下推理，**主动查询湿实验揭示新边**（= 隐图上的主动图探索 / active graph discovery）。跑真实 campaign，比 EVOLVEpro 更少轮次 / 更高终活性。（Paper 2 旗舰湿实验层）

**第 3 层的技术形式**（把 EVOLVEpro 的 active learning 抬高一层）：
- EVOLVEpro：每轮回归器挑 batch 变体测 → 更新回归器（无状态、无路径记忆）。
- 我们：agent 维护**可审计的进化状态轨迹**（已探明的 epistasis 结构 = 已揭示的子图）→ 推理"下一跳最有信息量 / 最可能可达的多突变体"→ 湿实验揭示新边 → 用实测逐跳校验并精化轨迹。**湿实验充当 GraphWalks 里那个程序 checker**（昂贵、噪声、批量），但过程监督结构不变。

---

## 4. 相对张老师给的文章，novelty 的精确表述（审稿人必问）

| 维度 | EVOLVEpro / MULTI-evolve | 我们 |
|---|---|---|
| 建模对象 | 单节点适应度回归 | WT→v* 的**可达路径**（序贯决策） |
| epistasis | FCNN 黑箱外推 | **显式状态**（visited 突变集）驱动可达性推理 |
| 推理 | 无（一次预测/排序） | **多跳状态追踪推理 + look-ahead** |
| 验证 | 回归精度指标 | **逐跳程序/实验可验证的过程监督** |
| 长上下文 | 不用 | 读全部历史变体+读数（你的强项） |
| 陷局部最优 | 会（贪心/排序） | 序贯规划显式规避 sign-epistasis 陷阱 |

**一句话 novelty**：把定向进化从"预测哪个突变好"（回归）升级为"在受 epistasis 约束的适应度图上规划一条可验证的可达路径"（序贯状态追踪推理）——这是范式升级，不是又一个更准的预测器。

---

## 5. 可证伪的机制预测（把故事变成可做的实验）

1. **序贯 > 回归的分界**：在低 ruggedness / 低阶 epistasis 景观上序贯规划 ≈ 回归；在高 ruggedness / 高阶 epistasis 景观上序贯规划**显著胜出**（回归卡局部最优）。→ 用 GB1/β-lactamase 的 ruggedness 分层验证。
2. **状态追踪先崩**：模型失败时，process-acc（逐跳可达性对不对）先于 answer-acc（终点对不对）崩塌 → 失败=状态丢失而非终点误判。
3. **可验证过程监督 > 自由 CoT > 答案级**：沿用 GraphWalks 的监督粒度阶梯消融，在 DMS 图上复现。
4. **长度/阶数外推**：训 ≤k 突变，测 >k 突变的可达路径（对应 GraphWalks 的 depth 外推）。
5. **前瞻性**（湿实验）：同实验预算下，路径规划 agent 的终活性 / 命中率 > EVOLVEpro。

---

## 6. 一篇 or 两篇（更新判断）

有了第 2 层（DMS 可验证图），**合成一篇的可行性大幅提高**：
- **合一篇（冲 Nature Methods / Science 系）**：GraphWalks=图2受控证明；DMS 可验证图=图3-4主 in-silico 结果打赢 baseline；630 湿实验=图5前瞻验证。故事完整、影响力最大，代价是等湿实验。
- **拆两篇**：Paper1（GraphWalks，ML 顶会，快）先落一作；Paper2（DMS+湿实验，跨学科顶刊）跟上。降风险。

**建议**：先按"Paper1 独立 + Paper2 以 DMS 层为地基先行做 in-silico"推进；DMS 层 in-silico 结果出来（能打赢 EVOLVEpro/MULTI-evolve）后，再决定 Paper2 是否等湿实验合成一篇大的。**第 2 层无论如何都要先做——它既是 Paper2 地基，也是最快能验证整个 thesis 的一步，且不花湿实验钱。**

---

## 7. 立即可做（不依赖湿实验、不依赖联网）
1. **建第一张可验证适应度图**：下 GB1（Wu 2016）或 ProteinGym 多突变 split，写 checker（节点=变体、边=Hamming-1、gold 路径=测好景观上单调可达路径），做成 GraphWalks 同构格式。
2. **打 baseline**：EVOLVEpro / MULTI-evolve / 贪心 / PLM zero-shot 在"找可达路径 / 到全局最优"任务上的成绩，作为对照。
3. **迁移 synth_engine**：把 BFS trace 合成器改成"epistasis-aware 可达路径 trace"合成器，逐跳对景观校验。
4. **神经元/细胞状态**：放 discussion 当"同一抽象的更广外延"，不做主线（与本实验室蛋白资产不对口）。

## 8. 风险 / 待办
- [ ] 校正张老师真实发表与在研方向（本机搜索引擎不可靠，未在线核实成功）。
- [ ] 查竞品避免撞车：MLDE/ftMLDE（Arnold/Yang）、self-driving lab / agentic science、ProteinGym leaderboard 上的 supervised 方法。差异化=**路径级可验证过程监督 + 序贯状态追踪**，不是更准的 fitness 回归。
- [ ] 确认组合完备 DMS 里 ruggedness/高阶 epistasis 足够强，能撑起"序贯>回归"的分界预测。
- 风险：若在真实 DMS 图上序贯规划相对回归优势不显著（景观没那么 rugged），则 thesis 削弱——**第 7 步就是最早的 go/no-go 检验点。**
