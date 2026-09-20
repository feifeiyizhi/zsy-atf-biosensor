# 生物方向调研 + 故事方案（把 GraphWalks 可验证状态追踪嫁接到蛋白进化）

日期 2026-08-15。作者视角：陈晨澍。导师：清华药学院 张数一（Zhang Shuyi）。
> 注：本机联网被封，张老师确切发表未在线核实；其方向由他给的 11 篇 zsyPDF + 630 机器人方案反推，请校正。

---

## 0. 张老师方向反推（从 zsyPDF 的选文口味）

选文高度收敛在一个主题群：**PLM/ML 引导的蛋白定向进化 + 进化引导的酶/工具工程 + 高通量湿实验自动化**。

| 文件 | 主题 | 给我们的关键词 |
|---|---|---|
| 04 EVOLVEpro (Science'25) | PLM embedding + few-shot 主动学习回归，**迭代**优化蛋白，每轮仅需 ~10 个实验点 | 迭代主动学习、模型在环 |
| 03 MULTI-evolve (Science'26) | PLM 集成 + **神经网络建模上位效应(epistasis)** 外推多突变体，单轮造 9 突变 | **上位效应 = 状态依赖** |
| 07 Evolutionary-scale enzymology (Science'25) | HT-MEK 高通量测 kcat/KM，解析 **sequence-catalysis landscape 的拓扑与可导航性**：rugged、多个全局邻域、**path-dependent、stepwise 突变仍 navigable** | **景观=图，进化=图上行走** |
| 06 EnzymeCAGE (Nat Catalysis'26) | 几何基础模型做酶检索/反应去孤儿/**生物合成通路重建** | 反应网络=图，多跳 |
| 09 NovaIscB (Nat Biotech'25) | 进化引导设计 IscB 表观编辑器 | 进化引导工程终点 |
| 10 trans-AT PKS (Science'24) | 共进化分析找模块插入位点，工程化 PKS | 共进化=边约束 |
| 08 HMD-AMP (Nat Biomed Eng'26) | PLM 挖进化远端抗菌肽，湿实验验证 91→74 有效 | PLM + 湿实验闭环 |
| 01 de novo design (Nature'26)、02 ESM world model、05 ESM compression | 综述 + PLM 基础设施 | 底座 |
| 11 630 蛋白进化机器人方案 | **自动化湿实验平台**：孔板上料→启动，四套难度脚本 | **湿实验闭环的物理载体** |

**一句话画像**：张老师做的是"用 ML/PLM 在崎岖的蛋白适应度景观上做定向进化，并用高通量自动化平台闭环验证"。核心痛点全在 **epistasis（上位效应）导致的路径依赖**——一个突变好不好，取决于已经积累了哪些突变（=状态）。

---

## 1. 核心桥梁：蛋白定向进化**本身就是一次 GraphWalk**

| GraphWalks（我们已有） | 蛋白定向进化（张老师战场） |
|---|---|
| 节点 = hash 节点 ID | 节点 = 蛋白序列/变体 |
| 边 = 有向边 | 边 = 单点突变 |
| 一次 walk = BFS 逐层前沿 | 一次进化 = stepwise 突变堆叠轨迹 |
| **状态** = visited 集合 + 当前 frontier | **状态** = 当前序列 + 已发现的上位效应结构（epistasis context） |
| 答案 = depth-d 的节点集 | 目标 = 到达高活性峰 / 满足多目标的多突变体 |
| **验证器 = 程序重算（免费、精确）** | **验证器 = 湿实验（昂贵、有噪声、真实）** |
| 图**已知、可观测** | 图（适应度景观）**未知、隐变量，只能实验查询** |

**这不是比喻，是同构。** 经典理论 Weinreich 2006（"Darwinian evolution can follow only very few mutational paths"）早就证明：sign epistasis 会把绝大多数突变路径变成"不可达"，进化只能沿极少数 accessible path 走——**这正是一个图上"哪些路径可达"的串行搜索问题**，和 GraphWalks 的"串行计算深度预算 + 状态追踪"论点是一回事。

**两者唯一的、本质的差别**：GraphWalks 的图和答案是**已知且免费可验证**的（in-silico oracle）；蛋白进化的图是**隐的**，唯一的 oracle 是**湿实验**。

---

## 2. 一句话 thesis（漂亮的故事）

> **定向进化是在一张隐适应度图上的、受上位效应约束的串行行走；其失败不是"预测单点突变好坏"的回归问题，而是"在实验预算下做可验证的多跳状态追踪推理"的问题。** 我们先在 GraphWalks（图已知、验证免费）上证明可验证过程监督能教会模型这种状态追踪推理，再把同一能力迁移到蛋白进化（图隐、湿实验做验证器），用长上下文推理模型显式维护一条"可审计的 epistasis 状态轨迹"，闭环驱动高通量进化平台。

这条故事线同时接住了：
1. **GraphWalks 论文**（in-silico、可验证、快）= 机制证明 / 受控 testbed。
2. **张老师的湿实验平台 + PLM 进化专长** = 真实部署 + 顶刊级 wet-lab 验证。
3. **比 EVOLVEpro/MULTI-evolve 深一层**：它们用回归器/NN 预测单点或外推多突变（黑箱、无显式状态、无多跳推理）；我们把它抬成"长上下文推理 + 可验证过程监督 + 实验在环"的**推理范式**。

---

## 3. 具体选题（推荐 Option A）

### Option A（首选，最贴实验室、最能打 baseline）：Epistasis-aware directed evolution as verifiable graph-walk reasoning
- **输入**：起始蛋白 + 迄今所有变体的实验读数（长上下文——正好是我们的强项）。
- **模型**：一个长上下文推理 agent，显式维护"已探明的 epistasis 图/状态轨迹"（frontier = 当前可扩展的高价值突变组合），**逐跳推理**下一批要造哪些多突变体。
- **可验证过程监督**：每一步预测的效应 vs 湿实验实测效应 = 湿实验充当"逐层 checker"（对应 GraphWalks 的程序 checker）；实测回来的结果反过来监督/精化下一跳。
- **湿实验验证（630 平台）**：跑一个真实 campaign（如 APEX / Cas 类 / 某酶），比 EVOLVEpro、MULTI-evolve **更少轮次 / 更高终活性 / 更少实验点**。
- **GraphWalks 在本文的角色**：受控消融——证明"可验证状态追踪推理"能力是真的、可训练的（图已知时打满分），再证明它迁到隐图 + 实验 oracle 仍有效（sim→real）。

### Option B（次选）：生物合成通路 / 代谢网络多跳推理（接 EnzymeCAGE）
- 图 = 反应网络；任务 = 通路重建 / 酶去孤儿 / 逆合成多跳。湿实验 = 表达候选酶测预测反应。可行但更慢、baseline 不如 A 直接。

### Option C（只作为 generalization/future，不做旗舰）：神经元/细胞状态轨迹
- 图 = 细胞状态转移图（单细胞轨迹）。和本实验室蛋白资产不对口，湿实验闭环难。你提到的"神经元状态"适合放在 discussion 里当"同一抽象的更广外延"，别当主线。

**判断**：你直觉里"蛋白进化那个点"比"神经元"更对——对这个实验室它是最紧的、最可 wet-lab、最能正面超基线的。

---

## 4. 相对张老师给的文章，novelty 在哪（审稿人必问）

- EVOLVEpro / MULTI-evolve：**回归器 + embedding** 预测突变效应，把 epistasis 当黑箱 NN 外推；**不推理轨迹、不维护显式可审计状态、不做可验证多跳规划**。
- 我们：把进化当**图上行走的推理问题**——长上下文读全部历史变体+读数 → 维护**可验证的 epistasis 状态轨迹** → 推理下一跳 → 湿实验验证 → 验证结果做过程监督。**这是范式升级（回归/预测 → 可验证过程推理），不是又一个预测器。**
- 理论锚点：Weinreich 2006（accessible path 极少）、NK rugged landscape、以及 GraphWalks 论文的串行深度/状态追踪三定理——把"为什么必须显式状态追踪推理"钉在理论上。

---

## 5. 一篇还是两篇？哪个更有意义？

**建议：一个统一 thesis 下的两篇，按顺序做。**

1. **Paper 1 = GraphWalks 可验证过程监督**（已在推进）。ML/长文顶会（NeurIPS/ICLR/ACL）。快、低风险、确立方法与你的一作，作为机制证明。
2. **Paper 2 = 蛋白进化旗舰**（wet-lab 闭环）。跨学科顶刊（Nat Methods / Nat Biotech / Science 系，配实验室湿实验）。高影响、贴导师、把 Paper 1 当受控 testbed 引用。

**"更有意义吗？"** 对药学院 + 真实影响 + 顶刊而言，**Paper 2 更有意义、更贴导师**；但它更慢、更依赖湿实验成功。Paper 1 是给 Paper 2 **降风险的机制证明**，两者不是二选一，是"先证机制、再落真实"。

**能否合成一篇？** 能，且对顶刊更有冲击力的写法是：GraphWalks 作为"in-silico 受控证明"当图 2/消融，蛋白进化 + 湿实验验证当主结果。代价是周期长、必须等湿实验出数。**推荐先按两篇推进，Paper 2 视湿实验进度再决定是否把 GraphWalks 并进去。**

---

## 6. 待办 / 风险
- [ ] 校正张老师真实发表与在研（本机无网，未在线核实）。
- [ ] 定 Paper 2 的靶蛋白 + 现成 DMS/campaign 数据集（先纯 in-silico 用公开 DMS 打 baseline，再上 630 平台做 prospective 验证）。
- [ ] 查竞品：LLM-for-science agent（Coscientist/ChemCrow 类）、MLDE/ftMLDE 主动学习、agentic self-driving lab——我们的差异化=可验证状态追踪 + 实验在环过程监督，别撞车。
- 风险：湿实验周期与成败决定 Paper 2 时间线；先用公开 DMS 数据把"图上推理 > 回归 baseline"证到位，再 prospective。
