# 拆解与规划：多小分子突变进化景观 + 潜空间映射

日期 2026-09-08。用途：把"前面的方向"与两个压缩包的**具体内容**拆解清楚，作为向张老师/师兄征求具体建议的底稿。
承接 `BIO_DIRECTION_v1/2/3.md` 与 `METHODOLOGY.md`。

---

## 0. 一句话定位

这是一个 **TetR 家族变构转录因子（aTF）生物传感器的定向进化 + 建模**项目。
目标（张老师问的"算法突破 / 构建 landscape"）拆成三层可交付：

1. **经验景观**：从多轮筛选 NGS 反推"基因型 → 适应度"图（每个小分子各一张）。
2. **潜空间景观（潜空间映射）**：用 PLM 嵌入 + 高斯过程（GOLLuM 思路）学一个**共享潜空间**，分出"共同进化景观"与"各小分子特异景观"。
3. **可验证路径推理**：在真实景观上做 GraphWalks 式"可达路径 + 逐步状态追踪"的可验证过程监督（这是论文的 ML 新意）。

用途闭环：**你 AI 预测 ↔ 师兄进化数据**互相验证 → 进而对**新目标小分子**做序列预测。

---

## 1. 蛋白与实验体系（从代码反推，需向师兄确认）

- **WT 蛋白**（由 `DATA_PREP.ipynb` 里的 `wt_whole_dna` 翻译）= 204 aa：
  `MQKKLTRSQQKHLDIINAAKEEFIEFGFLAANMDRITSSAEVSKRTLYRHFESKEVLFESVLTIINDSVNESISYHFDPNKSTEEQLTEIAYKEIDVLYKTYGIALARTIVMEFLRQPEMAKTLIQNIYSIRAITQWFRSAIEAKRLKDADPKLMTDVYVSLFQGLFFWPQVMHLDLEPHGEELSQKIETLTTIFLQSYGVAE`
- 结构判读：N 端 ~aa9–45 是 **HTH DNA 结合结构域**，C 端是**配体（小分子）结合/效应结构域** → 典型 **TetR-family regulator (TFR)**。这正好解释"DNA 结合能力 + 小分子响应"双属性。
- **文库**：**从完全随机库起始**，基于 **17 个位点**引入突变的多突变库。
- **本包的选择逻辑（师兄已确认）**：选择压力 = **DNA 结合能力差的留不下来 → 正向筛选出有 DNA 结合能力的变体**。`338lib1→lib5` = 逐轮选择；**lib0 = 起始的完全随机库（= 富集基线 / 第 0 轮）**。
- **参考序列（已实证澄清，关键）**：野生型全长 612bp、蛋白 204aa（掐头去尾 202aa）。存在**两个 DNA 版本、蛋白完全相同**：
  - `work/ref/wt_reference_libbackbone.fasta`（旧密码子版）——**实测 3000/3000 条 reads 全部匹配这条**，师兄的库就是建在它上面，**bowtie2 比对必须用它**。
  - `work/ref/wt_reference_altcodon_userpasted.fasta`（另一密码子版，与旧版差 167/612 nt，同义）——蛋白等价但 reads 不匹配，**不能用于比对**。
  - 影响面：仅"DNA 比对参考"这一步；蛋白层面两版等价，下游景观建模不受影响。
- **待测建模的数据（师兄，双向筛选）**：同一蛋白为起点，以**相似但不同家族的小分子**为目标，**双向**筛选 = 保留 DNA 结合 + 优化对目标小分子的响应。

> 仍需确认：17 个位点清单是否都在配体口袋？**起始随机库（lib0/随机池）是否单独测序**（决定富集基线用实测还是理论随机分布）？双向筛选数据的结构（几个小分子、每个几轮、读出定量还是二元）？

---

## 2. 旧包 `vfa_training.zip` 全拆解（= 前面的数据处理链）

路径：`/data/lab/A0724/0815paper_gw/vfa_training.zip`（完整可用；已抽出源码到 `work/vfa_src/`）。
除 `.venv/.idea`（无关）外，真身是一条 **NGS → 变体计数 → 富集标签 → PLM 嵌入 → MLP 二分类** 的链：

| 阶段                   | 文件                      | 输入 → 输出                                                                                                                                        | 做了什么                                                                                                                                                                                                                                                        |
| ---------------------- | ------------------------- | --------------------------------------------------------------------------------------------------------------------------------------------------- | --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| **A. 比对**      | `align_ngs_seq.sh`      | 双端`fq.gz` + 参考 `fasta` → `sorted_reads.bam` / `pileup.txt`                                                                             | fastp(质控+merge) → bowtie2(`--very-sensitive-local`) → samtools                                                                                                                                                                                            |
| **B. 计数+标签** | `DATA_PREP.ipynb`       | BAM → 每轮 unique 序列计数 →`seq_fold_change.csv` / `protein_seq_fold_change.csv` / `idx2*.pkl` / `df_balanced.csv`,`df_imbalanced.csv` | ①按 CIGAR 抽比对序列、去重计数(count≥3)；②合并 lib0–5、算各轮 percentage；③相对 lib0 的 fold-change；④**enrichment 启发式**打二分类标签 `enriched`；⑤DNA→蛋白翻译(截 355nt、起点 `AAAGCGTGCTGA`)、蛋白级再算一遍；⑥造平衡训练集/不平衡测试集 |
| **C. 特征**      | `FEATURE_EXTRACTION.py` | 序列 CSV →`esm2_embedding.h5`(embeddings+labels)                                                                                                 | ESM2(`facebook/esm2_t33_650M`) 或 Nucleotide Transformer(`InstaDeepAI/...`) → **mean-pooling** 末层隐状态                                                                                                                                            |
| **D. 训练**      | `train.py`              | h5 →`outputs/models/mlp_classifier.pt` + PR/AUC 曲线                                                                                             | MLP(256-128-64) 二分类 enriched，BCE + cosine LR + early-stop，报 AUC/PR                                                                                                                                                                                        |

产物（包内已带）：`esm2_embedding.h5`(196MB)、`mlp_classifier.pt`、`enrichment*.png`、`training result.png`、`VFA阴性数据建模.pdf`。

**旧链的局限（= 为什么还没回答张老师的"景观"问题）**

1. **只做了"富集/不富集"二分类**，丢掉了 0–5 轮的**轨迹信息** → 不是景观、不是适应度。
2. **mean-pooling 嵌入抹掉了"17 个位点里哪几个突变了"** → 没有位点级/上位（epistasis）结构。
3. **单一小分子、单一选择轴** → 没有跨小分子的共享/特异比较。
4. **enrichment 标签是手调启发式**，不是校准的适应度值。
5. **17 位点散布在 612nt 上**，短读长可能无法在单条 read 里联合定相（phasing）→ 未处理。

---

## 3. 新包 `lib0_lib5negative_selection.zip` 全拆解（师兄最新数据）

路径：`./lib0_lib5negative_selection.zip`（**已传完，882MB，`unzip -t` 完整**）。
内容 = **原始双端测序 FASTQ**（= 阶段 A 的最上游输入，比旧包更原始）：

```
lib0_lib5negative_selection/
  338lib1/  338lib1_1.fq.gz  338lib1_2.fq.gz  MD5.txt
  338lib2/  338lib2_1.fq.gz  338lib2_2.fq.gz  MD5.txt
  338lib3/  ...  338lib4/ ...  338lib5/ ...
```

- 平台：Illumina NovaSeq（`@A01426…`），**双端、read 长 ~250nt**；R1 起点即 `A(N)AGCGTGCTGACC…`，正好落在旧 pipeline 的锚点 `AAAGCGTGCTGA` 上 → 覆盖诱变区，格式与旧链一致。
- 抽样可见单条 read 已带多处突变（如 S→F、L→M 等），符合"多突变随机库"。

> **⚠️ 关键澄清（已实证）：`338lib1–lib5` 是相对于"起始完全随机库(lib0)"的逐轮正向筛选。** 参考序列有两版同义 DNA、蛋白完全相同——**实测 3000/3000 条 reads 全部匹配旧密码子版**（`work/ref/wt_reference_libbackbone.fasta`），比对必须用它；用户新贴的另一密码子版蛋白等价但 reads 不匹配，**不能用于比对**。缺 lib0=起始随机池，须确认师兄是否单独测序（否则富集基线用理论随机分布）。

**待办（可立即做）**：`unzip` 到 `work/` → 每个 lib 跑 `align_ngs_seq.sh`（先确认 fastp/bowtie2/samtools 可用、建 WT 参考 fasta）→ 逐轮基因型计数 → 接第 7 节 pipeline。

---

## 4. "前面的方向"拆解（BIO_DIRECTION v1–3 + METHODOLOGY）

核心思想：**把进化重述为"受上位效应约束的可达路径推理"，而非回归**，并用 GraphWalks 式**可验证过程监督**训练学生模型。要点：

- 三档湿实验（档0 零湿实验 / 档1 一块板确认 / 档2 小 campaign），默认**档1**。
- 关键技巧：用"已测全景观"当**模拟湿实验 oracle（reveal-on-query）**，把主动探索也搬进 in-silico，指标 = 达标所需查询数 vs EVOLVEpro。
- Baseline：PLM zero-shot / EVOLVEpro / MULTI-evolve / 贪心爬山。
- 可证伪预测：**高 ruggedness / 高阶上位景观上，序贯路径规划显著胜过回归/贪心**。

> 这套框架此前是围绕**公开 DMS（GB1 等）**写的。现在我们手上有**真实的、多轮、多小分子**的自有数据 → 可以把它落到真实项目上，比公开数据更有故事。

---

## 5. 两篇新文章如何接进来

- **Shen 2021 — "Reconstruction of evolving gene variants and fitness from short sequencing reads"**（23 页）
  → 升级**阶段 B**。从多轮短读长**重建全长变体单倍型 + 逐轮适应度**，解决"17 位点跨 612nt 无法联合定相"和"手调 fold-change 不是适应度"两个硬伤。**产物 = 景观节点的真实 fitness 值 f(v) 与轨迹**。
- **Ranković 2026 — GOLLuM: Large language models as uncertainty-calibrated optimizers**（Nat Mach Intell, 15 页）
  → 提供**潜空间映射 + AI 预测再验证**的方法骨架。LLM/PLM 嵌入 + **高斯过程头**、贝叶斯目标训练、**不确定性校准**；嵌入空间被重塑得"相似结果的实验聚类"（= 学到的潜空间景观），比传统 BO **少 40% 实验**。天然支持**多任务 GP** → 共享 vs 小分子特异景观，且后验均值+方差可直接**主动提序列**。

---

## 6. 核心提案：三种"景观"构建（由易到难，回答张老师）

**(A) 经验景观（数据驱动，每个小分子一张）**
对每个目标小分子 m：节点 = 17 位点上的氨基酸组合（基因型），值 = Shen 重建的**逐轮 log-enrichment 适应度**；Hamming-1 连边。输出：局部最优、可达路径比例、ruggedness 指标；UMAP/力导向可视化。

**(B) 潜空间景观（共享 + 特异）= "潜空间映射"**
多任务/多输出 GP 建在 ESM2 嵌入上（GOLLuM 思路）：**一个共享潜空间 + 每个小分子一个 GP 头**。

- **共享成分** = 对所有小分子都有利（主要维持 DNA 结合）的突变 → **共同进化景观**；
- **分子特异残差** = 小分子响应特异的突变 → **各自的进化景观**。
- 不确定性校准 → 直接做**主动序列提名**。

**(C) 可验证路径推理（GraphWalks 迁移）**
在 (A)/(B) 之上实例化"WT→目标多突变体的可达路径 + 逐步状态追踪"任务，程序 checker 逐跳查表校验 → 可验证过程监督（论文 ML 新意）。

**交叉验证设计（你 AI 预测 ↔ 师兄进化数据）**

- 留一验证：留出某小分子的**后期轮次**或**整个小分子**，用其余训练，预测富集突变/路径，与师兄真实进化序列比 **precision@k / 路径召回 / 突变共现**。
- **零样本迁移到新目标小分子**：用"共享+特异"模型排候选突变 → 提名 top-k 序列 → 一块板（档1）确认。

---

## 7. 数据处理落地步骤（"前面的数据处理先做起来"——具体 TODO）

先在**旧包**上把链跑顺、代码化，等新包传完直接切数据：

1. **把 notebook 代码化、参数化**：`count_variants.py`（BAM→逐轮基因型计数矩阵）+ `build_landscape.py`（计数→17位点突变表→适应度）。
2. **恢复 17 位点结构表示**（关键缺口）：每条蛋白变体对 WT(204aa) 做 diff → **位点级突变表**（哪几个位点、变成什么氨基酸）；这是做上位/景观的基础。
3. **用逐轮 log-enrichment 适应度替换二分类**（Enrich2 / Shen 思路）。
4. **环境核查**：`fastp / bowtie2 / samtools` 是否可用；构建参考 fasta（WT 序列已知）；python 用 uv。
5. **确认嵌入来源与维度**：`esm2_embedding.h5` 到底是 DNA(NT) 还是蛋白(ESM2)、多少维（`train.py` 里 `input_dim` 注释与实际需核对）。
6. **新包传完**：`unzip -t` 验证 → 判定 fq/bam/counts → 跑上面链 → 产出「基因型×轮次计数矩阵 + 17位点突变表 + 适应度」三件套。
   交付：`LANDSCAPE_PLAN.md`（本文件）+ 骨架脚本。

---

## 8. 需向张老师 / 师兄确认的问题（你征求建议时直接用）

1. **蛋白身份**：确认是否 TetR-family aTF？WT 是否就是这 204aa？**17 个位点清单**（是否配体口袋残基）？
2. **负筛逻辑**：lib0 是否 naive 起始？负筛去掉的是"丢失 DNA 结合"还是"组成型结合"？每轮压力如何递增？
3. **双向筛选数据结构**：几个目标小分子？每个几轮？读出是 FACS / 存活 / 报告基因？**定量还是二元**？双向的两个选择维度怎么施加？
4. **测序**：单条 read（或 merge 后 ~355nt 窗口）能否覆盖全部 17 位点？若不能 → 需要 Shen 式重建或长读长/条形码。
5. **目标 venue**：决定要不要做湿实验档1（一块板前瞻确认）。
6. **小分子结构**：能否拿到各目标小分子的结构/SMILES？（做 molecule-aware 特征、共享景观时可并入分子指纹，支撑"对新小分子预测序列"。）

---

## 附：文献夹 `zsyPDF/`（11 篇，可选背景）

01 de novo 设计(Nature'26)、02 ESM world model、03 MULTI-evolve(Science'26)、04 EVOLVEpro(Science'25)、
05 ESM 压缩、06 EnzymeCAGE、07 进化尺度酶学、ATGCAGAAAAAATTGACTCGTTCTCAGCAGAAACACTTGGACATCATCAACGCTGCTAAAGAAGAATTCATCGAATTCGGTTTCTTGGCTGCTAACATGGACCGTATCACTTCTTCTGCTGAAGTATCTAAACGTACTTTGTACCGTCACTTCGAATCTAAAGAAGTATTGTTCGAATCTGTATTGACTATCATCAACGACTCTGTAAACGAATCTATCTCTTACCACTTCGACCCGAACAAATCTACTGAAGAACAGTTGACTGAAATCGCTTACAAAGAAATCGACGTATTGTACAAAACTTACGGTATCGCTTTGGCTCGTACTATCGTAATGGAATTCTTGCGTCAGCCGGAAATGGCTAAAACTTTGATCCAGAACATCTACTCTATCCGTGCTATCACTCAGTGGTTCCGTTCTGCTATCGAAGCTAAACGTTTGAAAGACGCTGACCCGAAATTGATGACTGACGTATACGTATCTTTGTTCCAGGGTTTGTTCTTCTGGCCGCAGGTAATGCACTTGGACTTGGAACCGCACGGTGAAGAATTGTCTCAGAAAATCGAAACTTTGACTACTATCTTCTTGCAGTCTTACGGTGTAGCTGAATAA08 HMD-AMP、09 NovaIscB、10 trans-AT PKS、11 630 机器人平台方案。
其中 **03 MULTI-evolve / 04 EVOLVEpro** 是主 baseline。
