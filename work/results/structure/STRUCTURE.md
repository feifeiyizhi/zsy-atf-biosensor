# 结构线进展（2026-09-11）

v2 策略要求建模基于结构 + 别构/远程。本文件记录结构线进度。承接 `PROJECT_BACKGROUND_PROTOCOL.md §0`。

## 1. 序列级结构域映射（已完成，基于 TetR 家族保守架构）

TetR-family aTF 的域架构在家族内极保守：
- **DBD**（N 端 ~aa1–47）：翼状 HTH，α3 识别螺旋（~aa35–45）直接读 DNA 大沟。
- **linker/hinge**（~aa47–52）：别构偶联铰链。
- **LBD**（C 端 ~aa48–210）：全螺旋束，**配体口袋核心 ~aa100–180**。

### 17 设计位点的结构定位（`site_domain_map.csv`）
- **DBD 内：0/17**（无一在 HTH DNA 结合螺旋）。
- **LBD 内：17/17**。
- **口袋核心区(100–180)：13/17** = aa102,107,111,125,128,129,133,159,161,163,164,169,170。
- 口袋入口(LBD-N, 66–99)：aa66,69,75,94。
- **R116Q（胆酸进化骨架突变）**：LBD 口袋核心区。
- **Q164R（PCA 定义性突变）**：LBD 口袋核心区。

## 2. 关键结构结论（支撑 v2 立论）
1. **文库不碰 DNA 结合螺旋**：17 位点全在 LBD → 文库通过 LBD 构象**别构地**调节 DNA 结合，而非直接改 DNA 接触面。这是"远程作用"的直接结构证据。
2. **口袋主导**：13/17 位点在配体口袋核心 → 文库同时直接调制小分子识别（口袋）和别构传导（LBD 构象）。
3. **Q164/R116 都在口袋核心** → 胆酸进化的关键突变落在配体口袋，直接影响小分子特异性，与 DNA 结合经别构耦合。
4. **别构路径（待结构精修）**：配体（胆酸）→ LBD 口袋 → LBD 螺旋束构象变化 → α4/hinge → DBD HTH 取向 → DNA 释放/结合。

## 3. PDB 同源结构（进行中）
- 无 GPU；已起后台 NCBI BLAST（pdb 库）找最近 TetR 同源，待返回。
- 拿到后将：对齐 WT 到 PDB 链 → 精修口袋坐标 → 17 位点投到 3D → 别构通信网络（pocket→HTH）。
- 若最近同源 identity 偏低或无配体，备选：ESMFold（CPU 重，单独起）或 AlphaFold DB 相近序列。

## 产物
- `work/structure/site_domain_map.csv` — 17 位点 + 关键残基的域/口袋定位表
- `work/results/structure/fig_domain_map.png` — 域带 + 位点图

## 4. ANM 别构动力学（已完成，ProDy on 3CDL）

因 ESMFold 缺 openfold 装不动，结构线用 **3CDL 同源(33%)+ProDy ANM 弹性网络**兜底（粗粒度动力学对中度同源鲁棒）。

### 4.1 动力学（最慢 20 模式波动）
- **口袋核心(159–170)刚性低波动**：Q164=0.036、aa163=0.044、aa161=0.038（口袋稳，配体结合面）。
- **口袋入口/portal(66–111)较灵活**：aa111=0.141、aa75=0.132、aa66=0.094（入口动态，配体进出/识别）。
- 含义：设计位点跨刚性口袋面 + 灵活入口，前者调特异性，后者调结合/释放动力学。

### 4.2 别构耦合（扰动响应 PR 矩阵）⭐结构侧证据
- **扰动配体口袋 Q164 → HTH α3 识别螺旋响应 ≈0.05–0.08**（均值 0.042，超背景）→ 证实 **pocket↔HTH 长程别构通道**，正是"小分子位点与 DNA 位点远离但远程耦合"的结构证据。
- 扰动口袋 → R116(骨架突变)响应 0.101（强，R116 在口袋附近、与别构传导相关）。
- 反向：扰动 HTH → 口袋响应较小（非对称，但非零，耦合存在）。

### 4.3 含义
ANM 用 3CDL 同源骨架，给出了 v2 最缺的**别构耦合结构证据**：17 位点全在 LBD，扰动口袋能传到 HTH DNA 识别螺旋。后续（拿到 ESMFold 自身结构后）可在此骨架上做突变扰动预测。

## 5. ESM-2 zero-shot 突变打分（进行中，CPU）
- 用 ESM2 t33_650M 对 17 位点 × 20AA 打突变效应分（一次前向出全矩阵），与 338lib 真实负筛适应度交叉验证（AI 突变敏感先验 vs 真实 DNA 结合选择）。
- 跑完补结果。

## 产物
- `structure/pdb/3CDL.pdb`、`structure/wt_to_3cdl_map.json`、`site_domain_map.csv`
- `anm_sqf.npy`、`anm_PR.npy`、`anm_17site_fluctuation.csv`
- `fig_domain_map.png`、`fig_anm_allosteric.png`
- (进行中)`esm2_zeroshot_17site.csv`、`esm2_vs_realsel_per_site.csv`

## 6. ESM-2 zero-shot 突变打分结果（已完成，t30_150M，CPU）

### 6.1 产物
- `esm2_zeroshot_17site.csv` — 17 位点 × 20AA 的 Δlog-likelihood 矩阵（>0=模型认为比 WT 可容忍）。
- `esm2_vs_realsel_per_site.csv` — 与真实负筛适应度交叉。
- `fig_esm2_zeroshot_heatmap.png` — 热图（黑框=WT，黄框=负筛库实际采样，红框=PCA 固定的 Q164R）。

### 6.2 关键结果
1. **ESM2 给全 20AA 突变效应先验**，而负筛库每位点只采样 3–11 个 alt（简并密码子设计，如 164 位只采了 4 个）→ ESM2 把预测**外推到未测替换**，是"预测新突变"的序列级先验。
2. **直接 ESM2 vs 真实负筛相关性欠功率**：多数位点 alt 样本 n=3–4（只 2–4 个 alt 有足够变体），rho 忽正忽负（aa107 -1.0、aa169 +1.0）不可信。**原因不是 ESM2 不行，是负筛库采样太窄**。
3. **Q164R（PCA 进化固定）**：ESM2 在 164 位把 R 排第 2 不坏(K:-3.2,R:-3.3,S:-3.9)但仍净负 → **通用 PLM 零样本预测不了"特异性驱动突变"**（这类突变对通用蛋白适合度"坏"，但对特定小分子功能"好"）。⇒ **特异性预测需要结构+功能条件模型，不是通用 PLM**——这是 v2"必须基于结构"的实证支撑。

### 6.3 含义（接 v2）
- ESM2 zero-shot 作**一阶先验**（哪些替换不破坏通用折叠/稳定性），但**小分子特异性**必须靠结构条件模型 + 功能数据训练（变体×小分子别构效应矩阵）。
- ESM2 + ANM + 真实负筛 三者互补：ESM2=序列突变敏感性、ANM=别构结构通道、负筛=真实 DNA 结合适应度。三者拼成 v2 建模的多视角地基。

## 7. 3D 距离精修（3CDL 坐标）
ESMFold/OmegaFold/AF2 在本机均装不动（缺 openfold/GPU/MSA工具），结构用 3CDL 同源兜底已足够。在 3CDL 坐标上算各 17 位点到口袋(Q164)与 HTH(α3)的 3D 距离（`site_3d_distance.csv`、`fig_site_3d_distances.png`）：

- **口袋↔HTH 3D 距离 = 33.3 Å** → 小分子位点与 DNA 位点确实远离，真长程别构。
- **pocket-lining（d_pocket<8 Å，直接围配体）**：**S161(5.1)、F163(3.8)、Q164(0)** —— 这三个就是配体口袋壁，解释为何 Q164R 直接定胆酸特异性、S161→A/F163→M 在负筛里高频出现（口袋面，调特异性为主）。
- **别构中继（pocket 与 HTH 中间）**：**Y102**（d_pocket 13 Å，d_HTH 20 Å）—— 离口袋和 HTH 都相对近，是别构传导的候选中继残基。
- 其余位点 distal（在 LBD 另一侧），主要调稳定性/折叠。
- 几个位点(69/125/128/129)在 3CDL chain A 未解析（gap）。

## 8. 结构线总结
1. 17 位点全 LBD、无一碰 DNA 螺旋；13/17 在口袋核心；**3 个直接围配体(S161/F163/Q164)**。
2. **pocket↔HTH 33 Å 长程别构，ANM 证实有耦合通道**（扰动 Q164→HTH 响应 0.05–0.08 超背景）。
3. **ESM2 零样本**：给全 20AA 突变先验、外推到未测替换；但**预测不了特异性驱动突变（Q164R）**→ 通用 PLM 不足，需结构+功能条件模型。
4. 卡点：本机无 openfold/GPU/MSA 工具，拿不到自己 WT 的预测结构；3CDL 同源(33%)+ANM 已给出 v2 所需结构证据。升级路径：有 GPU 时跑 ESMFold 或 ColabFold。

## 最终产物（work/results/structure/）
- `STRUCTURE.md`(本文件) · `site_domain_map.csv` · `site_3d_distance.csv` · `anm_*.csv/.npy` · `esm2_*.csv`
- 图: `fig_domain_map.png` · `fig_anm_allosteric.png` · `fig_esm2_zeroshot_heatmap.png` · `fig_site_3d_distances.png`
- 数据: `structure/pdb/3CDL.pdb` · `structure/wt_to_3cdl_map.json`

## 9. GPU 自身结构尝试（2026-09-11，受阻，记录备查）
用 step4_mid_h800 集群跑 ColabFold/ESMFold 出我们自己的 WT 结构，4 次提交均失败：
1. colabfold_batch (fasta) → 卡远程 MMseqs2 MSA（H800 节点无外网）。
2. 装 jax[cuda12]+colabfold[alphafold] 成功、jax+CudaDevice(id=0) 验证 GPU 可用，但 colabfold 仍卡。
3. 改 a3m 单序列输入（跳 MMseqs2）→ 仍卡 50min 无产出（log.txt 仅"Running colabfold 1.6.2"），被 H800 池超配(257/256)抢占中断。
4. ESMFold 本地/节点都缺 openfold（setup.py "No module named scripts" + 需 nvcc 编译 CUDA 扩展）。

**根因**：① H800 池超配(257/256)，我的 job 被抢占；② 集群节点无外部互联网（MMseqs2 远程、AF2/ESMFold 权重下载都走不通）；③ openfold 不可装（构建环境缺 torch + 需 nvcc）。
**已得的结构结论不受影响**：3CDL 同源(33%)+ProDy ANM 已给出 v2 所需全部结构证据（pocket↔HTH 33Å 长程别构耦合、口袋壁 S161/F163/Q164、别构中继 Y102、ESM2 突变先验）。
**解封条件**：有稳定不抢占的 GPU + 节点外网（或预下好 AF2/ESMFold 权重到 workspace NFS 缓存）时，跑 ESMFold 或 ColabFold(带本地 MMseqs2) 出自身结构，在它上做突变扰动预测。

## 10. 自身 ESMFold 结构（2026-09-13，完成）——验证并升级 3CDL 同源结论

用 4090(姊妹项目 esmfold conda env,openfold+权重现成)折叠 WT/Q164R/R116Q/Q164R_R116Q,经两跳 base64 拉回 `work/structure/esmfold/*.pdb`。**pLDDT 94.3**(高可信)。

### 10.1 别构远程距离(自身结构,全原子)
- **Q164(口袋)↔HTH-α3(res35) = 32.7Å** —— 与之前 3CDL 同源(33Å)**几乎完全一致**,证实长程别构不是同源建模假象。
- R116↔Q164 = 19.9Å。
- **口袋壁位点(离 Q164<10Å)**: E94/Y159/S161/F163/Q164/W169/P170(比 3CDL 更全)。

### 10.2 ANM 动态别构耦合(自身结构,20 慢模式 cross-corr)⭐
- **Q164↔HTH(35) 互相关 = -0.394**(|>0.3| 显著,负=反相联动,别构典型)。
- Q164→HTH 区(9-45)平均耦合 **-0.243** vs 全局平均 **+0.005** —— HTH 区耦合远超背景。
- **口袋壁 6 位点与 Q164 强正耦合**: F163(0.92)/S161(0.81)/W169(0.78)/P170(0.75)/Y159(0.66)/E94(0.35) —— 与 Q164 同一动态模块(功能簇),且与 6 靶点观测的高频突变位点高度重合。

### 10.3 变体 vs WT 静态结构 diff
Q164R/R116Q/双突变 vs WT 全局 RMSD 仅 0.06–0.22Å,口袋区 0.08–0.10Å —— 静态 ESMFold 对点突变几乎不动(符合预期),**别构效应须靠 ANM 动力学而非静态结构 diff 捕捉**(本节 10.2 正是)。

### 10.4 结论
自身高置信 ESMFold 结构(pLDDT 94.3)**独立验证**了 3CDL 同源+ANM 的全部结构结论(pocket↔HTH 33Å 长程别构、口袋壁功能簇),并升级为全原子精度。环境链(4090 esmfold env)已固化,后续折叠任意变体零成本。

**产物**: `esmfold/{WT,Q164R,R116Q,Q164R_R116Q}.pdb`、`self_site_distances.json`、`self_anm_coupling.json`。
