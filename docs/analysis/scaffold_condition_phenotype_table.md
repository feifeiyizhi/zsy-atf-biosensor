# 骨架—条件—表型 对照表(28 克隆 provenance 修正)

> 2026-09-13 · 源数据 `mutation_activity_full.csv` + `ppt_all_data.json`(合作者组会 PPT 0119/0330/0622/0810 提取)+ `work/ref/wt_reference_libbackbone.fasta`(WT 612nt/204aa)+ 各克隆 seq 列实测位点。
> **本表是对 `zero_shot_cross_protein.md` 的 provenance 修正前置**:在重跑任何更大模型前,先把骨架/条件/表型钉死。

## 一、四个被 PPT 实锤的 provenance 问题(每个都直接影响零样本结论)

### 问题 1:"DC"=DCA,是独立第三种胆酸(交叉敏感性靶点,非选择条件)
合作者组会 20260119 PPT 反复出现:"Y159 mutations usually lose their DCA sensitivity"、"Y75...lose DCA"、"R116...remain DCA sensitivity"、"I128...lose DCA";slide 8 明列 **"7k-LCA, CDCA and DCA"**。DCA(脱氧胆酸)是**敏感性测试靶点**,不在 28 克隆的三个选择实验里(那是 7k-LCA/CDCA/23k-CDCA+PCA)。报告里任何 "DC" 指 DCA。

### 问题 2:"23k/PCA" 是两种分别测试,不是合并!csv 合并了
0622/0810 PPT slide 3-4 明写:**"A-D: 23k-CDCA / E-H: PCA"**——同一块 96 孔板,行 A-D=23k-CDCA 条件,行 E-H=PCA 条件。`mutation_activity_full.csv` 的 `exp=0810_23kPCA` 把两者合并成一条。拆开:
- **23k-CDCA**(行 A-D):A3, B1, B6, B12, C1, D1(6 克隆)
- **PCA**(行 E-H):E1, F2, F5, F9, G7, H5(6 克隆)

→ 零样本表里"23k/PCA"那一列实为两个独立条件混算,需拆成 23k-CDCA、PCA 两条。

### 问题 3:三个实验三个不同 origin/骨架,不是一个
PPT 明示:
- **7k-LCA 进化**:origin **Y159M**(0330 slide3 "origin Y159M")。但**序列实测:全部 7 个 7k-LCA 克隆 159=Y=WT,无一带 Y159M** → Y159M 在存活克隆里已 **M159Y 回复**(0119 slide10 "Seq4: R116Q, M159Y" 实锤)。故 28 克隆集里 7k-LCA 与 CDCA **实际共享同一骨架**(WT+R116Q+Q164R),只靠次级位点区分。
- **CDCA 进化**:origin **R116Q**(0330 slide4 "origin R116Q")。
- **23k-CDCA/PCA 进化(0810)**:origin **A11 = R116Q+Q136L+Q164R**(0810 slide4 "Compared to A11 from CDCA with: R116Q;Q136L;Q164R;H180H",后续克隆突变相对 A11 列)。

### 问题 4(最严重):csv 的 mutations 列用了两个不同参考系
- **7k-LCA / CDCA 克隆**:mutations 列是**相对 WT**(列出 R116Q/Q164R)。
- **0810 克隆**:mutations 列是**相对 A11**(A11 的 R116Q+Q136L+Q164R 不计,故 A3="none"、B6="E143D")。

序列实测验证:12 个 0810 克隆 seq 在 116/136/164 位 = Q/L/R(确带 A11 骨架),但 csv mutations 写 "none"/单突变。**同一张 csv 的 mutations 列在三个实验用了两个参考系** → 凡从 mutations/n_mut 列算的特征(位点 onehot、别构扰动位点集),0810 克隆被系统性少算 3 个 origin 突变 → 特征偏差。

## 二、统一(WT 相对)骨架—条件—表型对照表

> 0810 克隆已补回 A11 骨架(R116Q+Q136L+Q164R),与 7k-LCA/CDCA 同为 WT 相对。fc=fold-change=对该克隆**选择胆酸**的响应(来自 PPT slide3=7k-LCA 响应、slide4=CDCA 响应;0810 来自 0810 slide3)。

| clone | 条件 | origin/骨架 | 统一(WT相对)突变 | n_uni | fc | 结构组 |
|---|---|---|---|---|---|---|
| B9 | 7k-LCA | Y159M(已回复)+R116Q+Q164R | Q164R;R116Q | 2 | 23.81 | 7位 |
| B5 | 7k-LCA | 同上 | D176N;Q164R;R116Q | 3 | 23.11 | 7位 |
| A5 | 7k-LCA | 同上 | Q164R;R116Q | 2 | 16.76 | 7位 |
| C9 | 7k-LCA | 同上 | D176N;H174Y;Q164R;R116Q | 4 | 10.91 | 7位 |
| B2 | 7k-LCA | 同上 | Q164R;R116Q | 2 | 9.70 | 7位 |
| A8 | 7k-LCA | 同上 | A121T;D176N;Q164R;R116Q | 4 | 3.55 | 7位 |
| C11 | 7k-LCA | 同上 | E113K;Q164R;R116Q | 3 | 3.43 | 7位 |
| E9 | CDCA | R116Q+Q164R | Q164R;R116Q | 2 | 10.33 | 基准 |
| A11 | CDCA | R116Q+Q164R | Q164R;R116Q | 2 | 9.95 | 基准 |
| E10 | CDCA | R116Q+Q164R | I73N;Q164R;R116Q | 3 | 9.61 | 基准 |
| F8 | CDCA | R116Q+Q164R | I73N;R116Q | 2 | 8.84 | 基准 |
| F11 | CDCA | R116Q+Q164R | R116Q;Y129H | 2 | 8.68 | 基准 |
| B4 | CDCA | R116Q+Q164R | Q164R;R116Q | 2 | 8.52 | 基准 |
| B3 | CDCA | R116Q+Q164R | Q164R;R116Q | 2 | 8.44 | 基准 |
| A1 | CDCA | R116Q+Q164R | Q136L;Q164R;R116Q | 3 | 3.85 | 基准(CDCA里唯一带Q136L) |
| E2 | CDCA | R116Q+Q164R | R116Q | 1 | 5.48 | 基准 |
| F9 | PCA | A11(R116Q+Q136L+Q164R) | R116Q;Q136L;Q164R;L177S | 4 | 5.04 | C23 |
| B1 | 23k-CDCA | A11 | R116Q;Q136L;Q164R | 3 | 4.52 | C23 |
| D1 | 23k-CDCA | A11 | R116Q;Q136L;Q164R;E143D | 4 | 4.49 | C23 |
| H5 | PCA | A11 | R116Q;Q136L;Q164R | 3 | 4.48 | C23 |
| F2 | PCA | A11 | R116Q;Q136L;Q164R | 3 | 4.48 | C23 |
| C1 | 23k-CDCA | A11 | R116Q;Q136L;Q164R | 3 | 4.45 | C23 |
| B6 | 23k-CDCA | A11 | R116Q;Q136L;Q164R;E143D | 4 | 4.11 | C23 |
| A3 | 23k-CDCA | A11 | R116Q;Q136L;Q164R | 3 | 4.10 | C23 |
| F5 | PCA | A11 | R116Q;Q136L;Q164R | 3 | 3.97 | C23 |
| B12 | 23k-CDCA | A11 | R116Q;Q136L;Q164R | 3 | 3.94 | C23 |
| E1 | PCA | A11 | R116Q;Q136L;Q164R | 3 | 3.88 | C23 |
| G7 | PCA | A11 | R116Q;Q136L;Q164R | 3 | 3.83 | C23 |

### 统一后按条件汇总(正确参考系)

| 条件 | n | n_mut(WT相对)均值/范围 | fc 均值 | 骨架 |
|---|---|---|---|---|
| 7k-LCA | 7 | 2.86 / 2–4 | ~13.0 | R116Q+Q164R(Y159M 已回复) |
| CDCA | 9 | 2.11 / 1–3 | ~7.3 | R116Q+Q164R |
| 23k-CDCA | 6 | 3.33 / 3–4 | ~4.2 | A11(R116Q+Q136L+Q164R) |
| PCA | 6 | 3.17 / 3–4 | ~4.1 | A11 |

**两点**:
1. **WT 相对突变数梯度**:CDCA(2.11) < 7k-LCA(2.86) < PCA(3.17) < 23k-CDCA(3.33)。23k-CDCA 突变**最多**(因 A11 origin 自带 3 个),不是最少——纠正了之前"去骨架后 0.25 最少"的另一种问法(那问的是"本轮新增")。
2. **A11 origin(Q136L)与低响应强绑定**:12 个 A11 克隆 fc 全挤 3.83–5.04 紧簇;R116Q+Q164R 骨架克隆 fc 跨 3.43–23.81。CDCA 里唯一带 Q136L 的 A1 fc=3.85(也是 CDCA 最低档)。→ **Q136L/A11 系统性拉低响应**是 23k-CDCA/PCA 跨实验"退步"的实测落点,而非方法失败。

## 三、两个待查红旗(需师兄/PPT 原图确认)

1. **slide3/slide4 fold-change 矩阵在 JSON 里逐数字相同**(`0330_7kLCA_slide3`≡`0330_CDCA_slide3`、slide4 亦然)。推测:slide3=全体克隆对 7k-LCA 的响应、slide4=对 CDCA 的响应(同板克隆测两胆酸),csv 给 7k-LCA 克隆取 slide3、CDCA 克隆取 slide4(各自选择胆酸的响应)。但需确认这两 slide 的真实含义 + 原始荧光(分子/分母),JSON 里**只有 fold-change 比值,无原始荧光**。
2. **零样本特征构建脚本全盘找不到**——`clone_esm2_3B_emb.npy`(ESM2 3B 嵌入,seq 算,OK)+ allosteric_score/onehot 特征是在交互会话里建的、**没存盘** → 零样本结果(ρ=-0.37→+0.37)当前**不可从盘上复现**。若 allosteric_score/onehot 用了 csv 的 mutations 列(而非 seq),0810 克隆特征被少算 3 突变 → "23k/PCA 退步"可能部分是**特征偏差**而非纯 origin 偏移。**必须重建特征脚本、从 seq 重算位点,才能区分这两种解释**。

## 四、需师兄/张老师确认清单(盘上查不到,不编)

| 项 | 用途 | 现状 |
|---|---|---|
| 28 克隆原始荧光(fold-change 分子/分母) | 表型原始值、算误差棒 | PPT 只给比值 |
| 筛选方案:投入库大小、各轮±配体计数、FACS 分选方向/门槛、剂量、重复 | §四混杂审计、证明"胆酸确实改变输出" | PROJECT_BACKGROUND_PROTOCOL.md 有 338lib 负筛 protocol,但 6 胆酸正筛的 FACS 门槛/剂量/重复未在盘 |
| slide3/slide4 真实含义 | 确认 fc 轴=选择胆酸响应 | 推测,未确认 |
| 各克隆来源轮次(round) | "第二轮前已存在突变"标注 | 群体测序无 per-clone 轮次;TGS 逐轮测序是 PPT future plan(未做) |
| 6 胆酸 SMILES 立体(差 23K-CDCA+PCA) | 配体特征 | 4/6 已查 |

## 五、对零样本结论的方法学影响(诚实)

- **翻转信号本身**(7k-LCA↔CDCA -0.37→+0.37)**可能仍成立**:ESM2 嵌入来自 seq(全序列,参考系一致),翻转主要靠序列+结构协同;7k-LCA/CDCA 共享骨架,翻转比较的是次级位点驱动 → 物理因果链未被动摇。
- **但"23k/PCA 退步=origin 偏移"这个解释有混合**:可能是 origin 偏移(数据层)+ 特征少算(特征层)双重作用。拆开 23k-CDCA/PCA 后,6+6 样本太少,单条更不可靠。
- **结论**:零样本 ρ 表当前是 **PRELIMINARY(本就如此标注)**;在重建特征脚本(从 seq 重算位点,统一参考系)+ 拆 23k/PCA 重跑前,任何 ρ 数字(含翻转)的精确值不应被当终值引用。

## 六、下一步执行序(接你给的方向)

1. **重建零样本特征脚本并存盘**——从 seq 算突变位点(统一 WT 相对),不复用 csv mutations 列;存 `work/zero_shot_features.py` + 重跑 → 区分"origin 偏移"vs"特征偏差"对 23k 退步各贡献多少。(不需 GPU,我可直接做)
2. **拆 23k-CDCA / PCA 重跑跨实验**——4 条件(7k-LCA/CDCA/23k-CDCA/PCA)交叉,看翻转是否只属 7k-LCA↔CDCA。(不需 GPU)
3. **向师兄索取**:28 克隆原始荧光 + 6 胆酸 FACS 方案(门槛/剂量/重复)+ slide3/4 含义 + per-clone 轮次(若有)。(需师兄)
4. **同骨架配对交叉响应实验设计**(task5,见下)。

## 七、同骨架配对交叉响应实验设计(task5 草案)

**目的**:区分"只是 DNA 结合变了"vs"胆酸确实改变了 DNA 结合相关输出"。

**配对选择(同骨架、只差 1 个次级突变)**:
- CDCA 骨架(R116Q+Q164R):E9(Q164R+R116Q, fc10.33) vs E10(+I73N, fc9.61) vs A1(+Q136L, fc3.85)——同骨架,E10 加 7 位组位点 I73、A1 加 Q136L(低响应)。
- 7k-LCA 骨架(R116Q+Q164R):B9(fc23.81) vs A8(+A121T+D176N, fc3.55)——同骨架,加 7 位组位点后 fc 暴跌。
- A11 骨架:B1(A11, fc4.52) vs B6(+E143D, fc4.11) vs F9(+L177S, fc5.04)——C23 组次级。

**每对测**:
1. 无胆酸基线(DNA 结合荧光)
2. 目标胆酸剂量响应(多浓度,出 EC50/Hill)
3. **其他胆酸交叉响应**(7k-LCA/CDCA/DCA/UDCA/PCA/23k-CDCA 全矩阵)

**判据**:
- 若突变只改无胆酸基线 → "只是 DNA 结合变了"(非胆酸特异)。
- 若突变改的是"目标胆酸响应斜率/EC50"且不影响其他胆酸 → 胆酸特异传感器(强 claim)。
- 若广谱改多胆酸 → 广谱检测(另一种设计目标,需明示)。

**产物**:响应矩阵(克隆×胆酸)——这正是你说的"各传感器对各胆酸的响应矩阵";可验"特异 vs 广谱"是哪一种设计目标。

关联:[[zsy-atf-biosensor-landscape]]、`zero_shot_cross_protein.md`(本表修正其 provenance)、`PROJECT_BACKGROUND_PROTOCOL.md`(负筛 protocol)。
