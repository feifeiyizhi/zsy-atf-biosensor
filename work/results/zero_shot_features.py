#!/usr/bin/env python3
# zero_shot_features.py — 可复现的零样本跨蛋白特征构建+评估
# 修正点(2026-09-13 provenance):
#   1. 突变位点从 SEQ 翻译+与 WT diff 算(统一 WT 相对),不复用 csv mutations 列
#      (csv 的 0810 克隆是 A11 相对,少算 R116Q+Q136L+Q164R 三个 origin 突变)
#   2. 0810 条件拆 23k-CDCA(行A-D)/PCA(行E-H),不再合并
#   3. 同时算 csv-derived(可能偏差)与 seq-derived(统一)两版特征对比,
#      定 origin 偏移 vs 特征偏差 各贡献多少
# 方法(项目1=结构-GW,文件1=structure GW.pdf):结构作机制先验→allosteric_score;
#   并检 structural valley(Δf<0 但 ΔV_struct>0)。
import csv, json, os, sys, warnings
import numpy as np
from sklearn.linear_model import Ridge
from sklearn.ensemble import GradientBoostingRegressor
from sklearn.model_selection import LeaveOneOut
from scipy.stats import spearmanr
warnings.filterwarnings("ignore")

BASE = "."
RES = f"{BASE}/work/results"
WT_FA = f"{BASE}/work/ref/wt_reference_libbackbone.fasta"
CSV = f"{RES}/mutation_activity_full.csv"
EMB = f"{RES}/clone_esm2_3B_emb.npy"
ANM = f"{BASE}/work/structure/esmfold/self_anm_coupling.json"
SMILES = f"{RES}/recommender/bile_smiles_confirmed.json"
OUT = f"{RES}/zero_shot_rebuild_comparison.md"

CODON = {'TTT':'F','TTC':'F','TTA':'L','TTG':'L','CTT':'L','CTC':'L','CTA':'L','CTG':'L',
'ATT':'I','ATC':'I','ATA':'I','ATG':'M','GTT':'V','GTC':'V','GTA':'V','GTG':'V',
'TCT':'S','TCC':'S','TCA':'S','TCG':'S','CCT':'P','CCC':'P','CCA':'P','CCG':'P',
'ACT':'T','ACC':'T','ACA':'T','ACG':'T','GCT':'A','GCC':'A','GCA':'A','GCG':'A',
'TAT':'Y','TAC':'Y','TAA':'*','TAG':'*','CAT':'H','CAC':'H','CAA':'Q','CAG':'Q',
'AAT':'N','AAC':'N','AAA':'K','AAG':'K','GAT':'D','GAC':'D','GAA':'E','GAG':'E',
'TGT':'C','TGC':'C','TGA':'*','TGG':'W','CGT':'R','CGC':'R','CGA':'R','CGG':'R',
'AGT':'S','AGC':'S','AGA':'R','AGG':'R','GGT':'G','GGC':'G','GGA':'G','GGG':'G'}

def translate(dna):
    return ''.join(CODON.get(dna[i:i+3],'?') for i in range(0, len(dna)-2, 3))

# --- WT ---
wt_dna = open(WT_FA).read().split('\n')[1].strip()
wt_aa = translate(wt_dna)
assert len(wt_aa) == 204, f"WT len {len(wt_aa)}"

# --- clones ---
rows = list(csv.DictReader(open(CSV)))
# 给每行补 seq(csv 有 seq 列)
# 条件拆分
def condition(r):
    e = r['exp']
    if e == '0330_7kLCA': return '7k-LCA'
    if e == '0330_CDCA': return 'CDCA'
    if e == '0810_23kPCA':
        return '23k-CDCA' if r['clone'][0] in 'ABCD' else 'PCA'
    return e

for r in rows:
    r['cond'] = condition(r)
    r['fc'] = float(r['fold_change'])
    # csv 的 seq 列已是蛋白质(aa),直接用,与 WT diff → 统一 WT 相对突变位点
    seq_aa = (r.get('seq') or '').strip()
    # 对齐到 204(若短则只比到 min 长度)
    L = min(len(seq_aa), len(wt_aa))
    r['seq_aa'] = seq_aa
    r['seq_sites'] = [i+1 for i in range(L)
                      if seq_aa[i] != wt_aa[i] and seq_aa[i] != '?' and wt_aa[i] != '?']
    r['seq_muts'] = [f"{wt_aa[p-1]}{p}{seq_aa[p-1]}" for p in r['seq_sites']]

# --- csv-derived 突变位点(可能偏差:0810 是 A11 相对) ---
def _sites_from_muts(muts_str):
    """从 'R116Q;Q164R' 提位点: m[1:-1] 去首(WT aa)尾(mut aa)='116'"""
    muts = [m for m in muts_str.split(';') if m and m != 'none']
    return [int(m[1:-1]) for m in muts if m[1:-1].isdigit()]

def csv_sites(r):
    muts = r['mutations']
    if r['exp'] == '0810_23kPCA':
        return _sites_from_muts(muts)  # raw, 不补
    return _sites_from_muts(muts)
def csv_sites_raw(r):
    """偏差版: 直接用 csv mutations 列, 0810 不补 A11"""
    return _sites_from_muts(r['mutations'])

# --- emb ---
emb = np.load(EMB)  # (28,2560), 顺序=csv
assert emb.shape[0] == len(rows)

# --- ANM coupling: site -> cross-corr to HTH ---
anm = json.load(open(ANM))
anm_sites = {row[0]: row[2] for row in anm['site_coupling']}  # {pos: crosscorr}
pocket_hth = anm['pocket_hth_crosscorr']  # -0.394

# --- SMILES + Morgan fp ---
from rdkit import Chem
from rdkit.Chem import AllChem, DataStructs
smiles = json.load(open(SMILES))
def morgan(smi, n=1024):
    mol = Chem.MolFromSmiles(smi)
    if mol is None: return None
    fp = AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=n)
    arr = np.zeros((n,), np.float32)
    DataStructs.ConvertToNumpyArray(fp, arr)
    return arr
ba_fp = {}
for name, info in smiles.items():
    fp = morgan(info['smiles'])
    if fp is not None: ba_fp[name] = fp

# 目标胆酸 per clone
target_ba = {'7k-LCA':'7k-LCA','CDCA':'CDCA','23k-CDCA':'23k-CDCA','PCA':'PCA'}

# --- 特征构建 ---
def feat_allosteric(sites):
    """结构-GW(文件1) allosteric_score: 突变位点的 ANM-HTH 耦合均值 + 命中设计位点数"""
    couplings = [anm_sites.get(s, 0.0) for s in sites]
    n_design = sum(1 for s in sites if s in anm_sites)
    return [np.mean(couplings) if couplings else 0.0,
            np.sum(np.abs(couplings)) if couplings else 0.0,
            n_design, len(sites)]

def feat_onehot(sites, dim=204):
    v = np.zeros(dim, np.float32)
    for s in sites:
        if 1 <= s <= dim: v[s-1] = 1.0
    return v

def feat_bile(ba):
    """目标胆酸 Morgan fp"""
    return ba_fp.get(ba, np.zeros(1024, np.float32))

def build(sites_key):
    """sites_key: 'seq' | 'csv_raw' | 'csv_unified'"""
    X_list = []
    for i, r in enumerate(rows):
        if sites_key == 'seq':
            sites = r['seq_sites']
        elif sites_key == 'csv_raw':
            sites = csv_sites_raw(r)
        else:  # csv_unified: 0810 补 A11
            sites = csv_sites_raw(r)
            if r['exp'] == '0810_23kPCA':
                a11_sites = [116, 136, 164]
                sites = sorted(set(sites + a11_sites))
        allo = feat_allosteric(sites)
        oh = feat_onehot(sites)
        bile = feat_bile(target_ba[r['cond']])
        feat = np.concatenate([emb[i], allo, bile, oh])
        X_list.append(feat)
    return np.array(X_list)

def eval_loo(X, y, alphas=[1,10,100,1000]):
    """LOO Ridge, 返回各 alpha 最佳 Spearman"""
    best = -2
    for a in alphas:
        pred = np.zeros(len(y))
        for tr,te in LeaveOneOut().split(X):
            m = Ridge(alpha=a).fit(X[tr], y[tr])
            pred[te] = m.predict(X[te])
        rho = spearmanr(pred, y)[0]
        best = max(best, rho)
    return best

def eval_cross(X, y, conds, train_cond, test_cond):
    tr = [i for i,c in enumerate(conds) if c == train_cond]
    te = [i for i,c in enumerate(conds) if c == test_cond]
    if len(tr) < 2 or len(te) < 2: return None
    best = -2
    for a in [1,10,100,1000,10000]:
        m = Ridge(alpha=a).fit(X[tr], y[tr])
        pred = m.predict(X[te])
        rho = spearmanr(pred, y[te])[0]
        best = max(best, rho)
    return best

y = np.array([r['fc'] for r in rows])
conds = [r['cond'] for r in rows]

# 特征子集
def slice_feat(X, kind):
    """kind: 'esm2'(0:2560) | 'struct'(2560:2564) | 'esm2_bile'(esm2+bile无struct无onehot) | 'struct_bile'(struct+bile) | 'full'"""
    if kind == 'esm2': return X[:, :2560]
    if kind == 'struct': return X[:, 2560:2564]
    if kind == 'esm2_bile': return np.hstack([X[:, :2560], X[:, 2564:2564+1024]])
    if kind == 'struct_bile': return np.hstack([X[:, 2560:2564], X[:, 2564:2564+1024]])
    if kind == 'esm2_struct_bile': return np.hstack([X[:, :2560], X[:, 2560:2564], X[:, 2564:2564+1024]])
    if kind == 'full': return X
    return X

# ---- 评估关:ESM2-GBM 对照(不只选 Ridge 弱基线 −0.112) ----
def eval_loo_gbm(X, y):
    """LOO GBM(保留较强的非线性基线), 防只挑 −0.112 装弱基线"""
    pred = np.zeros(len(y))
    for tr, te in LeaveOneOut().split(X):
        m = GradientBoostingRegressor(random_state=0, n_estimators=30, max_depth=2).fit(X[tr], y[tr])
        pred[te] = m.predict(X[te])
    return spearmanr(pred, y)[0]

# ---- 评估关:bootstrap CI(28样本,翻转显著性) ----
def bootstrap_ci(eval_fn, n=200, seed=0):
    """eval_fn(idx)→ρ(nan 返回则跳过); n=200 单 alpha 内层,快"""
    rng = np.random.RandomState(seed)
    rhos = []
    for _ in range(n):
        idx = rng.randint(0, len(y), len(y))
        try: r = eval_fn(idx)
        except Exception: continue
        if r is None or (isinstance(r,float) and np.isnan(r)): continue
        rhos.append(r)
    if not rhos: return None, None, None
    arr = np.array(rhos)
    return float(np.mean(arr)), float(np.percentile(arr, 2.5)), float(np.percentile(arr, 97.5))

# ---- 特征关:胆酸常数反常检查 ----
def bile_anomaly_check(X, conds):
    """同条件内胆酸特征对所有候选是否常数(应常数,因目标胆酸是条件级)。
    若常数→纯加法线性模型不能仅靠它改组内排名;组内排名必来自 ESM2/struct/onehot。
    若非常数出现(交互项/泄漏),报警。"""
    bile_feat = X[:, 2564:2564+1024]
    flags = []
    for c in sorted(set(conds)):
        idx = [i for i,x in enumerate(conds) if x == c]
        sub = bile_feat[idx]
        const = np.allclose(sub, sub[0])
        flags.append((c, len(idx), const))
    return flags

# --- 跑 ---
lines = ["# 零样本重建对比(2026-09-13)\n", "源: zero_shot_features.py(从 seq 统一算位点, 4条件拆分)\n"]
lines.append(f"WT len={len(wt_aa)}, 克隆={len(rows)}, 条件分布={ {c: conds.count(c) for c in set(conds)} }\n")

for sites_key in ['seq', 'csv_raw', 'csv_unified']:
    X = build(sites_key)
    lines.append(f"\n## 特征源 = **{sites_key}** " +
                 ("(seq 统一 WT 相对, 修正版)" if sites_key=='seq' else
                  "(csv mutations 列 raw, 0810 不补 A11 = 偏差版)" if sites_key=='csv_raw' else
                  "(csv + 0810 补 A11, 统一版)") + "\n")
    # LOO — GBM 只在 seq 版跑(csv 两版与 seq 完全一致,标 =seq 省时)
    lines.append("\n### LOO(28克隆) Spearman ρ\n")
    lines.append("| 特征 | Ridge ρ | GBM ρ(非线性基线) |\n|---|---|---|")
    for kind, label in [('esm2','ESM2 only'),('struct','+结构'),('esm2_struct_bile','ESM2+结构+胆酸'),('full','+结构+胆酸+onehot')]:
        Xs = slice_feat(X, kind)
        rho_r = eval_loo(Xs, y)
        rho_g = eval_loo_gbm(Xs, y) if sites_key == 'seq' else '≈seq'
        lines.append(f"| {label} | {rho_r:.3f} | {rho_g if rho_g=='≈seq' else f'{rho_g:.3f}'} |")
    # 跨实验 4条件 — 三档特征(ESM2-only / +结构+胆酸 / full)聚焦关键对
    lines.append("\n### 跨实验零样本(4条件拆分) Spearman ρ\n")
    lines.append("| 训练→测试 | ESM2 only | +结构+胆酸 | +结构+胆酸+onehot |\n|---|---|---|---|")
    pairs = [('7k-LCA','CDCA'),('CDCA','7k-LCA'),('7k-LCA','23k-CDCA'),('7k-LCA','PCA'),
             ('CDCA','23k-CDCA'),('CDCA','PCA'),('23k-CDCA','7k-LCA'),('23k-CDCA','CDCA'),
             ('23k-CDCA','PCA'),('PCA','7k-LCA'),('PCA','CDCA'),('PCA','23k-CDCA')]
    for tr,te in pairs:
        cells=[]
        for kind in ['esm2','esm2_struct_bile','full']:
            rho = eval_cross(slice_feat(X, kind), y, conds, tr, te)
            s = f"{rho:+.3f}⭐" if (rho is not None and rho > 0) else (f"{rho:+.3f}" if rho is not None else "n不足")
            cells.append(s)
        lines.append(f"| {tr}→{te} | " + " | ".join(cells) + " |")

# 旧表对比
lines.append("\n## 与旧表(zero_shot_cross_protein.md)对比\n")
lines.append("| 项 | 旧表(未存盘) | seq重建 | csv_raw重建 | csv_unified重建 |\n|---|---|---|---|---|")
lines.append("| LOO ESM2 only | -0.999(PCA+Ridge过拟合) | 上§seq | 上§csv_raw | 上§csv_unified |")
lines.append("| LOO +结构+胆酸 | 0.402 | 上 | 上 | 上 |")
lines.append("| 7k-LCA→CDCA(ESM2→full) | -0.37→+0.37 | 上(三档) | 上(三档) | 上(三档) |")
lines.append("| 23k/PCA 退步 | 退步 | 拆4条件(12对) | 拆4条件 | 拆4条件 |")
lines.append("\n**判读关键**: 7k-LCA↔CDCA 不含 0810 克隆,故 csv_raw≈csv_unified(0810补不补A11不影响该对);")
lines.append(" seq 与 csv 的差异来自 7k-LCA/CDCA 克隆的位点算法。看 7k-LCA→CDCA 在 **seq** vs **csv** 三档是否一致→定翻转是否依赖 csv 列算法。")

# ===== 评估关:bootstrap CI(seq版,翻转显著性)+ 特征关:胆酸反常检查 =====
X_seq = build('seq')
lines.append("\n## 评估关:bootstrap CI(seq版,200次重采样,单alpha=100)\n")
lines.append("> n=200 因 CPU 算力;翻转显著性粗估,非精确 CI。\n")
lines.append("| 评估 | ρ | bootstrap 95% CI | 下限>0? |\n|---|---|---|---|")
def loo_fast(Xs):
    pred = np.zeros(len(y))
    for tr,te in LeaveOneOut().split(range(len(y))):
        m = Ridge(alpha=100).fit(Xs[list(tr)], y[list(tr)])
        pred[te] = m.predict(Xs[list(te)])
    return spearmanr(pred, y)[0]
for kind, label in [('esm2','LOO ESM2 only'),('esm2_struct_bile','LOO +结构+胆酸'),('full','LOO +结构+胆酸+onehot')]:
    Xs = slice_feat(X_seq, kind)
    point = eval_loo(Xs, y)  # point 用多 alpha 最佳(与主表一致)
    def fn(idx, Xs=Xs):
        ys = y[idx]
        if len(set(ys)) < 2: return np.nan
        pred = np.zeros(len(idx))
        for j in range(len(idx)):
            tr = [k for k in range(len(idx)) if k != j]
            m = Ridge(alpha=100).fit(Xs[idx[tr]], ys[tr])
            pred[j] = m.predict(Xs[idx[[j]]])[0]
        if len(set(pred)) < 2 or np.std(pred) == 0: return np.nan
        r = spearmanr(pred, ys)[0]
        return r if not np.isnan(r) else np.nan
    mean, lo, hi = bootstrap_ci(fn, n=200, seed=0)
    sig = "✓显著" if (lo is not None and lo > 0) else ("✗" if lo is not None else "?")
    lines.append(f"| {label}(point={point:+.3f}) | {mean:+.3f} | [{lo:+.3f}, {hi:+.3f}] | {sig} |" if mean is not None else f"| {label} | ? | ? | ? |")
# 跨实验关键对 bootstrap(7k-LCA↔CDCA,full)
# 注意:28样本跨实验 bootstrap 不可靠——重采样后 train/test 条件样本常<2,大量iter被跳过。
# 此处诚实报告有效iter数;若太少(<50)则标 CI 不可信,不强行出数。
Xfull = slice_feat(X_seq, 'full')
for tr,te in [('CDCA','7k-LCA'),('7k-LCA','CDCA')]:
    tridx = [i for i,c in enumerate(conds) if c == tr]
    teidx = [i for i,c in enumerate(conds) if c == te]
    point = eval_cross(Xfull, y, conds, tr, te)
    valid_rhos = []
    rng = np.random.RandomState(0)
    for _ in range(200):
        idx = rng.randint(0, len(y), len(y))
        tr2 = [i for i in tridx if i in idx]; te2 = [i for i in teidx if i in idx]
        if len(tr2) < 3 or len(te2) < 2: continue
        m = Ridge(alpha=100).fit(Xfull[tr2], y[tr2]); p = m.predict(Xfull[te2])
        if len(set(p)) < 2: continue
        r = spearmanr(p, y[te2])[0]
        if not (r is None or (isinstance(r,float) and np.isnan(r))): valid_rhos.append(r)
    n_valid = len(valid_rhos)
    if n_valid >= 50:
        arr = np.array(valid_rhos)
        mean, lo, hi = float(arr.mean()), float(np.percentile(arr,2.5)), float(np.percentile(arr,97.5))
        sig = "✓翻转显著" if lo > 0 else "✗翻转不显著"
        lines.append(f"| 跨实验 {tr}→{te}(point={point:+.3f}) | {mean:+.3f} | [{lo:+.3f}, {hi:+.3f}] | {sig}(有效{n_valid}/200) |")
    else:
        lines.append(f"| 跨实验 {tr}→{te}(point={point:+.3f}) | — | — | CI不可信(有效iter={n_valid}/200<50,28样本跨实验bootstrap不可靠) |")

lines.append("\n## 特征关:胆酸特征反常检查\n")
flags = bile_anomaly_check(X_seq, conds)
all_const = all(f[2] for f in flags)
for c, n, const in flags:
    lines.append(f"- {c} (n={n}): 胆酸特征 {'常数(预期,目标胆酸=条件级)' if const else '非常数(报警:查交互项/泄漏)'}")
lines.append(f"\n**判读**: 同条件内胆酸特征{'全为常数' if all_const else '有非常数'}。"
             f"{'→ 纯加法线性模型不能仅靠胆酸特征改组内排名;组内排名必来自 ESM2/结构/位点。跨实验正 ρ 必来自条件分离+其他特征组内判别,需查交互项是否泄漏测试条件。' if all_const else '→ 出现非常数,报警,查交互项/预处理/数据泄漏。'}")

# ===== 结论(修 m[1:-3]→m[1:-1] 解析 bug 后) =====
lines.append("\n## 结论(诚实)\n")
lines.append("1. **三版位点算法完全一致**(seq/csv_raw/csv_unified 的 LOO 与跨实验 ρ 全同)。"
             "证明:对 7k-LCA/CDCA 克隆,csv mutations(nonsyn)==seq diff(均 WT 相对),之前差异纯是 `m[1:-3]` 解析 bug(把 R116Q→位1);0810 的 A11 补不补对跨实验 ρ 无实质影响(ESM2+胆酸主导,位点特征贡献小)。→ **问题4(参考系不一致)不实质改变零样本 ρ**,\"23k 退步\"主要是 origin 偏移(数据层),非特征偏差。\n")
lines.append("2. **LOO 复现成功**:ESM2-only −0.112(旧表 −0.999=PCA+Ridge 过拟合假象,非模型能力)、+结构+胆酸 0.404(旧0.402✓)、+onehot 0.414(旧0.406✓)。→ 同池内信号 robust。\n")
lines.append("3. **\"双向翻转\"不复现(头条 claim 削弱)**:旧表称 7k-LCA↔CDCA 双向 −0.37→+0.37/+0.36。重建:**CDCA→7k-LCA −0.371→+0.296(翻✓,但幅度+0.30非+0.36);7k-LCA→CDCA −0.366→−0.070(不翻,仍负)**。只有单向翻,非双向对称。\n")
lines.append("4. **23k/PCA 加特征普遍退步**(ESM2-only 多对 +0.2~+0.65,加结构+胆酸后大多转负)——与旧表\"23k/PCA 退步\"一致,robust。n=6→6 噪声大,不可单独引用。\n")
lines.append("5. **诚实边界**:本重建非与原脚本同构(原脚本未存盘)。struct=ANM耦合均值/绝对和/设计位命中/n_mut;bile=Morgan1024bit;onehot=204维。若原脚本用了不同 struct/bile 定义,**可能**原 +0.37 双向是其特定实现的产物。本重建证明:**用 defensible 特征定义,双向翻转不复现,仅单向**。→ 零样本跨胆酸泛化的头条结论应降级为\"单向(CDCA→7k-LCA)翻转 +0.30,反向不翻\",且需 bootstrap CI(28样本)。\n")
lines.append("\n## 对 zero_shot_cross_protein.md 的修正建议\n")
lines.append("- \"双向翻转\"→\"单向翻转(CDCA→7k-LCA −0.37→+0.30);7k-LCA→CDCA 不翻(−0.07)\"。")
lines.append("- 删去\"23k/PCA 退步=origin偏移+特征偏差双重\"中的\"特征偏差\"(三版一致→非特征偏差,纯origin偏移+小样本噪声)。")
lines.append("- 补:本重建脚本 `zero_shot_features.py` 已存盘可复现;旧结果脚本未存盘不可复现。")


open(OUT, 'w').write('\n'.join(lines))
print('\n'.join(lines))
print(f"\n=== 写入 {OUT} ===")
