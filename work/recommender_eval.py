#!/usr/bin/env python3
"""推荐器评估 vs 随机对照. 正例=6靶点中任一观测>=5%的突变(已知有效). 评正例时LOTO排除观测它的靶点.
指标: AUC, precision@K, 富集倍数(模型top-K vs 随机top-K)."""
import pandas as pd, numpy as np, os, json, ast
from sklearn.metrics import roc_auc_score
R6="./work/results/per_target_six"
RSST="./work/results/structure"
OUT="./work/results/recommender"
TARGETS=["CDCA","PCA","7K-LCA","23K-CDCA","LCA","LLDCA"]
AAS=list("ACDEFGHIKLMNPQRSTVWY")
WT="MQKKLTRSQQKHLDIINAAKEEFIEFGFLAANMDRITSSAEVSKRTLYRHFESKEVLFESVLTIINDSVNESISYHFDPNKSTEEQLTEIAYKEIDVLYKTYGIALARTIVMEFLRQPEMAKTLIQNIYSIRAITQWFRSAIEAKRLKDADPKLMTDVYVSLFQGLFFWPQVMHLDLEPHGEELSQKIETLTTIFLQSYGVAE"
EXCLUDE={116,164}

# 正例: 6靶点观测>=5%的突变(除Q164R/R116Q)
positives=set(); pos_by_target={t:set() for t in TARGETS}
for t in TARGETS:
    pf=pd.read_csv(f"{R6}/{t}_position_freq.csv")
    for _,r in pf.iterrows():
        ta=r.top_alt
        if isinstance(ta,str):
            try: ta=ast.literal_eval(ta)
            except: ta=None
        if isinstance(ta,tuple) and len(ta)>=2 and r.nonwt_pct>=5:
            p=int(r.pos); a=ta[0]
            if p in EXCLUDE: continue
            mut=f"{WT[p-1]}{p}{a}"; positives.add(mut); pos_by_target[t].add((mut,r.nonwt_pct))
print(f"正例(任一靶点>=5%, 除Q164R/R116Q): {len(positives)}个")
for t in TARGETS: print(f"  {t}: {len(pos_by_target[t])}个")

# ESM2 + 别构(全位点)
esm=pd.read_csv(f"{RSST}/esm2_zeroshot_allpos.csv")
esm_map={(int(r.pos),r.alt):r.score for _,r in esm.iterrows()}
allo=pd.read_csv(f"{RSST}/allosteric_allpos.csv")
allo_map={int(r.wt_pos):r.allo_score for _,r in allo.iterrows()}

# 跨靶点频率: (pos,alt)->{target:freq}
pos_alts={}
for t in TARGETS:
    pf=pd.read_csv(f"{R6}/{t}_position_freq.csv")
    for _,r in pf.iterrows():
        ta=r.top_alt
        if isinstance(ta,str):
            try: ta=ast.literal_eval(ta)
            except: ta=None
        if isinstance(ta,tuple) and len(ta)>=2:
            pos_alts.setdefault((int(r.pos),ta[0]),{})[t]=r.nonwt_pct

def model_score(p,alt,exclude_target=None):
    """LOTO: 评该突变时排除exclude_target的频率(若该突变是被exclude_target观测的正例,则不算它)"""
    wt=WT[p-1]
    # borrow: 其他靶点(排除exclude_target)的最大频率
    ot_freqs=[fr for ot,fr in pos_alts.get((p,alt),{}).items() if ot!=exclude_target]
    borrow=max(ot_freqs)/100.0 if ot_freqs else 0
    esm_s=esm_map.get((p,alt),-10)
    # ESM2归一(该位点内)
    sub=[esm_map.get((p,a),-10) for a in AAS]; en=(esm_s-min(sub))/(max(sub)-min(sub)+1e-9)
    al=allo_map.get(p,0.0)
    return 0.45*borrow+0.30*en+0.25*al

# --- 评估1: AUC (正例 vs 随机负例), LOTO ---
print("\n"+"="*60); print("评估1: AUC (正例 vs 随机负例, LOTO)"); print("="*60)
all_muts=[f"{WT[p-1]}{p}{a}" for p in range(1,len(WT)+1) for a in AAS if a!=WT[p-1] and p not in EXCLUDE]
rng=np.random.default_rng(42); neg_pool=list(set(all_muts)-positives)
auc_runs=[]
for seed in range(10):
    negs=set(rng.choice(neg_pool,size=len(positives)*5,replace=False))
    labels=[]; scores=[]
    for mut in positives|negs:
        p=int(''.join(c for c in mut[1:] if c.isdigit())); alt=mut[-1]
        # 若是正例,LOTO排除观测它的靶点;负例无观测,不排除
        excl=None
        if mut in positives:
            excl=[t for t in TARGETS if mut in {m for m,_ in pos_by_target[t].__iter__()}]
            # 取唯一观测靶点(若多靶点都观测,留一:排除所有观测它的)
        s=model_score(p,alt,excl[0] if excl and len(excl)==1 else None)
        labels.append(1 if mut in positives else 0); scores.append(s)
    auc=roc_auc_score(labels,scores); auc_runs.append(auc)
print(f"模型 AUC = {np.mean(auc_runs):.3f} ± {np.std(auc_runs):.3f} (随机=0.5)")

# --- 评估2: precision@K + 富集 (模型top-K vs 随机top-K) ---
print("\n"+"="*60); print("评估2: precision@K + 富集倍数"); print("="*60)
# 给所有突变LOTO打分(正例排除其观测靶点)
scored=[]
for mut in all_muts:
    p=int(''.join(c for c in mut[1:] if c.isdigit())); alt=mut[-1]
    excl=None
    for t in TARGETS:
        if mut in {m for m,_ in pos_by_target[t]}:
            excl=t; break  # 排除第一个观测靶点(LOTO)
    s=model_score(p,alt,excl)
    scored.append((mut,s,mut in positives))
scored.sort(key=lambda x:-x[1])
prevalence=len(positives)/len(all_muts)
print(f"正例占比(随机期望precision): {prevalence:.4f} ({prevalence*100:.2f}%)")
print(f"{'K':>6} {'模型P@K':>10} {'随机P@K':>10} {'富集倍数':>10}")
for K in [10,20,50,100,200]:
    topk=scored[:K]; model_p=sum(1 for _,_,ispos in topk if ispos)/K
    random_p=prevalence
    print(f"{K:>6} {model_p*100:>9.1f}% {random_p*100:>9.2f}% {model_p/random_p:>9.1f}x")

# --- 评估3: 模型top-K的具体正例(召回哪些) ---
print("\n模型top-50里的正例:")
top50=scored[:50]
hit=[m for m,s,ispos in top50 if ispos]
print(f"  召回正例: {hit} ({len(hit)}/50)")
print(f"  随机top-50期望召回: {prevalence*50:.1f}个")

# 保存
pd.DataFrame(scored,columns=["mutation","score","is_positive"]).to_csv(f"{OUT}/all_scored_loto.csv",index=False)
print("\n产物已存 all_scored_loto.csv")
