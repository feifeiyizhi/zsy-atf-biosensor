#!/usr/bin/env python3
"""时序分割评估(最严格,零泄露): 训练=CDCA+PCA(09-11), 测试=7K-LCA/23K-CDCA/LCA/LLDCA(09-12).
测试正例=新4靶点观测>=5%的突变, 用旧2靶点(CDCA/PCA)的borrow+ESM2+别构打分.
对照随机. 这是真前向验证(预测新数据)."""
import pandas as pd, numpy as np, os, ast
from sklearn.metrics import roc_auc_score
R6="./work/results/per_target_six"
RSST="./work/results/structure"
OUT="./work/results/recommender"
TRAIN_T=["CDCA","PCA"]; TEST_T=["7K-LCA","23K-CDCA","LCA","LLDCA"]
AAS=list("ACDEFGHIKLMNPQRSTVWY")
WT="MQKKLTRSQQKHLDIINAAKEEFIEFGFLAANMDRITSSAEVSKRTLYRHFESKEVLFESVLTIINDSVNESISYHFDPNKSTEEQLTEIAYKEIDVLYKTYGIALARTIVMEFLRQPEMAKTLIQNIYSIRAITQWFRSAIEAKRLKDADPKLMTDVYVSLFQGLFFWPQVMHLDLEPHGEELSQKIETLTTIFLQSYGVAE"
EXCLUDE={116,164}

def get_muts(targets,thresh=5):
    """返回 {target: {(pos,alt):freq}} 和 全局 (pos,alt)->{target:freq}"""
    by_t={t:{} for t in targets}; pos_alts={}
    for t in targets:
        pf=pd.read_csv(f"{R6}/{t}_position_freq.csv")
        for _,r in pf.iterrows():
            ta=r.top_alt
            if isinstance(ta,str):
                try: ta=ast.literal_eval(ta)
                except: ta=None
            if isinstance(ta,tuple) and len(ta)>=2 and r.nonwt_pct>=thresh:
                p=int(r.pos); a=ta[0]
                if p in EXCLUDE: continue
                by_t[t][(p,a)]=r.nonwt_pct
                pos_alts.setdefault((p,a),{})[t]=r.nonwt_pct
    return by_t,pos_alts

# 训练集(旧2)的borrow信号
_,train_pos_alts=get_muts(TRAIN_T)
# 测试集(新4)的正例
test_by_t,_=get_muts(TEST_T,thresh=5)
test_positives=set()
for t in TEST_T:
    for (p,a),fr in test_by_t[t].items():
        test_positives.add(f"{WT[p-1]}{p}{a}")
print(f"训练靶点(旧): {TRAIN_T}")
print(f"测试靶点(新,刚测): {TEST_T}")
print(f"测试正例(新4靶点>=5%, 除Q164R/R116Q): {len(test_positives)}个")
for t in TEST_T:
    ms=sorted(f"{WT[p-1]}{p}{a}({fr:.0f}%)" for (p,a),fr in test_by_t[t].items())
    print(f"  {t}: {ms}")

# ESM2+别构(全位点)
esm=pd.read_csv(f"{RSST}/esm2_t36_3B_zeroshot_allpos.csv")
esm_map={(int(r.pos),r.alt):r.score for _,r in esm.iterrows()}
allo=pd.read_csv(f"{RSST}/allosteric_allpos.csv")
allo_map={int(r.wt_pos):r.allo_score for _,r in allo.iterrows()}

def score(p,alt):
    # borrow: 仅用训练靶点(旧2)的频率, 测试靶点完全不参与=零泄露
    borrow=max(train_pos_alts.get((p,alt),{}).get(t,0) for t in TRAIN_T)/100.0
    sub=[esm_map.get((p,a),-10) for a in AAS]
    en=(esm_map.get((p,alt),-10)-min(sub))/(max(sub)-min(sub)+1e-9)
    al=allo_map.get(p,0.0)
    return 0.45*borrow+0.30*en+0.25*al, borrow

# 评估: 测试正例 vs 随机负例
all_muts=[f"{WT[p-1]}{p}{a}" for p in range(1,len(WT)+1) for a in AAS if a!=WT[p-1] and p not in EXCLUDE]
rng=np.random.default_rng(42); neg_pool=list(set(all_muts)-test_positives)
auc_runs=[]; prec_runs={}
for seed in range(10):
    negs=set(rng.choice(neg_pool,size=len(test_positives)*5,replace=False))
    labels=[]; scores=[]
    for mut in test_positives|negs:
        p=int(''.join(c for c in mut[1:] if c.isdigit())); alt=mut[-1]
        s,_=score(p,alt); labels.append(1 if mut in test_positives else 0); scores.append(s)
    auc_runs.append(roc_auc_score(labels,scores))
# 全打分排序, precision@K
scored=[]
for mut in all_muts:
    p=int(''.join(c for c in mut[1:] if c.isdigit())); alt=mut[-1]
    s,borrow=score(p,alt); scored.append((mut,s,borrow,mut in test_positives))
scored.sort(key=lambda x:-x[2])
prevalence=len(test_positives)/len(all_muts)
print(f"\n{'='*60}\n时序评估(训练旧2→测试新4, 零泄露)\n{'='*60}")
print(f"模型AUC = {np.mean(auc_runs):.3f} ± {np.std(auc_runs):.3f} (随机=0.5)")
print(f"测试正例占比(随机期望): {prevalence*100:.2f}%")
print(f"{'K':>6} {'模型P@K':>10} {'随机P@K':>10} {'富集':>8}")
for K in [10,20,50,100,200]:
    topk=scored[:K]; mp=sum(1 for _,_,_,ispos in topk if ispos)/K
    print(f"{K:>6} {mp*100:>9.1f}% {prevalence*100:>9.2f}% {mp/prevalence:>7.1f}x")
print(f"\n模型top-50召回的正例: {[m for m,s,b,ispos in scored[:50] if ispos]}")
print(f"模型top-50里borrow>0的(有旧靶点转移信号): {[(m,round(b*100)) for m,s,b,_ in scored[:50] if b>0]}")
# 关键: 测试正例里有多少borrow>0(即可被旧靶点转移预测的)
test_borrow=[(m,score(int(''.join(c for c in m[1:] if c.isdigit())),m[-1])[1]) for m in test_positives]
predictable=[(m,b) for m,b in test_borrow if b>0]
print(f"\n测试正例中可被旧靶点转移预测的(borrow>0): {len(predictable)}/{len(test_positives)}: {predictable}")
pd.DataFrame(scored,columns=["mutation","score","borrow","is_positive"]).to_csv(f"{OUT}/temporal_eval_scored.csv",index=False)
print("\n产物: temporal_eval_scored.csv")
