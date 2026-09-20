#!/usr/bin/env python3
"""突变推荐模型v2: 预测'下一步突变哪个位置/AA增强胆酸响应'.
核心洞察: 靶点特异次级突变不跨靶点重复, 但聚集在LBD响应调谐热点位置
(64/73/93/143/152/176/177等). 位置可迁移,AA靶点特异.
推荐= 热点位置(scan已知调谐位置+别构路径位置) x ESM2可容忍AA x 别构相关."""
import pandas as pd, numpy as np, os, json, ast
R6="./work/results/per_target_six"
RS="./work/results"
RSST="./work/results/structure"
OUT="./work/results/recommender"; os.makedirs(OUT,exist_ok=True)
TARGETS=["CDCA","PCA","7K-LCA","23K-CDCA","LCA","LLDCA"]
AAS=list("ACDEFGHIKLMNPQRSTVWY")
WT="MQKKLTRSQQKHLDIINAAKEEFIEFGFLAANMDRITSSAEVSKRTLYRHFESKEVLFESVLTIINDSVNESISYHFDPNKSTEEQLTEIAYKEIDVLYKTYGIALARTIVMEFLRQPEMAKTLIQNIYSIRAITQWFRSAIEAKRLKDADPKLMTDVYVSLFQGLFFWPQVMHLDLEPHGEELSQKIETLTTIFLQSYGVAE"

# 1. 6靶点全突变频率(全位置)
full={}  # (target,pos,alt)->freq
hotspot_counts={}  # pos -> 有多少靶点在该位置出现>=10%突变
for t in TARGETS:
    pf=pd.read_csv(f"{R6}/{t}_position_freq.csv")
    for _,r in pf.iterrows():
        ta=r.top_alt
        if isinstance(ta,str):
            try: ta=ast.literal_eval(ta)
            except: ta=None
        if isinstance(ta,tuple) and len(ta)>=2 and r.nonwt_pct>=5:
            full[(t,int(r.pos),ta[0])]=r.nonwt_pct
            hotspot_counts[int(r.pos)]=hotspot_counts.get(int(r.pos),0)+1
# 热点位置: >=1靶点有次级突变(>=5%)
hotspots=sorted(hotspot_counts.keys())
print("响应调谐热点位置(>=1靶点有>=5%次级突变):",hotspots)
print("  各位置靶点数:",{p:hotspot_counts[p] for p in hotspots})

# 2. ESM2全位点零样本
esm=pd.read_csv(f"{RSST}/esm2_zeroshot_allpos.csv")
esm_map={(int(r.pos),r.alt):r.score for _,r in esm.iterrows()}

# 3. 别构相关(3D距口袋 + ANM). 用site_3d_distance的d_pocket(越近越相关), 对全位置用3CDL映射
dist=pd.read_csv(f"{RSST}/site_3d_distance.csv")
# 只有17位点有3D. 对全位置, 用ANM扰动响应(经3CDL映射)
mp=json.load(open("./work/structure/wt_to_3cdl_map.json"))
mp={int(k):v for k,v in mp.items()}
try:
    pr=np.load(f"{RSST}/anm_pocket_perturb_resp.npy")  # 3CDL残基->响应
except: pr=None
# d_pocket for 17 sites
dpock={int(r.wt_pos):float(r.d_pocket) for _,r in dist.iterrows() if pd.notna(r.d_pocket)}

def allo_score(p):
    # 17位点: 用d_pocket反比. 其他: 若有3CDL映射用ANM响应, 否则0
    if p in dpock: return 1.0/(1.0+dpock[p])
    return 0.0  # 非设计位点无结构3D, 用热点计数代替

def norm_pos_scores(d):
    if not d: return d
    v=np.array(list(d.values()),float); v=(v-v.min())/(v.max()-v.min()+1e-9)
    return dict(zip(d.keys(),v))

# --- 推荐: 对每靶点, 排序(位置,AA) ---
# 评分: hotspot(位置被多少靶点用作调谐) + ESM2容忍 + 别构相关 + (该位置-该靶点已固定则跳过)
print("\n"+"="*70); print("每靶点 下一步突变推荐(增强胆酸响应)"); print("="*70)
allrecs=[]
for t in TARGETS:
    # 该靶点已固定(>=10%)的位置-AA
    known=set()
    for (ot,pos,alt),fr in full.items():
        if ot==t and fr>=10: known.add((pos,alt))
    recs=[]
    for p in range(1,len(WT)+1):
        wt_aa=WT[p-1]
        # 位置热点分: 该位置被多少靶点用作调谐(归一)
        hs=hotspot_counts.get(p,0)
        # 该位置在各靶点出现的alt集合(已知调谐AA, 可借鉴)
        known_alts_at_p=set(alt2 for (ot2,p2,alt2) in full if p2==p)
        for alt in AAS:
            if alt==wt_aa: continue
            if (p,alt) in known: continue  # 已固定跳过
            esm_s=esm_map.get((p,alt),-10)  # ESM2容忍(越高越好)
            al=allo_score(p)
            # 综合: 热点位置(0.4) + ESM2容忍(0.3) + 别构(0.3) + 借鉴已知AA加成
            borrow=1.0 if alt in known_alts_at_p else 0.0
            score=0.40*hs+0.30*(esm_s+10)/10+0.30*al+0.15*borrow
            recs.append(dict(target=t,pos=p,wt=wt_aa,alt=alt,mutation=f"{wt_aa}{p}{alt}",
                score=round(score,3),hotspot_targets=hs,esm=round(esm_s,2),
                allosteric=round(al,2),borrowed_aa=(alt in known_alts_at_p)))
    df=pd.DataFrame(recs).sort_values("score",ascending=False)
    df.to_csv(f"{OUT}/{t}_recommendations_v2.csv",index=False)
    print(f"\n### {t}  (已固定: {sorted(WT[p-1]+str(p)+alt for (pp,alt) in known for p in [pp] if (p,alt) in known)}) ###")
    print(df.head(10)[["mutation","score","hotspot_targets","esm","allosteric","borrowed_aa"]].to_string(index=False))
    allrecs.append(df.head(10))
pd.concat(allrecs).to_csv(f"{OUT}/all_recommendations_v2_top10.csv",index=False)
print("\n产物:",os.listdir(OUT))
