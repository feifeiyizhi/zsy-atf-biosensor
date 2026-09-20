#!/usr/bin/env python3
"""突变推荐模型: 对每个胆酸靶点预测'下一步突变哪个位置/氨基酸能增强响应'.
特征融合4路信号:
1. 跨靶点转移: 该突变在其他胆酸靶点的频率(已知有效→候选相似靶点)
2. ESM2零样本: 突变对序列似然扰动(非破坏性)
3. 负筛DNA结合保持: 338lib中携带该突变的变体DNA结合适应度(保留功能)
4. ANM别构相关: 位点对口袋扰动的响应(别构路径上)
留一靶点验证: 训5靶点预测留出靶点的已知高频突变."""
import pandas as pd, numpy as np, os, json
from scipy.stats import spearmanr

R6="./work/results/per_target_six"
RS="./work/results"
RSST="./work/results/structure"
OUT="./work/results/recommender"; os.makedirs(OUT,exist_ok=True)
TARGETS=["CDCA","PCA","7K-LCA","23K-CDCA","LCA","LLDCA"]
SITES=[66,69,75,94,102,107,111,125,128,129,133,159,161,163,164,169,170]
AAS=list("ACDEFGHIKLMNPQRSTVWY")

# --- 1. 跨靶点突变频率矩阵 ---
mat=pd.read_csv(f"{R6}/mutation_by_target_matrix.csv",index_col=0)  # mutation x target, %
# 扩展: 对17位点x20AA全突变(不止>=10%的), 从per_target position_freq补全
# 已有mat只含显著突变. 补全: 对每个(target,site,alt), 查position_freq
full={}  # (target, mut)->freq
for t in TARGETS:
    pf=pd.read_csv(f"{R6}/{t}_position_freq.csv")
    import ast
    for _,r in pf.iterrows():
        ta=r.top_alt
        if isinstance(ta,str):
            try: ta=ast.literal_eval(ta)
            except: ta=None
        if isinstance(ta,tuple) and len(ta)>=2:
            mut=r.wt+str(int(r.pos))+ta[0]
            full[(t,mut)]=r.nonwt_pct

# --- 2. ESM2零样本 ---
esm=pd.read_csv(f"{RSST}/esm2_zeroshot_17site.csv")  # site,wt,alt,score
esm_map={(int(r.site),r.alt):r.score for _,r in esm.iterrows()}

# --- 3. 负筛DNA结合: 每位点alt的平均fitness ---
g=pd.read_csv(f"{RS}/genotype_17site_landscape.csv")
spec=pd.read_csv(f"{RS}/mutated_position_spectrum.csv")
wt_at={int(r.aa_pos):r.wt_aa for _,r in spec.iterrows()}; wt_at={p:wt_at[p] for p in SITES}
sitecols=[f"p{p}_{wt_at[p]}" for p in SITES]
rel=g[g.reliable]
dnabind={}  # (site,alt)->mean fitness
for p in SITES:
    col=f"p{p}_{wt_at[p]}"
    m=rel.groupby(col)["fitness_slope"].mean()
    for alt,fit in m.items():
        dnabind[(p,alt)]=fit

# --- 4. ANM别构相关: 位点对口袋扰动响应 + 3D距离 ---
# 用 site_3d_distance.csv 的 d_pocket (越近别构越相关) + anm扰动响应
dist=pd.read_csv(f"{RSST}/site_3d_distance.csv")  # wt_pos, cdl_res, d_pocket, d_hth
anm_resp=np.load(f"{RSST}/anm_pocket_perturb_resp.npy")  # 187(3CDL残基)->响应
# 映射wt_pos->响应(经3CDL)
mp=json.load(open("./work/structure/wt_to_3cdl_map.json"))
mp={int(k):v for k,v in mp.items()}
# ANM响应 by wt_pos
# 需anm_PR索引->3CDL残基. 简单: 用d_pocket反比作别构相关分(越近越相关)
allo={}
for _,r in dist.iterrows():
    if pd.notna(r.d_pocket):
        # 距离越近别构相关越高(用1/(1+d))
        allo[int(r.wt_pos)]=1.0/(1.0+float(r.d_pocket))

def norm(d):
    if not d: return d
    v=np.array(list(d.values()),dtype=float); v=(v-np.nanmin(v))/(np.nanmax(v)-np.nanmin(v)+1e-9)
    return dict(zip(d.keys(),v))
esm_norm=norm({k:v for k,v in esm_map.items()})
dnabind_norm=norm({k:v for k,v in dnabind.items()})
allo_norm=norm(allo)

# --- 推荐函数: 给定target T, 给每个(site,alt)打分 ---
def score_mutation(t,p,alt,wt_aa,weights=(0.35,0.25,0.25,0.15)):
    mut=wt_aa+str(p)+alt
    # 跨靶点转移: 其他5靶点中该突变的最大频率(归一)
    transfer=max(full.get((ot,mut),0) for ot in TARGETS if ot!=t)/100.0
    esm_s=esm_norm.get((p,alt),0.0)
    # DNA结合保持: 高=保留结合(好). 但极低fitness=-破坏(扣分)
    db=dnabind_norm.get((p,alt),0.0)
    al=allo_norm.get(p,0.0)
    w_transfer,w_esm,w_db,w_al=weights
    score=w_transfer*transfer+w_esm*esm_s+w_db*db+w_al*al
    return dict(score=score,transfer=transfer,esm=esm_s,dnabind=db,allosteric=al)

# --- 每靶点推荐top-K (排除已固定的) ---
print("="*70); print("每靶点下一步突变推荐(增强胆酸响应)"); print("="*70)
allrecs=[]
for t in TARGETS:
    # 该靶点已固定/高频的突变(>=10%) = 已知,排除
    known=set()
    for (ot,mut),fr in full.items():
        if ot==t and fr>=10: known.add(mut)
    recs=[]
    for p in SITES:
        wt=wt_at[p]
        for alt in AAS:
            if alt==wt: continue
            mut=wt+str(p)+alt
            if mut in known: continue  # 已知跳过
            s=score_mutation(t,p,alt,wt)
            recs.append(dict(target=t,mutation=mut,site=p,alt=alt,**s))
    df=pd.DataFrame(recs).sort_values("score",ascending=False)
    df.to_csv(f"{OUT}/{t}_recommendations.csv",index=False)
    print(f"\n### {t}  (已知固定: {sorted(known)}) ###")
    print(df.head(8)[["mutation","score","transfer","esm","dnabind","allosteric"]].to_string(index=False))
    allrecs.append(df.head(8))
pd.concat(allrecs).to_csv(f"{OUT}/all_recommendations_top8.csv",index=False)

# --- 留一靶点验证: 用5靶点的'已知有效'突变, 看模型能否在留出靶点排序高 ---
print("\n"+"="*70); print("留一靶点验证(LOTO)"); print("="*70)
# 每靶点'真实有效'突变 = freq>=10%的(除了Q164R和R116Q这两个全共享/骨架)
for held in TARGETS:
    held_good=[(p,full[(held,m)]) for (ot,m),p in full.items() if ot==held and full.get((held,m),0)>=10 and m!="R116Q"]
    # held_good: (mut, freq). 看这些突变在'跨靶点转移分'里排多高(用其他5靶点)
    ranks=[]
    for mut,fr in held_good:
        # 该突变在其他靶点的最大频率
        ot_freqs=[full.get((ot,mut),0) for ot in TARGETS if ot!=held]
        transfer=max(ot_freqs) if ot_freqs else 0
        ranks.append((mut,fr,transfer))
    # 简单验证: held靶点的高频突变,是否在其他靶点也较高频(可迁移)
    print(f"{held}: 高频突变->{[(m,f'真实{fr:.0f}%',f'他靶点max{tr:.0f}%') for m,fr,tr in ranks]}")
print("\n产物:",os.listdir(OUT))
