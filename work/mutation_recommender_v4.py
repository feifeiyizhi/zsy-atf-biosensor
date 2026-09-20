#!/usr/bin/env python3
"""推荐器v4: 全位点评测. 两种推荐:
A. 借鉴型(跨靶点转移): 其他靶点验证过的突变借到本靶点.
B. 全新型(别构路径): 没靶点试过、但别构相关+ESM2可容忍的位点的AA.
+ LOTO验证(留一靶点, 看能否召回该靶点的真实次级)."""
import pandas as pd, numpy as np, os, json, ast
from prody import parsePDB, ANM, calcPerturbResponse, LOGGER; LOGGER.verbosity='none'
R6="./work/results/per_target_six"
RSST="./work/results/structure"
OUT="./work/results/recommender"; os.makedirs(OUT,exist_ok=True)
TARGETS=["CDCA","PCA","7K-LCA","23K-CDCA","LCA","LLDCA"]
AAS=list("ACDEFGHIKLMNPQRSTVWY")
WT="MQKKLTRSQQKHLDIINAAKEEFIEFGFLAANMDRITSSAEVSKRTLYRHFESKEVLFESVLTIINDSVNESISYHFDPNKSTEEQLTEIAYKEIDVLYKTYGIALARTIVMEFLRQPEMAKTLIQNIYSIRAITQWFRSAIEAKRLKDADPKLMTDVYVSLFQGLFFWPQVMHLDLEPHGEELSQKIETLTTIFLQSYGVAE"
EXCLUDE={116,164}

# --- 全位点别构分(3CDL 3D距离 + ANM) ---
struct=parsePDB("./work/structure/pdb/3CDL.pdb")
ca=struct.select('protein and name CA and chain A')
coords=ca.getCoords(); resnums=ca.getResnums()
r2c={rn:i for i,rn in enumerate(resnums)}
mp=json.load(open("./work/structure/wt_to_3cdl_map.json")); mp={int(k):v for k,v in mp.items()}
POCKET=mp[164]; c_pocket=coords[r2c[POCKET]]
PR=np.load(f"{RSST}/anm_PR.npy")  # [187,187] 3CDL
allo_all={}  # wt_pos -> (d_pocket, anm_response)
for wp,cdlp in mp.items():
    if cdlp is None or cdlp not in r2c: continue
    d=float(np.linalg.norm(coords[r2c[cdlp]]-c_pocket))
    resp=float(PR[r2c[POCKET], r2c[cdlp]])  # 扰动口袋->该位响应
    allo_all[wp]=(d,resp)
# 归一: 别构相关 = 距近 + 响应高
d_arr=np.array([v[0] for v in allo_all.values()]); r_arr=np.array([v[1] for v in allo_all.values()])
d_norm=1/(1+d_arr); d_norm=(d_norm-d_norm.min())/(d_norm.max()-d_norm.min()+1e-9)
r_norm=(r_arr-r_arr.min())/(r_arr.max()-r_arr.min()+1e-9)
allo_score=dict(zip(allo_all.keys(), 0.5*d_norm+0.5*r_norm))
pd.DataFrame([{"wt_pos":p,"d_pocket":allo_all[p][0],"anm_resp":allo_all[p][1],"allo_score":allo_score[p]} for p in allo_all]).to_csv(f"{RSST}/allosteric_allpos.csv",index=False)

# --- 6靶点突变 + ESM2 ---
mut_by_target={t:{} for t in TARGETS}
pos_alts={}
for t in TARGETS:
    pf=pd.read_csv(f"{R6}/{t}_position_freq.csv")
    for _,r in pf.iterrows():
        ta=r.top_alt
        if isinstance(ta,str):
            try: ta=ast.literal_eval(ta)
            except: ta=None
        if isinstance(ta,tuple) and len(ta)>=2 and r.nonwt_pct>=5:
            mut_by_target[t][(int(r.pos),ta[0])]=r.nonwt_pct
            pos_alts.setdefault((int(r.pos),ta[0]),{})[t]=r.nonwt_pct
all_mutated_pos=set(p for (p,a) in pos_alts)  # 任意靶点突变过的位置
esm=pd.read_csv(f"{RSST}/esm2_t36_3B_zeroshot_allpos.csv")
esm_map={(int(r.pos),r.alt):r.score for _,r in esm.iterrows()}

def esm_norm_at(p):
    sub=[esm_map.get((p,a),-10) for a in AAS]
    arr=np.array(sub); return dict(zip(AAS,(arr-arr.min())/(arr.max()-arr.min()+1e-9)))

print("="*70); print("每靶点 推荐增强响应的下一步突变"); print("="*70)
allrecs=[]
for t in TARGETS:
    known=set(mut_by_target[t].keys())
    recs=[]
    for p in range(1,len(WT)+1):
        if p in EXCLUDE: continue
        wt_aa=WT[p-1]
        en=esm_norm_at(p); al=allo_score.get(p,0.0)
        for alt in AAS:
            if alt==wt_aa or (p,alt) in known: continue
            borrow=max(pos_alts.get((p,alt),{}).get(ot,0) for ot in TARGETS if ot!=t)/100.0 if (p,alt) in pos_alts else 0
            is_novel = p not in all_mutated_pos  # 全新位置(没靶点试过)
            # A借鉴型: borrow为主; B全新型: 别构+ESM2为主
            if borrow>0:
                score=0.50*borrow+0.25*en[alt]+0.25*al; typ="borrowed"
            elif is_novel and al>np.median(list(allo_score.values())):
                score=0.30*al+0.40*en[alt]+0.30*(1/ (1+allo_all.get(p,(0,0))[0])); typ="novel"
            else:
                score=0.25*al+0.50*en[alt]+0.25*(0.1 if p in all_mutated_pos else 0); typ="explore"
            recs.append(dict(target=t,pos=p,mutation=f"{wt_aa}{p}{alt}",score=round(score,3),type=typ,
                borrowed_from=f"{borrow*100:.0f}%" if borrow>0 else "-",esm=round(esm_map.get((p,alt),-10),2),
                allo=round(al,2),d_pocket=round(allo_all.get(p,(0,0))[0],1)))
    df=pd.DataFrame(recs).sort_values("score",ascending=False)
    df.to_csv(f"{OUT}/{t}_rec_v4.csv",index=False)
    known_str=sorted(f"{WT[p-1]}{p}{a}({fr:.0f}%)" for (p,a),fr in mut_by_target[t].items() if p not in EXCLUDE)
    print(f"\n### {t}  已固定次级: {known_str} ###")
    print(df.head(10)[["mutation","score","type","borrowed_from","esm","allo","d_pocket"]].to_string(index=False))
    allrecs.append(df.head(10))
pd.concat(allrecs).to_csv(f"{OUT}/all_rec_v4_top10.csv",index=False)

# --- LOTO验证: 留一靶点, 用其余5靶点的次级, 看能否召回留出靶点的真实次级 ---
print("\n"+"="*70); print("LOTO验证(留一靶点召回真实次级)"); print("="*70)
for held in TARGETS:
    held_muts=[(p,a,fr) for (p,a),fr in mut_by_target[held].items() if p not in EXCLUDE]
    recall=[]
    for p,a,fr in held_muts:
        # 该(pos,alt)在其他5靶点的最大频率(借鉴信号)
        ot_freq=max(pos_alts.get((p,a),{}).get(ot,0) for ot in TARGETS if ot!=held)
        recall.append((f"{WT[p-1]}{p}{a}",fr,ot_freq))
    n_recall=sum(1 for _,_,of in recall if of>=5)
    print(f"{held}: {len(recall)}个真实次级, 可借鉴召回{n_recall}个: "+", ".join(f"{m}(真实{fr:.0f}%/他靶{of:.0f}%)" for m,fr,of in recall))
print("\n产物:",os.listdir(OUT))
